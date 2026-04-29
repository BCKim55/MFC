!>
!! @file
!! @brief Contains module m_lso_filter

#:include 'macros.fpp'

!> @brief LSO variable-weight Gaussian filter applied to conserved variables at data-save steps.
!!
!! The filter approximates a Gaussian transfer function G(xi) = exp(-sigma^2 * xi^2 / 2) by composing N_passes applications of a
!! 9-point symmetric FIR stencil: H(xi) = a(1) + 2*[a(2)*cos(xi) + a(3)*cos(2*xi) + a(4)*cos(3*xi) + a(5)*cos(4*xi)]
!!
!! Weights a(1:5) vary per pass (computed by Block Coordinate Descent in the Python toolchain). The filter is applied independently
!! in x, then y, then z directions.
module m_lso_filter

    use m_derived_types
    use m_global_parameters
    use m_mpi_common
    use m_constants

    implicit none

    private

    public :: s_initialize_lso_filter_module, s_apply_lso_filter, s_copy_and_apply_lso_filter, s_lso_stride_sample, &
        & s_finalize_lso_filter_module, q_filt_vf, q_filt_ds_vf

    !> Temporary working buffer for one directional filter pass (interior cells only)
    real(wp), allocatable, dimension(:,:,:) :: lso_tmp
    $:GPU_DECLARE(create='[lso_tmp]')

    !> Filtered copy of the conserved variables (only allocated when lso_filter_wrt = T)
    type(scalar_field), allocatable :: q_filt_vf(:)

    !> Stride-sampled (coarsened) copy of q_filt_vf (only allocated when lso_filter_wrt .and. lso_down_sample_factor > 1)
    type(scalar_field), allocatable :: q_filt_ds_vf(:)

contains

    !> Allocate the temporary pass buffer and, when lso_filter_wrt is enabled, the filtered-copy array q_filt_vf with full
    !! ghost-cell bounds.
    impure subroutine s_initialize_lso_filter_module()

        integer :: i

        @:ALLOCATE(lso_tmp(0:m, 0:n, 0:p))

        if (lso_filter_wrt) then
            @:ALLOCATE(q_filt_vf(1:sys_size))
            do i = 1, sys_size
                @:ALLOCATE(q_filt_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                           & idwbuff(3)%beg:idwbuff(3)%end))
            end do
            do i = 1, sys_size
                @:ACC_SETUP_SFs(q_filt_vf(i))
            end do

            if (lso_down_sample_factor > 1) then
                @:ALLOCATE(q_filt_ds_vf(1:sys_size))
                do i = 1, sys_size
                    @:ALLOCATE(q_filt_ds_vf(i)%sf(0:m_lso_ds, 0:n_lso_ds, 0:p_lso_ds))
                end do
            end if
        end if

    end subroutine s_initialize_lso_filter_module

    !> Deallocate the temporary pass buffer and, when allocated, q_filt_vf and q_filt_ds_vf.
    impure subroutine s_finalize_lso_filter_module()

        integer :: i

        @:DEALLOCATE(lso_tmp)

        if (lso_filter_wrt) then
            do i = 1, sys_size
                @:DEALLOCATE(q_filt_vf(i)%sf)
            end do
            @:DEALLOCATE(q_filt_vf)

            if (lso_down_sample_factor > 1) then
                do i = 1, sys_size
                    @:DEALLOCATE(q_filt_ds_vf(i)%sf)
                end do
                @:DEALLOCATE(q_filt_ds_vf)
            end if
        end if

    end subroutine s_finalize_lso_filter_module

    !> Copy q_cons_vf into q_filt_vf (including ghost cells) on the device, then apply the LSO filter to q_filt_vf in place. Used
    !! when lso_filter_wrt = T so the original conserved-variable array is left unmodified for the primary (unfiltered) output
    !! write.
    impure subroutine s_copy_and_apply_lso_filter(q_cons_vf)

        type(scalar_field), intent(in) :: q_cons_vf(:)
        integer                        :: i, j, k, l

        do i = 1, sys_size
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        q_filt_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

        call s_apply_lso_filter(q_filt_vf)

    end subroutine s_copy_and_apply_lso_filter

    !> Apply the LSO variable-weight filter to all conserved variables.
    !!
    !! For each spatial direction with lso_n_passes > 0, the routine iterates over the prescribed number of passes, applying a
    !! 9-point symmetric stencil per pass. Ghost cells populated by the preceding halo exchange supply boundary stencil values. The
    !! filter is applied in-place on the interior (0:m, 0:n, 0:p).
    !!
    !! The outer loop runs over passes and the inner loop over variables so that an MPI halo exchange can be issued between passes,
    !! keeping the ghost cells consistent with the filtered interior of neighbouring ranks. Without this exchange, ghost cells at
    !! MPI rank boundaries remain frozen at pre-filter values for every pass, which causes Gibbs-like ringing near discontinuities
    !! (e.g. IBM ghost cells) that is absent in single-rank serial runs.
    !!
    !! @note  When time_stepper == 1 (forward Euler), stor == 1 and the live conserved
    !! variable array is modified in-place. For all Runge-Kutta time steppers (stor == 2) only the save-copy is affected and the
    !! simulation state is unchanged. IBM cases virtually always use an RK scheme, so this is safe in practice.
    impure subroutine s_apply_lso_filter(q_cons_vf)

        type(scalar_field), intent(inout) :: q_cons_vf(:)
        integer                           :: i, ipass, j, k, l
        real(wp)                          :: c0, c1, c2, c3, c4

        ! x-direction passes: outer loop over passes so all variables are updated before the
        ! inter-pass MPI exchange, which refreshes x-direction ghost cells at rank boundaries.
        do ipass = 1, lso_n_passes_x
            c0 = lso_a_x(1, ipass)
            c1 = lso_a_x(2, ipass)
            c2 = lso_a_x(3, ipass)
            c3 = lso_a_x(4, ipass)
            c4 = lso_a_x(5, ipass)
            do i = 1, sys_size
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            lso_tmp(j, k, l) = c0*real(q_cons_vf(i)%sf(j, k, l), wp) + c1*(real(q_cons_vf(i)%sf(j - 1, k, l), &
                                    & wp) + real(q_cons_vf(i)%sf(j + 1, k, l), wp)) + c2*(real(q_cons_vf(i)%sf(j - 2, k, l), &
                                    & wp) + real(q_cons_vf(i)%sf(j + 2, k, l), wp)) + c3*(real(q_cons_vf(i)%sf(j - 3, k, l), &
                                    & wp) + real(q_cons_vf(i)%sf(j + 3, k, l), wp)) + c4*(real(q_cons_vf(i)%sf(j - 4, k, l), &
                                    & wp) + real(q_cons_vf(i)%sf(j + 4, k, l), wp))
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_cons_vf(i)%sf(j, k, l) = real(lso_tmp(j, k, l), stp)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do
            if (ipass < lso_n_passes_x) call s_lso_filter_ghost_refresh(q_cons_vf, 1)
        end do

        ! y-direction passes (2D/3D only: n > 0)
        ! Exchange y ghost cells before starting y passes so that the ghost cells of each rank
        ! reflect the x-filtered interior values of its y-neighbours.
        if (n > 0) then
            call s_lso_filter_ghost_refresh(q_cons_vf, 2)
            do ipass = 1, lso_n_passes_y
                c0 = lso_a_y(1, ipass)
                c1 = lso_a_y(2, ipass)
                c2 = lso_a_y(3, ipass)
                c3 = lso_a_y(4, ipass)
                c4 = lso_a_y(5, ipass)
                do i = 1, sys_size
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                lso_tmp(j, k, l) = c0*real(q_cons_vf(i)%sf(j, k, l), wp) + c1*(real(q_cons_vf(i)%sf(j, k - 1, l), &
                                        & wp) + real(q_cons_vf(i)%sf(j, k + 1, l), wp)) + c2*(real(q_cons_vf(i)%sf(j, k - 2, l), &
                                        & wp) + real(q_cons_vf(i)%sf(j, k + 2, l), wp)) + c3*(real(q_cons_vf(i)%sf(j, k - 3, l), &
                                        & wp) + real(q_cons_vf(i)%sf(j, k + 3, l), wp)) + c4*(real(q_cons_vf(i)%sf(j, k - 4, l), &
                                        & wp) + real(q_cons_vf(i)%sf(j, k + 4, l), wp))
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                q_cons_vf(i)%sf(j, k, l) = real(lso_tmp(j, k, l), stp)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
                if (ipass < lso_n_passes_y) call s_lso_filter_ghost_refresh(q_cons_vf, 2)
            end do
        end if

        ! z-direction passes (3D only: p > 0)
        if (p > 0) then
            call s_lso_filter_ghost_refresh(q_cons_vf, 3)
            do ipass = 1, lso_n_passes_z
                c0 = lso_a_z(1, ipass)
                c1 = lso_a_z(2, ipass)
                c2 = lso_a_z(3, ipass)
                c3 = lso_a_z(4, ipass)
                c4 = lso_a_z(5, ipass)
                do i = 1, sys_size
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                lso_tmp(j, k, l) = c0*real(q_cons_vf(i)%sf(j, k, l), wp) + c1*(real(q_cons_vf(i)%sf(j, k, l - 1), &
                                        & wp) + real(q_cons_vf(i)%sf(j, k, l + 1), wp)) + c2*(real(q_cons_vf(i)%sf(j, k, l - 2), &
                                        & wp) + real(q_cons_vf(i)%sf(j, k, l + 2), wp)) + c3*(real(q_cons_vf(i)%sf(j, k, l - 3), &
                                        & wp) + real(q_cons_vf(i)%sf(j, k, l + 3), wp)) + c4*(real(q_cons_vf(i)%sf(j, k, l - 4), &
                                        & wp) + real(q_cons_vf(i)%sf(j, k, l + 4), wp))
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                q_cons_vf(i)%sf(j, k, l) = real(lso_tmp(j, k, l), stp)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
                if (ipass < lso_n_passes_z) call s_lso_filter_ghost_refresh(q_cons_vf, 3)
            end do
        end if

    end subroutine s_apply_lso_filter

    !> Stride-sample (decimate) q_src_vf into q_dst_vf with the given integer factor in all active directions.
    !!
    !! The filtered source data is already band-limited by the LSO filter, so a simple stride pick (every Nth
    !! cell, starting at index 0) is a valid decimation with no additional anti-aliasing needed. The destination
    !! array must be allocated to (0:m_lso_ds, 0:n_lso_ds, 0:p_lso_ds) by the caller.
    !!
    !! @param q_src_vf   Source conserved-variable array (full grid, interior cells 0:m, 0:n, 0:p on host).
    !! @param q_dst_vf   Destination array (coarsened, 0:m_lso_ds, etc.).
    !! @param fac        Stride factor (positive integer; 1 = identity copy).
    impure subroutine s_lso_stride_sample(q_src_vf, q_dst_vf, fac)

        type(scalar_field), intent(in)    :: q_src_vf(:)
        type(scalar_field), intent(inout) :: q_dst_vf(:)
        integer, intent(in)               :: fac
        integer                           :: i, j, k, l

        do i = 1, sys_size
            do l = 0, p_lso_ds
                do k = 0, n_lso_ds
                    do j = 0, m_lso_ds
                        q_dst_vf(i)%sf(j, k, l) = q_src_vf(i)%sf(j*fac, k*fac, l*fac)
                    end do
                end do
            end do
        end do

    end subroutine s_lso_stride_sample

    !> Refresh ghost cells for a given spatial direction between LSO filter passes.
    !!
    !! Handles two kinds of ghost cell sources:
    !!
    !! - MPI boundaries (bc >= 0): issues a sendrecv with the neighbouring rank so the ghost
    !!   cells reflect the latest filtered interior values of that rank.
    !!
    !! - Physical boundaries (bc < 0): re-applies the boundary ghost cell rule using the
    !!   current (filtered) state of the interior edge cells.  Currently implemented for
    !!   BC_GHOST_EXTRAP (-3), which is a zeroth-order extrapolation (constant copy of the
    !!   interior edge cell).  Other physical BCs (REFLECTIVE, PERIODIC single-rank, wall)
    !!   are left unchanged here; they use the pre-filter ghost values for all passes, which
    !!   is a small approximation acceptable for a save-time output filter.
    !!
    !! @param q_cons_vf  Conserved-variable array whose ghost cells are refreshed in-place.
    !! @param mpi_dir    Spatial direction: 1 = x, 2 = y, 3 = z.
    impure subroutine s_lso_filter_ghost_refresh(q_cons_vf, mpi_dir)

        type(scalar_field), intent(inout) :: q_cons_vf(:)
        integer, intent(in)               :: mpi_dir
        integer                           :: i, j, k, l, beg_bc, end_bc

        select case (mpi_dir)
        case (1)
            beg_bc = bc_x%beg
            end_bc = bc_x%end
        case (2)
            beg_bc = bc_y%beg
            end_bc = bc_y%end
        case (3)
            beg_bc = bc_z%beg
            end_bc = bc_z%end
        end select

        ! MPI rank boundaries
#ifdef MFC_MPI
        if (beg_bc >= 0) call s_mpi_sendrecv_variables_buffers(q_cons_vf, mpi_dir, -1, sys_size)
        if (end_bc >= 0) call s_mpi_sendrecv_variables_buffers(q_cons_vf, mpi_dir, 1, sys_size)
#endif

        ! Physical boundaries: BC_GHOST_EXTRAP — constant copy of the interior edge cell.
        ! Re-applying between passes ensures the stencil of the next pass reads the
        ! filtered boundary value rather than the original pre-filter value, which would
        ! otherwise produce a step-like discontinuity at the domain edge and cause ringing.
        select case (mpi_dir)
        case (1)
            if (beg_bc == BC_GHOST_EXTRAP) then
                do i = 1, sys_size
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                q_cons_vf(i)%sf(-j, k, l) = q_cons_vf(i)%sf(0, k, l)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
            end if
            if (end_bc == BC_GHOST_EXTRAP) then
                do i = 1, sys_size
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                q_cons_vf(i)%sf(m + j, k, l) = q_cons_vf(i)%sf(m, k, l)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
            end if
        case (2)
            if (beg_bc == BC_GHOST_EXTRAP) then
                do i = 1, sys_size
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 1, buff_size
                            do j = 0, m
                                q_cons_vf(i)%sf(j, -k, l) = q_cons_vf(i)%sf(j, 0, l)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
            end if
            if (end_bc == BC_GHOST_EXTRAP) then
                do i = 1, sys_size
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 1, buff_size
                            do j = 0, m
                                q_cons_vf(i)%sf(j, n + k, l) = q_cons_vf(i)%sf(j, n, l)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
            end if
        case (3)
            if (beg_bc == BC_GHOST_EXTRAP) then
                do i = 1, sys_size
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do k = 0, n
                        do j = 0, m
                            do l = 1, buff_size
                                q_cons_vf(i)%sf(j, k, -l) = q_cons_vf(i)%sf(j, k, 0)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
            end if
            if (end_bc == BC_GHOST_EXTRAP) then
                do i = 1, sys_size
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do k = 0, n
                        do j = 0, m
                            do l = 1, buff_size
                                q_cons_vf(i)%sf(j, k, p + l) = q_cons_vf(i)%sf(j, k, p)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
            end if
        end select

    end subroutine s_lso_filter_ghost_refresh

end module m_lso_filter
