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

    implicit none

    private

    public :: s_initialize_lso_filter_module, s_apply_lso_filter, s_copy_and_apply_lso_filter, s_finalize_lso_filter_module, &
        & q_filt_vf

    !> Temporary working buffer for one directional filter pass (interior cells only)
    real(wp), allocatable, dimension(:,:,:) :: lso_tmp
    $:GPU_DECLARE(create='[lso_tmp]')

    !> Filtered copy of the conserved variables (only allocated when lso_filter_wrt = T)
    type(scalar_field), allocatable :: q_filt_vf(:)

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
        end if

    end subroutine s_initialize_lso_filter_module

    !> Deallocate the temporary pass buffer and, when allocated, q_filt_vf.
    impure subroutine s_finalize_lso_filter_module()

        integer :: i

        @:DEALLOCATE(lso_tmp)

        if (lso_filter_wrt) then
            do i = 1, sys_size
                @:DEALLOCATE(q_filt_vf(i)%sf)
            end do
            @:DEALLOCATE(q_filt_vf)
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
    !! filter is applied in-place on the interior (0:m, 0:n, 0:p); ghost cells are not updated between passes - this is acceptable
    !! for a save-time output filter.
    !!
    !! @note  When time_stepper == 1 (forward Euler), stor == 1 and the live conserved
    !! variable array is modified in-place. For all Runge-Kutta time steppers (stor == 2) only the save-copy is affected and the
    !! simulation state is unchanged. IBM cases virtually always use an RK scheme, so this is safe in practice.
    impure subroutine s_apply_lso_filter(q_cons_vf)

        type(scalar_field), intent(inout) :: q_cons_vf(:)
        integer                           :: i, ipass, j, k, l
        real(wp)                          :: c0, c1, c2, c3, c4

        do i = 1, sys_size
            ! x-direction passes
            do ipass = 1, lso_n_passes_x
                c0 = lso_a_x(1, ipass)
                c1 = lso_a_x(2, ipass)
                c2 = lso_a_x(3, ipass)
                c3 = lso_a_x(4, ipass)
                c4 = lso_a_x(5, ipass)
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

            ! y-direction passes (2D/3D only: n > 0)
            if (n > 0) then
                do ipass = 1, lso_n_passes_y
                    c0 = lso_a_y(1, ipass)
                    c1 = lso_a_y(2, ipass)
                    c2 = lso_a_y(3, ipass)
                    c3 = lso_a_y(4, ipass)
                    c4 = lso_a_y(5, ipass)
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
            end if

            ! z-direction passes (3D only: p > 0)
            if (p > 0) then
                do ipass = 1, lso_n_passes_z
                    c0 = lso_a_z(1, ipass)
                    c1 = lso_a_z(2, ipass)
                    c2 = lso_a_z(3, ipass)
                    c3 = lso_a_z(4, ipass)
                    c4 = lso_a_z(5, ipass)
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
            end if
        end do

    end subroutine s_apply_lso_filter

end module m_lso_filter
