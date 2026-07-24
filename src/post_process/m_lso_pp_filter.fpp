!>
!! @file
!! @brief Contains module m_lso_pp_filter

#:include 'macros.fpp'

!> @brief Post_process-side additional Gaussian filter.
!!
!! Applies a second 9-point FIR pass with sigma_2 = sqrt(sigma_target^2 - sigma_in^2)
!! using the lso_pp_a_* weights from the Python BCD design. Stat product fields are
!! computed afterwards from the filtered state.
!!
!! IB-aware normalization: s_apply_lso_pp_filter_masked filters with a weight field w,
!! q <- filter(w*q)/filter(w). For original (unfiltered) input w is the binary gas mask
!! from ib_markers; for pre-filtered input w is the simulation-written stage-1 filtered
!! mask (lso_mask_<t>.dat), which composes the two-stage normalized convolution exactly:
!! filter2(w*qhat1)/filter2(w) = filter2(filter1(m*q))/filter2(filter1(m)).
module m_lso_pp_filter

    use m_derived_types
    use m_global_parameters
    use m_mpi_common
    use m_constants
    use m_variables_conversion, only: gammas
    use m_data_input, only: ib_markers

    implicit none

    private

    public :: s_initialize_lso_pp_filter_module, s_finalize_lso_pp_filter_module, s_apply_lso_pp_filter, &
        & s_compute_lso_pp_stat_fields, q_lso_pp_stat_vf, s_apply_lso_pp_filter_masked, s_lso_pp_mask_from_ib

    ! Scratch buffer for one directional pass.
    real(wp), allocatable :: lso_pp_tmp(:,:,:)

    ! Filtered stat fields (lso_stat_wrt = T).
    type(scalar_field), allocatable :: q_lso_pp_stat_vf(:)

contains

    impure subroutine s_initialize_lso_pp_filter_module()

        integer :: i

        allocate (lso_pp_tmp(0:m,0:n,0:p))

        if (lso_stat_wrt .and. n_lso_stat > 0) then
            allocate (q_lso_pp_stat_vf(1:n_lso_stat))
            do i = 1, n_lso_stat
                allocate (q_lso_pp_stat_vf(i)%sf(0:m,0:n,0:p))
            end do
        end if

    end subroutine s_initialize_lso_pp_filter_module

    impure subroutine s_finalize_lso_pp_filter_module()

        integer :: i

        if (allocated(lso_pp_tmp)) deallocate (lso_pp_tmp)

        if (allocated(q_lso_pp_stat_vf)) then
            do i = 1, n_lso_stat
                if (associated(q_lso_pp_stat_vf(i)%sf)) deallocate (q_lso_pp_stat_vf(i)%sf)
            end do
            deallocate (q_lso_pp_stat_vf)
        end if

    end subroutine s_finalize_lso_pp_filter_module

    !> Apply the post_process LSO filter to q_cons_vf in place.
    impure subroutine s_apply_lso_pp_filter(q_cons_vf)

        type(scalar_field), intent(inout) :: q_cons_vf(:)
        integer                           :: i, ipass, j, k, l, nv
        real(wp)                          :: c0, c1, c2, c3, c4

        nv = size(q_cons_vf)

        ! x-direction

        call s_lso_pp_filter_ghost_refresh(q_cons_vf, 1)
        do ipass = 1, lso_pp_n_passes_x
            c0 = lso_pp_a_x(1, ipass)
            c1 = lso_pp_a_x(2, ipass)
            c2 = lso_pp_a_x(3, ipass)
            c3 = lso_pp_a_x(4, ipass)
            c4 = lso_pp_a_x(5, ipass)
            do i = 1, nv
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            lso_pp_tmp(j, k, l) = c0*real(q_cons_vf(i)%sf(j, k, l), wp) + c1*(real(q_cons_vf(i)%sf(j - 1, k, l), &
                                       & wp) + real(q_cons_vf(i)%sf(j + 1, k, l), wp)) + c2*(real(q_cons_vf(i)%sf(j - 2, k, l), &
                                       & wp) + real(q_cons_vf(i)%sf(j + 2, k, l), wp)) + c3*(real(q_cons_vf(i)%sf(j - 3, k, l), &
                                       & wp) + real(q_cons_vf(i)%sf(j + 3, k, l), wp)) + c4*(real(q_cons_vf(i)%sf(j - 4, k, l), &
                                       & wp) + real(q_cons_vf(i)%sf(j + 4, k, l), wp))
                        end do
                    end do
                end do
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_cons_vf(i)%sf(j, k, l) = real(lso_pp_tmp(j, k, l), stp)
                        end do
                    end do
                end do
            end do
            if (ipass < lso_pp_n_passes_x) call s_lso_pp_filter_ghost_refresh(q_cons_vf, 1)
        end do

        ! y-direction (2D/3D)
        if (n > 0) then
            call s_lso_pp_filter_ghost_refresh(q_cons_vf, 2)
            do ipass = 1, lso_pp_n_passes_y
                c0 = lso_pp_a_y(1, ipass)
                c1 = lso_pp_a_y(2, ipass)
                c2 = lso_pp_a_y(3, ipass)
                c3 = lso_pp_a_y(4, ipass)
                c4 = lso_pp_a_y(5, ipass)
                do i = 1, nv
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                lso_pp_tmp(j, k, l) = c0*real(q_cons_vf(i)%sf(j, k, l), wp) + c1*(real(q_cons_vf(i)%sf(j, k - 1, &
                                           & l), wp) + real(q_cons_vf(i)%sf(j, k + 1, l), wp)) + c2*(real(q_cons_vf(i)%sf(j, &
                                           & k - 2, l), wp) + real(q_cons_vf(i)%sf(j, k + 2, l), &
                                           & wp)) + c3*(real(q_cons_vf(i)%sf(j, k - 3, l), wp) + real(q_cons_vf(i)%sf(j, k + 3, &
                                           & l), wp)) + c4*(real(q_cons_vf(i)%sf(j, k - 4, l), wp) + real(q_cons_vf(i)%sf(j, &
                                           & k + 4, l), wp))
                            end do
                        end do
                    end do
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                q_cons_vf(i)%sf(j, k, l) = real(lso_pp_tmp(j, k, l), stp)
                            end do
                        end do
                    end do
                end do
                if (ipass < lso_pp_n_passes_y) call s_lso_pp_filter_ghost_refresh(q_cons_vf, 2)
            end do
        end if

        ! z-direction (3D)
        if (p > 0) then
            call s_lso_pp_filter_ghost_refresh(q_cons_vf, 3)
            do ipass = 1, lso_pp_n_passes_z
                c0 = lso_pp_a_z(1, ipass)
                c1 = lso_pp_a_z(2, ipass)
                c2 = lso_pp_a_z(3, ipass)
                c3 = lso_pp_a_z(4, ipass)
                c4 = lso_pp_a_z(5, ipass)
                do i = 1, nv
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                lso_pp_tmp(j, k, l) = c0*real(q_cons_vf(i)%sf(j, k, l), wp) + c1*(real(q_cons_vf(i)%sf(j, k, &
                                           & l - 1), wp) + real(q_cons_vf(i)%sf(j, k, l + 1), wp)) + c2*(real(q_cons_vf(i)%sf(j, &
                                           & k, l - 2), wp) + real(q_cons_vf(i)%sf(j, k, l + 2), &
                                           & wp)) + c3*(real(q_cons_vf(i)%sf(j, k, l - 3), wp) + real(q_cons_vf(i)%sf(j, k, &
                                           & l + 3), wp)) + c4*(real(q_cons_vf(i)%sf(j, k, l - 4), wp) + real(q_cons_vf(i)%sf(j, &
                                           & k, l + 4), wp))
                            end do
                        end do
                    end do
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                q_cons_vf(i)%sf(j, k, l) = real(lso_pp_tmp(j, k, l), stp)
                            end do
                        end do
                    end do
                end do
                if (ipass < lso_pp_n_passes_z) call s_lso_pp_filter_ghost_refresh(q_cons_vf, 3)
            end do
        end if

    end subroutine s_apply_lso_pp_filter

    !> Fill the interior of w_vf with the binary gas mask from ib_markers: 1 in fluid, 0 inside an immersed body. Used when
    !! post_process filters ORIGINAL data.
    impure subroutine s_lso_pp_mask_from_ib(w_vf)

        type(scalar_field), intent(inout) :: w_vf(1:1)
        integer                           :: j, k, l

        do l = 0, p
            do k = 0, n
                do j = 0, m
                    if (ib_markers%sf(j, k, l) > 0) then
                        w_vf(1)%sf(j, k, l) = 0._stp
                    else
                        w_vf(1)%sf(j, k, l) = 1._stp
                    end if
                end do
            end do
        end do

    end subroutine s_lso_pp_mask_from_ib

    !> Mask-normalized filtering: q <- filter(w*q) / filter(w), applied in place. The solid interior takes the normalized fluid
    !! average (finite primitives everywhere); the floor only guards cells with almost no fluid within the filter reach.
    impure subroutine s_apply_lso_pp_filter_masked(q_cons_vf, w_vf)

        type(scalar_field), intent(inout) :: q_cons_vf(:), w_vf(1:1)
        integer                           :: i, j, k, l

        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_vf(i)%sf(j, k, l) = real(real(q_cons_vf(i)%sf(j, k, l), wp)*real(w_vf(1)%sf(j, k, l), wp), stp)
                    end do
                end do
            end do
        end do

        call s_apply_lso_pp_filter(q_cons_vf)  ! numerator   filter(w*q)
        call s_apply_lso_pp_filter(w_vf)  ! denominator filter(w)

        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_vf(i)%sf(j, k, l) = real(real(q_cons_vf(i)%sf(j, k, l), wp)/max(real(w_vf(1)%sf(j, k, l), wp), &
                                  & 1.0e-3_wp), stp)
                    end do
                end do
            end do
        end do

    end subroutine s_apply_lso_pp_filter_masked

    !> Refresh ghosts for direction mpi_dir (1=x, 2=y, 3=z). MPI ghosts come via s_mpi_sendrecv_variables_buffers; BC_GHOST_EXTRAP
    !! faces are re-extrapolated from the current edge cell.
    impure subroutine s_lso_pp_filter_ghost_refresh(q_cons_vf, mpi_dir)

        type(scalar_field), intent(inout) :: q_cons_vf(:)
        integer, intent(in)               :: mpi_dir
        integer                           :: i, j, k, l, beg_bc, end_bc, nv

        nv = size(q_cons_vf)

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

#ifdef MFC_MPI
        if (beg_bc >= 0) call s_mpi_sendrecv_variables_buffers(q_cons_vf, mpi_dir, -1, nv)
        if (end_bc >= 0) call s_mpi_sendrecv_variables_buffers(q_cons_vf, mpi_dir, 1, nv)
#endif

        select case (mpi_dir)
        case (1)
            if (beg_bc == BC_GHOST_EXTRAP) then
                do i = 1, nv
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                q_cons_vf(i)%sf(-j, k, l) = q_cons_vf(i)%sf(0, k, l)
                            end do
                        end do
                    end do
                end do
            end if
            if (end_bc == BC_GHOST_EXTRAP) then
                do i = 1, nv
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                q_cons_vf(i)%sf(m + j, k, l) = q_cons_vf(i)%sf(m, k, l)
                            end do
                        end do
                    end do
                end do
            end if
        case (2)
            if (beg_bc == BC_GHOST_EXTRAP) then
                do i = 1, nv
                    do l = 0, p
                        do k = 1, buff_size
                            do j = 0, m
                                q_cons_vf(i)%sf(j, -k, l) = q_cons_vf(i)%sf(j, 0, l)
                            end do
                        end do
                    end do
                end do
            end if
            if (end_bc == BC_GHOST_EXTRAP) then
                do i = 1, nv
                    do l = 0, p
                        do k = 1, buff_size
                            do j = 0, m
                                q_cons_vf(i)%sf(j, n + k, l) = q_cons_vf(i)%sf(j, n, l)
                            end do
                        end do
                    end do
                end do
            end if
        case (3)
            if (beg_bc == BC_GHOST_EXTRAP) then
                do i = 1, nv
                    do l = 1, buff_size
                        do k = 0, n
                            do j = 0, m
                                q_cons_vf(i)%sf(j, k, -l) = q_cons_vf(i)%sf(j, k, 0)
                            end do
                        end do
                    end do
                end do
            end if
            if (end_bc == BC_GHOST_EXTRAP) then
                do i = 1, nv
                    do l = 1, buff_size
                        do k = 0, n
                            do j = 0, m
                                q_cons_vf(i)%sf(j, k, p + l) = q_cons_vf(i)%sf(j, k, p)
                            end do
                        end do
                    end do
                end do
            end if
        end select

    end subroutine s_lso_pp_filter_ghost_refresh

    !> Build the LSO stat product fields from the post_process-filtered conserved state. No IB markers in post_process, so phi_p = 0
    !! and gas_mask = 1 everywhere. CPU loops.
    impure subroutine s_compute_lso_pp_stat_fields(q_cons_vf)

        type(scalar_field), intent(in) :: q_cons_vf(:)
        integer                        :: i, j, k, l
        real(wp)                       :: rho, rho_loc
        real(wp)                       :: mom1, mom2, mom3
        real(wp)                       :: u1, u2, u3
        real(wp)                       :: E_loc, ke, e_int, T_loc
        ! Gradient quantities (viscous pass)
        real(wp) :: rho_jm, rho_jp, rho_km, rho_kp, rho_lm, rho_lp
        real(wp) :: u1_jm, u1_jp, u1_km, u1_kp, u1_lm, u1_lp
        real(wp) :: u2_jm, u2_jp, u2_km, u2_kp, u2_lm, u2_lp
        real(wp) :: u3_jm, u3_jp, u3_km, u3_kp, u3_lm, u3_lp
        real(wp) :: T_jm, T_jp, T_km, T_kp, T_lm, T_lp
        real(wp) :: ddx, ddy, ddz
        real(wp) :: du1dx, du1dy, du1dz
        real(wp) :: du2dx, du2dy, du2dz
        real(wp) :: du3dx, du3dy, du3dz
        real(wp) :: dTdx, dTdy, dTdz
        real(wp) :: div_u
        real(wp) :: tau11, tau12, tau13, tau22, tau23, tau33
        real(wp) :: q1, q2, q3

        ! Pass 1: algebraic products. gas_mask = 1 everywhere (no IB).

        do l = 0, p
            do k = 0, n
                do j = 0, m
                    rho = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k, l), wp), sgm_eps)
                    mom1 = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, k, l), wp)
                    if (n > 0) then
                        mom2 = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k, l), wp)
                    else
                        mom2 = 0._wp
                    end if
                    if (p > 0) then
                        mom3 = real(q_cons_vf(eqn_idx%mom%beg + 2)%sf(j, k, l), wp)
                    else
                        mom3 = 0._wp
                    end if
                    E_loc = real(q_cons_vf(eqn_idx%E)%sf(j, k, l), wp)

                    u1 = mom1/rho
                    u2 = mom2/rho
                    u3 = mom3/rho
                    ke = 0.5_wp*(mom1**2 + mom2**2 + mom3**2)/rho
                    e_int = (E_loc - ke)/rho
                    T_loc = e_int/(gammas(1)*lso_R_gas)

                    ! phi_p, rho scalar, rhoke scalar, u_p (no IB in post_process filter)
                    q_lso_pp_stat_vf(lso_stat_phi_p_beg)%sf(j, k, l) = 0._stp
                    q_lso_pp_stat_vf(lso_stat_rho_beg)%sf(j, k, l) = real(rho, stp)
                    q_lso_pp_stat_vf(lso_stat_rhoke_beg)%sf(j, k, l) = real((mom1**2 + mom2**2 + mom3**2)/rho, stp)
                    q_lso_pp_stat_vf(lso_stat_up_beg)%sf(j, k, l) = 0._stp
                    if (n > 0) q_lso_pp_stat_vf(lso_stat_up_beg + 1)%sf(j, k, l) = 0._stp
                    if (p > 0) q_lso_pp_stat_vf(lso_stat_up_beg + 2)%sf(j, k, l) = 0._stp

                    ! rho*u
                    q_lso_pp_stat_vf(lso_stat_rhou_beg)%sf(j, k, l) = real(mom1, stp)
                    if (n > 0) q_lso_pp_stat_vf(lso_stat_rhou_beg + 1)%sf(j, k, l) = real(mom2, stp)
                    if (p > 0) q_lso_pp_stat_vf(lso_stat_rhou_beg + 2)%sf(j, k, l) = real(mom3, stp)

                    ! rho*u*u upper triangle
                    q_lso_pp_stat_vf(lso_stat_rhouu_beg)%sf(j, k, l) = real(mom1*u1, stp)
                    if (n > 0) then
                        q_lso_pp_stat_vf(lso_stat_rhouu_beg + 1)%sf(j, k, l) = real(mom1*u2, stp)
                        if (p > 0) then
                            q_lso_pp_stat_vf(lso_stat_rhouu_beg + 2)%sf(j, k, l) = real(mom1*u3, stp)
                            q_lso_pp_stat_vf(lso_stat_rhouu_beg + 3)%sf(j, k, l) = real(mom2*u2, stp)
                            q_lso_pp_stat_vf(lso_stat_rhouu_beg + 4)%sf(j, k, l) = real(mom2*u3, stp)
                            q_lso_pp_stat_vf(lso_stat_rhouu_beg + 5)%sf(j, k, l) = real(mom3*u3, stp)
                        else
                            q_lso_pp_stat_vf(lso_stat_rhouu_beg + 2)%sf(j, k, l) = real(mom2*u2, stp)
                        end if
                    end if

                    ! rho*u*|u|^2
                    q_lso_pp_stat_vf(lso_stat_rhouke_beg)%sf(j, k, l) = real(mom1*(u1**2 + u2**2 + u3**2), stp)
                    if (n > 0) q_lso_pp_stat_vf(lso_stat_rhouke_beg + 1)%sf(j, k, l) = real(mom2*(u1**2 + u2**2 + u3**2), stp)
                    if (p > 0) q_lso_pp_stat_vf(lso_stat_rhouke_beg + 2)%sf(j, k, l) = real(mom3*(u1**2 + u2**2 + u3**2), stp)

                    ! rho*u*T
                    q_lso_pp_stat_vf(lso_stat_rhouT_beg)%sf(j, k, l) = real(mom1*T_loc, stp)
                    if (n > 0) q_lso_pp_stat_vf(lso_stat_rhouT_beg + 1)%sf(j, k, l) = real(mom2*T_loc, stp)
                    if (p > 0) q_lso_pp_stat_vf(lso_stat_rhouT_beg + 2)%sf(j, k, l) = real(mom3*T_loc, stp)

                    ! Zero the viscous-pass slots.
                    do i = lso_stat_tau_beg, lso_stat_tau_end
                        q_lso_pp_stat_vf(i)%sf(j, k, l) = 0._stp
                    end do
                    do i = lso_stat_q_beg, lso_stat_q_end
                        q_lso_pp_stat_vf(i)%sf(j, k, l) = 0._stp
                    end do
                    do i = lso_stat_rhotau_u_beg, lso_stat_rhotau_u_end
                        q_lso_pp_stat_vf(i)%sf(j, k, l) = 0._stp
                    end do
                end do
            end do
        end do

        ! Pass 2: viscous stress, heat flux, viscous power flux from centred diffs.
        if (lso_mu > 0._wp) then
            if (n == 0) then
                ! 1D: x-gradients only.
                do j = 0, m
                    rho_loc = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, 0, 0), wp), sgm_eps)

                    rho_jm = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j - 1, 0, 0), wp), sgm_eps)
                    u1_jm = real(q_cons_vf(eqn_idx%mom%beg)%sf(j - 1, 0, 0), wp)/rho_jm
                    T_jm = (real(q_cons_vf(eqn_idx%E)%sf(j - 1, 0, 0), wp) - 0.5_wp*u1_jm**2*rho_jm)/(rho_jm*gammas(1)*lso_R_gas)

                    rho_jp = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j + 1, 0, 0), wp), sgm_eps)
                    u1_jp = real(q_cons_vf(eqn_idx%mom%beg)%sf(j + 1, 0, 0), wp)/rho_jp
                    T_jp = (real(q_cons_vf(eqn_idx%E)%sf(j + 1, 0, 0), wp) - 0.5_wp*u1_jp**2*rho_jp)/(rho_jp*gammas(1)*lso_R_gas)

                    ddx = x_cc(j + 1) - x_cc(j - 1)
                    du1dx = (u1_jp - u1_jm)/ddx
                    dTdx = (T_jp - T_jm)/ddx

                    tau11 = lso_mu*(2._wp*du1dx)
                    q1 = -lso_conductivity*dTdx

                    u1 = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, 0, 0), wp)/rho_loc

                    q_lso_pp_stat_vf(lso_stat_tau_beg)%sf(j, 0, 0) = real(tau11, stp)
                    q_lso_pp_stat_vf(lso_stat_q_beg)%sf(j, 0, 0) = real(q1, stp)
                    q_lso_pp_stat_vf(lso_stat_rhotau_u_beg)%sf(j, 0, 0) = real(rho_loc*tau11*u1, stp)
                end do
            else if (p == 0) then
                ! 2D: x- and y-gradients.
                do k = 0, n
                    do j = 0, m
                        rho_loc = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k, 0), wp), sgm_eps)

                        rho_jm = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j - 1, k, 0), wp), sgm_eps)
                        u1_jm = real(q_cons_vf(eqn_idx%mom%beg)%sf(j - 1, k, 0), wp)/rho_jm
                        u2_jm = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j - 1, k, 0), wp)/rho_jm
                        T_jm = (real(q_cons_vf(eqn_idx%E)%sf(j - 1, k, 0), &
                                & wp) - 0.5_wp*(u1_jm**2 + u2_jm**2)*rho_jm)/(rho_jm*gammas(1)*lso_R_gas)

                        rho_jp = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j + 1, k, 0), wp), sgm_eps)
                        u1_jp = real(q_cons_vf(eqn_idx%mom%beg)%sf(j + 1, k, 0), wp)/rho_jp
                        u2_jp = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j + 1, k, 0), wp)/rho_jp
                        T_jp = (real(q_cons_vf(eqn_idx%E)%sf(j + 1, k, 0), &
                                & wp) - 0.5_wp*(u1_jp**2 + u2_jp**2)*rho_jp)/(rho_jp*gammas(1)*lso_R_gas)

                        rho_km = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k - 1, 0), wp), sgm_eps)
                        u1_km = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, k - 1, 0), wp)/rho_km
                        u2_km = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k - 1, 0), wp)/rho_km
                        T_km = (real(q_cons_vf(eqn_idx%E)%sf(j, k - 1, 0), &
                                & wp) - 0.5_wp*(u1_km**2 + u2_km**2)*rho_km)/(rho_km*gammas(1)*lso_R_gas)

                        rho_kp = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k + 1, 0), wp), sgm_eps)
                        u1_kp = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, k + 1, 0), wp)/rho_kp
                        u2_kp = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k + 1, 0), wp)/rho_kp
                        T_kp = (real(q_cons_vf(eqn_idx%E)%sf(j, k + 1, 0), &
                                & wp) - 0.5_wp*(u1_kp**2 + u2_kp**2)*rho_kp)/(rho_kp*gammas(1)*lso_R_gas)

                        ddx = x_cc(j + 1) - x_cc(j - 1)
                        ddy = y_cc(k + 1) - y_cc(k - 1)

                        du1dx = (u1_jp - u1_jm)/ddx
                        du1dy = (u1_kp - u1_km)/ddy
                        du2dx = (u2_jp - u2_jm)/ddx
                        du2dy = (u2_kp - u2_km)/ddy
                        dTdx = (T_jp - T_jm)/ddx
                        dTdy = (T_kp - T_km)/ddy

                        div_u = du1dx + du2dy
                        tau11 = lso_mu*(2._wp*du1dx - (2._wp/3._wp)*div_u)
                        tau12 = lso_mu*(du1dy + du2dx)
                        tau22 = lso_mu*(2._wp*du2dy - (2._wp/3._wp)*div_u)
                        q1 = -lso_conductivity*dTdx
                        q2 = -lso_conductivity*dTdy

                        u1 = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, k, 0), wp)/rho_loc
                        u2 = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k, 0), wp)/rho_loc

                        q_lso_pp_stat_vf(lso_stat_tau_beg)%sf(j, k, 0) = real(tau11, stp)
                        q_lso_pp_stat_vf(lso_stat_tau_beg + 1)%sf(j, k, 0) = real(tau12, stp)
                        q_lso_pp_stat_vf(lso_stat_tau_beg + 2)%sf(j, k, 0) = real(tau22, stp)
                        q_lso_pp_stat_vf(lso_stat_q_beg)%sf(j, k, 0) = real(q1, stp)
                        q_lso_pp_stat_vf(lso_stat_q_beg + 1)%sf(j, k, 0) = real(q2, stp)
                        q_lso_pp_stat_vf(lso_stat_rhotau_u_beg)%sf(j, k, 0) = real(rho_loc*(tau11*u1 + tau12*u2), stp)
                        q_lso_pp_stat_vf(lso_stat_rhotau_u_beg + 1)%sf(j, k, 0) = real(rho_loc*(tau12*u1 + tau22*u2), stp)
                    end do
                end do
            else
                ! 3D: x-, y-, z-gradients.
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            rho_loc = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k, l), wp), sgm_eps)

                            rho_jm = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j - 1, k, l), wp), sgm_eps)
                            u1_jm = real(q_cons_vf(eqn_idx%mom%beg)%sf(j - 1, k, l), wp)/rho_jm
                            u2_jm = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j - 1, k, l), wp)/rho_jm
                            u3_jm = real(q_cons_vf(eqn_idx%mom%beg + 2)%sf(j - 1, k, l), wp)/rho_jm
                            T_jm = (real(q_cons_vf(eqn_idx%E)%sf(j - 1, k, l), &
                                    & wp) - 0.5_wp*(u1_jm**2 + u2_jm**2 + u3_jm**2)*rho_jm)/(rho_jm*gammas(1)*lso_R_gas)

                            rho_jp = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j + 1, k, l), wp), sgm_eps)
                            u1_jp = real(q_cons_vf(eqn_idx%mom%beg)%sf(j + 1, k, l), wp)/rho_jp
                            u2_jp = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j + 1, k, l), wp)/rho_jp
                            u3_jp = real(q_cons_vf(eqn_idx%mom%beg + 2)%sf(j + 1, k, l), wp)/rho_jp
                            T_jp = (real(q_cons_vf(eqn_idx%E)%sf(j + 1, k, l), &
                                    & wp) - 0.5_wp*(u1_jp**2 + u2_jp**2 + u3_jp**2)*rho_jp)/(rho_jp*gammas(1)*lso_R_gas)

                            rho_km = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k - 1, l), wp), sgm_eps)
                            u1_km = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, k - 1, l), wp)/rho_km
                            u2_km = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k - 1, l), wp)/rho_km
                            u3_km = real(q_cons_vf(eqn_idx%mom%beg + 2)%sf(j, k - 1, l), wp)/rho_km
                            T_km = (real(q_cons_vf(eqn_idx%E)%sf(j, k - 1, l), &
                                    & wp) - 0.5_wp*(u1_km**2 + u2_km**2 + u3_km**2)*rho_km)/(rho_km*gammas(1)*lso_R_gas)

                            rho_kp = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k + 1, l), wp), sgm_eps)
                            u1_kp = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, k + 1, l), wp)/rho_kp
                            u2_kp = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k + 1, l), wp)/rho_kp
                            u3_kp = real(q_cons_vf(eqn_idx%mom%beg + 2)%sf(j, k + 1, l), wp)/rho_kp
                            T_kp = (real(q_cons_vf(eqn_idx%E)%sf(j, k + 1, l), &
                                    & wp) - 0.5_wp*(u1_kp**2 + u2_kp**2 + u3_kp**2)*rho_kp)/(rho_kp*gammas(1)*lso_R_gas)

                            rho_lm = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k, l - 1), wp), sgm_eps)
                            u1_lm = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, k, l - 1), wp)/rho_lm
                            u2_lm = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k, l - 1), wp)/rho_lm
                            u3_lm = real(q_cons_vf(eqn_idx%mom%beg + 2)%sf(j, k, l - 1), wp)/rho_lm
                            T_lm = (real(q_cons_vf(eqn_idx%E)%sf(j, k, l - 1), &
                                    & wp) - 0.5_wp*(u1_lm**2 + u2_lm**2 + u3_lm**2)*rho_lm)/(rho_lm*gammas(1)*lso_R_gas)

                            rho_lp = max(real(q_cons_vf(eqn_idx%cont%beg)%sf(j, k, l + 1), wp), sgm_eps)
                            u1_lp = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, k, l + 1), wp)/rho_lp
                            u2_lp = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k, l + 1), wp)/rho_lp
                            u3_lp = real(q_cons_vf(eqn_idx%mom%beg + 2)%sf(j, k, l + 1), wp)/rho_lp
                            T_lp = (real(q_cons_vf(eqn_idx%E)%sf(j, k, l + 1), &
                                    & wp) - 0.5_wp*(u1_lp**2 + u2_lp**2 + u3_lp**2)*rho_lp)/(rho_lp*gammas(1)*lso_R_gas)

                            ddx = x_cc(j + 1) - x_cc(j - 1)
                            ddy = y_cc(k + 1) - y_cc(k - 1)
                            ddz = z_cc(l + 1) - z_cc(l - 1)

                            du1dx = (u1_jp - u1_jm)/ddx
                            du1dy = (u1_kp - u1_km)/ddy
                            du1dz = (u1_lp - u1_lm)/ddz
                            du2dx = (u2_jp - u2_jm)/ddx
                            du2dy = (u2_kp - u2_km)/ddy
                            du2dz = (u2_lp - u2_lm)/ddz
                            du3dx = (u3_jp - u3_jm)/ddx
                            du3dy = (u3_kp - u3_km)/ddy
                            du3dz = (u3_lp - u3_lm)/ddz
                            dTdx = (T_jp - T_jm)/ddx
                            dTdy = (T_kp - T_km)/ddy
                            dTdz = (T_lp - T_lm)/ddz

                            div_u = du1dx + du2dy + du3dz
                            tau11 = lso_mu*(2._wp*du1dx - (2._wp/3._wp)*div_u)
                            tau12 = lso_mu*(du1dy + du2dx)
                            tau13 = lso_mu*(du1dz + du3dx)
                            tau22 = lso_mu*(2._wp*du2dy - (2._wp/3._wp)*div_u)
                            tau23 = lso_mu*(du2dz + du3dy)
                            tau33 = lso_mu*(2._wp*du3dz - (2._wp/3._wp)*div_u)
                            q1 = -lso_conductivity*dTdx
                            q2 = -lso_conductivity*dTdy
                            q3 = -lso_conductivity*dTdz

                            u1 = real(q_cons_vf(eqn_idx%mom%beg)%sf(j, k, l), wp)/rho_loc
                            u2 = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k, l), wp)/rho_loc
                            u3 = real(q_cons_vf(eqn_idx%mom%beg + 2)%sf(j, k, l), wp)/rho_loc

                            q_lso_pp_stat_vf(lso_stat_tau_beg)%sf(j, k, l) = real(tau11, stp)
                            q_lso_pp_stat_vf(lso_stat_tau_beg + 1)%sf(j, k, l) = real(tau12, stp)
                            q_lso_pp_stat_vf(lso_stat_tau_beg + 2)%sf(j, k, l) = real(tau13, stp)
                            q_lso_pp_stat_vf(lso_stat_tau_beg + 3)%sf(j, k, l) = real(tau22, stp)
                            q_lso_pp_stat_vf(lso_stat_tau_beg + 4)%sf(j, k, l) = real(tau23, stp)
                            q_lso_pp_stat_vf(lso_stat_tau_beg + 5)%sf(j, k, l) = real(tau33, stp)
                            q_lso_pp_stat_vf(lso_stat_q_beg)%sf(j, k, l) = real(q1, stp)
                            q_lso_pp_stat_vf(lso_stat_q_beg + 1)%sf(j, k, l) = real(q2, stp)
                            q_lso_pp_stat_vf(lso_stat_q_beg + 2)%sf(j, k, l) = real(q3, stp)
                            q_lso_pp_stat_vf(lso_stat_rhotau_u_beg)%sf(j, k, l) = real(rho_loc*(tau11*u1 + tau12*u2 + tau13*u3), &
                                             & stp)
                            q_lso_pp_stat_vf(lso_stat_rhotau_u_beg + 1)%sf(j, k, &
                                             & l) = real(rho_loc*(tau12*u1 + tau22*u2 + tau23*u3), stp)
                            q_lso_pp_stat_vf(lso_stat_rhotau_u_beg + 2)%sf(j, k, &
                                             & l) = real(rho_loc*(tau13*u1 + tau23*u2 + tau33*u3), stp)
                        end do
                    end do
                end do
            end if
        end if

    end subroutine s_compute_lso_pp_stat_fields

end module m_lso_pp_filter
