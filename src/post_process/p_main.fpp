!>
!! @file
!! @brief Contains program p_main

!> Post-process raw simulation data into formatted database files (Silo-HDF5 or Binary)
program p_main

    use m_global_parameters
    use m_start_up
    use m_data_output, only: s_switch_output_dirs
    use m_data_input, only: q_cons_vf
    use m_lso_pp_filter, only: s_apply_lso_pp_filter, s_compute_lso_pp_stat_fields

    implicit none

    integer :: t_step, nt_step  !< Iterator for the main time-stepping loop
    !> Generic storage for the name(s) of the flow variable(s) that will be added to the formatted database file(s)
    character(LEN=name_len) :: varname
    real(wp)                :: pres
    real(wp)                :: c
    real(wp)                :: H
    real(wp)                :: start, finish

    call s_initialize_mpi_domain()

    call s_initialize_modules()

    if (cfl_dt) then
        t_step = n_start
        n_save = int(t_stop/t_save) + 1
    else
        ! Setting the time-step iterator to the first time step to be post-processed
        t_step = t_step_start
    end if

    ! Time-Marching Loop
    do
        ! If all time-steps are not ready to be post-processed and one rank is faster than another, the slower rank processing the
        ! last available step might be killed when the faster rank attempts to process the first missing step, before the slower
        ! rank finishes writing the last available step. To avoid this, we force synchronization here.
        call s_mpi_barrier()

        call cpu_time(start)

        ! Primary pass: read (LSO-filtered when lso_filter_wrt=T) and write.
        call s_perform_time_step(t_step)

        ! Post-process FIR filter: apply an additional Gaussian filter to the conserved variables in place,
        ! then reconvert to primitive so s_save_data emits the filtered cons AND prim fields (s_apply_lso_pp_filter
        ! only touches q_cons_vf). When lso_stat_wrt=T, compute the stat product fields now, while q_cons_vf still
        ! holds the filtered state (the secondary pass below reloads unfiltered data). Requires lso_filter_wrt=T so
        ! s_save_data routes output to silo_hdf5_lso/.
        if (lso_pp_filter) then
            call s_apply_lso_pp_filter(q_cons_vf)
            call s_reconvert_filtered_to_primitive()
            if (lso_stat_wrt) call s_compute_lso_pp_stat_fields(q_cons_vf)
        end if

        call s_save_data(t_step, varname, pres, c, H)

        ! Secondary pass: when lso_filter_wrt=T, also write the unfiltered data to silo_hdf5/. Grid files are NOT re-opened
        ! (grid_loaded flag skips them), avoiding a known MPI-IO hang from re-opening the same file with
        ! MPI_FILE_OPEN(MPI_COMM_WORLD, fp) a second time. NOTE: When lso_down_sample_factor > 1 the LSO files are on a coarser grid
        ! than the unfiltered data. The secondary pass cannot switch grid dimensions at runtime, so it is skipped. Run post_process
        ! again with lso_filter_wrt=.false. to process the full-resolution unfiltered data.
        if (lso_filter_wrt .and. lso_down_sample_factor <= 1) then
            lso_filter_wrt = .false.
            call s_reload_data(t_step)
            call s_switch_output_dirs(.false.)
            call s_save_data(t_step, varname, pres, c, H)
            call s_switch_output_dirs(.true.)
            lso_filter_wrt = .true.
        end if

        ! Third pass: write LSO statistical product fields to silo_hdf5_lso_stat/. In post_process filter mode the
        ! fields were computed in-process above (from the filtered state); otherwise they are read from the
        ! simulation-written lso_stat binary. Done here (after the secondary pass) so the lso_stat directory switch
        ! does not disturb the primary/secondary Silo writes.
        if (lso_stat_wrt) then
            if (lso_pp_filter) then
                call s_save_lso_pp_stat_data(t_step)
            else if (lso_filter_wrt) then
                call s_save_lso_stat_data(t_step)
            end if
        end if

        call cpu_time(finish)

        wall_time = abs(finish - start)

        if (cfl_dt) then
            nt_step = t_step - n_start + 1
        else
            nt_step = (t_step - t_step_start)/t_step_save + 1
        end if
        if (nt_step >= 2) then
            wall_time_avg = (wall_time + (nt_step - 2)*wall_time_avg)/(nt_step - 1)
        else
            wall_time_avg = 0._wp
        end if

        if (cfl_dt) then
            if (t_step == n_save - 1) then
                exit
            end if
        else
            ! Adjust time-step iterator to reach final step if needed, else exit
            if ((t_step_stop - t_step) < t_step_save .and. t_step_stop /= t_step) then
                t_step = t_step_stop - t_step_save
            else if (t_step == t_step_stop) then
                exit
            end if
        end if

        if (cfl_dt) then
            t_step = t_step + 1
        else
            ! Incrementing time-step iterator to next time-step to be post-processed
            t_step = t_step + t_step_save
        end if
    end do
    ! END: Time-Marching Loop

    close (11)

    call s_finalize_modules()
end program p_main
