! P3 AMR leaf-cell state and CUDA conservation diagnostic for S_N.
! It is enabled only by SNRT_P1_DIAGNOSTIC=1 at runtime.
module snrt_ramses_driver
  use amr_parameters, only: MAXLEVEL
  implicit none

  private
  public :: snrt_ramses_diagnose_level, snrt_ramses_advance_level

contains

  subroutine snrt_ramses_diagnose_level(ilevel)
    use amr_commons, only: levelmin, nstep_coarse, myid
    use snrt_state, only: snrt_state_sync_level, snrt_ndirection
    use snrt_cuda_interface, only: snrt_cuda_available, &
         snrt_cuda_angular_reduce_tf32
    use snrt_cuda_ledger_interface, only: snrt_cuda_weighted_sum_fp32
    use iso_c_binding, only: c_float
    implicit none

    integer, intent(in) :: ilevel
    integer, parameter :: nbin = 16
    character(len=16) :: env_value
    integer :: env_length, env_status
    integer :: i, idir, ibin, nrow, nnew, ierr
    integer, save :: last_checked(MAXLEVEL) = -1
    logical, save :: enabled_resolved = .false.
    logical, save :: enabled = .false.
    real(c_float) :: max_tensor_error, max_ledger_error, reference
    real(c_float), allocatable :: directional(:), weights(:), projection(:), &
         binned(:), scalar(:)

    if (.not. enabled_resolved) then
       env_value = ''
       call get_environment_variable('SNRT_P1_DIAGNOSTIC', env_value, &
            length=env_length, status=env_status)
       enabled = env_status == 0 .and. env_length == 1 .and. &
            env_value(1:1) == '1'
       enabled_resolved = .true.
       if (enabled .and. myid == 1) then
          write(*,'(A)') ' SNRT P3 diagnostic enabled: persistent AMR leaf state'
       endif
    endif
    if (.not. enabled) return
    if (ilevel < levelmin .or. ilevel > MAXLEVEL) return
    if (last_checked(ilevel) == nstep_coarse) return

    call snrt_state_sync_level(ilevel, nrow, nnew)
    last_checked(ilevel) = nstep_coarse
    if (nrow == 0) return
    if (.not. snrt_cuda_available()) then
       if (myid == 1) write(*,'(A)') ' SNRT P3 disabled: no CUDA device is visible'
       enabled = .false.
       return
    endif

    allocate(directional(nrow * snrt_ndirection), weights(snrt_ndirection), &
         projection(snrt_ndirection * nbin), binned(nrow * nbin), scalar(nrow))
    directional = 1.0_c_float
    weights = 1.0_c_float / real(snrt_ndirection, c_float)
    projection = 0.0_c_float
    do idir = 1, snrt_ndirection
       ibin = 1 + mod(idir - 1, nbin)
       projection((idir - 1) * nbin + ibin) = weights(idir)
    enddo

    call snrt_cuda_angular_reduce_tf32(directional, projection, binned, nrow, &
         snrt_ndirection, nbin, ierr)
    if (ierr /= 0) then
       if (myid == 1) write(*,'(A,I0)') ' SNRT P3 TF32 reduction failed, code=', ierr
       deallocate(directional, weights, projection, binned, scalar)
       enabled = .false.
       return
    endif
    call snrt_cuda_weighted_sum_fp32(directional, weights, scalar, nrow, &
         snrt_ndirection, ierr)
    if (ierr /= 0) then
       if (myid == 1) write(*,'(A,I0)') ' SNRT P3 ledger reduction failed, code=', ierr
       deallocate(directional, weights, projection, binned, scalar)
       enabled = .false.
       return
    endif

    reference = sum(weights)
    max_tensor_error = 0.0_c_float
    max_ledger_error = 0.0_c_float
    do i = 1, nrow
       max_tensor_error = max(max_tensor_error, &
            abs(sum(binned((i - 1) * nbin + 1:i * nbin)) - reference))
       max_ledger_error = max(max_ledger_error, abs(scalar(i) - reference))
    enddo
    if (myid == 1) write(*,'(A,I0,A,I0,A,I0,A,ES12.4,A,ES12.4)') &
         ' SNRT P3 leaf diagnostic level=', ilevel, ' rows=', nrow, &
         ' new_slots=', nnew, ' tensor_abs=', max_tensor_error, &
         ' ledger_abs=', max_ledger_error

    deallocate(directional, weights, projection, binned, scalar)
  end subroutine snrt_ramses_diagnose_level

  subroutine snrt_ramses_advance_level(ilevel)
    use amr_commons, only: levelmin, nstep_coarse, myid, dtnew, boxlen, &
         icoarse_min, icoarse_max
    use hydro_commons, only: uold
    use pm_commons, only: nsink, dMsmbh, xsink, eps_sink
    use snrt_state, only: snrt_ndirection, snrt_ngroups, snrt_intensity, &
         snrt_neutral_fraction, snrt_group_mean_energy_ev, &
         snrt_group_energy_fraction, snrt_group_cross_section_cm2, &
         snrt_state_get_slot
    use snrt_amr_topology, only: snrt_amr_build_same_level_neighbors
    use snrt_transport_step, only: snrt_transport_absorb_multigroup_prepared
    use snrt_angular_quadrature, only: snrt_angular_init
    use snrt_agn_locator, only: snrt_agn_find_local_leaf
    use snrt_agn_source, only: snrt_c_cgs, snrt_agn_photon_budget, &
         snrt_agn_deposit_isotropic
    use snrt_nlte_coupling, only: snrt_nlte_optical_depth_groups, &
         snrt_nlte_photo_source
    use snrt_cuda_interface, only: snrt_cuda_available
    use amr_parameters, only: dp, ndim
    use iso_c_binding, only: c_float
    use omp_lib, only: omp_get_wtime
    implicit none

    integer, intent(in) :: ilevel
    character(len=32) :: env_value
    integer :: env_length, env_status, read_status
    integer :: i, isink, igroup, ierr, nleaf, n_interface_face
    integer :: icell, islot, ilevel_found
    integer :: energy_index
    integer :: n_locator_calls, n_active_sources
    integer, allocatable :: leaf_cell(:), leaf_slot(:), neighbor(:,:)
    real(dp) :: direction_dp(snrt_ndirection,3), angular_weight(snrt_ndirection)
    real(dp) :: scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2
    real(dp) :: dt_s, dx_code, cell_volume_code, cdt_over_dx, scale_m
    real(dp) :: rho_code, neutral_fraction, ionized_fraction
    real(dp) :: tau_dp(snrt_ngroups), delta_accreted, epsilon_r
    real(dp) :: luminosity, emitted_photons, deposited_density
    real(dp) :: ionization_increment, heating_rate, heating_total
    real(dp) :: wall_start
    real(dp) :: wall_sub
    real(dp) :: t_setup, t_topology, t_nlte, t_source, t_transport, t_coupling
    real(dp) :: t_source_overhead, t_locator, t_budget, t_deposit
    real(c_float), allocatable :: optical_depth(:,:), neutral_hydrogen(:)
    real(c_float), allocatable :: absorbed_group(:,:)
    real(dp), allocatable :: accounted_mass_new(:)
    logical :: enabled, source_ok
    logical, save :: enabled_resolved = .false.
    real(dp), save :: reduced_c = 0.01d0
    integer, save :: level_filter = -1
    integer, save :: accounted_step = -huge(1)
    real(dp), allocatable, save :: accounted_mass(:)

    enabled = .false.
    if (.not. enabled_resolved) then
       env_value = ''
       call get_environment_variable('SNRT_RT_ENABLE', env_value, &
            length=env_length, status=env_status)
       enabled = env_status == 0 .and. env_length == 1 .and. &
            env_value(1:1) == '1'
       env_value = ''
       call get_environment_variable('SNRT_REDUCED_C', env_value, &
            length=env_length, status=env_status)
       if (env_status == 0 .and. env_length > 0) then
          read(env_value(1:env_length),*,iostat=read_status) reduced_c
          if (read_status /= 0 .or. reduced_c <= 0.0d0 .or. reduced_c > 1.0d0) &
               reduced_c = 0.01d0
       end if
       env_value = ''
       call get_environment_variable('SNRT_RT_LEVEL', env_value, &
            length=env_length, status=env_status)
       if (env_status == 0 .and. env_length > 0) then
          read(env_value(1:env_length),*,iostat=read_status) level_filter
          if (read_status /= 0 .or. level_filter < levelmin .or. &
               level_filter > MAXLEVEL) level_filter = -1
       end if
       enabled_resolved = .true.
       if (enabled .and. myid == 1) write(*,'(A,F8.4)') &
            ' SNRT S_N RT enabled; reduced speed factor=', reduced_c
       if (enabled .and. myid == 1 .and. level_filter > 0) &
            write(*,'(A,I0)') ' SNRT RT level filter=', level_filter
    else
       env_value = ''
       call get_environment_variable('SNRT_RT_ENABLE', env_value, &
            length=env_length, status=env_status)
       enabled = env_status == 0 .and. env_length == 1 .and. &
            env_value(1:1) == '1'
    end if
    if (.not. enabled) return
    if (ilevel < levelmin .or. dtnew(ilevel) <= 0.0d0) return
    if (level_filter > 0 .and. ilevel /= level_filter) return
    if (.not. snrt_cuda_available()) return
    if (.not. allocated(uold)) return

    wall_start = omp_get_wtime()
    call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)
    dt_s = dtnew(ilevel) * scale_t
    dx_code = boxlen / dble(icoarse_max - icoarse_min + 1) * 0.5d0**ilevel
    if (scale_l <= 0.0d0 .or. scale_nH <= 0.0d0 .or. dt_s <= 0.0d0 .or. &
         dx_code <= 0.0d0) return
    cell_volume_code = dx_code**ndim
    cdt_over_dx = snrt_c_cgs * reduced_c * dt_s / (dx_code * scale_l)
    energy_index = ndim + 2
    t_setup = omp_get_wtime() - wall_start

    wall_start = omp_get_wtime()
    call snrt_amr_build_same_level_neighbors(ilevel, leaf_cell, leaf_slot, &
         neighbor, nleaf, n_interface_face)
    t_topology = omp_get_wtime() - wall_start
    if (nleaf == 0) then
       deallocate(leaf_cell, leaf_slot, neighbor)
       return
    end if
    allocate(optical_depth(nleaf,snrt_ngroups), neutral_hydrogen(nleaf), &
         absorbed_group(nleaf,snrt_ngroups))
    call snrt_angular_init(direction_dp, angular_weight)

    wall_start = omp_get_wtime()
    do i = 1, nleaf
       icell = leaf_cell(i)
       islot = leaf_slot(i)
       rho_code = max(0.0d0, uold(icell,1))
       neutral_fraction = 1.0d0
       if (allocated(snrt_neutral_fraction)) neutral_fraction = &
            max(0.0d0, min(1.0d0, real(snrt_neutral_fraction(islot),dp)))
       neutral_hydrogen(i) = real(rho_code*neutral_fraction,c_float)
       call snrt_nlte_optical_depth_groups(rho_code*neutral_fraction, scale_nH, &
            dt_s, snrt_group_cross_section_cm2, tau_dp, ierr, reduced_c)
       if (ierr /= 0) then
          optical_depth(i,:) = 0.0_c_float
       else
          optical_depth(i,:) = real(max(tau_dp,0.0d0),c_float)
       end if
    end do
    t_nlte = omp_get_wtime() - wall_start

    wall_start = omp_get_wtime()
    t_locator = 0.0d0
    t_budget = 0.0d0
    t_deposit = 0.0d0
    n_locator_calls = 0
    n_active_sources = 0
    ! dMsmbh is cumulative within a coarse step.  Account for only the
    ! increment not already injected into the photon state on this rank.
    if (accounted_step /= nstep_coarse) then
       if (allocated(accounted_mass)) deallocate(accounted_mass)
       if (nsink > 0) then
          allocate(accounted_mass(nsink))
          accounted_mass = 0.0d0
       end if
       accounted_step = nstep_coarse
    else if (nsink > 0 .and. (.not. allocated(accounted_mass) .or. &
         size(accounted_mass) < nsink)) then
       allocate(accounted_mass_new(nsink))
       accounted_mass_new = 0.0d0
       if (allocated(accounted_mass)) accounted_mass_new(1:size(accounted_mass)) = &
            accounted_mass
       if (allocated(accounted_mass)) deallocate(accounted_mass)
       call move_alloc(accounted_mass_new, accounted_mass)
    end if
    if (nsink > 0 .and. allocated(dMsmbh) .and. allocated(xsink)) then
       scale_m = scale_d * scale_l**3
       do isink = 1, nsink
          delta_accreted = max(0.0d0, dMsmbh(isink) - accounted_mass(isink))
          if (delta_accreted <= 0.0d0) cycle
          wall_sub = omp_get_wtime()
          call snrt_agn_find_local_leaf(xsink(isink,1:ndim), icell, ilevel_found)
          t_locator = t_locator + omp_get_wtime() - wall_sub
          n_locator_calls = n_locator_calls + 1
          if (icell == 0 .or. ilevel_found /= ilevel) cycle
          islot = snrt_state_get_slot(icell)
          if (islot <= 0) cycle
          n_active_sources = n_active_sources + 1
          epsilon_r = 0.1d0
          if (allocated(eps_sink) .and. size(eps_sink) >= isink) &
               epsilon_r = max(1.0d-6, min(0.99d0, eps_sink(isink)))
          source_ok = .true.
          do igroup = 1, snrt_ngroups
             wall_sub = omp_get_wtime()
             call snrt_agn_photon_budget(delta_accreted, scale_m, dt_s, &
                  epsilon_r, 0.5d0*snrt_group_energy_fraction(igroup), &
                  snrt_group_mean_energy_ev(igroup), luminosity, emitted_photons)
             t_budget = t_budget + omp_get_wtime() - wall_sub
             wall_sub = omp_get_wtime()
             call snrt_agn_deposit_isotropic(snrt_intensity, islot, igroup, &
                  emitted_photons, cell_volume_code, scale_l, scale_nH, &
                  angular_weight, deposited_density, ierr)
             t_deposit = t_deposit + omp_get_wtime() - wall_sub
             if (ierr /= 0) source_ok = .false.
          end do
          if (source_ok) accounted_mass(isink) = accounted_mass(isink) + delta_accreted
       end do
    end if
    t_source = omp_get_wtime() - wall_start
    t_source_overhead = t_source - t_locator - t_budget - t_deposit

    wall_start = omp_get_wtime()
    call snrt_transport_absorb_multigroup_prepared(leaf_slot, neighbor, &
         cdt_over_dx, optical_depth, neutral_hydrogen, absorbed_group, ierr)
    t_transport = omp_get_wtime() - wall_start
    if (ierr /= 0) then
       if (myid == 1) write(*,'(A,I0,A,I0)') &
            ' SNRT RT transport failed, code=', ierr, ' level=', ilevel
       deallocate(leaf_cell, leaf_slot, neighbor, optical_depth, &
            neutral_hydrogen, absorbed_group)
       return
    end if

    wall_start = omp_get_wtime()
    do i = 1, nleaf
       icell = leaf_cell(i)
       islot = leaf_slot(i)
       rho_code = max(0.0d0, uold(icell,1))
       if (rho_code <= 0.0d0) cycle
       neutral_fraction = 1.0d0
       if (allocated(snrt_neutral_fraction)) neutral_fraction = &
            max(0.0d0, min(1.0d0, real(snrt_neutral_fraction(islot),dp)))
       ionized_fraction = 1.0d0 - neutral_fraction
       heating_total = 0.0d0
       do igroup = 1, snrt_ngroups
          call snrt_nlte_photo_source(real(absorbed_group(i,igroup),dp), &
               rho_code, ionized_fraction, scale_nH, dt_s, &
               snrt_group_mean_energy_ev(igroup), ionization_increment, &
               heating_rate)
          ionized_fraction = min(1.0d0, ionized_fraction + ionization_increment)
          heating_total = heating_total + heating_rate
       end do
       if (allocated(snrt_neutral_fraction)) snrt_neutral_fraction(islot) = &
            real(max(0.0d0, 1.0d0-ionized_fraction),c_float)
       if (size(uold,1) >= icell .and. size(uold,2) >= energy_index .and. &
            scale_d*scale_v**2 > 0.0d0) then
          uold(icell,energy_index) = uold(icell,energy_index) + heating_total*dt_s / &
               (scale_d*scale_v**2)
       end if
    end do
    t_coupling = omp_get_wtime() - wall_start

    if (myid == 1) then
       write(*,'(A,I0)') ' SNRT source internals level=', ilevel
       write(*,'(A,F10.3)') '   locator  : ', t_locator
       write(*,'(A,F10.3)') '   photon   : ', t_budget
       write(*,'(A,F10.3)') '   deposit  : ', t_deposit
       write(*,'(A,F10.3)') '   overhead : ', t_source_overhead
       write(*,'(A,I0)') '   locator calls: ', n_locator_calls
       write(*,'(A,I0)') '   active sources: ', n_active_sources
       write(*,'(A,I0,A,I0,6(A,F10.3,1X))') &
         ' SNRT stage timings level=', ilevel, ' leaves=', nleaf, &
         ' setup=', t_setup, ' topology=', t_topology, ' nlte=', t_nlte, &
         ' source=', t_source, ' transport=', t_transport, &
         ' coupling=', t_coupling
    endif

    deallocate(leaf_cell, leaf_slot, neighbor, optical_depth, &
         neutral_hydrogen, absorbed_group)
  end subroutine snrt_ramses_advance_level

end module snrt_ramses_driver
