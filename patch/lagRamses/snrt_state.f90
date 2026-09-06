! Persistent local S_N state, indexed by the RAMSES leaf-cell identifier.
module snrt_state
  use amr_parameters, only: MAXLEVEL, dp
  use snrt_spectral_contract, only: snrt_spectral_ngroups => snrt_ngroups, &
       snrt_spectral_contract_loaded, snrt_spectral_contract_runtime_allowed, &
       snrt_spectral_contract_status, snrt_spectral_contract_source_id, &
       snrt_spectral_contract_source_sha256, &
       snrt_spectral_contract_source_commit_binding, &
       snrt_spectral_contract_approval_id, &
       snrt_spectral_contract_group_edges_sha256, &
       snrt_spectral_contract_interval_convention, &
       snrt_spectral_contract_fraction_semantics, &
       snrt_spectral_contract_checkpoint_identity_matches
  use snrt_thermochemistry, only: snrt_secondary_source_id, &
       snrt_secondary_upstream_commit, snrt_secondary_manifest_sha256, &
       snrt_secondary_tables_loaded, snrt_secondary_loaded_source_id, &
       snrt_secondary_loaded_upstream_commit, snrt_secondary_loaded_manifest_sha256
  use iso_c_binding, only: c_float
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#include "amr_index.h"
  implicit none

  integer, parameter, public :: snrt_ndirection = 80
  ! The spectral module owns the canonical nine-group dimensions.  Keep this
  ! local public alias so existing transport code retains its state-facing
  ! interface while all source-dependent moments come from the loaded
  ! contract rather than from a second hard-coded table.
  integer, parameter, public :: snrt_ngroups = snrt_spectral_ngroups
  integer, public :: snrt_nslot = 0
  real(c_float), allocatable, public :: snrt_intensity(:, :, :)
  real(dp), allocatable, public :: snrt_neutral_fraction(:)
  ! H II, He II, and He III are the authoritative persistent chemistry state;
  ! snrt_neutral_fraction remains a compatibility mirror for H I opacity.
  real(dp), allocatable, public :: snrt_hydrogen_ii(:)
  real(dp), allocatable, public :: snrt_helium_ii(:)
  real(dp), allocatable, public :: snrt_helium_iii(:)

  integer :: snrt_capacity = 0
  integer, allocatable :: snrt_slot_of_cell(:)
  integer, allocatable :: snrt_cell_id(:)

  private
  public :: snrt_state_sync_level, snrt_state_get_slot, snrt_state_get_cell
  public :: snrt_state_checkpoint_write, snrt_state_checkpoint_read
  integer, parameter, public :: snrt_checkpoint_cell_width = 4 + snrt_ndirection*snrt_ngroups
  public :: snrt_state_pack_cell, snrt_state_restore_cell

contains

  subroutine snrt_state_pack_cell(icell, payload, ierr)
    integer, intent(in) :: icell
    real(dp), intent(out) :: payload(snrt_checkpoint_cell_width)
    integer, intent(out) :: ierr
    integer :: slot
    payload = 0.0_dp
    ierr = 0
    slot = snrt_state_get_slot(icell)
    if(slot == 0) return
    payload(1:4) = [1.0_dp, snrt_hydrogen_ii(slot), snrt_helium_ii(slot), snrt_helium_iii(slot)]
    payload(5:) = reshape(real(snrt_intensity(:,:,slot),dp),[snrt_ndirection*snrt_ngroups])
    call validate_cell_payload(payload,ierr)
  end subroutine

  subroutine validate_cell_payload(payload,ierr)
    real(dp), intent(in) :: payload(snrt_checkpoint_cell_width)
    integer, intent(out) :: ierr
    ierr = 10
    if(any(.not.ieee_is_finite(payload)).or.any(payload<0.0_dp)) return
    if(payload(1)==0.0_dp)then
       if(any(payload/=0.0_dp)) return
    else if(payload(1)==1.0_dp)then
       if(payload(2)>1.0_dp.or.sum(payload(3:4))>1.0_dp+1.0d-10) return
       if(any(payload(5:)>real(huge(0.0_c_float),dp))) return
    else
       return
    end if
    ierr = 0
  end subroutine

  subroutine snrt_state_restore_cell(icell,payload,ierr)
    integer, intent(in) :: icell
    real(dp), intent(in) :: payload(snrt_checkpoint_cell_width)
    integer, intent(out) :: ierr
    integer :: slot
    call validate_cell_payload(payload,ierr)
    if(ierr/=0) return
    call snrt_state_initialize()
    if(icell<1.or.icell>size(snrt_slot_of_cell))then
       ierr=11
       return
    end if
    if(payload(1)==0.0_dp) return
    slot=snrt_slot_of_cell(icell)
    if(slot==0)then
       call snrt_state_grow(snrt_nslot+1)
       snrt_nslot=snrt_nslot+1
       slot=snrt_nslot
       snrt_slot_of_cell(icell)=slot
       snrt_cell_id(slot)=icell
    end if
    snrt_hydrogen_ii(slot)=payload(2)
    snrt_neutral_fraction(slot)=1.0_dp-payload(2)
    snrt_helium_ii(slot)=payload(3)
    snrt_helium_iii(slot)=payload(4)
    snrt_intensity(:,:,slot)=reshape(real(payload(5:),c_float),[snrt_ndirection,snrt_ngroups])
  end subroutine

  function snrt_state_get_slot(icell) result(islot)
    integer, intent(in) :: icell
    integer :: islot

    islot = 0
    if (.not. allocated(snrt_slot_of_cell)) return
    if (icell < 1 .or. icell > size(snrt_slot_of_cell)) return
    islot = snrt_slot_of_cell(icell)
  end function snrt_state_get_slot

  function snrt_state_get_cell(islot) result(icell)
    integer, intent(in) :: islot
    integer :: icell

    icell = 0
    if (.not. allocated(snrt_cell_id)) return
    if (islot < 1 .or. islot > snrt_nslot) return
    icell = snrt_cell_id(islot)
  end function snrt_state_get_cell

  subroutine snrt_state_checkpoint_write(unit_id, ierr)
    integer, intent(in) :: unit_id
    integer, intent(out) :: ierr
    integer, parameter :: checkpoint_version = 6
    integer :: ios
    character(len=128) :: secondary_source_id, secondary_upstream_commit
    character(len=128) :: secondary_manifest_sha256

    ierr = 0
    if (.not. snrt_spectral_contract_loaded .or. &
         .not. snrt_spectral_contract_runtime_allowed) then
       ierr = 5
       return
    end if
    if (.not. snrt_secondary_tables_loaded .or. &
         trim(snrt_secondary_loaded_source_id) /= trim(snrt_secondary_source_id) .or. &
         trim(snrt_secondary_loaded_upstream_commit) /= trim(snrt_secondary_upstream_commit) .or. &
         trim(snrt_secondary_loaded_manifest_sha256) /= trim(snrt_secondary_manifest_sha256)) then
       ierr = 8
       return
    end if
    if (snrt_nslot > 0 .and. (.not. allocated(snrt_neutral_fraction) .or. &
         .not. allocated(snrt_hydrogen_ii) .or. .not. allocated(snrt_helium_ii) .or. &
         .not. allocated(snrt_helium_iii))) then
       ierr = 6
       return
    end if
    if (snrt_nslot > 0 .and. &
         (.not. all(ieee_is_finite(snrt_neutral_fraction(1:snrt_nslot))) .or. &
          .not. all(ieee_is_finite(snrt_hydrogen_ii(1:snrt_nslot))) .or. &
          .not. all(ieee_is_finite(snrt_helium_ii(1:snrt_nslot))) .or. &
          .not. all(ieee_is_finite(snrt_helium_iii(1:snrt_nslot))) .or. &
          any(snrt_neutral_fraction(1:snrt_nslot) < 0.0d0) .or. &
          any(snrt_neutral_fraction(1:snrt_nslot) > 1.0d0) .or. &
          any(snrt_hydrogen_ii(1:snrt_nslot) < 0.0d0) .or. &
          any(snrt_hydrogen_ii(1:snrt_nslot) > 1.0d0) .or. &
          any(snrt_helium_ii(1:snrt_nslot) < 0.0d0) .or. &
          any(snrt_helium_iii(1:snrt_nslot) < 0.0d0) .or. &
          any(snrt_helium_ii(1:snrt_nslot) + snrt_helium_iii(1:snrt_nslot) > &
             1.0d0 + 1.0d-10) .or. &
          any(abs(snrt_neutral_fraction(1:snrt_nslot) + &
             snrt_hydrogen_ii(1:snrt_nslot) - 1.0d0) > 1.0d-10))) then
       ierr = 7
       return
    end if
    ! Invalid photons must be rejected before even the header is published.
    if (snrt_nslot > 0) then
       if (.not. allocated(snrt_intensity)) then
          ierr = 10
          return
       end if
       if (size(snrt_intensity,3) < snrt_nslot) then
          ierr = 10
          return
       end if
       if (any(.not. ieee_is_finite(snrt_intensity(:,:,1:snrt_nslot))) .or. &
            any(snrt_intensity(:,:,1:snrt_nslot) < 0.0_c_float)) then
          ierr = 10
          return
       end if
    end if
    write(unit_id, iostat=ios) checkpoint_version, snrt_ndirection, &
         snrt_ngroups, snrt_nslot
    if (ios /= 0) then
       ierr = 1
       return
    end if
    write(unit_id, iostat=ios) snrt_spectral_contract_status, &
         snrt_spectral_contract_source_id, &
         snrt_spectral_contract_source_sha256, &
         snrt_spectral_contract_source_commit_binding, &
         snrt_spectral_contract_approval_id, &
         snrt_spectral_contract_group_edges_sha256, &
         snrt_spectral_contract_interval_convention, &
         snrt_spectral_contract_fraction_semantics
    if (ios /= 0) then
       ierr = 2
       return
    end if
    secondary_source_id = snrt_secondary_source_id
    secondary_upstream_commit = snrt_secondary_upstream_commit
    secondary_manifest_sha256 = snrt_secondary_manifest_sha256
    write(unit_id, iostat=ios) secondary_source_id, secondary_upstream_commit, &
         secondary_manifest_sha256
    if (ios /= 0) then
       ierr = 2
       return
    end if
    if (snrt_nslot <= 0) return
    write(unit_id, iostat=ios) snrt_cell_id(1:snrt_nslot)
    if (ios /= 0) then
       ierr = 3
       return
    end if
    write(unit_id, iostat=ios) snrt_intensity(:,:,1:snrt_nslot)
    if (ios /= 0) then
       ierr = 4
       return
    end if
    write(unit_id, iostat=ios) snrt_neutral_fraction(1:snrt_nslot)
    if (ios /= 0) then
       ierr = 5
       return
    end if
    write(unit_id, iostat=ios) snrt_hydrogen_ii(1:snrt_nslot), &
         snrt_helium_ii(1:snrt_nslot), snrt_helium_iii(1:snrt_nslot)
    if (ios /= 0) ierr = 6
  end subroutine snrt_state_checkpoint_write

  subroutine snrt_state_checkpoint_read(unit_id, ierr)
    integer, intent(in) :: unit_id
    integer, intent(out) :: ierr
    integer, parameter :: checkpoint_version = 6
    integer :: ios, version, ndirection_file, ngroups_file, nslot_file
    integer :: islot, icell
    character(len=64) :: checkpoint_status
    character(len=128) :: checkpoint_source_id, checkpoint_source_sha256
    character(len=128) :: checkpoint_source_commit_binding, checkpoint_approval_id
    character(len=128) :: checkpoint_edges_sha256, checkpoint_interval_convention
    character(len=64) :: checkpoint_fraction_semantics
    character(len=128) :: checkpoint_secondary_source_id
    character(len=128) :: checkpoint_secondary_upstream_commit
    character(len=128) :: checkpoint_secondary_manifest_sha256
    integer, allocatable :: saved_cell_id(:)
    real(c_float), allocatable :: saved_intensity(:,:,:)
    real(dp), allocatable :: saved_neutral(:)
    real(dp), allocatable :: saved_hydrogen_ii(:), saved_helium_ii(:), &
         saved_helium_iii(:)

    ierr = 0
    read(unit_id, iostat=ios) version, ndirection_file, ngroups_file, nslot_file
    if (ios /= 0) then
       ierr = 1
       return
    end if
    if (version /= checkpoint_version .or. ndirection_file /= snrt_ndirection .or. &
         ngroups_file /= snrt_ngroups .or. nslot_file < 0) then
       ierr = 2
       return
    end if
    read(unit_id, iostat=ios) checkpoint_status, checkpoint_source_id, &
         checkpoint_source_sha256, checkpoint_source_commit_binding, &
         checkpoint_approval_id, checkpoint_edges_sha256, &
         checkpoint_interval_convention, checkpoint_fraction_semantics
    if (ios /= 0) then
       ierr = 3
       return
    end if
    if (.not. snrt_spectral_contract_checkpoint_identity_matches( &
         checkpoint_source_id, checkpoint_source_sha256, &
         checkpoint_source_commit_binding, checkpoint_approval_id, &
         checkpoint_edges_sha256, checkpoint_interval_convention, &
         checkpoint_fraction_semantics, checkpoint_status)) then
       ierr = 4
       return
    end if
    read(unit_id, iostat=ios) checkpoint_secondary_source_id, &
         checkpoint_secondary_upstream_commit, checkpoint_secondary_manifest_sha256
    if (ios /= 0) then
       ierr = 4
       return
    end if
    if (.not. snrt_secondary_tables_loaded .or. &
         trim(checkpoint_secondary_source_id) /= trim(snrt_secondary_source_id) .or. &
         trim(checkpoint_secondary_upstream_commit) /= trim(snrt_secondary_upstream_commit) .or. &
         trim(checkpoint_secondary_manifest_sha256) /= trim(snrt_secondary_manifest_sha256) .or. &
         trim(snrt_secondary_loaded_source_id) /= trim(snrt_secondary_source_id) .or. &
         trim(snrt_secondary_loaded_upstream_commit) /= trim(snrt_secondary_upstream_commit) .or. &
         trim(snrt_secondary_loaded_manifest_sha256) /= trim(snrt_secondary_manifest_sha256)) then
       ierr = 5
       return
    end if
    if (nslot_file == 0) then
       snrt_nslot = 0
       if (allocated(snrt_slot_of_cell)) snrt_slot_of_cell = 0
       return
    end if

    allocate(saved_cell_id(nslot_file), &
         saved_intensity(snrt_ndirection,snrt_ngroups,nslot_file), &
         saved_neutral(nslot_file), saved_hydrogen_ii(nslot_file), &
         saved_helium_ii(nslot_file), saved_helium_iii(nslot_file))
    read(unit_id, iostat=ios) saved_cell_id
    if (ios /= 0) then
       ierr = 6
       deallocate(saved_cell_id, saved_intensity, saved_neutral, &
            saved_hydrogen_ii, saved_helium_ii, saved_helium_iii)
       return
    end if
    read(unit_id, iostat=ios) saved_intensity
    if (ios /= 0) then
       ierr = 7
       deallocate(saved_cell_id, saved_intensity, saved_neutral, &
            saved_hydrogen_ii, saved_helium_ii, saved_helium_iii)
      return
    end if
    read(unit_id, iostat=ios) saved_neutral
    if (ios /= 0) then
       ierr = 8
       deallocate(saved_cell_id, saved_intensity, saved_neutral, &
            saved_hydrogen_ii, saved_helium_ii, saved_helium_iii)
       return
    end if
    read(unit_id, iostat=ios) saved_hydrogen_ii, saved_helium_ii, saved_helium_iii
    if (ios /= 0) then
       ierr = 8
       deallocate(saved_cell_id, saved_intensity, saved_neutral, &
            saved_hydrogen_ii, saved_helium_ii, saved_helium_iii)
       return
    end if
    if (any(.not. ieee_is_finite(saved_intensity)) .or. &
         any(saved_intensity < 0.0_c_float)) then
       ierr = 10
       deallocate(saved_cell_id, saved_intensity, saved_neutral, &
            saved_hydrogen_ii, saved_helium_ii, saved_helium_iii)
       return
    end if
    if (.not. all(ieee_is_finite(saved_neutral)) .or. &
         .not. all(ieee_is_finite(saved_hydrogen_ii)) .or. &
         .not. all(ieee_is_finite(saved_helium_ii)) .or. &
         .not. all(ieee_is_finite(saved_helium_iii)) .or. &
         any(saved_neutral < 0.0d0) .or. any(saved_neutral > 1.0d0) .or. &
         any(saved_hydrogen_ii < 0.0d0) .or. any(saved_hydrogen_ii > 1.0d0) .or. &
         any(saved_helium_ii < 0.0d0) .or. any(saved_helium_iii < 0.0d0) .or. &
         any(saved_helium_ii + saved_helium_iii > 1.0d0 + 1.0d-10) .or. &
         any(abs(saved_neutral + saved_hydrogen_ii - 1.0d0) > 1.0d-10)) then
       ierr = 9
       deallocate(saved_cell_id, saved_intensity, saved_neutral, &
            saved_hydrogen_ii, saved_helium_ii, saved_helium_iii)
       return
    end if

    if (.not. allocated(snrt_slot_of_cell)) call snrt_state_initialize()
    if (.not. allocated(snrt_slot_of_cell)) then
       ierr = 9
       deallocate(saved_cell_id, saved_intensity, saved_neutral, &
            saved_hydrogen_ii, saved_helium_ii, saved_helium_iii)
       return
    end if
    call snrt_state_grow(nslot_file)
    if (snrt_nslot > 0) snrt_intensity(:,:,1:snrt_nslot) = 0.0_c_float
    snrt_slot_of_cell = 0
    snrt_nslot = nslot_file
    snrt_cell_id(1:nslot_file) = saved_cell_id
    snrt_intensity(:,:,1:nslot_file) = saved_intensity
    snrt_neutral_fraction(1:nslot_file) = 1.0d0 - saved_hydrogen_ii
    snrt_hydrogen_ii(1:nslot_file) = saved_hydrogen_ii
    snrt_helium_ii(1:nslot_file) = saved_helium_ii
    snrt_helium_iii(1:nslot_file) = saved_helium_iii
    do islot = 1, nslot_file
       icell = saved_cell_id(islot)
       if (icell >= 1 .and. icell <= size(snrt_slot_of_cell)) &
            snrt_slot_of_cell(icell) = islot
    end do
    deallocate(saved_cell_id, saved_intensity, saved_neutral, saved_hydrogen_ii, &
         saved_helium_ii, saved_helium_iii)
  end subroutine snrt_state_checkpoint_read

  subroutine snrt_state_sync_level(ilevel, nleaf, nnew)
    use amr_commons, only: active, ncoarse, twotondim, son
    implicit none

    integer, intent(in) :: ilevel
    integer, intent(out) :: nleaf, nnew
    integer :: i, ind, igrid, icell, islot

    call snrt_state_initialize()
    nleaf = 0
    nnew = 0
    do i = 1, active(ilevel)%ngrid
       igrid = active(ilevel)%igrid(i)
       do ind = 1, twotondim
          icell = ICELL_OF(igrid, ind)
          if (son(icell) /= 0) cycle
          nleaf = nleaf + 1
          islot = snrt_slot_of_cell(icell)
          if (islot == 0) then
             call snrt_state_grow(snrt_nslot + 1)
             snrt_nslot = snrt_nslot + 1
             islot = snrt_nslot
             snrt_slot_of_cell(icell) = islot
             snrt_cell_id(islot) = icell
             snrt_intensity(:, :, islot) = 0.0_c_float
             snrt_neutral_fraction(islot) = 1.0d0
             snrt_hydrogen_ii(islot) = 0.0d0
             snrt_helium_ii(islot) = 0.0d0
             snrt_helium_iii(islot) = 0.0d0
             nnew = nnew + 1
          endif
       enddo
    enddo
  end subroutine snrt_state_sync_level

  subroutine snrt_state_initialize()
    use amr_commons, only: ncoarse, ngridmax, twotondim
    implicit none

    integer :: cell_capacity

    ! Keep zero-sized payloads allocated on ranks with no local leaves.  The
    ! native transaction and AMR exchange paths are collective; an unallocated
    ! persistent payload on an empty rank would make the rank leave the path
    ! before its peers.  snrt_state_grow replaces these zero-sized arrays when
    ! the first local slot is created.
    if (.not. allocated(snrt_slot_of_cell)) then
       cell_capacity = ICELL_OF(ngridmax, twotondim)
       allocate(snrt_slot_of_cell(cell_capacity))
       snrt_slot_of_cell = 0
    end if
    if (.not. allocated(snrt_intensity)) &
         allocate(snrt_intensity(snrt_ndirection, snrt_ngroups, max(0,snrt_capacity)))
    if (.not. allocated(snrt_neutral_fraction)) &
         allocate(snrt_neutral_fraction(max(0,snrt_capacity)))
    if (.not. allocated(snrt_hydrogen_ii)) &
         allocate(snrt_hydrogen_ii(max(0,snrt_capacity)))
    if (.not. allocated(snrt_helium_ii)) &
         allocate(snrt_helium_ii(max(0,snrt_capacity)))
    if (.not. allocated(snrt_helium_iii)) &
         allocate(snrt_helium_iii(max(0,snrt_capacity)))
    if (snrt_capacity == 0) then
       snrt_intensity = 0.0_c_float
       snrt_neutral_fraction = 1.0d0
       snrt_hydrogen_ii = 0.0d0
       snrt_helium_ii = 0.0d0
       snrt_helium_iii = 0.0d0
    end if
  end subroutine snrt_state_initialize

  subroutine snrt_state_grow(required)
    implicit none

    integer, intent(in) :: required
    integer :: next_capacity
    integer, allocatable :: next_cell_id(:)
    real(c_float), allocatable :: next_intensity(:, :, :)
    real(dp), allocatable :: next_neutral_fraction(:)
    real(dp), allocatable :: next_hydrogen_ii(:), next_helium_ii(:), &
         next_helium_iii(:)

    if (required <= snrt_capacity) return
    next_capacity = max(required, max(1024, 2 * snrt_capacity))
    allocate(next_cell_id(next_capacity))
    allocate(next_intensity(snrt_ndirection, snrt_ngroups, next_capacity))
    allocate(next_neutral_fraction(next_capacity))
    allocate(next_hydrogen_ii(next_capacity), next_helium_ii(next_capacity), &
         next_helium_iii(next_capacity))
    next_cell_id = 0
    next_intensity = 0.0_c_float
    next_neutral_fraction = 1.0d0
    next_hydrogen_ii = 0.0d0
    next_helium_ii = 0.0d0
    next_helium_iii = 0.0d0
    if (snrt_nslot > 0) then
       next_cell_id(1:snrt_nslot) = snrt_cell_id(1:snrt_nslot)
       next_intensity(:, :, 1:snrt_nslot) = snrt_intensity(:, :, 1:snrt_nslot)
       if (allocated(snrt_neutral_fraction)) then
          next_neutral_fraction(1:snrt_nslot) = &
               snrt_neutral_fraction(1:snrt_nslot)
       end if
       if (allocated(snrt_hydrogen_ii)) next_hydrogen_ii(1:snrt_nslot) = &
            snrt_hydrogen_ii(1:snrt_nslot)
       if (allocated(snrt_helium_ii)) next_helium_ii(1:snrt_nslot) = &
            snrt_helium_ii(1:snrt_nslot)
       if (allocated(snrt_helium_iii)) next_helium_iii(1:snrt_nslot) = &
            snrt_helium_iii(1:snrt_nslot)
    endif
    call move_alloc(next_cell_id, snrt_cell_id)
    call move_alloc(next_intensity, snrt_intensity)
    call move_alloc(next_neutral_fraction, snrt_neutral_fraction)
    call move_alloc(next_hydrogen_ii, snrt_hydrogen_ii)
    call move_alloc(next_helium_ii, snrt_helium_ii)
    call move_alloc(next_helium_iii, snrt_helium_iii)
    snrt_capacity = next_capacity
  end subroutine snrt_state_grow

end module snrt_state
