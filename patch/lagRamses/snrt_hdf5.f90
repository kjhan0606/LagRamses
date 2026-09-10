! Radiation follows the same file grid ordering as the HDF5 hydro payload.
! RAMSES memory cell identifiers are reconstructed, never serialized here.
module snrt_hdf5
  use dust_composition_optics, only: d03_band_sha256
  use dust_iron_optics, only: fe_band_sha256
  use amr_commons
  use ramses_hdf5_io
  use snrt_state, only: primary_width=>snrt_checkpoint_cell_width, snrt_ndirection, &
       primary_pack=>snrt_state_pack_cell, primary_restore=>snrt_state_restore_cell, &
       snrt_checkpoint_number_width, snrt_state_clear_cell
  use snrt_agn_efficiency, only: snrt_agn_rt_requested
  use snrt_stellar_source, only: stellar_sed_enabled, stellar_sed_identity
  use snrt_spectral_contract, only: snrt_spectral_contract_source_sha256, &
       snrt_spectral_contract_source_commit_binding, snrt_spectral_contract_approval_id, &
       snrt_spectral_contract_group_edges_sha256, snrt_spectral_contract_status, &
       snrt_spectral_contract_fraction_semantics,snrt_band_enabled,snrt_band_model,snrt_band_kind, &
       snrt_d03_band_enabled,snrt_chimes_band_enabled,snrt_chimes_bank_sha256
  use snrt_thermochemistry, only: snrt_secondary_loaded_manifest_sha256
#ifdef DUST_LIVE
  use snrt_dust_contract
  use snrt_dust_live, only: snrt_dust_live_pack, snrt_dust_live_restore
#endif
#include "amr_index.h"
  implicit none
  private
  integer :: snrt_checkpoint_cell_width=primary_width
  integer :: snrt_checkpoint_file_width=primary_width
  logical :: legacy_number_only=.false.
  public :: snrt_hdf5_write, snrt_hdf5_read
contains
  subroutine snrt_state_pack_cell(icell,payload,ierr)
    integer, intent(in) :: icell
    real(dp), intent(out) :: payload(:)
    integer, intent(out) :: ierr
    call primary_pack(icell,payload(1:primary_width),ierr)
    if(ierr/=0)return
#ifdef DUST_LIVE
    if(snrt_checkpoint_cell_width>primary_width) &
         call snrt_dust_live_pack(icell,payload(primary_width+1:),ierr)
#endif
  end subroutine

  subroutine snrt_state_restore_cell(icell,file_payload,ierr,validate_only)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    integer, intent(in) :: icell
    real(dp), intent(in) :: file_payload(:)
    integer, intent(out) :: ierr
    logical, intent(in) :: validate_only
    real(dp) :: payload(snrt_checkpoint_cell_width)
#ifdef DUST_LIVE
    real(dp) :: saved_ir(snrt_checkpoint_cell_width-primary_width)
#endif
    integer, parameter :: old_primary_width=primary_width-snrt_checkpoint_number_width
    ! File strides remain in the declared disk layout. Insert zero correction
    ! between old primary photons and IR, never reinterpret IR as a shift.
    ierr=10
    if(size(file_payload)/=snrt_checkpoint_file_width)return
    payload=0.0_dp
    if(legacy_number_only)then
       payload(1:old_primary_width)=file_payload(1:old_primary_width)
       payload(primary_width+1:)=file_payload(old_primary_width+1:)
    else
       payload=file_payload
    end if
    if(any(.not.ieee_is_finite(payload)))return
    if(any(payload(primary_width+1:)<0.0_dp))return
    if(payload(1)==0.and.any(payload/=0))return
    call primary_restore(icell,payload(1:primary_width),ierr,validate_only=.true.)
    if(ierr/=0)return
#ifdef DUST_LIVE
    ! Prepare/validate the IR contract before a primary cell can be published.
    ! Packing may reserve storage, but preserves every existing IR value.
    if(snrt_checkpoint_cell_width>primary_width)then
       call snrt_dust_live_pack(icell,saved_ir,ierr)
       if(ierr/=0)return
    end if
#endif
    if(validate_only)return
    call primary_restore(icell,payload(1:primary_width),ierr)
    if(ierr/=0)return
    if(payload(1)==0.0_dp)call snrt_state_clear_cell(icell)
#ifdef DUST_LIVE
    if(snrt_checkpoint_cell_width>primary_width) &
         call snrt_dust_live_restore(icell,payload(primary_width+1:),ierr)
#endif
  end subroutine

  subroutine require_ok(status)
    integer, intent(in) :: status
    integer :: global_status, info
    include 'mpif.h'
    call MPI_Allreduce(status,global_status,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
    if(global_status/=0.or.info/=0)then
       if(myid==1)write(*,*) 'ERROR: SNRT HDF5 radiation checkpoint rejected ',global_status
       call MPI_Abort(MPI_COMM_WORLD,10,info)
    end if
  end subroutine

  subroutine identity(grp,writing)
    integer(HID_T), intent(in) :: grp
    logical, intent(in) :: writing
    character(len=128) :: values(7), loaded
    character(len=20) :: names(7)
    integer :: k,status,width,version,file_version,stellar_switch,count
    logical :: exists
    real(dp),allocatable :: stellar_values(:),saved_stellar(:)
#ifdef DUST_LIVE
    real(dp), allocatable :: dust_values(:), saved_dust(:)
#endif
    snrt_checkpoint_cell_width=primary_width
    version=1
#ifdef DUST_LIVE
    if(snrt_dust_contract_version>=3)then
       snrt_checkpoint_cell_width=primary_width+snrt_dust_contract_number_ir*snrt_ndirection
       version=2
    end if
#endif
    ! An older executable must not accept a checkpoint whose future stellar
    ! photons it cannot reproduce. Formats 3/4 add the live stellar source
    ! identity to primary-only / primary+IR formats 1/2 respectively.
    if(stellar_sed_enabled)version=version+2
    ! Formats 5..8 add primary energy correction to formats 1..4. Source,
    ! dust and stellar identities remain mandatory for legacy migration.
    version=version+4
    ! Spectral closure gets a disjoint identity, not additional cell state.
    ! No migration from a fixed-SED run merely by changing an environment flag.
    version=version+8*snrt_band_kind()
    legacy_number_only=.false.
    snrt_checkpoint_file_width=snrt_checkpoint_cell_width
    names=[character(len=20)::'source_sha256','source_commit','approval','edges_sha256', &
         'spectral_status','fraction_semantics','secondary_manifest']
    values=[character(len=128)::snrt_spectral_contract_source_sha256, &
         snrt_spectral_contract_source_commit_binding,snrt_spectral_contract_approval_id, &
         snrt_spectral_contract_group_edges_sha256,snrt_spectral_contract_status, &
         snrt_spectral_contract_fraction_semantics,snrt_secondary_loaded_manifest_sha256]
    if(writing)then
       call hdf5_write_attr_int(grp,'cell_width',snrt_checkpoint_cell_width)
       call hdf5_write_attr_int(grp,'format_version',version)
       if(snrt_band_enabled())call hdf5_write_attr_string(grp,'band_model',snrt_band_model)
       if(snrt_d03_band_enabled())call hdf5_write_attr_string(grp,'d03_node_sha256',d03_band_sha256)
       if(snrt_band_kind()==4)call hdf5_write_attr_string(grp,'fe_node_sha256',fe_band_sha256)
       if(snrt_chimes_band_enabled())call hdf5_write_attr_string(grp,'chimes_bank_sha256',snrt_chimes_bank_sha256)
       call hdf5_write_attr_string(grp,'primary_shift_units','photon CODE density * eV per direction')
       if(snrt_checkpoint_cell_width>primary_width) &
            call hdf5_write_attr_string(grp,'ir_energy_units','erg/cm3 per normalized direction')
    else
       call hdf5_read_attr_int_checked(grp,'format_version',file_version,status)
       call require_ok(status)
       call require_ok(merge(0,1,file_version==version.or. &
            (file_version==version-4.and..not.snrt_band_enabled())))
       legacy_number_only=file_version==version-4.and..not.snrt_band_enabled()
       if(snrt_band_enabled())then
          call hdf5_read_attr_string_checked(grp,'band_model',loaded,status)
          call require_ok(status)
          call require_ok(merge(0,1,trim(loaded)==snrt_band_model))
          if(snrt_d03_band_enabled())then
             call hdf5_read_attr_string_checked(grp,'d03_node_sha256',loaded,status)
             call require_ok(status)
             call require_ok(merge(0,1,trim(loaded)==d03_band_sha256))
          endif
          if(snrt_band_kind()==4)then
             call hdf5_read_attr_string_checked(grp,'fe_node_sha256',loaded,status)
             call require_ok(status)
             call require_ok(merge(0,1,trim(loaded)==fe_band_sha256))
          endif
          if(snrt_chimes_band_enabled())then
             call hdf5_read_attr_string_checked(grp,'chimes_bank_sha256',loaded,status)
             call require_ok(status)
             call require_ok(merge(0,1,trim(loaded)==snrt_chimes_bank_sha256))
          endif
       endif
       if(legacy_number_only)snrt_checkpoint_file_width=snrt_checkpoint_cell_width-snrt_checkpoint_number_width
       call hdf5_read_attr_int_checked(grp,'cell_width',width,status)
       call require_ok(status)
       call require_ok(merge(0,1,width==snrt_checkpoint_file_width))
       if(.not.legacy_number_only)then
          call hdf5_read_attr_string_checked(grp,'primary_shift_units',loaded,status)
          call require_ok(status)
          call require_ok(merge(0,1,trim(loaded)=='photon CODE density * eV per direction'))
       end if
    end if
    do k=1,size(names)
       if(writing)then
          call hdf5_write_attr_string(grp,trim(names(k)),trim(values(k)))
       else
          call hdf5_read_attr_string_checked(grp,trim(names(k)),loaded,status)
          call require_ok(status)
          call require_ok(merge(0,1,trim(loaded)==trim(values(k))))
       end if
    end do
    ! Persist the actual age/Z photon table, not just its pathname. Disabled
    ! legacy checkpoints remain readable, but cannot acquire a new source
    ! on restart without an explicit migration of the scientific model.
    stellar_switch=merge(1,0,stellar_sed_enabled)
    if(writing)then
       call hdf5_write_attr_int(grp,'stellar_sed_enabled',stellar_switch)
    else
       call h5aexists_f(grp,'stellar_sed_enabled',exists,status)
       call require_ok(abs(status))
       width=0
       if(exists)then
          call hdf5_read_attr_int_checked(grp,'stellar_sed_enabled',width,status)
          call require_ok(status)
       endif
       call require_ok(merge(0,1,width==stellar_switch))
    endif
    if(stellar_sed_enabled)then
       call stellar_sed_identity(stellar_values)
       if(writing)then
          count=0
          if(myid==1)count=size(stellar_values)
          call hdf5_write_dataset_1d_dp(grp,'stellar_sed_values',stellar_values,count, &
               0_i8b,int(size(stellar_values),i8b))
       else
          allocate(saved_stellar(size(stellar_values)))
          call hdf5_read_dataset_1d_dp_checked(grp,'stellar_sed_values',saved_stellar,size(saved_stellar), &
               0_i8b,int(size(saved_stellar),i8b),status)
          call require_ok(status)
          call require_ok(merge(0,1,all(saved_stellar==stellar_values)))
       endif
    endif
#ifdef DUST_LIVE
    ! Bind the actual opacity/thermal values as well as the field map. A
    ! changed constant heat capacity would reinterpret the saved dust energy.
    dust_values=[real(snrt_dust_contract_version,dp), &
         real(snrt_dust_contract_number_groups,dp),real(snrt_dust_contract_number_temperature,dp), &
         snrt_dust_contract_mass_per_h_g,snrt_dust_contract_heat_capacity_per_h_erg_k, &
         snrt_dust_contract_group_edges_ev,snrt_dust_contract_absorption_per_h_cm2, &
         snrt_dust_contract_absorption_mean_energy_ev,snrt_dust_contract_temperature_k, &
         snrt_dust_contract_emitted_power_per_h_erg_s]
    ! Keep the v2 attribute unchanged. Version 3 additionally binds the IR
    ! quadrature, so a restart cannot reinterpret radiation with new opacity.
    if (snrt_dust_contract_version >= 3) dust_values=[dust_values, &
         real(snrt_dust_contract_number_ir,dp),snrt_dust_contract_ir_background_k, &
         snrt_dust_contract_ir_energy_ev,snrt_dust_contract_ir_weight_ev, &
         snrt_dust_contract_ir_absorption_per_h_cm2]
    ! Bind U(T) itself; metadata alone must not reinterpret saved material energy.
    if (snrt_dust_contract_version >= 4) dust_values=[dust_values, &
         snrt_dust_contract_internal_energy_per_h_erg]
    ! Absorption-only identities remain byte-identical. The extra extent and
    ! model marker reject switching scattering on/off or changing sigma at restart.
    if (snrt_dust_contract_scattering_enabled) dust_values=[dust_values, &
         5.0_dp,1.0_dp,snrt_dust_contract_scattering_per_h_cm2]
    if (snrt_dust_contract_exchange_enabled) dust_values=[dust_values, &
         6.0_dp,3.0_dp,snrt_dust_contract_collision_area_per_h,snrt_dust_contract_accommodation]
    if(writing)then
       call hdf5_write_attr_1d_dp(grp,'dust_contract_values',dust_values,size(dust_values))
    else
       allocate(saved_dust(size(dust_values)))
       call hdf5_read_attr_1d_dp_checked(grp,'dust_contract_values',saved_dust,size(saved_dust),status)
       call require_ok(status)
       call require_ok(merge(0,1,all(saved_dust==dust_values)))
    end if
#endif
  end subroutine

  subroutine snrt_hdf5_write()
    integer :: lev,grid,i,ind,nlocal,status,err,info,counts(ncpu),base,pass
    integer(i8b) :: total,offset
    integer(HID_T) :: grp
    character(len=32) :: name
    real(dp), allocatable :: buffer(:)
    include 'mpif.h'
    if(.not.snrt_agn_rt_requested())return
    snrt_checkpoint_cell_width=primary_width
#ifdef DUST_LIVE
    if(snrt_dust_contract_version>=3) &
         snrt_checkpoint_cell_width=primary_width+snrt_dust_contract_number_ir*snrt_ndirection
#endif
    ! Reject invalid radiation anywhere before publishing even the SNRT header.
    do pass=1,2
    if(pass==2)then
       call hdf5_create_group('/snrt',grp)
       call identity(grp,.true.)
    end if
    do lev=1,nlevelmax
       nlocal=numbl(myid,lev)
       call MPI_Allgather(nlocal,1,MPI_INTEGER,counts,1,MPI_INTEGER,MPI_COMM_WORLD,info)
       call require_ok(abs(info))
       call require_ok(merge(0,1,all(counts>=0)))
       total=sum(int(counts,i8b))*twotondim*snrt_checkpoint_cell_width
       if(total==0)cycle
       offset=sum(int(counts(1:myid-1),i8b))*twotondim*snrt_checkpoint_cell_width
       allocate(buffer(max(1,nlocal*twotondim*snrt_checkpoint_cell_width)))
       grid=headl(myid,lev)
       status=0
       do i=1,nlocal
          if(grid<1.or.grid>size(next))then
             status=1
             exit
          end if
          do ind=1,twotondim
             base=((i-1)*twotondim+ind-1)*snrt_checkpoint_cell_width
             call snrt_state_pack_cell(ICELL_OF(grid,ind), &
                  buffer(base+1:base+snrt_checkpoint_cell_width),err)
             status=max(status,err)
          end do
          grid=next(grid)
       end do
       if(grid/=0)status=max(status,1)
       call require_ok(status)
       write(name,'("level_",I0)')lev
       if(pass==2)call hdf5_write_dataset_1d_dp(grp,trim(name),buffer, &
            nlocal*twotondim*snrt_checkpoint_cell_width,offset,total)
       deallocate(buffer)
    end do
    end do
    call hdf5_close_group(grp)
  end subroutine

  subroutine snrt_hdf5_read()
    integer :: lev,grid,ind,i,nlocal,status,err,info,counts(ncpu),base,fidx,first,count,pass,handled,last_fidx
    integer(i8b) :: total,offset,total_grids
    integer(HID_T) :: grp
    character(len=32) :: name
    real(dp), allocatable :: buffer(:)
    include 'mpif.h'
    if(.not.snrt_agn_rt_requested())return
    call h5gopen_f(hdf5_file_id,'/snrt',grp,status)
    call require_ok(abs(status))
    call identity(grp,.false.)
    ! First stream the complete checkpoint through validation on every rank.
    ! Only the second pass publishes cells, so even corruption in a later
    ! chunk/level cannot leave an earlier cell partially restored.
    do pass=1,2
    do lev=1,nlevelmax
       nlocal=numbl(myid,lev)
       if(varcpu_restart)then
          total_grids=int(varcpu_ngrid_file(lev),i8b)
          offset=0_i8b
       else
          call MPI_Allgather(nlocal,1,MPI_INTEGER,counts,1,MPI_INTEGER,MPI_COMM_WORLD,info)
          call require_ok(abs(info))
          total_grids=sum(int(counts,i8b))
          offset=sum(int(counts(1:myid-1),i8b))
       end if
       call require_ok(merge(0,1,total_grids>=0.and.nlocal>=0.and.int(nlocal,i8b)<=total_grids))
       if(total_grids==0)then
          call require_ok(merge(0,1,nlocal==0.and.headl(myid,lev)==0))
          cycle
       end if
       total=total_grids*twotondim*snrt_checkpoint_file_width
       write(name,'("level_",I0)')lev
       grid=headl(myid,lev)
       handled=0
       last_fidx=0
       if(varcpu_restart)then
          ! Reuse the restored AMR file-grid map, in bounded streaming chunks.
          first=0
          do while(int(first,i8b)<total_grids)
             count=int(min(64_i8b,total_grids-int(first,i8b)))
             allocate(buffer(count*twotondim*snrt_checkpoint_file_width))
             call hdf5_read_dataset_1d_dp_checked(grp,trim(name),buffer,size(buffer), &
                  int(first,i8b)*twotondim*snrt_checkpoint_file_width,total,status)
             call require_ok(status)
             status=0
             do while(grid>0)
                if(grid>size(varcpu_grid_file_idx).or.grid>size(next).or.handled>=nlocal)then
                   status=1
                   exit
                end if
                fidx=varcpu_grid_file_idx(grid)
                if(fidx>first+count)exit
                if(fidx<=first.or.fidx<=last_fidx)then
                   status=1
                   exit
                end if
                do ind=1,twotondim
                   base=((fidx-first-1)*twotondim+ind-1)*snrt_checkpoint_file_width
                   call snrt_state_restore_cell(ICELL_OF(grid,ind), &
                        buffer(base+1:base+snrt_checkpoint_file_width),err,validate_only=pass==1)
                   status=max(status,err)
                end do
                last_fidx=fidx
                grid=next(grid)
                handled=handled+1
             end do
             call require_ok(status)
             deallocate(buffer)
             first=first+count
          end do
          call require_ok(merge(0,1,handled==nlocal.and.grid==0))
       else
          allocate(buffer(max(1,nlocal*twotondim*snrt_checkpoint_file_width)))
          call hdf5_read_dataset_1d_dp_checked(grp,trim(name),buffer, &
               nlocal*twotondim*snrt_checkpoint_file_width, &
               offset*twotondim*snrt_checkpoint_file_width,total,status)
          call require_ok(status)
          status=0
          do i=1,nlocal
             if(grid<1.or.grid>size(next))then
                status=1
                exit
             end if
             do ind=1,twotondim
                base=((i-1)*twotondim+ind-1)*snrt_checkpoint_file_width
                call snrt_state_restore_cell(ICELL_OF(grid,ind), &
                     buffer(base+1:base+snrt_checkpoint_file_width),err,validate_only=pass==1)
                status=max(status,err)
             end do
             grid=next(grid)
          end do
          if(grid/=0)status=max(status,1)
          call require_ok(status)
          deallocate(buffer)
       end if
    end do
    end do
    call hdf5_close_group(grp)
    if(myid==1)write(*,*)'SNRT_HDF5_RADIATION_RESTORE_PASS'
  end subroutine
end module
