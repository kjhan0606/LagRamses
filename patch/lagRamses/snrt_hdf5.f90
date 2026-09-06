! Radiation follows the same file grid ordering as the HDF5 hydro payload.
! RAMSES memory cell identifiers are reconstructed, never serialized here.
module snrt_hdf5
  use amr_commons
  use ramses_hdf5_io
  use snrt_state, only: snrt_checkpoint_cell_width, snrt_state_pack_cell, snrt_state_restore_cell
  use snrt_agn_efficiency, only: snrt_agn_rt_requested
  use snrt_spectral_contract, only: snrt_spectral_contract_source_sha256, &
       snrt_spectral_contract_source_commit_binding, snrt_spectral_contract_approval_id, &
       snrt_spectral_contract_group_edges_sha256, snrt_spectral_contract_status, &
       snrt_spectral_contract_fraction_semantics
  use snrt_thermochemistry, only: snrt_secondary_loaded_manifest_sha256
#ifdef DUST_LIVE
  use snrt_dust_contract
#endif
#include "amr_index.h"
  implicit none
  private
  public :: snrt_hdf5_write, snrt_hdf5_read
contains
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
    integer :: k,status,width
#ifdef DUST_LIVE
    real(dp), allocatable :: dust_values(:), saved_dust(:)
#endif
    names=[character(len=20)::'source_sha256','source_commit','approval','edges_sha256', &
         'spectral_status','fraction_semantics','secondary_manifest']
    values=[character(len=128)::snrt_spectral_contract_source_sha256, &
         snrt_spectral_contract_source_commit_binding,snrt_spectral_contract_approval_id, &
         snrt_spectral_contract_group_edges_sha256,snrt_spectral_contract_status, &
         snrt_spectral_contract_fraction_semantics,snrt_secondary_loaded_manifest_sha256]
    if(writing)then
       call hdf5_write_attr_int(grp,'cell_width',snrt_checkpoint_cell_width)
       call hdf5_write_attr_int(grp,'format_version',1)
    else
       call hdf5_read_attr_int_checked(grp,'format_version',width,status)
       call require_ok(status)
       call require_ok(merge(0,1,width==1))
       call hdf5_read_attr_int_checked(grp,'cell_width',width,status)
       call require_ok(status)
       call require_ok(merge(0,1,width==snrt_checkpoint_cell_width))
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
    integer :: lev,grid,i,ind,nlocal,status,err,info,counts(ncpu),base
    integer(i8b) :: total,offset
    integer(HID_T) :: grp
    character(len=32) :: name
    real(dp), allocatable :: buffer(:)
    include 'mpif.h'
    if(.not.snrt_agn_rt_requested())return
    call hdf5_create_group('/snrt',grp)
    call identity(grp,.true.)
    do lev=1,nlevelmax
       nlocal=numbl(myid,lev)
       call MPI_Allgather(nlocal,1,MPI_INTEGER,counts,1,MPI_INTEGER,MPI_COMM_WORLD,info)
       total=sum(int(counts,i8b))*twotondim*snrt_checkpoint_cell_width
       if(total==0)cycle
       offset=sum(int(counts(1:myid-1),i8b))*twotondim*snrt_checkpoint_cell_width
       allocate(buffer(max(1,nlocal*twotondim*snrt_checkpoint_cell_width)))
       grid=headl(myid,lev)
       status=0
       do i=1,nlocal
          do ind=1,twotondim
             base=((i-1)*twotondim+ind-1)*snrt_checkpoint_cell_width
             call snrt_state_pack_cell(ICELL_OF(grid,ind), &
                  buffer(base+1:base+snrt_checkpoint_cell_width),err)
             status=max(status,err)
          end do
          grid=next(grid)
       end do
       call require_ok(status)
       write(name,'("level_",I0)')lev
       call hdf5_write_dataset_1d_dp(grp,trim(name),buffer, &
            nlocal*twotondim*snrt_checkpoint_cell_width,offset,total)
       deallocate(buffer)
    end do
    call hdf5_close_group(grp)
  end subroutine

  subroutine snrt_hdf5_read()
    integer :: lev,grid,ind,i,nlocal,status,err,info,counts(ncpu),base,fidx,first,count
    integer(i8b) :: total,offset,total_grids
    integer(HID_T) :: grp
    character(len=32) :: name
    real(dp), allocatable :: buffer(:)
    include 'mpif.h'
    if(.not.snrt_agn_rt_requested())return
    call h5gopen_f(hdf5_file_id,'/snrt',grp,status)
    call require_ok(abs(status))
    call identity(grp,.false.)
    do lev=1,nlevelmax
       nlocal=numbl(myid,lev)
       if(varcpu_restart)then
          total_grids=int(varcpu_ngrid_file(lev),i8b)
          offset=0_i8b
       else
          call MPI_Allgather(nlocal,1,MPI_INTEGER,counts,1,MPI_INTEGER,MPI_COMM_WORLD,info)
          total_grids=sum(int(counts,i8b))
          offset=sum(int(counts(1:myid-1),i8b))
       end if
       if(total_grids==0)cycle
       total=total_grids*twotondim*snrt_checkpoint_cell_width
       write(name,'("level_",I0)')lev
       grid=headl(myid,lev)
       if(varcpu_restart)then
          ! Reuse the restored AMR file-grid map, in bounded streaming chunks.
          first=0
          do while(int(first,i8b)<total_grids)
             count=int(min(64_i8b,total_grids-int(first,i8b)))
             allocate(buffer(count*twotondim*snrt_checkpoint_cell_width))
             call hdf5_read_dataset_1d_dp_checked(grp,trim(name),buffer,size(buffer), &
                  int(first,i8b)*twotondim*snrt_checkpoint_cell_width,total,status)
             call require_ok(status)
             status=0
             do while(grid>0)
                fidx=varcpu_grid_file_idx(grid)
                if(fidx>first+count)exit
                if(fidx<=first)then
                   status=1
                   exit
                end if
                do ind=1,twotondim
                   base=((fidx-first-1)*twotondim+ind-1)*snrt_checkpoint_cell_width
                   call snrt_state_restore_cell(ICELL_OF(grid,ind), &
                        buffer(base+1:base+snrt_checkpoint_cell_width),err)
                   status=max(status,err)
                end do
                grid=next(grid)
             end do
             call require_ok(status)
             deallocate(buffer)
             first=first+count
          end do
       else
          allocate(buffer(max(1,nlocal*twotondim*snrt_checkpoint_cell_width)))
          call hdf5_read_dataset_1d_dp_checked(grp,trim(name),buffer, &
               nlocal*twotondim*snrt_checkpoint_cell_width, &
               offset*twotondim*snrt_checkpoint_cell_width,total,status)
          call require_ok(status)
          status=0
          do i=1,nlocal
             do ind=1,twotondim
                base=((i-1)*twotondim+ind-1)*snrt_checkpoint_cell_width
                call snrt_state_restore_cell(ICELL_OF(grid,ind), &
                     buffer(base+1:base+snrt_checkpoint_cell_width),err)
                status=max(status,err)
             end do
             grid=next(grid)
          end do
          call require_ok(status)
          deallocate(buffer)
       end if
    end do
    call hdf5_close_group(grp)
    if(myid==1)write(*,*)'SNRT_HDF5_RADIATION_RESTORE_PASS'
  end subroutine
end module
