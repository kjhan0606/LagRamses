! Conservative state transfer at the actual grid creation/destruction hooks.
! This does not implement coarse-fine transport or load-balance migration.
module snrt_regrid
  use amr_commons, only: dp,ncoarse,ngridmax,myid,ncpu,active,father,headl,next
#ifndef WITHOUTMPI
  use mpi_mod
#endif
  use hydro_commons, only: uold
  use snrt_state, only: snrt_nslot,snrt_ndirection,snrt_checkpoint_cell_width, &
       snrt_state_get_slot,snrt_state_pack_cell,snrt_state_restore_cell, &
       snrt_state_clear_cell,validate_cell_payload
#ifdef DUST_LIVE
  use snrt_dust_contract, only: snrt_dust_contract_version,snrt_dust_contract_number_ir
  use snrt_dust_live, only: snrt_dust_live_pack,snrt_dust_live_restore
#endif
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#include "amr_index.h"
  implicit none
  private
  public :: snrt_regrid_refine,snrt_regrid_coarsen,snrt_regrid_check,snrt_regrid_upload
contains
  subroutine snrt_regrid_check(ierr)
    integer, intent(in) :: ierr
    integer :: info
    if(ierr==0)return
    write(*,'(A,I0,A,I0)')'SNRT regrid state transfer rejected: rank=',myid,' error=',ierr
#ifndef WITHOUTMPI
    ! Grid mutations are rank-local; MPI_FINALIZE here could strand peers.
    call MPI_ABORT(MPI_COMM_WORLD,10,info)
#endif
    error stop 'SNRT regrid state transfer failed'
  end subroutine

  pure integer function ir_width() result(n)
    n=0
#ifdef DUST_LIVE
    if(snrt_dust_contract_version>=3)n=snrt_dust_contract_number_ir*snrt_ndirection
#endif
  end function

  subroutine pack(icell,payload,ierr)
    integer, intent(in) :: icell
    real(dp), intent(out) :: payload(:)
    integer, intent(out) :: ierr
    call snrt_state_pack_cell(icell,payload(1:snrt_checkpoint_cell_width),ierr)
    if(ierr/=0)return
#ifdef DUST_LIVE
    if(ir_width()>0)call snrt_dust_live_pack(icell,payload(snrt_checkpoint_cell_width+1:),ierr)
#endif
  end subroutine

  subroutine restore(icell,payload,ierr)
    integer, intent(in) :: icell
    real(dp), intent(in) :: payload(:)
    integer, intent(out) :: ierr
    call snrt_state_restore_cell(icell,payload(1:snrt_checkpoint_cell_width),ierr)
    if(ierr/=0)return
#ifdef DUST_LIVE
    if(ir_width()>0)call snrt_dust_live_restore(icell,payload(snrt_checkpoint_cell_width+1:),ierr)
#endif
  end subroutine

  subroutine clear(icell,ierr)
    integer, intent(in) :: icell
    integer, intent(out) :: ierr
    real(dp) :: zero(ir_width())
    zero=0
    ierr=0
#ifdef DUST_LIVE
    if(ir_width()>0)call snrt_dust_live_restore(icell,zero,ierr)
#endif
    if(ierr==0)call snrt_state_clear_cell(icell)
  end subroutine

  subroutine cell_ids(parent,grid,children,ierr)
    integer, intent(in) :: parent,grid
    integer, intent(out) :: children(twotondim),ierr
    integer :: j
    ierr=1
    if(parent<1.or.parent>ICELL_OF(ngridmax,twotondim))return
    if(grid<1.or.grid>ngridmax)return
    do j=1,twotondim
       children(j)=ICELL_OF(grid,j)
    end do
    if(any(children==parent))return
    ierr=0
  end subroutine

  subroutine snrt_regrid_refine(parent,grid,ierr)
    integer, intent(in) :: parent,grid
    integer, intent(out) :: ierr
    real(dp) :: payload(snrt_checkpoint_cell_width+ir_width())
    integer :: children(twotondim),j
    ierr=0
    if(snrt_nslot==0)return
    call cell_ids(parent,grid,children,ierr)
    if(ierr/=0)return
    call pack(parent,payload,ierr)
    if(ierr/=0)return
    ! First-order, positivity-preserving prolongation of density variables.
    ! Each child has 1/8 the parent volume, so integrated photons/IR are kept.
    do j=1,twotondim
       if(payload(1)==0)then
          call clear(children(j),ierr)
       else
          call restore(children(j),payload,ierr)
       end if
       if(ierr/=0)return
    end do
  end subroutine

  subroutine snrt_regrid_coarsen(parent,grid,ierr)
    integer, intent(in) :: parent,grid
    integer, intent(out) :: ierr
    real(dp) :: merged(snrt_checkpoint_cell_width+ir_width())
    integer :: children(twotondim),j
    ierr=0
    if(snrt_nslot==0)return
    call cell_ids(parent,grid,children,ierr)
    if(ierr/=0)return
    call restricted_payload(parent,grid,merged,ierr)
    if(ierr/=0)return
    if(merged(1)==0)then
       call clear(parent,ierr)
    else
       call restore(parent,merged,ierr)
    end if
    if(ierr/=0)return
    do j=1,twotondim
       call clear(children(j),ierr)
       if(ierr/=0)return
    end do
  end subroutine

  subroutine restricted_payload(parent,grid,merged,ierr)
    integer, intent(in) :: parent,grid
    real(dp), intent(out) :: merged(:)
    integer, intent(out) :: ierr
    real(dp) :: payload(snrt_checkpoint_cell_width+ir_width(),twotondim)
    real(dp) :: weight(twotondim),norm
    integer :: children(twotondim),j
    ierr=0
    merged=0
    call cell_ids(parent,grid,children,ierr)
    if(ierr/=0)return
    do j=1,twotondim
       call pack(children(j),payload(:,j),ierr)
       if(ierr/=0)return
    end do
    ! Missing children cannot be interpreted as zero radiation if any child
    ! has evolved state. Detect this before modifying the parent or children.
    if(all(payload(1,:)==0))then
       call pack(parent,merged,ierr)
       if(ierr/=0)return
       if(any(merged(2:)/=0))then
          ierr=2
          return
       end if
       merged=0
       return
    end if
    ierr=2
    if(any(payload(1,:)/=1))return
    if(.not.allocated(uold))return
    if(maxval(children)>size(uold,1))return
    weight=uold(children,1)
    if(any(.not.ieee_is_finite(weight)).or.any(weight<0))return
    norm=maxval(weight)
    if(norm>0)then
       weight=weight/norm
       weight=weight/sum(weight)
    end if
    ! Photon/IR densities are volume averaged. Fractions instead preserve
    ! rho*x inventories; averaging fractions by volume loses ions when rho
    ! differs among children. Hydro owns the corresponding mass restriction.
    merged=sum(payload/real(twotondim,dp),dim=2)
    merged(1)=1
    do j=2,4
       merged(j)=min(1.0_dp,sum(payload(j,:)*weight))
    end do
    norm=max(1.0_dp,sum(merged(3:4)))
    merged(3:4)=merged(3:4)/norm
    call validate_cell_payload(merged(1:snrt_checkpoint_cell_width),ierr)
    if(ierr/=0)return
    if(any(.not.ieee_is_finite(merged)).or.any(merged<0))then
       ierr=2
       return
    end if
  end subroutine

  subroutine snrt_regrid_upload(ilevel,ierr)
    ! Restrict the completed fine level to its parents WITHOUT retiring the
    ! fine states. Parent cells may be virtual on the fine grid's rank:
    ! reverse communication sends each complete oct average to its owner.
    integer, intent(in) :: ilevel
    integer, intent(out) :: ierr
    integer :: nfield,ncache,width,i,j,g,parent,owner_grid,ind,nowned,info,global_error
    integer :: local_counts(2),global_counts(2)
    integer, allocatable :: parents(:),owned(:)
    real(dp), allocatable :: reduced(:,:),received(:,:),field(:),count_field(:)
    ierr=0
    if(ilevel<=1)return
    nfield=ICELL_OF(ngridmax,twotondim)
    ncache=active(ilevel)%ngrid
    width=snrt_checkpoint_cell_width+ir_width()
    allocate(parents(ncache),reduced(width,ncache),field(nfield),count_field(nfield))
    count_field=0
    do i=1,ncache
       g=active(ilevel)%igrid(i); parent=father(g); parents(i)=parent
       call restricted_payload(parent,g,reduced(:,i),j)
       ierr=max(ierr,j)
       if(j/=0)cycle
       count_field(parent)=count_field(parent)+1
    end do
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(ierr,global_error,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
    if(info/=0)call MPI_ABORT(MPI_COMM_WORLD,10,j)
    ierr=global_error
#endif
    if(ierr/=0)return
    if(ncpu>1)call make_virtual_reverse_dp(count_field,ilevel-1)
    ! Select owned parent cells through the persistent grid list, not a
    ! possibly absent active list below levelmin.
    nowned=0; owner_grid=headl(myid,ilevel-1)
    do while(owner_grid>0)
       do ind=1,twotondim
          parent=ICELL_OF(owner_grid,ind)
          if(count_field(parent)>0)nowned=nowned+1
       end do
       owner_grid=next(owner_grid)
    end do
    local_counts=[ncache,nowned]; global_counts=local_counts
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(local_counts,global_counts,2,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
    if(info/=0)call MPI_ABORT(MPI_COMM_WORLD,10,j)
#endif
    ! Every contributing oct must have exactly one receiving parent owner.
    if(global_counts(1)/=global_counts(2))then
       ierr=3
       return
    end if
    allocate(owned(nowned),received(width,nowned))
    i=0; owner_grid=headl(myid,ilevel-1)
    do while(owner_grid>0)
       do ind=1,twotondim
          parent=ICELL_OF(owner_grid,ind)
          if(count_field(parent)<=0)cycle
          i=i+1; owned(i)=parent
          if(count_field(parent)/=1)ierr=3
       end do
       owner_grid=next(owner_grid)
    end do
    do j=1,width
       field=0
       do i=1,ncache
          field(parents(i))=field(parents(i))+reduced(j,i)
       end do
       if(ncpu>1)call make_virtual_reverse_dp(field,ilevel-1)
       received(j,:)=field(owned)
    end do
    do i=1,nowned
       call validate_cell_payload(received(1:snrt_checkpoint_cell_width,i),j)
       ierr=max(ierr,j)
    end do
    if(any(.not.ieee_is_finite(received)).or.any(received<0))ierr=max(ierr,3)
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(ierr,global_error,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
    if(info/=0)call MPI_ABORT(MPI_COMM_WORLD,10,j)
    ierr=global_error
#endif
    if(ierr/=0)return
    do i=1,nowned
       if(received(1,i)==0)then
          call clear(owned(i),ierr)
       else
          call restore(owned(i),received(:,i),ierr)
       end if
       if(ierr/=0)return
    end do
  end subroutine
end module
