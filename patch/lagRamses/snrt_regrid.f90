! Conservative state transfer at the actual grid creation/destruction hooks.
! This does not implement coarse-fine transport or load-balance migration.
module snrt_regrid
  use amr_commons, only: dp,ncoarse,ngridmax,myid
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
  public :: snrt_regrid_refine,snrt_regrid_coarsen,snrt_regrid_check
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
    real(dp) :: payload(snrt_checkpoint_cell_width+ir_width(),twotondim)
    real(dp) :: merged(snrt_checkpoint_cell_width+ir_width()),weight(twotondim),norm
    integer :: children(twotondim),j
    ierr=0
    if(snrt_nslot==0)return
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
       call clear(parent,ierr)
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
    call restore(parent,merged,ierr)
    if(ierr/=0)return
    do j=1,twotondim
       call clear(children(j),ierr)
       if(ierr/=0)return
    end do
  end subroutine
end module
