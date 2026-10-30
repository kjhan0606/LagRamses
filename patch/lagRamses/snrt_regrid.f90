! Conservative state transfer at the actual grid creation/destruction hooks.
! Coarse-fine transport remains separate; migration copies the typed primary
! and IR payload before the old communication/ownership maps are retired.
module snrt_regrid
  use amr_commons, only: dp,ncoarse,ngridmax,myid,ncpu,active,father,headl,next, &
       emission,reception,cpu_map,cpu_map2
#ifndef WITHOUTMPI
  use mpi_mod
#endif
  use hydro_commons, only: uold
  use amr_parameters, only: snrt_moment_order
  use snrt_state, only: snrt_nslot,snrt_ndirection,snrt_checkpoint_cell_width, &
       snrt_state_get_slot,snrt_state_pack_cell,snrt_state_restore_cell, &
       snrt_state_clear_cell,validate_cell_payload,snrt_state_get_cell,snrt_state_rebind_cells, &
       snrt_state_payload_tag,snrt_state_is_moment,snrt_ngroups
  use snrt_moment_live, only: mn_live_initialize
  use snrt_amr_topology, only: snrt_halo_tile_exchange
#ifdef DUST_LIVE
  use snrt_dust_contract, only: snrt_dust_contract_version,snrt_dust_contract_number_ir
  use snrt_dust_live, only: snrt_dust_live_pack,snrt_dust_live_restore,snrt_dust_live_validate_payload
#endif
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#include "amr_index.h"
  implicit none
  private
  public :: snrt_regrid_refine,snrt_regrid_coarsen,snrt_regrid_check,snrt_regrid_upload
  public :: snrt_regrid_defrag
  public :: snrt_regrid_migrate
contains
  subroutine snrt_regrid_migrate(ilevel,ierr)
    ! Expanded reception grids still belong to their OLD owners here. Copy
    ! only grids becoming locally owned, using the existing matched peer
    ! grid order. Never allocate capacity*full_angular_payload scratch.
    integer, intent(in) :: ilevel
    integer, intent(out) :: ierr
#ifndef WITHOUTMPI
    integer :: sc(ncpu),rc(ncpu),expected(ncpu),so(ncpu),ro(ncpu),req(2*ncpu)
    integer :: cpu,i,j,g,k,width,info,global_error,npeer,global_peers,batch
    integer :: rounds,global_rounds,round,nreq,ns,nr,first,last,err,base
    integer, allocatable :: sends(:),receives(:)
    real(dp), allocatable :: sbuf(:),rbuf(:)
    integer(kind=8), parameter :: scratch_bytes=64_8*1024_8*1024_8
    ierr=0
    if(ncpu<=1)return
    width=twotondim*(snrt_checkpoint_cell_width+ir_width())
    sc=0;rc=0
    do cpu=1,ncpu
       do i=1,emission(cpu,ilevel)%ngrid
          g=emission(cpu,ilevel)%igrid(i)
          if(cpu_map(father(g))==myid.and.cpu_map2(father(g))==cpu.and.cpu/=myid)sc(cpu)=sc(cpu)+1
       end do
       do i=1,reception(cpu,ilevel)%ngrid
          g=reception(cpu,ilevel)%igrid(i)
          if(cpu_map(father(g))==cpu.and.cpu_map2(father(g))==myid.and.cpu/=myid)rc(cpu)=rc(cpu)+1
       end do
    end do
    call MPI_ALLTOALL(sc,1,MPI_INTEGER,expected,1,MPI_INTEGER,MPI_COMM_WORLD,info)
    ierr=max(abs(info),merge(0,1,all(rc==expected)))
    call MPI_ALLREDUCE(ierr,global_error,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
    ierr=max(global_error,abs(info));if(ierr/=0)return
    npeer=count(sc>0)+count(rc>0)
    call MPI_ALLREDUCE(npeer,global_peers,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
    ierr=abs(info);if(ierr/=0.or.global_peers==0)return
    batch=max(1,min(8,int(scratch_bytes/(8_8*int(width,8)*global_peers))))
    rounds=(max(maxval(sc),maxval(rc))+batch-1)/batch
    call MPI_ALLREDUCE(rounds,global_rounds,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
    ierr=abs(info);if(ierr/=0)return
    ns=0;nr=0
    do cpu=1,ncpu
       so(cpu)=ns;ro(cpu)=nr;ns=ns+sc(cpu);nr=nr+rc(cpu)
    end do
    allocate(sends(ns),receives(nr))
    do cpu=1,ncpu
       k=so(cpu)
       do i=1,emission(cpu,ilevel)%ngrid
          g=emission(cpu,ilevel)%igrid(i)
          if(cpu_map(father(g))/=myid.or.cpu_map2(father(g))/=cpu.or.cpu==myid)cycle
          k=k+1;sends(k)=g
       end do
       k=ro(cpu)
       do i=1,reception(cpu,ilevel)%ngrid
          g=reception(cpu,ilevel)%igrid(i)
          if(cpu_map(father(g))/=cpu.or.cpu_map2(father(g))/=myid.or.cpu==myid)cycle
          k=k+1;receives(k)=g
       end do
    end do
    allocate(sbuf(max(1,sum(min(sc,batch))*width)),rbuf(max(1,sum(min(rc,batch))*width)))
    do round=1,global_rounds
       first=(round-1)*batch+1;k=0
       do cpu=1,ncpu
          last=min(sc(cpu),first+batch-1)
          do i=first,last
             g=sends(so(cpu)+i)
             do j=1,twotondim
                call pack(ICELL_OF(g,j),sbuf(k+1:k+width/twotondim),err)
                ierr=max(ierr,err);k=k+width/twotondim
             end do
          end do
       end do
       call MPI_ALLREDUCE(ierr,global_error,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
       ierr=max(global_error,abs(info));if(ierr/=0)return
       nreq=0;k=0
       do cpu=1,ncpu
          nr=max(0,min(batch,rc(cpu)-first+1))*width
          if(nr==0)cycle
          nreq=nreq+1
          call MPI_IRECV(rbuf(k+1),nr,MPI_DOUBLE_PRECISION,cpu-1,937,MPI_COMM_WORLD,req(nreq),info)
          if(info/=0)call MPI_ABORT(MPI_COMM_WORLD,10,err)
          k=k+nr
       end do
       k=0
       do cpu=1,ncpu
          ns=max(0,min(batch,sc(cpu)-first+1))*width
          if(ns==0)cycle
          nreq=nreq+1
          call MPI_ISEND(sbuf(k+1),ns,MPI_DOUBLE_PRECISION,cpu-1,937,MPI_COMM_WORLD,req(nreq),info)
          if(info/=0)call MPI_ABORT(MPI_COMM_WORLD,10,err)
          k=k+ns
       end do
       call MPI_WAITALL(nreq,req,MPI_STATUSES_IGNORE,info)
       ierr=abs(info)
       call MPI_ALLREDUCE(ierr,global_error,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
       ierr=max(global_error,abs(info));if(ierr/=0)return
       ! Validate the whole received batch on every rank before publication.
       do base=1,2
          k=0
          do cpu=1,ncpu
             last=min(rc(cpu),first+batch-1)
             do i=first,last
                g=receives(ro(cpu)+i)
                do j=1,twotondim
                   if(base==1)then
                      call validate_restore(ICELL_OF(g,j),rbuf(k+1:k+width/twotondim),err)
                   else if(rbuf(k+1)==0)then
                      call clear(ICELL_OF(g,j),err)
                   else
                      call restore(ICELL_OF(g,j),rbuf(k+1:k+width/twotondim),err,prevalidated=.true.)
                   endif
                   ierr=max(ierr,err);k=k+width/twotondim
                end do
             end do
          end do
          call MPI_ALLREDUCE(ierr,global_error,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
          ierr=max(global_error,abs(info));if(ierr/=0)return
       end do
    end do
#else
    ierr=0
#endif
  end subroutine

  subroutine snrt_regrid_defrag(grid_map,ierr)
    integer, intent(in) :: grid_map(:)
    integer, intent(out) :: ierr
    integer, allocatable :: cells(:)
    logical, allocatable :: retired(:)
    integer :: slot,cell,grid,child,packed,unused
    ierr=0
    if(snrt_nslot==0)return
    allocate(cells(snrt_nslot),retired(snrt_nslot)); retired=.false.
    packed=0
    if(size(grid_map)>0)packed=maxval(grid_map)
    ierr=1
    if(any(grid_map<0).or.packed>ngridmax)return
    unused=0
    do slot=1,snrt_nslot
       cell=snrt_state_get_cell(slot)
       if(cell<1.or.cell>ICELL_OF(ngridmax,twotondim))return
       cells(slot)=cell
       if(cell<=ncoarse)cycle
       grid=IGRID_OF(cell); child=ICHILD_OF(cell)
       if(grid<=size(grid_map))then
          if(grid_map(grid)>0)then
             cells(slot)=ICELL_OF(grid_map(grid),child)
             cycle
          end if
       end if
       ! Retain zero-valued retired slots on free cells, as ordinary grid
       ! deletion already does. This preserves the raw checkpoint's unique,
       ! nonzero cell IDs and permits later pool reuse without slot growth.
       grid=packed+unused/twotondim+1; child=mod(unused,twotondim)+1
       if(grid>ngridmax)return
       cells(slot)=ICELL_OF(grid,child); retired(slot)=.true.; unused=unused+1
    end do
    call snrt_state_rebind_cells(cells,ierr)
    if(ierr/=0)return
    do slot=1,snrt_nslot
       if(.not.retired(slot))cycle
       call clear(cells(slot),ierr)
       if(ierr/=0)return
    end do
  end subroutine

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

  subroutine restore(icell,payload,ierr,prevalidated)
    integer, intent(in) :: icell
    real(dp), intent(in) :: payload(:)
    integer, intent(out) :: ierr
    logical,optional,intent(in) :: prevalidated
    logical :: trusted
    ierr=2
    if(size(payload)/=snrt_checkpoint_cell_width+ir_width())return
    if(any(.not.ieee_is_finite(payload)))return
    ! PREVALIDATED is legal only after this exact payload passed the earlier
    ! full validation phase (and MPI agreement where the operation is collective).
    trusted=.false.
    if(present(prevalidated))trusted=prevalidated
    ! Empty payloads use CLEAR at regrid call sites; never let the trusted path
    ! reinterpret an invalid zero-tag payload as a valid stored state.
    if(trusted.and.payload(1)==0.0_dp)trusted=.false.
    if(.not.trusted)then
       call validate_restore(icell,payload,ierr)
       if(ierr/=0)return
    endif
    call snrt_state_restore_cell(icell,payload(1:snrt_checkpoint_cell_width),ierr,prevalidated=trusted)
    if(ierr/=0)return
#ifdef DUST_LIVE
    if(ir_width()>0)call snrt_dust_live_restore(icell,payload(snrt_checkpoint_cell_width+1:),ierr, &
         prevalidated=trusted)
#endif
  end subroutine

  subroutine validate_restore(icell,payload,ierr)
    integer, intent(in) :: icell
    real(dp), intent(in) :: payload(:)
    integer, intent(out) :: ierr
    ierr=2
    if(size(payload)/=snrt_checkpoint_cell_width+ir_width())return
    if(any(.not.ieee_is_finite(payload)))return
    if(payload(1)==0.0_dp.and.any(payload/=0.0_dp))return
    call snrt_state_restore_cell(icell,payload(1:snrt_checkpoint_cell_width),ierr,validate_only=.true.)
#ifdef DUST_LIVE
    if(ierr==0.and.ir_width()>0)call snrt_dust_live_validate_payload(payload(snrt_checkpoint_cell_width+1:),ierr)
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
    ! All children receive the same payload, and cell_ids constructs each from
    ! this valid grid. Run the expensive M_N realizability check only once.
    call validate_restore(children(1),payload,ierr)
    if(ierr/=0)return
    ! First-order, positivity-preserving prolongation of density variables.
    ! Each child has 1/8 the parent volume, so integrated photons/IR are kept.
    do j=1,twotondim
       if(payload(1)==0)then
          call clear(children(j),ierr)
       else
          call restore(children(j),payload,ierr,prevalidated=.true.)
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
       call restore(parent,merged,ierr,prevalidated=.true.)
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
    if(any(payload(1,:)/=snrt_state_payload_tag()))return
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
    merged(1)=snrt_state_payload_tag()
    do j=2,4
       merged(j)=min(1.0_dp,sum(payload(j,:)*weight))
    end do
    norm=max(1.0_dp,sum(merged(3:4)))
    merged(3:4)=merged(3:4)/norm
    ! Number and signed correction densities use the same volume weights.
    ! Validate the eventual stored FP32 number as well as the FP64 average.
    call validate_restore(parent,merged,ierr)
  end subroutine

  subroutine snrt_regrid_upload(ilevel,ierr)
    ! Restrict the completed fine level to its parents WITHOUT retiring the
    ! fine states. Parent cells may be virtual on the fine grid's rank:
    ! reverse communication sends each complete oct average to its owner.
    integer, intent(in) :: ilevel
    integer, intent(out) :: ierr
    integer :: nfield,ncache,width,i,j,g,parent,owner_grid,ind,nowned,info,global_error
    integer :: first_column,column,ncolumn,tile_width
    integer(kind=8),parameter :: tile_extra_bytes=8_8*1024_8*1024_8
    integer :: local_counts(2),global_counts(2)
    integer, allocatable :: parents(:),owned(:)
    real(dp), allocatable :: reduced(:,:),received(:,:),field(:),count_field(:),tile(:,:)
    ierr=0
    if(ilevel<=1)return
    nfield=ICELL_OF(ngridmax,twotondim)
    ncache=active(ilevel)%ngrid
    width=snrt_checkpoint_cell_width+ir_width()
    allocate(parents(ncache),reduced(width,ncache),field(nfield),count_field(nfield))
    ! Validation below is read-only and cell independent. Initialize the shared
    ! M_N basis serially first so simultaneous validators never race lazy setup.
    if(snrt_state_is_moment())then
       call mn_live_initialize(snrt_moment_order,snrt_ngroups,ierr)
       if(ierr/=0)return
    endif
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
    ! The scalar reverse exchange performed one MPI message/handshake per
    ! photon or IR component. Pack the same parent densities into bounded
    ! tiles and use the existing reverse-add halo map. No reduction order
    ! changes within a component; the final validation/commit stays below.
    tile_width=max(16,min(width,int(tile_extra_bytes/(8_8*max(1_8,int(nfield,8))))))
    allocate(tile(nfield,tile_width))
    do first_column=1,width,tile_width
       ncolumn=min(tile_width,width-first_column+1)
       tile(:,1:ncolumn)=0
       do column=1,ncolumn
          j=first_column+column-1
          do i=1,ncache
             tile(parents(i),column)=tile(parents(i),column)+reduced(j,i)
          end do
       end do
       if(ncpu>1)then
          call snrt_halo_tile_exchange(tile(:,1:ncolumn),ilevel-1,j,reverse=.true.)
          if(j/=0)then
             ierr=j
             return
          end if
       end if
       do column=1,ncolumn
          received(first_column+column-1,:)=tile(owned,column)
       end do
    end do
    ! Regrid validation can require hundreds of safeguarded Newton iterations
    ! per group near the realizability boundary. Balance independent parents
    ! dynamically, and keep the global all-or-nothing commit after this barrier.
!$omp parallel do default(shared) private(i,j) reduction(max:ierr) schedule(dynamic,1) if(nowned>=4)
    do i=1,nowned
       call validate_restore(owned(i),received(:,i),j)
       ierr=max(ierr,j)
    end do
!$omp end parallel do
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
          call restore(owned(i),received(:,i),ierr,prevalidated=.true.)
       end if
       if(ierr/=0)return
    end do
  end subroutine
end module
