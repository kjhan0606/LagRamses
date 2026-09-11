module stellar_radioactive_runtime
  use amr_commons
  use hydro_commons
  use stellar_enrichment_config, only: configured_radioactive_model
  use stellar_radioactive_decay, only: radioactive_gas_decay
  use stellar_radioactive_transport, only: radioactive_limit_states
  use stellar_native_units, only: code_interval_to_age_gyr,seconds_per_gyr
#include "amr_index.h"
  implicit none
  private
  public :: radioactive_advance_level,radioactive_limit_faces,radioactive_limit_children
contains
  subroutine radioactive_advance_level(ilevel)
    ! Called BEFORE set_unew and finer recursion. Otherwise already-applied
    ! fine fluxes could remove material that a later coarse decay consumes.
    ! Only owned leaves evolve: restriction averages their updated subsets.
    integer,intent(in)::ilevel
    integer::n,k,i,ind,cell,status,bad,all_bad,info,e
    integer,allocatable::cells(:)
    real(dp),allocatable::stage(:,:)
    real(dp)::sl,st,sd,sv,snh,st2,dt_gyr,elements(11),parent(2)
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif
    if(trim(configured_radioactive_model)=='none')return
    call units(sl,st,sd,sv,snh,st2)
    call code_interval_to_age_gyr(dtnew(ilevel),st,aexp,dt_gyr,status)
    bad=0
    if(status/=0.or.iradioactive<1.or.iradioactive+1>nvar.or.ichem<1)bad=1
    n=active(ilevel)%ngrid*twotondim
    allocate(cells(n),stage(4,n));k=0
    if(bad==0)then
       do ind=1,twotondim
          do i=1,active(ilevel)%ngrid
             cell=ICELL_OF(active(ilevel)%igrid(i),ind)
             if(son(cell)/=0)cycle
             elements=uold(cell,ichem:ichem+10)
             parent=uold(cell,iradioactive:iradioactive+1)
             call radioactive_gas_decay(uold(cell,imetal),elements,parent,dt_gyr*seconds_per_gyr,status)
             if(status/=0)then
                if(bad==0)write(*,*)'Radioactive gas state rejected: level,cell,status ',ilevel,cell,status
                bad=1;cycle
             endif
             k=k+1;cells(k)=cell;stage(:,k)=[elements(7),elements(11),parent]
          enddo
       enddo
    endif
    all_bad=bad;info=0
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(bad,all_bad,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
#endif
    if(all_bad/=0.or.info/=0)then
       if(myid==1)write(*,*)'ERROR: radioactive gas decay rejected before level commit'
       call clean_stop;return
    endif
    do i=1,k
       cell=cells(i)
       uold(cell,ichem+6)=stage(1,i);uold(cell,ichem+10)=stage(2,i)
       uold(cell,iradioactive:iradioactive+1)=stage(3:4,i)
    enddo
    ! Communicate the same decayed states used by subsequent hydro/fine buffers.
    call make_virtual_fine_dp(uold(:,ichem+6),ilevel)
    call make_virtual_fine_dp(uold(:,ichem+10),ilevel)
    do e=iradioactive,iradioactive+1
       call make_virtual_fine_dp(uold(:,e),ilevel)
    enddo
  end subroutine

  subroutine radioactive_limit_faces(qin,qm,qp,ngrid)
    integer,intent(in)::ngrid
    real(dp),intent(in)::qin(nvector,iu1:iu2,ju1:ju2,ku1:ku2,nvar)
    real(dp),intent(inout)::qm(nvector,iu1:iu2,ju1:ju2,ku1:ku2,nvar,ndim)
    real(dp),intent(inout)::qp(nvector,iu1:iu2,ju1:ju2,ku1:ku2,nvar,ndim)
    real(dp)::center(14),faces(14,2*ndim)
    integer::idx(14),l,i,j,k,d,e,status
    if(trim(configured_radioactive_model)=='none')return
    idx=[imetal,(ichem+e-1,e=1,11),iradioactive,iradioactive+1]
    ! Trace routines initialize the interior stencil, not the outermost halo.
    do k=min(1,ku1+1),max(1,ku2-1)
       do j=min(1,ju1+1),max(1,ju2-1)
          do i=min(1,iu1+1),max(1,iu2-1)
             do l=1,ngrid
                center=qin(l,i,j,k,idx)
                do d=1,ndim
                   faces(:,2*d-1)=qm(l,i,j,k,idx,d)
                   faces(:,2*d)=qp(l,i,j,k,idx,d)
                enddo
                call radioactive_limit_states(center,faces,status)
                if(status/=0)then
                   write(*,*)'ERROR: invalid radioactive reconstruction center/states'
                   call clean_stop;return
                endif
                do d=1,ndim
                   qm(l,i,j,k,idx,d)=faces(:,2*d-1)
                   qp(l,i,j,k,idx,d)=faces(:,2*d)
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine

  subroutine radioactive_limit_children(u1,u2,nn)
    integer,intent(in)::nn
    real(dp),intent(in)::u1(nvector,0:twondim,nvar)
    real(dp),intent(inout)::u2(nvector,twotondim,nvar)
    real(dp)::children(14,twotondim)
    integer::idx(14),l,e,status
    if(trim(configured_radioactive_model)=='none')return
    idx=[imetal,(ichem+e-1,e=1,11),iradioactive,iradioactive+1]
    do l=1,nn
       children=transpose(u2(l,:,idx))
       call radioactive_limit_states(u1(l,0,idx),children,status)
       if(status/=0)then
          write(*,*)'ERROR: invalid radioactive prolongation center/states'
          call clean_stop;return
       endif
       u2(l,:,idx)=transpose(children)
    enddo
  end subroutine
end module stellar_radioactive_runtime
