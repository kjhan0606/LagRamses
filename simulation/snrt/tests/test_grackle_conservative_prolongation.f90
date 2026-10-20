program check_prolongation
  use amr_parameters
  use hydro_parameters
  implicit none
  real(dp)::parent(nvar),u1(nvector,0:twondim,nvar)
  real(dp)::u2(nvector,twotondim,nvar),unlimited(nvector,twotondim,nvar)
  real(dp)::internal,momentum,thermal,min_thermal,average(nvar)
  integer::test_case,j,ind
  character(40)::mode
  call get_command_argument(1,mode)
  if(trim(mode)=='--reject-thermal-parent'.or.trim(mode)=='--reject-thermal-margin')then
     call reject_parent_from_worker(trim(mode)=='--reject-thermal-margin')
     error stop 'invalid parent was accepted'
  endif
  if(trim(mode)=='--reproduce-mass-closure')then
     call reproduce_mass_closure(.true.)
     stop
  endif
  if(trim(mode)=='--verify-mass-closure')then
     call reproduce_mass_closure(.false.)
     stop
  endif
  if(trim(mode)=='--reproduce-face-closure'.or.trim(mode)=='--verify-face-closure')then
     call check_face_closure(trim(mode)=='--reproduce-face-closure')
     stop
  endif
  if(trim(mode)=='--verify-inherited-face-advection')then
     call check_face_closure(.false.,9.95d-9)
     call check_face_closure(.false.,-9.95d-9)
     call check_face_closure(.false.,1d-10)
     write(*,*)'GRACKLE_NATIVE_INHERITED_FACE_ADVECTION_PASS'
     stop
  endif
  if(len_trim(mode)>0)error stop 'unknown native regression mode'
  imetal=6;ichem=7;igrackle=19
  interpol_var=0;interpol_type=1
  do test_case=1,4
     internal=1d-6
     if(test_case==2)internal=1d-2
     if(test_case==4)internal=2d-13
     parent=1d-3;parent(1)=1d0;parent(2:ndim+1)=0d0
     parent(2)=1d0;parent(ndim+2)=0.5d0+internal
     u1=0d0
     do j=0,twondim
        u1(1,j,:)=parent
     enddo
     if(test_case/=3)then
        do j=1,2
           momentum=1d0+merge(-0.01d0,0.01d0,j==1)
           u1(1,j,2)=momentum
           u1(1,j,ndim+2)=0.5d0*momentum**2+internal
        enddo
     endif
     cooling_method='original'
     call interpol_hydro(u1,unlimited,1)
     if(test_case==1.or.test_case==4)then
        min_thermal=huge(1d0)
        do ind=1,twotondim
           thermal=unlimited(1,ind,ndim+2)-0.5d0*sum(unlimited(1,ind,2:ndim+1)**2)/unlimited(1,ind,1)
           min_thermal=min(min_thermal,thermal)
        enddo
        if(min_thermal>=0d0)error stop 'fixture did not reproduce negative thermal energy'
     endif
     cooling_method='grackle_neq'
     call interpol_hydro(u1,u2,1)
     average=sum(u2(1,:,:),dim=1)/real(twotondim,dp)
     if(maxval(abs(average-parent))>128d0*spacing(1d0))error stop 'parent average not conserved'
     do ind=1,twotondim
        thermal=u2(1,ind,ndim+2)-0.5d0*sum(u2(1,ind,2:ndim+1)**2)/u2(1,ind,1)
        if(thermal<=0d0)error stop 'nonpositive limited child'
     enddo
     if(test_case==2.or.test_case==3)then
        if(any(u2(1,:,:)/=unlimited(1,:,:)))error stop 'admissible stencil changed'
     endif
     write(*,*)'PROLONGATION_PASS case=',test_case,' max_budget_error=',maxval(abs(average-parent))
  enddo
contains
  subroutine reject_parent_from_worker(positive_below_margin)
    use omp_lib, only: omp_get_thread_num,omp_get_num_threads
    logical,intent(in)::positive_below_margin
    real(dp)::stencil(nvector,0:twondim,nvar),children(nvector,twotondim,nvar)
    integer::neighbor
    imetal=6;ichem=7;igrackle=19
    interpol_var=0;interpol_type=1;cooling_method='grackle_neq'
!$omp parallel num_threads(2) private(stencil,children,neighbor)
    if(omp_get_num_threads()/=2)error stop 'negative test requires two OpenMP threads'
    if(omp_get_thread_num()==1)then
       stencil=0d0
       do neighbor=0,twondim
          stencil(1,neighbor,1)=1d0
          stencil(1,neighbor,2)=1d0
          stencil(1,neighbor,ndim+2)=.4d0
          if(positive_below_margin)stencil(1,neighbor,ndim+2)=.5d0+64d0*spacing(.5d0)
       enddo
       write(*,*)'GRACKLE_THERMAL_REJECT_WORKER=',omp_get_thread_num()
       call interpol_hydro(stencil,children,1)
    endif
!$omp end parallel
  end subroutine reject_parent_from_worker

  subroutine check_face_closure(expect_defect,inherited_residual)
    logical,intent(in)::expect_defect
    real(dp),intent(in),optional::inherited_residual
    real(dp)::states(nvar,4),newstate(nvar),before(nvar),after(nvar),fractions(3,4)
    real(dp)::density(4),ion(4),mol(4),heion(4),velocity,dt,residual,max_flux_error,max_cell_error
    real(dp)::carrier_sum
    real(dp)::u(nvector,iu1:iu2,ju1:ju2,ku1:ku2,nvar)
    real(dp)::grav(nvector,iu1:iu2,ju1:ju2,ku1:ku2,ndim)
    real(dp)::outer(nvector,2,0:3,0:3,nvar,ndim)
    real(dp)::flux(nvector,if1:if2,jf1:jf2,kf1:kf2,nvar,ndim)
    real(dp)::tmp(nvector,if1:if2,jf1:jf2,kf1:kf2,2,ndim)
    integer::case_id,tile,ii,jj,kk,gid,direction,ends(3),zero_mass_faces,gradient_axis,position(3)
    if(ndim/=3.or.nvar<27)error stop 'face regression needs the 3D 27-field Grackle build'
    imetal=6;ichem=7;igrackle=19
    cooling_method='grackle_neq';scheme='muscl';slope_type=1;gamma=5d0/3d0;difmag=0d0
    carrier_sum=1d0
    if(present(inherited_residual))carrier_sum=1d0+inherited_residual
    fractions(:,1)=[.70d0,.20d0,.10d0];fractions(:,2)=[.65d0,.23d0,.12d0]
    fractions(:,3)=[.55d0,.24d0,.21d0];fractions(:,4)=[.60d0,.28d0,.12d0]
    ion=[.1d0,.3d0,.6d0,.2d0];mol=[.01d0,.02d0,.03d0,.005d0]
    heion=[.3d0,.1d0,.2d0,.05d0];dt=.01d0
    do gradient_axis=1,ndim
    density=[1d0,2d0,4d0,2d0]
    do case_id=1,3
       riemann='hllc'
       if(case_id==3)then
          riemann='llf';density=1d0
       endif
       states=0d0
       do gid=1,4
          velocity=merge(1d0,-1d0,case_id==1)
          if(case_id==3)velocity=merge(1d0,-1d0,gid<=2)
          states(1,gid)=density(gid);states(gradient_axis+1,gid)=density(gid)*velocity
          states(ndim+2,gid)=density(gid)*(10d0+.5d0*velocity**2)
          states(ichem:ichem+1,gid)=density(gid)*fractions(1:2,gid)
          states(imetal,gid)=density(gid)*fractions(3,gid)
          states(ichem+2,gid)=.3d0*states(imetal,gid)
          states(ichem+10,gid)=.2d0*states(imetal,gid)
          states(igrackle,gid)=states(ichem,gid)*(1d0-ion(gid)-mol(gid))
          states(igrackle+1,gid)=states(ichem,gid)*ion(gid)
          states(igrackle+6,gid)=states(ichem,gid)*mol(gid)
          states(igrackle+3,gid)=states(ichem+1,gid)*(1d0-heion(gid))
          states(igrackle+4,gid)=states(ichem+1,gid)*heion(gid)
          states(igrackle+8,gid)=states(igrackle+1,gid)+states(igrackle+4,gid)/4d0
          states(imetal:nvar,gid)=states(imetal:nvar,gid)*carrier_sum
       enddo
       before=4d0*sum(states,dim=2);after=0d0
       max_flux_error=0d0;max_cell_error=0d0;zero_mass_faces=0
       do tile=0,1
          u=0d0;grav=0d0;outer=0d0;flux=0d0;tmp=0d0
          do kk=ku1,ku2
             do jj=ju1,ju2
                do ii=iu1,iu2
                   position=[ii,jj,kk]
                   gid=modulo(position(gradient_axis)-1+2*tile,4)+1
                   u(1,ii,jj,kk,:)=states(:,gid)
                enddo
             enddo
          enddo
          call unsplit(u,grav,flux,tmp,1d0,1d0,1d0,dt,1,outer)
          do direction=1,ndim
             ends=2;ends(direction)=3
             do kk=1,ends(3)
                do jj=1,ends(2)
                   do ii=1,ends(1)
                      residual=flux(1,ii,jj,kk,ichem,direction)+flux(1,ii,jj,kk,ichem+1,direction)+ &
                           flux(1,ii,jj,kk,imetal,direction)-carrier_sum*flux(1,ii,jj,kk,1,direction)
                      max_flux_error=max(max_flux_error,abs(residual))
                      if(case_id==3.and.direction==gradient_axis.and.abs(flux(1,ii,jj,kk,1,direction))<1d-15) &
                           zero_mass_faces=zero_mass_faces+1
                   enddo
                enddo
             enddo
          enddo
          do kk=1,2
             do jj=1,2
                do ii=1,2
                   position=[ii,jj,kk]
                   gid=modulo(position(gradient_axis)-1+2*tile,4)+1
                   newstate=states(:,gid)+flux(1,ii,jj,kk,:,1)-flux(1,ii+1,jj,kk,:,1)+ &
                        flux(1,ii,jj,kk,:,2)-flux(1,ii,jj+1,kk,:,2)+ &
                        flux(1,ii,jj,kk,:,3)-flux(1,ii,jj,kk+1,:,3)
                   after=after+newstate
                   residual=newstate(ichem)+newstate(ichem+1)+newstate(imetal)-carrier_sum*newstate(1)
                   max_cell_error=max(max_cell_error,abs(residual))
                   if(present(inherited_residual))then
                      residual=(newstate(ichem)+newstate(ichem+1)+newstate(imetal)-newstate(1))/newstate(1)
                      if(abs(residual)>1d-8)error stop 'transport amplified admitted relative residual past gate'
                   endif
                   if(newstate(1)<=0d0.or.any(newstate(ichem:ichem+10)<0d0).or. &
                        any(newstate(igrackle:igrackle+8)<0d0).or.newstate(imetal)<0d0) &
                        error stop 'negative periodic update'
                   if(.not.expect_defect)then
                      if(abs(sum(newstate(igrackle+[0,1,2,6,7]))-newstate(ichem))> &
                           256d0*spacing(maxval(density)).or. &
                           abs(sum(newstate(igrackle+3:igrackle+5))-newstate(ichem+1))> &
                           256d0*spacing(maxval(density)))error stop 'periodic species budget is not closed'
                      residual=newstate(igrackle+1)-newstate(igrackle+2)+newstate(igrackle+7)/2d0+ &
                           newstate(igrackle+4)/4d0+newstate(igrackle+5)/2d0-newstate(igrackle+8)
                      if(abs(residual)>256d0*spacing(maxval(density)))error stop 'periodic charge is not closed'
                   endif
                enddo
             enddo
          enddo
       enddo
       if(maxval(abs(after-before))>256d0*spacing(maxval(abs(before)))) &
            error stop 'periodic transport changed an independent global inventory'
       if(case_id==3.and.zero_mass_faces==0)error stop 'zero-mass-flux LLF case was not exercised'
       write(*,'(A,2I2,3ES25.16,I4)')'NATIVE_FACE axis/case/flux/cell/global/zero_mass_faces ',gradient_axis,case_id, &
            max_flux_error,max_cell_error,maxval(abs(after-before)),zero_mass_faces
       if(expect_defect)then
          if(case_id<=2.and.max_flux_error<=1d-12)error stop 'expected face mass defect was not reproduced'
       else
          if(max_flux_error>256d0*spacing(1d0).or.max_cell_error>256d0*spacing(maxval(density))) &
               error stop 'face or periodic-cell mass identity is not closed'
       endif
    enddo
    enddo
    if(expect_defect)then
       write(*,*)'GRACKLE_NATIVE_FACE_CLOSURE_REPRODUCED: not a repair pass'
    else
       write(*,*)'GRACKLE_NATIVE_FACE_CLOSURE_PASS'
    endif
  end subroutine check_face_closure

  subroutine reproduce_mass_closure(expect_defect)
    logical,intent(in)::expect_defect
    real(dp)::stencil(nvector,0:twondim,nvar),child(nvector,twotondim,nvar)
    real(dp)::central(nvar),mean_state(nvar),residual,max_residual
    integer::neighbor,octant
    if(nvar<27)error stop 'mass-closure reproduction needs the 27-field Grackle build'
    imetal=6;ichem=7;igrackle=19
    cooling_method='grackle_neq';interpol_var=0;interpol_type=1
    central=0d0
    central(1)=2d0;central(ndim+2)=10d0
    central(imetal)=.2d0;central(ichem)=1.2d0;central(ichem+1)=.6d0
    central(igrackle)=central(ichem);central(igrackle+3)=central(ichem+1)
    stencil=0d0
    do neighbor=0,twondim
       stencil(1,neighbor,:)=central
    enddo
    stencil(1,1,1)=1d0
    stencil(1,1,ichem)=.7d0;stencil(1,1,ichem+1)=.2d0;stencil(1,1,imetal)=.1d0
    stencil(1,2,1)=4d0
    stencil(1,2,ichem)=2.8d0;stencil(1,2,ichem+1)=.8d0;stencil(1,2,imetal)=.4d0
    do neighbor=1,2
       stencil(1,neighbor,igrackle)=stencil(1,neighbor,ichem)
       stencil(1,neighbor,igrackle+3)=stencil(1,neighbor,ichem+1)
    enddo
    do neighbor=0,twondim
       residual=stencil(1,neighbor,ichem)+stencil(1,neighbor,ichem+1)+ &
            stencil(1,neighbor,imetal)-stencil(1,neighbor,1)
       if(abs(residual)>64d0*spacing(stencil(1,neighbor,1)))error stop 'fixture parent is not closed'
    enddo
    call interpol_hydro(stencil,child,1)
    mean_state=sum(child(1,:,:),dim=1)/real(twotondim,dp)
    if(maxval(abs(mean_state-central))>128d0*spacing(maxval(abs(central)))) &
         error stop 'prolongation lost an independent parent integral'
    max_residual=0d0
    do octant=1,twotondim
       residual=child(1,octant,ichem)+child(1,octant,ichem+1)+ &
            child(1,octant,imetal)-child(1,octant,1)
       max_residual=max(max_residual,abs(residual))
       if(.not.expect_defect)then
          if(any(child(1,octant,ichem:ichem+10)<0d0).or.child(1,octant,imetal)<0d0.or. &
               any(child(1,octant,igrackle:igrackle+8)<0d0))error stop 'negative child composition'
          if(maxval(abs(child(1,octant,ichem:ichem+10)/child(1,octant,1)- &
               central(ichem:ichem+10)/central(1)))>64d0*spacing(1d0)) &
               error stop 'child abundance differs from parent'
          if(abs(child(1,octant,igrackle)-child(1,octant,ichem))>64d0*spacing(central(1)).or. &
               abs(child(1,octant,igrackle+3)-child(1,octant,ichem+1))>64d0*spacing(central(1))) &
               error stop 'child species budget differs from element budget'
       endif
       write(*,'(A,I2,5ES25.16)')'NATIVE_CHILD rho/H/He/Z/residual ',octant, &
            child(1,octant,1),child(1,octant,ichem),child(1,octant,ichem+1), &
            child(1,octant,imetal),residual
    enddo
    if(expect_defect)then
       if(max_residual<=1d-8*maxval(child(1,:,1)))error stop 'expected mass-closure defect was not reproduced'
       write(*,'(A,ES25.16)')'GRACKLE_NATIVE_MASS_CLOSURE_REPRODUCED max_absolute_residual=',max_residual
       ! This expected-defect mode is a source reproduction, not a repair pass.
    else
       if(max_residual>64d0*spacing(maxval(child(1,:,1))))error stop 'child mass identity is not closed'
       write(*,'(A,ES25.16)')'GRACKLE_NATIVE_MASS_CLOSURE_PASS max_absolute_residual=',max_residual
       ! Inherited residuals must survive conservative prolongation, not be
       ! hidden by normalization. An inadmissible restart still needs diagnosis.
       stencil(1,0,ichem)=stencil(1,0,ichem)-1d-6
       stencil(1,0,igrackle)=stencil(1,0,igrackle)-1d-6
       central=stencil(1,0,:)
       call interpol_hydro(stencil,child,1)
       mean_state=sum(child(1,:,:),dim=1)/real(twotondim,dp)
       if(maxval(abs(mean_state-central))>128d0*spacing(maxval(abs(central)))) &
            error stop 'inconsistent parent inventory was altered'
       do octant=1,twotondim
          residual=(child(1,octant,ichem)+child(1,octant,ichem+1)+ &
               child(1,octant,imetal)-child(1,octant,1))/child(1,octant,1)
          if(abs(residual+5d-7)>128d0*spacing(1d0))error stop 'inherited closure defect was hidden'
       enddo
       write(*,*)'GRACKLE_NATIVE_INHERITED_RESIDUAL_PRESERVED'
    endif
  end subroutine reproduce_mass_closure
end program
