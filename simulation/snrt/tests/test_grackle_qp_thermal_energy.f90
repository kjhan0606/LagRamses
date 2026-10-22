program test_grackle_qp_thermal_energy
  use amr_parameters, only: dp,ndim
  use grackle_neq_runtime, only: grackle_neq_thermal_energy
  implicit none
  integer,parameter::qp=selected_real_kind(33,4931)
  real(dp)::rho(4),momentum(3,4),energy(4),thermal_dp,thermal_grackle,thermal_scale
  real(qp)::thermal_exact
  integer::i,near_margin_count
  if(ndim/=3)error stop 'Grackle thermal helper regression requires 3D'
  rho=[6.0224011572532203d-2,3.1771656350967686d-1, &
       2.3974987490560998d-2,2.0619790522295375d-1]
  momentum(:,1)=[2.9694834334669629d-4,-2.9898351544043049d-5,-5.9883391717737419d-5]
  momentum(:,2)=[-8.4576341883053094d-4,-1.1312328652722783d-3,-1.7432285066087646d-3]
  momentum(:,3)=[-2.4274468670122109d-5,-1.0518456302331167d-5,9.2086595984272390d-5]
  momentum(:,4)=[7.8543022308416848d-4,5.1513466580278204d-4,-9.0783119697291498d-4]
  energy=[7.6927996180911619d-7,7.9219180875411757d-6, &
          1.9144595833022912d-7,4.1378253566591302d-6]
  near_margin_count=0
  do i=1,4
     thermal_dp=energy(i)-0.5d0*sum(momentum(:,i)**2)/rho(i)
     thermal_scale=spacing(max(abs(energy(i)),abs(0.5d0*sum(momentum(:,i)**2)/rho(i))))
     thermal_exact=real(energy(i),qp)-sum(real(momentum(:,i),qp)**2)/(2._qp*real(rho(i),qp))
     thermal_grackle=grackle_neq_thermal_energy(rho(i),momentum(:,i),energy(i),0d0)
     if(thermal_exact<=0._qp.or.thermal_grackle<=0d0) &
          error stop 'positive conserved-state residual did not reach Grackle as positive energy'
     if(thermal_exact>=8._qp*real(thermal_scale,qp).or.abs(thermal_dp)>=8d0*thermal_scale) &
          error stop 'regression fixture is outside the 8-ULP near-cancellation band'
     if(thermal_grackle/=real(thermal_exact,dp)) &
          error stop 'Grackle thermal helper disagrees with extended-precision residual'
     near_margin_count=near_margin_count+1
  enddo
  if(near_margin_count/=4)error stop 'near-cancellation fixture count mismatch'
  write(*,'(A)')'GRACKLE_QP_THERMAL_COOLING_INPUT_PASS cases=4 dp_near_margin_cases=4'
  write(*,'(A)')'GRACKLE_QP_THERMAL_CONSERVED_INPUT_PASS'
end program test_grackle_qp_thermal_energy
