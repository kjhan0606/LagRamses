program test_exact
  use eunha_cooling_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  real(8), parameter :: kb=1.3806504d-16,mh=1.6726219d-24,c=kb/((5d0/3-1)*mh)
  real(8) :: temp(4),raw(4,2,2),mu(4,2,2),rates(4),rho
  real(8) :: nh(1),z(1),t(1),change(1),whole,half,expected,dt
  integer :: mode,i,j,k
  character(1024) :: real_table
  temp=[100d0,1000d0,10000d0,100000d0]
  mu=1d0;nh=sqrt(10d0);z=sqrt(0.1d0)
  do mode=1,6
     mu=1d0
     if(mode==5) then
        do k=1,4
           mu(k,:,:)=1.2d0-.05d0*(k-1)
        enddo
     endif
     do i=1,2
        do j=1,2
           rho=10d0**(i-1)*mh/(.76d0*(1-.02d0*10d0**(j-2)))
           do k=1,4
              select case(mode)
              case(1); rates(k)=c*10d0**(i+j-2) ! different corner normalizations
              case(2); rates(k)=-c*10d0**(i+j-2)
              case(3); rates(k)=c*(temp(k)-3000d0)*1d-3
              case(4); rates(k)=c*temp(k)*1d-3
              case(5); rates(k)=c*temp(k)/mu(k,j,i)*1d-3
              case(6); rates(k)=c*(temp(k)-3000d0*(i+j-1))*1d-3
              end select
              raw(k,j,i)=-rates(k)*rho
           enddo
        enddo
     enddo
     open(20,file='synthetic.bin',access='stream',form='unformatted',status='replace')
     write(20) 2,2,4
     write(20) 0d0,0d0,1d0,-1d0,0d0,2d0,5d0
     write(20) log10(temp),raw,mu
     close(20)
     call eunha_load_table('synthetic.bin')
     t=5000d0;dt=100d0
     select case(mode)
     case(1);expected=t(1)-30.25d0*dt
     case(2);expected=t(1)+30.25d0*dt
     case(3);expected=3000d0+(t(1)-3000d0)*exp(-1d-3*dt)
     case(4);expected=t(1)*exp(-1d-3*dt)
     case(5);expected=t(1)*exp(-1d-3*dt)
     case(6);expected=6000d0+(t(1)-6000d0)*exp(-1d-3*dt)
     end select
     call eunha_solve(nh,t,z,dt,change,1)
     whole=t(1)+change(1)
     if(abs(whole-expected)>1d-10*abs(expected)) error stop 'analytic mismatch'
     call eunha_solve(nh,t,z,dt/2,change,1)
     t=t+change
     call eunha_solve(nh,t,z,dt/2,change,1)
     half=t(1)+change(1)
     if(abs(half-whole)>1d-10*abs(whole)) error stop 'split mismatch'
     if(mode==3) then
        t=5000d0
        call eunha_solve(nh,t,z,1d10,change,1)
        if(abs(t(1)+change(1)-3000d0)>1d-8) error stop 'cooling equilibrium'
        t=2000d0
        call eunha_solve(nh,t,z,1d10,change,1)
        if(abs(t(1)+change(1)-3000d0)>1d-8) error stop 'heating equilibrium'
     endif
     t=5000d0
     call eunha_solve(nh,t,z,0d0,change,1)
     if(change(1)/=0) error stop 'zero timestep'
     if(mode<=2) then
        call eunha_solve(nh,t,z,1d10,change,1)
        if(mode==1.and.abs(t(1)+change(1)-100d0)>1d-8) error stop 'lower boundary'
        if(mode==2.and.abs(t(1)+change(1)-100000d0)>1d-8) error stop 'upper boundary'
     endif
  enddo
  print *, 'PASS: corner interpolation, analytic cooling/heating, equilibrium, split, dt=0'
  call get_command_argument(1,real_table)
  if(len_trim(real_table)>0) then
     call eunha_load_multi_z(trim(real_table))
     do i=0,2
        call eunha_interp_redshift(3d0*i)
        do j=0,3
           nh=10d0**(j-3);z=0.1d0
           do k=0,3
              t=10d0**(k+2)
              call eunha_solve(nh,t,z,1d10,change,1)
              whole=t(1)+change(1)
              if(.not.ieee_is_finite(whole).or.whole<=0) error stop 'real table state'
              call eunha_solve(nh,t,z,5d9,change,1)
              t=t+change
              call eunha_solve(nh,t,z,5d9,change,1)
              if(abs(t(1)+change(1)-whole)>1d-10*abs(whole)) error stop 'real table split'
           enddo
        enddo
     enddo
     print *, 'PASS: existing multi-z table, 48 states, redshift interpolation and split'
  endif
end program
