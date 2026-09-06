! Compact numeric-fixture driver for the actual native module, no RAMSES stubs.
program snrt_dust_ir_smoke
  use snrt_dust_ir
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  type(dust_ir_table) :: table
  type(dust_ir_diagnostics) :: diag, before_diag
  real(dust_dp), allocatable :: freq(:), fw(:), sigma(:), nodes(:), direction(:,:), weight(:)
  real(dust_dp), allocatable :: density(:), primary(:), energy(:,:,:), photons(:,:), temperature(:)
  real(dust_dp), allocatable :: before(:,:,:), before_p(:,:), before_t(:), integral(:,:)
  integer, allocatable :: neighbor(:,:)
  real(dust_dp) :: dx,dt,c,cmb,tol,escaped,absorbed,injected,stored,balance
  integer :: ng,nt,nd,nc,steps,i,g,d,it,unit,ierr
  character(len=2048) :: filename
  call get_command_argument(1,filename)
  open(newunit=unit,file=trim(filename),status='old',action='read')
  read(unit,*) ng,nt,nd,nc,steps
  read(unit,*) dx,dt,c,cmb,tol
  allocate(freq(ng),fw(ng),sigma(ng),nodes(nt),direction(3,nd),weight(nd),neighbor(6,nc))
  allocate(density(nc),primary(nc),energy(ng,nd,nc),photons(ng,nc),temperature(nc),integral(ng,nc))
  read(unit,*) freq
  read(unit,*) fw
  read(unit,*) sigma
  read(unit,*) nodes
  read(unit,*) direction
  read(unit,*) weight
  read(unit,*) neighbor
  read(unit,*) density
  read(unit,*) primary
  close(unit)
  call snrt_dust_ir_initialize(table,freq,fw,sigma,nodes,cmb,ierr)
  if (ierr/=dust_ok) stop 1
  energy=0; photons=0; temperature=0; escaped=0; absorbed=0; injected=0
  do it=1,steps
     call snrt_dust_ir_advance(table,direction,weight,neighbor,dx,dt,c,density,primary, &
          energy,temperature,photons,diag,ierr,tol,128)
     if (ierr/=dust_ok) then
        write(*,*) 'NATIVE_DUST_STEP_ERROR',ierr,it
        stop 2
     end if
     escaped=escaped+diag%escaped_erg; absorbed=absorbed+diag%absorbed_erg
     injected=injected+diag%primary_erg
  end do
  before=energy; before_p=photons; before_t=temperature; before_diag=diag
  ! Force a late nonlinear failure, and separately a CFL admission failure.
  call snrt_dust_ir_advance(table,direction,weight,neighbor,dx,dt,c,density,primary, &
       energy,temperature,photons,diag,ierr,tol,1)
  if (ierr/=dust_err_convergence) stop 3
  call assert_unchanged()
  call snrt_dust_ir_advance(table,direction,weight,neighbor,dx,100*dt,c,density,primary, &
       energy,temperature,photons,diag,ierr,tol,128)
  if (ierr/=dust_err_cfl) stop 4
  call assert_unchanged()
  ! A reciprocal-link violation is caught without a partial state commit.
  i=neighbor(2,1)
  if (i<=0) stop 5
  neighbor(1,i)=0
  call snrt_dust_ir_advance(table,direction,weight,neighbor,dx,dt,c,density,primary, &
       energy,temperature,photons,diag,ierr,tol,128)
  if (ierr/=dust_err_state) stop 6
  call assert_unchanged()
  neighbor(1,i)=1
  ! Save physical result before zero/weak checks, which have their own state.
  integral=0
  do d=1,nd
     integral=integral+energy(:,d,:)*weight(d)
  end do
  stored=sum(integral)*dx**3
  balance=abs(stored+escaped-injected)/injected
  if (balance>tol) stop 7
  write(*,'(a)') 'NATIVE_DUST_IR_OK'
  write(*,'(es26.16e3)') stored,escaped,absorbed,balance
  write(*,'(es26.16e3)') temperature
  write(*,'(es26.16e3)') integral
  write(*,'(es26.16e3)') photons
  energy=0; photons=0; temperature=0; primary=0
  call snrt_dust_ir_advance(table,direction,weight,neighbor,dx,dt,c,density,primary, &
       energy,temperature,photons,diag,ierr,tol,128)
  if (ierr/=dust_ok .or. any(energy/=0) .or. any(photons/=0)) stop 8
  density=1; primary=1d-42
  call snrt_dust_ir_advance(table,direction,weight,neighbor,dx,dt,c,density,primary, &
       energy,temperature,photons,diag,ierr,tol,128)
  if (ierr/=dust_ok .or. maxval(energy)<=0 .or. maxval(photons)<=0) stop 9
contains
  subroutine assert_unchanged()
    if (any(transfer(energy,[0_int64],size(energy))/=transfer(before,[0_int64],size(before)))) stop 11
    if (any(transfer(photons,[0_int64],size(photons))/=transfer(before_p,[0_int64],size(before_p)))) stop 12
    if (any(transfer(temperature,[0_int64],size(temperature))/=transfer(before_t,[0_int64],size(before_t)))) stop 13
    ! Compare fields, not derived-type padding bytes.
    if (any(transfer([diag%escaped_erg,diag%absorbed_erg,diag%primary_erg, &
         diag%balance_relative,diag%local_relative],[0_int64],5)/= &
         transfer([before_diag%escaped_erg,before_diag%absorbed_erg,before_diag%primary_erg, &
         before_diag%balance_relative,before_diag%local_relative],[0_int64],5))) stop 14
    if (diag%iterations/=before_diag%iterations) stop 15
  end subroutine
end program
