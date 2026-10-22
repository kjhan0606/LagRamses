!###########################################################
!###########################################################
! Eunha Exact Cooling Integrator (Fortran port)
!
! Method: analytic segment integration of a piecewise-linear energy/rate curve
!         Zhu, Smith & Hernquist (2017), MNRAS 470, 1017
!
! Table: Grackle equilibrium net-rate binary
!        (nH, Z_solar, T) with mean molecular weight
!
! Units: RAMSES T/mu [K], nH [cm^-3], dt [s]
!###########################################################
!###########################################################
module eunha_cooling_mod
  use amr_parameters, only: dp
  implicit none
  private

  ! Physical constants (CGS)
  real(dp), parameter :: k_B       = 1.3806504d-16   ! Boltzmann [erg/K]
  real(dp), parameter :: m_H       = 1.6726219d-24   ! proton mass [g]
  real(dp), parameter :: GAMMA_AD  = 5.0d0/3.0d0     ! adiabatic index
  real(dp), parameter :: X_H       = 0.76d0          ! primordial H fraction
  real(dp), parameter :: Y_He      = 0.24d0          ! primordial He fraction
  real(dp), parameter :: Z_SOLAR   = 0.02d0          ! solar metallicity
  real(dp), parameter :: A_metal   = 16.0d0          ! mean metal atomic weight
  real(dp), parameter :: T_ion_lo  = 6000.0d0        ! H ionization start
  real(dp), parameter :: T_ion_hi  = 20000.0d0       ! fully ionized
  real(dp), parameter :: PI_VAL    = 3.14159265358979323846d0
  real(dp), parameter :: EPSILON   = 1.0d-30
  real(dp), parameter :: REL_TOL   = 1.0d-8

  ! Table arrays — (N_temp, M_metal, L_dens) for column-major efficiency
  real(dp), allocatable :: tab_density(:)        ! (L) nH [cm^-3]
  real(dp), allocatable :: tab_metallicity(:)    ! (M) Z_solar
  real(dp), allocatable :: tab_temperature(:)    ! (N) T [K]
  real(dp), allocatable :: tab_net_rate(:,:,:)   ! (N,M,L) erg/g/s
  real(dp), allocatable :: tab_mmw(:,:,:)        ! (N,M,L) mean mol weight
  integer  :: tab_L, tab_M, tab_N
  logical  :: tab_loaded = .false.

  ! Multi-z table storage
  integer  :: tab_Nz = 0                            ! number of z snapshots (0 = single-z)
  real(dp), allocatable :: tab_z_values(:)           ! (Nz) redshift array
  real(dp), allocatable :: tab_net_rate_z(:,:,:,:)   ! (N, M, L, Nz) all z snapshots [erg/g/s]
  real(dp), allocatable :: tab_mmw_z(:,:,:,:)        ! (N, M, L, Nz) all z snapshots
  real(dp) :: tab_current_z = -1.0d0                 ! last interpolated z (cache)

  ! Comparison accumulators (OMP atomic safe)
  real(dp) :: acc_sum_rdiff  = 0.0d0
  real(dp) :: acc_max_rdiff  = 0.0d0
  real(dp) :: acc_sum_dT_orig = 0.0d0
  real(dp) :: acc_sum_dT_exact = 0.0d0
  integer(8) :: acc_ncell = 0
  integer(8) :: acc_nsign = 0

  public :: eunha_load_table, eunha_load_multi_z, eunha_interp_redshift
  public :: eunha_solve, eunha_report, eunha_reset_stats
  public :: acc_sum_rdiff, acc_max_rdiff, acc_sum_dT_orig, acc_sum_dT_exact
  public :: acc_ncell, acc_nsign

contains

!###########################################################
! Load Grackle binary table
!
! Binary format (field-by-field, no struct padding):
!   3 ints:  L, M, N  (density, metallicity, temperature)
!   7 doubles: redshift, log10(dens_min), log10(dens_max),
!              log10(Z_min), log10(Z_max),
!              log10(T_min), log10(T_max)
!   N doubles: log10(temperature) grid
!   L*M*N doubles: net_rate [erg/cm^3/s]
!   L*M*N doubles: mmw [dimensionless]
!###########################################################
subroutine eunha_load_table(filename)
  implicit none
  character(len=*), intent(in) :: filename
  integer :: iunit, ierr, i, j, k
  integer :: nd, nm, nt
  real(dp) :: z_val, log_dmin, log_dmax, log_zmin, log_zmax, log_tmin, log_tmax
  real(dp) :: dstep, zstep
  real(dp), allocatable :: log_T(:)
  real(dp), allocatable :: raw_rate(:,:,:), raw_mmw(:,:,:)
  real(dp) :: nH_val, Z_metal, rho
  real(dp) :: mf_X, mf_Z, f_scale

  iunit = 42
  open(unit=iunit, file=filename, status='old', access='stream', &
       form='unformatted', iostat=ierr)
  if(ierr /= 0) then
     write(*,*) 'ERROR: eunha_load_table: cannot open ', trim(filename)
     return
  end if

  ! Read dimensions
  read(iunit) nd, nm, nt
  tab_L = nd
  tab_M = nm
  tab_N = nt

  ! Read header doubles
  read(iunit) z_val, log_dmin, log_dmax, log_zmin, log_zmax, log_tmin, log_tmax

  ! Read temperature grid (log10 values)
  allocate(log_T(nt))
  read(iunit) log_T(1:nt)

  ! Read net_rate and mmw — raw C layout: (density, metallicity, temperature)
  ! In file: fastest index = temperature, matching Fortran (N,M,L)
  allocate(raw_rate(nt, nm, nd))
  allocate(raw_mmw(nt, nm, nd))
  read(iunit) raw_rate
  read(iunit) raw_mmw
  close(iunit)

  ! Allocate table arrays
  if(allocated(tab_density))     deallocate(tab_density)
  if(allocated(tab_metallicity)) deallocate(tab_metallicity)
  if(allocated(tab_temperature)) deallocate(tab_temperature)
  if(allocated(tab_net_rate))    deallocate(tab_net_rate)
  if(allocated(tab_mmw))        deallocate(tab_mmw)

  allocate(tab_density(nd))
  allocate(tab_metallicity(nm))
  allocate(tab_temperature(nt))
  allocate(tab_net_rate(nt, nm, nd))
  allocate(tab_mmw(nt, nm, nd))

  ! Build density grid (log-spaced from log10 values)
  dstep = (log_dmax - log_dmin) / dble(nd - 1)
  do i = 1, nd
     tab_density(i) = 10.0d0 ** (log_dmin + dble(i-1) * dstep)
  end do

  ! Build metallicity grid (log-spaced, in Z_solar units)
  zstep = (log_zmax - log_zmin) / dble(nm - 1)
  do j = 1, nm
     tab_metallicity(j) = 10.0d0 ** (log_zmin + dble(j-1) * zstep)
  end do

  ! Temperature grid (convert from log10)
  do k = 1, nt
     tab_temperature(k) = 10.0d0 ** log_T(k)
  end do
  deallocate(log_T)

  ! Copy mmw table (dimensionless, no conversion)
  tab_mmw(:,:,:) = raw_mmw(:,:,:)
  deallocate(raw_mmw)

  ! Convert net_rate: erg/cm^3/s -> erg/g/s
  ! rho = nH * m_H / X_H  (with metallicity-dependent X_H)
  ! Negate sign: Grackle convention (negative=cooling) -> RAMSES (positive=cooling)
  do i = 1, nd
     nH_val = tab_density(i)
     do j = 1, nm
        Z_metal = tab_metallicity(j) * Z_SOLAR
        if(Z_metal < 0.0d0) Z_metal = 0.0d0
        if(Z_metal > 1.0d0) Z_metal = 1.0d0
        f_scale = (1.0d0 - Z_metal) / (X_H + Y_He)
        mf_X = X_H * f_scale
        rho = nH_val * m_H / mf_X
        do k = 1, nt
           tab_net_rate(k, j, i) = -raw_rate(k, j, i) / rho
        end do
     end do
  end do
  deallocate(raw_rate)


  tab_loaded = .true.

  write(*,'(A,A)') ' Eunha cooling table loaded: ', trim(filename)
  write(*,'(A,I4,A,I4,A,I4)') '   Grid: L=', nd, ' M=', nm, ' N=', nt
  write(*,'(A,F6.2)') '   Redshift: ', z_val
  write(*,'(A,ES10.3,A,ES10.3)') '   Density: ', tab_density(1), ' - ', tab_density(nd)
  write(*,'(A,ES10.3,A,ES10.3)') '   Z_solar: ', tab_metallicity(1), ' - ', tab_metallicity(nm)
  write(*,'(A,ES10.3,A,ES10.3)') '   Temperature: ', tab_temperature(1), ' - ', tab_temperature(nt)

end subroutine eunha_load_table

!###########################################################
! T -> u (specific internal energy)
! Eq. (3): u = kT / ((gamma-1) * mu * m_H)
!###########################################################
function T_to_u(T, mu) result(u)
  implicit none
  real(dp), intent(in) :: T, mu
  real(dp) :: u
  u = k_B * T / ((GAMMA_AD - 1.0d0) * mu * m_H)
end function T_to_u

!###########################################################
! Simple mean molecular weight (smooth ionization transition)
! Used as fallback; primary is Grackle mmw table
!###########################################################
function calc_mu_simple(T, Z_met) result(mu)
  implicit none
  real(dp), intent(in) :: T, Z_met
  real(dp) :: mu
  real(dp) :: Zabs, f_scale, mf_X, mf_Y, mf_Z
  real(dp) :: denom_n, mu_n, denom_i, mu_i, f

  Zabs = Z_met * Z_SOLAR
  Zabs = max(0.0d0, min(1.0d0, Zabs))
  f_scale = (1.0d0 - Zabs) / (X_H + Y_He)
  mf_X = X_H * f_scale
  mf_Y = Y_He * f_scale
  mf_Z = Zabs

  ! Neutral
  denom_n = mf_X + mf_Y / 4.0d0 + mf_Z / A_metal
  mu_n = 1.0d0 / denom_n
  ! Ionized
  denom_i = 2.0d0 * mf_X + 0.75d0 * mf_Y + mf_Z / 2.0d0
  mu_i = 1.0d0 / denom_i

  if(T <= T_ion_lo) then
     mu = mu_n
  else if(T >= T_ion_hi) then
     mu = mu_i
  else
     f = (T - T_ion_lo) / (T_ion_hi - T_ion_lo)
     f = 0.5d0 * (1.0d0 - cos(PI_VAL * f))
     mu = mu_n * (1.0d0 - f) + mu_i * f
  end if
end function calc_mu_simple

!###########################################################
! Interpolate mmw from Grackle table (trilinear in log-space)
! For use during table building (at exact grid points)
!###########################################################
function interp_mmw_cell(nH, Z_met, T) result(mu)
  implicit none
  real(dp), intent(in) :: nH, Z_met, T
  real(dp) :: mu
  integer :: i, j, k
  real(dp) :: f_n, f_z, f_t
  real(dp) :: c000, c001, c010, c011, c100, c101, c110, c111
  real(dp) :: c00, c01, c10, c11, c0, c1

  call find_index_log(tab_density, tab_L, nH, i, f_n)
  call find_index_log(tab_metallicity, tab_M, Z_met, j, f_z)
  call find_index_log(tab_temperature, tab_N, T, k, f_t)

  c000 = tab_mmw(k,   j,   i  )
  c001 = tab_mmw(k+1, j,   i  )
  c010 = tab_mmw(k,   j+1, i  )
  c011 = tab_mmw(k+1, j+1, i  )
  c100 = tab_mmw(k,   j,   i+1)
  c101 = tab_mmw(k+1, j,   i+1)
  c110 = tab_mmw(k,   j+1, i+1)
  c111 = tab_mmw(k+1, j+1, i+1)

  c00 = c000 + f_t * (c001 - c000)
  c01 = c010 + f_t * (c011 - c010)
  c10 = c100 + f_t * (c101 - c100)
  c11 = c110 + f_t * (c111 - c110)

  c0 = c00 + f_z * (c01 - c00)
  c1 = c10 + f_z * (c11 - c10)

  mu = c0 + f_n * (c1 - c0)
  if(mu < 0.5d0) mu = 0.5d0
  if(mu > 1.5d0) mu = 1.5d0
end function interp_mmw_cell

!###########################################################
! Trilinear interpolation of net cooling rate
!###########################################################
function interp_net_rate(nH, Z_met, T) result(rate)
  implicit none
  real(dp), intent(in) :: nH, Z_met, T
  real(dp) :: rate
  integer :: i, j, k
  real(dp) :: f_n, f_z, f_t
  real(dp) :: c000, c001, c010, c011, c100, c101, c110, c111
  real(dp) :: c00, c01, c10, c11, c0, c1

  call find_index_log(tab_density, tab_L, nH, i, f_n)
  call find_index_log(tab_metallicity, tab_M, Z_met, j, f_z)
  call find_index_log(tab_temperature, tab_N, T, k, f_t)

  c000 = tab_net_rate(k,   j,   i  )
  c001 = tab_net_rate(k+1, j,   i  )
  c010 = tab_net_rate(k,   j+1, i  )
  c011 = tab_net_rate(k+1, j+1, i  )
  c100 = tab_net_rate(k,   j,   i+1)
  c101 = tab_net_rate(k+1, j,   i+1)
  c110 = tab_net_rate(k,   j+1, i+1)
  c111 = tab_net_rate(k+1, j+1, i+1)

  c00 = c000 + f_t * (c001 - c000)
  c01 = c010 + f_t * (c011 - c010)
  c10 = c100 + f_t * (c101 - c100)
  c11 = c110 + f_t * (c111 - c110)

  c0 = c00 + f_z * (c01 - c00)
  c1 = c10 + f_z * (c11 - c10)

  rate = c0 + f_n * (c1 - c0)
end function interp_net_rate

!###########################################################
! Find index in log-spaced array (binary search + fine-tune)
! Returns idx in [1, n-1] and fraction f in [0,1]
!###########################################################
subroutine find_index_log(arr, n, val, idx, frac)
  implicit none
  real(dp), intent(in)  :: arr(:)
  integer,  intent(in)  :: n
  real(dp), intent(in)  :: val
  integer,  intent(out) :: idx
  real(dp), intent(out) :: frac
  real(dp) :: log_val, log_lo, log_hi
  integer :: lo, hi, mid

  if(val <= arr(1)) then
     idx = 1; frac = 0.0d0; return
  end if
  if(val >= arr(n)) then
     idx = n - 1; frac = 1.0d0; return
  end if

  ! Binary search
  lo = 1; hi = n
  do while(hi - lo > 1)
     mid = (lo + hi) / 2
     if(arr(mid) > val) then
        hi = mid
     else
        lo = mid
     end if
  end do
  idx = lo

  ! Log-space fraction
  log_val = log10(max(val, 1.0d-30))
  log_lo  = log10(max(arr(idx),   1.0d-30))
  log_hi  = log10(max(arr(idx+1), 1.0d-30))
  if(abs(log_hi - log_lo) > 1.0d-30) then
     frac = (log_val - log_lo) / (log_hi - log_lo)
  else
     frac = 0.0d0
  end if
  frac = max(0.0d0, min(1.0d0, frac))
end subroutine find_index_log

!###########################################################
! T/mu -> T conversion using Grackle mmw table (iterative)
!###########################################################
function T2_to_T(T2, nH, Z_met) result(T)
  implicit none
  real(dp), intent(in) :: T2, nH, Z_met
  real(dp) :: T
  real(dp) :: mu
  integer :: iter

  ! Initial guess with simple mu
  mu = calc_mu_simple(T2 * 0.6d0, Z_met)
  T  = T2 * mu

  ! Iterate with Grackle mmw
  do iter = 1, 3
     mu = interp_mmw_cell(nH, Z_met, T)
     T  = T2 * mu
  end do
end function T2_to_T

!###########################################################
! Main exact integration for a single cell
! Input: nH [cm^-3], Z [Z_solar], specific internal energy [erg/g], dt [s]
! Output: specific internal energy [erg/g]
!###########################################################
! Integrate du/dt=-Lambda(u) exactly for the piecewise-linear, bilinearly
! interpolated energy/rate curve. No interpolation of separately normalized
! Y functions, no corner-dependent equilibrium, and no chemical state solve.
function exact_energy(nH, Z_met, u_initial, dt) result(u)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  real(dp), intent(in) :: nH, Z_met, u_initial, dt
  real(dp) :: u, energy(tab_N), rate(tab_N), remaining
  real(dp) :: slope, current, target, endpoint, travel, exponent, factor,fn,fz,mu
  integer :: k, step,i,j
  if(.not.tab_loaded) error stop 'Eunha: cooling table not loaded'
  if(.not.ieee_is_finite(dt).or.dt<0) error stop 'Eunha: invalid timestep'
  if(.not.ieee_is_finite(nH).or.nH<=0.or..not.ieee_is_finite(Z_met).or.Z_met<0.or. &
     .not.ieee_is_finite(u_initial).or.u_initial<=0) error stop 'Eunha: invalid cell state'
  u=u_initial
  if(dt==0) return
  call find_index_log(tab_density,tab_L,nH,i,fn)
  call find_index_log(tab_metallicity,tab_M,Z_met,j,fz)
  do k=1,tab_N
     mu=(1-fn)*((1-fz)*tab_mmw(k,j,i)+fz*tab_mmw(k,j+1,i)) &
          +fn*((1-fz)*tab_mmw(k,j,i+1)+fz*tab_mmw(k,j+1,i+1))
     if(.not.ieee_is_finite(mu).or.mu<=0) error stop 'Eunha: invalid table molecular weight'
     energy(k)=T_to_u(tab_temperature(k),mu)
     rate(k)=(1-fn)*((1-fz)*tab_net_rate(k,j,i)+fz*tab_net_rate(k,j+1,i)) &
          +fn*((1-fz)*tab_net_rate(k,j,i+1)+fz*tab_net_rate(k,j+1,i+1))
     if(.not.ieee_is_finite(rate(k))) error stop 'Eunha: invalid table rate'
  enddo
  if(any(energy(2:)<=energy(:tab_N-1))) error stop 'Eunha: nonmonotone energy table'
  u=max(energy(1),min(u,energy(tab_N)))
  remaining=dt
  do step=1,tab_N+1
     k=1
     do while(k<tab_N-1)
        if(u<energy(k+1)) exit
        k=k+1
     enddo
     current=rate(k)+(rate(k+1)-rate(k))*(u-energy(k))/(energy(k+1)-energy(k))
     if(abs(current)<EPSILON) return
     ! At a knot select the segment in the direction of evolution.
     if(current>0.and.k>1.and.u==energy(k)) k=k-1
     slope=(rate(k+1)-rate(k))/(energy(k+1)-energy(k))
     current=rate(k)+slope*(u-energy(k))
     if(current>0) then
        target=energy(k); endpoint=rate(k)
     else
        target=energy(k+1); endpoint=rate(k+1)
     endif
     if(target==u) return ! outward motion at the table boundary
     travel=huge(1.0_dp)
     if(current*endpoint>0) then
        if(abs(slope)*(energy(k+1)-energy(k))<=1d-12*abs(current)) then
           travel=abs((target-u)/current)
        else
           travel=log(abs(current/endpoint))/slope
        endif
     endif
     if(remaining>=travel) then
        u=target
        remaining=remaining-travel
        if(remaining<=0) return
     else
        exponent=-slope*remaining
        if(abs(exponent)<1d-5) then
           factor=remaining*(1+exponent/2+exponent**2/6+exponent**3/24+exponent**4/120)
        else
           factor=(1-exp(exponent))/slope
        endif
        u=u-current*factor
        ! A zero of Lambda is approached asymptotically, never crossed.
        if(current*endpoint<=0.and.slope/=0) then
           target=energy(k)-rate(k)/slope
        endif
        u=max(min(u,max(u_initial,target)),min(u_initial,target))
        u=max(energy(k),min(u,energy(k+1)))
        return
     endif
  enddo
  error stop 'Eunha: segment traversal failed'
end function exact_energy

subroutine eunha_solve(nH, T2, Zsolar, dt, delta_T2, n)
  implicit none
  integer, intent(in) :: n
  real(dp), intent(in) :: dt,nH(n),T2(n),Zsolar(n)
  real(dp), intent(out) :: delta_T2(n)
  integer :: ic
  real(dp), parameter :: specific_heat=k_B/((GAMMA_AD-1)*m_H)
  do ic=1,n
     delta_T2(ic)=exact_energy(nH(ic),Zsolar(ic),specific_heat*T2(ic),dt)/specific_heat-T2(ic)
  enddo
end subroutine eunha_solve

!###########################################################
! Reset comparison statistics
!###########################################################
subroutine eunha_reset_stats()
  implicit none
  acc_sum_rdiff   = 0.0d0
  acc_max_rdiff   = 0.0d0
  acc_sum_dT_orig = 0.0d0
  acc_sum_dT_exact= 0.0d0
  acc_ncell       = 0
  acc_nsign       = 0
end subroutine eunha_reset_stats

!###########################################################
! Report comparison statistics (MPI allreduce + print)
!###########################################################
subroutine eunha_report(nstep)
  use amr_commons, only: myid
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer, intent(in) :: nstep
  real(dp) :: g_sum_rdiff, g_max_rdiff, g_sum_dT_orig, g_sum_dT_exact
  integer(8) :: g_ncell, g_nsign
  integer :: ierr
  real(dp) :: mean_rdiff, sign_pct

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(acc_sum_rdiff,  g_sum_rdiff,  1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(acc_max_rdiff,  g_max_rdiff,  1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(acc_sum_dT_orig, g_sum_dT_orig, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(acc_sum_dT_exact,g_sum_dT_exact,1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(acc_ncell,      g_ncell,      1, MPI_INTEGER8,          MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(acc_nsign,      g_nsign,      1, MPI_INTEGER8,          MPI_SUM, MPI_COMM_WORLD, ierr)
#else
  g_sum_rdiff   = acc_sum_rdiff
  g_max_rdiff   = acc_max_rdiff
  g_sum_dT_orig = acc_sum_dT_orig
  g_sum_dT_exact= acc_sum_dT_exact
  g_ncell       = acc_ncell
  g_nsign       = acc_nsign
#endif

  if(myid == 1 .and. g_ncell > 0) then
     mean_rdiff = g_sum_rdiff / dble(g_ncell)
     sign_pct   = 100.0d0 * dble(g_nsign) / dble(g_ncell)
     write(*,'(A,I6,A,I12)') &
          ' Cooling compare step=', nstep, '  ncell=', g_ncell
     write(*,'(A,ES10.3,A,ES10.3)') &
          '   mean_rel_diff=', mean_rdiff, '  max_rel_diff=', g_max_rdiff
     write(*,'(A,ES10.3,A,ES10.3,A,F5.1,A)') &
          '   <dT2_original>=', g_sum_dT_orig/dble(g_ncell), &
          '  <dT2_exact>=', g_sum_dT_exact/dble(g_ncell), &
          '  sign_mismatch=', sign_pct, '%'
  end if

  ! Reset for next step
  call eunha_reset_stats()

end subroutine eunha_report

!###########################################################
! Load multi-z or single-z Grackle binary table
!
! Auto-detects format:
!   Multi-z: first int32 is n_z (small number, <= 200)
!            followed by L, M, N ints
!   Single-z: first int32 is L (=60), next M (=40), N (=256)
!
! Multi-z binary format:
!   int32:           n_z
!   int32:           L, M, N
!   float64[n_z]:    z_values (sorted ascending)
!   float64[7]:      log_dmin, log_dmax, log_zmin, log_zmax, log_tmin, log_tmax, unused
!   float64[N]:      log10(temperature)
!   For iz = 1..n_z:
!     float64[L*M*N]:  net_rate [erg/cm^3/s]
!     float64[L*M*N]:  mmw
!###########################################################
subroutine eunha_load_multi_z(filename)
  use amr_commons, only: myid
  implicit none
  character(len=*), intent(in) :: filename
  integer :: iunit, ierr, iz, i, j, k
  integer :: first_int, nd, nm, nt, nz
  real(dp) :: log_dmin, log_dmax, log_zmin, log_zmax, log_tmin, log_tmax, unused
  real(dp) :: dstep, zstep
  real(dp), allocatable :: log_T(:)
  real(dp), allocatable :: raw_rate(:,:,:), raw_mmw(:,:,:)
  real(dp) :: nH_val, Z_metal, rho, mf_X, f_scale
  logical  :: is_multi_z

  iunit = 42
  open(unit=iunit, file=filename, status='old', access='stream', &
       form='unformatted', iostat=ierr)
  if(ierr /= 0) then
     write(*,*) 'ERROR: eunha_load_multi_z: cannot open ', trim(filename)
     return
  end if

  ! Read first int to detect format
  read(iunit) first_int

  ! Heuristic: single-z files have L=60 as first int
  ! Multi-z files have n_z as first int (small, typically 5-200)
  ! We distinguish by reading 3 more ints and checking consistency
  if(first_int <= 200) then
     ! Could be multi-z: read L, M, N
     read(iunit) nd, nm, nt
     ! Verify: if nd=60, nm=40, nt=256 this is multi-z with n_z=first_int
     ! If first_int=60, nd=40, nt=256, this is single-z
     if(nd >= 10 .and. nm >= 10 .and. nt >= 10 .and. first_int /= nd) then
        is_multi_z = .true.
        nz = first_int
     else
        ! Reinterpret: first_int=L, nd=M, nm=N (single-z)
        is_multi_z = .false.
     end if
  else
     is_multi_z = .false.
  end if

  close(iunit)

  if(.not. is_multi_z) then
     ! Delegate to single-z loader
     if(myid==1) write(*,'(A)') ' Detected single-z Grackle table format'
     call eunha_load_table(filename)
     tab_Nz = 0
     return
  end if

  ! --- Multi-z loading ---
  if(myid==1) write(*,'(A,I3,A)') ' Detected multi-z Grackle table (', nz, ' snapshots)'

  open(unit=iunit, file=filename, status='old', access='stream', &
       form='unformatted', iostat=ierr)

  ! Skip first_int (n_z) and L,M,N — already read
  read(iunit) nz
  read(iunit) nd, nm, nt
  tab_L = nd
  tab_M = nm
  tab_N = nt
  tab_Nz = nz

  ! Read z_values
  if(allocated(tab_z_values)) deallocate(tab_z_values)
  allocate(tab_z_values(nz))
  read(iunit) tab_z_values(1:nz)

  ! The bundled multi-z file retains the single-z seven-double header:
  ! placeholder redshift first, then the six log-domain bounds.
  read(iunit) unused, log_dmin, log_dmax, log_zmin, log_zmax, log_tmin, log_tmax

  ! Read temperature grid
  allocate(log_T(nt))
  read(iunit) log_T(1:nt)

  ! Allocate working arrays (same as single-z)
  if(allocated(tab_density))     deallocate(tab_density)
  if(allocated(tab_metallicity)) deallocate(tab_metallicity)
  if(allocated(tab_temperature)) deallocate(tab_temperature)
  if(allocated(tab_net_rate))    deallocate(tab_net_rate)
  if(allocated(tab_mmw))        deallocate(tab_mmw)

  allocate(tab_density(nd))
  allocate(tab_metallicity(nm))
  allocate(tab_temperature(nt))
  allocate(tab_net_rate(nt, nm, nd))
  allocate(tab_mmw(nt, nm, nd))

  ! Build density grid
  dstep = (log_dmax - log_dmin) / dble(nd - 1)
  do i = 1, nd
     tab_density(i) = 10.0d0 ** (log_dmin + dble(i-1) * dstep)
  end do

  ! Build metallicity grid
  zstep = (log_zmax - log_zmin) / dble(nm - 1)
  do j = 1, nm
     tab_metallicity(j) = 10.0d0 ** (log_zmin + dble(j-1) * zstep)
  end do

  ! Temperature grid
  do k = 1, nt
     tab_temperature(k) = 10.0d0 ** log_T(k)
  end do
  deallocate(log_T)

  ! Allocate multi-z storage (store as erg/g/s after conversion)
  if(allocated(tab_net_rate_z)) deallocate(tab_net_rate_z)
  if(allocated(tab_mmw_z))     deallocate(tab_mmw_z)
  allocate(tab_net_rate_z(nt, nm, nd, nz))
  allocate(tab_mmw_z(nt, nm, nd, nz))

  ! Read each z snapshot
  allocate(raw_rate(nt, nm, nd))
  allocate(raw_mmw(nt, nm, nd))

  do iz = 1, nz
     read(iunit) raw_rate
     read(iunit) raw_mmw

     ! Store mmw directly
     tab_mmw_z(:,:,:,iz) = raw_mmw(:,:,:)

     ! Convert net_rate: erg/cm^3/s -> erg/g/s
     ! Negate sign: Grackle convention (negative=cooling) -> RAMSES (positive=cooling)
     do i = 1, nd
        nH_val = tab_density(i)
        do j = 1, nm
           Z_metal = tab_metallicity(j) * Z_SOLAR
           if(Z_metal < 0.0d0) Z_metal = 0.0d0
           if(Z_metal > 1.0d0) Z_metal = 1.0d0
           f_scale = (1.0d0 - Z_metal) / (X_H + Y_He)
           mf_X = X_H * f_scale
           rho = nH_val * m_H / mf_X
           do k = 1, nt
              tab_net_rate_z(k, j, i, iz) = -raw_rate(k, j, i) / rho
           end do
        end do
     end do
  end do

  deallocate(raw_rate, raw_mmw)
  close(iunit)

  ! Initialize working tables to first z snapshot
  tab_current_z = -1.0d0
  call eunha_interp_redshift(tab_z_values(1))

  tab_loaded = .true.

  if(myid==1) then
     write(*,'(A,A)') ' Eunha multi-z cooling table loaded: ', trim(filename)
     write(*,'(A,I4,A,I4,A,I4,A,I3)') '   Grid: L=', nd, ' M=', nm, ' N=', nt, '  Nz=', nz
     write(*,'(A,F6.2,A,F6.2)') '   Redshift range: ', tab_z_values(1), ' - ', tab_z_values(nz)
     write(*,'(A,ES10.3,A,ES10.3)') '   Density: ', tab_density(1), ' - ', tab_density(nd)
     write(*,'(A,ES10.3,A,ES10.3)') '   Z_solar: ', tab_metallicity(1), ' - ', tab_metallicity(nm)
     write(*,'(A,ES10.3,A,ES10.3)') '   Temperature: ', tab_temperature(1), ' - ', tab_temperature(nt)
  end if

end subroutine eunha_load_multi_z

!###########################################################
! Interpolate multi-z table to target redshift
! Updates tab_net_rate, tab_mmw, and rebuilds Y-function
!###########################################################
subroutine eunha_interp_redshift(z_target)
  use amr_commons, only: myid
  implicit none
  real(dp), intent(in) :: z_target
  integer :: iz, i, j, k
  real(dp) :: w, z_lo, z_hi

  ! Skip if no multi-z data loaded
  if(tab_Nz <= 0) return

  ! Skip if z hasn't changed enough
  if(abs(z_target - tab_current_z) < 0.01d0) return

  ! Find bracket: tab_z_values(iz) <= z_target < tab_z_values(iz+1)
  if(z_target <= tab_z_values(1)) then
     iz = 1
     w  = 0.0d0
  else if(z_target >= tab_z_values(tab_Nz)) then
     iz = tab_Nz - 1
     w  = 1.0d0
  else
     ! Linear search (Nz is small, ~21)
     iz = 1
     do i = 1, tab_Nz - 1
        if(tab_z_values(i+1) > z_target) then
           iz = i
           exit
        end if
     end do
     z_lo = tab_z_values(iz)
     z_hi = tab_z_values(iz + 1)
     w = (z_target - z_lo) / (z_hi - z_lo)
  end if

  ! Linear interpolation of net_rate and mmw
  !$omp parallel do collapse(3) private(k,j,i)
  do i = 1, tab_L
     do j = 1, tab_M
        do k = 1, tab_N
           tab_net_rate(k,j,i) = (1.0d0-w)*tab_net_rate_z(k,j,i,iz) &
                               + w*tab_net_rate_z(k,j,i,iz+1)
           tab_mmw(k,j,i)      = (1.0d0-w)*tab_mmw_z(k,j,i,iz) &
                               + w*tab_mmw_z(k,j,i,iz+1)
        end do
     end do
  end do


  if(myid==1) write(*,'(A,F7.3,A,I2,A,F5.3)') &
       ' Eunha: interpolated to z=', z_target, '  bracket iz=', iz, '  w=', w

  tab_current_z = z_target

end subroutine eunha_interp_redshift

end module eunha_cooling_mod
