module pbh_commons
  !=====================================================================
  ! Evaporating-PBH dark matter: global tables and runtime parameters.
  !
  ! Design (paper appendix A): ONE global mass function, all DM particles
  ! share the mixed-mass factor w(a) = 1 - f + f*g(a), heating from the
  ! cumulative kernel Qtilde(a) = int q*g dt. Everything here is exact
  ! table arithmetic; no per-particle state, no evaporation CFL.
  !
  ! This module is deliberately dependency-light (only amr_parameters)
  ! so it can be unit-tested standalone. All run context (aexp, myid,
  ! nrestart, ...) is passed in as arguments by pbh_evap_fine.
  !=====================================================================
  use amr_parameters, only: dp
  implicit none

  ! ---- &PBH_PARAMS namelist parameters (read in read_params) ----
  logical            :: use_pbh              = .false.
  character(LEN=512) :: pbh_table_file       = 'pbh_evap_table.dat'
  real(dp)           :: pbh_fraction         = 1.0d0   ! f_PBH of the DM
  real(dp)           :: pbh_boost            = 1.0d0   ! linear-response boost, q only
  character(LEN=16)  :: pbh_energy_sink      = 'local_heat' ! local_heat | removed
  real(dp)           :: pbh_bkg_warn         = 1.0d-3  ! warn if f*(1-g(z=0)) exceeds
  logical            :: pbh_check_provenance = .true.

  ! ---- loaded table (t stored in seconds) ----
  integer :: pbh_n = 0
  real(dp), allocatable :: tab_a(:), tab_la(:), tab_t(:), tab_g(:)
  real(dp), allocatable :: tab_q(:), tab_lam(:), tab_qc(:)
  integer(kind=8)    :: pbh_cksum = 0        ! adler32 of the data section
  character(LEN=256) :: pbh_model_line = ''  ! "# model: ..." header line
  logical :: pbh_table_loaded = .false.

  ! ---- run normalisation (persisted through pbh_provenance.txt) ----
  logical  :: pbh_ready  = .false.
  real(dp) :: pbh_anorm  = 1.0d0   ! scale factor where g is normalised (sim start)
  real(dp) :: pbh_gnorm  = 1.0d0   ! g_table(pbh_anorm)
  real(dp) :: pbh_qcnorm = 0.0d0   ! Qtilde_table(pbh_anorm)

  ! ---- diagnostics (accumulated by pbh_evap_fine; einj in erg) ----
  real(dp)        :: pbh_einj_loc     = 0.0d0
  real(dp)        :: pbh_einj_tot     = 0.0d0
  integer(kind=8) :: pbh_nfallback_loc = 0

  ! ---- per-level epoch bookkeeping: the STARTING scale factor of the
  !      current level step, recorded by amr_step (via pbh_mark_level)
  !      before the fine-level recursion advances aexp. Every step ratio
  !      then spans the exact [a_start, a_end] of its own step, including
  !      the very first step and the first step after a restart ----
  real(dp), allocatable :: pbh_aold(:)

  ! ---- remote-deposit buffers for the deterministic thread-safe deposit
  !      scheme (owner-writes + deferred buffer, no atomics): deposits whose
  !      target cell is not owned by the particle's own grid (parent-cell
  !      fallback, boundary-drift neighbours) are queued per thread and
  !      applied serially after the parallel loop, sorted by (cell, id) so
  !      the result is bitwise independent of the thread count ----
  type :: pbh_rbuf_t
     integer :: n = 0
     integer,         allocatable :: icell(:)
     integer(kind=8), allocatable :: pid(:)
     real(dp),        allocatable :: de(:)
  end type pbh_rbuf_t
  type(pbh_rbuf_t), allocatable :: pbh_rbuf(:)

  real(dp), parameter :: GYR2S = 3.155760d16  ! matches make_pbh_tables.py

  character(LEN=*), parameter :: pbh_prov_file = 'pbh_provenance.txt'

contains

  !=====================================================================
  subroutine pbh_read_table(myid)
    ! Load the 7-column table {a z t[Gyr] g q Lambda Qtilde}. Every rank
    ! reads the file directly (same pattern as the namelist itself).
    integer, intent(in) :: myid
    integer :: u, ios, n, i, s1, s2, j, b
    integer(kind=8) :: ck_hdr
    character(LEN=512) :: line
    logical :: ok

    if(pbh_table_loaded) return

    open(newunit=u, file=trim(pbh_table_file), form='formatted', &
         & status='old', action='read', iostat=ios)
    if(ios /= 0) then
       if(myid==1) write(*,*) 'PBH ERROR: cannot open pbh_table_file: ', &
            & trim(pbh_table_file)
       call clean_stop
    end if

    ! -- first pass: count data lines, grab header info --
    n = 0
    ck_hdr = -1
    do
       read(u,'(A)',iostat=ios) line
       if(ios /= 0) exit
       if(line(1:1) == '#') then
          if(index(line,'model:') > 0) pbh_model_line = trim(line)
          j = index(line,'checksum_adler32:')
          if(j > 0) call pbh_parse_hex(line(j+17:), ck_hdr)
       else if(len_trim(line) > 0) then
          n = n + 1
       end if
    end do
    if(n < 4) then
       if(myid==1) write(*,*) 'PBH ERROR: table too short: ', n, ' rows'
       call clean_stop
    end if

    pbh_n = n
    allocate(tab_a(n), tab_la(n), tab_t(n), tab_g(n))
    allocate(tab_q(n), tab_lam(n), tab_qc(n))

    ! -- second pass: parse data + Adler-32 over data bytes (line+'\n') --
    rewind(u)
    i = 0; s1 = 1; s2 = 0
    do
       read(u,'(A)',iostat=ios) line
       if(ios /= 0) exit
       if(line(1:1) == '#' .or. len_trim(line) == 0) cycle
       i = i + 1
       read(line,*,iostat=ios) tab_a(i), tab_la(i), tab_t(i), tab_g(i), &
            & tab_q(i), tab_lam(i), tab_qc(i)   ! tab_la holds z, unused; reset below
       if(ios /= 0) then
          if(myid==1) write(*,*) 'PBH ERROR: bad table row ', i
          call clean_stop
       end if
       do j = 1, len_trim(line)
          b = iachar(line(j:j))
          s1 = mod(s1 + b, 65521)
          s2 = mod(s2 + s1, 65521)
       end do
       s1 = mod(s1 + 10, 65521)    ! newline byte
       s2 = mod(s2 + s1, 65521)
    end do
    close(u)
    pbh_cksum = ior(ishft(int(s2,kind=8),16), int(s1,kind=8))

    ! -- units, derived grids, sanity checks --
    tab_t  = tab_t * GYR2S          ! Gyr -> s
    tab_la = log(tab_a)
    ok = .true.
    do i = 2, n
       if(tab_a(i) <= tab_a(i-1)) ok = .false.
       if(tab_t(i) <= tab_t(i-1)) ok = .false.
       if(tab_g(i) >  tab_g(i-1)*(1.0d0+1.0d-12)) ok = .false.
       if(tab_qc(i) < tab_qc(i-1)) ok = .false.
    end do
    if(.not.ok .or. minval(tab_g) < 0.0d0 .or. maxval(tab_g) > 1.0d0+1.0d-12) then
       if(myid==1) write(*,*) 'PBH ERROR: table not monotone/sane (a,t up; g down; Qc up)'
       call clean_stop
    end if
    if(ck_hdr >= 0 .and. ck_hdr /= pbh_cksum) then
       if(myid==1) write(*,'(A,Z8.8,A,Z8.8)') ' PBH ERROR: table checksum mismatch, header=0x', &
            & ck_hdr, ' computed=0x', pbh_cksum
       call clean_stop
    end if

    pbh_table_loaded = .true.
    if(myid==1) then
       write(*,'(A)') ' PBH: loaded table '//trim(pbh_table_file)
       if(len_trim(pbh_model_line)>0) write(*,'(A)') ' PBH: '//trim(pbh_model_line)
       write(*,'(A,I6,2(A,ES11.4),A,Z8.8)') ' PBH: rows=', pbh_n, &
            & '  a=[', tab_a(1), ',', tab_a(pbh_n), ']  adler32=0x', pbh_cksum
    end if
  end subroutine pbh_read_table

  !=====================================================================
  subroutine pbh_parse_hex(str, val)
    ! Parse the first 0x... token in str into val (Adler-32 fits int8).
    character(LEN=*), intent(in)  :: str
    integer(kind=8),  intent(out) :: val
    integer :: i, j, d
    character :: c
    val = -1
    i = index(str, '0x')
    if(i == 0) return
    val = 0
    do j = i+2, len_trim(str)
       c = str(j:j)
       select case(c)
       case('0':'9'); d = iachar(c) - iachar('0')
       case('a':'f'); d = iachar(c) - iachar('a') + 10
       case('A':'F'); d = iachar(c) - iachar('A') + 10
       case default;  exit
       end select
       val = val*16 + d
    end do
  end subroutine pbh_parse_hex

  !=====================================================================
  function pbh_interp_la(arr, a) result(v)
    ! Linear interpolation of arr on the log-a grid, clamped at the ends.
    real(dp), intent(in) :: arr(:), a
    real(dp) :: v, la, w
    integer  :: j
    la = log(max(a, 1.0d-30))
    if(la <= tab_la(1)) then
       v = arr(1); return
    else if(la >= tab_la(pbh_n)) then
       v = arr(pbh_n); return
    end if
    j = pbh_bisect(tab_la, la)
    w = (la - tab_la(j)) / (tab_la(j+1) - tab_la(j))
    v = arr(j)*(1.0d0-w) + arr(j+1)*w
  end function pbh_interp_la

  function pbh_t_of_a(a) result(t)
    real(dp), intent(in) :: a
    real(dp) :: t
    t = pbh_interp_la(tab_t, a)
  end function pbh_t_of_a

  function pbh_a_of_t(t) result(a)
    ! Inverse map: linear in (t, ln a), clamped.
    real(dp), intent(in) :: t
    real(dp) :: a, w
    integer  :: j
    if(t <= tab_t(1)) then
       a = tab_a(1); return
    else if(t >= tab_t(pbh_n)) then
       a = tab_a(pbh_n); return
    end if
    j = pbh_bisect(tab_t, t)
    w = (t - tab_t(j)) / (tab_t(j+1) - tab_t(j))
    a = exp(tab_la(j)*(1.0d0-w) + tab_la(j+1)*w)
  end function pbh_a_of_t

  function pbh_bisect(arr, x) result(j)
    ! Largest j with arr(j) <= x, for strictly increasing arr; 1<=j<=n-1.
    real(dp), intent(in) :: arr(:), x
    integer :: j, lo, hi, mid
    lo = 1; hi = pbh_n
    do while(hi - lo > 1)
       mid = (lo + hi)/2
       if(arr(mid) <= x) then
          lo = mid
       else
          hi = mid
       end if
    end do
    j = lo
  end function pbh_bisect

  !=====================================================================
  subroutine pbh_lazy_init(aexp_now, aexp_ini, nrestart, myid, cosmo)
    ! One-time run initialisation at the first pbh_evap_fine call.
    ! Fresh run: normalise g/Qtilde at the simulation start and persist
    ! the normalisation in pbh_provenance.txt. Restart: read it back so
    ! the evaporation history embedded in mp stays consistent.
    real(dp), intent(in) :: aexp_now, aexp_ini
    integer,  intent(in) :: nrestart, myid
    logical,  intent(in) :: cosmo
    real(dp) :: a0, gz0, drift
    logical  :: have_prov

    if(pbh_ready) return
    if(.not.cosmo) then
       if(myid==1) write(*,*) 'PBH ERROR: use_pbh requires cosmo=.true.'
       call clean_stop
    end if

    if(nrestart > 0 .and. pbh_check_provenance) then
       call pbh_read_provenance(myid, have_prov)
       if(.not.have_prov) then
          if(myid==1) write(*,*) 'PBH ERROR: restart but no/invalid ', &
               & pbh_prov_file, ' (set pbh_check_provenance=.false. to override)'
          call clean_stop
       end if
    else
       a0 = aexp_ini
       if(a0 <= 0.0d0 .or. a0 > 1.0d0) a0 = aexp_now
       pbh_anorm  = a0
       pbh_gnorm  = pbh_interp_la(tab_g,  a0)
       pbh_qcnorm = pbh_interp_la(tab_qc, a0)
       if(nrestart > 0 .and. myid==1) write(*,*) &
            & 'PBH WARNING: provenance check disabled, renormalising at aexp_ini'
       if(nrestart == 0 .and. myid==1) call pbh_write_provenance
    end if

    ! background-consistency guard (paper appendix A)
    gz0   = pbh_interp_la(tab_g, 1.0d0) / pbh_gnorm
    drift = pbh_fraction * (1.0d0 - min(gz0, 1.0d0))
    if(myid==1) then
       write(*,'(A,ES11.4,A,ES11.4,A,ES11.4)') ' PBH: a_norm=', pbh_anorm, &
            & '  g_norm=', pbh_gnorm, '  g(z=0)/g_norm=', gz0
       write(*,'(A,ES11.4)') ' PBH: expected matter-density drift f*(1-g(z=0)) =', drift
       if(drift > pbh_bkg_warn) write(*,'(A,ES9.2,A)') &
            & ' PBH WARNING: drift exceeds ', pbh_bkg_warn, &
            & ' ; fixed-Omega_m background is a poor approximation for this setup'
       write(*,'(A,ES11.4,A,ES11.4,A)') ' PBH: f_PBH=', pbh_fraction, &
            & '  boost=', pbh_boost, '  sink='//trim(pbh_energy_sink)
    end if
    pbh_ready = .true.
  end subroutine pbh_lazy_init

  !=====================================================================
  subroutine pbh_rbuf_push(b, icell, pid, de)
    ! Append one remote deposit to a thread-private buffer (geometric
    ! growth; only the owning thread ever touches b).
    type(pbh_rbuf_t), intent(inout) :: b
    integer,          intent(in)    :: icell
    integer(kind=8),  intent(in)    :: pid
    real(dp),         intent(in)    :: de
    integer :: cap
    integer,         allocatable :: ti(:)
    integer(kind=8), allocatable :: tp(:)
    real(dp),        allocatable :: td(:)
    if(.not.allocated(b%icell)) then
       allocate(b%icell(64), b%pid(64), b%de(64))
    else if(b%n == size(b%icell)) then
       cap = 2*size(b%icell)
       allocate(ti(cap)); ti(1:b%n)=b%icell(1:b%n); call move_alloc(ti,b%icell)
       allocate(tp(cap)); tp(1:b%n)=b%pid(1:b%n);   call move_alloc(tp,b%pid)
       allocate(td(cap)); td(1:b%n)=b%de(1:b%n);    call move_alloc(td,b%de)
    end if
    b%n = b%n + 1
    b%icell(b%n) = icell
    b%pid(b%n)   = pid
    b%de(b%n)    = de
  end subroutine pbh_rbuf_push

  !=====================================================================
  subroutine pbh_mark_level(ilevel, nlev, a)
    ! Called from amr_step at the START of every level step (before the
    ! fine-level recursion), while aexp still holds the step-start value.
    integer,  intent(in) :: ilevel, nlev
    real(dp), intent(in) :: a
    if(.not.allocated(pbh_aold)) then
       allocate(pbh_aold(nlev))
       pbh_aold = -1.0d0
    end if
    pbh_aold(ilevel) = a
  end subroutine pbh_mark_level

  !=====================================================================
  subroutine pbh_step(ilevel, nlev, a_end, dtsec, ratio, dQ, w0)
    ! Step-global factors for the level step that ENDS at a_end (aexp is
    ! already advanced when pbh_evap_fine runs). The step start a0 was
    ! recorded by pbh_mark_level at step entry, so the interval is exact
    ! for every step; if the mark is missing (pathological) fall back to
    ! the table time map a0 = a(t(a_end)-dtsec).
    !   ratio = w(a1)/w(a0) with w = 1-f+f*g   (exact mixed-mass update)
    !   dQ    = boost * [Qtilde(a1)-Qtilde(a0)]/g_norm   [erg per gram of
    !           initial PBH-share mass]
    !   w0    = w(a0), needed for m_initial = mp/w0
    integer,  intent(in)  :: ilevel, nlev
    real(dp), intent(in)  :: a_end, dtsec
    real(dp), intent(out) :: ratio, dQ, w0
    real(dp) :: a0, a1, g0, g1, w1, f
    if(.not.allocated(pbh_aold)) then
       allocate(pbh_aold(nlev))
       pbh_aold = -1.0d0
    end if
    f  = pbh_fraction
    a1 = a_end
    if(pbh_aold(ilevel) > 0.0d0) then
       a0 = min(pbh_aold(ilevel), a1)
    else
       a0 = min(pbh_a_of_t(pbh_t_of_a(a1) - dtsec), a1)
    end if
    g0 = pbh_interp_la(tab_g, a0) / pbh_gnorm
    g1 = pbh_interp_la(tab_g, a1) / pbh_gnorm
    w0 = 1.0d0 - f + f*g0
    w1 = 1.0d0 - f + f*g1
    if(w0 <= 1.0d-30) then
       ratio = 1.0d0
       dQ    = 0.0d0
       w0    = 1.0d0
    else
       ratio = w1/w0
       dQ    = pbh_boost * (pbh_interp_la(tab_qc, a1) - pbh_interp_la(tab_qc, a0)) &
             &           / pbh_gnorm
    end if
  end subroutine pbh_step

  !=====================================================================
  subroutine pbh_write_provenance
    integer :: u
    open(newunit=u, file=pbh_prov_file, form='formatted', status='replace')
    write(u,'(A)') '# PBH run provenance (written by the lagRamses PBH patch)'
    write(u,'(A)') 'table_file= '//trim(pbh_table_file)
    write(u,'(A,Z8.8)')   'checksum_adler32= 0x', pbh_cksum
    write(u,'(A,ES23.15)') 'pbh_fraction= ', pbh_fraction
    write(u,'(A,ES23.15)') 'pbh_boost= ', pbh_boost
    write(u,'(A,ES23.15)') 'a_norm= ', pbh_anorm
    write(u,'(A,ES23.15)') 'g_norm= ', pbh_gnorm
    write(u,'(A,ES23.15)') 'qc_norm= ', pbh_qcnorm
    close(u)
    write(*,'(A)') ' PBH: wrote '//pbh_prov_file
  end subroutine pbh_write_provenance

  subroutine pbh_read_provenance(myid, ok)
    integer, intent(in)  :: myid
    logical, intent(out) :: ok
    integer :: u, ios, j
    integer(kind=8) :: ck
    real(dp) :: f, b, an, gn, qn
    character(LEN=512) :: line
    ok = .false.
    ck = -1; f = -1.0d0; b = -1.0d0; an = -1.0d0; gn = -1.0d0; qn = 0.0d0
    open(newunit=u, file=pbh_prov_file, form='formatted', status='old', &
         & action='read', iostat=ios)
    if(ios /= 0) return
    do
       read(u,'(A)',iostat=ios) line
       if(ios /= 0) exit
       j = index(line,'=')
       if(j == 0) cycle
       if(index(line,'checksum_adler32') > 0) then
          call pbh_parse_hex(line(j+1:), ck)
       else if(index(line,'pbh_fraction') > 0) then
          read(line(j+1:),*,iostat=ios) f
       else if(index(line,'pbh_boost') > 0) then
          read(line(j+1:),*,iostat=ios) b
       else if(index(line,'a_norm') > 0) then
          read(line(j+1:),*,iostat=ios) an
       else if(index(line,'g_norm') > 0) then
          read(line(j+1:),*,iostat=ios) gn
       else if(index(line,'qc_norm') > 0) then
          read(line(j+1:),*,iostat=ios) qn
       end if
    end do
    close(u)
    if(ck /= pbh_cksum) then
       if(myid==1) write(*,'(A,Z8.8,A,Z8.8)') ' PBH ERROR: provenance table checksum 0x', &
            & ck, ' /= loaded table 0x', pbh_cksum
       return
    end if
    if(abs(f - pbh_fraction) > 1.0d-12*max(1.0d0,abs(f)) .or. &
       abs(b - pbh_boost)    > 1.0d-12*max(1.0d0,abs(b))) then
       if(myid==1) write(*,*) 'PBH ERROR: provenance pbh_fraction/pbh_boost differ from namelist'
       return
    end if
    if(an <= 0.0d0 .or. gn <= 0.0d0) return
    pbh_anorm  = an
    pbh_gnorm  = gn
    pbh_qcnorm = qn
    ok = .true.
    if(myid==1) write(*,'(A)') ' PBH: provenance verified ('//pbh_prov_file//')'
  end subroutine pbh_read_provenance

end module pbh_commons
