program fof_backend_equivalence_driver
  use AA, only: dp, vrel_merge, xsink, msink, vsink, rmerge
  implicit none
  include 'mpif.h'

  integer, parameter :: nsink = 9
  real(dp), parameter :: fof_link = 0.07d0
  integer :: ierr, rank, nrank, nbad, nbad_global
  integer :: nseq, nomp, i, j
  integer, allocatable :: pseq(:), gseq(:), pomp(:), gomp(:)
  integer, parameter :: expected(nsink) = [1, 1, 1, 2, 2, 2, 3, 3, 4]
  integer, external :: Do_FoF_Sequential, Omp_Do_FoF_Tree

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nrank, ierr)

  if (nrank /= 2) then
    if (rank == 0) write(*, '(A,I0)') 'FAIL: expected 2 MPI ranks; got ', nrank
    call MPI_Abort(MPI_COMM_WORLD, 2, ierr)
  endif

  allocate(xsink(nsink, 3), msink(nsink), vsink(nsink, 3))
  allocate(pseq(nsink), gseq(nsink), pomp(nsink), gomp(nsink))
  xsink = 0.5d0
  xsink(:, 1) = [0.02d0, 0.98d0, 0.06d0, 0.20d0, 0.26d0, &
                 0.32d0, 0.60d0, 0.66d0, 0.85d0]
  msink = 1.0d0
  vsink = 0.0d0
  vrel_merge = .false.
  rmerge = 1.0d0

  do i = 1, nsink
    pseq(i) = i
    pomp(i) = i
    gseq(i) = 0
    gomp(i) = 0
  enddo

  nseq = Do_FoF_Sequential(nsink, pseq, gseq, fof_link**2, 1.0d0)
  nomp = Omp_Do_FoF_Tree(nsink, pomp, gomp, fof_link**2, 1.0d0)

  nbad = 0
  if (nseq /= 4 .or. nomp /= 4) nbad = nbad + 1
  do i = 1, nsink
    if (gseq(i) <= 0 .or. gomp(i) <= 0) nbad = nbad + 1
    do j = 1, nsink
      if ((gseq(i) == gseq(j)) .neqv. (expected(i) == expected(j))) nbad = nbad + 1
      if ((gomp(i) == gomp(j)) .neqv. (expected(i) == expected(j))) nbad = nbad + 1
      if ((gseq(i) == gseq(j)) .neqv. (gomp(i) == gomp(j))) nbad = nbad + 1
    enddo
  enddo
  call MPI_Allreduce(nbad, nbad_global, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

  if (nbad_global /= 0) then
    if (rank == 0) then
      write(*, '(A,I0,A,I0)') 'FAIL: sequential groups=', nseq, ' OMP+MPI groups=', nomp
      write(*, '(A,*(1X,I0))') 'expected: ', expected
      write(*, '(A,*(1X,I0))') 'sequential:', gseq
      write(*, '(A,*(1X,I0))') 'OMP+MPI:   ', gomp
      flush(6)
    endif
    call MPI_Finalize(ierr)
    stop 1
  endif

  if (rank == 0) then
    write(*, '(A)') 'PASS: sequential and OMP+MPI FoF partitions match the expected four components.'
    write(*, '(A)') 'Coverage: periodic wrap, transitive links, thread/rank-crossing component, isolated sink.'
  endif
  call MPI_Finalize(ierr)
end program fof_backend_equivalence_driver
