! [RESIZABLE] Verify the post-MPI_Init CPU mask used by each Slurm task.
program mpi_affinity_probe
  use mpi
  implicit none

  integer :: ierr, rank, nproc, unit, ios, env_status, localid, expected_nproc
  character(len=1024) :: output_dir, output_file, line, cpu_mask
  character(len=64) :: localid_text, expected_text

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)

  call get_command_argument(1,output_dir)
  expected_nproc=32
  expected_text=''
  call get_command_argument(2,expected_text)
  if(len_trim(expected_text)>0)then
     read(expected_text,*,iostat=ios)expected_nproc
     if(ios/=0 .or. expected_nproc<1)then
        call MPI_Abort(MPI_COMM_WORLD,17,ierr)
     endif
  endif
  if(len_trim(output_dir)==0 .or. nproc/=expected_nproc)then
     call MPI_Abort(MPI_COMM_WORLD,10,ierr)
  endif

  localid_text=''
  call get_environment_variable('SLURM_LOCALID',localid_text, &
       status=env_status)
  if(env_status/=0 .or. len_trim(localid_text)==0)then
     call MPI_Abort(MPI_COMM_WORLD,14,ierr)
  endif
  read(localid_text,*,iostat=ios)localid
  if(ios/=0)call MPI_Abort(MPI_COMM_WORLD,15,ierr)
  if(localid<0 .or. localid>=nproc)then
     call MPI_Abort(MPI_COMM_WORLD,16,ierr)
  endif

  cpu_mask=''
  open(newunit=unit,file='/proc/self/status',status='old',action='read', &
       iostat=ios)
  if(ios/=0)call MPI_Abort(MPI_COMM_WORLD,11,ierr)
  do
     read(unit,'(A)',iostat=ios)line
     if(ios/=0)exit
     if(index(line,'Cpus_allowed_list:')==1)then
        cpu_mask=adjustl(line(len('Cpus_allowed_list:')+1:))
        exit
     endif
  enddo
  close(unit)
  if(len_trim(cpu_mask)==0)call MPI_Abort(MPI_COMM_WORLD,12,ierr)

  write(output_file,'(A,"/rank_",I4.4,".txt")')trim(output_dir),rank
  open(newunit=unit,file=trim(output_file),status='replace',action='write', &
       iostat=ios)
  if(ios/=0)call MPI_Abort(MPI_COMM_WORLD,13,ierr)
  write(unit,'(I0,1X,I0,1X,A)')rank,localid,trim(cpu_mask)
  close(unit)

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call MPI_Finalize(ierr)
end program mpi_affinity_probe
