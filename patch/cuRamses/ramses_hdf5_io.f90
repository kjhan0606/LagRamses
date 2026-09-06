!###########################################################################
! ramses_hdf5_io.f90 — HDF5 wrapper module for cuRAMSES
!
! Provides helper routines for parallel HDF5 file create/open/close,
! attribute and dataset I/O.
! All routines are guarded by #ifdef HDF5.
!###########################################################################
#ifdef HDF5
module ramses_hdf5_io
  use amr_parameters, only: dp, i8b
  use hdf5
  implicit none

  integer(HID_T) :: hdf5_file_id    ! current file handle
  integer(HID_T) :: hdf5_plist_id   ! file access property list (MPI-IO)

  ! ----- Phase D: file pool (LRU cache of open HDF5 file handles) -----
  ! Future-proofing for multi-file outputs (one file per writer-group).
  ! Single-shared-file callers (current production path) bypass the pool
  ! entirely; nothing changes for them. Callers that want pooling use
  !   call file_pool_init(max_open)
  !   call file_pool_get(filename, comm, file_id)
  !   ... (handle is owned by the pool; do not close it directly)
  !   call file_pool_close_all()
  integer, parameter :: FP_NAME_LEN    = 512
  integer, parameter :: FP_MAX_DEFAULT = 8
  integer, parameter :: FP_MAX_LIMIT   = 1024
  integer                    :: fp_cap   = 0
  integer                    :: fp_count = 0
  integer(i8b)               :: fp_tick  = 0
  character(len=FP_NAME_LEN), allocatable :: fp_name(:)
  integer(HID_T),             allocatable :: fp_file_id(:)
  integer(i8b),               allocatable :: fp_lru(:)
  ! Diagnostic counters
  integer(i8b) :: fp_hits = 0, fp_misses = 0, fp_evictions = 0

contains

  !=========================================================================
  ! File-level operations
  !=========================================================================
  subroutine hdf5_create_parallel(filename, comm)
    implicit none
    include 'mpif.h'
    character(len=*), intent(in) :: filename
    integer, intent(in) :: comm
    integer :: ierr

    call h5open_f(ierr)
    ! Create file access property list for parallel I/O
    call h5pcreate_f(H5P_FILE_ACCESS_F, hdf5_plist_id, ierr)
    call h5pset_fapl_mpio_f(hdf5_plist_id, comm, MPI_INFO_NULL, ierr)
    ! Create file
    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, hdf5_file_id, ierr, &
         access_prp=hdf5_plist_id)
    call h5pclose_f(hdf5_plist_id, ierr)
  end subroutine hdf5_create_parallel

  subroutine hdf5_open_parallel(filename, comm)
    implicit none
    include 'mpif.h'
    character(len=*), intent(in) :: filename
    integer, intent(in) :: comm
    integer :: ierr

    call h5open_f(ierr)
    call h5pcreate_f(H5P_FILE_ACCESS_F, hdf5_plist_id, ierr)
    call h5pset_fapl_mpio_f(hdf5_plist_id, comm, MPI_INFO_NULL, ierr)
    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, hdf5_file_id, ierr, &
         access_prp=hdf5_plist_id)
    call h5pclose_f(hdf5_plist_id, ierr)
  end subroutine hdf5_open_parallel

  subroutine hdf5_close_file()
    implicit none
    integer :: ierr
    call h5fclose_f(hdf5_file_id, ierr)
    call h5close_f(ierr)
  end subroutine hdf5_close_file

  !=========================================================================
  ! Group operations
  !=========================================================================
  subroutine hdf5_create_group(grp_name, grp_id)
    implicit none
    character(len=*), intent(in) :: grp_name
    integer(HID_T), intent(out) :: grp_id
    integer :: ierr
    call h5gcreate_f(hdf5_file_id, trim(grp_name), grp_id, ierr)
  end subroutine hdf5_create_group

  subroutine hdf5_open_group(grp_name, grp_id)
    implicit none
    character(len=*), intent(in) :: grp_name
    integer(HID_T), intent(out) :: grp_id
    integer :: ierr
    call h5gopen_f(hdf5_file_id, trim(grp_name), grp_id, ierr)
  end subroutine hdf5_open_group

  subroutine hdf5_close_group(grp_id)
    implicit none
    integer(HID_T), intent(in) :: grp_id
    integer :: ierr
    call h5gclose_f(grp_id, ierr)
  end subroutine hdf5_close_group

  !=========================================================================
  ! Attribute write helpers
  !=========================================================================
  subroutine hdf5_write_attr_int(loc_id, name, val)
    implicit none
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: val
    integer(HID_T) :: space_id, attr_id
    integer(HSIZE_T), dimension(1) :: dims = (/1/)
    integer :: ierr
    call h5screate_simple_f(1, dims, space_id, ierr)
    call h5acreate_f(loc_id, trim(name), H5T_NATIVE_INTEGER, space_id, attr_id, ierr)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, val, dims, ierr)
    call h5aclose_f(attr_id, ierr)
    call h5sclose_f(space_id, ierr)
  end subroutine hdf5_write_attr_int

  subroutine hdf5_write_attr_int8(loc_id, name, val)
    implicit none
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: name
    integer(i8b), intent(in) :: val
    integer(HID_T) :: space_id, attr_id, type_id
    integer(HSIZE_T), dimension(1) :: dims = (/1/)
    integer :: ierr
    call h5screate_simple_f(1, dims, space_id, ierr)
    call h5tcopy_f(H5T_NATIVE_INTEGER, type_id, ierr)
    call h5tset_size_f(type_id, int(8, SIZE_T), ierr)
    call h5acreate_f(loc_id, trim(name), type_id, space_id, attr_id, ierr)
    call h5awrite_f(attr_id, type_id, val, dims, ierr)
    call h5aclose_f(attr_id, ierr)
    call h5tclose_f(type_id, ierr)
    call h5sclose_f(space_id, ierr)
  end subroutine hdf5_write_attr_int8

  subroutine hdf5_write_attr_dp(loc_id, name, val)
    implicit none
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: val
    integer(HID_T) :: space_id, attr_id
    integer(HSIZE_T), dimension(1) :: dims = (/1/)
    integer :: ierr
    call h5screate_simple_f(1, dims, space_id, ierr)
    call h5acreate_f(loc_id, trim(name), H5T_NATIVE_DOUBLE, space_id, attr_id, ierr)
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, val, dims, ierr)
    call h5aclose_f(attr_id, ierr)
    call h5sclose_f(space_id, ierr)
  end subroutine hdf5_write_attr_dp

  subroutine hdf5_write_attr_string(loc_id, name, val)
    implicit none
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: val
    integer(HID_T) :: space_id, attr_id, type_id
    integer(HSIZE_T), dimension(1) :: dims = (/1/)
    integer(SIZE_T) :: slen
    integer :: ierr
    slen = len_trim(val)
    call h5screate_simple_f(1, dims, space_id, ierr)
    call h5tcopy_f(H5T_NATIVE_CHARACTER, type_id, ierr)
    call h5tset_size_f(type_id, slen, ierr)
    call h5acreate_f(loc_id, trim(name), type_id, space_id, attr_id, ierr)
    call h5awrite_f(attr_id, type_id, trim(val), dims, ierr)
    call h5aclose_f(attr_id, ierr)
    call h5tclose_f(type_id, ierr)
    call h5sclose_f(space_id, ierr)
  end subroutine hdf5_write_attr_string

  !=========================================================================
  ! Attribute read helpers
  !=========================================================================
  subroutine hdf5_read_attr_int(loc_id, name, val)
    implicit none
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: name
    integer, intent(out) :: val
    integer(HID_T) :: attr_id
    integer(HSIZE_T), dimension(1) :: dims = (/1/)
    integer :: ierr
    call h5aopen_f(loc_id, trim(name), attr_id, ierr)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, val, dims, ierr)
    call h5aclose_f(attr_id, ierr)
  end subroutine hdf5_read_attr_int

  subroutine hdf5_read_attr_int8(loc_id, name, val)
    implicit none
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: name
    integer(i8b), intent(out) :: val
    integer(HID_T) :: attr_id, type_id
    integer(HSIZE_T), dimension(1) :: dims = (/1/)
    integer :: ierr
    call h5aopen_f(loc_id, trim(name), attr_id, ierr)
    call h5tcopy_f(H5T_NATIVE_INTEGER, type_id, ierr)
    call h5tset_size_f(type_id, int(8, SIZE_T), ierr)
    call h5aread_f(attr_id, type_id, val, dims, ierr)
    call h5aclose_f(attr_id, ierr)
    call h5tclose_f(type_id, ierr)
  end subroutine hdf5_read_attr_int8

  subroutine hdf5_read_attr_int_checked(loc_id, name, val, status)
    implicit none
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: name
    integer, intent(out) :: val
    integer, intent(out) :: status
    integer(HID_T) :: attr_id, space_id
    integer(HSIZE_T), dimension(1) :: dims = (/1/), maxdims
    integer :: ierr, rank
    logical :: have_attr, have_space

    status = 1
    have_attr = .false.
    have_space = .false.
    call h5aopen_f(loc_id, trim(name), attr_id, ierr)
    if(ierr /= 0) then
       status = 10
       goto 930
    end if
    have_attr = .true.
    call h5aget_space_f(attr_id, space_id, ierr)
    if(ierr /= 0) then
       status = 11
       goto 930
    end if
    have_space = .true.
    call h5sget_simple_extent_ndims_f(space_id, rank, ierr)
    if(ierr /= 0 .or. (rank /= 0 .and. rank /= 1)) then
       status = 12
       goto 930
    end if
    if(rank == 1) then
       call h5sget_simple_extent_dims_f(space_id, dims, maxdims, ierr)
       if(ierr < 0 .or. ierr /= rank .or. dims(1) /= 1_HSIZE_T) then
          status = 13
          goto 930
       end if
    end if
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, val, dims, ierr)
    if(ierr /= 0) then
       status = 14
       goto 930
    end if
    status = 0
930 continue
    if(have_space) call h5sclose_f(space_id, ierr)
    if(have_attr) call h5aclose_f(attr_id, ierr)
  end subroutine hdf5_read_attr_int_checked

  subroutine hdf5_read_attr_int8_checked(loc_id, name, val, status)
    implicit none
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: name
    integer(i8b), intent(out) :: val
    integer, intent(out) :: status
    integer(HID_T) :: attr_id, space_id, type_id
    integer(HSIZE_T), dimension(1) :: dims = (/1/), maxdims
    integer :: ierr, rank
    logical :: have_attr, have_space, have_type

    status = 1
    have_attr = .false.
    have_space = .false.
    have_type = .false.
    call h5aopen_f(loc_id, trim(name), attr_id, ierr)
    if(ierr /= 0) then
       status = 10
       goto 940
    end if
    have_attr = .true.
    call h5aget_space_f(attr_id, space_id, ierr)
    if(ierr /= 0) then
       status = 11
       goto 940
    end if
    have_space = .true.
    call h5sget_simple_extent_ndims_f(space_id, rank, ierr)
    if(ierr /= 0 .or. (rank /= 0 .and. rank /= 1)) then
       status = 12
       goto 940
    end if
    if(rank == 1) then
       call h5sget_simple_extent_dims_f(space_id, dims, maxdims, ierr)
       if(ierr < 0 .or. ierr /= rank .or. dims(1) /= 1_HSIZE_T) then
          status = 13
          goto 940
       end if
    end if
    call h5tcopy_f(H5T_NATIVE_INTEGER, type_id, ierr)
    if(ierr /= 0) then
       status = 14
       goto 940
    end if
    have_type = .true.
    call h5tset_size_f(type_id, int(8, SIZE_T), ierr)
    if(ierr /= 0) then
       status = 15
       goto 940
    end if
    call h5aread_f(attr_id, type_id, val, dims, ierr)
    if(ierr /= 0) then
       status = 16
       goto 940
    end if
    status = 0
940 continue
    if(have_type) call h5tclose_f(type_id, ierr)
    if(have_space) call h5sclose_f(space_id, ierr)
    if(have_attr) call h5aclose_f(attr_id, ierr)
  end subroutine hdf5_read_attr_int8_checked

  subroutine hdf5_read_attr_dp(loc_id, name, val)
    implicit none
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: val
    integer(HID_T) :: attr_id
    integer(HSIZE_T), dimension(1) :: dims = (/1/)
    integer :: ierr
    call h5aopen_f(loc_id, trim(name), attr_id, ierr)
    call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, val, dims, ierr)
    call h5aclose_f(attr_id, ierr)
  end subroutine hdf5_read_attr_dp

  subroutine hdf5_read_attr_dp_checked(loc_id, name, val, status)
    implicit none
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: val
    integer, intent(out) :: status
    integer(HID_T) :: attr_id, space_id
    integer(HSIZE_T), dimension(1) :: dims = (/1/), maxdims
    integer :: ierr, rank
    logical :: have_attr, have_space

    status = 1
    have_attr = .false.
    have_space = .false.
    call h5aopen_f(loc_id, trim(name), attr_id, ierr)
    if(ierr /= 0) then
       status = 10
       goto 950
    end if
    have_attr = .true.
    call h5aget_space_f(attr_id, space_id, ierr)
    if(ierr /= 0) then
       status = 11
       goto 950
    end if
    have_space = .true.
    call h5sget_simple_extent_ndims_f(space_id, rank, ierr)
    if(ierr /= 0 .or. (rank /= 0 .and. rank /= 1)) then
       status = 12
       goto 950
    end if
    if(rank == 1) then
       call h5sget_simple_extent_dims_f(space_id, dims, maxdims, ierr)
       if(ierr < 0 .or. ierr /= rank .or. dims(1) /= 1_HSIZE_T) then
          status = 13
          goto 950
       end if
    end if
    call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, val, dims, ierr)
    if(ierr /= 0) then
       status = 14
       goto 950
    end if
    status = 0
950 continue
    if(have_space) call h5sclose_f(space_id, ierr)
    if(have_attr) call h5aclose_f(attr_id, ierr)
  end subroutine hdf5_read_attr_dp_checked

  subroutine hdf5_read_attr_string(loc_id, name, val)
    implicit none
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: name
    character(len=*), intent(out) :: val
    integer(HID_T) :: attr_id, type_id
    integer(HSIZE_T), dimension(1) :: dims = (/1/)
    integer(SIZE_T) :: slen
    integer :: ierr
    slen = len(val)
    call h5aopen_f(loc_id, trim(name), attr_id, ierr)
    call h5tcopy_f(H5T_NATIVE_CHARACTER, type_id, ierr)
    call h5tset_size_f(type_id, slen, ierr)
    val = ''
    call h5aread_f(attr_id, type_id, val, dims, ierr)
    call h5aclose_f(attr_id, ierr)
    call h5tclose_f(type_id, ierr)
  end subroutine hdf5_read_attr_string

  subroutine hdf5_read_attr_string_checked(loc_id, name, val, status)
    implicit none
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: name
    character(len=*), intent(out) :: val
    integer, intent(out) :: status
    integer(HID_T) :: attr_id, space_id, type_id
    integer(HSIZE_T), dimension(1) :: dims = (/1/), maxdims
    integer(SIZE_T) :: slen
    integer :: ierr, rank
    logical :: have_attr, have_space, have_type

    status = 1
    have_attr = .false.
    have_space = .false.
    have_type = .false.
    val = ''
    call h5aopen_f(loc_id, trim(name), attr_id, ierr)
    if(ierr /= 0) then
       status = 10
       goto 960
    end if
    have_attr = .true.
    call h5aget_space_f(attr_id, space_id, ierr)
    if(ierr /= 0) then
       status = 11
       goto 960
    end if
    have_space = .true.
    call h5sget_simple_extent_ndims_f(space_id, rank, ierr)
    if(ierr /= 0 .or. (rank /= 0 .and. rank /= 1)) then
       status = 12
       goto 960
    end if
    if(rank == 1) then
       call h5sget_simple_extent_dims_f(space_id, dims, maxdims, ierr)
       if(ierr < 0 .or. ierr /= rank .or. dims(1) /= 1_HSIZE_T) then
          status = 13
          goto 960
       end if
    end if
    slen = len(val)
    call h5tcopy_f(H5T_NATIVE_CHARACTER, type_id, ierr)
    if(ierr /= 0) then
       status = 14
       goto 960
    end if
    have_type = .true.
    call h5tset_size_f(type_id, slen, ierr)
    if(ierr /= 0) then
       status = 15
       goto 960
    end if
    call h5aread_f(attr_id, type_id, val, dims, ierr)
    if(ierr /= 0) then
       status = 16
       goto 960
    end if
    status = 0
960 continue
    if(have_type) call h5tclose_f(type_id, ierr)
    if(have_space) call h5sclose_f(space_id, ierr)
    if(have_attr) call h5aclose_f(attr_id, ierr)
  end subroutine hdf5_read_attr_string_checked

  !=========================================================================
  ! 1D array attribute write/read
  !=========================================================================
  subroutine hdf5_write_attr_1d_dp(loc_id, name, arr, n)
    implicit none
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    real(dp), intent(in) :: arr(n)
    integer(HID_T) :: space_id, attr_id
    integer(HSIZE_T), dimension(1) :: dims
    integer :: ierr
    dims(1) = n
    call h5screate_simple_f(1, dims, space_id, ierr)
    call h5acreate_f(loc_id, trim(name), H5T_NATIVE_DOUBLE, space_id, attr_id, ierr)
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, arr, dims, ierr)
    call h5aclose_f(attr_id, ierr)
    call h5sclose_f(space_id, ierr)
  end subroutine hdf5_write_attr_1d_dp

  subroutine hdf5_read_attr_1d_dp(loc_id, name, arr, n)
    implicit none
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    real(dp), intent(out) :: arr(n)
    integer(HID_T) :: attr_id
    integer(HSIZE_T), dimension(1) :: dims
    integer :: ierr
    dims(1) = n
    call h5aopen_f(loc_id, trim(name), attr_id, ierr)
    call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, arr, dims, ierr)
    call h5aclose_f(attr_id, ierr)
  end subroutine hdf5_read_attr_1d_dp

  subroutine hdf5_read_attr_1d_dp_checked(loc_id, name, arr, n, status)
    implicit none
    integer(HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    real(dp), intent(out) :: arr(n)
    integer, intent(out) :: status
    integer(HID_T) :: attr_id, space_id
    integer(HSIZE_T), dimension(1) :: dims, maxdims
    integer :: ierr, rank
    logical :: have_attr, have_space

    status = 1
    have_attr = .false.
    have_space = .false.
    if(n < 1) then
       status = 2
       return
    end if
    call h5aopen_f(loc_id, trim(name), attr_id, ierr)
    if(ierr /= 0) then
       status = 10
       goto 970
    end if
    have_attr = .true.
    call h5aget_space_f(attr_id, space_id, ierr)
    if(ierr /= 0) then
       status = 11
       goto 970
    end if
    have_space = .true.
    call h5sget_simple_extent_ndims_f(space_id, rank, ierr)
    if(ierr /= 0 .or. rank /= 1) then
       status = 12
       goto 970
    end if
    dims(1) = 0_HSIZE_T
    maxdims(1) = 0_HSIZE_T
    call h5sget_simple_extent_dims_f(space_id, dims, maxdims, ierr)
    if(ierr < 0 .or. ierr /= rank .or. dims(1) /= int(n, HSIZE_T)) then
       status = 13
       goto 970
    end if
    call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, arr, dims, ierr)
    if(ierr /= 0) then
       status = 14
       goto 970
    end if
    status = 0
970 continue
    if(have_space) call h5sclose_f(space_id, ierr)
    if(have_attr) call h5aclose_f(attr_id, ierr)
  end subroutine hdf5_read_attr_1d_dp_checked

  !=========================================================================
  ! Dataset write: all ranks collectively write with hyperslab selection
  !=========================================================================
  subroutine hdf5_write_dataset_1d_dp(grp_id, name, data, nlocal, offset_global, ntotal)
    ! Write a 1D double-precision dataset with hyperslab parallel write
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: nlocal
    integer(i8b), intent(in) :: offset_global, ntotal
    real(dp), intent(in) :: data(nlocal)
    integer(HID_T) :: dspace_id, dset_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: gdims, ldims, offset, count
    integer :: ierr

    gdims(1) = ntotal
    call h5screate_simple_f(1, gdims, dspace_id, ierr)
    call h5dcreate_f(grp_id, trim(name), H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)

    ! Select hyperslab in file space
    offset(1) = offset_global
    count(1) = nlocal
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, ierr)

    ! Create memory space
    ldims(1) = nlocal
    call h5screate_simple_f(1, ldims, memspace_id, ierr)

    ! Collective write
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, ldims, ierr, &
         mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)

    call h5pclose_f(plist_id, ierr)
    call h5sclose_f(memspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(dspace_id, ierr)
  end subroutine hdf5_write_dataset_1d_dp

  subroutine hdf5_write_dataset_1d_int(grp_id, name, data, nlocal, offset_global, ntotal)
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: nlocal
    integer(i8b), intent(in) :: offset_global, ntotal
    integer, intent(in) :: data(nlocal)
    integer(HID_T) :: dspace_id, dset_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: gdims, ldims, offset, count
    integer :: ierr

    gdims(1) = ntotal
    call h5screate_simple_f(1, gdims, dspace_id, ierr)
    call h5dcreate_f(grp_id, trim(name), H5T_NATIVE_INTEGER, dspace_id, dset_id, ierr)

    offset(1) = offset_global
    count(1) = nlocal
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, ierr)

    ldims(1) = nlocal
    call h5screate_simple_f(1, ldims, memspace_id, ierr)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, ldims, ierr, &
         mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)

    call h5pclose_f(plist_id, ierr)
    call h5sclose_f(memspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(dspace_id, ierr)
  end subroutine hdf5_write_dataset_1d_int

  subroutine hdf5_write_dataset_1d_int8(grp_id, name, data, nlocal, offset_global, ntotal)
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: nlocal
    integer(i8b), intent(in) :: offset_global, ntotal
    integer(i8b), intent(in) :: data(nlocal)
    integer(HID_T) :: dspace_id, dset_id, memspace_id, plist_id, type_id
    integer(HSIZE_T), dimension(1) :: gdims, ldims, offset, count
    integer :: ierr

    gdims(1) = ntotal
    call h5screate_simple_f(1, gdims, dspace_id, ierr)
    call h5tcopy_f(H5T_NATIVE_INTEGER, type_id, ierr)
    call h5tset_size_f(type_id, int(8, SIZE_T), ierr)
    call h5dcreate_f(grp_id, trim(name), type_id, dspace_id, dset_id, ierr)

    offset(1) = offset_global
    count(1) = nlocal
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, ierr)

    ldims(1) = nlocal
    call h5screate_simple_f(1, ldims, memspace_id, ierr)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dwrite_f(dset_id, type_id, data, ldims, ierr, &
         mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)

    call h5pclose_f(plist_id, ierr)
    call h5tclose_f(type_id, ierr)
    call h5sclose_f(memspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(dspace_id, ierr)
  end subroutine hdf5_write_dataset_1d_int8

  !=========================================================================
  ! Dataset read: all ranks collectively read with hyperslab selection
  !=========================================================================
  subroutine hdf5_read_dataset_1d_dp(grp_id, name, data, nlocal, offset_global)
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: nlocal
    integer(i8b), intent(in) :: offset_global
    real(dp), intent(out) :: data(nlocal)
    integer(HID_T) :: dset_id, dspace_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: ldims, offset, count
    integer :: ierr

    call h5dopen_f(grp_id, trim(name), dset_id, ierr)
    call h5dget_space_f(dset_id, dspace_id, ierr)

    offset(1) = offset_global
    count(1) = nlocal
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, ierr)

    ldims(1) = nlocal
    call h5screate_simple_f(1, ldims, memspace_id, ierr)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, ldims, ierr, &
         mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)

    call h5pclose_f(plist_id, ierr)
    call h5sclose_f(memspace_id, ierr)
    call h5sclose_f(dspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
  end subroutine hdf5_read_dataset_1d_dp

  subroutine hdf5_read_dataset_1d_dp_checked(grp_id, name, data, nlocal, &
       offset_global, expected_total, status)
    ! Checked counterpart used by restart-critical state.  The historical
    ! helper above intentionally has no status return because many old
    ! callers predate fail-closed restart handling.  Do not use that relaxed
    ! path for serialized stellar release cursors.
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: nlocal
    integer(i8b), intent(in) :: offset_global, expected_total
    real(dp), intent(out) :: data(nlocal)
    integer, intent(out) :: status
    integer(HID_T) :: dset_id, dspace_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: dims, maxdims, ldims, offset, count
    integer :: ierr, rank
    logical :: have_dset, have_dspace, have_memspace, have_plist

    status = 1
    have_dset = .false.
    have_dspace = .false.
    have_memspace = .false.
    have_plist = .false.

    if(nlocal < 0 .or. offset_global < 0_i8b .or. expected_total < 0_i8b) then
       status = 2
       return
    end if
    if(offset_global > expected_total) then
       status = 3
       return
    end if
    if(int(nlocal, i8b) > expected_total-offset_global) then
       status = 4
       return
    end if

    call h5dopen_f(grp_id, trim(name), dset_id, ierr)
    if(ierr /= 0) then
       status = 10
       goto 900
    end if
    have_dset = .true.
    call h5dget_space_f(dset_id, dspace_id, ierr)
    if(ierr /= 0) then
       status = 11
       goto 900
    end if
    have_dspace = .true.
    call h5sget_simple_extent_ndims_f(dspace_id, rank, ierr)
    if(ierr /= 0 .or. rank /= 1) then
       status = 12
       goto 900
    end if
    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, ierr)
    ! HDF5's Fortran wrapper returns the dataspace rank on success (not zero)
    ! and -1 on failure; the dimensions themselves are written to dims.
    if(ierr < 0 .or. ierr /= rank) then
       status = 13
       goto 900
    end if
    if(dims(1) /= expected_total) then
       status = 14
       goto 900
    end if

    offset(1) = int(offset_global, HSIZE_T)
    count(1) = int(nlocal, HSIZE_T)
    if(nlocal > 0) then
       call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, ierr)
    else
       call h5sselect_none_f(dspace_id, ierr)
    end if
    if(ierr /= 0) then
       status = 20
       goto 900
    end if

    ldims(1) = int(nlocal, HSIZE_T)
    call h5screate_simple_f(1, ldims, memspace_id, ierr)
    if(ierr /= 0) then
       status = 21
       goto 900
    end if
    have_memspace = .true.
    if(nlocal == 0) then
       call h5sselect_none_f(memspace_id, ierr)
       if(ierr /= 0) then
          status = 22
          goto 900
       end if
    end if

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    if(ierr /= 0) then
       status = 23
       goto 900
    end if
    have_plist = .true.
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    if(ierr /= 0) then
       status = 24
       goto 900
    end if
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, ldims, ierr, &
         mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)
    if(ierr /= 0) then
       status = 25
       goto 900
    end if

    status = 0
900 continue
    if(have_plist) call h5pclose_f(plist_id, ierr)
    if(have_memspace) call h5sclose_f(memspace_id, ierr)
    if(have_dspace) call h5sclose_f(dspace_id, ierr)
    if(have_dset) call h5dclose_f(dset_id, ierr)
  end subroutine hdf5_read_dataset_1d_dp_checked

  subroutine hdf5_read_dataset_1d_int(grp_id, name, data, nlocal, offset_global)
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: nlocal
    integer(i8b), intent(in) :: offset_global
    integer, intent(out) :: data(nlocal)
    integer(HID_T) :: dset_id, dspace_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: ldims, offset, count
    integer :: ierr

    call h5dopen_f(grp_id, trim(name), dset_id, ierr)
    call h5dget_space_f(dset_id, dspace_id, ierr)

    offset(1) = offset_global
    count(1) = nlocal
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, ierr)

    ldims(1) = nlocal
    call h5screate_simple_f(1, ldims, memspace_id, ierr)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data, ldims, ierr, &
         mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)

    call h5pclose_f(plist_id, ierr)
    call h5sclose_f(memspace_id, ierr)
    call h5sclose_f(dspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
  end subroutine hdf5_read_dataset_1d_int

  subroutine hdf5_read_dataset_1d_int_checked(grp_id, name, data, nlocal, &
       offset_global, expected_total, status)
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: nlocal
    integer(i8b), intent(in) :: offset_global, expected_total
    integer, intent(out) :: data(nlocal)
    integer, intent(out) :: status
    integer(HID_T) :: dset_id, dspace_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: dims, maxdims, ldims, offset, count
    integer :: ierr, rank
    logical :: have_dset, have_dspace, have_memspace, have_plist

    status = 1
    have_dset = .false.
    have_dspace = .false.
    have_memspace = .false.
    have_plist = .false.

    if(nlocal < 0 .or. offset_global < 0_i8b .or. expected_total < 0_i8b) then
       status = 2
       return
    end if
    if(offset_global > expected_total) then
       status = 3
       return
    end if
    if(int(nlocal, i8b) > expected_total-offset_global) then
       status = 4
       return
    end if

    call h5dopen_f(grp_id, trim(name), dset_id, ierr)
    if(ierr /= 0) then
       status = 10
       goto 900
    end if
    have_dset = .true.
    call h5dget_space_f(dset_id, dspace_id, ierr)
    if(ierr /= 0) then
       status = 11
       goto 900
    end if
    have_dspace = .true.
    call h5sget_simple_extent_ndims_f(dspace_id, rank, ierr)
    if(ierr /= 0 .or. rank /= 1) then
       status = 12
       goto 900
    end if
    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, ierr)
    if(ierr < 0 .or. ierr /= rank) then
       status = 13
       goto 900
    end if
    if(dims(1) /= expected_total) then
       status = 14
       goto 900
    end if

    offset(1) = int(offset_global, HSIZE_T)
    count(1) = int(nlocal, HSIZE_T)
    if(nlocal > 0) then
       call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, ierr)
    else
       call h5sselect_none_f(dspace_id, ierr)
    end if
    if(ierr /= 0) then
       status = 20
       goto 900
    end if

    ldims(1) = int(nlocal, HSIZE_T)
    call h5screate_simple_f(1, ldims, memspace_id, ierr)
    if(ierr /= 0) then
       status = 21
       goto 900
    end if
    have_memspace = .true.
    if(nlocal == 0) then
       call h5sselect_none_f(memspace_id, ierr)
       if(ierr /= 0) then
          status = 22
          goto 900
       end if
    end if

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    if(ierr /= 0) then
       status = 23
       goto 900
    end if
    have_plist = .true.
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    if(ierr /= 0) then
       status = 24
       goto 900
    end if
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data, ldims, ierr, &
         mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)
    if(ierr /= 0) then
       status = 25
       goto 900
    end if

    status = 0
900 continue
    if(have_plist) call h5pclose_f(plist_id, ierr)
    if(have_memspace) call h5sclose_f(memspace_id, ierr)
    if(have_dspace) call h5sclose_f(dspace_id, ierr)
    if(have_dset) call h5dclose_f(dset_id, ierr)
  end subroutine hdf5_read_dataset_1d_int_checked

  subroutine hdf5_read_dataset_1d_int8(grp_id, name, data, nlocal, offset_global)
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: nlocal
    integer(i8b), intent(in) :: offset_global
    integer(i8b), intent(out) :: data(nlocal)
    integer(HID_T) :: dset_id, dspace_id, memspace_id, plist_id, type_id
    integer(HSIZE_T), dimension(1) :: ldims, offset, count
    integer :: ierr

    call h5dopen_f(grp_id, trim(name), dset_id, ierr)
    call h5dget_space_f(dset_id, dspace_id, ierr)

    offset(1) = offset_global
    count(1) = nlocal
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, ierr)

    ldims(1) = nlocal
    call h5screate_simple_f(1, ldims, memspace_id, ierr)

    call h5tcopy_f(H5T_NATIVE_INTEGER, type_id, ierr)
    call h5tset_size_f(type_id, int(8, SIZE_T), ierr)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dread_f(dset_id, type_id, data, ldims, ierr, &
         mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)

    call h5pclose_f(plist_id, ierr)
    call h5tclose_f(type_id, ierr)
    call h5sclose_f(memspace_id, ierr)
    call h5sclose_f(dspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
  end subroutine hdf5_read_dataset_1d_int8

  subroutine hdf5_read_dataset_1d_int8_checked(grp_id, name, data, nlocal, &
       offset_global, expected_total, status)
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: nlocal
    integer(i8b), intent(in) :: offset_global, expected_total
    integer(i8b), intent(out) :: data(nlocal)
    integer, intent(out) :: status
    integer(HID_T) :: dset_id, dspace_id, memspace_id, plist_id, type_id
    integer(HSIZE_T), dimension(1) :: dims, maxdims, ldims, offset, count
    integer :: ierr, rank
    logical :: have_dset, have_dspace, have_memspace, have_plist, have_type

    status = 1
    have_dset = .false.
    have_dspace = .false.
    have_memspace = .false.
    have_plist = .false.
    have_type = .false.

    if(nlocal < 0 .or. offset_global < 0_i8b .or. expected_total < 0_i8b) then
       status = 2
       return
    end if
    if(offset_global > expected_total) then
       status = 3
       return
    end if
    if(int(nlocal, i8b) > expected_total-offset_global) then
       status = 4
       return
    end if

    call h5dopen_f(grp_id, trim(name), dset_id, ierr)
    if(ierr /= 0) then
       status = 10
       goto 910
    end if
    have_dset = .true.
    call h5dget_space_f(dset_id, dspace_id, ierr)
    if(ierr /= 0) then
       status = 11
       goto 910
    end if
    have_dspace = .true.
    call h5sget_simple_extent_ndims_f(dspace_id, rank, ierr)
    if(ierr /= 0 .or. rank /= 1) then
       status = 12
       goto 910
    end if
    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, ierr)
    if(ierr < 0 .or. ierr /= rank) then
       status = 13
       goto 910
    end if
    if(dims(1) /= expected_total) then
       status = 14
       goto 910
    end if

    offset(1) = int(offset_global, HSIZE_T)
    count(1) = int(nlocal, HSIZE_T)
    if(nlocal > 0) then
       call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, ierr)
    else
       call h5sselect_none_f(dspace_id, ierr)
    end if
    if(ierr /= 0) then
       status = 20
       goto 910
    end if

    ldims(1) = int(nlocal, HSIZE_T)
    call h5screate_simple_f(1, ldims, memspace_id, ierr)
    if(ierr /= 0) then
       status = 21
       goto 910
    end if
    have_memspace = .true.
    if(nlocal == 0) then
       call h5sselect_none_f(memspace_id, ierr)
       if(ierr /= 0) then
          status = 22
          goto 910
       end if
    end if

    call h5tcopy_f(H5T_NATIVE_INTEGER, type_id, ierr)
    if(ierr /= 0) then
       status = 23
       goto 910
    end if
    have_type = .true.
    call h5tset_size_f(type_id, int(8, SIZE_T), ierr)
    if(ierr /= 0) then
       status = 24
       goto 910
    end if
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    if(ierr /= 0) then
       status = 25
       goto 910
    end if
    have_plist = .true.
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    if(ierr /= 0) then
       status = 26
       goto 910
    end if
    call h5dread_f(dset_id, type_id, data, ldims, ierr, &
         mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)
    if(ierr /= 0) then
       status = 27
       goto 910
    end if

    status = 0
910 continue
    if(have_plist) call h5pclose_f(plist_id, ierr)
    if(have_type) call h5tclose_f(type_id, ierr)
    if(have_memspace) call h5sclose_f(memspace_id, ierr)
    if(have_dspace) call h5sclose_f(dspace_id, ierr)
    if(have_dset) call h5dclose_f(dset_id, ierr)
  end subroutine hdf5_read_dataset_1d_int8_checked

  !=========================================================================
  ! Rank-0-only dataset write/read (for sinks, small metadata)
  !=========================================================================
  subroutine hdf5_write_dataset_serial_dp(grp_id, name, data, n, myid)
    implicit none
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n, myid
    real(dp), intent(in) :: data(n)
    integer(HID_T) :: dspace_id, dset_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: gdims, ldims, offset, count
    integer :: ierr

    gdims(1) = n
    call h5screate_simple_f(1, gdims, dspace_id, ierr)
    call h5dcreate_f(grp_id, trim(name), H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)

    ! Only rank 1 (myid==1) writes; others select none
    if(myid == 1) then
       offset(1) = 0; count(1) = n
    else
       offset(1) = 0; count(1) = 0
    end if
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, ierr)

    ldims(1) = count(1)
    call h5screate_simple_f(1, ldims, memspace_id, ierr)
    if(myid /= 1) call h5sselect_none_f(memspace_id, ierr)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, ldims, ierr, &
         mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)

    call h5pclose_f(plist_id, ierr)
    call h5sclose_f(memspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(dspace_id, ierr)
  end subroutine hdf5_write_dataset_serial_dp

  subroutine hdf5_write_dataset_serial_int(grp_id, name, data, n, myid)
    implicit none
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n, myid
    integer, intent(in) :: data(n)
    integer(HID_T) :: dspace_id, dset_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: gdims, ldims, offset, count
    integer :: ierr

    gdims(1) = n
    call h5screate_simple_f(1, gdims, dspace_id, ierr)
    call h5dcreate_f(grp_id, trim(name), H5T_NATIVE_INTEGER, dspace_id, dset_id, ierr)

    if(myid == 1) then
       offset(1) = 0; count(1) = n
    else
       offset(1) = 0; count(1) = 0
    end if
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, ierr)

    ldims(1) = count(1)
    call h5screate_simple_f(1, ldims, memspace_id, ierr)
    if(myid /= 1) call h5sselect_none_f(memspace_id, ierr)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, ldims, ierr, &
         mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)

    call h5pclose_f(plist_id, ierr)
    call h5sclose_f(memspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(dspace_id, ierr)
  end subroutine hdf5_write_dataset_serial_int

  subroutine hdf5_read_dataset_all_dp(grp_id, name, data, n)
    ! Replicate dataset to ALL ranks via rank-0 streamed read + MPI_Bcast.
    ! Was: every rank collectively re-read the full dataset → O(N²) GPFS traffic.
    ! Now: rank 0 reads n/nproc-sized chunks sequentially, broadcasting each chunk
    !      so the GPFS volume sees ~1× n bytes total instead of nproc × n bytes.
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    real(dp), intent(out) :: data(n)
    integer(HID_T) :: dset_id, dspace_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: ldims, fooffset, count
    integer :: ierr, myrank, nproc, ichunk, nchunk, chunk_size, this_chunk, istart

    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

    if (n <= 0) return

    if (n < nproc) then
       chunk_size = n
    else
       chunk_size = (n + nproc - 1) / nproc
    end if
    nchunk = (n + chunk_size - 1) / chunk_size

    call h5dopen_f(grp_id, trim(name), dset_id, ierr)
    if (myrank == 0) call h5dget_space_f(dset_id, dspace_id, ierr)
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierr)

    do ichunk = 1, nchunk
       istart = (ichunk-1) * chunk_size + 1
       this_chunk = min(chunk_size, n - istart + 1)
       if (this_chunk <= 0) exit
       if (myrank == 0) then
          fooffset(1) = int(istart - 1, HSIZE_T)
          count(1) = int(this_chunk, HSIZE_T)
          call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, fooffset, count, ierr)
          ldims(1) = int(this_chunk, HSIZE_T)
          call h5screate_simple_f(1, ldims, memspace_id, ierr)
          call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data(istart:istart+this_chunk-1), &
               ldims, ierr, mem_space_id=memspace_id, file_space_id=dspace_id, &
               xfer_prp=plist_id)
          call h5sclose_f(memspace_id, ierr)
       end if
       call MPI_Bcast(data(istart), this_chunk, MPI_DOUBLE_PRECISION, 0, &
            MPI_COMM_WORLD, ierr)
    end do

    call h5pclose_f(plist_id, ierr)
    if (myrank == 0) call h5sclose_f(dspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
  end subroutine hdf5_read_dataset_all_dp

  subroutine hdf5_read_dataset_all_int(grp_id, name, data, n)
    ! Integer counterpart of hdf5_read_dataset_all_dp (rank-0 chunk read + Bcast).
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    integer, intent(out) :: data(n)
    integer(HID_T) :: dset_id, dspace_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: ldims, fooffset, count
    integer :: ierr, myrank, nproc, ichunk, nchunk, chunk_size, this_chunk, istart

    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

    if (n <= 0) return

    if (n < nproc) then
       chunk_size = n
    else
       chunk_size = (n + nproc - 1) / nproc
    end if
    nchunk = (n + chunk_size - 1) / chunk_size

    call h5dopen_f(grp_id, trim(name), dset_id, ierr)
    if (myrank == 0) call h5dget_space_f(dset_id, dspace_id, ierr)
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierr)

    do ichunk = 1, nchunk
       istart = (ichunk-1) * chunk_size + 1
       this_chunk = min(chunk_size, n - istart + 1)
       if (this_chunk <= 0) exit
       if (myrank == 0) then
          fooffset(1) = int(istart - 1, HSIZE_T)
          count(1) = int(this_chunk, HSIZE_T)
          call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, fooffset, count, ierr)
          ldims(1) = int(this_chunk, HSIZE_T)
          call h5screate_simple_f(1, ldims, memspace_id, ierr)
          call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data(istart:istart+this_chunk-1), &
               ldims, ierr, mem_space_id=memspace_id, file_space_id=dspace_id, &
               xfer_prp=plist_id)
          call h5sclose_f(memspace_id, ierr)
       end if
       call MPI_Bcast(data(istart), this_chunk, MPI_INTEGER, 0, &
            MPI_COMM_WORLD, ierr)
    end do

    call h5pclose_f(plist_id, ierr)
    if (myrank == 0) call h5sclose_f(dspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
  end subroutine hdf5_read_dataset_all_int

  subroutine hdf5_read_dataset_all_int_checked(grp_id, name, data, n, status)
    ! Checked rank-0 read + broadcast for restart metadata.  Unlike the
    ! historical helper above, every HDF5/MPI stage is reduced to all ranks
    ! before any rank proceeds, so a missing or malformed count array cannot
    ! leave the other ranks waiting in a broadcast with untrusted offsets.
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    integer, intent(out) :: data(n)
    integer, intent(out) :: status
    integer(HID_T) :: dset_id, dspace_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: dims, maxdims, ldims
    integer(HSIZE_T), dimension(1) :: offset, count
    integer :: ierr, myrank, local_status, global_status, rank
    logical :: have_dset, have_dspace, have_memspace, have_plist

    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
    status = 0
    local_status = 0
    have_dset = .false.
    have_dspace = .false.
    have_memspace = .false.
    have_plist = .false.

    if(n <= 0) local_status = 2
    if(local_status == 0) then
       call h5dopen_f(grp_id, trim(name), dset_id, ierr)
       if(ierr /= 0) then
          local_status = 10
       else
          have_dset = .true.
       end if
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    if(global_status /= 0) then
       status = global_status
       goto 920
    end if

    if(myrank == 0) then
       call h5dget_space_f(dset_id, dspace_id, ierr)
       if(ierr /= 0) then
          local_status = 11
       else
          have_dspace = .true.
       end if
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    if(global_status /= 0) then
       status = global_status
       goto 920
    end if

    if(myrank == 0) then
       call h5sget_simple_extent_ndims_f(dspace_id, rank, ierr)
       if(ierr /= 0 .or. rank /= 1) then
          local_status = 12
       else
          call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, ierr)
          if(ierr < 0 .or. ierr /= rank) then
             local_status = 13
          else if(dims(1) /= int(n, HSIZE_T)) then
             local_status = 14
          end if
       end if
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    if(global_status /= 0) then
       status = global_status
       goto 920
    end if

    if(myrank == 0) then
       ldims(1) = int(n, HSIZE_T)
       offset(1) = 0_HSIZE_T
       count(1) = int(n, HSIZE_T)
       call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, ierr)
       if(ierr /= 0) then
          local_status = 20
       else
          call h5screate_simple_f(1, ldims, memspace_id, ierr)
          if(ierr /= 0) then
             local_status = 21
          else
             have_memspace = .true.
          end if
       end if
    end if
    if(myrank == 0 .and. local_status == 0) then
       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
       if(ierr /= 0) then
          local_status = 22
       else
          have_plist = .true.
          call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierr)
          if(ierr /= 0) local_status = 23
       end if
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    if(global_status /= 0) then
       status = global_status
       goto 920
    end if

    if(myrank == 0) then
       call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data, ldims, ierr, &
            mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)
       if(ierr /= 0) local_status = 24
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    if(global_status /= 0) then
       status = global_status
       goto 920
    end if
    call MPI_Bcast(data, n, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if(ierr /= 0) status = 25

920 continue
    if(myrank == 0) then
       if(have_plist) call h5pclose_f(plist_id, ierr)
       if(have_memspace) call h5sclose_f(memspace_id, ierr)
       if(have_dspace) call h5sclose_f(dspace_id, ierr)
    end if
    if(have_dset) call h5dclose_f(dset_id, ierr)
  end subroutine hdf5_read_dataset_all_int_checked

  subroutine hdf5_read_dataset_all_dp_checked(grp_id, name, data, n, status)
    ! Checked rank-0 read + broadcast for restart metadata stored as doubles.
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    real(dp), intent(out) :: data(n)
    integer, intent(out) :: status
    integer(HID_T) :: dset_id, dspace_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: dims, maxdims, ldims
    integer(HSIZE_T), dimension(1) :: offset, count
    integer :: ierr, myrank, local_status, global_status, rank
    logical :: have_dset, have_dspace, have_memspace, have_plist

    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
    status = 0
    local_status = 0
    have_dset = .false.
    have_dspace = .false.
    have_memspace = .false.
    have_plist = .false.

    if(n <= 0) local_status = 2
    if(local_status == 0) then
       call h5dopen_f(grp_id, trim(name), dset_id, ierr)
       if(ierr /= 0) then
          local_status = 10
       else
          have_dset = .true.
       end if
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    if(global_status /= 0) then
       status = global_status
       goto 950
    end if

    if(myrank == 0) then
       call h5dget_space_f(dset_id, dspace_id, ierr)
       if(ierr /= 0) then
          local_status = 11
       else
          have_dspace = .true.
       end if
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    if(global_status /= 0) then
       status = global_status
       goto 950
    end if

    if(myrank == 0) then
       call h5sget_simple_extent_ndims_f(dspace_id, rank, ierr)
       if(ierr /= 0 .or. rank /= 1) then
          local_status = 12
       else
          call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, ierr)
          if(ierr < 0 .or. ierr /= rank) then
             local_status = 13
          else if(dims(1) /= int(n, HSIZE_T)) then
             local_status = 14
          end if
       end if
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    if(global_status /= 0) then
       status = global_status
       goto 950
    end if

    if(myrank == 0) then
       ldims(1) = int(n, HSIZE_T)
       offset(1) = 0_HSIZE_T
       count(1) = int(n, HSIZE_T)
       call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, ierr)
       if(ierr /= 0) then
          local_status = 20
       else
          call h5screate_simple_f(1, ldims, memspace_id, ierr)
          if(ierr /= 0) then
             local_status = 21
          else
             have_memspace = .true.
          end if
       end if
    end if
    if(myrank == 0 .and. local_status == 0) then
       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
       if(ierr /= 0) then
          local_status = 22
       else
          have_plist = .true.
          call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierr)
          if(ierr /= 0) local_status = 23
       end if
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    if(global_status /= 0) then
       status = global_status
       goto 950
    end if

    if(myrank == 0) then
       call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, ldims, ierr, &
            mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)
       if(ierr /= 0) local_status = 24
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    if(global_status /= 0) then
       status = global_status
       goto 950
    end if
    call MPI_Bcast(data, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if(ierr /= 0) status = 25

950 continue
    if(myrank == 0) then
       if(have_plist) call h5pclose_f(plist_id, ierr)
       if(have_memspace) call h5sclose_f(memspace_id, ierr)
       if(have_dspace) call h5sclose_f(dspace_id, ierr)
    end if
    if(have_dset) call h5dclose_f(dset_id, ierr)
  end subroutine hdf5_read_dataset_all_dp_checked

  subroutine hdf5_read_dataset_chunk_dp(grp_id, name, data, n, offset_global)
    ! Replicate a chunk [offset_global+1 .. offset_global+n] of dataset to all
    ! ranks via rank-0 read + MPI_Bcast. Used by streaming varcpu restore to
    ! cap per-rank transient buffers (vs. full-N allocation in _all variants).
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    integer(i8b), intent(in) :: offset_global
    real(dp), intent(out) :: data(n)
    integer(HID_T) :: dset_id, dspace_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: ldims, fooffset, count
    integer :: ierr, myrank

    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
    if (n <= 0) return

    call h5dopen_f(grp_id, trim(name), dset_id, ierr)
    if (myrank == 0) then
       call h5dget_space_f(dset_id, dspace_id, ierr)
       fooffset(1) = int(offset_global, HSIZE_T)
       count(1) = int(n, HSIZE_T)
       call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, fooffset, count, ierr)
       ldims(1) = int(n, HSIZE_T)
       call h5screate_simple_f(1, ldims, memspace_id, ierr)
       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierr)
       call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, ldims, ierr, &
            mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)
       call h5pclose_f(plist_id, ierr)
       call h5sclose_f(memspace_id, ierr)
       call h5sclose_f(dspace_id, ierr)
    end if
    call h5dclose_f(dset_id, ierr)
    call MPI_Bcast(data, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  end subroutine hdf5_read_dataset_chunk_dp

  subroutine hdf5_read_dataset_chunk_int(grp_id, name, data, n, offset_global)
    ! Integer counterpart of hdf5_read_dataset_chunk_dp.
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    integer(i8b), intent(in) :: offset_global
    integer, intent(out) :: data(n)
    integer(HID_T) :: dset_id, dspace_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: ldims, fooffset, count
    integer :: ierr, myrank

    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
    if (n <= 0) return

    call h5dopen_f(grp_id, trim(name), dset_id, ierr)
    if (myrank == 0) then
       call h5dget_space_f(dset_id, dspace_id, ierr)
       fooffset(1) = int(offset_global, HSIZE_T)
       count(1) = int(n, HSIZE_T)
       call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, fooffset, count, ierr)
       ldims(1) = int(n, HSIZE_T)
       call h5screate_simple_f(1, ldims, memspace_id, ierr)
       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierr)
       call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data, ldims, ierr, &
            mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)
       call h5pclose_f(plist_id, ierr)
       call h5sclose_f(memspace_id, ierr)
       call h5sclose_f(dspace_id, ierr)
    end if
    call h5dclose_f(dset_id, ierr)
    call MPI_Bcast(data, n, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  end subroutine hdf5_read_dataset_chunk_int

  subroutine hdf5_read_dataset_collective_dp(grp_id, name, data, n, offset_global)
    ! Tier 2 parallel restore primitive: each rank reads its own disjoint
    ! hyperslab [offset_global+1 .. offset_global+n] in a single collective
    ! H5Dread. No post-read Bcast: each rank gets only its slice.
    ! Ranks with n==0 must still participate (H5Sselect_none) — collective.
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    integer(i8b), intent(in) :: offset_global
    real(dp), intent(out) :: data(*)
    integer(HID_T) :: dset_id, dspace_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: ldims, fooffset, count
    integer :: ierr

    call h5dopen_f(grp_id, trim(name), dset_id, ierr)
    call h5dget_space_f(dset_id, dspace_id, ierr)
    if (n > 0) then
       fooffset(1) = int(offset_global, HSIZE_T)
       count(1) = int(n, HSIZE_T)
       call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, fooffset, count, ierr)
       ldims(1) = int(n, HSIZE_T)
       call h5screate_simple_f(1, ldims, memspace_id, ierr)
    else
       call h5sselect_none_f(dspace_id, ierr)
       ldims(1) = 0_HSIZE_T
       call h5screate_simple_f(1, ldims, memspace_id, ierr)
       call h5sselect_none_f(memspace_id, ierr)
    end if
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, ldims, ierr, &
         mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)
    call h5pclose_f(plist_id, ierr)
    call h5sclose_f(memspace_id, ierr)
    call h5sclose_f(dspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
  end subroutine hdf5_read_dataset_collective_dp

  subroutine hdf5_read_dataset_collective_int(grp_id, name, data, n, offset_global)
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    integer(i8b), intent(in) :: offset_global
    integer, intent(out) :: data(*)
    integer(HID_T) :: dset_id, dspace_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: ldims, fooffset, count
    integer :: ierr

    call h5dopen_f(grp_id, trim(name), dset_id, ierr)
    call h5dget_space_f(dset_id, dspace_id, ierr)
    if (n > 0) then
       fooffset(1) = int(offset_global, HSIZE_T)
       count(1) = int(n, HSIZE_T)
       call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, fooffset, count, ierr)
       ldims(1) = int(n, HSIZE_T)
       call h5screate_simple_f(1, ldims, memspace_id, ierr)
    else
       call h5sselect_none_f(dspace_id, ierr)
       ldims(1) = 0_HSIZE_T
       call h5screate_simple_f(1, ldims, memspace_id, ierr)
       call h5sselect_none_f(memspace_id, ierr)
    end if
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data, ldims, ierr, &
         mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)
    call h5pclose_f(plist_id, ierr)
    call h5sclose_f(memspace_id, ierr)
    call h5sclose_f(dspace_id, ierr)
    call h5dclose_f(dset_id, ierr)
  end subroutine hdf5_read_dataset_collective_int

  subroutine hdf5_read_dataset_collective_dp_checked(grp_id, name, data, n, &
       offset_global, expected_total, status)
    ! Fail-closed counterpart of the legacy collective reader.  Every rank
    ! validates the one-dimensional dataset and its own hyperslab before the
    ! collective H5Dread.  This is used for AMR payloads whose offsets are
    ! derived from checkpoint count arrays; a malformed extent must therefore
    ! stop all ranks before any untrusted selection is issued.
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    integer(i8b), intent(in) :: offset_global, expected_total
    real(dp), intent(out) :: data(*)
    integer, intent(out) :: status
    integer(HID_T) :: dset_id, dspace_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: dims, maxdims, ldims, offset, count
    integer :: ierr, rank, local_status, global_status
    logical :: have_dset, have_dspace, have_memspace, have_plist

    status = 0
    local_status = 0
    have_dset = .false.
    have_dspace = .false.
    have_memspace = .false.
    have_plist = .false.

    if(n < 0 .or. offset_global < 0_i8b .or. expected_total < 0_i8b) then
       local_status = 2
    else if(offset_global > expected_total) then
       local_status = 3
    else if(int(n, i8b) > expected_total-offset_global) then
       local_status = 4
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    if(global_status /= 0) then
       status = global_status
       goto 1010
    end if

    call h5dopen_f(grp_id, trim(name), dset_id, ierr)
    if(ierr /= 0) then
       local_status = 10
    else
       have_dset = .true.
       local_status = 0
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    if(global_status /= 0) then
       status = global_status
       goto 1010
    end if

    call h5dget_space_f(dset_id, dspace_id, ierr)
    if(ierr /= 0) then
       local_status = 11
    else
       have_dspace = .true.
       local_status = 0
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    if(global_status /= 0) then
       status = global_status
       goto 1010
    end if

    call h5sget_simple_extent_ndims_f(dspace_id, rank, ierr)
    if(ierr /= 0 .or. rank /= 1) then
       local_status = 12
    else
       call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, ierr)
       if(ierr < 0 .or. ierr /= rank) then
          local_status = 13
       else if(dims(1) /= int(expected_total, HSIZE_T)) then
          local_status = 14
       else
          local_status = 0
       end if
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    if(global_status /= 0) then
       status = global_status
       goto 1010
    end if

    if(n > 0) then
       offset(1) = int(offset_global, HSIZE_T)
       count(1) = int(n, HSIZE_T)
       call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, ierr)
    else
       call h5sselect_none_f(dspace_id, ierr)
    end if
    if(ierr /= 0) then
       local_status = 20
    else
       local_status = 0
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    if(global_status /= 0) then
       status = global_status
       goto 1010
    end if

    ldims(1) = int(n, HSIZE_T)
    call h5screate_simple_f(1, ldims, memspace_id, ierr)
    if(ierr /= 0) then
       local_status = 21
    else
       have_memspace = .true.
       local_status = 0
    end if
    if(n == 0 .and. local_status == 0) then
       call h5sselect_none_f(memspace_id, ierr)
       if(ierr /= 0) local_status = 22
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    if(global_status /= 0) then
       status = global_status
       goto 1010
    end if

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    if(ierr /= 0) then
       local_status = 23
    else
       have_plist = .true.
       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
       if(ierr /= 0) then
          local_status = 24
       else
          local_status = 0
       end if
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    if(global_status /= 0) then
       status = global_status
       goto 1010
    end if

    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, ldims, ierr, &
         mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)
    if(ierr /= 0) then
       local_status = 25
    else
       local_status = 0
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    status = global_status

1010 continue
    if(have_plist) call h5pclose_f(plist_id, ierr)
    if(have_memspace) call h5sclose_f(memspace_id, ierr)
    if(have_dspace) call h5sclose_f(dspace_id, ierr)
    if(have_dset) call h5dclose_f(dset_id, ierr)
  end subroutine hdf5_read_dataset_collective_dp_checked

  subroutine hdf5_read_dataset_collective_int_checked(grp_id, name, data, n, &
       offset_global, expected_total, status)
    ! Integer version of hdf5_read_dataset_collective_dp_checked.
    implicit none
    include 'mpif.h'
    integer(HID_T), intent(in) :: grp_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    integer(i8b), intent(in) :: offset_global, expected_total
    integer, intent(out) :: data(*)
    integer, intent(out) :: status
    integer(HID_T) :: dset_id, dspace_id, memspace_id, plist_id
    integer(HSIZE_T), dimension(1) :: dims, maxdims, ldims, offset, count
    integer :: ierr, rank, local_status, global_status
    logical :: have_dset, have_dspace, have_memspace, have_plist

    status = 0
    local_status = 0
    have_dset = .false.
    have_dspace = .false.
    have_memspace = .false.
    have_plist = .false.

    if(n < 0 .or. offset_global < 0_i8b .or. expected_total < 0_i8b) then
       local_status = 2
    else if(offset_global > expected_total) then
       local_status = 3
    else if(int(n, i8b) > expected_total-offset_global) then
       local_status = 4
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    if(global_status /= 0) then
       status = global_status
       goto 1110
    end if

    call h5dopen_f(grp_id, trim(name), dset_id, ierr)
    if(ierr /= 0) then
       local_status = 10
    else
       have_dset = .true.
       local_status = 0
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    if(global_status /= 0) then
       status = global_status
       goto 1110
    end if

    call h5dget_space_f(dset_id, dspace_id, ierr)
    if(ierr /= 0) then
       local_status = 11
    else
       have_dspace = .true.
       local_status = 0
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    if(global_status /= 0) then
       status = global_status
       goto 1110
    end if

    call h5sget_simple_extent_ndims_f(dspace_id, rank, ierr)
    if(ierr /= 0 .or. rank /= 1) then
       local_status = 12
    else
       call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, ierr)
       if(ierr < 0 .or. ierr /= rank) then
          local_status = 13
       else if(dims(1) /= int(expected_total, HSIZE_T)) then
          local_status = 14
       else
          local_status = 0
       end if
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    if(global_status /= 0) then
       status = global_status
       goto 1110
    end if

    if(n > 0) then
       offset(1) = int(offset_global, HSIZE_T)
       count(1) = int(n, HSIZE_T)
       call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, ierr)
    else
       call h5sselect_none_f(dspace_id, ierr)
    end if
    if(ierr /= 0) then
       local_status = 20
    else
       local_status = 0
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    if(global_status /= 0) then
       status = global_status
       goto 1110
    end if

    ldims(1) = int(n, HSIZE_T)
    call h5screate_simple_f(1, ldims, memspace_id, ierr)
    if(ierr /= 0) then
       local_status = 21
    else
       have_memspace = .true.
       local_status = 0
    end if
    if(n == 0 .and. local_status == 0) then
       call h5sselect_none_f(memspace_id, ierr)
       if(ierr /= 0) local_status = 22
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    if(global_status /= 0) then
       status = global_status
       goto 1110
    end if

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    if(ierr /= 0) then
       local_status = 23
    else
       have_plist = .true.
       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
       if(ierr /= 0) then
          local_status = 24
       else
          local_status = 0
       end if
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    if(global_status /= 0) then
       status = global_status
       goto 1110
    end if

    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data, ldims, ierr, &
         mem_space_id=memspace_id, file_space_id=dspace_id, xfer_prp=plist_id)
    if(ierr /= 0) then
       local_status = 25
    else
       local_status = 0
    end if
    call MPI_Allreduce(local_status, global_status, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    status = global_status

1110 continue
    if(have_plist) call h5pclose_f(plist_id, ierr)
    if(have_memspace) call h5sclose_f(memspace_id, ierr)
    if(have_dspace) call h5sclose_f(dspace_id, ierr)
    if(have_dset) call h5dclose_f(dset_id, ierr)
  end subroutine hdf5_read_dataset_collective_int_checked

  !=========================================================================
  ! Phase D: file pool (LRU cache of open HDF5 file handles)
  !
  ! Rationale: at exascale, restore may need to read from O(10^3-10^5)
  ! per-writer files. Keeping all of them open exhausts file descriptors;
  ! open/close per access wastes O(ms) HDF5 setup per dataset. An LRU pool
  ! of size K caps open descriptors at K while amortising open overhead
  ! over many dataset accesses to the same file.
  !
  ! Implementation: parallel array of (name, file_id, last_use_tick).
  ! get() does linear search (K small, typically 8-64) → returns cached
  ! file_id on hit; on miss opens with parallel MPI-IO, evicting the
  ! lowest-tick slot if full. Tick is monotonically incremented on every
  ! get() to define LRU order.
  !=========================================================================
  subroutine file_pool_init(max_open)
    implicit none
    integer, intent(in), optional :: max_open
    integer :: cap, ierr
    if(present(max_open)) then
       cap = max(1, min(max_open, FP_MAX_LIMIT))
    else
       cap = FP_MAX_DEFAULT
    end if
    ! Idempotent: re-init flushes whatever was there
    if(allocated(fp_file_id)) call file_pool_close_all()
    fp_cap = cap
    allocate(fp_name(cap), fp_file_id(cap), fp_lru(cap))
    fp_name = ''
    fp_file_id = -1_HID_T
    fp_lru = 0_i8b
    fp_count = 0
    fp_tick = 0_i8b
    fp_hits = 0_i8b
    fp_misses = 0_i8b
    fp_evictions = 0_i8b
    call h5open_f(ierr)
  end subroutine file_pool_init

  subroutine file_pool_get(filename, comm, file_id)
    implicit none
    include 'mpif.h'
    character(len=*), intent(in)  :: filename
    integer,          intent(in)  :: comm
    integer(HID_T),   intent(out) :: file_id
    integer :: i, victim, ierr
    integer(HID_T) :: plist_id

    ! Auto-init if caller forgot (uses default capacity)
    if(.not. allocated(fp_file_id)) call file_pool_init()

    fp_tick = fp_tick + 1_i8b

    ! Cache lookup
    do i = 1, fp_count
       if(trim(fp_name(i)) == trim(filename)) then
          file_id = fp_file_id(i)
          fp_lru(i) = fp_tick
          fp_hits = fp_hits + 1_i8b
          return
       end if
    end do

    fp_misses = fp_misses + 1_i8b

    ! Evict LRU slot if pool is full
    if(fp_count >= fp_cap) then
       victim = 1
       do i = 2, fp_count
          if(fp_lru(i) < fp_lru(victim)) victim = i
       end do
       call h5fclose_f(fp_file_id(victim), ierr)
       fp_evictions = fp_evictions + 1_i8b
       ! Compact: move the last live slot into the victim's hole
       if(victim /= fp_count) then
          fp_name(victim)    = fp_name(fp_count)
          fp_file_id(victim) = fp_file_id(fp_count)
          fp_lru(victim)     = fp_lru(fp_count)
       end if
       fp_name(fp_count)    = ''
       fp_file_id(fp_count) = -1_HID_T
       fp_lru(fp_count)     = 0_i8b
       fp_count = fp_count - 1
    end if

    ! Open new
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
    call h5pset_fapl_mpio_f(plist_id, comm, MPI_INFO_NULL, ierr)
    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, ierr, &
         access_prp=plist_id)
    call h5pclose_f(plist_id, ierr)

    fp_count = fp_count + 1
    fp_name(fp_count)    = trim(filename)
    fp_file_id(fp_count) = file_id
    fp_lru(fp_count)     = fp_tick
  end subroutine file_pool_get

  subroutine file_pool_close_all()
    implicit none
    integer :: i, ierr
    if(.not. allocated(fp_file_id)) return
    do i = 1, fp_count
       call h5fclose_f(fp_file_id(i), ierr)
    end do
    deallocate(fp_name, fp_file_id, fp_lru)
    fp_cap = 0
    fp_count = 0
    fp_tick = 0_i8b
  end subroutine file_pool_close_all

  subroutine file_pool_stats(n_open, n_hits, n_misses, n_evictions, cap)
    implicit none
    integer,      intent(out) :: n_open, cap
    integer(i8b), intent(out) :: n_hits, n_misses, n_evictions
    n_open      = fp_count
    cap         = fp_cap
    n_hits      = fp_hits
    n_misses    = fp_misses
    n_evictions = fp_evictions
  end subroutine file_pool_stats

  ! Convenience wrapper: open via pool AND publish handle into module's
  ! hdf5_file_id global so callers using the single-file API keep working.
  ! Pair with hdf5_release_pooled() (does NOT close — pool owns the handle).
  subroutine hdf5_open_pooled(filename, comm)
    implicit none
    character(len=*), intent(in) :: filename
    integer,          intent(in) :: comm
    integer(HID_T) :: fid
    call file_pool_get(filename, comm, fid)
    hdf5_file_id = fid
  end subroutine hdf5_open_pooled

  subroutine hdf5_release_pooled()
    implicit none
    hdf5_file_id = -1_HID_T
  end subroutine hdf5_release_pooled

  subroutine hdf5_suppress_errors()
    implicit none
    integer :: ierr
    call h5eset_auto_f(0, ierr)
  end subroutine hdf5_suppress_errors

  subroutine hdf5_restore_errors()
    implicit none
    integer :: ierr
    call h5eset_auto_f(1, ierr)
  end subroutine hdf5_restore_errors

end module ramses_hdf5_io
#endif
