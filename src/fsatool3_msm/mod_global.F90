module mod_global
  use util
  implicit none
  integer, parameter :: max_str_len = 256
  integer,parameter :: maxnlevel=20
  real*8 :: kelvinarray(maxnlevel)
  integer :: numFileLevel(maxnlevel)
  integer :: ndim,nsnap,procid,procnum,ntraj,nss,ncv, trajtype,nlevel, ncluster, clustercycle
  integer :: cutnumber, lagstep, nstate
  integer, dimension(:), allocatable :: snap2cluster, levelindex, ifcutcluster, trajIndex
  real*8, dimension(:,:), pointer :: traj,cvs
  real*8, dimension(:), allocatable :: refcoor,snappot
  real*8 :: lagtime,kelvin
  real*8,parameter :: gasconst=1.9872065d-3
  logical :: ifreadcoor, ifOnlyRunClustering
  character(max_str_len),save :: resultdir, tpmmethod, clustermethod, lumpingmethod, inputfile
  integer, dimension(maxnlevel), target :: startstate, endstate
contains

subroutine mod_global_mpiinit()
  use mpi
  integer :: ierr

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,procID,ierr)
  call mpi_comm_size(mpi_comm_world,procNum,ierr)
end subroutine

subroutine mod_global_readtraj()
  use mpi
  use mpi_shared_memory
  use netcdf_func
  integer, parameter :: maxnfile=30
  character(max_str_len) :: cvsfiles(maxnfile),coorfiles(maxnfile)
  integer :: numFile
  integer, allocatable :: numFramePerFile(:)
  integer :: numFrame, frameStride, natom
  integer :: iofile, ierr, shared_cv_win, shared_traj_win
  namelist /trajs/ ncv, framestride, kelvinarray, trajtype, cvsfiles, coorfiles, numFileLevel

  cvsfiles=""; framestride = 1; kelvinarray=0

  if ( procid == 0 ) then
    call loginfo("Read Collective Variables File")
    call getfreeunit(iofile)
    open(unit=iofile, file=inputfile, action="read")
    read(iofile, nml=trajs, iostat=ierr)
    if ( ierr < 0 ) call errormsg("error in reading trajs namelist")
    close(iofile)
    call get_file_number(cvsfiles, numFile)
    allocate(numFramePerFile(numFile))
    call getTrajFrame(cvsfiles, numFile, frameStride, 1, numFramePerFile, nsnap)
  endif

  call mpi_bcast(trajtype, 1, mpi_integer, 0, mpi_comm_world, ierr)
  call mpi_bcast(numFile, 1, mpi_integer,0,mpi_comm_world,ierr)
  call mpi_bcast(cvsfiles, numFile*max_str_len, mpi_char, 0, mpi_comm_world, ierr)
  if(trajtype /= 0) call mpi_bcast(coorFiles, numFile*max_str_len, mpi_char, 0, mpi_comm_world, ierr)
  call mpi_bcast(ncv, 1, mpi_int, 0, mpi_comm_world, ierr)
  call mpi_bcast(framestride,1,mpi_integer,0,mpi_comm_world,ierr)
  call mpi_bcast(nsnap,1,mpi_integer,0,mpi_comm_world,ierr)

  allocate(trajIndex(nsnap))

  call init_shared_memory(mpi_comm_world)
  call allocate_shared_memory_for_2dim_array([ncv, nsnap], cvs, shared_cv_win)

  if (trajtype /= 0) then
    call allocate_shared_memory_for_2dim_array([ndim, nsnap], traj, shared_traj_win)
  endif

  if(procid == 0) then
    call find_level_from_kelvinarray()
    write(*, "(a, I8)") "numer of files         = ", numFile
    write(*, "(a, I8)") "number of cv           = ", ncv
    write(*, "(a, I8)") "dimension of traj      = ", ndim
    write(*, "(a, I8)") "Trajectory type        = ", trajtype
    write(*, "(a, I8)") "number of total Frames = ", nsnap
    write(*, "(a, I8)") "Frame stride           = ", frameStride
    write(*, "(a, I8)") "Number of temperature  = ", nlevel
    write(*, "(a, 10f8.3)") "Kelvin array           = ", kelvinarray(1:nlevel)
  endif

  if(shared_id == 0) then
    call read_cv_file(cvsfiles, numFile, ncv, cvs, frameStride, trajIndex)
    if(trajtype == 1) then
      call read_pdb_file(coorFiles, numFile, ndim, traj, frameStride)
    else if(trajtype == 2) then
      call read_NCcoorfile_into_coor(coorFiles, numFile, ndim, traj, 1)
    endif
  endif
  call mpi_win_fence(0, shared_cv_win, ierr)
  if(trajtype /= 0) call mpi_win_fence(0, shared_traj_win, ierr)
  if(procid==0) call loginfo()
end subroutine

subroutine read_pdb_file(files, numFile, ndim, trajs, frameStride)
  character(max_str_len), intent(in) :: files(:)
  integer, intent(in) :: numFile, ndim, frameStride
  real*8, intent(out) :: trajs(:, :)

  integer :: snapIndex = 0
  integer :: frameLine, iofile
  integer :: numFrames(numFile)
  integer :: i, j

  frameLine = ceiling( dble(ndim) / 10)

  do i = 1, numFile
    call getfreeunit(iofile)
    open(unit = iofile, file=files(i), action="read")
    read(iofile, *)
    do  j = 1, frameLine * (frameStride - 1) ! skip using frameStride
      read(iofile, *)
    enddo
    snapIndex = snapIndex + 1
    if(snapIndex > nsnap) exit
    read(iofile, "(10f8.3)") trajs(:, snapIndex)
  enddo
end subroutine


subroutine read_cv_file(files, numFile, ncv, cvs, frameStride, trajIndex)
  character(max_str_len), intent(in) :: files(:)
  integer, intent(in) :: numFile, ncv, frameStride
  integer, intent(out) :: trajIndex(:)
  real*8, intent(out) :: cvs(:, :)

  real*8 ::  cv_temp(ncv)
  integer :: snapIndex = 1
  integer :: frameLine, iofile
  integer :: numFrames(numFile)
  integer :: i, j, k, ierr

  do i = 1, numFile
    call getfreeunit(iofile)
    open(unit = iofile, file=files(i), action="read")
    outer: do 
      do  j = 1, (frameStride - 1) ! skip using frameStride
        read(iofile, *, iostat=ierr)
        if(ierr < 0) exit outer
      enddo
      read(iofile, *, iostat=ierr) trajIndex(snapIndex), cvs(:, snapIndex)
      if (ierr < 0) exit
      snapIndex = snapIndex + 1
    enddo outer
  enddo
end subroutine

subroutine getTrajFrame(files, numFile, frameStride, frameLine, numFrames, nsnap)
  character(max_str_len), intent(in) :: files(:)
  integer, intent(in) :: numFile, frameStride, frameLine
  integer, intent(out) :: numFrames(:)

  logical :: ifExist
  integer :: i, j, ierr, iofile
  integer :: nsnap

  nsnap = 0
  do i = 1, numFile
    inquire (file=files(i), exist=ifExist)
    if(.not. ifExist) call errormsg("traj file" // trim(files(i)) // "does not exist")
    call getfreeunit(iofile)
    open(file=files(i), unit=iofile, action="read")
    outer: do
      inner: do j = 1, frameStride * frameLine 
        read(iofile, *, iostat=ierr)
        if (ierr < 0) exit outer
      enddo inner
      nsnap = nsnap + 1
      numFrames(i) = numFrames(i) + 1
    enddo outer
    close(iofile)
  enddo
end subroutine

subroutine get_file_number(files, numFile)
  character(max_str_len), intent(in) :: files(:)
  integer, intent(out) :: numFile
  integer :: i

  i = 1
  do 
    if(files(i) == "" ) exit
    i = i + 1
  enddo
  numFile = i - 1
end subroutine

subroutine find_level_from_kelvinarray()
  integer :: i
  nlevel = 0
  do i = 1, maxnlevel
    if(kelvinarray(i) > 0) nlevel = nlevel + 1
  enddo
  allocate(levelindex(nsnap))
end subroutine

subroutine partition_process(procnum, nsnap, displs, counts, left, right)
  integer, intent(in) :: procnum, nsnap
  integer, intent(out) :: displs(:), counts(:)
  integer, intent(out) :: left, right
  integer :: i, length
  do i = 0, procnum-1
    left = i * nsnap / procnum + 1
    right = (i+1) * nsnap / procnum
    length = right - left + 1
    displs(i+1) = left - 1
    counts(i+1) = length
  enddo
  left = procid * nsnap / procnum + 1
  right = (procid+1) * nsnap / procnum
end subroutine

end module