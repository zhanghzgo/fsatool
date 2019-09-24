module mod_global
  use util
  implicit none
  integer, parameter :: max_str_len = 256
  integer,parameter :: maxnlevel=20
  real*8 :: kelvinarray(maxnlevel)
  integer :: numFileLevel(maxnlevel)
  integer :: ndim,nsnap,procid,procnum,ntraj,nss,ncv,ss2cv(20),trajtype,nlevel, ncluster, clustercycle
  integer :: cutnumber, lagstep, nstate
  integer,dimension(:),allocatable :: trajindex,snap2cluster,levelindex, ifcutcluster
  real*8,dimension(:,:),allocatable :: traj,cvs
  real*8,dimension(:),allocatable :: refcoor,snappot
  real*8 :: lagtime,kelvin
  real*8,parameter :: gasconst=1.9872065d-3
  logical :: ifreadcoor
  character(max_str_len),save :: resultdir, tpmmethod, clustermethod, lumpingmethod, inputfile
  integer, dimension(maxnlevel), target :: startstate, endstate
contains

subroutine mod_global_mpiinit()
  include 'mpif.h'
  character*(MPI_MAX_PROCESSOR_NAME) :: procName
  integer :: nameLen,ierr

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,procID,ierr)
  call mpi_comm_size(mpi_comm_world,procNum,ierr)
  call mpi_get_processor_name(procName,nameLen,ierr)
  !print "(i5,a4,i5,a4,a30)",procID," of ",procNum," at ",procName
  call mpi_barrier(mpi_comm_world,ierr)
end subroutine

subroutine mod_global_readtraj()
  include 'mpif.h'
  integer,parameter :: maxnfile=30
  integer :: narcfile,arcnums(maxnfile),readmode
  character*100 :: cvsfiles(maxnfile),coorfiles(maxnfile)
  integer :: iofile,iproc,ifile,index,iarc,isnap,ierr,iatom,ioarc,iocvs,i,j,k,narc,fragment,istart,istop
  integer :: status,nline,itraj,framestride,tpiarc, accumfile
  integer,dimension(0:procnum-1) :: tpstarts,tpstops,tpcountarray,tpposarray
  integer,dimension(:),allocatable :: tptrajindex
  real*8,dimension(:,:),allocatable :: tptraj,tpcvs
  real*8,dimension(:),allocatable :: tpsnappot
  real*8 :: centerpos(3),tpcoor(ndim),tpcvvec(ncv),tppot,tptime
  character*100 :: tpfilename,tpline
  logical :: ifexist
  namelist /trajs/ ncv, ss2cv, framestride, kelvinarray, trajtype,cvsfiles, coorfiles, numFileLevel

  cvsfiles=""; framestride = 1
  if ( procid == 0 ) then
    call loginfo("Read Collective Variables File")
    call getfreeunit(iofile); open(unit=iofile,file=inputfile,action="read")
    read(iofile,nml=trajs,iostat=ierr)
    if ( ierr < 0 ) call errormsg("error in reading trajs namelist")
    close(iofile)
    if ( trajtype == 1 ) then
      nline = int(ndim/10) + 1; if ( mod(ndim,10) == 0 ) nline = nline - 1
    else if ( trajtype == 2 ) then
      nline=1
    end if

    do i=1,ncv
        if ( ss2cv(i) == 0 ) exit
    end do
    nss = i-1

    if ( procid == 0 ) then
      print "(a,i4,a,i4,a,i4,a,f6.1)","ncv ",ncv,",  nss ",nss,",  framestride ",framestride
    end if

    do ifile = 1, maxnfile
      if ( cvsfiles(ifile) == "" ) exit
    end do
    narcfile=ifile-1; arcnums=0
    do ifile = 1, narcfile
      cvsfiles(ifile) = cvsfiles(ifile)
      coorfiles(ifile) = coorfiles(ifile)
      if ( trajtype == 1 ) then
          tpfilename = coorfiles(ifile)
      else if ( trajtype == 2 ) then
          tpfilename = cvsfiles(ifile)
      end if
      inquire(file=tpfilename,exist=ifexist)

      if ( ifexist .eqv. .false. ) call errormsg("traj file "//trim(tpfilename)//" does not exist!")
      call getfreeunit(ioarc); open(file=tpfilename,unit=ioarc,action="read")
      if ( trajtype == 1 ) read(ioarc,*)
      do
        do i=1,framestride*nline
          read(ioarc,"(a)",iostat=ierr) tpline
          if ( ierr < 0 ) exit
        end do
        if ( ierr < 0 ) exit
        arcnums(ifile) = arcnums(ifile) + 1
        ! if ( mod(arcnums(ifile),10000)==0 ) print *,"reading snapshot ",arcnums(ifile)
      end do
      close(ioarc)
    end do
    narc=sum(arcnums); nsnap = narc

    nlevel = 0
    do i = 1, maxnlevel
      if(kelvinarray(i) > 0.0d0) nlevel = nlevel + 1
    enddo

    allocate(levelindex(nsnap),trajindex(nsnap),cvs(ncv,nsnap),snappot(nsnap))
    levelindex=0; trajindex=0; cvs=0.0d0; snappot=0.0d0
    print "(a,i5,5x,a,i15)","narcfile ",narcfile,",  total structure no. ",narc

     if (nlevel/=1) then
       j = 1; i = 1; accumfile=0
       do ifile = 1, narcfile
         if ((ifile - accumfile) > numFileLevel(i)) then
           i = i + 1
           accumfile = accumfile + numFileLevel(i)
         endif
         levelindex(j:j + arcnums(ifile)-1) = i
         j = j + arcnums(ifile)
       enddo
     else
       levelindex= 1
     endif
  end if ! procid == 0

  call mpi_bcast(trajtype,1,mpi_integer,0,mpi_comm_world,ierr)
  call mpi_bcast(nline,1,mpi_integer,0,mpi_comm_world,ierr)
  call mpi_bcast(nsnap,1,mpi_integer,0,mpi_comm_world,ierr)
  call mpi_bcast(narcfile,1,mpi_integer,0,mpi_comm_world,ierr)
  call mpi_bcast(arcnums(1:narcfile),narcfile,mpi_integer,0,mpi_comm_world,ierr)
  call mpi_bcast(framestride,1,mpi_integer,0,mpi_comm_world,ierr)
  call mpi_bcast(ndim, 1, mpi_int, 0, mpi_comm_world, ierr)
  ! call mpi_bcast(ncv, 1, mpi_int, 0, mpi_comm_world, ierr)

  ! broadcast for msm cluster
  ! call mpi_bcast(ncluster, 1, mpi_int, 0, mpi_comm_world, ierr)
  ! call mpi_bcast(clustermethod, len(clustermethod), mpi_char, 0, mpi_comm_world, ierr)
  ! call mpi_bcast(clustercycle, 1, mpi_int, 0, mpi_comm_world, ierr)

  do i=1,narcfile
    call mpi_bcast(cvsfiles(i),100,mpi_character,0,mpi_comm_world,ierr)
    if ( trajtype == 1 ) call mpi_bcast(coorfiles(i),100,mpi_character,0,mpi_comm_world,ierr)
  end do
  call mpi_bcast(nss,1,mpi_integer,0,mpi_comm_world,ierr)
  call mpi_bcast(ss2cv(1:nss),nss,mpi_integer,0,mpi_comm_world,ierr)

  tpstarts=0
  tpstops=0
  fragment=narcfile/procNum
  if ( fragment < 1 ) fragment=1

  do i=0,procnum-1
      tpstarts(i)=i*fragment+1
      tpstops(i)=tpstarts(i)+fragment-1; 
      if ( tpstarts(i) > narcfile ) then
        tpstarts(i) = 0; tpstops(i) = 0
        exit
      else if ( tpstops(i) > narcfile) then
        tpstops(i) = narcfile
        exit
      else if ( i == (procnum-1) .and. (tpstops(i) < narcfile) ) then !this is a problm
        tpstops(i) = narcfile
      end if
  end do
  istart = tpstarts(procid)
  istop = tpstops(procid)

  if ( istart > 0 ) then
      i = istop - istart + 1
      print "(a,i4,a,i4,a,i4,a,i4)","proc ",procid,", arcfile ranges from ",istart," to ",istop,", total ",i
  end if
  call mpi_barrier(mpi_comm_world,ierr)

  if ( allocated(traj) ) deallocate(traj)
  allocate(traj(ndim,nsnap))
  traj=0.0d0
  iarc = 0
  call getfreeunit(ioarc)
  iocvs = 10*ioarc
  if ( istart > 0 ) then 
    i = sum(arcnums(istart:istop))
  else
    i = 1
  end if

  allocate(tptraj(ndim,i), tptrajindex(i),tpcvs(ncv,i),tpsnappot(i))
  tptraj=0.0d0; tptrajindex=0; tpcvs=0.0d0; tpsnappot=0.0d0
  do ifile=istart,istop
    if ( (ifile > narcfile) .or. (ifile==0) ) exit
    if ( trajtype == 1 ) then
      open(unit=ioarc,file=coorfiles(ifile),action="read")
      read(ioarc,*)
      open(unit=iocvs,file=cvsfiles(ifile),action="read")
    else if ( trajtype == 2 ) then
      open(unit=ioarc,file=cvsfiles(ifile),action="read")
    end if
    print "(a,i4,a,i4,a,a,a,i8)","proc ",procid," is reading file ",ifile," : ",trim(cvsfiles(ifile)),", ndata ",arcnums(ifile)
    tpiarc = 0
    do isnap=1,arcnums(ifile)
      iarc = iarc + 1
      tpiarc = tpiarc + 1
      do iatom=1,nline*(framestride-1)
        read(ioarc,*)
      end do
      if ( ierr < 0 ) exit
      if ( mod(iarc,10000) == 0 ) print *,"reading snapshot ",iarc
      if ( trajtype == 1 ) then
        centerpos(1:3)=0.0d0; tpcoor=0.0d0
        read(ioarc,"(10f8.3)") tpcoor
        do iatom=1,ndim/3
          centerpos(1:3) = centerpos(1:3) + tpcoor(iatom*3-2:iatom*3)
        end do
        centerpos(1:3) = centerpos(1:3)/dble(ndim/3)
        do iatom=1,ndim/3
          tpcoor(iatom*3-2:iatom*3) = tpcoor(iatom*3-2:iatom*3) - centerpos(1:3)
        end do
        ifexist = .false.
        do 
           read(iocvs, *) itraj, tpcvvec(1:ncv)
          !read(iocvs,"(i8,i3,f8.2,f11.3,10f8.3)") itraj,tpiproc,tptime,tppot,tpcvvec(1:ncv)
          !if ( abs(tptime-dble(tpiarc*framestride*coortimestep)/1000.0d0) < 1.0d-8 ) then
          !    ifexist = .true.; exit
          !end if
        end do
        if ( ifexist .eqv. .false. ) call errormsg("cvs data error!")
      else if ( trajtype == 2 ) then
        read(ioarc, *) itraj, tpcvvec(1:ncv)
        !read(ioarc,"(i8,i3,f8.2,f11.3,10f8.3)") itraj,tpiproc,tptime,tppot,tpcvvec(1:ncv)
        do i=1,ndim
          tpcoor(i) = tpcvvec(ss2cv(i))
        end do
      end if ! trajtype = 1
      tptraj(1:ndim,iarc) = tpcoor
      tptrajindex(iarc) = itraj
      tpcvs(1:ncv,iarc)=tpcvvec(1:ncv)
      !tpsnappot(iarc)=tppot
    end do
    close(ioarc); close(iocvs)
  end do


  if ( (istart>0) .and. (iarc /= sum(arcnums(istart:istop))) ) call errormsg("narc is not correct!")

  tpposarray = 0; tpcountarray = 0
  do iproc=0,procnum-1
    if ( tpstarts(iproc) == 0 ) exit
    tpcountarray(iproc) = ndim*sum(arcnums(tpstarts(iproc):tpstops(iproc)))
    if ( iproc > 0 ) tpposarray(iproc) = (sum(arcnums(1:tpstops(iproc-1))))*ndim
  end do
  i = 1
  if ( istart > 0 ) i = sum(arcnums(tpstarts(procid):tpstops(procid)))

  call mpi_gatherv(tptraj(1:ndim,1:i), tpcountarray(procid), mpi_double_precision, &
      traj, tpcountarray, tpposarray, mpi_double_precision,0,mpi_comm_world, ierr)
  deallocate(tptraj)

  tpposarray = 0; tpcountarray = 0
  do iproc=0,procnum-1
    if ( tpstarts(iproc) == 0 ) exit
    tpcountarray(iproc) = ncv*sum(arcnums(tpstarts(iproc):tpstops(iproc)))
    if ( iproc > 0 ) tpposarray(iproc) = (sum(arcnums(1:tpstops(iproc-1))))*ncv
  end do
  i = 1
  if ( istart > 0 ) i = sum(arcnums(tpstarts(procid):tpstops(procid)))
  call mpi_gatherv(tpcvs(1:ncv,1:i), tpcountarray(procid), mpi_double_precision, &
      cvs, tpcountarray, tpposarray, mpi_double_precision,0,mpi_comm_world, ierr)
  deallocate(tpcvs)


  tpposarray = 0; tpcountarray = 0
  do iproc=0,procnum-1
    if ( tpstarts(iproc) == 0 ) exit
    tpcountarray(iproc) = sum(arcnums(tpstarts(iproc):tpstops(iproc)))
    if ( iproc > 0 ) tpposarray(iproc) = (sum(arcnums(1:tpstops(iproc-1))))
  end do
  i = 1
  if ( istart > 0 ) i = sum(arcnums(tpstarts(procid):tpstops(procid)))
  call mpi_gatherv(tptrajindex(1:i), tpcountarray(procid), mpi_integer, &
      trajindex, tpcountarray, tpposarray, mpi_integer,0,mpi_comm_world, ierr)
  deallocate(tptrajindex)
  call mpi_gatherv(tpsnappot(1:i), tpcountarray(procid), mpi_double_precision, &
      snappot, tpcountarray, tpposarray, mpi_double_precision,0,mpi_comm_world, ierr)
  deallocate(tpsnappot)
  ! call mpi_bcast(traj, nsnap*ndim, mpi_double_precision, 0, mpi_comm_world, ierr)

  if ( procid == 0 ) then
    ntraj = maxval(trajindex)
    allocate(tptrajindex(ntraj)); tptrajindex=-1
    do i=1,nsnap
        j = trajindex(i)
        if ( tptrajindex(j) == -1 ) tptrajindex(j)=0
    end do
    i = ntraj + sum(tptrajindex)
    ntraj = i
    deallocate(tptrajindex)
    print "(a,i4,5x,a,i8,5x,a,i12)","nlevel ",nlevel,"ntraj ",ntraj,"nsnap ",nsnap
    call loginfo()
  end if
end subroutine
end module