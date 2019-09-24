module fileio
  use util
  use mod_global
  use cluster, only: maptocluster, centersnaps, numberincluster
  use coarse_grain, only: maptomacro
  use markov, only: reducedtcmindex, tcmdim
  use math, only: math_lRMSD
  implicit none
  integer, allocatable :: clusterindex(:), clusters(:), stateindex(:)
  integer, allocatable :: macrocenter, states_to_cluster(:)
contains

  subroutine fileio_default_parameters()
    use mod_global, only: clustermethod, clustercycle, cutnumber, &
                          tpmmethod, lumpingmethod, &
                          nstate, ifreadcoor, trajtype
    clustermethod="kmeans"
    clustercycle=600
    cutnumber=0
    tpmmethod = "reversible"
    lumpingmethod = "pccaplus"
    nstate = 8
    ifreadcoor = .true.
    trajtype = 2
  end subroutine

  subroutine fileio_init_parameters()
    include 'mpif.h'
    logical :: file_exist
    integer :: ierr, system, iofile, i, j
    namelist /msm/ ndim, ncluster, clustermethod, clustercycle,cutnumber, &
         lagstep,ifreadcoor, lumpingmethod, tpmmethod, nstate, startstate, endstate

    call fileio_default_parameters()
    if(procid == 0) call loginfo("Initialize parameters")
    call getfreeunit(iofile); open(unit=iofile,file=trim(inputfile),action="read")
    read(iofile,nml=msm,iostat=ierr)
    if ( ierr < 0 ) call errormsg("error in reading the msm namelist!")
    close(iofile)
    if ( procid == 0 ) then
       if (ifreadcoor .eqv. .true.) then
          ierr = system("rm -rf "//trim(resultdir))
          ierr = system("mkdir -p "//trim(resultdir))
       endif
       write(*, "(a,a)") "result directory: ",trim(resultdir)
    end if
    if ( procid == 0 ) then 
      call loginfo()
    endif
  end subroutine fileio_init_parameters

  subroutine fileio_writeclusterinfo()
    integer :: iocluster,index,isnap,icluster,i,j,totalcount, ndim

    if(procid > 0) return
    allocate(clusterindex(ncluster+1),clusters(nsnap))
    clusterindex = 0
    totalcount = 0
    do i = 1, ncluster
       do j = 1, nsnap
          if (maptocluster(j) == i) then
             totalcount = totalcount + 1
             clusters(totalcount) = j ! clusters map cluster index to snap
          endif
       enddo
       clusterindex(i+1) = totalcount
    enddo

    call loginfo("Writing cluster information")
    write(*, "(a)") "Cluster information has been putted into:" // trim(resultdir)
    write(*, "(a)") "cluster info has two files clusterindex.txt and clusters.txt"

    call getfreeunit(iocluster); open(unit=iocluster,file=trim(resultdir)//"/clusterindex.txt",action="write")
    write(iocluster,"(10i10)") ncluster,nsnap
    do icluster=1,ncluster
       write(iocluster,"(i6,3i10,10f10.5)") icluster,clusterindex(icluster+1),clusterindex(icluster+1)-clusterindex(icluster), &
            centersnaps(icluster), traj(:, centersnaps(icluster))
    end do
    close(iocluster)
    
    call getfreeunit(iocluster); open(unit=iocluster,file=trim(resultdir)//"/cluster_forcheck.txt",action="write")
    do i = 1, nsnap
       write(iocluster,"(i6,3i10,10f10.5)") maptocluster(i), trajindex(i)
    end do
    close(iocluster)

    call getfreeunit(iocluster); open(unit=iocluster,file=trim(resultdir)//"/clusters.txt",action="write")
    do icluster=1,ncluster
       do i=clusterindex(icluster)+1,clusterindex(icluster+1)
          isnap = clusters(i)
          write(iocluster,"(i10,2i8,10f10.5)") isnap,icluster,trajindex(isnap), traj(:, isnap)
       end do
    end do
    close(iocluster)
    call loginfo()
  end subroutine fileio_writeclusterinfo

  subroutine fileio_readclusterinfo()
    integer :: iocluster,isnap,icluster,i,j,k,ierr

    if(procid > 0) return

    call loginfo("Reading cluster information")
    print "(a)", "Reading Cluster information from:" // trim(resultdir)
    print "(a)", "Reading two files clusterindex.txt and clusters.txt"

    call getfreeunit(iocluster); open(unit=iocluster,file=trim(resultdir)//"/clusterindex.txt",action="read")
    read(iocluster,"(10i10)") ncluster,nsnap
    print "(5(a,i7,2x))","read cluster file, cluster number: ",ncluster," snap number:",nsnap
    if ( allocated(clusters) ) deallocate(clusters)
    if ( allocated(clusterindex) ) deallocate(clusterindex)
    if ( allocated(centersnaps) ) deallocate(centersnaps)
    allocate(clusters(nsnap),clusterindex(ncluster+1),centersnaps(ncluster)); clusters=0; clusterindex=0; centersnaps=0
    if ( allocated(trajindex) ) deallocate(trajindex,cvs)
    if (allocated(numberincluster)) deallocate(numberincluster)
    if (allocated(maptocluster)) deallocate(maptocluster)
    allocate(trajindex(nsnap),cvs(ncv,nsnap),numberincluster(ncluster), maptocluster(nsnap))
    trajindex=0; cvs=0.0d0
    do icluster=1,ncluster
       read(iocluster,"(i6,3i10)") isnap,clusterindex(icluster+1),numberincluster(icluster),centersnaps(icluster)
    end do
    close(iocluster)

    call getfreeunit(iocluster); open(unit=iocluster,file=trim(resultdir)//"/clusters.txt",action="read")
    do icluster=1,ncluster
       do i=clusterindex(icluster)+1,clusterindex(icluster+1)
          read(iocluster,"(i10,2i8,10f10.6)") isnap,maptocluster(isnap),trajindex(isnap), cvs(1:ncv,isnap)
          clusters(i) = isnap
       end do
    end do
    close(iocluster)
    call loginfo()
  end subroutine fileio_readclusterinfo

  subroutine fileio_writestateinfo()

    integer :: i, j, k, numstate, index, istate, totalcount, isnap, iostate
    integer :: numcluster_state(nstate), temp(tcmdim), state_center_snap(nstate), helparray(tcmdim)
    integer ::  state_center_cluster(nstate), numsnap_state(nstate)
    real*8 :: total, bignum

    if(procid > 0) return
    allocate(states_to_cluster(tcmdim))

    totalcount = 0
    do i = 1, nstate
       do j = 1, tcmdim
          if(maptomacro(j) == i) then
             totalcount = totalcount + 1
             states_to_cluster(totalcount) = j !map same continuous state index to cluster index
          endif
       enddo
       numcluster_state(i) = count(maptomacro==i)
    enddo

    allocate(stateindex(nstate+1)); stateindex=0
    do i = 1, nstate
       stateindex(i+1) = stateindex(i) + numcluster_state(i)
    enddo

    numsnap_state = 0
    do i = 1, tcmdim
       numsnap_state(maptomacro(i)) = numsnap_state(maptomacro(i)) + numberincluster(reducedtcmindex(i))
    enddo

    helparray = (/(i, i=1, tcmdim)/)
    ! find state center correspond to cluster center
    do i = 1, nstate
       numstate = numcluster_state(i)
       temp(1:numstate) = centersnaps(reducedtcmindex(pack(helparray, maptomacro==i)))
       bignum = 1e10
       do j = 1, numstate
          total = 0.0d0
          do k = 1, numstate
             total = total + sum((traj(:, temp(j)) - traj(:, temp(k))) ** 2)
          enddo
          if (total < bignum) then
             bignum = total
             state_center_snap(i) = temp(j)
            !  state_center_cluster(i) = states(stateindex(i)+j)
             state_center_cluster(i) = maptocluster(temp(j))
          endif
       enddo
    enddo
    call getfreeunit(iostate); open(unit=iostate,file=trim(resultdir)//"/stateindex.txt",action="write")
    write(iostate,"(i6,2i10)") nstate, tcmdim, sum(numsnap_state)
    do i = 1, nstate
       write(iostate,"(5i10,i10,10f8.3)") i, stateindex(i+1),numcluster_state(i), numsnap_state(i), &
            state_center_cluster(i), state_center_snap(i), cvs(1:ndim, state_center_snap(i))
    enddo
    close(iostate)

    call getfreeunit(iostate); open(unit=iostate,file=trim(resultdir)//"/states.txt",action="write")
    index=0
    do istate=1,nstate
       do i=stateindex(istate)+1,stateindex(istate+1)
          index = index + 1
          isnap = centersnaps(reducedtcmindex(states_to_cluster(i)))
         !  write(iostate,"(3i8, 1x, i10,10f8.3)") istate,states(i),reducedtcmindex(states(i)), isnap, cvs(1:ndim, isnap)
          write(iostate,"(2i8, 1x, i10,10f8.3)") istate,states_to_cluster(i), isnap, cvs(1:ndim, isnap)
       end do
    end do
    close(iostate)
  end subroutine fileio_writestateinfo
end module fileio
