! Chapter 2 Application of Markov State Models to Simulate Long Timescale Dynamics of Biological Macromolecules Lin-Tai Da*, Fu Kit Sheong*, Daniel-Adriano Silva*, and Xuhui Huang

module cluster
  use util
  use mod_global
  implicit none
  integer,dimension(:),allocatable :: maptocluster, centersnaps, numberincluster
contains

  subroutine cluster_analysis()
    integer :: stat, i, j, ierr
    real*8 :: total
    real*8 :: average(ndim), standard_deviation(ndim)
    real*8 :: normalize_traj(ndim, nsnap)

    !  first normalize the trajectory (ndim * nsnap)
    if (procid == 0) then
      average = sum(traj, dim=2)/nsnap
      do i = 1, ndim
        total = 0.0d0
        do j = 1, nsnap
          total = total + (traj(i,j) - average(i)) * (traj(i,j)-average(i))
        enddo
          standard_deviation(i) = sqrt(total/dble((nsnap-1)))
      enddo
      do i = 1, ndim
        do j = 1, nsnap
          normalize_traj(i, j) = (traj(i,j) - average(i))/(standard_deviation(i))
        enddo
      enddo
    endif

    if (trim(clustermethod) == "kmeans") then
       if(procid == 0) call loginfo("Running Kmeans++ Cluster Algorithm")
       call cluster_kmeans(normalize_traj, ndim, ncluster, nsnap, clustercycle, stat)
       if(procid == 0) then
          write(*, "(a, I4,2x,a)")"kmeans has ", ncluster, "clusters"
          if (stat >= clustercycle) then
              write(*, "(a,I4,2x,a)")"kmeans has run", stat, "iterations, maybe increase clustercycle larger"
          else
              write(*, "(a,I4,2x,a)")"kmeans has run", stat, "iterations and has converged"
          end if
          call loginfo()
       endif
    elseif (trim(clustermethod)=="kmedoids") then
       call loginfo("Running Kmedoids Cluster Algorithm")
       call cluster_kmedoids(normalize_traj, ndim, ncluster, nsnap, clustercycle, stat)
       if (stat == clustercycle) then
          write(*, "(a)") "Kmedoids has not convegerd, maybe increase clustercycle larger"
       else
          write(*, "(a, I4, a)") "Kmeans has run", stat, " times, and converged"
       endif
       call loginfo()
    else
       call errormsg("clustermethod must choose be kmeans or kmedoids")
    endif
  end subroutine cluster_analysis

  subroutine mod_cluster(inputfile, resultfile)
    use mod_global, only: ncluster, clustercycle, clustermethod, traj,  ndim, nsnap
    include "mpif.h"
    integer :: iofile, ierr, i
    character(256) :: inputfile, resultfile, datafile
    namelist /cluster/ ncluster, clustercycle, clustermethod, datafile, ndim, nsnap

    clustercycle = 400
    if (procid == 0) then
      call getfreeunit(iofile)
      open(unit = iofile, file=trim(inputfile), action="read")
      read(iofile, nml=cluster, iostat=ierr)
      if (ierr < 0 ) call errormsg("error in reading the cluster namelists")
      close(iofile)
      allocate(traj(ndim, nsnap))
      call getfreeunit(iofile)
      open(unit = iofile, file=trim(datafile), action="read")
      read(iofile, *) traj
      close(iofile)
    endif
    call mpi_bcast(ncluster, 1, mpi_int, 0, mpi_comm_world, ierr)
    call mpi_bcast(ndim, 1, mpi_int, 0, mpi_comm_world, ierr)
    call mpi_bcast(nsnap, 1, mpi_int, 0, mpi_comm_world, ierr)
    call mpi_bcast(clustermethod, len(clustermethod), mpi_char, 0, mpi_comm_world, ierr)
    call mpi_bcast(clustercycle, 1, mpi_int, 0, mpi_comm_world, ierr)

    call cluster_analysis()

    if (procid == 0) then
      call getfreeunit(iofile)
      open(unit = iofile, file=trim(resultfile), action="write")
      do i = 1, nsnap
          write(iofile, *) maptocluster(i), centersnaps(maptocluster(i))
      enddo
      close(iofile)
    endif
  end subroutine

  subroutine cluster_kmeans_init(coor_center, coor, ndim, ncluster, nsnap)
    integer :: i, j, ndim, ncluster, nsnap
    real*8 :: weights(nsnap), coor(ndim, nsnap), coor_center(ndim, ncluster)
    real*8 :: ran, tot, dist2, best
    weights = 1e20

    write(*, "(a)") "Using kmeans++ to select initial center "
    call random_seed()
    call random_number(ran)
    j = int(ran*nsnap) + 1
    coor_center(:, 1) = coor(:, j)
    do i = 2, ncluster
       tot = 0
       do j = 1, nsnap
          best = weights(j)
          dist2 = sum((coor(:, j) - coor_center(:, i-1)) ** 2)
          if (dist2 < best) then
             best = dist2
             weights(j) = best
          endif
          tot = tot + best
       enddo
       call random_number(ran)
       ran = ran * tot
       tot = 0
       do j =1, nsnap
          tot = tot + weights(j)
          if (tot > ran) exit
       enddo
       coor_center(:, i) = coor(:, j)
    enddo
  end subroutine

  subroutine cluster_kmeans(coor, ndim, ncluster, nsnap, niter, stat)
    include "mpif.h"
    integer :: ncluster, nsnap, niter, ndim
    integer :: i, j, k, last_center, now_center, stat
    real*8 :: dist2, best, dist, maxnum, minnum, tolerance, center_shift
    real*8 :: coor(ndim, nsnap), coor_center(ndim, ncluster), temp_dist(ncluster)
    real*8 :: old_coor_center(ndim, ncluster)
    real*8, allocatable :: sub_coor(:, :)
    integer, allocatable :: sub_maptocluster(:)
    real*8 :: subsums(ndim, ncluster), sums(ndim, ncluster)
    integer :: subnumbers(ncluster), numbers(ncluster)
    logical :: continue_run
    integer, DIMENSION(procnum) :: displs, counts
    integer :: left, right, length
    integer :: ierr

    tolerance = 1e-5
     if (procid == 0) then
      allocate(maptocluster(nsnap), centersnaps(ncluster), numberincluster(ncluster))
      call cluster_kmeans_init(coor_center, coor, ndim, ncluster, nsnap)
    endif

    ! find the start and end index of each processor
    do i = 0, procnum-1
      left = i * nsnap /procnum + 1
      right = (i+1) * nsnap/ procnum
      length = right - left + 1
      displs(i+1) = left - 1
      counts(i+1) = length
    enddo

    ! save the left and right index for the current processor
    left = procid * nsnap / procnum + 1
    right = (procid+1) * nsnap / procnum
    length = right - left + 1
    allocate(sub_coor(ndim, length), sub_maptocluster(length))


    ! broadcast the coordinate and scatter the data
    call mpi_bcast(coor_center, ncluster*ndim, mpi_double, 0, mpi_comm_world, ierr)
    call mpi_scatterv(coor, ndim*counts, ndim*displs, mpi_double, sub_coor, length*ndim, mpi_double, 0, mpi_comm_world, ierr)
    call mpi_barrier(mpi_comm_world, ierr)

    old_coor_center = coor_center
    center_shift = 0d0
    j = 0 

    do while (j < niter)
      j = j + 1
      call assignToKmeansCluster(sub_coor, coor_center, sub_maptocluster, subSums, subnumbers)
      call mpi_gatherv(sub_maptocluster, length, mpi_int, maptocluster, counts, displs, mpi_int, 0, mpi_comm_world, ierr)
      call mpi_reduce(subsums, sums, ncluster*ndim, mpi_double, mpi_sum, 0, mpi_comm_world, ierr)
      call mpi_reduce(subnumbers, numberincluster, ncluster, mpi_int, mpi_sum, 0, mpi_comm_world, ierr)
      if (procid == 0) then
        center_shift = 0d0
        do i = 1, ncluster
          coor_center(:, i) = sums(:, i)/numberincluster(i)
          center_shift = center_shift + sqrt(sum( (coor_center(:, i) - old_coor_center(:, i)) ** 2 ))
        enddo
        if(mod(j, 20) == 0) then
          write(*, "(a, I5, a, F8.5)") "kmeans has run", j ," iterations, the error is ", center_shift**2
        endif
      endif
      call mpi_bcast(coor_center, ndim*ncluster, mpi_double, 0, mpi_comm_world, ierr)
      call mpi_bcast(center_shift, 1, mpi_double, 0, mpi_comm_world, ierr)
      old_coor_center = coor_center
      if(center_shift**2 < tolerance) exit
    enddo

    ! finalize the kmeans cluster
    if(procid == 0) then
      temp_dist = 10000; stat = j
      do i = 1, nsnap
        dist = sum((coor(:, i) - coor_center(:, maptocluster(i)))**2)
        if (dist < temp_dist(maptocluster(i))) then
            centersnaps(maptocluster(i)) = i
            temp_dist(maptocluster(i)) = dist
        end if
      enddo
    endif
  end subroutine
  
  subroutine assignToKmeansCluster(data, center, assignment, sums, number)
    real*8 :: data(:, :), center(: ,:), sums(:, :)
    integer :: number(:)
    integer :: assignment(:), i, j
    real*8 :: min_dist
    real*8 :: dist

    sums = 0d0
    number = 0
    do i = 1, size(data, 2)
        min_dist = huge(1.0)
        do j = 1, size(center, 2)
            dist = sum((data(:, i) - center(:, j))**2)
            if (dist < min_dist) then
                min_dist = dist
                assignment(i) = j
            endif
        enddo
        sums(:, assignment(i)) = sums(:, assignment(i)) + data(:, i)
        number(assignment(i)) = number(assignment(i)) + 1
    enddo
  end subroutine

  subroutine cluster_kmedoids(coor, ndim, ncluster, nsnap, niter, stat)
    use util
    integer :: ncluster, nsnap, niter, ndim, stat
    real*8, allocatable :: distmat(:, :), percent
    integer, allocatable :: nowcentersnap(:), oneclusterindex(:)
    real*8 :: coor(ndim, nsnap), coor_center(ndim, ncluster)
    integer :: temparray(nsnap), temparray2(nsnap)
    integer :: i, j

    if(allocated(maptocluster)) deallocate(maptocluster, centersnaps, numberincluster)
    allocate(maptocluster(nsnap), centersnaps(ncluster), numberincluster(ncluster), &
        nowcentersnap(ncluster))

    if (nsnap<=10000) then ! calculate distance matrix
      allocate(distmat(nsnap, nsnap)); distmat = 0
      do i = 1, nsnap
          do j = i+1, nsnap
            distmat(i, j) = sum((coor(:, i) - coor(:, j)) ** 2)
            distmat(j, i) = distmat(i, j)
          enddo
      enddo

    ! find initialize cluster center from kmeans++
      call cluster_kmeans_init(coor_center, coor, ndim, ncluster, nsnap)

      temparray = (/(i, i=1,nsnap)/); stat=niter
      do i=1, niter
          maptocluster = minloc(distmat(centersnaps, :), dim=1)
          do j=1, ncluster
            allocate(oneclusterindex(count(maptocluster==j)))
            oneclusterindex = pack(temparray, maptocluster==j)
            nowcentersnap(j) = oneclusterindex(minloc(sum(distmat(oneclusterindex, oneclusterindex), dim=2), dim=1))
            deallocate(oneclusterindex)
          enddo
          call isort(nowcentersnap, temparray2, 1, nsnap)
          call isort(centersnaps, temparray2, 1, nsnap)
          if(all(nowcentersnap == centersnaps)) then
             stat = i
             exit
          else
             centersnaps = nowcentersnap
          endif
       enddo
       deallocate(nowcentersnap, distmat)
    else
       write(*, *)"nsnap is too large, unable to construct distance matrix, Using CLARANS method"
       percent = 0.001d0
       call  kmedoids_largedata(coor, ndim, ncluster, nsnap, niter, percent, stat)
    end if

    do i = 1, ncluster
      numberincluster(i) = count(maptocluster == i)
    enddo

  end subroutine cluster_kmedoids

  subroutine kmedoids_largedata(coor, ndim, ncluster, nsnap, niter, percent, stat)
    integer :: ncluster, nsnap, niter, ndim, stat
    real*8 :: percent
    real*8 :: coor(ndim, nsnap), coor_center(ndim, ncluster)

    integer :: i, k, nowstep, nTestPoints, swapIndex
    real*8 :: score, testScore
    real*8, allocatable :: testPoints(:, :)
    integer, allocatable :: testPointsIndex(:)
    logical :: isBetter

    nowstep = 1; nTestPoints = int(nsnap*percent); stat=niter; isBetter=.true.

    call cluster_kmeans_init(coor_center, coor, ndim, ncluster, nsnap)
    call assignToCluster(coor, coor_center, nsnap, ndim, ncluster)

    do while(nowstep < niter)
      if(isBetter .eqv. .false.) then
        stat = nowstep
        exit
      endif
      isBetter = .false.
      do k = 1, ncluster
        call generateNeighbors(ntestPoints, testPointsIndex, testPoints, coor, nsnap, ndim, k)
        call dist_score(coor, nsnap, ndim, k, coor_center(k, :), score)
        swapIndex = 0
        do i = 1, nTestPoints
          call dist_score(coor, nsnap, ndim, k, testPoints(i, :), testScore)
          if (testScore < score) then
            score = testScore
            swapIndex = i
          endif
        enddo
        if(swapIndex /= 0) then
          isBetter = .true.
          coor_center(k, :) = testPoints(swapIndex, :)
          centersnaps(k) = testPointsIndex(swapIndex)
        endif
        deallocate(testPoints, testPointsIndex)
      enddo
      call assignToCluster(coor, coor_center, nsnap, ndim, ncluster)
      nowstep = nowstep + 1
    enddo
    do i = 1, ncluster
      numberincluster(i) = count(maptocluster == i)
    enddo
  end subroutine

  subroutine generateNeighbors(nPoints, pointsIndex, points, coor, nsnap,ndim, index)
    integer ::  nPoints, nsnap, ndim, index
    real*8 ::  coor(ndim, nsnap)
    integer, allocatable :: subIndex(:), pointsIndex(:)
    real*8, allocatable :: randArray(:), points(:, :)

    integer :: i, j, count
    count = 0

    do i = 1, nsnap
        if(maptocluster(i) == index) then
            count= count+ 1
        endif
    enddo
    allocate(subIndex(count), randArray(count))

    j = 1
    do i =1, nsnap
        if(maptocluster(i) == index) then
            subIndex(j) = i
            j = j + 1
        endif
    enddo

    nPoints = min(count, nPoints)
    allocate(pointsIndex(nPoints), points(ndim, nPoints))

    call random_seed()
    call RANDOM_NUMBER(randArray)
    do i = 1, nPoints
      j = minloc(randArray, dim=1)
      randArray(j) = huge(1.0)
      pointsIndex(i) = subIndex(j)
      points(:, i) = coor(:, subIndex(j))
    enddo
  end subroutine

  subroutine dist_score(coor, nsnap, ndim, index, testPoint, score)
    integer :: ndim, nsnap, index
    real*8 :: coor(ndim, nsnap), testPoint(ndim), score
    integer :: i

    score = 0d0
    do i =1, nsnap
      if(maptocluster(i) == index) then
        score = score + dot_product(coor(:, i) - testPoint(:), coor(:,i) - testPoint(:))
      endif
    enddo
  end subroutine

  subroutine assignToCluster(coor, coor_center, nsnap, ndim, ncluster)
    integer :: nsnap, ndim, ncluster
    real*8 :: coor(ndim, nsnap), coor_center(ndim, ncluster)
    integer :: i, j
    real*8 :: minDist, temp

    do i = 1, nsnap
      minDist = 10e9
      do  j = 1, ncluster
        temp = dot_product(coor(:, i) - coor_center(:, j), coor(:, i) - coor_center(:, j))
        if (temp < minDist) then
          minDist = temp
          maptocluster(i) = j
        endif
      enddo
    enddo
  end subroutine

end module
