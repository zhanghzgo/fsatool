! Chapter 2 Application of Markov State Models to Simulate Long Timescale Dynamics of Biological Macromolecules Lin-Tai Da*, Fu Kit Sheong*, Daniel-Adriano Silva*, and Xuhui Huang

module cluster
  use math, only: math_lRMSD, math_euclidean_distance, math_normalization
  use mpi_shared_memory
  use util
  use mod_global
  use netcdf_func
  ! use coor_cluster
  implicit none
  integer,dimension(:),allocatable :: maptocluster, centersnaps, numberincluster
  ! real*8, pointer :: shared_array(:, :)
contains

  subroutine cluster_analysis()
    use mpi
    integer :: stat, i, j, ierr
    real*8, pointer :: normalized_cvs(:, :)
    integer :: win_normalized_shared

    call allocate_shared_memory_for_2dim_array([ndim, nsnap], normalized_cvs, win_normalized_shared)

    if (trajtype == 0) then
      if (shared_id == 0) then
        call math_normalization(cvs, normalized_cvs, ndim, nsnap)
      endif
    endif
    call mpi_win_fence(0, win_normalized_shared, ierr)

    if (trim(clustermethod) == "kmeans") then
       if(procid == 0) call loginfo("Running Kmeans++ Cluster Algorithm")
       call cluster_kmeans(normalized_cvs, ndim, ncluster, nsnap, clustercycle, stat, math_euclidean_distance)
    elseif (trim(clustermethod)=="kmedoids") then
      if(procid == 0) call loginfo("Running Kmedoids Cluster Algorithm")
      call cluster_kmedoids(normalized_cvs, ndim, ncluster, nsnap, clustercycle, stat, math_euclidean_distance)
    elseif (trim(clustermethod) == "coordinate") then
      if(procid == 0) call loginfo("Running Kmedoids Coordinate Cluster Algorithm")
      call cluster_kmedoids(traj, ndim, ncluster, nsnap, clustercycle, stat, math_lRMSD)
    else
       if(procid == 0) call errormsg("clustermethod must choose be kmeans, kmedoids or coordinate")
    endif

    if(procid == 0) then
      write(*, "(a, I4,2x,a)") trim(clustermethod) // " has ", ncluster, "clusters"
      if (stat >= clustercycle) then
          write(*, "(a,I4,2x,a)")"clustering has run", stat, "iterations, maybe increase clustercycle larger"
      else
          write(*, "(a,I4,2x,a)")"clustering  has run", stat, "iterations and converged"
      end if
      call loginfo()
    endif
    call mpi_win_fence(0, win_normalized_shared, ierr)
    call mpi_win_free(win_normalized_shared, ierr)
    ! call deallocate_shared_memory()
  end subroutine cluster_analysis

  subroutine mod_cluster(inputfile, resultfile)
    use mpi

    logical :: file_exist
    integer :: iofile, ierr, i, num_file, shared_win
    real*8 :: start, finish
    character(256) :: inputfile, resultfile
    character(256), dimension(100) :: datafile
    namelist /cluster/ ncluster, clustercycle, clustermethod, datafile, ndim, nsnap, trajtype

    clustercycle = 400
    trajtype = 0
    datafile = ""

    if (procid == 0) then
      call cpu_time(start)
      call getfreeunit(iofile)
      open(unit = iofile, file=trim(inputfile), action="read")
      read(iofile, nml=cluster, iostat=ierr)
      if (ierr < 0 ) call errormsg("error in reading the cluster namelists")
      close(iofile)
     
      i = 1
      do 
        if(trim(datafile(i)) /= "") then
           i = i + 1 
        else
            exit
        endif
      enddo
      num_file = i-1
    endif

    call mpi_bcast(ncluster, 1, mpi_int, 0, mpi_comm_world, ierr)
    call mpi_bcast(ndim, 1, mpi_int, 0, mpi_comm_world, ierr)
    call mpi_bcast(nsnap, 1, mpi_int, 0, mpi_comm_world, ierr)
    call mpi_bcast(clustermethod, len(clustermethod), mpi_char, 0, mpi_comm_world, ierr)
    call mpi_bcast(clustercycle, 1, mpi_int, 0, mpi_comm_world, ierr)
    call mpi_bcast(trajtype, 1, mpi_int, 0, mpi_comm_world, ierr)
    call mpi_bcast(num_file, 1, mpi_int, 0, mpi_comm_world, ierr)
    call mpi_bcast(datafile, 256*num_file, mpi_char, 0, mpi_comm_world, ierr)

    if(ncluster < procnum) STOP "procnum  must be less than the ncluster"
    if(nsnap < procnum) STOP "procnum  must be less than the nsnap"

    call init_shared_memory(mpi_comm_world)

    if (trajtype == 0) then
      call allocate_shared_memory_for_2dim_array((/ndim, nsnap/), cvs, shared_win)
    else
      call allocate_shared_memory_for_2dim_array((/ndim, nsnap/), traj, shared_win)
    endif

    ! print*, procid, ndim, nsnap, size(cvs, 1), size(cvs, 2)

    if (shared_id == 0) then
      if (trajtype == 0) then
          inquire(file=datafile(1), exist=file_exist)
          if (.not. file_exist) STOP "Data file not exist"
          call getfreeunit(iofile)
          open(unit = iofile, file=trim(datafile(1)), action="read")
          read(iofile, *) cvs
          close(iofile)
      else if(trajtype == 2) then
        call loginfo("Running Coordinate Cluster Algorithm")
        call read_NCcoorfile_into_coor(datafile, num_file, ndim, traj, 1)
      endif
    endif
    call mpi_win_fence(0, shared_win, ierr)
    call cluster_analysis()
    if (procid == 0) then
      call getfreeunit(iofile)
      open(unit = iofile, file=trim(resultfile), action="write")
      do i = 1, nsnap
          write(iofile, *) maptocluster(i), centersnaps(maptocluster(i))
      enddo
      close(iofile)
      resultDir="./"
      if(trajtype==2) call write_coor_cluster_result(traj, ndim, nsnap, [(i,i=1,ncluster)], resultDir)
      call cpu_time(finish)
      write(*, "(a, f8.2, a)") "The cluster program has run ", finish-start, " seconds"
    endif
  end subroutine

  subroutine cluster_kmeans_init(coor, nsnap, ndim, ncluster, dist_function, clusterCenterIndex)
    use mpi
    integer, intent(in) :: nsnap, ndim, ncluster
    real*8, intent(in) :: coor(ndim, nsnap)
    integer,intent(out)  :: clusterCenterIndex(ncluster)
    external :: dist_function

    integer :: i, j, left_snap, right_snap, snap_length, ierr
    real*8 :: weights(nsnap)
    real*8 :: ran, tot, dist2, sub_tot
    real*8 , allocatable :: sub_dist2s(:), sub_weights(:)
    integer, dimension(procnum) :: snap_displs, snap_counts

    if(procid == 0) then
      write(*, "(a)") "Using kmeans++ to select initial center "
    endif

    call partition_process(procnum, nsnap, snap_displs, snap_counts, left_snap, right_snap)
    snap_length = right_snap - left_snap + 1
    allocate(sub_dist2s(snap_length), sub_weights(snap_length))

    if(procid == 0) then
      call init_random_seed()
      call random_number(ran)
      j = int(ran*nsnap) + 1
    endif

    call mpi_bcast(j, 1, mpi_integer, 0, mpi_comm_world, ierr)
    clusterCenterIndex(1) = j
    sub_weights = huge(1.0)

    do i = 2, ncluster
      sub_tot = 0
      do j = left_snap, right_snap
        call dist_function(coor(:, j), coor(:, clusterCenterIndex(i-1)), dist2, ndim)
        if (dist2 < sub_weights(j-left_snap+1)) sub_weights(j-left_snap+1) = dist2
        sub_tot = sub_tot + sub_weights(j-left_snap+1)
      enddo
      call mpi_gatherv(sub_weights, snap_length, mpi_double, weights, &
                      snap_counts, snap_displs, mpi_double, 0, mpi_comm_world, ierr)
      call mpi_reduce(sub_tot, tot, 1, mpi_double, mpi_sum, 0, mpi_comm_world, ierr)
      
      if (procid == 0) then
        call random_number(ran)
        ran = ran * tot
        tot = 0
        do j = 1, nsnap
          tot = tot + weights(j)
          if (tot > ran) exit
        enddo
      endif
      call mpi_bcast(j, 1, mpi_integer, 0, mpi_comm_world, ierr)
      clusterCenterIndex(i) = j
      if(procid == 0) then
        if( mod (i, 20) == 0) write(*, *) "Kmeans++ init subroutine generates ", i, " centers"
      endif
    enddo
    deallocate(sub_weights, sub_dist2s)
  end subroutine

  subroutine cluster_kmeans(coor, ndim, ncluster, nsnap, niter, stat, dist_function)
    use mpi
    integer,intent(in) :: ncluster, nsnap, niter, ndim
    external :: dist_function
    integer :: i, j, k, last_center, now_center, stat
    real*8 :: dist2, maxnum, minnum, tolerance, center_shift
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

    allocate(maptocluster(nsnap), centersnaps(ncluster), numberincluster(ncluster))
    ! allocate(centersnaps(ncluster))
    call cluster_kmeans_init(coor, nsnap, ndim, ncluster, dist_function, centersnaps)

    if (procid == 0) then
      do i = 1, ncluster
        coor_center(:, i) = coor(:, centersnaps(i))
      enddo
    endif

    call partition_process(procnum, nsnap, displs, counts, left, right)

    length = right - left + 1
    allocate(sub_coor(ndim, length), sub_maptocluster(length))

    ! broadcast the coordinate and scatter the data
    call mpi_bcast(coor_center, ncluster*ndim, mpi_double, 0, mpi_comm_world, ierr)
    call mpi_scatterv(coor, ndim*counts, ndim*displs, mpi_double, sub_coor, length*ndim, mpi_double, 0, mpi_comm_world, ierr)
    ! call mpi_barrier(mpi_comm_world, ierr)

    old_coor_center = coor_center
    j = 0 

    do while (j < niter)
      j = j + 1
      call assignToKmeansCluster(sub_coor, coor_center, sub_maptocluster, subSums, subnumbers, dist_function)
      call mpi_gatherv(sub_maptocluster, length, mpi_int, maptocluster, counts, displs, mpi_int, 0, mpi_comm_world, ierr)
      call mpi_reduce(subsums, sums, ncluster*ndim, mpi_double, mpi_sum, 0, mpi_comm_world, ierr)
      call mpi_reduce(subnumbers, numberincluster, ncluster, mpi_int, mpi_sum, 0, mpi_comm_world, ierr)
      if (procid == 0) then
        center_shift = 0d0
        do i = 1, ncluster
          coor_center(:, i) = sums(:, i) / numberincluster(i)
          call dist_function(coor_center(:, i), old_coor_center(:, i), dist2, ndim)
          center_shift = center_shift + sqrt(dist2)
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
      temp_dist = 1e9; stat = j
      do i = 1, nsnap
        call dist_function(coor(:, i), coor_center(:, maptocluster(i)), dist2, ndim)
        if (dist2 < temp_dist(maptocluster(i))) then
            centersnaps(maptocluster(i)) = i
            temp_dist(maptocluster(i)) = dist2
        end if
      enddo
    endif
  end subroutine

  subroutine assignToKmeansCluster(data, center, assignment, sums, number, dist_function)
    real*8 :: data(:, :), center(: ,:), sums(:, :)
    external :: dist_function
    integer :: number(:)
    integer :: assignment(:), i, j
    real*8 :: min_dist
    real*8 :: dist2

    sums = 0d0
    number = 0
    do i = 1, size(data, 2)
        min_dist = 1e9
        do j = 1, size(center, 2)
            call dist_function(data(:, i), center(:, j), dist2, ndim)
            if (dist2 < min_dist) then
                min_dist = dist2
                assignment(i) = j
            endif
        enddo
        sums(:, assignment(i)) = sums(:, assignment(i)) + data(:, i)
        number(assignment(i)) = number(assignment(i)) + 1
    enddo
  end subroutine

  subroutine cluster_kmedoids(coor, ndim, ncluster, nsnap, niter, stat, dist_function)
    use mpi
    integer, intent(in) :: ncluster, nsnap, niter, ndim 
    real*8,  intent(in) :: coor(ndim, nsnap)
    integer, intent(out) :: stat
    external :: dist_function

    integer :: i
    real*8 ::  percent

    allocate(maptocluster(nsnap), centersnaps(ncluster), numberincluster(ncluster))

    if (nsnap<5000) then ! calculate distance matrix
      call kmedoids_smalldata(coor, ndim, ncluster, nsnap, niter, stat, dist_function)
    else
      if (procid == 0) then
        write(*, *)"nsnap is too large, unable to construct distance matrix, Using CLARANS method"
      endif
      percent = 100 / dble(nsnap)
      if(procid == 0) then 
        print*, "Maximum Generate number of test point is :" , int(percent * nsnap)
      endif
      call  kmedoids_largedata(coor, ndim, ncluster, nsnap, niter, percent, stat, dist_function)
    end if
    do i = 1, ncluster
      numberincluster(i) = count(maptocluster == i)
    enddo
  end subroutine cluster_kmedoids

  subroutine kmedoids_smalldata(coor, ndim, ncluster, nsnap, niter, stat, dist_function)
    use util
    use mpi
    integer, intent(in):: ndim, ncluster, nsnap, niter
    integer, intent(out) :: stat
    real*8, intent(in) :: coor(ndim, nsnap)
    external :: dist_function

    real*8, pointer :: distmat(:, :)
    real*8 :: min_value
    integer :: shared_distmat_win, left, right, ierr
    integer, dimension(procnum) :: counts, displs
    integer, allocatable :: oneclusterindex(:)
    integer :: i, j, k, temparray(nsnap), temparray2(nsnap), nowcentersnap(ncluster)

    call partition_process(procnum, nsnap, displs, counts, left, right)
    call allocate_shared_memory_for_2dim_array([nsnap, nsnap], distmat, shared_distmat_win)

    do i = left, right
      do j = i + 1, nsnap
        call dist_function(coor(:, i), coor(:, j), distmat(i, j), ndim)
        distmat(j, i) = distmat(i, j)
      enddo
    enddo

    call mpi_win_fence(0, shared_distmat_win, ierr)
    call cluster_kmeans_init(coor, nsnap, ndim, ncluster, dist_function, centersnaps)
    if(procid ==  0) then
      temparray = (/(i, i=1,nsnap)/)
      stat=niter
      do i=1, niter
        maptocluster = minloc(distmat(centersnaps, :), dim=1)
        do j=1, ncluster
          allocate(oneclusterindex(count(maptocluster==j)))
          oneclusterindex = pack(temparray, maptocluster==j)
          nowcentersnap(j) = oneclusterindex(minloc(sum(distmat(oneclusterindex, oneclusterindex), dim=2), dim=1))
          deallocate(oneclusterindex)
        enddo
        call isort(nowcentersnap, temparray2, 1, ncluster)
        call isort(centersnaps, temparray2, 1, ncluster)
        if(all(nowcentersnap == centersnaps)) then
            stat = i
            exit
        else
            centersnaps = nowcentersnap
        endif
      enddo
    endif
    ! deallocate(distmat)
    call mpi_win_fence(0, shared_distmat_win, ierr)
    call mpi_win_free(shared_distmat_win, ierr)
  end subroutine

  subroutine kmedoids_largedata(coor, ndim, ncluster, nsnap, niter, percent, stat, dist_function)
    use mpi
    integer, intent(in) :: ncluster, nsnap, ndim, niter
    real*8,  intent(in) :: coor(ndim, nsnap), percent
    integer, intent(out) :: stat
    external :: dist_function

    integer :: newCenterSnaps(ncluster)
    integer :: i, ierr
    real*8 :: time

    call cluster_kmeans_init(coor, nsnap, ndim, ncluster, dist_function, centersnaps)

    i = 0; stat=0
    do while(i < niter)
      if(procid == 0) write(*, *) "kmedoids has run ", i, "/", niter, " steps" 
      call assignToCluster(coor, ncluster, nsnap, ndim, dist_function)
      call calculateClusterCentroids(coor, ncluster, nsnap, ndim, newCenterSnaps, percent, dist_function)
      if (all(centersnaps == newCenterSnaps)) exit
      centersnaps = newCenterSnaps
      i = i + 1
      stat = i
      call mpi_barrier(mpi_comm_world, ierr)
    enddo
  end subroutine
  
  subroutine assignToCluster(coor, ncluster, nsnap, ndim, dist_function)
    use mpi
    integer, intent(in) :: ncluster, nsnap, ndim
    real*8, intent(in) :: coor(ndim, nsnap)
    external :: dist_function

    integer, allocatable :: sub_maptocluster(:)
    integer :: i, j, ierr
    integer :: left_snap, right_snap, snap_length
    integer, dimension(procnum) :: snap_displs, snap_counts
    real*8 :: min_value, dist2s


    ! partition snap
    call partition_process(procnum, nsnap, snap_displs, snap_counts, left_snap, right_snap)
    snap_length = right_snap - left_snap + 1

    allocate(sub_maptocluster(snap_length))

    do i = left_snap, right_snap
      min_value = huge(1.0)
      do j = 1, ncluster
        call dist_function(coor(:, i), coor(:, centersnaps(j)), dist2s, ndim)
        if(dist2s < min_value)  then
          min_value = dist2s
          sub_maptocluster(i-left_snap+1) = j
        endif
      enddo
    enddo
    call mpi_allgatherv(sub_maptocluster, snap_length, mpi_integer, maptocluster, snap_counts, &
                        snap_displs, mpi_integer, mpi_comm_world, ierr)

    deallocate(sub_maptocluster)
  end subroutine

  subroutine calculateClusterCentroids(coor, ncluster, nsnap, ndim, newCentersnaps, percent, dist_function)
    use mpi
    real*8, intent(in) :: coor(ndim, nsnap), percent
    integer, intent(in) :: ncluster, nsnap, ndim
    integer, intent(out) :: newCenterSnaps(ncluster)
    external :: dist_function
    integer, allocatable :: sub_centersnaps(:)
    integer :: cluster_displs(ncluster), cluster_counts(ncluster)
    integer :: i, j, ierr
    integer :: left_cluster, right_cluster, cluster_length

    call partition_process(procnum, ncluster, cluster_displs, cluster_counts, left_cluster, right_cluster)
    cluster_length = right_cluster - left_cluster + 1

    allocate(sub_centersnaps(cluster_length))

    do i = left_cluster, right_cluster
      call calculateSingleClusterCentroids(i, coor, nsnap, sub_centersnaps(i-left_cluster+1), ndim, percent, dist_function)
    enddo

    call mpi_allgatherv(sub_centersnaps, cluster_length, mpi_integer, newCenterSnaps, &
                        cluster_counts, cluster_displs, mpi_integer, mpi_comm_world, ierr)
    
    deallocate(sub_centersnaps)
  end subroutine

  subroutine calculateSingleClusterCentroids(icluster, coor, nsnap, icentersnap, ndim, percent, dist_function)
    integer, intent(in) :: icluster,  nsnap, ndim
    real*8, intent(in) :: percent
    real*8, intent(in) :: coor(ndim, nsnap)
    integer, intent(out) :: icentersnap
    external :: dist_function

    real*8 :: testscore, score
    integer, allocatable :: testPointsIndex(:)
    integer :: j, swapIndex
    integer :: nTestPoints

    nTestPoints = int(percent*nsnap)
    allocate(testPointsIndex(nTestPoints))
    ! do i = 1, ncluster
    call generateNeighbors(nTestPoints, icluster, testPointsIndex, coor, nsnap, ndim)
    call dist_score(coor, nsnap, ndim, icluster, centersnaps(icluster), score, dist_function)
    swapIndex = 0

    do j = 1, nTestPoints
      call dist_score(coor, nsnap, ndim, icluster, testPointsIndex(j), testscore, dist_function)
      if(testScore < score) then
        score = testScore
        swapIndex = j
      endif
    enddo
    if(swapIndex /= 0) then
      icentersnap = testPointsIndex(swapIndex)
    else
      icentersnap = centersnaps(icluster)
    endif
    deallocate(testpointsIndex)
  end subroutine

  subroutine generateNeighbors(nPoints, clusterIndex, pointsIndex, coor, nsnap, ndim)
    integer, intent(inout) :: nPoints
    integer, intent(in) ::  nsnap, ndim, clusterIndex
    real*8,  intent(in) ::  coor(ndim, nsnap)
    integer, intent(out) :: pointsIndex(nPoints)

    integer, allocatable :: subIndex(:)
    real*8 :: temp_real
    integer :: i, j, count, temp_int

    count = 0
    do i = 1, nsnap
      if (maptocluster(i) == clusterindex) then
        count= count+ 1
      endif
    enddo
    allocate(subIndex(count))
    if (count < nPoints) nPoints = count

    j = 1
    do i = 1, nsnap
        if (maptocluster(i) == clusterindex) then
            subIndex(j) = i
            j = j + 1
        endif
    enddo

    ! shuffle array
    call init_random_seed()
    do i = count, 1, -1
      call RANDOM_NUMBER(temp_real)
      j = int(ceiling(temp_real * i))
      temp_int = subIndex(i)
      subIndex(i) = subIndex(j)
      subIndex(j) = temp_int
    enddo
    pointsIndex(1:nPoints) = subIndex(1:nPoints)
   deallocate(subIndex)

  end subroutine

  subroutine dist_score(coor, nsnap, ndim, clusterindex, testPointIndex, score, dist_function)
    real*8, intent(in) :: coor(ndim, nsnap)
    integer, intent(in) :: nsnap, ndim
    real*8, intent(out) :: score
    external :: dist_function
    integer :: i, testPointIndex, clusterindex
    real*8 :: dist2

    score = 0d0
    do i = 1, nsnap
      if(maptocluster(i) == clusterindex) then
        call dist_function(coor(:, i), coor(:, testPointIndex), dist2, ndim)
        score = score + dist2
      endif
    enddo
  end subroutine

  subroutine write_coor_cluster_result(traj, ndim, nsnap, clusterIndex, resultDir)
    real*8, intent(in) :: traj(ndim, nsnap)
    integer, intent(in) :: ndim, nsnap
    integer, intent(in) :: clusterIndex(:)
    character(max_str_len) :: temp_char, resultDIr
    integer :: w_ncid, w_coorDVID, traj_index_counts
    integer :: i, j

    do i = 1, size(clusterIndex)
      write(temp_char, "(I6)") clusterIndex(i)
      temp_char = trim(resultDir)//"cluster_"//trim(adjustl(temp_char)) // ".nc"
      write(*, *) "writing cluster ", clusterIndex(i), "into " // trim(temp_char)
      call create_NCfile(temp_char, w_ncid, w_coorDVID)
      traj_index_counts = 1
      do j = 1, nsnap
         if(maptocluster(j) == clusterIndex(i)) then
            call write_NCcoor(w_ncid, w_coorDVID, reshape(traj(:, j), (/3, ndim/3/)), traj_index_counts)
            traj_index_counts = traj_index_counts + 1
         endif
      enddo
      call close_NCfile(w_ncid)
   enddo
  end subroutine

end module