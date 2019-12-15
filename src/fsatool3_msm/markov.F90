module markov
  use util
  use mod_global
  use cluster
  use coarse_grain, only: coarse_grain_analysis, maptomacro
  use math, only: math_solve_eigenproblem_tpm
  implicit none
  integer :: tcmdim
  integer,allocatable :: reducedtcmindex(:)
  real*8,dimension(:,:),allocatable :: tpm,tcm,eigenvec,macrotpm, impliedtime
  real*8,dimension(:),allocatable :: limitpi,eigenval,macropi
  real*8 :: lagfactor
contains


  subroutine markov_check(inputfile, resultfile)
    use cluster, only: maptocluster
    use mod_global, only:  nsnap, ncluster, trajindex, ifcutcluster
    character(256) :: inputfile, resultfile, datafile, checkmethod
    integer :: i, iofile, ierr, lagstart, lagend, nits
    namelist /check/ checkmethod, datafile, lagstart, lagend, nits, lagfactor

    nits = 5; lagfactor = 1.0d0
    call getfreeunit(iofile)
    open(unit=iofile, file=trim(inputfile), action="read")
    read(iofile, nml=check, iostat=ierr)
    if(ierr < 0) call errormsg("error in reading check namelist")
    close(iofile)

    allocate(maptocluster(nsnap), trajindex(nsnap), ifcutcluster(nsnap))
    ifcutcluster = 0

    call getfreeunit(iofile)
    open(unit=iofile, file=trim(datafile), action="read")
    read(iofile, *) nsnap, ncluster
    do i = 1, nsnap
      read(iofile, *) maptocluster(i), trajindex(i)
    enddo
    close(iofile)

    if (trim(checkmethod) == "timescales") then
      call markov_task1_impliedtimescale_vs_lagtime(resultfile, nits, lagstart, lagend)
    ! else if (trim(checkmethod) == "cktest")  then
    !   call makrov_task2_cktest()
    endif
  end subroutine

  subroutine markov_analysis()
    use tram, only:tram_analysis
    integer :: i, j, iofile
    integer,allocatable :: trimmed_cluster(:)
    allocate(ifcutcluster(ncluster));  ifcutcluster = 0

    if(procid > 0) return
    if (cutnumber /= 0) then
      do i = 1, ncluster
          if (numberincluster(i) < cutnumber) ifcutcluster(i) = 1
      enddo
    endif

    if(nlevel == 1) then
      call loginfo("Build MSM from clusters")
      call markov_buildtcmtpm_fromclusters()
      allocate(eigenval(tcmdim), eigenvec(tcmdim, tcmdim), limitpi(tcmdim))
      call math_solve_eigenproblem_tpm(tpm, eigenval, eigenvec, limitpi)

      ! save reduced tcm matrix
      call getfreeunit(iofile); open(unit=iofile, file=trim(resultdir)//"/tcm_reduce.txt", action="write")
      do i = 1, tcmdim
        write(iofile, "(10F6.0)") tcm(i, :)
      enddo
      close(iofile)

      ! save tpm matrix
      call getfreeunit(iofile); open(unit=iofile, file=trim(resultdir)//"/tpm.txt", action="write")
      write(iofile, "(10E12.5)") limitpi
      do i = 1, tcmdim
        write(iofile, "(10E12.5)") tpm(i, 1:tcmdim)
      enddo
      close(iofile)
    else
      call tram_analysis(tpm, reducedtcmindex)
      tcmdim = size(reducedtcmindex, 1)
    endif

    if((ncluster-tcmdim) > 0) then
      allocate(trimmed_cluster(ncluster-tcmdim))
      j = 1
      do i = 1, ncluster
          if (.not. any(i == reducedtcmindex(1:tcmdim))) then
            trimmed_cluster(j) = i
            j = j+1
          endif
      enddo
    endif
    write(*, "(a, I8, a, I5, a)") "System has ", nsnap, " microstates, corresponding  to", ncluster, " cluster number"
    if(cutnumber /= 0) then
      write(*,  "(a, I4)") "remove the cluster that has number of states which less than", cutnumber
    endif
    write(*, "(a)")"build tpm using the "//trim(tpmmethod)// "  method."
    write(*, "(a,i4)")"new dimension of TCM after reducing ",tcmdim
    write(*, "(a)") "cluster index which has been trimmed:"
    write(*, "(20I4)") trimmed_cluster
    call loginfo()

    if (lumpingmethod == "bace") then
       call coarse_grain_analysis(tcm, nstate)
    else
       call coarse_grain_analysis(tpm, nstate)
    end if

    ! call loginfo("Build MSM from Marco States")
    ! call markov_build_coarsetpm()
    ! write(*,"(a, I3, a, I5, a)") "The system has", nstate, " states, correspond to", tcmdim, " clusters."
    ! write(*, "(a)") "Marco tpm:"
    ! do i = 1,nstate; write(*, "(100f10.6)") macrotpm(i,:); enddo
    ! write(*, "(a)") "Pi of marco states:"
    ! write(*, "(10f8.3)")macropi

    ! call getfreeunit(iofile); open(unit=iofile, file=trim(resultdir)//"/macroinfo", action="write")
    !   write(iofile, "(100ES15.8)") macropi
    !   write(iofile, *)
    !   do i = 1,nstate; write(iofile, "(100ES15.8)") macrotpm(i,:); enddo
    ! close(iofile)
  end subroutine markov_analysis

  subroutine markov_buildtcmtpm_fromclusters()
    integer :: i,j,icluster,jcluster,isnap, nums_scc
    integer :: scc_label_array(ncluster)
    logical :: scc_filter(ncluster)
    real*8 :: temp,tpsum,tptcm(ncluster,ncluster)

    if ( allocated(tcm) ) deallocate(tcm);
    allocate(tcm(ncluster,ncluster)); tcm=0.0d0; tcmdim = 0

    do isnap=1,nsnap-lagstep
      i=maptocluster(isnap); j=maptocluster(isnap+lagstep)
      if ( (ifcutcluster(i)==0) .and. (ifcutcluster(j)==0) .and. (trajindex(isnap)==trajindex(isnap+lagstep))) then
          tcm(i,j)=tcm(i,j)+1.0d0
      end if
    end do

    tptcm = tcm
    call find_scc(tptcm, ncluster, scc_label_array, nums_scc)
    ! scc_filter is the logical array to keep the states at the largest scc
    do i = 1, ncluster
      if(scc_label_array(i) == 1) tcmdim = tcmdim + 1
    enddo

    if(allocated(tcm)) deallocate(tcm)
    if (allocated(reducedtcmindex)) deallocate(reducedtcmindex)
    if (allocated(tpm)) deallocate(tpm)
    allocate(reducedtcmindex(tcmdim), tcm(tcmdim, tcmdim), tpm(tcmdim, tcmdim))

    j = 1
    scc_filter = .true.
    do i = 1, ncluster
      if(scc_label_array(i) == 1) then
        reducedtcmindex(j) = i
        j = j+1
      else
        scc_filter(i) = .false.
      endif
    enddo

    j = 1;
    do i =1, ncluster
      if(scc_filter(i)) then
        tcm(j, :) = pack(tptcm(i, :), scc_filter)
        j = j+1
      endif
    enddo

    tpm=0d0
    select case (trim(tpmmethod))
    case("normal")
      tpsum = sum(tcm)
      do icluster=1,tcmdim
          temp = sum(tcm(icluster,:)); tpm(icluster,:)=tcm(icluster,:)/temp
      end do
    case("non-reversible")
      do icluster=1,tcmdim
          do jcluster=icluster+1,tcmdim
            temp = (tcm(icluster,jcluster)+tcm(jcluster,icluster))*0.5d0
            tcm(icluster,jcluster)=temp; tcm(jcluster,icluster)=temp
          end do
      end do
      tpsum = sum(tcm)
      do icluster=1,tcmdim
          temp = sum(tcm(icluster,:)); tpm(icluster,:)=tcm(icluster,:)/temp
          !limitpi(icluster) = temp/tpsum
      end do
    case("reversible")
      call markov_symmetrizetcm_and_normalizetpm()
    case default
      call markov_symmetrizetcm_and_normalizetpm()
   end select
  end subroutine markov_buildtcmtpm_fromclusters

  subroutine markov_build_coarsetpm()
    integer :: i,icluster,jcluster, index

    if ( allocated(macrotpm) ) deallocate(macrotpm)
    if ( allocated(macropi) ) deallocate(macropi)
    allocate(macrotpm(nstate,nstate),macropi(nstate)); macrotpm=0.0d0; macropi=0.0d0
    do icluster=1,tcmdim
      do jcluster=1,tcmdim
          macrotpm(maptomacro(icluster),maptomacro(jcluster)) = macrotpm(maptomacro(icluster), maptomacro(jcluster)) &
              + tcm(icluster,jcluster)
       end do
    enddo
    forall(i=1:nstate)
      macrotpm(i, :) = macrotpm(i,:)/sum(macrotpm(i, :))
    end forall
    do i = 1, tcmdim
      index = maptomacro(i)
      macropi(index) = numberincluster(reducedtcmindex(i)) + macropi(index)
    enddo
    macropi = macropi/sum(macropi)
  end subroutine markov_build_coarsetpm


  subroutine markov_timescale(tplagtime, tpimpliedtime, nlagtime)
    real*8,intent(in) :: tplagtime
    integer, intent(in) :: nlagtime
    real*8, intent(inout) :: tpimpliedtime(:)
    allocate(eigenval(tcmdim), eigenvec(tcmdim, tcmdim), limitpi(tcmdim))
    call math_solve_eigenproblem_tpm(tpm, eigenval, eigenvec, limitpi)
    tpimpliedtime(2:nlagtime) = (-tplagtime)/log(abs(eigenval(2:nlagtime)))
    deallocate(eigenval, eigenvec, limitpi)
  end subroutine markov_timescale

  !Book:An Introduction to Markov State Models and Their Application to Long Timescale Molecular Simulation
  ! Chapter 4.6.1, page 50

  subroutine markov_symmetrizetcm_and_normalizetpm()
    integer :: i,j,icycle
    real*8 :: temp,xivec(tcmdim),civec(tcmdim),qlog,lastqlog
    real*8,dimension(:,:),allocatable :: tpmatrix

    allocate(tpmatrix(tcmdim,tcmdim)); tpmatrix=0.0d0
    do i=1,tcmdim
      do j=i,tcmdim
          call random_number(temp)
          tpmatrix(i,j)=temp
          tpmatrix(j,i)=temp
      end do
    end do

    do i=1,tcmdim
      civec(i)=sum(tcm(i,:))
    end do

    do
      icycle = icycle + 1
      do i=1,tcmdim;  xivec(i)=sum(tpmatrix(i,:)); end do
      qlog=0.0d0
      do i=1,tcmdim
          do j=i,tcmdim
            temp = (civec(i)/xivec(i)) + (civec(j)/xivec(j))
            tpmatrix(i,j)=(tcm(i,j)+tcm(j,i))/temp; tpmatrix(j,i)=tpmatrix(i,j)
            if ( tpmatrix(i,j) > 0.0d0 ) qlog = qlog + tcm(i,j)*log(tpmatrix(i,j)/xivec(i))
            if ( (j > i) .and. (tpmatrix(j,i) > 0.0d0) ) qlog = qlog + tcm(j,i)*log(tpmatrix(j,i)/xivec(j))
          end do
      end do
      if ( abs(qlog-lastqlog) < 1d-10) exit
      lastqlog=qlog
    end do
    do i=1,tcmdim
      xivec(i)=sum(tpmatrix(i,:)); tpm(i,:) = tpmatrix(i,:)/xivec(i)
    end do
    deallocate(tpmatrix)
  end subroutine

  subroutine markov_task1_impliedtimescale_vs_lagtime(resultfile, nits, lag_start, lag_end)
    use mod_global, only: procid, procnum
    include "mpif.h"
    character(256) :: resultfile, ierr
    integer :: i, j, ioimplied, nits, lag_size, lag_start, lag_end, lag_temp
    integer, allocatable :: lag_array(:), sub_lagarray(:)
    integer,dimension(0:procnum) ::  job_index(procnum+1)
    integer :: left, right, length
    integer :: displs(procnum), counts(procnum)
    real*8, allocatable :: sub_impliedtime(:, :)

    if(procid == 0) then
      call loginfo("Task1: calculate the impliedtimescale on lagtime")
      i = 1
      lag_temp = lag_start

      if (lagfactor /= 1) then
        do while(lag_temp <= lag_end)
          lag_temp = ceiling(lag_temp * lagfactor)
          i = i + 1
        enddo

        lag_size = i-1
        allocate(lag_array(lag_size))
        lag_array(1) = lag_start
        do  i = 2, lag_size
          lag_array(i) = ceiling(lagfactor*lag_array(i-1))
        enddo
      else
        lag_size = lag_end - lag_start + 1
        allocate(lag_array(lag_size))
        do i = lag_start, lag_end
          lag_array(i-lag_start + 1) = i
        enddo
      endif

      allocate(impliedtime(nits, lag_size))
      call getfreeunit(ioimplied); open(unit=ioimplied,file=trim(resultfile),action="write")
      
      write(*, "(a, 10I5)") "lagstep:", lag_array
      write(ioimplied, *) "#lagstep #implied_timescale(nits-1)"
    endif

    call mpi_bcast(lag_size, 1, mpi_int, 0, mpi_comm_world, ierr)

    do i = 0, procnum - 1
      left = i * lag_size /procnum + 1
      right = (i+1) * lag_size / procnum
      length = right - left + 1
      if (right < left) length = 0
      displs(i+1) = left - 1
      counts(i+1) = length
    enddo

    left = procid * lag_size / procnum + 1
    right = (procid + 1 ) * lag_size / procnum
    length = right - left + 1

    if(length <= 0 ) length = 0
    allocate(sub_lagarray(length))
    allocate(sub_impliedtime(nits, length))
    
    call mpi_scatterv(lag_array, counts, displs, mpi_int, sub_lagarray, length, mpi_int, 0, mpi_comm_world, ierr)
    do i  = 1, size(sub_lagarray, 1)
      lagstep = dble(sub_lagarray(i))
      call markov_buildtcmtpm_fromclusters()
      call markov_timescale(dble(lagstep), sub_impliedtime(:, i), nits)
      write(*, "(a, I6, a)") "finished calculating", lagstep, " lagstep"
    enddo

    call mpi_gatherv(sub_impliedtime, length*nits, mpi_double,& 
      impliedtime, counts*nits, displs*nits, mpi_double, 0, mpi_comm_world, ierr)

    if(procid == 0) then
      do i = 1, lag_size
        lagstep = lag_array(i)
        write(ioimplied,"(i4, 1x, 20f12.6)") lagstep, log(impliedtime(2:nits, i))
      enddo
      close(ioimplied)
      call loginfo()
    endif
  endsubroutine markov_task1_impliedtimescale_vs_lagtime
end module markov
