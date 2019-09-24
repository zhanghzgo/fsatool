module tpt
  use util
  use mod_global, only: procid, nstate, ncluster, startstate, endstate, resultdir
  use fileio, only: stateindex, maptomacro => states_to_cluster
  use math, only: math_tpm_pi, math_inv
  use markov, only: tpm, pi=>limitpi
  implicit none
  integer :: num_path
  integer, parameter :: max_state_num = 10
  integer, allocatable :: Astate(:), Bstate(:)
  integer, allocatable :: strong_pathway(:)
  integer, pointer :: sendstate(:), recv_state(:)
  integer :: num_state_inA, num_state_inB
  real*8 :: totalflux,kab
  real*8,dimension(:,:),allocatable :: mfpt, netflux, flux, macro_net_flux
  real*8,dimension(:),allocatable :: forward_committor, backward_committor
  type :: pathwayinfo
     real*8 :: flux
     integer :: nodes
     integer, allocatable :: pathway(:)
     real*8 :: fluxratio
  end type pathwayinfo
  TYPE(pathwayinfo), allocatable, dimension(:) :: pathways
contains

  subroutine tpt_readfile(infile)
    character(256) :: infile, tpmfile, statefile
    integer :: i, j, k, prev_temp, iofile, ierr
    integer :: temp 
    namelist /tpt/ tpmfile, statefile,  startstate, endstate

    if (procid  > 0) return
    call getfreeunit(iofile)
    open(unit=iofile, file=trim(infile), action="read")
    read(iofile, nml=tpt,iostat=ierr)
    if (ierr<0) call errormsg("error in reading the tpt namelist")
    close(iofile)

    call getfreeunit(iofile)
    open(unit=iofile, file=trim(statefile), action="read")
    i = 0
    do
      read(iofile, *, iostat=ierr) nstate
      i = i + 1
      if (ierr < 0) exit
    enddo
    ncluster = i - 1
    close(iofile)

    allocate(tpm(ncluster, ncluster), pi(ncluster), &
    stateindex(nstate+1), maptomacro(ncluster))
    if(allocated(netflux)) deallocate(netflux)
    allocate(netflux(ncluster, ncluster))

    call getfreeunit(iofile)
    open(unit=iofile, file=trim(tpmfile), action="read")
    read(iofile, *) pi
    do i = 1, ncluster
      read(iofile, *) tpm(i, :)
    enddo
    close(iofile)

    call getfreeunit(iofile)
    open(unit=iofile, file=trim(statefile), action="read")
    prev_temp = 1; j = 1; stateindex(j) = 1
    do i = 1, ncluster
      read(iofile, *) temp, maptomacro(i)
      if (temp /= prev_temp) then 
        prev_temp = temp
        j = j + 1
        stateindex(j) = i
      endif
    enddo
    stateindex(j+1) = ncluster+1
    close(iofile)
  end subroutine

  subroutine mod_tpt(inputfile, outputfile)
    character(256) :: inputfile, outputfile
    if(procid > 0) return
    startstate = 0; endstate = 0
    call tpt_readfile(inputfile)
    call assign_states()
    allocate(flux(size(tpm, 1), size(tpm, 1)), macro_net_flux(nstate, nstate))
    resultdir = "./"
    call tpt_analyze(outputfile)
  end subroutine mod_tpt

  subroutine assign_states()
    integer :: i, j, k, iofile, count
    count = 0
    do i = 1,  max_state_num
      if (startstate(i) /= 0) then
        count = count + stateindex(startstate(i)+1) - stateindex(startstate(i))
      else
        exit
      endif
    enddo
    allocate(Astate(count))
    sendstate => startstate(1:i-1)
    num_state_inA = i - 1

    count = 0
    do i = 1, max_state_num
      if (endstate(i) /= 0) then
        count = count + stateindex(endstate(i)+1) - stateindex(endstate(i))
      else
        exit
      endif
    enddo
    allocate(Bstate(count))
    recv_state => endstate(1:i-1)
    num_state_inB = i - 1

    j = 1
    do k = 1, max_state_num
      if(startstate(k) /= 0) then
        do i = stateindex(startstate(k)), stateindex(startstate(k)+1) - 1
          Astate(j) = maptomacro(i)
          j = j + 1
        enddo
      endif
    enddo

    j = 1
    do k = 1, max_state_num
      if(endstate(k) /= 0) then
        do i = stateindex(endstate(k)), stateindex(endstate(k)+1)-1
          Bstate(j) = maptomacro(i)
          j = j + 1
        enddo
      endif
    enddo
  end subroutine

  subroutine tpt_analyze(outputfile)
    character*80, optional :: outputfile
    integer :: i, iofile
    character*80 :: resultfile = "macro_net_flux.txt"
    if (present(outputfile)) resultfile = outputfile
    call loginfo("TPT Analysis")
    call tpt_flux(tpm, size(tpm, 1), Astate, Bstate, flux, reversible=.true., pi=pi)
    call coarse_grain_flux(flux, macro_net_flux)
    call tpt_extract_pathways(macro_net_flux, num_path, sendstate, recv_state, totalflux)
    call printinfo(num_path)
    call getfreeunit(iofile)
    open(unit=iofile, file= trim(resultdir) //"/" // outputfile, action="write")
    do i = 1, size(macro_net_flux ,1)
      write(iofile, "(100f12.9)") macro_net_flux(i, :)
    enddo
    close(iofile)
    call loginfo()
  end subroutine

  subroutine tpt_msm_analysis()
    integer :: i
    stateindex(1) = 1
    do i = 2, nstate + 1
      stateindex(i) = stateindex(i) + 1
    enddo
    call assign_states()
    allocate(flux(size(tpm, 1), size(tpm, 1)), macro_net_flux(nstate, nstate))
    call tpt_analyze()
  end subroutine

  subroutine printinfo(num_path)
    integer :: num_path 
    integer :: i
    real*8 :: accumratio
    real*8,allocatable :: temp_array(:)
    integer, allocatable :: sortindex(:)

    write(*, "(a, I8)") "number of micro state: ", ncluster
    write(*, "(a, I5)")"number of macro state: ", nstate
    write(*, "(a, 10I4)") "int state: ", startstate(1:num_state_inA)
    write(*, "(a, 10I4)")"end state: ", endstate(1:num_state_inB)
    write(*, "(a, 50f12.8)")"forward committor: ",forward_committor(:)
    write(*, "(a, 500f12.8)")"backward committor: ",backward_committor(:)
    write(*, "(a, 50f10.6)")"pi: ", pi
    write(*, "(a, f12.8)")"totalflux: ", totalflux
    write(*, "(a, f12.8)")"kab: ", kab
    write(*, "(A, I3, A)")"There are total ", num_path, " pathways in the network"
    allocate(temp_array(num_path), sortindex(num_path))
    do i = 1, num_path
        temp_array(i) = pathways(i) % fluxratio
        sortindex(i) = i
    enddo
   ! numpath is not so big, so using simple sort to return sortindex
    accumratio = 0d0
    call get_feature_index(temp_array, num_path, sortindex)
    write(*, "(a, 5x, a, 2x, a, 2x, a)") "flux", "ratio", "accumratio", "pathway"
    do i = 1, num_path
      accumratio = accumratio + pathways(sortindex(i))%fluxratio
      write(*, "(E12.4, 2f8.4, 3x, 20I3)") pathways(sortindex(i))%flux, pathways(sortindex(i)) %fluxratio, &
      & accumratio, pathways(sortindex(i))%pathway
    enddo
  end subroutine

  subroutine tpt_flux(p, num_state, initstate, endstate, flux, reversible, pi)
    integer,intent(in) :: num_state
    real*8, intent(in):: p(num_state, num_state)
    real*8, intent(out) :: flux(num_state, num_state)
    logical, optional :: reversible
    real*8, dimension(:), optional :: pi
    integer, dimension(:) :: initstate, endstate
    real*8 :: lhs(num_state, num_state), rhs(num_state)
    real*8 :: subpi(num_state)
    integer :: i, j, ipiv(num_state), info

    lhs=p; rhs = 0.0d0
    ! calculate forward committor
    forall(i=1:num_state) lhs(i,i) = lhs(i,i) - 1.0d0
    do i = 1, size(initstate, 1)
      lhs(initstate(i), :)=0.0d0
      lhs(:, initstate(i)) = 0.0d0
      lhs(initstate(i), initstate(i)) = 1.0d0
    enddo
    do i =1, size(endstate, 1)
      lhs(endstate(i),:) = 0.0d0
      lhs(endstate(i), endstate(i)) = 1.0d0
      rhs(endstate(i)) = 1.0d0
    enddo

    call DGESV(num_state, 1, lhs, num_state, ipiv, rhs, num_state, info)
    forward_committor = rhs
    if (present(reversible) .and. reversible) then
      backward_committor = 1.0d0 - forward_committor
      subpi = pi
    else
      ! calculate backward committor
      call math_tpm_pi(p, num_state, subpi)
      rhs = 0.0d0;
      do i = 1, num_state
         do j = 1, num_state
            lhs(i,j) = subpi(j)* p(j,i) / subpi(i)
         enddo
      enddo
      forall(i=1:num_state) lhs(i,i) = lhs(i,i) - 1.0d0
      do i = 1, size(initstate, 1)
         lhs(initstate(i), :)=0.0d0
         lhs(initstate(i), initstate(i)) = 1.0d0
         rhs(initstate(i)) = 1.0d0
      enddo
      do i = 1, size(endstate, 1)
         lhs(endstate(i),:) = 0.0d0
         lhs(:, endstate(i)) = 0.0d0
         lhs(endstate(i), endstate(i)) = 1.0d0
      enddo
      call DGESV(num_state, 1, lhs, num_state, ipiv, rhs, num_state, info)
      backward_committor = rhs
    endif

    flux = 0.0d0
    forall(i=1:num_state, j=1:num_state, i/=j) flux(i,j)=subpi(i)*p(i,j)*backward_committor(i)*forward_committor(j)
  end subroutine 

  subroutine tpt_extract_pathways(netflux, numpath, initstate, endstate, totalflux)
    real*8 :: minflux, totalflux
    real*8,dimension(:, :) :: netflux
    integer :: numpath
    integer, dimension(:), pointer :: initstate, endstate
    real*8, allocatable :: netfluxcopy(:, :)
    integer :: i, nodes, index, j, k

   numpath = 0; allocate(pathways(100))
    do k = 1, size(initstate, 1)
      do j = 1, size(endstate, 1)
         netfluxcopy = netflux
         do
            call dijkstra(netfluxcopy, size(netfluxcopy, 1), initstate(k), endstate(j))
            nodes = size(strong_pathway, 1)
            if(nodes == 1) exit
            numpath = numpath + 1
            minflux = 1e10
            do i = 1, nodes-1
               if (minflux > netfluxcopy(strong_pathway(i), strong_pathway(i+1)))  then
                  minflux = netfluxcopy(strong_pathway(i), strong_pathway(i+1))
                  index = i
               endif
            enddo
            do i = 1, nodes-1
               netfluxcopy(strong_pathway(i), strong_pathway(i+1)) = &
                     netfluxcopy(strong_pathway(i), strong_pathway(i+1)) - &
                     minflux
            enddo
            netfluxcopy(strong_pathway(index), strong_pathway(index+1)) = 0.0d0
            pathways(numpath)%flux = minflux
            pathways(numpath)%nodes = nodes
            pathways(numpath)%pathway = strong_pathway
            pathways(numpath)%fluxratio =  minflux/ totalflux
         enddo
      enddo
   enddo
   if (numpath > 100) write(*,"(A)") "numpath is larger than allocated number of pathway, allocate pathway larger"
  end subroutine tpt_extract_pathways

  subroutine dijkstra(network, num_state, initstate, endstate)
    integer :: num_state, initstate, endstate
    real*8,intent(in) :: network(num_state, num_state)
    real*8 ::  dist(num_state), temp
    integer :: calculated_state(num_state), laststate(num_state)
    integer :: j, newstate, oldstate, templist(num_state)

    calculated_state = 0; dist = -huge(1.0); oldstate = endstate
    newstate = initstate; calculated_state(newstate) = 1
    do
       if(newstate == endstate) exit ! if find end, then exit
       if(newstate == oldstate) then ! if can't find other node next to intistate
          if (allocated(strong_pathway)) deallocate(strong_pathway)
          allocate(strong_pathway(1))
          return
       endif
       oldstate = newstate
       do  j = 1, num_state
          if(network(oldstate, j) > 0d0 .and. calculated_state(j) == 0) then
             if(network(oldstate, j) > dist(j))  then
                dist(j) = network(oldstate, j)
                laststate(j) = oldstate
             endif
          endif
       enddo
       temp = 0
       do j = 1, num_state
          if(calculated_state(j) == 0) then
             if(dist(j) > temp) then
                temp = dist(j)
                newstate = j
             endif
          endif
       enddo
       calculated_state(newstate) = 1
    enddo
    newstate = endstate; j = 1; templist(j) = endstate
    do
       if(newstate == initstate) exit
       newstate = laststate(newstate)
       j = j + 1
       templist(j) = newstate
    enddo
    if (allocated(strong_pathway)) deallocate(strong_pathway)
    allocate(strong_pathway(j)); strong_pathway = templist(j:1:-1)
  end subroutine dijkstra

  subroutine tpt_mfpt(tpm, num_state)
    integer :: num_state
    real*8 :: tpm(num_state, num_state)
    real*8 :: fundamental_matrix(num_state, num_state)
    real*8 :: temp_matrix(num_state, num_state)
    integer :: i, j

    allocate(mfpt(num_state, num_state))
    if (allocated(pi)) deallocate(pi)
    allocate(pi(num_state))
    call math_tpm_pi(tpm, num_state, pi)
    temp_matrix = - tpm
    forall(i=1:num_state)
       temp_matrix(i, i) = 1 + temp_matrix(i, i)
       temp_matrix(i, :) = temp_matrix(i, :) + pi(:)
    end forall
    call math_inv(temp_matrix, fundamental_matrix)
    mfpt = 0
    forall(i=1:num_state, j=1:num_state,i/=j)
       mfpt(i, j) = (fundamental_matrix(j,j) - fundamental_matrix(i, j)) / pi(j)
    end forall
    deallocate(pi)
  end subroutine tpt_mfpt

  subroutine coarse_grain_flux(flux, macro_net_flux)
    real*8, intent(in), dimension(:, :) :: flux
    real*8, intent(out), dimension(nstate, nstate) ::  macro_net_flux
    real*8, DIMENSION(nstate, nstate) :: macro_flux 
    real*8, dimension(nstate) :: pfold, macro_pi, qminus
    integer :: i, j, m, n
    real*8 :: flux_temp, tempsum, temp, tempsum2

    macro_flux = 0.d0
    do i = 1, nstate
      do j = 1, nstate
        flux_temp = 0d0
        if(i/=j) then
          do m = stateindex(i), stateindex(i+1) - 1
            do n = stateindex(j), stateindex(j+1) - 1
              flux_temp = flux_temp + flux(maptomacro(m), maptomacro(n))
            enddo
          enddo
          macro_flux(i, j) = flux_temp
        endif
      enddo
    enddo

    forall(i=1:nstate, j=1:nstate, i/=j) 
      macro_net_flux(i,j)=max(macro_flux(i,j)-macro_flux(j,i), 0.0d0)
    endforall

    totalflux = 0d0
    do i = 1, max_state_num
      if(endstate(i)/= 0) then
        totalflux = sum(macro_net_flux(:, endstate(i)))
      else
        exit
      endif
    enddo

    do i = 1, nstate
      temp = 0d0; tempsum=0d0; tempsum2 = 0d0
      do j = stateindex(i), stateindex(i+1) - 1
        temp = temp + pi(maptomacro(j))
      enddo
      macro_pi(i) = temp
      do j = stateindex(i), stateindex(i + 1) - 1
        tempsum = tempsum + forward_committor(maptomacro(j)) * pi(maptomacro(j))
        tempsum2 = tempsum2 + backward_committor(maptomacro(j)) * pi(maptomacro(j))
      enddo
      pfold(i) = tempsum/temp
      qminus(i) = tempsum2/temp
    enddo

    if (allocated(forward_committor)) then
      DEALLOCATE(forward_committor, pi, backward_committor)
      allocate(forward_committor(nstate), pi(nstate), backward_committor(nstate))
      forward_committor = pfold
      pi = macro_pi
      backward_committor = qminus
    endif
    kab = totalflux / sum(pi * backward_committor)
  end subroutine

  subroutine get_feature_index(tparray, features, sortindex)
    integer :: features
    real*8 :: temp, tparray(features)
    integer :: tempint, sortindex(features)
    integer :: i, j

    do i = 1, features
      sortindex(i) = i
    enddo
    do i = 1, features
      do j = i+1, features
          if (tparray(j) > tparray(i)) then
            temp = tparray(i)
            tparray(i) = tparray(j)
            tparray(j)  = temp

            tempint = sortindex(i)
            sortindex(i) = sortindex(j)
            sortindex(j) = tempint
          endif
      enddo
    enddo
  end subroutine
end module
