module coarse_grain
  use util
  use mod_global, only: nstate, lumpingmethod
  use math
  implicit none
  integer, allocatable :: maptomacro(:)
  real*8, allocatable :: eigenvector(:,:)
  real*8, allocatable :: membership(:, :)
  real*8, allocatable :: macrotpm(:, :)
contains

! Improved coarse-graining of Markov state models via explicit consideration of statistical uncertainty
! Gregory R. Bowman, J. Chem. Phys. 137, 134111 (2012);

  subroutine coarse_grain_analysis(tcmortpm, nstate)
    integer :: nstate, reduced_nstate
    integer :: i
    integer :: nmicro
    real*8, intent(in) :: tcmortpm(:, :)
    nmicro = size(tcmortpm, 1)

    if(lumpingmethod == "bace") then
       call loginfo("BACE Coarse Grain")
       call coarse_grain_bace(tcmortpm, nstate)
    elseif(lumpingmethod == "pcca") then
       call loginfo("PCCA Coarse Grain")
       call coarse_grain_pcca(tcmortpm, nstate)
    elseif (lumpingmethod == "pccaplus") then
       call loginfo("PCCAplus Coarse Grain")
       call coarse_grain_pccaplus(tcmortpm, nstate)
#ifdef DEBUG
       write(*, "(a)") "membership functions:"
       do i = 1, nmicro; print"(40f10.5)", membership(i, :);  enddo
#endif
       call coarse_grain_pccaplus_tpm(tcmortpm, membership, nstate)
      ! write(*, "(a)") "PCCA Coarse Grain TPM:"
      ! do i = 1, nstate; print"(40f10.5)", macrotpm(i, :);  enddo
    else
       call errormsg("lumpingmethod must be bace, pcca, or pccaplus")
    endif

    ! remap the number maptomacro to 1~nstate
    i = 1; reduced_nstate = nstate;
    do while(.true.)
       if(i > reduced_nstate) exit
       if (.not. any(maptomacro == i)) then
          where(maptomacro>i) maptomacro = maptomacro-1
          reduced_nstate = reduced_nstate - 1;
       endif
       i = i +1
    enddo

    write(*, "(A)") "cluster index belong to macro states index:"
    write(*, "(20I4)") maptomacro

    if (nstate > reduced_nstate)  then
       write(*, "(A, I4, 1x, A)") "After lumping only have", reduced_nstate, "macro states."
       write(*, "(A, I4, A, I3, A)")"Lumping algorithm can't coarse grain microstates to", &
       & nstate, " states, but has ", reduced_nstate, " states"
    endif
    nstate = reduced_nstate
    call loginfo()
  end subroutine coarse_grain_analysis

  subroutine mod_lumping(inputfile, resultfile)
   use mod_global, only: ncluster, nstate, lumpingmethod
   real*8, allocatable :: tpm(:, :), macrotpm(:, :), pi(:)
   character(256) :: inputfile, resultfile, datafile
   integer :: iofile, i, ierr
   namelist/lumping/ datafile, ncluster, nstate , lumpingmethod

   call getfreeunit(iofile)
   open(unit=iofile, file=trim(inputfile), action="read")
   read(iofile, nml=lumping, iostat=ierr)
   if(ierr<0) call errormsg("error in reading lumping namelist")
   close(iofile)

   allocate(tpm(ncluster, ncluster), macrotpm(nstate, nstate), pi(ncluster))

   call getfreeunit(iofile)
   open(unit=iofile, file=trim(datafile), action="read")
   ! read(iofile, *) pi
   do i = 1, ncluster
      read(iofile, *) tpm(i, :)
   enddo
   close(iofile)
   call coarse_grain_analysis(tpm, nstate)

   macrotpm = 0d0
   call getfreeunit(iofile)

   ! if(lumpingmethod == "bace") then
   !    do i=1,ncluster
   !       do j=1,ncluster
   !          macrotpm(maptomacro(i),maptomacro(j)) = macrotpm(maptomacro(i), maptomacro(j)) &
   !            + tpm(i,j)
   !       end do
   !    enddo
   !    forall(i=1:nstate)
   !       macrotpm(i, :) = macrotpm(i,:)/sum(macrotpm(i, :))
   !    end forall
   !    open(unit=iofile, file=trim(resultfile), action="write")
   !    do i = 1, nstate
   !       write(iofile, *) macrotpm(i, :)
   !    enddo
   ! endif
   open(unit=iofile, file=trim(resultfile), action="write")
   do i = 1, ncluster
      write(iofile, *) maptomacro(i)
   enddo
   close(iofile)
  
  end subroutine

  subroutine coarse_grain_pccaplus_tpm(tpm, membership, nstate)
    real*8, intent(in) :: tpm(:, :) , membership(:, :)
    integer, intent(in) :: nstate
    real*8 :: inverse_matrix(nstate, nstate)
    call math_inv(matmul(transpose(membership), membership), inverse_matrix)
    macrotpm=matmul(inverse_matrix, matmul(matmul(transpose(membership), tpm), membership))
  end subroutine coarse_grain_pccaplus_tpm

  subroutine coarse_grain_bace(tcm, nstate)
    real*8, intent(in) :: tcm(:, :)
    integer, intent(in) :: nstate
    integer :: i,j,istate,jstate,iter,imin,jmin,stateleft, tcmdim, minxminy(2), helparray(nstate)
    real*8,dimension(:), allocatable :: tpqvec, sumci
    real*8 :: spiq,spjq
    real*8,dimension(:,:),allocatable :: tptcm, tpm, bayesmat
    integer,allocatable :: allowedcalc(:,:), statesold(:)

    tcmdim = size(tcm,1)
    allocate(tpqvec(tcmdim), maptomacro(tcmdim), sumci(tcmdim), statesold(tcmdim))
    allocate(bayesmat(tcmdim, tcmdim), tpm(tcmdim, tcmdim), tptcm(tcmdim, tcmdim),allowedcalc(tcmdim, tcmdim))
    tptcm=tcm; bayesmat=0.0d0; stateleft=tcmdim; imin=0

    allowedcalc = 0
    where (tcm>1) allowedcalc = 1

    do i=1,tcmdim
       statesold(i) = i
       maptomacro(i)=i
       tptcm(i, 1:tcmdim) = tptcm(i, 1:tcmdim) + 1.0d0/dble(tcmdim)
       sumci(i) = sum(tptcm(i,1:tcmdim))
       tpm(i,1:tcmdim) = tptcm(i,1:tcmdim)/sumci(i)
    end do

    do iter=tcmdim,nstate+1,-1
       if (iter==tcmdim) then
          do istate=1,stateleft
             do jstate=istate+1, stateleft
                if (allowedcalc(istate, jstate) > 0) then
                   tpqvec(1:stateleft) = (tptcm(istate,1:stateleft)+tptcm(jstate,1:stateleft))/(sumci(istate)+sumci(jstate))
                   spiq = dot_product(tptcm(istate, 1:stateleft), log(tpm(istate,1:stateleft)/tpqvec(1:stateleft)))
                   spjq = dot_product(tptcm(jstate, 1:stateleft), log(tpm(jstate,1:stateleft)/tpqvec(1:stateleft)))
                   bayesmat(istate,jstate) = 1.0d0/(spiq + spjq)
                endif
             enddo
          enddo
       else
          do istate=1,stateleft
             if (allowedcalc(imin, istate)>0) then
                tpqvec(1:stateleft) = (tptcm(imin,1:stateleft)+tptcm(istate,1:stateleft))/(sumci(imin)+sumci(istate))
                spiq = dot_product(tptcm(imin, 1:stateleft), log(tpm(imin,1:stateleft)/tpqvec(1:stateleft)))
                spjq = dot_product(tptcm(istate, 1:stateleft), log(tpm(istate,1:stateleft)/tpqvec(1:stateleft)))
                bayesmat(imin,istate) = 1.0/(spiq + spjq)
             endif
          enddo
       endif

       minxminy = maxloc(bayesmat(1:stateleft, 1:stateleft));
       if (minxminy(1) < minxminy(2)) then
          imin = minxminy(1); jmin=minxminy(2)
       else
          imin = minxminy(2); jmin=minxminy(1)
       endif

       tptcm(imin, 1:stateleft) = tptcm(imin, 1:stateleft) + tptcm(jmin, 1:stateleft)
       tptcm(1:stateleft, imin) = tptcm(1:stateleft, imin) + tptcm(1:stateleft, jmin)
       bayesmat(imin, 1:stateleft) = 0
       bayesmat(1:stateleft, imin) = 0

       where (maptomacro == maptomacro(statesold(jmin)))
          maptomacro = maptomacro(statesold(imin))
       end where

       do jstate=jmin,stateleft-1
          tptcm(jstate,1:stateleft) = tptcm(jstate+1,1:stateleft)
          tptcm(1:stateleft,jstate) = tptcm(1:stateleft,jstate+1)
          bayesmat(jstate,1:stateleft) = bayesmat(jstate+1,1:stateleft)
          bayesmat(1:stateleft,jstate) = bayesmat(1:stateleft,jstate+1)
       end do

       do jstate = jmin+1, stateleft
          statesold(jstate-1) = statesold(jstate)
       enddo

       stateleft = stateleft - 1
       do jstate=1,stateleft
          sumci(jstate) = sum(tptcm(jstate,1:stateleft))
          tpm(jstate,1:stateleft) = tptcm(jstate,1:stateleft)/sumci(jstate)
       end do

       allowedcalc(imin, :) = 0
       do jstate=1, stateleft
          if(tptcm(imin, jstate)>1 .and. jstate /= imin) then
             allowedcalc(imin, jstate) = 1
          endif
       enddo
    enddo
    ! renumber the maptomacro to 1:n
    helparray=0; j=1
    do i =1, tcmdim
       if(any(helparray == maptomacro(i))) then
          maptomacro(i) = minloc(helparray, dim=1, mask=(helparray>=maptomacro(i)))
       else
          helparray(j) = maptomacro(i)
          maptomacro(i) = j
          j = j + 1
       endif
    enddo
    deallocate(tpqvec, sumci, statesold, bayesmat, tpm, tptcm, allowedcalc)
  end subroutine coarse_grain_bace

  subroutine coarse_grain_pcca(tpm, nstate)
    ! coarse_grain by pcca method
    real*8, dimension(:,:), intent(in) :: tpm
    integer, intent(in) :: nstate
    real*8, allocatable :: temparray(:), temparray2(:)
    integer :: nmicro,i,j,maxpos

    nmicro = size(tpm, 1)
    if (allocated(maptomacro)) deallocate(maptomacro)
    if (allocated(eigenvector)) deallocate(eigenvector)
    allocate(temparray(nmicro), eigenvector(nmicro, nmicro), temparray2(nmicro), maptomacro(nmicro))
    call math_solve_eigenproblem_tpm(tpm, temparray, eigenvector, temparray2)
    maptomacro=1
    do i = 1, nstate - 1
       temparray =  eigenvector(:, i+1)
       do j = 1, i
          temparray2(j) = maxval(pack(temparray, maptomacro==j))-minval(pack(temparray, maptomacro==j))
       enddo
       maxpos = maxloc(temparray2(1:j-1), dim=1)
       where((maptomacro == maxpos) .and. (temparray>0.0d0))
          maptomacro = i+1
       end where
    enddo
    deallocate(temparray, eigenvector, temparray2)
  end subroutine coarse_grain_pcca

  subroutine coarse_grain_pccaplus(tpm, nstate)
    real*8, intent(in) :: tpm(:, :)
    integer, intent(in) :: nstate
    integer :: nmicro, i, j, icount, index, nparameter, numres, ifault
    real*8 :: temp, dist2
    integer, allocatable :: stateindex(:)
    real*8, allocatable :: a(:, :), ainv(:, :)
    real*8, allocatable :: pi(:), eigenvalue(:), temparray(:)
    real*8, allocatable :: reducedeigvec(:,:), step(:), mincoor(:), xinit(:)

    nmicro = size(tpm, 1)
    allocate(pi(nmicro), eigenvalue(nmicro), eigenvector(nmicro, nmicro),&
         stateindex(nstate), temparray(nstate), a(nstate, nstate), ainv(nstate, nstate), &
         reducedeigvec(nmicro, nstate))
    call math_solve_eigenproblem_tpm(tpm, eigenvalue, eigenvector, pi)

    !normalize the eigenvalue by the weight pi
    do i = 1, nmicro
       temp = dot_product(eigenvector(:, i)**2, pi)
       eigenvector(:, i) = eigenvector(:, i) / sqrt(temp)
    enddo

    eigenvector(:, 1) = abs(eigenvector(:, 1))
    reducedeigvec = eigenvector(:, 1:nstate)

    ! Need to find the first index, select the most spread index
    ! https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
    stateindex = 0; dist2 = 0
    do i = 1, nmicro
       temp = dot_product(reducedeigvec(i,:), reducedeigvec(i, :))
       if (temp > dist2) then
          dist2 = temp
          stateindex(1) = i
       endif
    enddo
    index = stateindex(1)
    temparray = reducedeigvec(index, :)

    do i=1, nmicro
       reducedeigvec(i, :) = reducedeigvec(i, :) - temparray
    enddo

    do i = 1, nstate - 1
       dist2 = 0.0d0
       temparray = reducedeigvec(stateindex(i), :)
       do j = 1, nmicro
          reducedeigvec(j,:) = reducedeigvec(j, :) - dot_product(reducedeigvec(j, :), temparray)*temparray
          temp = dot_product(reducedeigvec(j, :), reducedeigvec(j, :))
          if ((temp > dist2) .and. (.not. any(stateindex == j))) then
             dist2  = temp
             stateindex(i+1) = j
          endif
       enddo
       temp = sqrt(dot_product(reducedeigvec(stateindex(i+1), :), reducedeigvec(stateindex(i+1), :)))
       reducedeigvec = reducedeigvec / temp
    enddo

    a =  eigenvector(stateindex, 1:nstate)
    call math_inv(a, ainv) ! get the inverse matrix of a
    nparameter = (nstate-1)*(nstate-1)
    allocate(mincoor(nparameter), step(nparameter), xinit(nparameter))
    xinit = pack(transpose(ainv(2:nstate, 2:nstate)), .true.)

    step = 0.05
    ! call Neled-Mead method to find mincoor

    membership = matmul(eigenvector(:, 1:nstate), ainv)
    call math_nelmin(object_func, nparameter, xinit, mincoor, temp, 0.00005d0, step, 10, 2000000, icount, numres, ifault)
    print*, "call Neled-Mead method to minimize the func"
    if (ifault == 2) then
       print*, "iteration terminated because kcount was exceed without convergence, maybe increase max steps"
    else
       print*, "icount, numres:", icount, numres
    endif

    ! rescale the mincoor to ainv
    ainv(2:nstate, 2:nstate) = reshape(mincoor, (/nstate-1, nstate-1/), order=[2,1])
    do i = 2, nstate
       ainv(i, 1) = -sum(ainv(i,2:nstate))
    enddo
    ainv(1, :) = -minval(matmul(eigenvector(:, 2:nstate), ainv(2:nstate, :)), dim=1)
    ainv = ainv / sum(ainv(1,:))

    ! At last, we can get membership function, normalize the membership function
    membership = matmul(eigenvector(:, 1:nstate), ainv)
    where (membership<0) membership = 0.0d0
    where (membership>1) membership = 1.0d0
    do i = 1, nmicro
       membership(i, :) = membership(i, :) / sum(membership(i, :))
    enddo

    ! get the maptomacro function
    maptomacro = maxloc(membership, dim=2)
    deallocate(pi, eigenvalue, eigenvector, stateindex, &
         temparray, a, ainv, reducedeigvec)
  end subroutine coarse_grain_pccaplus

  real*8 function object_func(ainv, nparameter)
    integer :: nparameter
    real*8,intent(in) :: ainv(nparameter)
    real*8,allocatable :: expandainv(:, :)
    integer :: n, i, j

    i = int(sqrt(dble(nparameter)))
    n = i + 1
    allocate(expandainv(n, n))
    expandainv(2:n, 2:n) = reshape(ainv, (/i, i/), order=[2,1])
    do i = 2, n
       expandainv(i, 1) = -sum(expandainv(i,2:n))
    enddo
    expandainv(1, :) = -minval(matmul(eigenvector(:, 2:n), expandainv(2:n, :)), dim=1)
    expandainv = expandainv / sum(expandainv(1,:))
    object_func = 0.0d0
    do i= 1, n
       do j=1, n
          object_func = object_func + (expandainv(j,i)**2) / expandainv(1,i)
       enddo
    enddo
    object_func = -object_func
    deallocate(expandainv)
    return
  end function object_func

end module coarse_grain
