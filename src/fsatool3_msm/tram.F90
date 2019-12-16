module tram
  use util
  use mod_global
  use cluster, only: maptocluster
  implicit none
  PRIVATE
   integer,parameter :: maxiter=5000
   real*8, allocatable :: nk(:, :), leveltcm(:, :, :), biask(:, :), unbiased_energy(:)
   integer, allocatable :: snapindex(:)
   logical, allocatable :: unbiased_keep(:), ifkeepArray(:)
   real*8, allocatable ::  pi(:), fk(:, :), logvk(:, :), logrk(:, :)
  PUBLIC :: tram_analysis, mod_tram
contains

! Multiensemble Markov models of molecular thermodynamics and kinetics.
! Hao Wua, Fabian Paula, Christoph Wehmeyer, and Frank Noe. PNAS, 2016, 25, 3221-3230
   subroutine mod_tram(inputfile, resultfile)
      use  mod_global, only: ncluster, nlevel, kelvinarray, trajindex, snappot, levelindex, lagstep, ifcutcluster 
      use cluster, only: maptocluster
      character(256) :: inputfile, resultfile, datafile
      integer :: iofile, ierr, i, index, index2, index3
      real*8 :: temppot
      real*8, allocatable :: tpm(: ,:)
      integer, allocatable :: reducedtcmindex(:)
      namelist /tram/ kelvinarray, nlevel, lagstep, datafile

      if(procid > 0) return
      resultdir = "."
      call getfreeunit(iofile)
      open(unit=iofile, file=trim(inputfile), action="read")
      read(iofile, nml=tram, iostat=ierr)
      if(ierr < 0) call errormsg("error in reading tram namelist")
      close(iofile)

      call  getfreeunit(iofile)
      open(iofile, file=trim(datafile), action="read")
      read(iofile, *) nsnap, ncluster
      allocate(maptocluster(nsnap), trajindex(nsnap), snappot(nsnap), levelindex(nsnap), ifcutcluster(nsnap))
      ifcutcluster=0
      do i = 1, nsnap
         read(iofile, *) index, index2, index3, temppot
         maptocluster(index) = index2
         trajindex(index) = index3
         snappot(index) = temppot
         levelindex(index) = ceiling(dble(index)*nlevel/nsnap)
      enddo
      close(iofile)
      call tram_analysis(tpm, reducedtcmindex)
      deallocate(logvk, fk, logrk, unbiased_energy, nk, leveltcm, biask, snapindex)
   end subroutine

   subroutine tram_analysis(tpm, reducedTcmIndex)
      use mpi
      real*8, allocatable :: tpm(:, :)
      integer, allocatable :: reducedTcmIndex(:)
      integer :: i, iofile

      if ( procid > 1 ) return
      call loginfo("Perfom Tram Analysis")
      write(*, "(a,i5,5x,a,100f10.3)")"nlevel ",nlevel,"kelvinarray ",kelvinarray(1:nlevel)

      call checkConnect(tpm, reducedTcmIndex)
      call markov_tram(leveltcm, nk, biask, nsnap, ncluster, nlevel)
      call tram_result(tpm, reducedTcmIndex)

      write(*,"(a)") "the pi and tpm has been write in the info/tram_tpm.txt file"
      write(*, "(a, I6, a)") "After tram analysis, the system has ", size(tpm, 1), " clusters"
      call getfreeunit(iofile)
      open(unit=iofile, file=trim(resultdir)//"/tram_tpm.txt", action="write")
      write(iofile, "(10E13.6)") pi
      do i = 1, size(tpm, 1)
         write(iofile, "(10E13.6)") tpm(i, :)
      enddo
      close(iofile)
      write(*, "(10E12.4)") pi
      call loginfo()
   end subroutine

   subroutine checkConnect(tpm, reducedTcmIndex)
      implicit none
      real*8, allocatable :: tpm(:, :)
      integer, allocatable :: reducedTcmIndex(:)
      integer :: ilevel, isnap, iclu, i, j, count
      real*8 :: betak
      real*8, DIMENSION(1:ncluster, 1:ncluster):: temparray2d
      logical, DIMENSION(1:ncluster, 1:nlevel) :: levelKeepSet ! .true. for keep, .false. for remove

      allocate(nk(ncluster, nlevel), leveltcm(ncluster, ncluster, nlevel), biask(nsnap, nlevel), &
               snapindex(nsnap), unbiased_keep(ncluster), ifKeepArray(ncluster))

      temparray2d = 0d0; nk=0d0; leveltcm = 0d0; biask=0d0

      !find connected sets and project connected sets
      do isnap = 1, nsnap
         nk(maptocluster(isnap), levelindex(isnap)) = nk(maptocluster(isnap), levelindex(isnap)) + 1.0d0
      enddo

      do ilevel=1,nlevel
         do isnap=1,nsnap
            betak = 1d0/(gasconst*kelvinarray(ilevel))
            biask(isnap,ilevel) = biask(isnap,ilevel)+snappot(isnap)*( betak - 1d0/(gasconst*kelvinarray(1)) )
         end do
         do isnap=1,nsnap-lagstep
            i=maptocluster(isnap); j=maptocluster(isnap+lagstep)
            if ( (ifcutcluster(i)==1) .or. (ifcutcluster(j)==1) ) cycle
            if ( trajindex(isnap) /= trajindex(isnap+lagstep) ) cycle
            if ( levelindex(isnap) == ilevel ) leveltcm(i,j,ilevel)=leveltcm(i,j,ilevel)+1.0d0
         end do
      end do

      call connectState(leveltcm)

      levelKeepSet = .false.
      do ilevel = 1, nlevel
         do iclu = 1, ncluster
            if((nk(iclu, ilevel) > 0d0) .and. ifKeepArray(iclu)) then
               levelKeepSet(iclu, ilevel) = .true.
            endif
         enddo
      enddo

      unbiased_keep(:) = levelKeepSet(:, 1)
      do ilevel = 1, nlevel
         do iclu = 1, ncluster
            if(levelKeepSet(iclu, ilevel) .eqv. .false.) then
               nk(iclu ,ilevel) = 0d0
               leveltcm(iclu, :, ilevel) = 0d0
               leveltcm(:, iclu, ilevel) = 0d0
           endif
         enddo
      enddo

      ! change the nsnap array, if the snap has been removed, then set the index to -1
      do i = 1, nsnap
         j = levelindex(i)
         if (levelKeepSet(maptocluster(i), j) .eqv. .true.) then
            snapindex(i) = maptocluster(i)
         else
            snapindex(i) = -1
         endif
      enddo

      ! find the  number of state in unbiased replica 
      count = 0;
      do i = 1, ncluster
         if(levelKeepSet(i, 1)) count = count+1;
      enddo
      allocate(tpm(count, count), pi(count), reducedTcmIndex(count))

      j = 1
      do i = 1, ncluster
         if(levelKeepSet(i, 1)) then
            reducedTcmIndex(j) = i
            j = j+ 1
         endif
      enddo
   end subroutine

   subroutine markov_tram(subTcm, subNk, subBiask, tpnsnap, tpncluster, tpnlevel)
      implicit none
      real*8, INTENT(IN) :: subTcm(:, : ,:), subNk(:, :), subBiask(:, :)
      integer, INTENT(IN) :: tpnsnap, tpncluster, tpnlevel
      integer :: iclu,jclu,ilevel,iter,isnap,i,j
      real*8 :: temparray(tpncluster), temparray2(tpnlevel), tpsum, tpsum2, tramerror, errormark, tempreal
      real*8,dimension(1:tpncluster,1:tpnlevel) :: newlogvk,newfk, stat_vectors, old_stat_vectors
      real*8 :: betak, delta, thermal_energies(tpnlevel), new_thermal_energies(tpnlevel)
      real*8 :: inf= huge(0d0), diff_energy

      allocate(logvk(tpncluster, tpnlevel), fk(tpncluster, tpnlevel), &
         logrk(tpncluster, tpnlevel), unbiased_energy(tpncluster))

      fk = 0d0; iter=0; logrk=0d0; tramerror=0d0; errormark=1d-15;thermal_energies=0d0; old_stat_vectors=0d0
      do ilevel = 1, tpnlevel
         do iclu = 1, tpncluster
            tpsum = 0d0
            do jclu= 1, tpncluster
               tpsum = tpsum + subTcm(iclu, jclu, ilevel) + subTcm(jclu, iclu, ilevel)
            enddo
            if (tpsum == 0d0) then
               logvk(iclu, ilevel) = -inf
            else
               logvk(iclu, ilevel) = log(0.5d0 * tpsum)
            endif
         enddo
      enddo

      do while((iter <= maxiter)) 
      if(mod(iter, 100) == 0 .and. iter /= 0) write(*,"(a,I5,2x,a, E12.4)") "tram has run", iter, "iterations, &
      the error is", diff_energy
         iter = iter + 1
         !update logvk. Eq(18)
         do ilevel=1,tpnlevel
            do iclu=1,tpncluster
               if(nk(iclu, ilevel) == 0) then
                  newlogvk(iclu, ilevel) = -INF
                  cycle
               endif
               i = 1
               do jclu=1,tpncluster
                  tpsum = leveltcm(iclu,jclu,ilevel) + leveltcm(jclu,iclu,ilevel)
                  if(jclu == iclu) then
                     if (leveltcm(iclu, jclu, ilevel) /=0 ) then
                        temparray(i) = log(leveltcm(iclu, jclu, ilevel))
                        i = i+1
                        cycle
                     endif
                  endif
                  if ( tpsum == 0d0 ) cycle
                  tpsum2 = logSumPair(logvk(jclu, ilevel)+fk(jclu,ilevel)-fk(iclu, ilevel)-logvk(iclu, ilevel), 0d0)
                  temparray(i) = log(tpsum) - tpsum2
                  i = i + 1
               end do
               newlogvk(iclu,ilevel) = logSumExp(temparray(1:i-1))
            end do
         end do

         !update R_i^k
         do ilevel = 1, tpnlevel
            do iclu = 1, tpncluster
               if(nk(iclu, ilevel) == 0) then
                  logrk(iclu, ilevel) = -inf
                  cycle
               endif
               i = 1
               tpsum = 0d0
               do jclu = 1, tpncluster
                  tpsum = leveltcm(jclu, iclu, ilevel) + tpsum
                  if(iclu == jclu) then
                     if (leveltcm(iclu, jclu, ilevel) /= 0d0) then
                        temparray(i) = log(leveltcm(iclu, jclu, ilevel))
                        temparray(i) = temparray(i) + fk(iclu, ilevel)
                        i = i+1
                        cycle
                     endif
                  endif
                  if((leveltcm(iclu,jclu,ilevel) + leveltcm(jclu, iclu, ilevel)) == 0d0) cycle
                  tpsum2 = logsumPair(newlogvk(jclu, ilevel) - fk(iclu, ilevel), newlogvk(iclu, ilevel) - fk(jclu, ilevel))
                  temparray(i) = log(leveltcm(iclu, jclu, ilevel) + leveltcm(jclu, iclu, ilevel))  &
                                       + newlogvk(jclu, ilevel) - tpsum2
                  i = i + 1
               enddo
               tempreal = nk(iclu ,ilevel) - tpsum
               if(tempreal > 0d0) then
                  tempreal = log(tempreal) + fk(iclu, ilevel)
               else
                  tempreal = -inf
               endif
               logrk(iclu, ilevel) = logsumpair(logsumexp(temparray(1:i-1)), tempreal)
            enddo
         enddo

         !calculation of  f_i^k, eq(19)
         newfk = inf
         do isnap = 1, tpnsnap
            if (snapindex(isnap) < 0) cycle
            i = 1
            j = snapindex(isnap)
            do ilevel = 1, tpnlevel
               if(logrk(j, ilevel) == -inf) cycle
               temparray2(i) = logrk(j, ilevel) - biask(isnap, ilevel)
               i = i+1
            enddo
            tpsum = logSumExp(temparray2(1:i-1))
            do ilevel=1, tpnlevel
               newfk(j, ilevel) = -logSumPair(-newfk(j, ilevel), -(tpsum+biask(isnap, ilevel)))
            enddo
         enddo


         tempreal = minval(newfk)
         newfk = newfk - tempreal

         do ilevel = 1,tpnlevel
            new_thermal_energies(ilevel)= -logsumexp(-newfk(:,ilevel))
         enddo

         tempreal = 0
         do ilevel = 1, tpnlevel
            stat_vectors(:, ilevel) = exp(thermal_energies(ilevel) - newfk(:, ilevel))
            tempreal = max(tempreal, maxval(abs(old_stat_vectors(:, ilevel) - stat_vectors(:, ilevel))))
         enddo
         diff_energy = max(tempreal, maxval(abs(thermal_energies-new_thermal_energies)))

         if (diff_energy <  errormark) exit
         fk = newfk
         logvk = newlogvk
         thermal_energies = new_thermal_energies
         old_stat_vectors = stat_vectors
      enddo

      do ilevel = 1,tpnlevel
         new_thermal_energies(ilevel)= -logsumexp(-newfk(:,ilevel))
      enddo

      ! get conf energies(temparray)
      unbiased_energy = inf
      do isnap = 1, tpnsnap
         if(snapindex(isnap) < 0) cycle
         i = 1; j=snapindex(isnap)
         do ilevel = 1,tpnlevel
            if(logrk(j, ilevel)==-inf) cycle
            temparray2(i) = logrk(j, ilevel) - biask(isnap, ilevel)
            i = i+1
         enddo
         tpsum = logSumExp(temparray2(1:i-1))
         unbiased_energy(j) = -logsumPair(-unbiased_energy(j), -tpsum)
      enddo

      ! normalize
      tempreal = -logsumexp(-unbiased_energy)
      fk = fk-tempreal
      new_thermal_energies = new_thermal_energies -tempreal
      unbiased_energy = unbiased_energy -tempreal
   end subroutine

   subroutine tram_result(tpm, reducedTcmIndex)
      implicit none
      real*8, intent(inout), allocatable :: tpm(:, :)
      integer, intent(inout), allocatable :: reducedTcmIndex(:)
      logical :: mask(ncluster, ncluster)
      integer :: ilevel, i, j, isnap, count, nums_scc, unbiased_count
      real*8 :: tpsum, tempdouble, temp_tpm(ncluster, ncluster), row_sum(ncluster)
      real*8 :: max_sum
      integer,ALLOCATABLE :: label(:)

      ! the p_ij, equ(16)
      temp_tpm = 0d0; row_sum=0d0
      do i=1,ncluster
         do j=1,ncluster
            tpsum = leveltcm(i, j, 1) + leveltcm(j, i , 1)
            if(tpsum == 0d0) cycle
            if (i == j) then
               temp_tpm(i, i) = 0.5 * tpsum * exp(-logvk(i, 1))
            else
               tempdouble = logSumPair(logvk(j, 1)- unbiased_energy(i), &
                                       logvk(i, 1) - unbiased_energy(j))
               temp_tpm(i, j) = tpsum * exp(-(unbiased_energy(j) + tempdouble))
            endif
            row_sum(i) = row_sum(i) + temp_tpm(i,j)
         end do
      end do

      max_sum = 0d0
      do i = 1, ncluster
         if(row_sum(i) > max_sum)  max_sum=row_sum(i)
      enddo
      do i = 1, ncluster
         do j = 1, ncluster
            if(i == j) then
               temp_tpm(i,j) = (temp_tpm(i,j) + max_sum - row_sum(i)) / max_sum
            else
               temp_tpm(i,j) = temp_tpm(i,j) / max_sum
            endif
         enddo
      enddo
     
      mask = .true.; unbiased_count = 0
      do i = 1, ncluster
         if(.not. unbiased_keep(i)) then
            mask(i, :) = .false.
            mask(:, i) = .false.
         else
            unbiased_count = unbiased_count + 1;
         endif
      enddo

      count = 0
      allocate(label(size(tpm, 1)))

      call find_undirected_scc(reshape(pack(temp_tpm, mask), [unbiased_count, unbiased_count]),&
       size(tpm , 1), label, nums_scc)

      do i = 1, size(tpm, 1)
         if(label(i) /= 1) then
            count = count + 1
            unbiased_keep(reducedTcmIndex(i)) = .false.
            mask(reducedTcmIndex(i), :) = .false.
            mask(:, reducedTcmIndex(i)) = .false.
         endif
      enddo

      if (count > 0) then
         unbiased_count = unbiased_count - count
         DEALLOCATE(tpm, reducedTcmIndex, pi)
         allocate(tpm(unbiased_count, unbiased_count), pi(unbiased_count), &
                  reducedTcmIndex(unbiased_count))
      endif

      tpm = reshape(pack(temp_tpm, mask), [unbiased_count, unbiased_count])
      pi = exp(pack(-unbiased_energy, unbiased_keep)) /  &
         sum(exp(pack(-unbiased_energy, unbiased_keep)))
      reducedTcmIndex = pack([(i, i=1, ncluster)], unbiased_keep)
      deallocate(label)
   end subroutine

   function logSumExp(array) result(value)
      implicit none
      real*8, intent(in) :: array(:)
      real*8 :: maxvalue, value

      maxvalue = maxval(array)
      value = maxvalue + log(sum(exp(array - maxvalue)))
   end function

   function logSumPair(a, b) result(value)
      implicit none
      real*8, intent(in) :: a,b
      real*8 :: value
      if( b > a) then
         value = b + log(1d0 + exp(a-b))
      else
         value = a+log(1d0 + exp(b-a))
      endif
   end function

   ! using the method of reversible pathways to check the states 

   subroutine connectState(tcms)
      implicit none
      real*8, intent(in) ::  tcms(:, :, :)
      real*8, allocatable :: tcm_temp(:, :)
      integer, allocatable :: label(:), filter(:), label_index(:)
      integer :: i,j,k,states, levels, nums_scc

      states = size(tcms, 1)
      levels = size(tcms, 3)
      allocate(label(states), tcm_temp(states,states))

      label_index = [(i,i=1,states)]
      ifKeepArray = .true.
      do i = 1, levels
         call find_scc(tcms(:,:,i), states, label, nums_scc)
         do j = 1, nums_scc
            filter = pack(label_index, label==j)
            do k = 1, size(filter, 1)-1
               tcm_temp(filter(k), filter(k+1)) = 1.0d0
            enddo
         enddo
      enddo
      call find_undirected_scc(tcm_temp(:, :), states, label, nums_scc)
      do i = 1, states
         if(label(i) /= 1) ifKeepArray(i) = .false.
      enddo
      DEALLOCATE(label, tcm_temp)
   end subroutine connectState
end module