! References
!
! Walking freely in the energy and temperature space by the modified replica exchange molecular dynamics method.
! Changjun Chen* and Yanzhao Huang
! J. Comput. Chem., 2016, 37: 1565-1575.
!
! Replica-exchange molecular dynamics method for protein folding.
! Yuji Sugita, Yuko Okamoto*
! Chem. Phys. Lett. 1999, 314: 141-151



module temd
  use moldata; use colvar
  implicit none

  integer,private :: montensample=0,nkelvin,oremdstep,moviecv
  real*8,private:: kelvinup,kelvindown
  integer,private :: exchangestep
  integer,private :: nkperproc,nexcperproc
  integer,dimension(:),allocatable,private :: allcoorindex
  real*8,dimension(:),allocatable,private :: kelvinarray,kelvinpot
  real*8,dimension(:,:,:),allocatable,private :: allcoor,allvel,allforce
contains

subroutine temd_init()
  integer :: ierr,iofile,iatom,ires,ikelvin,i,j,k,iproc,system
  real*8 :: kelvinbinratio=1.0d0
  logical :: ifexist
  real*8 :: temp
  namelist /temdss/nkelvin,kelvindown,kelvinup,kelvinbinratio,exchangestep,montensample,oremdstep,moviecv

  if ( procid == 0 ) write(iosiminfo,"(a)") "Perform REMD simulation"

  read(iomdin,nml=temdss,iostat=ierr)
  if ( (ierr < 0) .or. (nkelvin == 0) ) then
     call errormsg("temd namelist does not exist!")
  else if ( nkelvin == 1 ) then
     return
  end if
  if ( procnum > nkelvin ) call errormsg("procnum can not be larger than nkelvin")
  if ( montensample == 0 ) then
     nkperproc = nkelvin/procnum; nexcperproc = exchangestep/nkperproc
  end if  

  allocate(kelvinarray(0:nkelvin-1)); kelvinarray=0.0d0
  if ( kelvinbinratio == 1.0d0 ) then
     temp = log(kelvinup/kelvindown)
     do i=0,nkelvin-1
        kelvinarray(i)=kelvindown*exp(temp*dble(i)/dble(nkelvin-1))
     end do
  else
     temp = (kelvinup - kelvindown)*(1.0d0-kelvinbinratio)/(1.0d0-kelvinbinratio**(dble(nkelvin-1)))
     kelvinarray(0) = kelvindown
     do i=1,nkelvin-1
        kelvinarray(i) = kelvinarray(i-1) + temp*(kelvinbinratio**(dble(i-1)))
     end do
     kelvinarray(nkelvin-1) = kelvinup
  end if

  if ( procid == 0 ) then

     if ( montensample == 0 ) then
        write(iosiminfo,"(a)") "Standard REMD simulation"
        if ( procnum*nkperproc /= nkelvin ) call errormsg("nkelvin must be an integral multiple of procnum in standard REMD")
        if ( nkperproc*nexcperproc /= exchangestep ) then
          call errormsg("exchangestep must be an integral multiple of nkperproc ",int1=exchangestep,int2=nkperproc)
        end if
     else
        write(iosiminfo,"(a)") "Non-standard REMD simulation (TEMD)"
     end if
     write(iosiminfo,"(3(a,i6,5x))") "prcnum ",procnum,"nkelvin ",nkelvin
     write(iosiminfo, "(a,i5,5x,a,i10)") "samplestep ",samplestep,"excexcstep ",exchangestep*samplestep
     write(iosiminfo, "(a,f10.5,5x,a,i10)") "kelvin bin ratio ",kelvinbinratio
     write(iosiminfo, "(a,100f10.3)") "Temperature range: ",kelvinarray(0:nkelvin-1)
     write(iosiminfo, "(a,i5)") "save movie for cv ",moviecv
     write(iosiminfo, "(a,i10)") "OREMD step: ",oremdstep
     if ( montensample > 0 ) then
        write(iosiminfo, "(a,i10)") "monte carlo step number: ",montensample
     else
        write(iosiminfo,"(2(a,i6,5x))") "nkperproc ",nkperproc,"nexcperproc ",nexcperproc
     end if
     flush(iosiminfo)

     allocate(replicankelvin(0:procnum-1)); replicankelvin=0
     allocate(kelvinpot(0:nkelvin-1)); kelvinpot=0.0d0
  end if

  if ( montensample > 0 ) then
     nowikelvin = procid; kelvin0=kelvinarray(nowikelvin)
     if ( procid == 0 ) then
        do iproc=0,procnum-1
           replicankelvin(iproc) = iproc
        end do
     end if
  else
     nowikelvin = procid*nkperproc; kelvin0 = kelvinarray(nowikelvin)
     if ( procid == 0 ) then
        do iproc=0,procnum-1
           replicankelvin(iproc) = iproc*nkperproc
        end do
        if ( procnum < nkelvin ) then
           allocate(allcoor(3,numboxatom,0:nkelvin-1),allvel(3,numboxatom,0:nkelvin-1))
           allocate(allforce(3,numboxatom,0:nkelvin-1),allcoorindex(0:nkelvin-1))
           allcoor=0.0d0; allvel=0.0d0; allforce=0.0d0; allcoorindex=0
           do ikelvin=0,nkelvin-1
              allcoorindex(ikelvin)=ikelvin
           end do
        end if
     end if
  end if
end subroutine

subroutine temd_updatecoor()
  integer :: i,j,k,iss,ierr,irep,jrep,index,index2,ikelvin,jkelvin,status,iproc
  integer,dimension(0:nkelvin-1) :: kelvinreplica
  integer,dimension(0:procnum-1) :: tpindexarray,tpcountarray
  integer :: tpikelvin1,tpikelvin2,iadd,stkelvin
  integer,save :: itrial,ioexchangerate,nsample,iokelvinreplica,iokelvinhis,ioavgepot
  integer,allocatable,dimension(:),save :: jumparray,numepot,jumpsum,monteindex
  real*8 :: pab,xrand,temp,tpkelvin,temp2,cvvec(ncv)
  real*8,dimension(0:nkelvin-1) :: tpvec,exchangerate,avgepot
  real*8,dimension(:,:),allocatable,save :: prekelvinpot
  real*8,dimension(:),allocatable,save :: totepot
  real*8 :: coeffa,coeffb,coeffc,tpcoeffc,kslope,stpot
  character*100 :: tpfilename
  logical :: ifexchange,ifexceed,tplogical
  logical,save :: ifinitialdone

!  for OREMD
  real*8 :: tpkelvinprob(0:nkelvin-1),tpkelvinproberror
  real*8,dimension(:,:),allocatable,save :: repscaledfree,repkelvincounts,tprepkelvincounts

  if ( (montensample>0) .and. (procnum<nkelvin) .and. (ifinitialdone .eqv. .false.) ) then
     call temd_prepare(ifinitialdone)
     if ( (ifinitialdone .eqv. .true.) .and. (procid==0) ) then
        if ( procnum < nkelvin ) then
           if ( montensample > 0 ) then
              allocate(prekelvinpot(montensample,0:nkelvin-1),monteindex(0:nkelvin-1)); prekelvinpot=0.0d0; monteindex=0
              prekelvinpot(1,0:nkelvin-1)=kelvinpot(0:nkelvin-1)
           end if
        end if
     end if
     return
  end if

  if ( .not. allocated(totepot) ) then
     allocate(totepot(0:nkelvin-1),jumparray(0:nkelvin-1),numepot(0:nkelvin-1)); totepot=0.0d0; jumparray=0; numepot=0
     allocate(jumpsum(0:nkelvin-1)); jumpsum=0
     if ( procid == 0 ) then
        allocate(repscaledfree(0:nkelvin-1,0:procnum-1),tprepkelvincounts(0:nkelvin-1,0:procnum-1))
        allocate(repkelvincounts(0:nkelvin-1,0:procnum-1))
        repscaledfree=1.0d0; tprepkelvincounts=0.0d0; repkelvincounts=0.0d0
     end if
  end if

  nsample = nsample + 1

  if ( mod(nowstep,samplestep) == 0 ) then
     call moldata_download_crd(nowcoor)
     call colvar_calandrec(cvvec,tpmoviecv=moviecv)
  end if

  if ( nkelvin == 1 ) return

  if ( procnum < nkelvin ) then
     if ( (montensample==0) .and. (mod(nsample,nexcperproc) == 0) ) then
        call temd_collectpot()
        if ( nowikelvin < (procid*nkperproc+nkperproc-1) ) then
           nowikelvin = nowikelvin + 1; kelvin0=kelvinarray(nowikelvin)
           if ( procid == 0 ) then
              do iproc=0,procnum-1
                 replicankelvin(iproc) = replicankelvin(iproc) + 1
              end do
           end if
        end if
     else if ( montensample > 0 ) then
        call temd_collectpot()
        if ( procid == 0 ) then
           do iproc=0,procnum-1
              ikelvin=replicankelvin(iproc)
              monteindex(ikelvin) = monteindex(ikelvin) + 1
              if ( monteindex(ikelvin) > montensample ) monteindex(ikelvin)=1
              prekelvinpot(monteindex(ikelvin),ikelvin)=kelvinpot(ikelvin)
           end do
        end if
     end if
  else if ( mod(nsample,exchangestep) == 0 ) then
     call temd_collectpot()
  end if

  if ( mod(nsample,exchangestep) == 0 ) then

     if ( procid == 0 ) then
        kelvinreplica=-1
        do iproc=0,procnum-1
           ikelvin=replicankelvin(iproc); kelvinreplica(ikelvin)=iproc
        end do

        do iproc=0,procnum-1
           ikelvin=replicankelvin(iproc)
           repkelvincounts(ikelvin,iproc) = repkelvincounts(ikelvin,iproc) + 1.0d0
           tprepkelvincounts(ikelvin,iproc) = tprepkelvincounts(ikelvin,iproc) + 1.0d0
        end do

        if ( (oremdstep>0) .and. (mod(nsample,oremdstep) == 0) ) then
           do iproc=0,procnum-1
              temp = sum(tprepkelvincounts(:,iproc)); tpkelvinprob(:) = tprepkelvincounts(:,iproc)/temp
              do ikelvin=0,nkelvin-1
                 if ( tpkelvinprob(ikelvin) < 0.01d0 ) tpkelvinprob(ikelvin)=0.01d0
                 temp = (1.0d0/dble(nkelvin))/tpkelvinprob(ikelvin)
                 repscaledfree(ikelvin,iproc)=log(temp)+repscaledfree(ikelvin,iproc)
              end do
           end do
           tprepkelvincounts=0.0d0
        end if

        if ( (montensample > 0) .and. (procnum < nkelvin) ) then
           do ikelvin=0,nkelvin-1
              if ( kelvinreplica(ikelvin) > -1 ) cycle
              if ( prekelvinpot(montensample,ikelvin) == 0.0d0 ) then
                 i=monteindex(ikelvin)
              else
                 i=montensample
              end if
              call random_number(xrand); temp = 1.0d0/dble(i); j = int(xrand/temp)+1
              kelvinpot(ikelvin)=prekelvinpot(j,ikelvin)
           end do
        end if

        if ( (montensample == 0) .and. (procnum < nkelvin) ) then
           do ikelvin=0,nkelvin-1
              numepot(ikelvin) = numepot(ikelvin) + 1
              totepot(ikelvin) = totepot(ikelvin) + kelvinpot(ikelvin)
           end do
        else
           do iproc=0,procnum-1
              ikelvin=replicankelvin(iproc)
              numepot(ikelvin) = numepot(ikelvin) + 1
              totepot(ikelvin) = totepot(ikelvin) + kelvinpot(ikelvin)
           end do
        end if

        iadd = -1; tpikelvin1=nkelvin-1; tpikelvin2=1;
        do ikelvin=tpikelvin1,tpikelvin2,iadd
           jkelvin=ikelvin-1
           irep=kelvinreplica(ikelvin); jrep=kelvinreplica(jkelvin)
           if ( (montensample > 0) .and. (procnum < nkelvin) .and. (irep < 0) .and. (jrep < 0) ) cycle

           pab = - (kelvinpot(ikelvin) - kelvinpot(jkelvin))/(gasconst*kelvinarray(jkelvin))  &
                 - (kelvinpot(jkelvin) - kelvinpot(ikelvin))/(gasconst*kelvinarray(ikelvin))

           if ( oremdstep > 0 ) then
              if ( (irep>=0) .and. (jrep>=0) ) then
                 pab = pab + repscaledfree(jkelvin,irep) + repscaledfree(ikelvin,jrep)  &
                           - repscaledfree(jkelvin,jrep) - repscaledfree(ikelvin,irep)
              else if ( irep >= 0 ) then
                 pab = pab + repscaledfree(jkelvin,irep) - repscaledfree(ikelvin,irep)
              else if ( jrep >= 0 ) then
                 pab = pab + repscaledfree(ikelvin,jrep) - repscaledfree(jkelvin,jrep)
              end if
           end if

           pab = exp(pab)

           jumpsum(ikelvin)=jumpsum(ikelvin)+1
           ifexchange=.false.
           if ( pab > 1.0d0 ) then
              ifexchange=.true.
           else
              call random_number(xrand)
              if ( pab > xrand ) ifexchange=.true.
           end if

           if ( ifexchange .eqv. .true. ) then
              jumparray(ikelvin)=jumparray(ikelvin) + 1
              temp = kelvinpot(ikelvin); kelvinpot(ikelvin)=kelvinpot(jkelvin); kelvinpot(jkelvin)=temp
              j=kelvinreplica(ikelvin); kelvinreplica(ikelvin)=kelvinreplica(ikelvin-1); kelvinreplica(ikelvin-1)=j

              if ( (montensample == 0 ) .and. (procnum < nkelvin) ) then
                 j=allcoorindex(ikelvin); allcoorindex(ikelvin)=allcoorindex(ikelvin-1); allcoorindex(ikelvin-1)=j
              end if
           end if
        end do

        if ( (montensample == 0 ) .and. (procnum < nkelvin) ) then
           do iproc=0,procnum-1
              replicankelvin(iproc) = nkperproc*iproc
           end do
        else
           do ikelvin=0,nkelvin-1
              j = kelvinreplica(ikelvin)
              if ( j > -1 ) replicankelvin(j)=ikelvin
           end do
        end if
     end if

     call temd_distributekelvin()
     kelvin0=kelvinarray(nowikelvin)

  end if

  if ( (procid == 0) .and. (mod(nsample,exchangestep*100)==0) ) then
     if ( ioexchangerate == 0 ) then
        call getfreeunit(iokelvinhis); open(file="procinfo/kelvin_histogram.txt", unit=iokelvinhis, action="write")
        call getfreeunit(ioavgepot); open(file="procinfo/kelvin_avgepot.txt", unit=ioavgepot, action="write")
        call getfreeunit(ioexchangerate); open(file="procinfo/kelvin_exchangerate.txt", unit=ioexchangerate, action="write")
     end if

     tpkelvinproberror=0.0d0
     do iproc=0,procnum-1
        temp = sum(repkelvincounts(:,iproc))
        tpkelvinprob(:) = repkelvincounts(:,iproc)/temp - (1.0d0/dble(nkelvin))
        tpkelvinproberror = tpkelvinproberror + dot_product(tpkelvinprob,tpkelvinprob)
     end do
     tpkelvinproberror = sqrt(tpkelvinproberror)
     temp = sum(numepot(:))
     write(iokelvinhis,"(f10.1,50f8.4)") deltatime*dble(nowstep),tpkelvinproberror,dble(numepot(:))/temp
     flush(iokelvinhis)
     do ikelvin=0,nkelvin-1
        if ( (ikelvin>0) .and. (jumpsum(ikelvin) > 0) ) exchangerate(ikelvin)=dble(jumparray(ikelvin))/dble(jumpsum(ikelvin))
        if ( numepot(ikelvin) > 0 ) avgepot(ikelvin) = totepot(ikelvin)/dble(numepot(ikelvin)) 
     end do
     write(ioavgepot,"(f10.1,50f8.3)") deltatime*dble(nowstep),avgepot(:)
     flush(ioavgepot)
     temp = sum(exchangerate)/dble(nkelvin-1)
     write(ioexchangerate,"(f10.1,50f8.4)") deltatime*dble(nowstep),temp,exchangerate(1:nkelvin-1)
     flush(ioexchangerate)
  end if

end subroutine

subroutine temd_prepare(ifinitialdone)
  logical,intent(out) :: ifinitialdone
  integer,dimension(0:procnum-1) :: tpcountarray,tpposarray
  integer :: i,ierr,ikelvin,index
  real*8 :: tpkelvin
  integer,save :: nsample
  logical,save :: ifnotfirst

  nsample = nsample + 1
  if ( mod(nsample,exchangestep) /= 0 ) return

  if ( ifnotfirst .eqv. .false.) then
     ifnotfirst = .true.
  else
     if ( procid == 0 ) then
        tpcountarray(0:procnum-1)=1
        do i=0,procnum-1; tpposarray(i)=i; end do
     end if
     call mpi_gatherv(nowikelvin, 1, mpi_integer, &
             replicankelvin, tpcountarray, tpposarray, mpi_integer, 0, mpi_comm_world, ierr)

! epot: instant energy
     if ( procid == 0 ) then
        do i=0,procnum-1
           ikelvin=replicankelvin(i)
           tpposarray(i)=ikelvin
        end do
     end if
     call mpi_gatherv(epot, 1, mpi_double_precision, &
          kelvinpot, tpcountarray, tpposarray, mpi_double_precision,0, mpi_comm_world, ierr)

     nowikelvin = nowikelvin + procnum
     if ( nowikelvin > nkelvin-1 ) then
        nowikelvin = nowikelvin - procnum
        if ( procid == 0 ) ifinitialdone=.true.
     end if

     call mpi_bcast(ifinitialdone,1,mpi_logical,0, mpi_comm_world,ierr)
  end if

  kelvin0=kelvinarray(nowikelvin)
end subroutine

subroutine temd_collectpot()
  integer :: ikelvin,i,ierr,index
  integer,dimension(0:procnum-1) :: tpcountarray,tpposarray
  real*8 :: tpkelvin

! epot: instant energy
  if ( procid == 0 ) then
     tpcountarray(0:procnum-1)=1
     do i=0,procnum-1
        ikelvin=replicankelvin(i)
        tpposarray(i)=ikelvin
     end do
  end if
  call mpi_gatherv(epot, 1, mpi_double_precision, &
       kelvinpot, tpcountarray, tpposarray,mpi_double_precision,0,mpi_comm_world, ierr)

  if ( (montensample == 0) .and. (procnum<nkelvin) ) then
     call gpu_download_crd(nowcoor); call gpu_download_vel(nowvel); call gpu_download_frc(nowforce)

     if ( procid == 0 ) then
        tpcountarray(0:procnum-1)=3*numboxatom
        do i=0,procnum-1
           ikelvin = replicankelvin(i)
           index=allcoorindex(ikelvin)
           tpposarray(i)=index*3*numboxatom
        end do
     end if
     call mpi_gatherv(nowcoor, 3*numboxatom, mpi_double_precision, &
          allcoor, tpcountarray, tpposarray, mpi_double_precision,0,mpi_comm_world, ierr)
     call mpi_gatherv(nowvel, 3*numboxatom, mpi_double_precision, &
          allvel, tpcountarray, tpposarray, mpi_double_precision,0,mpi_comm_world, ierr)
     call mpi_gatherv(nowforce, 3*numboxatom, mpi_double_precision, &
          allforce, tpcountarray, tpposarray, mpi_double_precision,0,mpi_comm_world, ierr)
  end if
end subroutine

subroutine temd_distributekelvin()
  use mdin_ctrl_dat_mod, only : using_pme_potential
  integer,dimension(0:procnum-1) :: tpcountarray,tpposarray
  integer :: i,index,ikelvin,ierr

  if ( procid == 0 ) then
     tpcountarray(0:procnum-1)=1
     do i=0,procnum-1
        tpposarray(i)=i
     end do
  end if
  call mpi_scatterv(replicankelvin, tpcountarray, tpposarray, mpi_integer, &
             nowikelvin, 1, mpi_integer, 0, mpi_comm_world, ierr)

  if ( (montensample == 0) .and. (procnum<nkelvin) ) then
     if ( procid == 0 ) then
        tpcountarray(0:procnum-1)=3*numboxatom
        do i=0,procnum-1
           ikelvin = replicankelvin(i)
           index=allcoorindex(ikelvin)
           tpposarray(i)=index*3*numboxatom
        end do
     end if
     call mpi_scatterv(allcoor, tpcountarray, tpposarray, mpi_double_precision, &
          nowcoor, 3*numboxatom, mpi_double_precision,0,mpi_comm_world, ierr)
     call mpi_scatterv(allvel, tpcountarray, tpposarray, mpi_double_precision, &
          nowvel, 3*numboxatom, mpi_double_precision,0,mpi_comm_world, ierr)
     call mpi_scatterv(allforce, tpcountarray, tpposarray, mpi_double_precision, &
          nowforce, 3*numboxatom, mpi_double_precision,0,mpi_comm_world, ierr)
     if (using_pme_potential) call gpu_force_new_neighborlist()
     call gpu_upload_crd(nowcoor); call gpu_upload_vel(nowvel); call gpu_upload_frc(nowforce)
  end if
end subroutine

end module

subroutine ccj_interface_mpiinit(tpmytaskid,tpnumtasks)
  use moldata
  implicit none
  integer,intent(out) :: tpmytaskid,tpnumtasks
  integer :: ierr
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,procID,ierr)
  call mpi_comm_size(mpi_comm_world,procNum,ierr)
  tpmytaskid = procid; tpnumtasks = procnum
end subroutine

subroutine ccj_interface_getmolecules(tpnmol,tpmolnum,tpnumatom,tpcharge,tpnumres,tpresname,tpresindex,tpatomname,tpmass,tppbcbox)
  use mdin_ctrl_dat_mod, only  : ntb
  use moldata
  implicit none
  integer,intent(in) :: tpnmol,tpmolnum(tpnmol),tpnumatom,tpnumres,tpresindex(tpnumres+1)
  real*8,intent(inout) :: tpcharge(tpnumatom),tpmass(tpnumatom)
  real*8,intent(in) :: tppbcbox(3)
  integer :: ierr,i,ires,jres,imol,iatom
  real*8 :: temp
  character*4,intent(in) :: tpresname(tpnumres),tpatomname(tpnumatom)
  logical :: ifmdec
  namelist / task / iftemd, ifflood, ifsmd, ifmdec, ifeds

  iftemd=.false.; ifflood=.false.; ifsmd=.false.; ifmdec=.false.
  rewind(iomdin); read(iomdin,nml=task,iostat=ierr)
  if ( procid == 0 ) then
     ierr = system("rm -r procinfo"); ierr = system("mkdir procinfo")
     call getfreeunit(iosiminfo); open(unit=iosiminfo,file="siminfo.txt",action="write")
     flush(iosiminfo)
  end if

  nbox = ntb
  do i=1,tpnumres
     if ( tpresname(i) == "WAT " ) exit
  end do
  numboxres=tpnumres; numres=i-1; numatom=tpresindex(numres+1)-1

  numboxatom = tpnumatom; allocate(nowcoor(3,numboxatom),refcoor(3,numboxatom)); nowcoor=0.0d0; refcoor=0d0
  allocate(nowvel(3,numboxatom),nowforce(3,numboxatom)); nowvel=0.0d0; nowforce=0.0d0
  allocate(atomname(numboxatom)); atomname(1:numboxatom)(1:4)=tpatomname(1:numboxatom)(1:4)
  allocate(resname(numboxres),resindex(numboxres+1)); resname(1:numboxres)=tpresname(1:numboxres)
  resindex(1:numboxres+1)=tpresindex(1:numboxres+1)-1

  allocate(atommass(numboxatom)); atommass(1:numboxatom) = tpmass(1:numboxatom)

  allocate(resmolindex(numboxres)); resmolindex=0
  iatom=0; jres=1; imol=1
  do ires=1,numboxres
     iatom = iatom + resindex(ires+1) - resindex(ires)
     if ( iatom == sum(tpmolnum(1:imol)) ) then
        resmolindex(jres:ires) = imol
        jres = ires+1; imol = imol + 1
     end if
  end do
  nummol = tpnmol; allocate(molindex(nummol+1)); molindex=0
  do ires=2,numboxres
     if ( resmolindex(ires) /= resmolindex(ires-1) ) molindex(resmolindex(ires)) = ires-1
  end do
  molindex(nummol+1) = numboxres

  if ( ifmdec .eqv. .true. ) call moldata_setccharge_formdec(tpcharge)
  allocate(atomcharge(numboxatom)); atomcharge(1:numboxatom)=tpcharge(1:numboxatom)/18.2223

  if ( nbox > 0 ) boxsize = tppbcbox
end subroutine

subroutine ccj_interface_init(tpkelvin,tpdt,tpfreq,tpcoor)
  use nfe_lib_mod, only : remember_atom_names
  use nfe_colvar_mod, only : colvar_nlread,colvar_bootstrap,cv_min,cv_max
  use prmtop_dat_mod, only : atm_igraph, atm_mass
  use file_io_mod
  use file_io_dat_mod, only : mdin_name
  use mdin_ctrl_dat_mod, only : ntwx
  use moldata; use colvar
  use temd
  use flood
  use smd
  use eds
  implicit none

  integer,intent(in) :: tpfreq
  real*8,intent(in) :: tpdt,tpcoor(3,numboxatom),tpkelvin
  integer :: ierr,icv,ifind

  call remember_atom_names(atm_igraph)
  icv = 0; rewind(iomdin)
  do
     call nmlsrc('colvar', iomdin, ifind)
     if (ifind == 0 ) exit
     read(iomdin,'(a80)')
     icv = icv + 1
  end do
  nfencv = icv; allocate(nfecvs(nfencv),cvlowrange(nfencv),cvuprange(nfencv))
  icv = 1; rewind(iomdin)
  do icv=1,nfencv
     call colvar_nlread(iomdin,nfecvs(icv))
     cvlowrange(icv)=cv_min; cvuprange(icv)=cv_max
  end do
  rewind(iomdin)
  do icv = 1, nfencv; call colvar_bootstrap(nfecvs(icv), icv, atm_mass); end do

  deltatime = tpdt; samplestep=tpfreq; savetrajstep=ntwx; kelvin0 = tpkelvin; nowcoor = tpcoor
  call moldata_init(); call colvar_init()

  if ( iftemd .or. ifflood .or. ifsmd .or. ifeds ) ntwx=0

  if ( iftemd .eqv. .true. ) call temd_init()
  if ( ifflood .eqv. .true. ) call flood_init()
  if ( ifsmd .eqv. .true. ) call smd_init()
  if ( ifeds .eqv. .true. ) call eds_init()
  call random_seed(); call mpi_barrier(mpi_comm_world,ierr)
end subroutine

subroutine ccj_interface_updateforce(tpnstep,tpifpotupdate,tpepot,tpkelvin)
  use moldata
  use flood
  use smd
  use eds
  implicit none
  integer,intent(in) :: tpnstep
  logical,intent(in) :: tpifpotupdate
  real*8,intent(in) :: tpepot
  real*8,intent(inout) :: tpkelvin
  real*8 :: temp
  character*100 :: filename,tpfilename

  nowstep = tpnstep+1; nowtime = dble(nowstep)*deltatime*1d-3; epot = tpepot
  if ( tpifpotupdate .eqv. .false. ) epot=0.0d0

  if ( ifflood .eqv. .true. ) then
     call flood_updateforce()
  else if ( ifsmd .eqv. .true. ) then
     call smd_updateforce()
  end if

  if ( kelvin0 /= tpkelvin ) then
     temp = sqrt(kelvin0/tpkelvin); tpkelvin=kelvin0
     call gpu_scale_velocities(temp)
  end if
end subroutine

subroutine ccj_interface_simulate(tpnstep,tpifpotupdate,tpepot,tpkelvin)
  use moldata
  use temd
  use eds
  implicit none
  integer,intent(in) :: tpnstep
  logical,intent(in) :: tpifpotupdate
  real*8,intent(in) :: tpepot
  real*8,intent(inout) :: tpkelvin

  integer :: i,j,k,ierr
  real*8 :: temp
  character*100 :: filename,tpfilename

  nowstep = tpnstep; nowtime = dble(nowstep)*deltatime*1d-3
  if ( iftemd .eqv. .true. ) then
     epot = tpepot; if ( tpifpotupdate .eqv. .false. ) return
     call temd_updatecoor()
  else if ( ifeds .eqv. .true. ) then
     call eds_updatecoor()
  end if

  if ( kelvin0 /= tpkelvin ) then
     temp = sqrt(kelvin0/tpkelvin); tpkelvin=kelvin0
     call gpu_scale_velocities(temp)
  end if
end subroutine

