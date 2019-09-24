! References
!
! Combining the biased and unbiased sampling strategy into one convenient free energy calculation method.
! Haomiao Zhang, Qiankun Gong, Haozhe Zhang and Changjun Chen*
! J. Comput. Chem. 2019 (accepted)
!
! Adaptively biased molecular dynamics for free energy calculations.
! Volodymyr Babin, Christopher Roland, and Celeste Sagui*
! J. Chem. Phys. 2008, 128, 134101

module subspace
  use moldata; use colvar
  implicit none
  integer,parameter :: maxnss=10
  integer :: ssntotgrid,nss,ssngrid(maxnss)
  integer,dimension(maxnss) :: ss2cv
  real*8,dimension(maxnss) :: sslowrange,ssuprange,ssresolution,sscenter,ssrad
  real*8,dimension(:),allocatable :: freegriddata,hisgriddata
  real*8 :: ssslope
  character*3 :: ssshape
contains

subroutine subspace_init()
  use nfe_colvar_mod, only : colvar_is_periodic
  integer :: iss,ierr,system,iofile

  do iss=1,maxnss
     if ( ss2cv(iss) == 0 ) exit
  end do
  nss = iss - 1
  ssntotgrid=1; do iss=1,nss; ssntotgrid = ssntotgrid*ssngrid(iss); end do
  allocate(freegriddata(ssntotgrid),hisgriddata(ssntotgrid)); freegriddata(:)=0.0d0; hisgriddata(:)=0.0d0

  do iss=1,nss
     sslowrange(iss) = cvlowrange(ss2cv(iss))
     ssuprange(iss) = cvuprange(ss2cv(iss))
     if ( (ss2cv(iss) <= nfencv) .and. (colvar_is_periodic(nfecvs(ss2cv(iss))) .eqv. .true.) ) then
        sslowrange(iss) = -pi; ssuprange(iss)=pi
        ssresolution(iss) = 2d0*pi/dble(ssngrid(iss))
        ssuprange(iss) = sslowrange(iss) + dble(ssngrid(iss)-1)*ssresolution(iss)
     end if
  end do

  ssresolution(:)=(ssuprange(:)-sslowrange(:))/dble(ssngrid(:)-1)
  if ( procid == 0 ) then
     write(iosiminfo,"(a)") "Subspace info"
     do iss=1,nss
        write(iosiminfo,"(a,2i5,5x,a,2f8.3,5x,a,i5,5x,a,f8.3)")  &
              "ss axis ",iss,ss2cv(iss),"range ",sslowrange(iss),ssuprange(iss),"ngrid ",ssngrid(iss),"res ",ssresolution(iss)
     end do
     if ( ssslope /= 0d0 ) then
        write(iosiminfo,"(a,2f10.3)") "ss center ",sscenter(1:nss)
        write(iosiminfo,"(a,2f10.3)") "ss radius ",ssrad(1:nss)
        write(iosiminfo,"(a,f10.3,5x,a,a3)") "ss slope  ",ssslope,"ss shape: ",ssshape
     end if
     flush(iosiminfo)
  end if
end subroutine

subroutine subspace_writesurfacedata(freefilename,tpgriddata)
  character(len=*),intent(in) :: freefilename
  real*8,intent(in),optional :: tpgriddata(ssntotgrid)
  integer :: iofree,igrid,tpindex,idim,stack(nss),stride
  real*8 :: tpvec(nss)
  call getfreeunit(iofree); open(file=freefilename, unit=iofree, action="write")
  write(iofree,"(i6)") nss
  do igrid=1,nss
     write(iofree,"(2i6,2f10.3)") igrid,ssngrid(igrid),sslowrange(igrid),ssuprange(igrid)
  end do
  do igrid=1,ssntotgrid
     stride=ssntotgrid; tpindex=igrid
     do idim=nss,1,-1
        stride=stride/ssngrid(idim)
        stack(idim)=(tpindex-1)/stride + 1
        tpindex = tpindex - (stack(idim)-1)*stride
     end do
     tpvec(1:nss) = sslowrange(1:nss) + dble(stack(1:nss)-1)*ssresolution(1:nss)
     if ( .not. present(tpgriddata) ) then
        write(iofree,"(100f10.4)") freegriddata(igrid),0.0d0,tpvec(1:nss)
     else
        write(iofree,"(100f10.4)") freegriddata(igrid),tpgriddata(igrid),tpvec(1:nss)
     end if
  end do
  close(iofree)
end subroutine

subroutine ccj_subspace_readsurfacedata(tpfilename,tpgriddata)
  character(len=*),intent(in) :: tpfilename
  real*8,intent(out),optional :: tpgriddata(ssntotgrid)
  integer :: iofile,i,igrid,tpnss,tpssngrid(maxnss)
  real*8 :: tpvec(nss),temp
  logical :: ifexist

  inquire(file=tpfilename, exist=ifexist)
  if ( ifexist .eqv. .false. ) call errormsg("surface data file does not exist!"//trim(tpfilename))
  call getfreeunit(iofile); open(unit=iofile,file=tpfilename,action="read")
  do i=1,nss+1; read(iofile,*); end do
  if ( .not. present(tpgriddata) ) then
     do igrid=1,ssntotgrid
        read(iofile,"(100f10.4)") freegriddata(igrid), temp, tpvec(1:nss)
     end do
  else
     do igrid=1,ssntotgrid
        read(iofile,"(100f10.4)") freegriddata(igrid), tpgriddata(igrid), tpvec(1:nss)
     end do
  end if
  close(iofile)
end subroutine

subroutine subspace_getgridvec(gridindex,gridvec)
  integer,intent(in) :: gridindex
  integer,intent(out) :: gridvec(nss)
  integer :: idim,stride,tpindex

  stride=ssntotgrid; tpindex=gridindex
  do idim=nss,1,-1
     stride=stride/ssngrid(idim)
     gridvec(idim)=(tpindex-1)/stride + 1
     if ( (gridvec(idim)<1) .or. (gridvec(idim)>ssngrid(idim)) ) call errormsg("gridindex error!")
     tpindex = tpindex - (gridvec(idim)-1)*stride
  end do
end subroutine

subroutine subspace_getgridindex(gridvec,gridindex)
  integer,intent(in) :: gridvec(nss)
  integer,intent(out) :: gridindex
  integer :: i,tpstride

  gridindex=0
  do i=1,nss
     if ( (gridvec(i)<1) .or. (gridvec(i)>ssngrid(i)) ) then
        call errormsg("ssubspace: illegal gridindex ",int1=gridvec(i),int2=ssngrid(i))
     end if
  end do
  tpstride=1; gridindex=gridvec(1)
  do i=2,nss
     tpstride=tpstride*ssngrid(i-1)
     gridindex = gridindex + tpstride*(gridvec(i)-1)
  end do
end subroutine

subroutine subspace_getowngridandbias(tpsspos, owngrid, localbias, closestgrid)
  real*8,intent(in) :: tpsspos(nss)
  integer,intent(out) :: owngrid(nss)
  real*8,intent(out) :: localbias(nss)
  integer,intent(out),optional :: closestgrid
  integer :: iss,igrid,stride,tpgridindex,stack(nss),tpgridvec(nss)

  owngrid(:)=0; localbias(:)=0.0d0
  do iss=1,nss
     owngrid(iss) = int((tpsspos(iss) - sslowrange(iss))/ssresolution(iss)) + 1
     localbias(iss) = (tpsspos(iss) - sslowrange(iss))/ssresolution(iss) - dble(owngrid(iss) - 1)
  end do
end subroutine

subroutine subspace_getsurfaceinfo(tpsspos,tpsurfacedata,surfacevalue,surfacedev)
  use nfe_colvar_mod, only : colvar_is_periodic
  real*8,intent(in) :: tpsspos(nss),tpsurfacedata(ssntotgrid)
  real*8,intent(out) :: surfacevalue
  real*8,intent(out),optional :: surfacedev(nss)
  integer :: tpindex,iss,i,stride,igrid,owngrid(nss)
  integer :: stack(nss+1),tpstack(nss),tpssntotgrid
  real*8 :: cubicvec(4,nss+1),cubicdev(4,nss+1,nss),temp,localbias(nss)

  call subspace_getowngridandbias(tpsspos, owngrid, localbias)
  surfacevalue=0.0d0; cubicvec(:,:)=0.0d0; stack(:)=1
  if ( present(surfacedev) ) surfacedev=0.0d0
  tpssntotgrid=4**nss
  do igrid=1,tpssntotgrid
     stride=tpssntotgrid; tpindex=igrid
     do iss=nss,1,-1
        stride=stride/4
        stack(iss)=(tpindex-1)/stride + 1
        tpindex = tpindex - (stack(iss)-1)*stride
     end do

     tpstack(1:nss)=owngrid(1:nss)+stack(1:nss)-2
     do iss=1,nss
        if ( (ss2cv(iss) <= nfencv) .and. (colvar_is_periodic(nfecvs(ss2cv(iss))) .eqv. .true.) ) then
           if ( tpstack(iss) < 1 ) then
              tpstack(iss) = tpstack(iss) + ssngrid(iss)
           else if ( tpstack(iss) > ssngrid(iss) ) then
              tpstack(iss) = tpstack(iss) - ssngrid(iss)
           end if
        else
           if ( tpstack(iss) == 0 ) then
              tpstack(iss) = 1
           else if ( (tpstack(iss)==(ssngrid(iss)+1)) .or. (tpstack(iss)==(ssngrid(iss)+2)) ) then
              tpstack(iss) = ssngrid(iss)
           end if
        end if
     end do

     call subspace_getgridindex(tpstack,tpindex)
     cubicvec(stack(1),1)=tpsurfacedata(tpindex)

     do iss=1,nss
        if ( stack(iss) == 4 ) then
           if ( present(surfacedev) ) then
              call getcubicbspline(localbias(iss),cubicvec(1:4,iss), &
                           cubicvec(stack(iss+1),iss+1),cubicdev(stack(iss+1),iss+1,iss) )
              do i=1,iss-1
                 call getcubicbspline(localbias(iss),cubicdev(1:4,iss,i),cubicdev(stack(iss+1),iss+1,i))
                 cubicdev(:,iss,i)=0.0d0
              end do
           else
              call getcubicbspline(localbias(iss),cubicvec(1:4,iss),cubicvec(stack(iss+1),iss+1))
           end if

           if ( iss == nss ) then
              surfacevalue=cubicvec(stack(iss+1),iss+1)
              if ( present(surfacedev) ) surfacedev(1:nss)=cubicdev(stack(iss+1),iss+1,1:nss)
           end if
           cubicvec(:,iss)=0.0d0;
        else
           exit
        end if
     end do
  end do

  if ( present(surfacedev) ) surfacedev(1:nss)=surfacedev(1:nss)/ssresolution(1:nss)
end subroutine

subroutine getcubicbspline(bias,gridvec,splinedvalue,splineddev)
  implicit none
  real*8,intent(in) :: bias,gridvec(4)
  real*8,intent(out) :: splinedvalue
  real*8,intent(out),optional :: splineddev
  real*8 :: xi

  splinedvalue=0.0d0; if ( present(splineddev) ) splineddev=0.0d0

  xi = 1.0d0 + bias
  splinedvalue = splinedvalue + gridvec(1)*((2.0d0-xi)**3)/6.0d0
  if ( present(splineddev) ) splineddev = splineddev - gridvec(1)*0.5d0*((2.0d0-xi)**2)

  xi = bias
  splinedvalue = splinedvalue + gridvec(2)*( (xi**2)*(xi-2.0d0)/2.0d0 + 2.0d0/3.0d0 )
  if ( present(splineddev) ) splineddev = splineddev + gridvec(2)*( 1.5d0*xi*xi - 2.0d0*xi )

  xi = bias - 1.0d0
  splinedvalue = splinedvalue + gridvec(3)*( (xi**2)*(-xi-2.0d0)/2.0d0 + 2.0d0/3.0d0 )
  if ( present(splineddev) ) splineddev = splineddev + gridvec(3)*( -1.5d0*xi*xi - 2.0d0*xi  )

  xi = bias - 2.0d0
  splinedvalue = splinedvalue + gridvec(4)*((2.0d0+xi)**3)/6.0d0
  if ( present(splineddev) ) splineddev = splineddev + gridvec(4)*0.5d0*((2.0d0+xi)**2)

end subroutine

subroutine subspace_hisgriddata2freegriddata(tphisgriddata,tpfreegriddata)
  real*8,intent(in) :: tphisgriddata(ssntotgrid)
  real*8,intent(out) :: tpfreegriddata(ssntotgrid)
  integer :: igrid
  real*8 :: surfaceheight,temp

  temp = sum(tphisgriddata); tpfreegriddata = tphisgriddata/temp
  do igrid=1,ssntotgrid
     if ( tpfreegriddata(igrid) > 1.0d-10 ) then
        tpfreegriddata(igrid) = -gasconst*kelvin0*log(tpfreegriddata(igrid))
     else
        tpfreegriddata(igrid) = 0.0d0
     end if
  end do

  surfaceheight=0.0d0
  do igrid=1,ssntotgrid
     if ( tpfreegriddata(igrid) == 0.0d0 ) cycle
     if ( surfaceheight == 0.0d0 ) surfaceheight=tpfreegriddata(igrid)
     if ( tpfreegriddata(igrid) > surfaceheight ) surfaceheight=tpfreegriddata(igrid)
  end do
  do igrid=1,ssntotgrid
     if ( tpfreegriddata(igrid) /= 0.0d0 ) tpfreegriddata(igrid) = tpfreegriddata(igrid)-surfaceheight
  end do
end subroutine

end module

module flood
  use moldata
  use subspace
  implicit none

  integer,private :: moviecv
  integer,private :: nkelvin,exchangestep,biasrep,biasexc
  real*8,private :: floodingtime,biasweight
  real*8,dimension(:),allocatable,private :: iniaddgriddata,addbiasgriddata,prebiasgriddata,totbiasgriddata
  real*8,dimension(:),allocatable,private :: kelvinarray,biasarray,wellkelvinarray
  real*8,dimension(:,:),allocatable,private :: kelvinrepdata
  real*8,private :: wellkelvin,kelvindown,kelvinup
  integer,private :: savesurfacefreq
  logical,private :: iffloodeds
contains

subroutine flood_init()
  use eds, only : eds_init
  include 'mpif.h'
  integer :: ierr,iofile,i,j,k,system
  real*8 :: biasratio,tpbias,temp
  logical :: ifexist
  integer,dimension(0:procnum-1) :: tpcountarray,tpposarray
  namelist /flood/ floodingtime, savesurfacefreq, wellkelvin, moviecv, &
                   ss2cv, ssngrid, &
                   kelvindown, kelvinup, biasratio, exchangestep, &
                   ssshape, sscenter, ssrad, ssslope, biasrep, biasexc, iffloodeds

  if ( procid == 0 ) write(iosiminfo,"(a)") "Perform mixing REMD simulation"

  savesurfacefreq=1000000; wellkelvin=5000.0; biasratio=1.0d0 
  read(iomdin,nml=flood,iostat=ierr); if ( ierr < 0 ) call errormsg("namelist flood error!")
  if ( floodingtime == 0d0 ) wellkelvin=0d0
  call subspace_init()
  allocate(addbiasgriddata(ssntotgrid),iniaddgriddata(ssntotgrid)); addbiasgriddata=0.0d0; iniaddgriddata=0.0d0

  if ( procid == 0 ) then
     write(iosiminfo,*)
     write(iosiminfo,"(a,f10.3)") "floodingtime: ",floodingtime
     write(iosiminfo,"(a,i10,5x,a,i5)") "savesurfacefreq ",savesurfacefreq,"save movie for cv ",moviecv
     write(iosiminfo,*)
     flush(iosiminfo)
     allocate(totbiasgriddata(ssntotgrid),prebiasgriddata(ssntotgrid))
     totbiasgriddata=0.0d0; prebiasgriddata=0.0d0; addbiasgriddata=0.0d0

     inquire(file="./surfacedata_ref.txt", exist=ifexist)
     if ( ifexist .eqv. .true. ) then
        write(iosiminfo,"(a)") "reading reference biasgriddata"
        call ccj_subspace_readsurfacedata("./surfacedata_ref.txt",tpgriddata=addbiasgriddata)
     end if
     if ( ssslope /= 0d0 ) then
        call flood_setInitialBoundPotential(); 
        call subspace_writesurfacedata("procinfo/inibiasdata.txt",tpgriddata=iniaddgriddata)
     end if
     addbiasgriddata=addbiasgriddata+iniaddgriddata; prebiasgriddata=addbiasgriddata

     allocate(replicankelvin(0:procnum-1)); do i=0,procnum-1; replicankelvin(i)=i; end do
     tpcountarray=1; do i=0,procnum-1; tpposarray(i)=i; end do
  end if

  call mpi_scatterv(replicankelvin, tpcountarray, tpposarray, mpi_integer, &
          nowikelvin, 1, mpi_integer, 0, mpi_comm_world, ierr)
  call mpi_bcast(iniaddgriddata, ssntotgrid , mpi_double_precision, 0, mpi_comm_world, ierr)
  call mpi_bcast(addbiasgriddata, ssntotgrid , mpi_double_precision, 0, mpi_comm_world, ierr)

  if ( procid == 0 ) then
     allocate(kelvinarray(0:procnum-1),biasarray(0:procnum-1)); kelvinarray=kelvindown; biasarray=1.0d0; nkelvin=procnum
     allocate(kelvinrepdata(ssntotgrid,0:procnum-1)); kelvinrepdata=0.0d0
     allocate(wellkelvinarray(0:procnum-1)); wellkelvinarray=wellkelvin

     write(iosiminfo,"(a,i6,5x,a,i10)") "biased replicas: ",biasrep,"copy steps: ",biasexc*exchangestep
     if ( biasrep > procnum ) call errormsg("biasrep cannot be larger than procnum")

     wellkelvinarray = 0d0; biasarray=0d0
     if ( biasrep > 0 ) then
        temp = wellkelvin/dble(biasrep)
        do i=1,biasrep
           wellkelvinarray(procnum-1-biasrep+i) = dble(i)*temp
        end do
        if ( biasratio < 0d0 ) then
           do i=1,biasrep
              biasarray(procnum-1-biasrep+i) = 1.0d0
           end do
        else if ( biasratio == 0d0 ) then
           biasarray(:)=0.0d0
        else if ( biasratio == 1.0d0 ) then
           temp = 1.0d0/dble(biasrep)
           do i=1,biasrep
              biasarray(procnum-1-biasrep+i) = dble(i)*temp
           end do
        else if ( biasratio /= 1.0d0 ) then
           biasarray=0d0; tpbias=0.0d0
           temp = (1.0d0 - 0.0d0)*(1.0d0-biasratio)/(1.0d0-biasratio**(dble(biasrep)))
           do i=1,biasrep
              tpbias = tpbias + temp*(biasratio**(dble(i-1)))
              biasarray(procnum-1-biasrep+i) = tpbias
           end do
           biasarray(procnum-1) = 1.0d0
        end if
     end if
     do i=0,procnum-1
        if ( floodingtime > 0d0 ) then
           if ( (biasarray(i) /= 1.0d0) .and. (biasarray(i) /= 0.0d0) ) then
              call errormsg("error biasweight when flooding > 0 ",real1=biasarray(i))
           end if
        else if ( floodingtime == 0d0 ) then
           if ( wellkelvinarray(i) /= 0.0d0 ) call errormsg("error wellkelvin when floodingtime==0 ",real1=wellkelvinarray(i))
        else
           call errormsg("error flooding time! ",real1=floodingtime)
        end if
     end do 
     kelvinarray = kelvindown; j=-1 
     if ( kelvinup > kelvindown ) then
        if ( (biasexc == 0) .or. (floodingtime==0d0) ) then
           if (biasrep < (procnum-1) ) then
              j = procnum-1-biasrep
           end if
        else
           if (biasrep < (procnum-2) ) then
              j = procnum-1-biasrep-1
           end if
        end if
        temp = log(kelvinup/kelvindown)
        do i=0,j
           kelvinarray(i)=kelvindown*exp(temp*dble(i)/dble(j))
        end do
        do i=j+1,procnum-1
           kelvinarray(i) = kelvinup
        end do
     end if

     write(iosiminfo,"(a,i10)") "exchange step: ",exchangestep
     write(iosiminfo,"(a,100f10.3)") "Temperatures: ",kelvinarray(:)
     write(iosiminfo,"(a,100f10.3)") "Bias weights: ",biasarray(:)
     write(iosiminfo,"(a,100f10.3)") "Well-kelvin:  ",wellkelvinarray(:)
     write(iosiminfo,"(a,l10)") "iffloodeds ",iffloodeds
     flush(iosiminfo)
     tpcountarray(0:procnum-1)=1
     do i=0,procnum-1; tpposarray(i)=replicankelvin(i); end do
  end if

  call mpi_scatterv(biasarray, tpcountarray, tpposarray, mpi_double_precision, &
            biasweight, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
  call mpi_scatterv(kelvinarray, tpcountarray, tpposarray, mpi_double_precision, &
            kelvin0, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
  call mpi_scatterv(wellkelvinarray, tpcountarray, tpposarray, mpi_double_precision, &
            wellkelvin, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
  call mpi_bcast(nkelvin,1,mpi_integer,0, mpi_comm_world,ierr)

  if ( iffloodeds .eqv. .true. ) call eds_init()
end subroutine

subroutine flood_updateforce()
  use eds, only : eds_updatecoor,eds_flood_readbox
  use nfe_colvar_mod, only : colvar_force,colvar_is_periodic
  use mdin_ctrl_dat_mod, only : using_pme_potential,ntpr
  use energy_records_mod, only : gb_pot_ene_rec
  integer,dimension(:),allocatable,save :: jumpsum,jumparray
  real*8,dimension(:),allocatable,save :: lastfreegriddata
  integer,save :: icalfree, iomovie,ioexchangerate,ionegrec
  real*8,save :: toprepkelvin
  integer :: i,iss,icv,iproc,igrid,ierr,iatom,tpgridvec(nss),ikelvin,jkelvin,irep,jrep,status
  real*8 :: surfacedev(nss),cvvec(ncv),sspos(nss),tpheight,tpvec(nss)
  real*8 :: tpgriddata(ssntotgrid),itererror,tpavg1,tpavg2,temp,pab,xrand,tprepsspos(nss,0:procnum-1)
  real*8 :: tpbiasioni,tpbiasionj,tpbiasjonj,tpbiasjoni,tpkelvin,tpepot,tpenmr(3)
  logical :: ifexceed,ifexchange,ifreadeds
  character*100 :: tpfilename
  integer,dimension(0:procnum-1) :: tpcountarray,tpposarray,kelvinreplica
  real*8,dimension(0:procnum-1) :: tpreppot,tprepbias,exchangerate
  real*8,dimension(3,numboxatom) :: tpcoor,tpvel,tpfrc
  type(gb_pot_ene_rec)          :: pot_ene

  call moldata_download_crd(nowcoor)

  if ( iffloodeds .eqv. .true. ) call eds_updatecoor()
  call colvar_calandrec(cvvec,tpmoviecv=moviecv)

  ifexceed=.false.
  do iss=1,nss
     icv = ss2cv(iss)
     if ( icv > 0 ) then
        sspos(iss) = cvvec(icv)
     else if ( icv == -1 ) then
        if ( numatom /= numboxatom ) call errormsg("icv cannot be -1 in explicit solvent!")
        call gpu_download_frc(tpfrc); sspos(iss) = epot; ntpr=1
        if ( epot == 0.0d0 ) then
           call gpu_gb_ene(pot_ene, tpenmr); sspos(iss)=pot_ene%total
        end if
     else
        call errormsg("ss2cv index error ",int1=iss, int2=icv)
     end if 
     if ( ( (ss2cv(iss) > nfencv) .or. (colvar_is_periodic(nfecvs(ss2cv(iss))) .eqv. .false.) ) &
               .and. ((sspos(iss) < sslowrange(iss)) .or. (sspos(iss) > ssuprange(iss)))  ) then
        ifexceed=.true.
        if ( iffloodeds .eqv. .false. ) then
           write (iosiminfo,"(a,i12,i4,3f15.3)") "ss exceed ",nowstep,iss,sspos(iss), &
                   sslowrange(iss),ssuprange(iss)
           stop
        end if
     end if
  end do

  if ( ( biasweight > 0d0 ) .and. ( iffloodeds .eqv. .false. ) ) then
     call subspace_getsurfaceinfo(sspos,addbiasgriddata,tpheight,surfacedev=surfacedev)
     nowforce=0.0d0; surfacedev=surfacedev*biasweight; tpheight=tpheight*biasweight
     do iss = 1, nss
        icv = ss2cv(iss)
        if ( icv > 0 ) then
           temp = -surfacedev(iss)
           if ( icv <= nfencv ) then
              call colvar_force(nfecvs(icv), nowcoor, temp, nowforce)
           else
              call colvar_fsainfo(icv-nfencv, nowcoor, xrand, tpgrad=nowforce, tpfactor=temp);
           end if
        else if ( icv == -1 ) then
           nowforce = nowforce + surfacedev(iss)*tpfrc
        else
           call errormsg("error ss2cv")
        end if
     end do
     call moldata_upload_frc_add(nowforce)
  end if

  if ( ifexceed .eqv. .false. ) then
     call subspace_getowngridandbias(sspos,tpgridvec,tpvec)
     call subspace_getgridindex(tpgridvec,igrid)

     if ( (floodingtime == 0.0d0) .or. (biasweight == 0d0) .or. (iffloodeds .eqv. .true.)) then
        hisgriddata(igrid) = hisgriddata(igrid) + 1.0d0
     else
        call flood_updatebiaspotentialgriddata(sspos,tpheight=tpheight)
     end if
  end if

  if ( (exchangestep>0) .and. (mod(nowstep,exchangestep)==0) .and. (nkelvin>1) ) then
     if ( procid == 0 ) then
        tpcountarray(0:procnum-1)=1
        do iproc=0,procnum-1; tpposarray(iproc)=iproc; end do
     end if
     call mpi_gatherv(epot, 1, mpi_double_precision, &
             tpreppot, tpcountarray, tpposarray,mpi_double_precision,0,mpi_comm_world, ierr)

     if ( procid == 0 ) then
        tpcountarray(0:procnum-1)=nss
        do iproc=0,procnum-1; tpposarray(iproc)=iproc*nss; end do
     end if
     call mpi_gatherv(sspos, nss, mpi_double_precision, &
             tprepsspos, tpcountarray, tpposarray, mpi_double_precision,0,mpi_comm_world, ierr)

     if ( procid == 0 ) then
        tpcountarray(0:procnum-1)=ssntotgrid
        do iproc=0,procnum-1; tpposarray(iproc)=replicankelvin(iproc)*ssntotgrid; end do
     end if
     if ( floodingtime > 0.0d0 ) then
        if ( (biasweight == 0.0d0) .or. (iffloodeds .eqv. .true.) ) then
           tpgriddata=hisgriddata
        else
           tpgriddata=addbiasgriddata
        end if
        call mpi_gatherv(tpgriddata, ssntotgrid, mpi_double_precision, &
             kelvinrepdata, tpcountarray, tpposarray, mpi_double_precision,0,mpi_comm_world, ierr)
     else if ( floodingtime == 0.0d0 ) then
        call mpi_gatherv(hisgriddata, ssntotgrid, mpi_double_precision, &
             kelvinrepdata, tpcountarray, tpposarray, mpi_double_precision,0,mpi_comm_world, ierr)
     end if

     if ( procid == 0 ) then

        if ( .not. allocated(jumpsum) ) then
           allocate(jumpsum(0:nkelvin-1),jumparray(0:nkelvin-1)); jumparray=0; jumpsum=0
        end if

        kelvinreplica=-1
        do iproc=0,procnum-1
           ikelvin=replicankelvin(iproc); kelvinreplica(ikelvin)=iproc
           if ( floodingtime==0.0d0 ) call subspace_getsurfaceinfo(tprepsspos(1:nss,iproc),addbiasgriddata,tprepbias(iproc))
        end do

        do ikelvin=0,nkelvin-2
           jkelvin = ikelvin + 1
           irep=kelvinreplica(ikelvin); jrep=kelvinreplica(jkelvin)

           pab = - (tpreppot(irep) - tpreppot(jrep))/(gasconst*kelvinarray(jkelvin))  &
                 - (tpreppot(jrep) - tpreppot(irep))/(gasconst*kelvinarray(ikelvin))

           if ( (floodingtime > 0.0d0) .and. (biasarray(ikelvin) > 0.0d0) ) then

              call subspace_getsurfaceinfo(tprepsspos(1:nss,irep),kelvinrepdata(:,ikelvin),tpbiasioni)
              call subspace_getsurfaceinfo(tprepsspos(1:nss,jrep),kelvinrepdata(:,jkelvin),tpbiasjonj)
              call subspace_getsurfaceinfo(tprepsspos(1:nss,irep),kelvinrepdata(:,jkelvin),tpbiasionj)
              call subspace_getsurfaceinfo(tprepsspos(1:nss,jrep),kelvinrepdata(:,ikelvin),tpbiasjoni)

              pab = pab - (tpbiasionj - tpbiasjonj)/(gasconst*kelvinarray(jkelvin)) &
                        - (tpbiasjoni - tpbiasioni)/(gasconst*kelvinarray(ikelvin))

           else if ( floodingtime == 0.0d0 ) then

             pab = pab - biasarray(jkelvin)*(tprepbias(irep) - tprepbias(jrep))/(gasconst*kelvinarray(jkelvin)) &
                       - biasarray(ikelvin)*(tprepbias(jrep) - tprepbias(irep))/(gasconst*kelvinarray(ikelvin))

           end if

           pab = exp(pab)

           if ( (floodingtime>0d0) .and. (biasrep>0) .and. (biasexc>0) ) then
              if ( ikelvin == (nkelvin-1-biasrep) ) then 
                 pab = 0.0d0
              else if ( ikelvin == (nkelvin-1-biasrep-1) ) then
                 pab = exp(-(tpreppot(jrep)-tpreppot(irep))/(gasconst*kelvinarray(ikelvin)))
              end if
           end if

           ifexchange=.false.
           if ( pab >= 1.0d0 ) then
              ifexchange=.true.
           else
              call random_number(xrand)
              if ( pab > xrand ) ifexchange=.true.
           end if

           jumpsum(jkelvin) = jumpsum(jkelvin) + 1
           if ( ifexchange .eqv. .true. ) then
              jumparray(jkelvin) = jumparray(jkelvin) + 1 
              kelvinreplica(ikelvin)=jrep; kelvinreplica(jkelvin)=irep
           end if
        end do

        if ( (mod(nowstep,exchangestep*100)==0) ) then
           if ( ioexchangerate == 0 ) then
              call getfreeunit(ioexchangerate)
              open(file="procinfo/kelvin_exchangerate.txt", unit=ioexchangerate, action="write")
           end if
           exchangerate=0d0
           do ikelvin=1,nkelvin-1
              if ( jumpsum(ikelvin) > 0 ) exchangerate(ikelvin)=dble(jumparray(ikelvin))/dble(jumpsum(ikelvin))
           end do
           temp = sum(exchangerate)/dble(nkelvin-1)
           write(ioexchangerate,"(f10.1,50f8.4)") deltatime*dble(nowstep)/1000d0,temp,exchangerate(1:nkelvin-1)
           flush(ioexchangerate)
        end if

        do ikelvin=0,nkelvin-1
           ierr = kelvinreplica(ikelvin); replicankelvin(ierr)=ikelvin
        end do

     end if

     if ( procid == 0 ) then
        tpcountarray(0:procnum-1)=1
        do iproc=0,procnum-1; tpposarray(iproc)=replicankelvin(iproc); end do
     end if
     call mpi_scatterv(biasarray, tpcountarray, tpposarray,mpi_double_precision, &
             biasweight, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
     temp = kelvin0
     call mpi_scatterv(kelvinarray, tpcountarray, tpposarray,mpi_double_precision, &
             kelvin0, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
     if ( (iffloodeds .eqv. .true.) .and. (procid == (procnum-1)) ) kelvin0=temp
     if ( floodingtime > 0.0d0 ) then
        call mpi_scatterv(wellkelvinarray, tpcountarray, tpposarray,mpi_double_precision, &
             wellkelvin, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
     end if

     if ( procid == 0 ) then
        tpcountarray(0:procnum-1)=ssntotgrid
        do iproc=0,procnum-1; tpposarray(iproc)=replicankelvin(iproc)*ssntotgrid; end do
     end if
     if ( floodingtime > 0d0 ) then
        call mpi_scatterv(kelvinrepdata, tpcountarray, tpposarray, mpi_double_precision, &
             tpgriddata, ssntotgrid, mpi_double_precision, 0, mpi_comm_world, ierr)
        if ( (biasweight == 0.0d0) .or. (iffloodeds .eqv. .true.) ) then
           hisgriddata=tpgriddata
        else
           addbiasgriddata=tpgriddata
        end if
     else if ( floodingtime == 0d0 ) then
        call mpi_scatterv(kelvinrepdata, tpcountarray, tpposarray,mpi_double_precision, &
             hisgriddata, ssntotgrid, mpi_double_precision, 0,mpi_comm_world, ierr)
     end if

     if ( procid == 0 ) then
        tpcountarray(0:procnum-1)=1
        do iproc=0,procnum-1; tpposarray(iproc)=iproc; end do
     end if
     call mpi_scatterv(replicankelvin, tpcountarray, tpposarray, mpi_integer, &
             nowikelvin, 1, mpi_integer, 0, mpi_comm_world, ierr)

     if ( (biasrep>0) .and. (biasexc>0) .and. (mod(nowstep,exchangestep*biasexc)==0) ) then
        call colvar_calandrec(cvvec,ifonlycv=.true.)

        i=procnum; call mpi_bcast(kelvinreplica(0:i-1), i, mpi_integer, 0, mpi_comm_world, ierr)
        jkelvin = nkelvin - 1 - biasrep; ikelvin = jkelvin + 1
        irep=kelvinreplica(ikelvin); jrep=kelvinreplica(jkelvin)

        if ( iffloodeds .eqv. .false. ) then
           if ( procid == irep ) then
              call gpu_download_crd(tpcoor); call gpu_download_vel(tpvel); call gpu_download_frc(tpfrc)
              i=jrep; call mpi_send(tpcoor,3*numboxatom,mpi_double_precision,i,91,mpi_comm_world,ierr)
              i=jrep; call mpi_send(tpvel,3*numboxatom,mpi_double_precision,i,92,mpi_comm_world,ierr)
              i=jrep; call mpi_send(tpfrc,3*numboxatom,mpi_double_precision,i,93,mpi_comm_world,ierr)
           else if ( procid == jrep ) then
              i=irep; call mpi_recv(nowcoor,3*numboxatom,mpi_double_precision,i,91,mpi_comm_world,status,ierr)
              i=irep; call mpi_recv(nowvel,3*numboxatom,mpi_double_precision,i,92,mpi_comm_world,status,ierr)
              i=irep; call mpi_recv(nowforce,3*numboxatom,mpi_double_precision,i,93,mpi_comm_world,status,ierr)
              call gpu_upload_crd(nowcoor); call gpu_upload_vel(nowvel); call gpu_upload_frc(nowforce)
              if (using_pme_potential) call gpu_force_new_neighborlist()
           end if
        else
           if ( irep /= (procnum-1) ) call errormsg("floodeds error!")
           if (procid == (procnum-1)) call eds_flood_readbox(ifreadeds,tpcoor,tpvel,tpfrc)
           call mpi_bcast(ifreadeds, 1, mpi_logical, procnum-1, mpi_comm_world, ierr)
           if ( ifreadeds .eqv. .true. ) then
              if ( procid == irep ) then
                 i=jrep; call mpi_send(tpcoor,3*numboxatom,mpi_double_precision,i,91,mpi_comm_world,ierr)
                 i=jrep; call mpi_send(tpvel,3*numboxatom,mpi_double_precision,i,92,mpi_comm_world,ierr)
                 i=jrep; call mpi_send(tpfrc,3*numboxatom,mpi_double_precision,i,93,mpi_comm_world,ierr)
              else if ( procid == jrep ) then
                 i=irep; call mpi_recv(nowcoor,3*numboxatom,mpi_double_precision,i,91,mpi_comm_world,status,ierr)
                 i=irep; call mpi_recv(nowvel,3*numboxatom,mpi_double_precision,i,92,mpi_comm_world,status,ierr)
                i=irep; call mpi_recv(nowforce,3*numboxatom,mpi_double_precision,i,93,mpi_comm_world,status,ierr)
                 call gpu_upload_crd(nowcoor); call gpu_upload_vel(nowvel); call gpu_upload_frc(nowforce)
                 if (using_pme_potential) call gpu_force_new_neighborlist()
              end if
           end if
        end if
     end if
  end if

  if ( mod(nowstep,savesurfacefreq) == 0 ) then
     if ( floodingtime > 0.0d0 ) then
        if ( exchangestep == 0 ) then
           call mpi_reduce(addbiasgriddata, totbiasgriddata, ssntotgrid,  &
                 mpi_double_precision, MPI_SUM, 0, mpi_comm_world, ierr)
           if ( procid == 0 ) then
              addbiasgriddata = totbiasgriddata - dble(procnum-1)*prebiasgriddata
              prebiasgriddata = addbiasgriddata
           end if
           call mpi_bcast(addbiasgriddata, ssntotgrid , mpi_double_precision, 0, mpi_comm_world, ierr)
           if ( procid == 0 ) then
              call flood_generatefreesurface()
              call subspace_writesurfacedata("procinfo/surfacedata.txt",tpgriddata=addbiasgriddata-iniaddgriddata)
           end if
        else
           call flood_generatefreesurface()
           write(tpfilename,"(i8)") nowikelvin; tpfilename=adjustl(tpfilename)
           tpfilename="procinfo/surfacedata_"//trim(tpfilename)//".txt"
           if ( biasweight == 0.0d0 ) then
              temp = sum(hisgriddata); tpgriddata=hisgriddata/temp
              call subspace_writesurfacedata(tpfilename,tpgriddata=tpgriddata)
           else
              call subspace_writesurfacedata(tpfilename,tpgriddata=addbiasgriddata-iniaddgriddata)
           end if
        end if

     else if ( floodingtime == 0.0d0 ) then

        if ( exchangestep == 0 ) then
           call mpi_reduce(hisgriddata, totbiasgriddata, ssntotgrid,  &
                 mpi_double_precision, MPI_SUM, 0, mpi_comm_world, ierr)
           if ( procid == 0 ) then
              hisgriddata = totbiasgriddata - dble(procnum-1)*prebiasgriddata
              prebiasgriddata = hisgriddata
           end if
           call mpi_bcast(hisgriddata, ssntotgrid , mpi_double_precision, 0, mpi_comm_world, ierr)
           if ( procid == 0 ) then
              call flood_generatefreesurface()
              call subspace_writesurfacedata("procinfo/surfacedata.txt",tpgriddata=(addbiasgriddata-iniaddgriddata)*biasweight)
           end if
        else
           call flood_generatefreesurface()
           write(tpfilename,"(i8)") nowikelvin; tpfilename=adjustl(tpfilename)
           tpfilename="procinfo/surfacedata_"//trim(tpfilename)//".txt"
           if ( biasweight == 0.0d0 ) then
              temp = sum(hisgriddata); tpgriddata=hisgriddata/temp
              call subspace_writesurfacedata(tpfilename,tpgriddata=tpgriddata)
           else
              call subspace_writesurfacedata(tpfilename,tpgriddata=(addbiasgriddata-iniaddgriddata)*biasweight)
           end if
        end if

     end if

!     icalfree = icalfree + 1
!     if ( (nowikelvin == 0) .and. (mod(icalfree,100) == 0) ) then
!        write(tpfilename,"(i8)") nint(nowstep*deltatime/1000.0d0); tpfilename=adjustl(tpfilename)
!        tpfilename="procinfo/surfacedata_"//trim(tpfilename)//"ns.txt"
!        call subspace_writesurfacedata(tpfilename,tpgriddata=addbiasgriddata-iniaddgriddata)
!     end if
     call mpi_barrier(mpi_comm_world,ierr)
  end if
end subroutine

subroutine flood_setInitialBoundPotential()
  integer :: i,j,tpindex,iss,stride,igrid
  integer :: stack(nss)
  real*8 :: tpvec(nss),tpratiovec(nss),tpmarkdis,tpdis
  logical :: ifvalidgrid

  if ( procid == 0 ) write(iosiminfo,"(a)") "set up initial bias potential"
  iniaddgriddata=0.0d0
  do igrid=1,ssntotgrid
     stride=ssntotgrid; tpindex=igrid
     do iss=nss,1,-1
        stride=stride/ssngrid(iss)
        stack(iss)=(tpindex-1)/stride + 1
        tpindex = tpindex - (stack(iss)-1)*stride
     end do

     ifvalidgrid=.true.
     do iss=1,nss
        if ( (stack(iss) < 1) .or. (stack(iss)>ssngrid(iss)) ) then
           ifvalidgrid=.false.; exit
        end if
     end do
     if ( ifvalidgrid .eqv. .false. ) cycle
     tpvec(1:nss) = sslowrange(1:nss) + (stack(1:nss)-1)*ssresolution(1:nss)

     tpmarkdis=0.0d0; tpvec(1:nss) = tpvec(1:nss) - sscenter(1:nss)
     j=1; tpdis=ssrad(1)
     do i=2,nss
        if ( ssrad(i) > tpdis ) then
           j = i; tpdis=ssrad(i)
        end if
     end do
     if ( ssshape(1:3) == "cir" ) then
        tpratiovec(1:nss) = ssrad(1:nss)/tpdis
        tpvec(1:nss) = (tpvec(1:nss)/tpratiovec(1:nss))**2
        tpdis = sqrt(sum(tpvec))
        if ( tpdis > ssrad(j) ) tpmarkdis = tpdis - ssrad(j)
     else if ( ssshape(1:3) == "rec" ) then
        do iss=1,nss
           if ( abs(tpvec(iss)) < ssrad(iss) ) cycle
           tpmarkdis = tpmarkdis + (abs(tpvec(iss)) - ssrad(iss))**2
        end do
        tpmarkdis = sqrt(tpmarkdis)
     else
!        write(iosiminfo,"(a)") "no ssshape for iniaddgriddata"
     end if
     iniaddgriddata(igrid) =  ssslope*tpmarkdis
  end do
end subroutine

subroutine flood_updatebiaspotentialgriddata(tpsspos,tpheight)
  use nfe_colvar_mod, only : colvar_value, colvar_force,colvar_is_periodic 
  real*8,intent(in) :: tpsspos(nss)
  real*8,intent(in),optional :: tpheight
  integer :: owngrid(nss),tpindex,iss,stride,igrid
  integer :: stack(nss),tpstack(nss),tpssntotgrid
  real*8 :: xi,biweightkernel,tpconst,localbias(nss),tpfactor
  logical :: ifvalidgrid

  if ( floodingtime == 0.0d0 ) call errormsg("floodingtime cannot be zero!")
  if ( wellkelvin < 0d0 ) call errormsg("when the bias potential is updating, wellkelvin cannot be a negative number!")
 
  tpconst=deltatime*gasconst*kelvin0/abs(floodingtime)
  if ( wellkelvin > 0.0d0 ) then
     if ( .not. present(tpheight) ) call errormsg("parameter tpheight must be presented for well-tempered ABMD!")
     tpconst=exp(-tpheight/(gasconst*wellkelvin))*tpconst
  end if

  call subspace_getowngridandbias(tpsspos,owngrid,localbias)
  tpssntotgrid=4**nss
  do igrid=1,tpssntotgrid
     stride=tpssntotgrid; tpindex=igrid
     do iss=nss,1,-1
        stride=stride/4
        stack(iss)=(tpindex-1)/stride + 1
        tpindex = tpindex - (stack(iss)-1)*stride
     end do

     tpstack(1:nss)=owngrid(1:nss)+stack(1:nss)-2
     ifvalidgrid=.true.
     do iss=1,nss
        if ( colvar_is_periodic(nfecvs(ss2cv(iss))) .eqv. .true. ) then
           if ( tpstack(iss) < 1 ) then
              tpstack(iss) = tpstack(iss) + ssngrid(iss)
           else if ( tpstack(iss) > ssngrid(iss) ) then
              tpstack(iss) = tpstack(iss) - ssngrid(iss)
           end if
        else
           if ( (tpstack(iss) < 1) .or. (tpstack(iss)>ssngrid(iss)) ) then
              ifvalidgrid=.false.; exit
           end if
        end if
     end do
     if ( ifvalidgrid .eqv. .true. ) then
        call subspace_getgridindex(tpstack,tpindex); tpfactor=tpconst
        do iss=1,nss
           select case(stack(iss))
           case (1)
              xi=1.0d0+localbias(iss)
           case (2)
              xi=localbias(iss)
           case (3)
              xi=1.0d0 - localbias(iss)
           case (4)
              xi=2.0d0 - localbias(iss)
           case default
              call errormsg("xi error!")
           end select
           biweightkernel = ((1-((xi**2)/4.0d0))**2)*48.0d0/41.0d0
           tpfactor = tpfactor*biweightkernel
        end do
        addbiasgriddata(tpindex) = addbiasgriddata(tpindex) + tpfactor
     end if
  end do
end subroutine

subroutine flood_generatefreesurface()
  integer :: igrid,stride,stack(nss),idim,tpindex
  real*8 :: surfaceheight,temp,tpvec(nss),tpgriddata(ssntotgrid),tpheight

  tpgriddata = addbiasgriddata - iniaddgriddata
  do igrid=1,ssntotgrid
     stride=ssntotgrid; tpindex=igrid
     do idim=nss,1,-1
        stride=stride/ssngrid(idim)
        stack(idim)=(tpindex-1)/stride + 1
        tpindex = tpindex - (stack(idim)-1)*stride
     end do
     tpvec(1:nss) = sslowrange(1:nss) + dble(stack(1:nss)-1)*ssresolution(1:nss)
     call subspace_getsurfaceinfo(tpvec,tpgriddata*biasweight,tpheight); freegriddata(igrid)=-tpheight
  end do

  if ( floodingtime > 0.0d0 ) then

     if ( (biasweight == 0.0d0) .or. (iffloodeds .eqv. .true.) ) then
        call subspace_hisgriddata2freegriddata(hisgriddata,freegriddata)
     else if ( biasweight == 1d0 ) then
        if ( wellkelvin > 0.0d0 ) freegriddata = (1.0d0+(kelvin0/wellkelvin))*freegriddata
     else 
        call errormsg("error biasweight ",real1=biasweight) 
     end if

  else if ( floodingtime == 0.0d0 ) then

     call subspace_hisgriddata2freegriddata(hisgriddata,tpgriddata)
     freegriddata = tpgriddata + freegriddata

  end if
end subroutine

end module 
