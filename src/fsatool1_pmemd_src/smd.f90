! Written by Chen Changjun, Huazhong University of Science and Technology, 2018

! Reference 1
! C. Jarzynski
! Nonequilibrium equality for free energy differences.
! Phys. Rev. Lett., 1997, 78, 2690-2693

! Reference 2
! Park S , Khaliliaraghi F , Tajkhorshid E , et al. 
! Free energy calculation from steered molecular dynamics simulations using Jarzynski's equality.
! Journal of Chemical Physics, 2003, 119(6):3559-3566.

! Reference 3
! Chen L Y , Bastien D A , Espejel H E.
! Determination of equilibrium free energy from nonequilibrium work measurements.
! Physical Chemistry Chemical Physics, 2010, 12(25):6579.

module smd
  use moldata; use colvar
  implicit none
  public :: smd_init, smd_updateforce

  private
  integer :: smdpathfrag,smdnpath,smdleadcv,smdrelaxstep,reversefrag,speed2frag,smdweightfreq
  real*8,allocatable,dimension(:) :: smdcvspeed,smdinicv,smdendcv,smdcv,refcvspeed,refcvspeed2
  real*8,allocatable,dimension(:) :: smdk,smdkv,smdkini,smdkend,refsmdkv
  real*8,allocatable,dimension(:,:) :: fragcv
contains

subroutine smd_init()
  integer :: ierr,icv,ifrag,smdtotstep,tpstep
  real*8 :: cvvec(ncv),smdspeed,smdspeed2
  namelist / smd / smdleadcv, smdrelaxstep, reversefrag, smdspeed, smdinicv, smdendcv, smdk, smdkend, &
              smdpathfrag, smdnpath, smdspeed2, speed2frag, smdweightfreq

  allocate(smdcvspeed(ncv),smdinicv(ncv),smdendcv(ncv),smdcv(ncv))
  smdspeed=0d0; smdspeed2=0d0; smdinicv=0d0; smdendcv=0d0; smdcv=0d0; reversefrag=0
  allocate(refcvspeed(ncv),refcvspeed2(ncv)); refcvspeed=0d0; refcvspeed2=0d0
  allocate(smdk(ncv),smdkini(ncv),smdkend(ncv),smdkv(ncv),refsmdkv(ncv))
  smdkini=0d0; smdk=0d0; smdkend=-1d0; smdkv=0d0; refsmdkv=0d0
 
  read(iomdin,nml=smd,iostat=ierr); if ( ierr < 0 ) call errormsg("namelist smd error! ")
  if ( (speed2frag<0) .or. (speed2frag>smdpathfrag) ) then
     call errormsg("speed2frag exceeds the range ",int1=speed2frag,int2=smdpathfrag)
  else if ( speed2frag == 0 ) then
     speed2frag = smdpathfrag
     if ( smdspeed2 /= 0d0 ) call errormsg("smdspeed2 should be zero when speed2frag is zero")
  else
     if ( smdspeed2 == 0d0 ) call errormsg("smdspeed2 cannot be zero with non-negative speed2frag")
  end if
  if ( (smdendcv(smdleadcv)-smdinicv(smdleadcv))*smdspeed < 0d0 ) call errormsg("smdspeed error!",real1=smdspeed)

  allocate(fragcv(ncv,0:smdpathfrag)); fragcv(1:ncv,0) = smdinicv(1:ncv)
  cvvec(1:ncv) = (smdendcv(1:ncv) - smdinicv(1:ncv))/dble(smdpathfrag)
  do ifrag=1,smdpathfrag
     fragcv(1:ncv,ifrag) = smdinicv(1:ncv) + dble(ifrag)*cvvec(1:ncv)
  end do

  smdtotstep = nint((fragcv(smdleadcv,speed2frag) - fragcv(smdleadcv,0))/(smdspeed*1.0d-3*deltatime))
  refcvspeed(1:ncv) = (fragcv(1:ncv,speed2frag) - fragcv(1:ncv,0))/dble(smdtotstep*deltatime*1d-3)
  refcvspeed(smdleadcv) = smdspeed; tpstep=0; smdtotstep = abs(smdtotstep)
  
  if ( smdspeed2 /= 0d0 ) then
     tpstep = nint((fragcv(smdleadcv,smdpathfrag) - fragcv(smdleadcv,speed2frag))/(smdspeed2*1.0d-3*deltatime))
     refcvspeed2(1:ncv) = (fragcv(1:ncv,smdpathfrag) - fragcv(1:ncv,speed2frag))/dble(tpstep*deltatime*1d-3)
     refcvspeed2(smdleadcv) = smdspeed2; tpstep = abs(tpstep)
  end if
  smdcv = smdinicv; smdcvspeed = refcvspeed
  smdkini = smdk; if ( smdkend(1) < 0d0 ) smdkend = smdk
  refsmdkv = (smdkend - smdkini)/(dble(smdtotstep+tpstep)*deltatime*1d-3); smdkv = refsmdkv

  if ( procid == 0 ) then
     write(iosiminfo,*); write(iosiminfo,"(a)") "Perform SMD simulation"
     write(iosiminfo,"(a,i5,a,i8,a,i5,a,i5,a,i5)") "path number ",smdnpath, &
             ", relax step ",smdrelaxstep, &
             ", reverse frag ",reversefrag,", speed2frag ",speed2frag,", smdweightfreq ",smdweightfreq

     call colvar_calandrec(cvvec,ifnorecord=.true.)
     write(iosiminfo,"(a,i4,5x,a,100f10.3)") "leading cv index ",smdleadcv,"initial cv ",cvvec(:)
     write(iosiminfo,"(a)") "first stage"
     do icv=1,ncv
        write(iosiminfo,"(a,i3,a,2f8.3,2(a,f8.3))") "cv ",icv,": ini-end ",smdinicv(icv), &
            fragcv(icv,speed2frag), ", speed ",refcvspeed(icv),", k ",smdk(icv)
     end do
     if ( reversefrag > 0 ) smdtotstep = smdtotstep*2
     write(iosiminfo,"(a,f8.3,a,i3)") "simulation time ",smdtotstep*deltatime/1d3," ns"

     if ( smdspeed2 /= 0d0 ) then
        write(iosiminfo,"(a)") "second stage"
        do icv=1,ncv
           write(iosiminfo,"(a,i3,a,2f8.3,2(a,f8.3))") "cv ",icv,": ini-end ",fragcv(icv,speed2frag), &
               smdendcv(icv), ", speed ",refcvspeed2(icv),", k ",smdk(icv)
        end do
        if ( reversefrag > 0 ) tpstep = tpstep*2
        write(iosiminfo,"(a,f8.3,a)") "simulation time ",tpstep*deltatime/1d3," ns"
     end if
     write(iosiminfo,"(a,f8.3,a)") "total simulatin time ",(smdtotstep+tpstep+smdrelaxstep)*deltatime/1d3," ns"
     
     write(iosiminfo,*)
     write(iosiminfo,"(a,100f8.3)") "initial spring constant ",smdkini
     write(iosiminfo,"(a,100f8.3)") "  final spring constant ",smdkend
     write(iosiminfo,*)
     flush(iosiminfo)
  end if
end subroutine

subroutine smd_updateforce()
  use mdin_ctrl_dat_mod, only : using_pme_potential
  use nfe_colvar_mod, only : colvar_force
  integer :: icv,ierr,iproc,i,k,imol,ires,tpstep,iatom
  integer,dimension(0:procnum-1) :: tpcountarray,tpposarray
  real*8 :: cvvec(ncv),tpcoeff,tppotvec(ncv),temp,tptime,smdpot,tptpfrc(3,numboxatom)
  real*8,dimension(0:smdpathfrag) :: tppathvec,pathfreef,pathfreeb,pathfreefb
  integer,save :: ioweight,fragdelta,iopathfree,ifrag,iopaths,npath,smdstartstep,prefragmark,fragmark
  real*8,save :: smdwork,smdtotwork
  real*8,dimension(:),allocatable,save :: pathworkf,pathworkb,precvforce,pathexpsumf2,pathexpsumb2
  real*8,dimension(:),allocatable,save :: pathworksumf,pathworksumb,pathexpsumf,pathexpsumb,prekforce
  real*8,dimension(:,:),allocatable,save :: tpcoor,tpvel,tpfrc,tpinicoor,tpinifrc
  character*100 :: tpfilename
  character*100,save :: trajfilename
  logical,save :: ifpathend,ifrelaxed

  if ( .not. allocated (pathworkf) ) then

     smdwork = 0d0; npath=0; ifrag=1; ifpathend=.false.; ifrelaxed=.false.; smdstartstep=nowstep; fragdelta=0
     smdtotwork = 0d0

     if ( smdweightfreq > 0 ) call getfreeunit(ioweight)
     write(tpfilename,"(i4)") npath+procid+1; tpfilename=adjustl(tpfilename)
     tpfilename="procinfo/weight_"//trim(tpfilename)//".txt"
     open(file=tpfilename, unit=ioweight, action="write")
     write(tpfilename,"(i4)") npath+procid+1; tpfilename=adjustl(tpfilename)
     trajfilename="procinfo/traj_"//trim(tpfilename)//".mdcrd"
     call moldata_write_traj(tpfilename=trajfilename)

     allocate(pathworkf(0:smdpathfrag),pathworkb(0:smdpathfrag),precvforce(ncv),prekforce(ncv))
     pathworkf=0d0; pathworkb=0d0; precvforce=0d0; fragmark=0; prefragmark=0; prekforce=0d0
     allocate(tpcoor(3,numboxatom),tpvel(3,numboxatom),tpfrc(3,numboxatom))
     tpcoor=0d0; tpvel=0d0; tpfrc=0d0
     allocate(tpinicoor(3,numboxatom),tpinifrc(3,numboxatom))
     tpinicoor=0d0; tpinifrc=0d0

     if ( procid == 0 ) then
        allocate(pathworksumf(0:smdpathfrag),pathworksumb(0:smdpathfrag))
        pathworksumf=0d0; pathworksumb=0d0
        allocate(pathexpsumf(0:smdpathfrag),pathexpsumb(0:smdpathfrag))
        pathexpsumf=0d0; pathexpsumb=0d0
        allocate(pathexpsumf2(0:smdpathfrag),pathexpsumb2(0:smdpathfrag))
        pathexpsumf2=0d0; pathexpsumb2=0d0
     end if
  end if

  call moldata_download_crd(nowcoor)
  call colvar_calandrec(cvvec,ifonlycv=.true.)

  nowforce=0.0d0
  do icv = 1, ncv
     tpcoeff = -smdk(icv)*(cvvec(icv)-smdcv(icv))
     if ( icv <= nfencv ) then
        call colvar_force(nfecvs(icv), nowcoor, tpcoeff, nowforce)
     else
        call colvar_fsainfo(icv-nfencv, nowcoor, temp, tpgrad=nowforce, tpfactor=tpcoeff);
     end if
     if ( smdcvspeed(icv) /= 0d0 ) then
        if (precvforce(icv) /= 0d0) then
           temp = 0.5d0*deltatime*1d-3*(tpcoeff + precvforce(icv))*smdcvspeed(icv)
           smdwork = smdwork + temp
           smdtotwork = smdtotwork + temp
        end if
     end if
     precvforce(icv) = tpcoeff
     if ( smdkv(icv) /= 0d0 ) then
        tpcoeff = 0.5d0*((cvvec(icv)-smdcv(icv))**2)
        if ( prekforce(icv) /= 0d0 ) then
           temp = 0.5d0*deltatime*1d-3*(tpcoeff + prekforce(icv))*smdkv(icv)
           smdwork = smdwork + temp
           smdtotwork = smdtotwork + temp
        end if
        prekforce(icv) = tpcoeff
     end if
  end do

  call moldata_upload_frc_add(nowforce)

  tpstep = nowstep - smdstartstep
  if ( (ifrelaxed .eqv. .false.) .and. (tpstep <= smdrelaxstep) ) then
     if ( tpstep == smdrelaxstep ) then
        ifrelaxed=.true.; smdstartstep=nowstep
        call gpu_download_crd(tpinicoor); call gpu_download_frc(tpinifrc)
        precvforce=0d0; prekforce=0d0; smdwork=0d0; smdtotwork=0d0
     end if
     return
  else if ( mod(tpstep,samplestep) == 0 ) then
     call mpi_barrier(mpi_comm_world,ierr)
  end if

  if ( ( smdweightfreq > 0 ) .and. (mod(tpstep,smdweightfreq) == 0) ) then
     tppotvec(1:ncv) = 0.5d0*smdk(1:ncv)*((cvvec(1:ncv)-smdcv(1:ncv))**2); smdpot = sum(tppotvec)
     temp = sqrt(dot_product(precvforce,precvforce))
     write(ioweight, "(f6.2,100f8.3)") dble(tpstep)*deltatime*1d-3, smdcv(smdleadcv), cvvec(smdleadcv), &
           precvforce(smdleadcv),temp,smdpot,smdtotwork
     flush(ioweight)
  end if
  if ( mod(tpstep,savetrajstep) == 0 ) call moldata_write_traj(tpfilename=trajfilename)

  smdcv(1:ncv) = smdcv(1:ncv) + smdcvspeed(1:ncv)*1d-3*deltatime
  smdk(1:ncv) = smdk(1:ncv) + smdkv(1:ncv)*1d-3*deltatime
 
  if ( (smdcvspeed(smdleadcv)*refcvspeed(smdleadcv)>0d0) .and. &
       ( ( (refcvspeed(smdleadcv)>0d0) .and. (smdcv(smdleadcv)>=fragcv(smdleadcv,ifrag)) ) .or.  &
         ( (refcvspeed(smdleadcv)<0d0) .and. (smdcv(smdleadcv)<=fragcv(smdleadcv,ifrag)) ) )   ) then

     if ( pathworkf(ifrag) /= 0d0 ) call errormsg("pathworkf(ifrag) should be zero! ",int1=ifrag)
     fragdelta = fragdelta + 1; pathworkf(ifrag) = smdwork; smdwork=0d0
     if ( (reversefrag == 0) .and. (ifrag == smdpathfrag) ) then
        ifpathend=.true.
     else if ( (fragdelta == reversefrag) .or. (ifrag == smdpathfrag) ) then
        smdcvspeed(1:ncv) = -refcvspeed(1:ncv)
        if ( ifrag > speed2frag ) smdcvspeed(1:ncv) = -refcvspeed2(1:ncv)
        if (fragmark/=0) prefragmark=fragmark; fragmark=ifrag
        fragdelta=0; precvforce=0d0; prekforce=0d0
        call gpu_download_crd(tpcoor); call gpu_download_vel(tpvel); call gpu_download_frc(tpfrc)
        tpcoor = nowcoor; tpvel = nowvel; tpfrc = nowforce
        smdkv(1:ncv) = -refsmdkv(1:ncv)
     else
        ifrag = ifrag + 1
        if ( ifrag > speed2frag ) smdcvspeed(1:ncv) = refcvspeed2(1:ncv)
     end if

  else if ( (smdcvspeed(smdleadcv)*refcvspeed(smdleadcv)<0d0) .and. &
            ( ( (refcvspeed(smdleadcv)>0d0) .and. (smdcv(smdleadcv)<=fragcv(smdleadcv,ifrag-1)) ) .or.  &
              ( (refcvspeed(smdleadcv)<0d0) .and. (smdcv(smdleadcv)>=fragcv(smdleadcv,ifrag-1)) ) )   ) then

     if ( pathworkb(ifrag) /= 0d0 ) call errormsg("pathworkb(ifrag) should be zero! ",int1=ifrag)
     pathworkb(ifrag) = smdwork; smdwork=0d0
     if ( ifrag == (prefragmark+1) ) then
        if ( fragmark == smdpathfrag ) then
           ifpathend=.true.
        else
           precvforce=0d0; prekforce=0d0; smdcv(1:ncv)=fragcv(1:ncv,fragmark); ifrag=fragmark+1
           smdcvspeed(1:ncv) = refcvspeed(1:ncv)
           if ( ifrag > speed2frag ) smdcvspeed(1:ncv) = refcvspeed2(1:ncv)
           call gpu_upload_crd(tpcoor); call gpu_upload_vel(tpvel); call gpu_upload_frc(tpfrc)
           if (using_pme_potential) call gpu_force_new_neighborlist()
           smdwork = 0d0
           do i=1,fragmark
              smdwork = smdwork + 0.5d0*( pathworkf(i) - pathworkb(i))
           end do
           smdtotwork = smdwork; smdwork=0d0
           smdkv(1:ncv) = refsmdkv(1:ncv)
        end if
     else
        ifrag = ifrag - 1
        if ( ifrag <= speed2frag ) smdcvspeed(1:ncv) = -refcvspeed(1:ncv)
     end if

  end if
 
  if ( ifpathend .eqv. .true. ) then

     do i=1,smdpathfrag
        if ( pathworkf(i) == 0d0 ) call errormsg("pathworkf has zero elements! ",int1=procid, int2=i)
        if ( (reversefrag>0) .and. (pathworkb(i) == 0d0) ) call errormsg("pathworkb has zero elements! ",int1=procid, int2=i)
     end do

     if ( iopaths == 0 ) then
        call getfreeunit(iopaths); write(tpfilename,"(i4)") procid; tpfilename=adjustl(tpfilename)
        tpfilename="procinfo/smd_paths_"//trim(tpfilename)//".txt"
        open(file=tpfilename, unit=iopaths, action="write")
     end if

     tppathvec = pathworkf
     do i=2,smdpathfrag
        tppathvec(i) = tppathvec(i) + tppathvec(i-1)
     end do
     write(iopaths,"(a14,i5,a)") "forward path  ",npath+procid+1,", work"
     write(iopaths,"(10f10.3)") tppathvec
     pathworkf = tppathvec

     tppathvec = pathworkb
     do i=2,smdpathfrag
        tppathvec(i) = -tppathvec(i) + tppathvec(i-1)
     end do
     write(iopaths,"(a14,i5,a)") "backward path  ",npath+procid+1,", work"
     write(iopaths,"(10f10.3)") tppathvec
     flush(iopaths)
     pathworkb = tppathvec

     call mpi_reduce(pathworkf, tppathvec, smdpathfrag+1,  &
              mpi_double_precision, MPI_SUM, 0, mpi_comm_world, ierr)
     if ( procid == 0 ) then
        pathworksumf = pathworksumf + tppathvec
     end if

     pathfreef = exp(-pathworkf/(gasconst*kelvin0))
     call mpi_reduce( pathfreef, tppathvec, smdpathfrag+1,  &
              mpi_double_precision, MPI_SUM, 0, mpi_comm_world, ierr)
     if ( procid == 0 ) then
        pathexpsumf = pathexpsumf + tppathvec
     end if

     pathfreef = exp(-0.5d0*pathworkf/(gasconst*kelvin0))
     call mpi_reduce( pathfreef, tppathvec, smdpathfrag+1,  &
              mpi_double_precision, MPI_SUM, 0, mpi_comm_world, ierr)
     if ( procid == 0 ) then
        pathexpsumf2 = pathexpsumf2 + tppathvec
     end if

     if ( reversefrag > 0 ) then
     call mpi_reduce(pathworkb, tppathvec, smdpathfrag+1,  &
              mpi_double_precision, MPI_SUM, 0, mpi_comm_world, ierr)
     if ( procid == 0 ) then
        pathworksumb = pathworksumb + tppathvec
     end if

     pathfreeb = exp(-pathworkb/(gasconst*kelvin0))
     call mpi_reduce( pathfreeb, tppathvec, smdpathfrag+1,  &
              mpi_double_precision, MPI_SUM, 0, mpi_comm_world, ierr)
     if ( procid == 0 ) then
        pathexpsumb = pathexpsumb + tppathvec
     end if

     pathfreeb = exp(0.5d0*pathworkb/(gasconst*kelvin0))
     call mpi_reduce( pathfreeb, tppathvec, smdpathfrag+1,  &
              mpi_double_precision, MPI_SUM, 0, mpi_comm_world, ierr)
     if ( procid == 0 ) then
        pathexpsumb2 = pathexpsumb2 + tppathvec
     end if
     end if

     npath = npath + procnum

     if ( procid == 0 ) then

        pathfreef=0d0; pathfreeb=0d0; pathfreefb=0d0

        pathworkf = pathworksumf/dble(npath)
        pathfreef = -gasconst*kelvin0*log(pathexpsumf/dble(npath))
        pathfreef = pathfreef - minval(pathfreef)

        if ( reversefrag > 0 ) then
           pathworkb = pathworksumb/dble(npath)
           pathfreeb = -gasconst*kelvin0*log(pathexpsumb/dble(npath))
           pathfreeb = pathfreeb - minval(pathfreeb)
           pathfreefb = -gasconst*kelvin0*log( (pathexpsumf2/dble(npath))/(pathexpsumb2/dble(npath)))
           pathfreefb = pathfreefb - minval(pathfreefb)
        end if

        call getfreeunit(iopathfree); open(unit=iopathfree,file="procinfo/smd_pathfree.txt",action="write")
        if ( fragcv(smdleadcv,smdpathfrag) > fragcv(smdleadcv,0) ) then
           do ifrag=0,smdpathfrag
              write(iopathfree,"(100f10.3)") dble(ifrag)/dble(smdpathfrag), fragcv(smdleadcv,ifrag), &
                   pathworkf(ifrag), pathfreef(ifrag), pathworkb(ifrag), pathfreeb(ifrag), pathfreefb(ifrag)
           end do
        else if ( fragcv(smdleadcv,smdpathfrag) < fragcv(smdleadcv,0) ) then
           do ifrag=smdpathfrag,0,-1
              write(iopathfree,"(100f10.3)") dble(ifrag)/dble(smdpathfrag), fragcv(smdleadcv,ifrag), &
                   pathworkf(ifrag), pathfreef(ifrag), pathworkb(ifrag), pathfreeb(ifrag), pathfreefb(ifrag)
           end do
        end if
        close(iopathfree)

     end if

     call mpi_barrier(mpi_comm_world,ierr)
     if ( npath >= smdnpath ) call errormsg("smd simulation has been finished!")

     smdwork = 0d0; precvforce=0d0; prekforce=0d0;
     pathworkf=0d0; pathworkb=0d0; fragdelta=0; fragmark=0; prefragmark=0
     ifpathend=.false.; ifrelaxed=.false.; smdstartstep=nowstep; smdtotwork=0d0; smdk=smdkini
     smdcv = fragcv(1:ncv,0);  smdcvspeed(1:ncv) = refcvspeed(1:ncv); ifrag = 1
     if (using_pme_potential) call gpu_force_new_neighborlist()
     call gpu_upload_crd(tpinicoor); call gpu_upload_frc(tpinifrc); nowcoor=tpinicoor; nowforce=tpinifrc
     smdk(1:ncv) = smdkini(1:ncv); smdkv(1:ncv) = refsmdkv(1:ncv)

     if ( smdweightfreq > 0 ) then
        close(ioweight)
        write(tpfilename,"(i4)") npath+procid+1; tpfilename=adjustl(tpfilename)
        tpfilename="procinfo/weight_"//trim(tpfilename)//".txt"
        open(file=tpfilename, unit=ioweight, action="write")
     end if
     write(tpfilename,"(i4)") npath+procid+1; tpfilename=adjustl(tpfilename)
     trajfilename="procinfo/traj_"//trim(tpfilename)//".mdcrd"
     call moldata_write_traj(tpfilename=trajfilename) 
  end if
end subroutine

end module 
