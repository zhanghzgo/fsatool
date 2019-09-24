! Reference
! Protein Folding Pathways Revealed by Essential Dynamics Sampling
! D. Narzi, I. Daidone, A. Amadei and A.D. Nola
! J. Chem. Theory Comput. 2008, 4, 1940-1948

module eds
  use moldata; use colvar
  implicit none
  public :: eds_init, eds_updatecoor, eds_flood_readbox

  private
  integer :: edsnpc,edsnatom,step0low,step0up,edsfixstep,edsmdstep,edsfixstep2,edscv,edscv2,edsnbox
  real*8 :: edsdislow,edsdisup,shell,expandkelvin0,cpdislen
  integer,allocatable,dimension(:) :: edsatoms
  real*8,allocatable,dimension(:) :: xi0,edsrefcoor
  real*8,allocatable,dimension(:,:) :: edsvec

  integer :: edsnpc2,edsnatom2,edsnsnap
  real*8 :: edsdislow2,edsdisup2,cvxi0,cvxi02
  integer,allocatable,dimension(:) :: edsatoms2
  real*8,allocatable,dimension(:) :: xi02,edsrefcoor2
  real*8,allocatable,dimension(:,:) :: edsvec2

  integer :: edsreadbox,edswritebox
  real*8 :: edsboxreadtime
  real*8,allocatable,dimension(:) :: edsboxtime
  real*8,allocatable,dimension(:,:,:) :: edscrdbox,edsvelbox,edsfrcbox
contains

subroutine eds_init()
  namelist / eds / edsatoms, edsdislow, edsdisup, pcstart, pcend, eigenfile, edsfixstep, edsmdstep, &
                   edsatoms2,  pcstart2, pcend2, eigenfile2, edsfixstep2, edsdislow2,edsdisup2, shell, &
                   edscv, edscv2, cvxi0, cvxi02, step0low, step0up, expandkelvin0, cpdislen,edsnbox
  integer :: tpndim,idim,ierr,iatom,tpnatom,tpatoms(numatom),pcstart,pcend,iofile, pcstart2, pcend2, tpndim2
  real*8 :: tpvec(numatom*3),tpdis,tpdis2,temp,cvvec(ncv)
  character*100 :: eigenfile,eigenfile2

  allocate(edsatoms(numatom),edsatoms2(numatom)); edsnatom=0; edsnatom2=0; edsatoms=0; edsatoms2=0; shell=0d0
  edsmdstep=0; edsfixstep=1; edsfixstep2=1; edsnpc=0; edsnpc2=0; pcstart=0; pcend=0; pcstart2=0; pcend2=0
  edscv=0; edscv2=0; edsdislow=0d0; edsdisup=0d0; step0low=0; step0up=0; expandkelvin0=0d0
  cpdislen=10d0; edsnbox=10
  read(iomdin,nml=eds,iostat=ierr); if ( ierr < 0 ) call errormsg("namelist eds error! ")
  if ( step0low == 0 ) then
     call errormsg("step0low should not be zero!")
  else if ( step0up == 0 ) then
     step0up = step0low
  end if
  if ( edsatoms(1) /= 0 ) then
     do iatom=1,numatom
        if ( edsatoms(iatom) == 0 ) exit
     end do
     edsnatom = iatom - 1
     tpatoms(1:edsnatom) = edsatoms(1:edsnatom)
     deallocate(edsatoms); allocate(edsatoms(1:edsnatom)); edsatoms(1:edsnatom) = tpatoms(1:edsnatom)
     allocate(edsrefcoor(3*edsnatom)); edsrefcoor=0d0
     tpndim=edsnatom*3
  end if

  if ( edsnatom > 0 ) then
  edsnpc = pcend - pcstart + 1
  allocate(edsvec(tpndim,edsnpc),xi0(edsnpc)); edsvec=0d0; xi0=0d0 
  call getfreeunit(iofile); open(unit=iofile,file=eigenfile,action="read")
  read(iofile,"(10x,10000f8.3)") tpvec(:)
  if ( edsnatom > 0 ) then
     edsrefcoor(1:3*edsnatom)=tpvec(1:3*edsnatom)
     tpvec(1:3)=0.0d0; do idim=1,edsnatom; tpvec(1:3)=tpvec(1:3)+edsrefcoor(3*idim-2:3*idim); end do
     tpvec(1:3)=tpvec(1:3)/dble(edsnatom)
     do idim=1,edsnatom; edsrefcoor(3*idim-2:3*idim)=edsrefcoor(3*idim-2:3*idim)-tpvec(1:3); end do
  end if
  read(iofile,"(10x,10000f8.3)") tpvec(:)
  do idim=pcstart, pcend 
     xi0(idim-pcstart+1) = tpvec(idim)
  end do
  read(iofile,*) 
  do idim=1,pcend
     read(iofile,"(4x,10000f8.3)") tpvec
     temp = sqrt(dot_product(tpvec(1:tpndim),tpvec(1:tpndim)))
     if ( abs(1.0d0-temp) > 0.01d0 ) call errormsg("eds axis is not normalized ",int1=idim,real1=temp)
     if ( idim >= pcstart) edsvec(1:tpndim,idim-pcstart+1)=tpvec(1:tpndim)
  end do
  close(iofile)
  end if

  if ( pcstart2 > 0 ) then
     do iatom=1,numatom
        if ( edsatoms2(iatom) == 0 ) exit
     end do
     edsnatom2 = iatom - 1
     if ( edsnatom2 > 0 ) then
     edsnpc2 = pcend2 - pcstart2 + 1; if ( edsnpc2 <= 0 ) call errormsg("wrong edsnpc2!")
     tpatoms=0; tpatoms(1:edsnatom2) = edsatoms2(1:edsnatom2)
     deallocate(edsatoms2); allocate(edsatoms2(1:edsnatom2)); edsatoms2(1:edsnatom2) = tpatoms(1:edsnatom2)
     allocate(edsrefcoor2(3*edsnatom2)); edsrefcoor2=0d0
     tpndim2=edsnatom2*3

     allocate(edsvec2(tpndim2,edsnpc2),xi02(edsnpc2)); edsvec2=0d0; xi02=0d0

     open(unit=iofile,file=eigenfile2,action="read")
     read(iofile,"(10x,10000f8.3)") tpvec(:)
     edsrefcoor2(1:3*edsnatom2)=tpvec(1:3*edsnatom2)
     tpvec(1:3)=0.0d0; do idim=1,edsnatom2; tpvec(1:3)=tpvec(1:3)+edsrefcoor2(3*idim-2:3*idim); end do
     tpvec(1:3)=tpvec(1:3)/dble(edsnatom2)
     do idim=1,edsnatom2; edsrefcoor2(3*idim-2:3*idim)=edsrefcoor2(3*idim-2:3*idim)-tpvec(1:3); end do

     read(iofile,"(10x,10000f8.3)") tpvec(:)
     do idim=pcstart2, pcend2
        xi02(idim-pcstart2+1) = tpvec(idim)
     end do
     read(iofile,*)
     do idim=1,pcend2
        read(iofile,"(4x,10000f8.3)") tpvec
        temp = sqrt(dot_product(tpvec(1:tpndim2),tpvec(1:tpndim2)))
        if ( abs(1.0d0-temp) > 0.01d0 ) call errormsg("eds2 axis is not normalized ",int1=idim,real1=temp)
        if ( idim >= pcstart2) edsvec2(1:tpndim2,idim-pcstart2+1)=tpvec(1:tpndim2)
     end do
     close(iofile)
     end if
  else
     deallocate(edsatoms2); edsnpc2=0; edsnatom2=0; edsfixstep2=0
  end if

  if ( edsnatom2 == 0 ) then
     call eds_pcspace_info(tpedsdis=tpdis)
  else
     call eds_pcspace_info(tpedsdis=tpdis,tpedsdis2=tpdis2)
  end if

  if ( procid == 0 ) then
     write(iosiminfo,*) 
     write(iosiminfo,"(a,f5.2,a,i10)") "Performing essential dyanmics sampling"
     write(iosiminfo,"(2(a,i8),4(a,f8.3))") "step0low ",step0low,", step0up ",step0up, &
                 ", shell ",shell,", expandkelvin0 ",expandkelvin0
     write(iosiminfo,"(a,i5,3x,a,f10.3)") "edsnbox ",edsnbox,", cpdislen ",cpdislen
     write(iosiminfo,"(3(a,i8))") "edsmdstep ",edsmdstep, &
                                  ", edsfixstep ",edsfixstep,", edsfixstep2 ",edsfixstep2
     write(iosiminfo,*)

     if ( pcstart > 0 ) then
        if ( edsnatom > 0 ) then
           write(iosiminfo,"(a,i4)") "Princial Component Space 1, edsnatom=",edsnatom
!           write(iosiminfo,"(a,1000f8.3)") "xi0:  ",xi0(1:edsnpc)
        else
           write(iosiminfo,"(a)") "Collective Variable Space 1"
           write(iosiminfo,"(a,i4,5x,a,f8.3)") "cv ",edscv," cvxi0 ",cvxi0
           if ( edscv <= 0 ) call errormsg("error edscv ",int1=edscv)
        end if
        write(iosiminfo,"(a,i10,3x,a,i10)") "pc start:  ",pcstart,", pc end: ",pcend
        write(iosiminfo,"(3(a,f8.3))") "dislow: ",edsdislow,", disup: ",edsdisup,", ini dis ",tpdis
     end if

     if ( pcstart2 > 0 ) then
        write(iosiminfo,*)
        if ( edsnatom2 > 0 ) then
           write(iosiminfo,"(a,i4)") "Princial Component Space 2, edsnatom2=",edsnatom2
!           write(iosiminfo,"(a,1000f8.3)") "xi02:  ",xi02(1:edsnpc2)
        else
           write(iosiminfo,"(a)") "Collective Variable Space 2"
           write(iosiminfo,"(a,i4,5x,a,f8.3)") "cv ",edscv2," cvxi0 ",cvxi02
           if ( edscv2 <= 0 ) call errormsg("error edscv ",int1=edscv2)
        end if
        write(iosiminfo,"(a,i10,3x,a,i10)") "pc start2  ",pcstart2,", pc end2: ",pcend2
        write(iosiminfo,"(3(a,f8.3))") "dislow2: ",edsdislow2,", disup2: ",edsdisup2,", ini dis ",tpdis2
     end if

     write(iosiminfo,*); flush(iosiminfo)
  end if

! should be modified here
! if ( nbox > 0 ) then
!     call gpu_setup_shuttle_info(cvnatom+edsnatom, 0, cvatoms+edsatoms)
! end if
end subroutine

subroutine eds_updatecoor()
  use mdin_ctrl_dat_mod, only : using_pme_potential
  integer,save :: ioedsinfo,iotraj,edsnpath, edspathstep,iedsfixstep, iedsmdstep, iedsfixstep2
  integer,save :: step0, lastupdatestep, ioedswrite
  integer :: idim,iatom,ierr,edssimstate
  real*8 :: tptime,tpdis,tpdis2,cvvec(ncv),tpcoor(3,numboxatom),tpbcastvec(8)
  real*8,save :: inwalldis,outwalldis,inwalldis2,outwalldis2,predis,predis2,inikelvin0,lastcpdis
  real*8,allocatable,dimension(:,:),save :: precoor,inicoor,inifrc
  character*100 :: tpfilename
  logical :: ifend,ifcpeds

  if ( ifeds .eqv. .true. ) call moldata_download_crd(nowcoor)

  if ( .not. allocated(precoor) ) then
     allocate(precoor(3,numboxatom)); precoor=nowcoor

     if ( procid == (procnum-1) ) then
        if ( ifeds .eqv. .false. ) then
           call getfreeunit(ioedsinfo); open(file="procinfo/eds_info.txt", unit=ioedsinfo, action="write")
           if ( edsnbox == 0 ) call errormsg("edsnbox is zero!")
           edsreadbox=1; edswritebox=1
           allocate(edsboxtime(edsnbox)); edsboxtime=0d0
           allocate(edscrdbox(3,numboxatom,edsnbox)); edscrdbox=0d0
           allocate(edsvelbox(3,numboxatom,edsnbox)); edsvelbox=0d0
           allocate(edsfrcbox(3,numboxatom,edsnbox)); edsfrcbox=0d0
        end if

        if ( edsdisup == 0d0 ) then
           allocate(inicoor(3,numboxatom),inifrc(3,numboxatom)); inicoor=0d0; inifrc=0d0
           call gpu_download_crd(inicoor); call gpu_download_frc(inifrc)
        else if ( expandkelvin0 > 0d0 ) then
           inikelvin0 = kelvin0
        end if
     end if

     if ( edsnpc2 == 0 ) then
        call eds_pcspace_info(tpedsdis=inwalldis)
        if ( ifeds .eqv. .true. ) call colvar_calandrec(cvvec,tpaddvar=(/inwalldis/),ifonlycv=.true.)
     else
        call eds_pcspace_info(tpedsdis=inwalldis,tpedsdis2=inwalldis2)
        if ( ifeds .eqv. .true. ) call colvar_calandrec(cvvec,tpaddvar=(/inwalldis,inwalldis2/),ifonlycv=.true.)
     end if
     outwalldis=inwalldis; outwalldis2=inwalldis2
     predis=inwalldis; predis2=inwalldis2; lastcpdis=predis

     call getfreeunit(iotraj); open(file="procinfo/eds_path_1.traj", unit=iotraj, action="write")
     write(iotraj,"(a12)") "default_name"; edsnpath=1; edspathstep=0
     tpcoor=nowcoor; call moldata_moveto_comorcenterbox(tpcoor)
     write(iotraj,"(10f8.3)") tpcoor; if (nbox>0) write(iotraj,"(3f8.3)") boxsize(1:3); flush(iotraj)

     iedsmdstep=0; iedsfixstep = edsfixstep; iedsfixstep2 = edsfixstep2; lastupdatestep=nowstep; step0=step0low
     if ( edsmdstep == 0 ) iedsfixstep = 0
     if ( edsfixstep == 0 ) iedsfixstep2 = 0

     if ( (edsmdstep+edsfixstep+edsfixstep2) == 0 ) &
         call errormsg("error edsmdstep, edsfixstep or edsfixstep2 ",int1=edsmdstep,int2=edsfixstep,int3=edsfixstep2)
  end if

  if ( edsnpc2 == 0 ) then
      call eds_pcspace_info(tpedsdis=tpdis)
  else
      call eds_pcspace_info(tpedsdis=tpdis, tpedsdis2=tpdis2)
  end if

  if ( ifeds .eqv. .false. ) then
     tpbcastvec(1:8)=(/tpdis,predis,inwalldis,outwalldis,tpdis2,predis2,inwalldis2,outwalldis2/)
     call mpi_bcast(tpbcastvec(1:8), 8, mpi_double_precision, procnum-1, mpi_comm_world, ierr)
     tpdis =tpbcastvec(1); predis =tpbcastvec(2); inwalldis =tpbcastvec(3); outwalldis =tpbcastvec(4)
     tpdis2=tpbcastvec(5); predis2=tpbcastvec(6); inwalldis2=tpbcastvec(7); outwalldis2=tpbcastvec(8)
     ifcpeds = .false.
  end if

  if ( (iedsmdstep < edsmdstep) .and. (edsmdstep > 0) ) then

     iedsmdstep = iedsmdstep + 1; edssimstate=0
     if ( iedsmdstep == edsmdstep ) then
        if ( edsfixstep > 0 ) then
           iedsfixstep = 0
        else if ( edsfixstep2 > 0 ) then
           iedsfixstep2 = 0
        else
           iedsmdstep = 0
        end if
     end if

  else if ( (iedsfixstep < edsfixstep) .and. (edsfixstep > 0) .and. (edsnpc > 0) ) then

     if ( step0 > 0 ) then
     if ( tpdis > outwalldis ) then
        if ( tpdis > predis ) then
           call eds_pcspace_fix(precoor,tprefdis=predis); tpdis=predis
           if ( procid == (procnum-1) ) call moldata_upload_crd(nowcoor)
        end if
     else if ( tpdis < inwalldis ) then
        inwalldis = tpdis; lastupdatestep=nowstep
        if ( outwalldis > (inwalldis+shell) ) outwalldis = inwalldis + shell
        if ( (ifeds .eqv. .false.) .and. (procid==(procnum-1)).and. (edsboxtime(edswritebox) == 0d0) ) then
           if ( (abs(lastcpdis-tpdis) > cpdislen) ) then
              ifcpeds=.true.; lastcpdis=tpdis
           end if
        end if
     end if
     else
     if ( tpdis < inwalldis ) then
        if ( (tpdis < predis) .and. (expandkelvin0 == 0d0) ) then
           call eds_pcspace_fix(precoor,tprefdis=predis); tpdis=predis
           if ( procid == (procnum-1) ) call moldata_upload_crd(nowcoor)
        end if
     else if ( tpdis > outwalldis ) then
        outwalldis = tpdis; lastupdatestep=nowstep
        if ( outwalldis > (inwalldis+shell) ) inwalldis = outwalldis - shell
        if ( (ifeds .eqv. .false.) .and. (procid==(procnum-1)).and. (edsboxtime(edswritebox) == 0d0) ) then
           if ( (abs(lastcpdis-tpdis) > cpdislen) ) then
              ifcpeds=.true.; lastcpdis=tpdis
           end if
        end if
     end if
     end if

     iedsfixstep = iedsfixstep + 1; edssimstate=1
     if ( iedsfixstep == edsfixstep ) then
        if ( edsfixstep2 > 0 ) then
           iedsfixstep2 = 0
        else if ( edsmdstep > 0 ) then
           iedsmdstep = 0
        else
           iedsfixstep = 0
        end if
     end if

  else if ( (iedsfixstep2 < edsfixstep2) .and. (edsfixstep2 > 0) .and. (edsnpc2 > 0) ) then

     if ( step0 > 0 ) then
     if ( tpdis2 > outwalldis2 ) then
        if ( tpdis2 > predis2 ) then
           call eds_pcspace_fix(precoor,tprefdis2=predis2); tpdis2=predis2
           if ( procid == (procnum-1) ) call moldata_upload_crd(nowcoor)
        end if
     else if ( tpdis2 < inwalldis2 ) then
        inwalldis2 = tpdis2; lastupdatestep=nowstep
        if ( outwalldis2 > (inwalldis2+shell) ) outwalldis2 = inwalldis2 + shell
        if ( (ifeds .eqv. .false.) .and. (procid==(procnum-1)).and. (edsboxtime(edswritebox) == 0d0) ) then
           if ( (abs(lastcpdis-tpdis2) > cpdislen) ) then
              ifcpeds=.true.; lastcpdis=tpdis2
           end if
        end if
     end if
     else
     if ( tpdis2 < inwalldis2 ) then
        if ( (tpdis2 < predis2) .and. (expandkelvin0 == 0d0) ) then
           call eds_pcspace_fix(precoor,tprefdis2=predis2); tpdis2=predis2
           if ( procid == (procnum-1) ) call moldata_upload_crd(nowcoor)
        end if
     else if ( tpdis2 > outwalldis2 ) then
        outwalldis2 = tpdis2; lastupdatestep=nowstep
        if ( outwalldis2 > (inwalldis2+shell) ) inwalldis2 = outwalldis2 - shell
        if ( (ifeds .eqv. .false.) .and. (procid==(procnum-1)).and. (edsboxtime(edswritebox) == 0d0) ) then
           if ( (abs(lastcpdis-tpdis2) > cpdislen) ) then
              ifcpeds=.true.; lastcpdis=tpdis2
           end if
        end if
     end if
     end if
     iedsfixstep2 = iedsfixstep2 + 1; edssimstate=2
     if ( iedsfixstep2 == edsfixstep2 ) then
        if ( edsmdstep > 0 ) then
           iedsmdstep = 0
        else if ( edsfixstep > 0 ) then
           iedsfixstep = 0
        else
           iedsfixstep2 = 0
        end if
     end if

  else

     call errormsg("error eds state ",int1=iedsmdstep,int2=iedsfixstep,int3=iedsfixstep2)

  end if
  precoor=nowcoor; predis=tpdis; predis2=tpdis2

  if ( mod(nowstep,samplestep) == 0 ) then

     tptime = dble(edspathstep)*deltatime*1d-3
     if ( ifeds .eqv. .true. ) then
        if ( edsnpc2 == 0 ) then
           call colvar_calandrec(cvvec,tpaddvar=(/dble(edsnpath),tptime,tpdis,kelvin0/),ifonlycv=.true.)
        else
           call colvar_calandrec(cvvec,tpaddvar=(/dble(edsnpath),tptime,tpdis,tpdis2,kelvin0/),ifonlycv=.true.)
        end if
     else if ( procid == (procnum-1) ) then
        if ( edsnpc2 == 0 ) then
           write(ioedsinfo,"(f10.3,f8.3,i6,2f10.3)") nowtime,tptime,edsnpath,tpdis,kelvin0
        else
           write(ioedsinfo,"(f10.3,f8.3,i6,3f10.3)") nowtime,tptime,edsnpath,tpdis,tpdis2,kelvin0
        end if
        flush(ioedsinfo)
     end if

     if ( mod(savetrajstep,samplestep) /= 0 ) call errormsg("step error ",int1=savetrajstep,int2=samplestep)
     edspathstep = edspathstep + samplestep

     ifend=.false.
     if ( step0 > 0 ) then
        if ( (nowstep-lastupdatestep) > step0 ) then
           ifend=.true.
        else if ( tpdis < edsdislow ) then
           if ( (edsnpc2 == 0) .or. (tpdis2 < edsdislow2) ) ifend=.true.
        end if
     else
        if ( (nowstep-lastupdatestep) > -step0 ) then
           ifend=.true.
        else if ( tpdis > edsdisup ) then
           if ( (edsnpc2 == 0) .or. (tpdis2 > edsdisup2) ) ifend=.true.
        end if
     end if
     if ( .not. ifeds ) call mpi_bcast(ifend, 1, mpi_logical, procnum-1, mpi_comm_world, ierr)

     if ( ((mod(nowstep,savetrajstep) == 0) .or. (ifend .eqv. .true.)) .and. (procid==(procnum-1)) ) then
        tpcoor=nowcoor; call moldata_moveto_comorcenterbox(tpcoor)
        write(iotraj,"(10f8.3)") tpcoor
        if (nbox>0) write(iotraj,"(3f8.3)") boxsize(1:3); flush(iotraj)
     end if

     if ( ifend .eqv. .true. ) then
        edsnpath=edsnpath+1
        if ( procid == (procnum-1) ) then
           close(iotraj)
           if ( (ifeds .eqv. .false.) .and. (edsboxtime(edswritebox) == 0d0) ) then
              if ( edssimstate == 1 ) then
                 ifcpeds=.true.; lastcpdis=tpdis
              else if ( edssimstate == 2 ) then
                 ifcpeds=.true.; lastcpdis=tpdis2
              end if
           end if
           if ( edsdisup == 0d0 ) then
              nowcoor=inicoor; nowforce=inifrc; precoor=nowcoor
              if (using_pme_potential) call gpu_force_new_neighborlist();
              call gpu_upload_crd(nowcoor); call gpu_upload_frc(nowforce)
           end if
           write(tpfilename,"(i4)") edsnpath; tpfilename=adjustl(tpfilename)
           tpfilename="procinfo/eds_path_"//trim(tpfilename)//".traj"
           open(file=tpfilename, unit=iotraj, action="write")
           write(iotraj,"(a12)") "default_name"; edspathstep=0
           tpcoor=nowcoor; call moldata_moveto_comorcenterbox(tpcoor)
           write(iotraj,"(10f8.3)") tpcoor; if (nbox>0) write(iotraj,"(3f8.3)") boxsize(1:3); flush(iotraj)
        end if

        iedsmdstep = 0; iedsfixstep = edsfixstep; iedsfixstep2 = edsfixstep2; lastupdatestep=nowstep
        if ( edsmdstep == 0 ) iedsfixstep = 0; if ( edsfixstep == 0 ) iedsfixstep2 = 0
        if ( edsdisup > edsdislow ) then
           if ( step0 > 0 ) then
              step0 = -step0up
              if ( (expandkelvin0 > 0d0) .and. (procid==(procnum-1)) ) kelvin0 = expandkelvin0
           else
              step0 = step0low
              if ( (expandkelvin0 > 0d0) .and. (procid==(procnum-1)) ) kelvin0 = inikelvin0
           end if
        else
           if ( edsnpc2 == 0 ) then
              call eds_pcspace_info(tpedsdis=inwalldis)
              if ( ifeds .eqv. .true. ) call colvar_calandrec(cvvec,tpaddvar=(/inwalldis/),ifonlycv=.true.)
           else
              call eds_pcspace_info(tpedsdis=inwalldis,tpedsdis2=inwalldis2)
              if ( ifeds .eqv. .true. ) &
                 call colvar_calandrec(cvvec,tpaddvar=(/inwalldis,inwalldis2/),ifonlycv=.true.)
           end if
           outwalldis=inwalldis; outwalldis2=inwalldis2
           predis=inwalldis; predis2=inwalldis2
        end if
     end if
  end if

  if ( ifcpeds .eqv. .true. ) then
     call colvar_calandrec(cvvec,ifnorecord=.true.)
     if ( ioedswrite == 0 ) then
        call getfreeunit(ioedswrite); open(unit=ioedswrite,file="procinfo/eds_for_flood.txt",action="write")
     end if
     write(ioedswrite,"(2f12.6,10f8.3)") nowtime,edsboxreadtime,tpdis,tpdis2,cvvec(1:ncv)
     flush(ioedswrite)

     edsboxtime(edswritebox)=nowtime
     edscrdbox(:,:,edswritebox) = nowcoor
     call gpu_download_vel(edsvelbox(:,:,edswritebox))
     call gpu_download_frc(edsfrcbox(:,:,edswritebox))
     edswritebox = edswritebox + 1; if ( edswritebox > edsnbox ) edswritebox=1
  end if
end subroutine

subroutine eds_flood_readbox(ifread,tpcrd,tpvel,tpfrc)
  logical,intent(out) :: ifread
  real*8,intent(out),dimension(3,numboxatom) :: tpcrd,tpvel,tpfrc

  ifread=.false.
  if ( (procid==(procnum-1)) .and. (edsboxtime(edsreadbox) > 0d0) ) then
     edsboxreadtime = edsboxtime(edsreadbox)
     edsboxtime(edsreadbox) = 0d0; ifread=.true.
     tpcrd = edscrdbox(:,:,edsreadbox)
     tpvel = edsvelbox(:,:,edsreadbox)
     tpfrc = edsfrcbox(:,:,edsreadbox)
     edsreadbox = edsreadbox + 1; if ( edsreadbox > edsnbox ) edsreadbox=1
  end if
end subroutine

subroutine eds_pcspace_info(tpprecoor,tpedsdis,tpedsdis2,tpedsgrad,tpedsgrad2)
  use nfe_colvar_mod, only : colvar_value, colvar_force
  real*8,intent(in),optional :: tpprecoor(3,numboxatom)
  real*8,intent(out),optional :: tpedsdis,tpedsdis2,tpedsgrad(3,numboxatom),tpedsgrad2(3,numboxatom)
  integer,save :: istart,istop,istart2,istop2
  logical,save :: ifmpiprepared,ifmpiprepared2
  real*8,allocatable,dimension(:,:),save :: edsvecsum,edsvecsum2
  integer :: i,j,icv,ipc,iatom,ierr
  real*8 :: tpcoor(3*max(edsnatom,edsnatom2)),temp,tpedsvec(max(edsnpc,edsnpc2)),tpcvvec(ncv),tpvalue
  real*8 :: tprotmat(3,3),tpgrad(3*max(edsnatom,edsnatom2)),tpedsdissq,tpedsdissq2

  if ( present(tpedsdis) ) then
     tpedsdis=0d0

     if ( edsnatom > 0 ) then

        if ( .not. allocated(edsvecsum) ) then
           allocate(edsvecsum(3,edsnpc)); edsvecsum=0d0
           do ipc=1,edsnpc
              do i=1,edsnatom
                 edsvecsum(1:3,ipc) = edsvecsum(1:3,ipc) + edsvec(3*i-2:3*i,ipc)
              end do
           end do
           edsvecsum = edsvecsum/dble(edsnatom)
        end if
        do i=1, edsnatom
           if ( present(tpprecoor) ) then
              tpcoor(3*i-2:3*i) = tpprecoor(1:3, edsatoms(i))
           else
              tpcoor(3*i-2:3*i) = nowcoor(1:3, edsatoms(i))
           end if
        end do
        call moldata_lRMSD(edsnatom*3,tpcoor,edsrefcoor,temp,ifrotatecoor1=.true.,ifmovecoor2=.false.,rotatemat=tprotmat)
        if ( ifeds .eqv. .false. ) then
           if ( ifmpiprepared .eqv. .false. ) then
              j=(edsnpc/procnum); ifmpiprepared=.true.
              istart=procid*j+1; istop=istart+j-1; if ( procid == (procnum-1) ) istop=edsnpc
           end if
           call mpi_bcast(tpcoor, edsnatom*3, mpi_double_precision, procnum-1, mpi_comm_world, ierr)
           call mpi_bcast(tprotmat, 9, mpi_double_precision, procnum-1, mpi_comm_world, ierr)
        else
           istart=1; istop=edsnpc
        end if
        do ipc=istart,istop
           tpedsvec(ipc) = dot_product(edsvec(1:3*edsnatom,ipc),tpcoor(1:3*edsnatom)) - xi0(ipc)
        end do
        tpedsdissq = dot_product(tpedsvec(istart:istop),tpedsvec(istart:istop))
        tpedsdis = sqrt(tpedsdissq)
        if ( ifeds .eqv. .false. ) then
           call mpi_reduce(tpedsdissq,temp,1,mpi_double_precision, MPI_SUM, procnum-1, mpi_comm_world, ierr)
           tpedsdissq = temp; tpedsdis = sqrt(tpedsdissq) 
           call mpi_bcast(tpedsdis, 1, mpi_double_precision, procnum-1, mpi_comm_world, ierr)
        end if

        if ( present(tpedsgrad) ) then
           tpgrad=0d0
           do ipc=istart,istop
              temp = tpedsvec(ipc)/tpedsdis
              do i=1,edsnatom
                 do j=1,3
                    tpgrad(3*i-2:3*i) = tpgrad(3*i-2:3*i) + temp*  &
                          ( edsvec(3*i-3+j,ipc)*tprotmat(j,1:3) - edsvecsum(j,ipc)*tprotmat(j,1:3) )
                 end do
              end do
           end do
           if ( ifeds .eqv. .false. ) then
              call mpi_reduce(tpgrad(1:edsnatom*3), tpcoor(1:edsnatom*3), 3*edsnatom, &
                     mpi_double_precision, MPI_SUM, procnum-1, mpi_comm_world, ierr)
              tpgrad(1:edsnatom*3) = tpcoor(1:edsnatom*3)
           end if
           tpedsgrad=0d0
           do i=1,edsnatom
              tpedsgrad(1:3,edsatoms(i)) = tpgrad(3*i-2:3*i)
           end do
        end if


     else if ( edsnatom == 0 ) then

        if ( present(tpprecoor) ) then
           if ( edscv <= nfencv ) then
              tpvalue= colvar_value(nfecvs(edscv), tpprecoor)
           else
              call colvar_fsainfo(edscv-nfencv, tpprecoor, tpvalue)
           end if
        else
           if ( edscv <= nfencv ) then
              tpvalue = colvar_value(nfecvs(edscv), nowcoor)
           else
              call colvar_fsainfo(edscv-nfencv, nowcoor, tpvalue)
           end if
        end if
        tpedsdis = abs(tpvalue - cvxi0)
        if ( present(tpedsgrad) ) then
           tpedsgrad=0d0; nowforce=0d0
           if ( present(tpprecoor) ) then
              if ( edscv <= nfencv ) then
                 call colvar_force(nfecvs(edscv), tpprecoor, -1.0d0, nowforce)
              else
                 call colvar_fsainfo(edscv-nfencv, tpprecoor, temp, tpgrad=nowforce, tpfactor=-1d0)
              end if
           else
              if ( edscv <= nfencv ) then
                 call colvar_force(nfecvs(edscv), nowcoor, -1.0d0, nowforce)
              else
                 call colvar_fsainfo(edscv-nfencv, nowcoor, temp, tpgrad=nowforce, tpfactor=-1d0)
              end if
           end if
           temp = (tpvalue-cvxi0)/tpedsdis
           do i=1,cvnatom
              iatom = cvatoms(i)
              tpedsgrad(1:3,iatom) = tpedsgrad(1:3,iatom) - temp*nowforce(1:3,iatom)
           end do
        end if

     end if

  end if


  if ( present(tpedsdis2) ) then
     tpedsdis2=0d0
 
     if ( edsnatom2 > 0 ) then

       if ( .not. allocated(edsvecsum2) ) then
           allocate(edsvecsum2(3,edsnpc2)); edsvecsum2=0d0
           do ipc=1,edsnpc2
              do i=1,edsnatom2
                 edsvecsum2(1:3,ipc) = edsvecsum2(1:3,ipc) + edsvec2(3*i-2:3*i,ipc)
              end do
           end do
           edsvecsum2 = edsvecsum2/dble(edsnatom2)
        end if
        do i=1, edsnatom2
           if ( present(tpprecoor) ) then
              tpcoor(3*i-2:3*i) = tpprecoor(1:3, edsatoms2(i))
           else
              tpcoor(3*i-2:3*i) = nowcoor(1:3, edsatoms2(i))
           end if
        end do
        call moldata_lRMSD(edsnatom2*3,tpcoor,edsrefcoor2,temp,ifrotatecoor1=.true.,ifmovecoor2=.false.,rotatemat=tprotmat)
        if ( ifeds .eqv. .false. ) then
           if ( ifmpiprepared2 .eqv. .false. ) then
              j=(edsnpc2/procnum); ifmpiprepared2=.true.
              istart2=procid*j+1; istop2=istart2+j-1; if ( procid == (procnum-1) ) istop2=edsnpc2
           end if
           call mpi_bcast(tpcoor, edsnatom2*3, mpi_double_precision, procnum-1, mpi_comm_world, ierr)
           call mpi_bcast(tprotmat, 9, mpi_double_precision, procnum-1, mpi_comm_world, ierr)
        else
           istart2=1; istop2=edsnpc2
        end if
        do ipc=istart2,istop2
           tpedsvec(ipc) = dot_product(edsvec2(1:3*edsnatom2,ipc),tpcoor(1:3*edsnatom2)) - xi02(ipc)
        end do
        tpedsdissq2 = dot_product(tpedsvec(istart2:istop2),tpedsvec(istart2:istop2))
        tpedsdis2 = sqrt(tpedsdissq2)
        if ( ifeds .eqv. .false. ) then
           call mpi_reduce(tpedsdissq2,temp,1,mpi_double_precision, MPI_SUM, procnum-1, mpi_comm_world, ierr)
           tpedsdissq2 = temp; tpedsdis2 = sqrt(tpedsdissq2)
           call mpi_bcast(tpedsdis2, 1, mpi_double_precision, procnum-1, mpi_comm_world, ierr)
        end if
        if ( present(tpedsgrad2) ) then
           tpgrad=0d0
           do ipc=istart2,istop2
              temp = tpedsvec(ipc)/tpedsdis2
              do i=1,edsnatom2
                 do j=1,3
                    tpgrad(3*i-2:3*i) = tpgrad(3*i-2:3*i) + temp*  &
                           ( edsvec2(3*i-3+j,ipc)*tprotmat(j,1:3) - edsvecsum2(j,ipc)*tprotmat(j,1:3) )
                 end do
              end do
           end do
           if ( ifeds .eqv. .false. ) then
              call mpi_reduce(tpgrad(1:edsnatom2*3), tpcoor(1:edsnatom2*3), 3*edsnatom2, &
                    mpi_double_precision, MPI_SUM, procnum-1, mpi_comm_world, ierr)
              tpgrad(1:edsnatom2*3) = tpcoor(1:edsnatom2*3)
           end if
           tpedsgrad2=0d0
           do i=1,edsnatom2
              tpedsgrad2(1:3,edsatoms2(i)) = tpgrad(3*i-2:3*i)
           end do
        end if

     else if ( edsnatom2 == 0 ) then

        if ( present(tpprecoor) ) then
           if ( edscv2 <= nfencv ) then
              tpvalue= colvar_value(nfecvs(edscv2), tpprecoor)
           else
              call colvar_fsainfo(edscv2-nfencv, tpprecoor, tpvalue)
           end if
        else
           if ( edscv2 <= nfencv ) then
              tpvalue = colvar_value(nfecvs(edscv2), nowcoor)
           else
              call colvar_fsainfo(edscv2-nfencv, nowcoor, tpvalue)
           end if
        end if
        tpedsdis2 = abs(tpvalue - cvxi02)
        if ( present(tpedsgrad2) ) then
           tpedsgrad2=0d0; nowforce=0d0
           if ( present(tpprecoor) ) then
              if ( edscv2 <= nfencv ) then
                 call colvar_force(nfecvs(edscv2), tpprecoor, -1.0d0, nowforce)
              else
                 call colvar_fsainfo(edscv2-nfencv, tpprecoor, temp, tpgrad=nowforce, tpfactor=-1d0)
              end if
           else
              if ( edscv2 <= nfencv ) then
                 call colvar_force(nfecvs(edscv2), nowcoor, -1.0d0, nowforce)
              else
                 call colvar_fsainfo(edscv2-nfencv, nowcoor, temp, tpgrad=nowforce, tpfactor=-1d0)
              end if
           end if
           temp = (tpvalue-cvxi02)/tpedsdis2
           do i=1,cvnatom
              iatom = cvatoms(i)
              tpedsgrad2(1:3,iatom) = tpedsgrad2(1:3,iatom) - temp*nowforce(1:3,iatom)
           end do
        end if

     end if

  end if
end subroutine

subroutine eds_pcspace_fix(tpprecoor,tprefdis,tprefdis2)
  real*8,intent(in) :: tpprecoor(3,numboxatom)
  real*8,intent(in),optional :: tprefdis,tprefdis2
  real*8 :: dlambda(3,numboxatom),dlambdaold(3,numboxatom)
  integer :: i,j,k,ipc,icv,iatom,maxiter,iter,ierr
  real*8 :: tpcstvalue,sor,eps,sigma,temp,tempmulti,denominator
  logical :: ifconvergent

  maxiter = 10; sor=1.2d0; eps = 1d-1; ifconvergent=.false.; iter=0

  if ( (.not. present(tprefdis)) .and. (.not. present(tprefdis2)) ) call errormsg("eds_pcspace_fix error: no refdis or refdis2!")
  if ( (present(tprefdis)) .and. (present(tprefdis2)) ) call errormsg("eds_pcspace_fix error: double refdis!")

  if ( present(tprefdis) ) then
     call eds_pcspace_info(tpedsdis=tpcstvalue,tpedsgrad=dlambdaold,tpprecoor=tpprecoor)
  else
     call eds_pcspace_info(tpedsdis2=tpcstvalue,tpedsgrad2=dlambdaold,tpprecoor=tpprecoor)
  end if

  do while ( (ifconvergent .eqv. .false.) .and. (iter < maxiter) )
     iter = iter + 1
     ifconvergent = .true.

     if ( present(tprefdis) ) then
        call eds_pcspace_info(tpedsdis=tpcstvalue,tpedsgrad=dlambda)
        sigma = tprefdis - tpcstvalue
     else
        call eds_pcspace_info(tpedsdis2=tpcstvalue,tpedsgrad2=dlambda)
        sigma = tprefdis2 - tpcstvalue
     end if

     if ( ifeds .eqv. .false. ) call mpi_bcast(sigma, 1, mpi_double_precision, procnum-1, mpi_comm_world, ierr)

     if ( (abs(sigma) < eps) .and. (iter>1) ) cycle
     ifconvergent = .false.; denominator=0.0d0

     if ( (present(tprefdis)) .and. (edsnatom > 0) ) then

        do i=1,edsnatom
           iatom=edsatoms(i)
           do j=1,3
              denominator = denominator + dlambda(j,iatom)*dlambdaold(j,iatom)/(2.0d0*atommass(iatom))
           end do
        end do
        tempmulti = sor*sigma/((deltatime**2)*denominator)
        temp = tempmulti*(deltatime**2)
        do i=1,edsnatom
           iatom=edsatoms(i)
           nowcoor(1:3,iatom) = nowcoor(1:3,iatom) +  temp*dlambdaold(1:3,iatom)/(2.0d0*atommass(iatom))
        end do

     else if ( (present(tprefdis2)) .and. (edsnatom2 > 0) ) then

        do i=1,edsnatom2
           iatom=edsatoms2(i)
           do j=1,3
              denominator = denominator + dlambda(j,iatom)*dlambdaold(j,iatom)/(2.0d0*atommass(iatom))
           end do
        end do
        tempmulti = sor*sigma/((deltatime**2)*denominator)
        temp = tempmulti*(deltatime**2)
        do i=1,edsnatom2
           iatom=edsatoms2(i)
           nowcoor(1:3,iatom) = nowcoor(1:3,iatom) +  temp*dlambdaold(1:3,iatom)/(2.0d0*atommass(iatom))
        end do

     else if ( ( (present(tprefdis)) .and. (edsnatom == 0) ) .or. &
               ( (present(tprefdis2)) .and. (edsnatom2 == 0) )  ) then
        do i=1,cvnatom
           iatom=cvatoms(i)
           do j=1,3
              denominator = denominator + dlambda(j,iatom)*dlambdaold(j,iatom)/(2.0d0*atommass(iatom))
           end do
        end do
        tempmulti = sor*sigma/((deltatime**2)*denominator)
        temp = tempmulti*(deltatime**2)
        do i=1,cvnatom
           iatom=cvatoms(i)
           nowcoor(1:3,iatom) = nowcoor(1:3,iatom) + temp*dlambdaold(1:3,iatom)/(2.0d0*atommass(iatom))
        end do

     else
        call errormsg("wrong input parameters!")
     end if
  end do
end subroutine

end module
