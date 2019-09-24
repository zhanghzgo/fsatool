! Written by Chen Changjun, Huazhong University of Science and Technology, 2017

module moldata
  implicit none
  include 'mpif.h'
  integer,parameter :: iomdin=5
  real*8,dimension(:),allocatable :: atomcharge,atommass
  real*8,parameter :: gasconst=0.0019872065d0, pi = 3.141592653589793d0
  integer :: procid,procnum,samplestep,iosiminfo,numboxatom,numboxres,savetrajstep
  integer :: nowstep,numatom,numres,nowikelvin,nummol,nbox,numalpha,numback,numheavy
  real*8 :: kelvin0,epot,deltatime,boxsize(3),nowtime
  real*8,dimension(:,:),allocatable :: nowcoor,nowvel,nowforce,refcoor
  integer,dimension(:),allocatable :: resindex,calpha,resmolindex,molindex,backbone
  character*4,dimension(:),allocatable :: atomname,resname

  integer,dimension(:),allocatable :: replicankelvin
  logical :: iftemd,ifflood,ifsmd,ifeds
contains

subroutine moldata_init()
  use file_io_mod
  integer :: imol,ires,iatom,i,j,ierr,tpnumres,tpnumatom,tparray(numboxatom)
  real*8 :: tpchargesum,tpvec(3),temp
  logical :: iffind

  allocate(calpha(numres)); calpha=0; numalpha=0
  do ires=1,numres
     iffind = .false.
     do iatom=resindex(ires)+1,resindex(ires+1)
        if ( (atomname(iatom) == "CA  ") .or. (atomname(iatom) == "P   ") ) then
           iffind=.true.; numalpha=numalpha+1; calpha(ires) = iatom; exit
        end if
     end do
     if ( iffind .eqv. .false. ) exit
  end do

  allocate(backbone(numres*3)); backbone=0; numback=0
  do ires=1,numres
     iffind = .false.
     do iatom=resindex(ires)+1,resindex(ires+1)
        if ( atomname(iatom) == "N   " )  then
           backbone(3*ires-2) = iatom; iffind=.true.
        else if ( atomname(iatom) == "CA  " )  then
           backbone(3*ires-1) = iatom; iffind=.true.
        else if ( atomname(iatom) == "C   " )  then
           backbone(3*ires) = iatom; iffind=.true.
        end if
     end do
     if ( iffind .eqv. .false. ) exit
  end do
  numback = 3*(ires-1)

  numheavy=0
  do ires=1,numres
     do iatom=resindex(ires)+1,resindex(ires+1)
        if ( (atomname(iatom)(1:1) == "N") .or. (atomname(iatom)(1:1) == "C") .or. &
             (atomname(iatom)(1:1) == "O") .or. (atomname(iatom)(1:1) == "S") .or. &
             (atomname(iatom)(1:1) == "P") )  then
           numheavy = numheavy + 1
        end if
     end do
  end do

  if ( procid == 0 ) then

  tpnumatom=0; tpnumres=0
  do ires=1,numres
     tpchargesum = sum(atomcharge(resindex(ires)+1:resindex(ires+1)))
     if ( abs(tpchargesum) < 1.0d-6 ) tpchargesum=0.0d0
     tpnumatom = resindex(ires+1)-resindex(ires); if ( tpnumatom > 1 ) tpnumres = ires
    write(iosiminfo,"(a,i3,3x,a,i5,a5,5x,a,i5,5x,a,f10.3)") "mol ",resmolindex(ires),"res ",ires,resname(ires), &
           "atom number ",tpnumatom,"charge ",tpchargesum
!     do iatom=resindex(ires)+1,resindex(ires+1)
!        write(iosiminfo,"(i5,a5,2f8.3,3(f10.3,a1))") iatom,atomname(iatom),atomcharge(iatom),atommass(iatom), &
!                        (nowcoor(i,iatom),",",i=1,3)
!     end do
  end do

  write(iosiminfo,*)
  write(iosiminfo,"(a)") "Biomolecule Sequence: "
  write(iosiminfo,"(20(a4))") resname(1:numres)
  write(iosiminfo,"(a,i7,5x,a,i5)") "tot atom ",numboxatom,", tot res ",numboxres
  write(iosiminfo,"(3(a,i5,3x),a,f10.3)") "biomolecules info: nmol ",nummol,"nres ",numres,  &
             "natom ",numatom,"tot charge ",sum(atomcharge)

  if ( nbox > 0 ) write(iosiminfo,"(a,3f8.3)") "periodic box size ",boxsize(1:3)

  write(iosiminfo,*); flush(iosiminfo)

  inquire(file="native.crd", exist=iffind)
  if ( iffind .eqv. .true. ) then
     call moldata_read_crd("native.crd",refcoor)
     tpvec=0d0
     do iatom=1, numboxatom
        tpvec(1:3) = tpvec(1:3) + refcoor(1:3,iatom)
     end do
     tpvec(1:3)=tpvec(1:3)/dble(numboxatom)
     do iatom=1, numboxatom
        refcoor(1:3,iatom) = refcoor(1:3,iatom) - tpvec(1:3)
     end do
     if ( sum(abs(refcoor(1:3,1))) < 1d-5 ) call errormsg("refcoor is not prepared!")
  end if

  end if

  call mpi_bcast(refcoor,3*numboxatom,mpi_double_precision,0, mpi_comm_world,ierr)
  call mpi_bcast(numatom,1,mpi_integer,0, mpi_comm_world,ierr)
  call mpi_bcast(numres,1,mpi_integer,0, mpi_comm_world,ierr)

end subroutine

subroutine moldata_upload_crd(tpcoor)
  real*8,intent(in) :: tpcoor(3,numboxatom)
  if ( nbox == 0 ) then
     call gpu_upload_crd(tpcoor)
  else
     call gpu_shuttle_post_data(tpcoor, 0)
  end if
end subroutine

subroutine moldata_download_crd(tpcoor)
  real*8,intent(out) :: tpcoor(3,numboxatom)
  if ( nbox == 0 ) then
     call gpu_download_crd(tpcoor)
  else
     call gpu_shuttle_retrieve_data(tpcoor, 0)
  end if
end subroutine

subroutine moldata_upload_frc_add(tpfrc)
  real*8,intent(in) :: tpfrc(3,numboxatom)
  if ( nbox == 0 ) then
     call gpu_upload_frc_add(tpfrc)
  else
     call gpu_shuttle_post_data(tpfrc, 4)
  end if
end subroutine

subroutine moldata_read_crd(tpfilename,tpcoor)
  character(len=*),intent(in) :: tpfilename
  real*8,intent(out) :: tpcoor(3,numboxatom)
  integer :: iotraj

  call getfreeunit(iotraj); open(file=tpfilename, unit=iotraj, action="read")
  read(iotraj,*)
  read(iotraj,"(10f8.3)") tpcoor
  close(iotraj)
end subroutine

subroutine moldata_read_traj_ascii(tpfilename,tpcoor)
  character(len=*),intent(in) :: tpfilename
  real*8,intent(out) :: tpcoor(3,numboxatom)
  integer,save :: iotraj
  character*100,save :: prefilename

  if ( trim(tpfilename(:)) /= trim(prefilename(:)) ) then
     if (iotraj > 0) close(iotraj)
     call getfreeunit(iotraj); open(file=tpfilename, unit=iotraj, action="read")
     read(iotraj,*)
     prefilename=tpfilename
  end if
  read(iotraj,"(10f8.3)") tpcoor
  if (nbox>0) read(iotraj,*)
end subroutine

subroutine moldata_read_traj_netcdf(tpfilename,tpindex,tpcrd)
  use netcdf
  character(len=*),intent(in) :: tpfilename
  integer,intent(in) :: tpindex
  real*8,intent(out) :: tpcrd(:,:)
  character*100,save :: prefilename
  integer,save :: ncid, tpcoorid, tpframenum, tpatomnum,tptimeid
  integer :: status, dimid
  real*8 :: tptime(1)
  logical,save :: ifopened

  if ( trim(tpfilename(:)) /= trim(prefilename(:)) ) then
     if (ifopened .eqv. .true.) status = nf90_close(ncid)
     ifopened = .true.; prefilename=tpfilename
     status = nf90_open(tpfilename, nf90_nowrite, ncid)
     status = nf90_inq_dimid(ncid, "atom", dimid)
     status = nf90_inquire_dimension(ncid, dimid, len=tpatomnum)
     status = nf90_inq_dimid(ncid, "frame", dimid)
     status = nf90_inquire_dimension(ncid, dimid, len=tpframenum)
     if ( (tpindex > tpframenum) .or. (size(tpcrd,2) /= tpatomnum ) ) then
        print *,"index error! ",tpindex; stop
     end if
     status = nf90_inq_varid(ncid, "coordinates", tpcoorid)
     status = nf90_inq_varid(ncid, "time", tptimeid)
  end if
  status = nf90_get_var(ncid, tpcoorid, tpcrd, start=(/1, 1, tpindex/), count=(/3, tpatomnum, 1/))
!  status = nf90_get_var(ncid, tptimeid, tptime, start=(/tpindex/), count=(/1/))
!  status  = nf90_close(ncid)
end subroutine

subroutine moldata_write_traj(tpfilename,tpcoor)
  use mdin_ctrl_dat_mod, only: ioutfm
  character*100,intent(in),optional :: tpfilename
  real*8,intent(in),optional :: tpcoor(3,numboxatom)
  character*100 :: tpfilename2
  real*8 :: tpcoor2(3,numboxatom)

  if ( present(tpcoor) ) then
     tpcoor2 = tpcoor
  else
     call gpu_download_crd(tpcoor2)
  end if
  call moldata_moveto_comorcenterbox(tpcoor2)
  if ( present(tpfilename) ) then
     tpfilename2=trim(tpfilename)
  else
     write(tpfilename2,"(i4)") procid; tpfilename2=adjustl(tpfilename2)
     tpfilename2="procinfo/traj_"//trim(tpfilename2)//".mdcrd"
  end if
  if ( ioutfm == 0 ) then
     call moldata_write_traj_ascii(tpfilename2,tpcoor2)
  else
     call moldata_write_traj_netcdf(tpfilename2,tpcoor2)
  end if
end subroutine

subroutine moldata_write_traj_ascii(tpfilename,tpcoor)
  character(len=*),intent(in) :: tpfilename
  real*8,intent(in) :: tpcoor(3,numboxatom)
  integer,save :: iotraj
  character*100,save :: prefilename

  if ( trim(tpfilename(:)) /= trim(prefilename(:)) ) then
     if (iotraj > 0) close(iotraj)
     call getfreeunit(iotraj); open(file=tpfilename, unit=iotraj, action="write")
     write(iotraj,"(a12)") "default_name"
     prefilename=tpfilename
  end if
  write(iotraj,"(10f8.3)") tpcoor
  if (nbox>0) write(iotraj,"(3f8.3)") boxsize(1:3); flush(iotraj)
end subroutine

subroutine moldata_write_traj_netcdf(tpfilename,tpcrd)
  use netcdf
  use axis_optimize_mod, only : axis_flipback_ords
  use pbc_mod, only : pbc_alpha,pbc_beta,pbc_gamma, pbc_box
  use mdin_ctrl_dat_mod, only : ntwprt
  character(len=*),intent(in) :: tpfilename
  real*8,intent(in) :: tpcrd(:,:)
  character*100,save :: prefilename
  integer,save :: tpinistep, tpframenum, tpatomnum, tpcoorid, tptimeid, ncid, tpframeid,tpcelllengthsid,tpcellanglesid
  logical,save :: ifopened
  integer :: status,tpatomid,tpspatialid,tpaxisid,tpcellspatialid,tplabelid,tpcellangularid
  integer :: tpcellspatialvid,tpcellangularvid,ord1,ord2,ord3
  character*1 :: tpspatial(3)

  tpatomnum = numboxatom
  if ( ntwprt > 0 ) tpatomnum = ntwprt

  if ( trim(tpfilename(:)) /= trim(prefilename(:)) ) then
     if (ifopened .eqv. .true.) status = nf90_close(ncid)
     ifopened = .true.; prefilename=tpfilename; tpframenum=0; tpinistep=nowstep
     status = nf90_create(tpfilename, nf90_clobber, ncid)
     status = nf90_def_dim(ncid, "frame", NF90_UNLIMITED, tpframeid)
     status = nf90_def_dim(ncid, "spatial", 3, tpaxisid)
     status = nf90_def_dim(ncid, "atom", tpatomnum, tpatomid)

     if ( nbox > 0 ) then
     status = nf90_def_dim(ncid, "cell_spatial", 3, tpcellspatialid)
     status = nf90_def_dim(ncid, "label", 5, tplabelid)
     status = nf90_def_dim(ncid, "cell_angular", 3, tpcellangularid)
     end if

     status = nf90_def_var(ncid, "time", nf90_float, tpframeid, tptimeid)
     status = nf90_put_att(ncid, tptimeid, "units", "picosecond")
     status = nf90_def_var(ncid, "coordinates", nf90_float, (/tpaxisid,tpatomid,tpframeid/), tpcoorid)
     status = nf90_put_att(ncid, tpcoorid, "units", "angstrom")

     if ( nbox > 0 ) then
     status = nf90_def_var(ncid, "cell_lengths", nf90_float, (/tpcellspatialid,tpframeid/), tpcelllengthsid)
     status = nf90_put_att(ncid, tpcelllengthsid, "units", "angstrom")
     status = nf90_def_var(ncid, "cell_angles", nf90_float, (/tpcellangularid,tpframeid/), tpcellanglesid)
     status = nf90_put_att(ncid, tpcellanglesid, "units", "degree")
     end if

     status = nf90_def_var(ncid, "spatial", nf90_char, tpaxisid, tpspatialid)

     if ( nbox > 0 ) then
     status = nf90_def_var(ncid, "cell_spatial", nf90_char, (/tpcellspatialid/), tpcellspatialvid)
     status = nf90_def_var(ncid, "cell_angular", nf90_char, (/tplabelid,tpcellangularid/), tpcellangularvid)
     end if

     status = nf90_put_att(ncid, nf90_global, "title", "default_name")
     status = nf90_put_att(ncid, nf90_global, "application", "AMBER")
     status = nf90_put_att(ncid, nf90_global, "program", "pmemd")
     status = nf90_put_att(ncid, nf90_global, "programVersion", "16.0")
     status = nf90_put_att(ncid, nf90_global, "Conventions", "AMBER")
     status = nf90_put_att(ncid, nf90_global, "ConventionVersion", "1.0")
     status = nf90_enddef(ncid)

     status = nf90_put_var(ncid, tpspatialid, (/"x","y","z"/),start=(/1/), count=(/3/))

     if ( nbox > 0 ) then   
     status = nf90_put_var(ncid, tpcellangularvid, (/"a","l","p","h","a", &
                 "b","e","t","a"," ", "g","a","m","m","a" /), &
                     start=(/1,1/), count=(/5,3/))
     status = nf90_put_var(ncid, tpcellspatialvid, (/"a","b","c"/),start=(/1/), count=(/3/))
     end if
  end if
  tpframenum = tpframenum + 1
  status = nf90_redef(ncid)
  status = nf90_def_dim(ncid, "frame", tpframenum, tpframeid) 
  status = nf90_enddef(ncid)
  status = nf90_put_var(ncid, tptimeid, (/dble(nowstep-tpinistep)*deltatime/), start=(/tpframenum/), count=(/1/))
  status = nf90_put_var(ncid, tpcoorid, tpcrd, start=(/1, 1, tpframenum/), count=(/3, tpatomnum, 1/))

  if ( nbox > 0 ) then
  ord1 = axis_flipback_ords(1)
  ord2 = axis_flipback_ords(2)
  ord3 = axis_flipback_ords(3)
  status = nf90_put_var(ncid, tpcelllengthsid, (/ pbc_box(ord1), pbc_box(ord2), pbc_box(ord3) /), &
           start=(/1,tpframenum/), count=(/3,1/))
  status = nf90_put_var(ncid, tpcellanglesid, (/pbc_alpha, pbc_beta, pbc_gamma/), &
           start=(/1,tpframenum/), count=(/3,1/))
  end if

  status = nf90_sync(ncid)
end subroutine

subroutine moldata_writecrd(tpfilename,tpcoor)
  character(len=*),intent(in) :: tpfilename
  real*8,intent(in) :: tpcoor(3,numboxatom)
  integer :: iofile

  call getfreeunit(iofile); open(file=tpfilename, unit=iofile, action="write")
  write(iofile,"(a12)") "default_name"
  write(iofile,"(i6)") numboxatom
  write(iofile,"(6f12.7)") tpcoor(1:3,1:numboxatom)
  write(iofile,"(6f12.7)") boxsize(1:3),90.0d0,90.0d0,90.0d0
  close(iofile)
end subroutine

subroutine errormsg(msg,real1,real2,real3,int1,int2,int3)
  implicit none
  character(len=*),intent(in) :: msg
  real*8,intent(in),optional :: real1,real2,real3
  integer,intent(in),optional :: int1,int2,int3
  integer :: ierr
  write (iosiminfo, "(a,a)") "Error: ",msg
  if ( present(real1) ) write(iosiminfo,"(a,f10.3)") "real1= ",real1
  if ( present(real2) ) write(iosiminfo,"(a,f10.3)") "real2= ",real2
  if ( present(real3) ) write(iosiminfo,"(a,f10.3)") "real3= ",real3
  if ( present(int1) ) write(iosiminfo,"(a,i10)") "int1= ",int1
  if ( present(int2) ) write(iosiminfo,"(a,i10)") "int2= ",int2
  if ( present(int3) ) write(iosiminfo,"(a,i10)") "int3= ",int3
  flush(iosiminfo)
  call mpi_finalize(ierr)
  stop
end subroutine

subroutine getfreeunit(tpunit)
  implicit none
  integer tpunit
  logical ifused

  tpunit = 10; ifused = .true.
  do
     if ( ifused .eqv. .false. ) exit
     tpunit = tpunit + 1
     if (tpunit > 99) call errormsg ("unit has been used ",int1=tpunit)
     inquire (unit=tpunit,opened=ifused)
  end do
end subroutine

subroutine moldata_getpotandfrc(tpnatom,tpepot,crd,frc)
  use runmd_mod, only : runmd_interface_getpotandfrc
  implicit none
  integer,intent(in) :: tpnatom
  real*8,intent(out) :: tpepot
  real*8,dimension(3,tpnatom) :: crd, frc

  call runmd_interface_getpotandfrc(tpnatom,nowstep,crd,frc,tpepot)
end subroutine

subroutine rotateVector(vector, direction, angle)
  implicit none
  real*8,intent(in) :: direction(3),angle
  real*8,intent(inout) :: vector(3)
  real*8 :: matrix(3,3)
  real*8 n1,n2,n3,ct,st,temp

    n1 = direction(1); n2 = direction(2); n3 = direction(3)
    ct=cos(angle); st=sin(angle)

    matrix(1,1)=n1**2 + (1.0-n1**2)*ct
    matrix(1,2)=n1*n2*(1.0 - ct) + n3*st
    matrix(1,3)=n1*n3*(1.0 - ct) - n2*st

    matrix(2,1)=n1*n2*(1.0 - ct) - n3*st
    matrix(2,2)=n2**2 + (1.0-n2**2)*ct
    matrix(2,3)=n2*n3*(1.0 - ct) + n1*st

    matrix(3,1)=n1*n3*(1.0 - ct) + n2*st
    matrix(3,2)=n2*n3*(1.0 - ct) - n1*st
    matrix(3,3)=n3**2 + (1.0-n3**2)*ct

    call MatrixDotVector(3,matrix, vector)
end subroutine

subroutine MatrixDotVector(ndim,matrix,vector)
  implicit none
    integer, intent(in) :: ndim
    real*8,intent(in):: matrix(ndim,ndim)
    real*8,intent(inout):: vector(ndim)
    real*8 :: temp,tempvector(ndim)
    integer i,j
    tempvector=vector
    do i=1,ndim
       temp=0.0
       do j=1,ndim
          temp = temp + matrix(i,j)*tempvector(j)
       end do
       vector(i) = temp
     end do
end subroutine 

subroutine cross_product ( A, B, C)
   real*8, dimension(1:3), intent(in) :: A, B
   real*8, dimension(1:3), intent(out) :: C

   C(1) = A(2) * B(3) - A(3) * B(2)
   C(2) = A(3) * B(1) - A(1) * B(3)
   C(3) = A(1) * B(2) - A(2) * B(1)
end subroutine

subroutine moldata_moveto_comorcenterbox(tpcoor)
  real*8,intent(inout) :: tpcoor(3,numboxatom)
  integer :: imol,ires,i,j,k,tpnbox
  real*8 :: tpvec(3),tpdis,tpmass,tpbias,tpcenter(3)

  if ( nbox > 0 ) then
     do imol=1,nummol
        call moldata_get_com(tpcoor,tpvec,tpmols=(/imol/))
        if ( imol == 1 ) then
!           tpcenter(1:3) = tpvec(1:3)
            tpcenter(1:3) = boxsize(1:3)*0.5d0
            cycle
        end if
        do j=1,3
           tpdis = tpvec(j) - tpcenter(j)
           tpbias = (-anint(tpdis/boxsize(j)))*boxsize(j)
           do ires=molindex(imol)+1,molindex(imol+1)
              do k=resindex(ires)+1,resindex(ires+1)
                 tpcoor(j,k)=tpcoor(j,k) + tpbias
              end do
           end do
        end do
     end do
  else
     tpvec=0d0
     do i=1,numboxatom
        tpvec(1:3) = tpvec(1:3) + tpcoor(1:3,i)
     end do
     tpvec(1:3) = tpvec(1:3)/dble(numboxatom)
     do i=1,numboxatom
        tpcoor(1:3,i) = tpcoor(1:3,i) - tpvec(1:3)
     end do
  end if
end subroutine

subroutine moldata_setccharge_formdec(tpcharge)
  real*8,intent(inout) :: tpcharge(numboxatom)
  integer :: ires,iatom,iofile,imol,istart,istop
  real*8 :: temp,tpscale,mdeceps
  character*8 :: tpresstr
  character*50 :: tpstring

  mdeceps=1.78d0
  if ( procid == 0 ) then
     write(iosiminfo,"(a,f10.3)") "USE MDEC, eps: ",mdeceps
     call getfreeunit(iofile); open(unit=iofile,file="parmmd_scale_neutralgroup.in",action="write")
  end if

  tpscale = (1.0d0/sqrt(mdeceps))
  do imol=1,nummol
     if ( resname(molindex(imol)+1) == "WAT " ) exit
     if ( procid == 0 ) write(iofile,"(a,i8)") "# mol ",imol
     do ires=molindex(imol)+1,molindex(imol+1)
         istart = resindex(ires)+1; istop = resindex(ires+1)
         temp = sum(tpcharge(istart:istop))
         if ( abs(temp) > 0.01d0 ) then
            tpcharge(istart:istop) = tpcharge(istart:istop)*tpscale
         else
            if ( procid == 0 ) then
               write(tpresstr,"(i8)") ires; tpresstr = adjustl(tpresstr)
               do iatom=istart,istop
                  write(tpstring,"(a,f8.4)") "change charge :"//trim(tpresstr)//"@"//atomname(iatom), &
                      tpcharge(iatom)*1.334d0/18.2223
                  write(iofile,"(a)") tpstring
               end do
            end if
         end if
     end do
  end do
  if ( procid == 0 ) then
     write(iofile,"(a)") "outparm pbsa_com_prmtop"
     write(iofile,"(a)") "quit"
     close(iofile)
  end if
end subroutine

subroutine moldata_get_com(tpcoor,tppos,tpmols,tpres,tpatoms,ifcalpha)
  real*8,intent(in) :: tpcoor(3,numatom)
  real*8,intent(out) :: tppos(3)
  integer,intent(in),dimension(:),optional :: tpmols,tpatoms,tpres
  logical,intent(in),optional :: ifcalpha
  integer :: imol,ires,i,k,tpnmol,istart,istop,tpnatom,tpatomarray(numboxatom)
  real*8 :: tpvec(3),tpmass

  if ( present(tpmols) ) then
     call moldata_getatoms(tpnatom,tpatomarray,tpmols=tpmols, ifcalpha=ifcalpha)
  else if ( present(tpres) ) then
     call moldata_getatoms(tpnatom,tpatomarray,tpres=tpres, ifcalpha=ifcalpha)
  else if ( present(tpatoms) ) then
     tpnatom = size(tpatoms); tpatomarray(1:tpnatom)=tpatoms(1:tpnatom)
  end if

  tpvec(1:3)=0.0d0; tpmass=0d0
  do i=1,tpnatom
     k = tpatomarray(i)
     tpvec(1:3)=tpvec(1:3)+tpcoor(1:3,k)*atommass(k); tpmass = tpmass + atommass(k)
  end do
  tppos(1:3)=tpvec(1:3)/tpmass
end subroutine

subroutine moldata_getatoms(tpnatom,tpatoms,tpmols,tpres,ifcalpha)
  integer,intent(out) :: tpnatom
  integer,dimension(numboxatom),intent(out) :: tpatoms
  integer,dimension(:),intent(in),optional :: tpmols,tpres
  logical,intent(in),optional :: ifcalpha
  integer :: i,k,iatom,ires,imol,istart,istop

  tpatoms=0; tpnatom=0
  if ( present(tpmols) ) then
     do i=1,size(tpmols)
        imol = tpmols(i)
        do ires=molindex(imol)+1,molindex(imol+1)
           istart = resindex(ires)+1; istop = resindex(ires+1)
           if ( (present(ifcalpha)) .and. (ifcalpha .eqv. .true.) ) then
              istart=calpha(ires); istop=calpha(ires)
           end if
           do k=istart,istop
              tpnatom = tpnatom + 1; tpatoms(tpnatom) = k
           end do
        end do
     end do
  else if ( present(tpres) ) then
     do i=1,size(tpres)
        ires = tpres(i)
        istart = resindex(ires)+1; istop = resindex(ires+1)
        if ( (present(ifcalpha)) .and. (ifcalpha .eqv. .true.) ) then
           istart=calpha(ires); istop=calpha(ires)
        end if
        do k = istart,istop
           tpnatom = tpnatom + 1; tpatoms(tpnatom) = k
        end do
     end do
  else
     call errormsg("wrong in moldata_getatoms!")
  end if
end subroutine

subroutine moldata_anchor_and_rotate(tpcoor,ifwritecrd,anchormol,anchorres,anchoratom,rotatemol,rotateres,rotateatom,ifcalpha)
  real*8,intent(inout) :: tpcoor(3,numboxatom)
  logical,intent(in),optional :: ifwritecrd,ifcalpha
  integer,intent(in),dimension(:),optional :: anchormol,anchorres,anchoratom,rotatemol,rotateres,rotateatom
  integer :: iatom,ires,tpnanchor,tpanchoratoms(numboxatom),tpnrotate,tprotateatoms(numboxatom)
  real*8 :: tpdis,temp,tpcenter(3),tpvec(3),tpdirect(3),tpangle
 
  tpcenter=0d0; tpnanchor=0; tpanchoratoms=0
  if ( present(anchormol) ) then
     call moldata_getatoms(tpnanchor,tpanchoratoms,tpmols=anchormol,ifcalpha=ifcalpha)
  else if ( present(anchorres) ) then
     call moldata_getatoms(tpnanchor,tpanchoratoms,tpres=anchorres,ifcalpha=ifcalpha)
  else if ( present(anchoratom) ) then
     tpnanchor=size(anchoratom); tpanchoratoms(1:tpnanchor)=anchoratom(1:tpnanchor)
  end if
  call moldata_get_com(tpcoor,tpcenter,tpatoms=tpanchoratoms)

  tpvec=0d0; tpnrotate=0; tprotateatoms=0
  if ( present(rotatemol) ) then
     call moldata_getatoms(tpnrotate,tprotateatoms,tpmols=rotatemol,ifcalpha=ifcalpha)
  else if ( present(rotateres) ) then
     call moldata_getatoms(tpnrotate,tprotateatoms,tpres=rotateres,ifcalpha=ifcalpha)
  else if ( present(rotateatom) ) then
     tpnrotate=size(rotateatom); tprotateatoms(1:tpnrotate)=rotateatom(1:tpnrotate)
  end if
  call moldata_get_com(tpcoor,tpvec,tpatoms=tprotateatoms)

  tpdirect(1:3) = tpvec(1:3) - tpcenter(1:3)
  temp = sqrt(dot_product(tpdirect,tpdirect)); tpdirect=tpdirect/temp
  do iatom=1,numboxatom
     tpcoor(1:3,iatom) = tpcoor(1:3,iatom) - tpcenter(1:3)
  end do
  call cross_product((/1d0,0d0,0d0/),tpdirect,tpvec)
  temp = sqrt(dot_product(tpvec,tpvec)); tpvec=tpvec/temp
  temp = dot_product((/1d0,0d0,0d0/),tpdirect); tpangle=acos(temp)
  do iatom=1,numatom
     call rotateVector(tpcoor(1:3,iatom), tpvec, tpangle)
  end do

  call moldata_get_com(tpcoor,tpcenter,tpatoms=tpanchoratoms)
  call moldata_get_com(tpcoor,tpvec,tpatoms=tprotateatoms)
  write(iosiminfo,"(2(a,i4,3f8.3,5x))") "anchor atoms ",tpnanchor,tpcenter,"rotate atoms ",tpnrotate,tpvec
  flush(iosiminfo)

  if ( (present(ifwritecrd)) .and. (ifwritecrd .eqv. .true.) ) call moldata_writecrd("tp.crd",tpcoor)
end subroutine

subroutine moldata_lRMSD(ndim,coor1,coor2,rmsd,ifrotatecoor1,rotatemat,tpgrad,ifmovecoor2)
  integer,intent(in) :: ndim
  real*8,intent(inout),dimension(ndim) :: coor1,coor2
  real*8,intent(out) :: rmsd
  logical,intent(in),optional :: ifrotatecoor1,ifmovecoor2
  real*8,intent(out),optional :: rotatemat(3,3)
  real*8,intent(out),optional :: tpgrad(ndim)
  integer :: inum,i,j,idim,info,tpnatom
  real*8 :: tpvec(3),temp,umatrix(3,3),vmatrix(3,3),wvector(3),tpmatrix(3,3),tpvec2(ndim)
  real*8 :: originmatrix(3,3),determinant,dsign,tpcoor(ndim),work(15),rotmatdotdx(3)
  
  if ( mod(ndim,3) /= 0 ) call errormsg("wrong ndim in lRMSD!")
  tpvec(1:3)=0.0d0; tpnatom=ndim/3
  do i=1,tpnatom
     tpvec(1:3)=tpvec(1:3)+coor1(3*i-2:3*i)
  end do
  tpvec(1:3)=tpvec(1:3)/dble(tpnatom)
  do i=1,tpnatom
     coor1(3*i-2:3*i)=coor1(3*i-2:3*i)-tpvec(1:3)
  end do

  if ( (present(ifmovecoor2)) .and. (ifmovecoor2 .eqv. .true.) ) then
     tpvec(1:3)=0.0d0
     do i=1,tpnatom
        tpvec(1:3)=tpvec(1:3)+coor2(3*i-2:3*i)
     end do
     tpvec(1:3)=tpvec(1:3)/dble(tpnatom)
     do i=1,tpnatom
        coor2(3*i-2:3*i)=coor2(3*i-2:3*i)-tpvec(1:3)
     end do
  end if

  originmatrix(:,:)=0.0d0
  do i=1,3
     do j=1,3
        do idim=1,ndim,3
           originmatrix(i,j) = originmatrix(i,j) + coor1(idim+i-1)*coor2(idim+j-1)
        end do
     end do
  end do
  call getMatrixDeterminant(3,originmatrix,determinant); dsign=sign(1.0d0, determinant)

  wvector=0d0; umatrix=0d0; vmatrix=0d0; work=0d0
  call dgesvd("A", "A", 3, 3, originmatrix, 3, wvector, umatrix, 3, vmatrix, 3, work, 15, info)

  tpvec(1:3)=(/1.0d0,1.0d0,dsign/)
  do i=1,3
     do j=1,3
        temp=0.0d0
        do inum=1,3
           temp = temp + vmatrix(inum,i)*tpvec(inum)*umatrix(j,inum)
        end do
        tpmatrix(i,j)=temp
     end do
  end do
  rmsd=0.0d0; tpcoor=0d0
  do idim=1,ndim,3
     do i=1,3
        tpcoor(idim+i-1) = dot_product(tpmatrix(i,1:3),coor1(idim:idim+2))
        rmsd = rmsd + (tpcoor(idim+i-1) - coor2(idim+i-1))**2 
     end do
  end do
  rmsd=sqrt(rmsd/dble(tpnatom))
  if ( (present(ifrotatecoor1)) .and. (ifrotatecoor1 .eqv. .true.)) coor1=tpcoor
  if ( present(rotatemat) .eqv. .true. ) rotatemat = tpmatrix

  if ( present(tpgrad) ) then
     tpvec2(1:ndim) = 2.0d0*(tpcoor(1:ndim) - coor2(1:ndim))
     tpvec=0d0; tpgrad=0d0
     do i=1,tpnatom
        tpvec(1:3) = tpvec(1:3) + tpvec2(3*i-2:3*i)*2.0d0
     end do
     do j=1,3
        rotmatdotdx(j) = dot_product(tpmatrix(1:3,j), tpvec(1:3))/dble(tpnatom)
     end do
     do i=1,tpnatom
        do j=1,3
           tpgrad(3*(i-1)+j) = dot_product(tpmatrix(1:3,j),tpvec2(3*i-2:3*i)) - rotmatdotdx(j)
        end do
     end do
     tpgrad = 0.5d0/(rmsd*dble(tpnatom))*tpgrad
  end if

contains
subroutine getMatrixDeterminant(ndim,amatrix,determinant)
  implicit none
  integer,intent(in):: ndim
  real*8,intent(in),dimension(ndim,ndim):: amatrix
  real*8,intent(out) :: determinant
  integer:: idim,jdim
  real*8:: mij,tpmatrix(ndim,ndim)

  tpmatrix=amatrix
  do jdim=1,ndim
     do idim=jdim+1,ndim
        if ( tpmatrix(jdim,jdim) == 0.0d0 ) cycle
        mij=tpmatrix(idim,jdim)/tpmatrix(jdim,jdim)
        tpmatrix(idim,jdim:ndim) = tpmatrix(idim,jdim:ndim) - mij*tpmatrix(jdim,jdim:ndim)
     end do
  end do
  determinant=1.0d0
  do idim=1,ndim
     determinant=determinant*tpmatrix(idim,idim)
  end do
end subroutine
end subroutine

end module
