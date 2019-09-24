! Written by Chen Changjun, Huazhong University of Science and Technology, 2017

module colvar
  use moldata
  use nfe_colvar_mod, only : colvar_t
  implicit none

  integer :: ncv,nfencv,cvnatom
  integer,dimension(:),allocatable :: cvatoms
  real*8,dimension(:),allocatable :: cvlowrange,cvuprange
  type(colvar_t), allocatable :: nfecvs(:)
  type fsacvs_t
     integer :: type = -1, ni,nr,ngroup,cv_index
     integer, pointer :: i(:) => null(), groups(:,:) => null()
     double precision, pointer :: r(:) => null()
     character*30 :: cv_file
  end type
  type(fsacvs_t), allocatable :: fsacvs(:)
  integer,parameter :: CV_COM_POS=1, CV_COM_DIS=2, CV_COM_POSAXIS=3, CV_PCA=4, CV_RMSD=5, CV_DISSUM=6
contains

subroutine colvar_init()
  use file_io_mod
  use nfe_colvar_mod, only : colvar_force
  logical :: iffind

  integer,parameter :: maxcvatom=500
  integer,dimension(maxcvatom) :: cv_i
  real*8,dimension(maxcvatom*3) :: cv_r
  integer :: i, j, tpni, iatom, icv, cv_ni, cv_nr,ifind,cv_index, tparray(numboxatom)
  character*30 :: cv_type,cv_file
  real*8 :: cv_min, cv_max, tpvec(3)
  namelist / fsacolvar / cv_type, cv_i, cv_r, cv_index, cv_file, cv_min, cv_max

  icv = 0; rewind(iomdin)
  do
     call nmlsrc('fsacolvar', iomdin, ifind)
     if (ifind == 0 ) exit
     read(iomdin,'(a80)')
     icv = icv + 1
  end do
  ncv = nfencv + icv;  allocate(fsacvs(ncv-nfencv))
  cv_r(1:nfencv)=cvlowrange(1:nfencv); deallocate(cvlowrange)
  allocate(cvlowrange(ncv)); cvlowrange(1:nfencv)=cv_r(1:nfencv)
  cv_r(1:nfencv)=cvuprange(1:nfencv); deallocate(cvuprange)
  allocate(cvuprange(ncv)); cvuprange(1:nfencv)=cv_r(1:nfencv)

  icv = 1; rewind(iomdin);
  do icv=1,ncv-nfencv
     call nmlsrc('fsacolvar', iomdin, ifind); cv_i(:)=-10000; cv_r(:)=-10000.0d0
     read(iomdin,nml=fsacolvar)
     do i=1,maxcvatom
        if ( cv_i(i) == -10000 ) exit
     end do
     fsacvs(icv)%ni=i-1; allocate(fsacvs(icv)%i(fsacvs(icv)%ni))
     do i=1,maxcvatom*3
        if ( cv_r(i) == -10000.0d0 ) exit
     end do
     fsacvs(icv)%nr=i-1; allocate(fsacvs(icv)%r(fsacvs(icv)%nr))
     fsacvs(icv)%i(1:fsacvs(icv)%ni) = cv_i(1:fsacvs(icv)%ni)
     fsacvs(icv)%r(1:fsacvs(icv)%nr) = cv_r(1:fsacvs(icv)%nr)
     if ( cv_type == 'COM_POS' ) then
        fsacvs(icv)%type = CV_COM_POS
     else if ( cv_type == 'COM_DIS' ) then
        fsacvs(icv)%type = CV_COM_DIS
     else if ( cv_type == 'COM_POSAXIS' ) then
        fsacvs(icv)%type = CV_COM_POSAXIS
     else if ( cv_type == 'PCA' ) then
        fsacvs(icv)%type = CV_PCA; fsacvs(icv)%cv_index=cv_index; fsacvs(icv)%cv_file=cv_file
        cvlowrange(icv+nfencv) = cv_min; cvuprange(icv+nfencv)= cv_max
     else if ( cv_type == 'RMSD' ) then
        fsacvs(icv)%type = CV_RMSD
        if (fsacvs(icv)%nr /= fsacvs(icv)%ni*3 ) call errormsg("RMSD definition error ",int1=fsacvs(icv)%ni,int2=fsacvs(icv)%nr)
        tpvec=0d0
        do i=1,fsacvs(icv)%ni
           tpvec(1:3) = tpvec(1:3) + fsacvs(icv)%r(3*i-2:3*i)
        end do
        tpvec(1:3)=tpvec(1:3)/dble(fsacvs(icv)%ni)
        do i=1,fsacvs(icv)%ni
           fsacvs(icv)%r(3*i-2:3*i) = fsacvs(icv)%r(3*i-2:3*i) - tpvec(1:3)
        end do
     else if ( cv_type == 'DISSUM' ) then
        if ( mod(fsacvs(icv)%ni,2) /= 0 ) call errormsg("DISSUM: atom number is not even!")
        fsacvs(icv)%type = CV_DISSUM
     else
        call errormsg("colvar type error :"//trim(cv_type))
     end if
     j=0
     do i=1,fsacvs(icv)%ni
        iatom = fsacvs(icv)%i(i)
        if ( (iatom == 0) .or. ((i==fsacvs(icv)%ni)) .and. (iatom>0) ) j=j+1 
     end do
     fsacvs(icv)%ngroup = j; ! if ( j == 0 ) call errormsg("group num error!")
     allocate(fsacvs(icv)%groups(2,j)); fsacvs(icv)%groups=0
     j=0; fsacvs(icv)%groups(1,1)=1
     do i=1,fsacvs(icv)%ni
        iatom = fsacvs(icv)%i(i)
        if ( (iatom == 0) .or. ((i==fsacvs(icv)%ni)) .and. (iatom>0) ) then
           j = j + 1
           if ( (i==fsacvs(icv)%ni) .and. (iatom>0) ) then
              fsacvs(icv)%groups(2,j) = i; exit
           else if ( iatom == 0 ) then
              fsacvs(icv)%groups(2,j) = i-1; if ( i==fsacvs(icv)%ni ) exit
              fsacvs(icv)%groups(1,j+1) = i+1
           end if
        end if
     end do
  end do
  rewind(iomdin); flush(iosiminfo)

  cvnatom=0; tparray=0
  do icv = 1, ncv
     if ( icv <= nfencv ) then
        tpni = size(nfecvs(icv)%i)
     else
        tpni = fsacvs(icv-nfencv)%ni
     end if
     do i=1,tpni
        if ( icv <= nfencv ) then
           iatom = nfecvs(icv)%i(i)
        else
           iatom = fsacvs(icv-nfencv)%i(i)
        end if
        if ( iatom == 0 ) cycle
        iffind=.false.
        do j=1,cvnatom
           if ( iatom == tparray(j) ) then
              iffind=.true.; exit
           end if
        end do
        if ( iffind .eqv. .true. ) then
           cycle
        else
           cvnatom = cvnatom + 1
           tparray(cvnatom) = iatom
        end if
     end do
  end do
  allocate(cvatoms(cvnatom)); cvatoms(1:cvnatom) = tparray(1:cvnatom)
  if ( procid == 0 ) then

     write(iosiminfo,"(a,i5,5x,a,i5,5x,a,i5)") "colvar number ",ncv," nfe cv ",nfencv,"fsa cv ",ncv-nfencv
     write(iosiminfo,*); write(iosiminfo,"(a)" ) "fsa collective variables info"
     do i=nfencv+1,ncv
        icv = i - nfencv
        write(iosiminfo,"(4(a,i4,3x))") "cv ",icv,", group number ",fsacvs(icv)%ngroup, &
              ", indices ",fsacvs(icv)%ni, ", paramters ",fsacvs(icv)%nr
        write(iosiminfo,"(5x,20(a,i3,a,i4,a2,3x))") ("group ",j,", atoms ", &
             fsacvs(icv)%groups(2,j) - fsacvs(icv)%groups(1,j) + 1,", ", j=1,fsacvs(icv)%ngroup)
     end do
     write(iosiminfo,*); flush(iosiminfo)

     write(iosiminfo,"(a,i6)") "cvnatom ",cvnatom
     write(iosiminfo,"(a)") "cvatoms "
     write(iosiminfo,"(18i6)") cvatoms(1:cvnatom)
     write(iosiminfo,*); flush(iosiminfo)
  end if
  call gpu_setup_shuttle_info(cvnatom, 0, cvatoms)

!  additional operations
!  call colvar_get_com(nowcoor,tpvec,tpmols=(/1,2/),ifcalpha=.true.); write(iosiminfo,"(a,10f8.3)") "com1 ",tpvec(:)
!  call colvar_get_com(nowcoor,tpvec,tpmols=(/3,4/),ifcalpha=.true.); write(iosiminfo,"(a,10f8.3)") "com2 ",tpvec(:)
!  call colvar_anchor_and_rotate(nowcoor,ifwritecrd=.true.,anchormol=(/1,2/),rotatemol=(/3,4/),ifcalpha=.true.)
!  call colvar_printcv_calpharmsd((/1,51,0/), skip=3)
!  call colvar_printcv_calpharmsd((/52,102,0/), skip=3)
!  call colvar_printcv_backbonermsd((/1,51,0/))
!  call colvar_printcv_backbonermsd((/52,102,0/))

!  call colvar_printcv_heavyrmsd((/1,numres,0/))
!  call colvar_printcv_backbonenhb(1,numres,1)

!  call colvar_printcv_contacts((/1,numres,0/),4.0d0,tpstride=18)
end subroutine


subroutine colvar_printcv_contacts(tpgs1,tpmarkdis,tpstride,tpgs2,ifcalpha,tpncontact,tpcontactatoms,tpcoor,tpcontactnatom)
  integer,dimension(:),intent(in) :: tpgs1
  real*8,intent(in) :: tpmarkdis
  integer,intent(in),optional :: tpstride
  integer,intent(in),dimension(:),optional :: tpgs2
  logical,intent(in),optional :: ifcalpha
  integer,intent(out),optional :: tpncontact,tpcontactnatom
  integer,intent(out),optional :: tpcontactatoms(numatom)
  real*8,intent(in),optional :: tpcoor(3,numboxatom)
  integer,parameter :: maxncontact=1000
  integer :: istride,igroup,jgroup,tpngroup1,tpngroup2,tpgroups1(2,numres),tpgroups2(2,numres)
  integer :: ncontact,iatom,jatom,i,j,k,m,n,ires,jres,ibond,tpcontacts(4,maxncontact)
  logical :: iffind

  istride=1; tpngroup1 = size(tpgs1)/3; tpgroups1=0
  if ( present(tpstride) ) istride=tpstride
  do igroup=1,tpngroup1
     tpgroups1(1,igroup) = tpgs1(3*igroup-2); tpgroups1(2,igroup) = tpgs1(3*igroup-1)
  end do

  ncontact=0; tpcontacts=0
  if ( .not. present(tpgs2) ) then

     do igroup=1,tpngroup1
        do ires=tpgroups1(1,igroup),tpgroups1(2,igroup)
           do jgroup=igroup,tpngroup1
              do jres=tpgroups1(1,jgroup),tpgroups1(2,jgroup)
                 if ( (jres-ires) >= istride ) call findrescontact()
              end do
           end do
        end do
     end do

  else

     tpngroup2 = size(tpgs2)/3; tpgroups2=0
     do igroup=1,tpngroup2
        tpgroups2(1,igroup) = tpgs2(3*igroup-2); tpgroups2(2,igroup) = tpgs2(3*igroup-1)
     end do

     do igroup=1,tpngroup1
        do ires=tpgroups1(1,igroup),tpgroups1(2,igroup),istride
           do jgroup=1,tpngroup2
              do jres=tpgroups2(1,jgroup),tpgroups2(2,jgroup),istride
                 call findrescontact()
              end do
           end do
        end do
     end do

  end if

  if ( procid == 0 ) then
  write(iosiminfo,*)
  write(iosiminfo,"(a)") "&colvar"
  write(iosiminfo,"(a)") "  cv_type = 'N_OF_BONDS'"
  write(iosiminfo,"(a,i3,a,i2)") "  cv_ni= ",ncontact*2,", cv_nr= ",1
  write(iosiminfo,"(a,1000(i4,a2))") "  cv_i = ",((tpcontacts(iatom,ibond),", ",iatom=1,2),ibond=1,ncontact)
  write(iosiminfo,"(a,f4.1))") "  cv_r = ",tpmarkdis
  write(iosiminfo,"(a)") "/"
  end if

  if ( present(tpncontact) ) tpncontact=ncontact
  if ( present(tpcontactatoms) ) then
     tpcontactatoms=0; k=0; m=0
     do i=1,ncontact
        do j=1,2
           iffind = .false.
           do k=1,m
              if ( tpcontactatoms(k) == tpcontacts(j,i) ) then
                 iffind = .true.; exit
              end if
           end do
           if ( iffind .eqv. .false. ) then
              m = m + 1; tpcontactatoms(m) = tpcontacts(j,i)
              do k=1,m-1
                 if ( tpcontactatoms(m) < tpcontactatoms(k) ) then
                    do n=m,k+1,-1
                       tpcontactatoms(n) = tpcontactatoms(n-1)
                    end do
                    tpcontactatoms(k) = tpcontacts(j,i)
                    exit
                 end if
              end do
           end if 
        end do
     end do
     if ( procid == 0 ) then
        write(iosiminfo,"(a,i6)") "number of contact atoms ",m
        write(iosiminfo,"(18i6)") tpcontactatoms(1:m)
     end if
     if ( present(tpcontactnatom) ) tpcontactnatom = m
  end if

contains
  subroutine findrescontact()
    real*8 :: tpdis,tpmindis,tpvec(3)

    tpmindis = 10000.0d0
    do i=resindex(ires)+1,resindex(ires+1)
       if ( (present(ifcalpha)) .and. (atomname(i)(1:4)/="CA  ") ) cycle
       do j=resindex(jres)+1,resindex(jres+1)
          if ( (present(ifcalpha)) .and. (atomname(j)(1:4)/="CA  ") ) cycle
          if ( present(tpcoor) ) then
             tpvec(1:3) = tpcoor(1:3,i) - tpcoor(1:3,j)
          else
             tpvec(1:3) = nowcoor(1:3,i) - nowcoor(1:3,j)
          end if
          tpdis = sqrt(dot_product(tpvec,tpvec))
          if ( tpdis < tpmindis ) then
             tpmindis = tpdis; iatom = i; jatom=j
          end if
       end do
    end do
    if ( tpmindis < tpmarkdis ) then
       ncontact = ncontact + 1
       if ( ncontact > maxncontact ) call errormsg("too much contacts ",int1=ncontact,int2=maxncontact)
       tpcontacts(1:4,ncontact) = (/iatom,jatom,ires,jres/)
!       write(iosiminfo,"(a,i3,a1,3x,a4,i3,a1,a4,a3,a4,i3,a1,a4)")  &
!          "contacts ",ncontact,":",resname(tpcontacts(3,ncontact)), &
!          tpcontacts(3,ncontact),":",atomname(tpcontacts(1,ncontact))," - ",resname(tpcontacts(4,ncontact)), &
!          tpcontacts(4,ncontact),":",atomname(tpcontacts(2,ncontact))
    end if
  end subroutine
end subroutine

subroutine colvar_printcv_calpharmsd(tpgs, skip)
  integer,intent(in),dimension(:) :: tpgs
  integer,intent(in),optional :: skip
  integer :: igroup,tpngroup,tpnres,ires,index,tpgroups(2,numres),tpskip
  character*150:: tpformat

  tpskip=1; if ( present(skip) ) tpskip=skip
  tpngroup = size(tpgs)/3; tpnres=0
  do igroup=1,tpngroup
     tpgroups(1,igroup) = tpgs(3*igroup-2); tpgroups(2,igroup) = tpgs(3*igroup-1)
     tpnres = tpnres + (tpgroups(2,igroup) - tpgroups(1,igroup) + 1)/tpskip
  end do
  write(iosiminfo,*)
  write(iosiminfo,"(a)") "&colvar"
  write(iosiminfo,"(a)") "  cv_type = 'MULTI_RMSD'"
  write(iosiminfo,"(a,i3,a,i4)") "  cv_ni= ",tpnres+tpngroup,", cv_nr=",(tpnres*3)

  tpformat(1:3)="(a,"; index=4
  do igroup=1,tpngroup 
     write(tpformat(index+(igroup-1)*16:index+(igroup-1)*16+15),"(i2,a14)") &
          (tpgroups(2,igroup)-tpgroups(1,igroup)+1)/tpskip,"(i4,a2),'0, ',"
  end do
  tpformat(index+tpngroup*16-1:index+tpngroup*16-1) = ")"
  write(iosiminfo,fmt=tpformat) "  cv_i = ",((calpha(ires),", ",  &
          ires=tpgroups(1,igroup),tpgroups(2,igroup),tpskip),igroup=1,tpngroup)

  write(iosiminfo,"(a)") "  cv_r = "
  do igroup=1,tpngroup
      write(iosiminfo,"(6(f10.3,a2))") ((nowcoor(index,calpha(ires)),", ",index=1,3), &
                 ires=tpgroups(1,igroup),tpgroups(2,igroup),tpskip)
  end do
  write(iosiminfo,"(a)") "/"
  flush(iosiminfo)
end subroutine

subroutine colvar_printcv_backbonermsd(tpgs)
  integer,intent(in),dimension(:) :: tpgs
  integer :: igroup,tpngroup,tpnres,ires,index,tpgroups(2,numres),i
  character*150:: tpformat

  tpngroup = size(tpgs)/3; tpnres=0
  do igroup=1,tpngroup
     tpgroups(1,igroup) = tpgs(3*igroup-2); tpgroups(2,igroup) = tpgs(3*igroup-1)
     tpnres = tpnres + tpgroups(2,igroup) - tpgroups(1,igroup) + 1
  end do
  write(iosiminfo,*)
  write(iosiminfo,"(a)") "&colvar"
  write(iosiminfo,"(a)") "  cv_type = 'MULTI_RMSD'"
  write(iosiminfo,"(a,i3,a,i4)") "  cv_ni= ",tpnres*3+tpngroup,", cv_nr=",(tpnres*9)

  tpformat(1:3)="(a,"; index=4
  do igroup=1,tpngroup 
     write(tpformat(index+(igroup-1)*17:index+(igroup-1)*17+16),"(i3,a14)") &
          (tpgroups(2,igroup)-tpgroups(1,igroup)+1)*3,"(i4,a2),'0, ',"
  end do
  tpformat(index+tpngroup*17-1:index+tpngroup*17-1) = ")"
  write(iosiminfo,fmt=tpformat) "  cv_i = ",(( (backbone(i),", ", i=3*ires-2,3*ires),  &
          ires=tpgroups(1,igroup),tpgroups(2,igroup)),igroup=1,tpngroup)

  write(iosiminfo,"(a)") "  cv_r = "
  do igroup=1,tpngroup
      write(iosiminfo,"(6(f10.3,a2))") (((nowcoor(index,backbone(i)),", ",index=1,3),i=3*ires-2,3*ires), &
                 ires=tpgroups(1,igroup),tpgroups(2,igroup))
  end do
  write(iosiminfo,"(a)") "/"
  flush(iosiminfo)
end subroutine

subroutine colvar_printcv_heavyrmsd(tpgs)
  integer,intent(in),dimension(:) :: tpgs
  integer :: igroup,tpngroup,tpnres,ires,index,i,iatom,tpnheavy,tpheavys(numheavy),tpheavyindex(0:numheavy)
  character*150:: tpformat

  tpngroup = size(tpgs)/3; tpnres=0; tpnheavy=0; tpheavyindex=0
  do igroup=1,tpngroup
     do ires = tpgs(3*igroup-2), tpgs(3*igroup-1)
        do iatom=resindex(ires)+1, resindex(ires+1)
           if ( (atomname(iatom)(1:1) == "N") .or. (atomname(iatom)(1:1) == "C") .or. &
                (atomname(iatom)(1:1) == "O") .or. (atomname(iatom)(1:1) == "S") .or. &
                (atomname(iatom)(1:1) == "P") )  then
              tpnheavy = tpnheavy + 1; tpheavys(tpnheavy) = iatom
           end if
        end do
     end do
     tpheavyindex(igroup) = tpnheavy
  end do

  write(iosiminfo,*)
  write(iosiminfo,"(a)") "&colvar"
  write(iosiminfo,"(a)") "  cv_type = 'MULTI_RMSD'"
  write(iosiminfo,"(a,i3,a,i4)") "  cv_ni= ",tpnheavy+tpngroup,", cv_nr=",(tpnheavy*3)

  tpformat(1:3)="(a,"; index=4
  do igroup=1,tpngroup
     write(tpformat(index+(igroup-1)*17:index+(igroup-1)*17+16),"(i3,a14)") &
          tpheavyindex(igroup)-tpheavyindex(igroup-1),"(i4,a2),'0, ',"
  end do
  tpformat(index+tpngroup*17-1:index+tpngroup*17-1) = ")"
  write(iosiminfo,fmt=tpformat) "  cv_i = ", ( tpheavys(i),", ", i=1,tpnheavy )

  write(iosiminfo,"(a)") "  cv_r = "
  do igroup=1,tpngroup
      write(iosiminfo,"(6(f10.3,a2))") &
      ( (nowcoor(index,tpheavys(i)),", ",index=1,3),i=tpheavyindex(igroup-1)+1, tpheavyindex(igroup) )
  end do
  write(iosiminfo,"(a)") "/"
  flush(iosiminfo)
end subroutine

subroutine colvar_printcv_backbonenhb(istart,istop,iskip)
  integer,intent(in) :: istart,istop,iskip
  integer :: nhb,iatom,jatom,ires,jres,io,ih,io2,ih2,hbonds(4,numres*4),ibond

  nhb=0
  do ires=istart,istop
!     if ( (resname(ires)=="ACE ") .or. (resname(ires)=="NME ") .or. (resname(ires)=="NH2 ") ) cycle 
     do iatom=resindex(ires)+1,resindex(ires+1)
        io=0; ih=0
        if ( atomname(iatom)(1:4) == "O   " ) then
           io = iatom
        else if ( atomname(iatom)(1:4) == "H   " ) then
           ih = iatom
        end if
        if ( (io == 0) .and. (ih == 0) ) cycle
     
        do jres = ires + iskip+1, istop
           do jatom = resindex(jres)+1,resindex(jres+1)
              io2=0; ih2=0
              if ( atomname(jatom)(1:4) == "O   " ) then
                 io2 = jatom
              else if ( atomname(jatom)(1:4) == "H   " ) then
                 ih2 = jatom
              end if
              if ( (io2 == 0) .and. (ih2 == 0) ) cycle
              if ( (io>0) .and. (ih2>0) ) then
                 nhb = nhb + 1
                 hbonds(1:4,nhb) = (/io,ih2,ires,jres/)
              else if ( (ih>0) .and. (io2>0) ) then
                 nhb = nhb + 1
                 hbonds(1:4,nhb) = (/ih,io2,ires,jres/)
              end if
              write(iosiminfo,"(a,i3,a1,3x,a4,i3,a1,a4,a3,a4,i3,a1,a4)")  &
                 "hbond ",nhb,":",resname(hbonds(3,nhb)),hbonds(3,nhb),":", &
                 atomname(hbonds(1,nhb))," - ",resname(hbonds(4,nhb)),hbonds(4,nhb),":",atomname(hbonds(2,nhb))

              if ( nhb >= 4*numres ) call errormsg("hbonds array is too small!",int1=4*numres)
           end do
        end do
     end do
  end do

  write(iosiminfo,*)
  write(iosiminfo,"(a)") "&colvar"
  write(iosiminfo,"(a)") "  cv_type = 'N_OF_BONDS'"
  write(iosiminfo,"(a,i3,a,i2)") "  cv_ni= ",nhb*2,", cv_nr= ",1
  write(iosiminfo,"(a,100(i4,a2))") "  cv_i = ",((hbonds(iatom,ibond),", ",iatom=1,2),ibond=1,nhb)
  write(iosiminfo,"(a,f4.1))") "  cv_r = ",3.0d0
  write(iosiminfo,"(a)") "/"

  return
  write(iosiminfo,"(a)") "every backbone h-bond distance"
  do ibond = 1, nhb
     write(iosiminfo,*)
     write(iosiminfo,"(a)") "&colvar"
     write(iosiminfo,"(a)") "  cv_type = 'DISTANCE'"
     write(iosiminfo,"(a,i3,a,i4)") "  cv_ni= ",2
     write(iosiminfo,"(a,i4,a2,i4)") "  cv_i = ",hbonds(1,ibond),", ",hbonds(2,ibond)
     write(iosiminfo,"(a)") "/"
  end do
end subroutine

subroutine colvar_printcv_chain2chaindis(chainarray1,chainarray2)
  integer,intent(in),dimension(:) :: chainarray1,chainarray2
  integer :: i,ires,igroup,tpngroup1,tpngroup2,tpgroup1(2,numres),tpgroup2(2,numres),tpnres1,tpnres2

  tpngroup1=size(chainarray1,1); tpngroup2=size(chainarray2,1); tpnres1=0; tpnres2=0
  do i=1,tpngroup1
     igroup = chainarray1(i)
     tpgroup1(1,i) = molindex(igroup)+1; tpgroup1(2,i) = molindex(igroup+1)
     tpnres1 = tpnres1 + tpgroup1(2,i) - tpgroup1(1,i) + 1
  end do
  do i=1,tpngroup2
     igroup = chainarray2(i)
     tpgroup2(1,i) = molindex(igroup)+1; tpgroup2(2,i) = molindex(igroup+1)
     tpnres2 = tpnres2 + tpgroup2(2,i) - tpgroup2(1,i) + 1
  end do

  write(iosiminfo,*)
  write(iosiminfo,"(a)") "&colvar"
  write(iosiminfo,"(a)") "  cv_type = 'COM_DISTANCE'"
  write(iosiminfo,"(a,i6)") "  cv_ni= ",tpnres1+tpnres2+1
  write(iosiminfo,"(a,100(i4,a2))") "  cv_i = ",((calpha(ires),",",ires=tpgroup1(1,igroup),tpgroup1(2,igroup)), &
                             igroup=1,tpngroup1)
  write(iosiminfo,"(a,100(i4,a2))") "      0, ",((calpha(ires),",",ires=tpgroup2(1,igroup),tpgroup2(2,igroup)), &
                             igroup=1,tpngroup2)
  write(iosiminfo,"(a)") "/"
  flush(iosiminfo)
end subroutine

subroutine colvar_calandrec(tpcvvec,tpmoviecv,ifonlycv,ifnorecord,tpaddvar)
  use nfe_colvar_mod, only : colvar_value
  real*8,intent(out) :: tpcvvec(ncv)
  integer,intent(in),optional :: tpmoviecv
  logical,intent(in),optional :: ifonlycv,ifnorecord
  real*8,intent(in),optional :: tpaddvar(:)
  integer,save :: iomovie,ioprocinfo,tptotntraj,iotraj
  integer,dimension(:),allocatable,save :: iolevel,tplastnkelvin,tptrajindex,ioleveltraj
  real*8,dimension(:,:,:),allocatable,save :: tpallcoor
  real*8,save :: moviecvvalue,lasttptime
  integer :: i,j,k,icv,ierr,iatom,ires,ikelvin
  integer,dimension(0:procnum-1) :: tpcountarray,tpposarray,tprepnkelvin
  real*8 :: temp,tpcoor(3,numboxatom),tpreppot(0:procnum-1),tprepcv(1:ncv,0:procnum-1)
  character*100 :: tpfilename
  logical :: ifrecord

  do icv = 1, nfencv
     tpcvvec(icv) = colvar_value(nfecvs(icv), nowcoor)
  end do
  do icv = nfencv+1, ncv
     call colvar_fsainfo(icv-nfencv,nowcoor,tpcvvec(icv))
  end do
  if ( (present(ifnorecord)) .and. (ifnorecord .eqv. .true.) ) return

  if ( mod(nowstep,samplestep) == 0 ) then

! save cvs on each proc0
     if ( ioprocinfo == 0 ) then
        call getfreeunit(ioprocinfo); write(tpfilename,"(i4)") procid; tpfilename=adjustl(tpfilename)
        tpfilename="procinfo/procinfo_"//trim(tpfilename)//".txt"
        open(file=tpfilename, unit=ioprocinfo, action="write")
     end if
     if ( present(tpaddvar) ) then
        write(ioprocinfo,"(f10.3,i3,f11.3,100f8.3)") nowtime,nowikelvin,epot,tpcvvec(1:ncv),tpaddvar(:)
     else
        write(ioprocinfo,"(f10.3,i3,f11.3,100f8.3)") nowtime,nowikelvin,epot,tpcvvec(1:ncv)
     end if
     flush(ioprocinfo)

goto 1001
! save cvs on proc0
     if ( procid == 0 ) then
        tpcountarray(0:procnum-1)=1
        do i=0,procnum-1; tpposarray(i)=i; end do
     end if
     call mpi_gatherv(nowikelvin, 1, mpi_integer, &
          tprepnkelvin, tpcountarray, tpposarray,mpi_integer,0,mpi_comm_world, ierr)
     call mpi_gatherv(epot, 1, mpi_double_precision, &
          tpreppot, tpcountarray, tpposarray,mpi_double_precision,0,mpi_comm_world, ierr)

     if ( procid == 0 ) then
        tpcountarray(0:procnum-1)=ncv
        do i=0,procnum-1; tpposarray(i)=i*ncv; end do
     end if
     call mpi_gatherv(tpcvvec, ncv, mpi_double_precision, &
          tprepcv, tpcountarray, tpposarray, mpi_double_precision,0,mpi_comm_world, ierr)

     if ( procid == 0 ) then
        if ( .not. allocated(iolevel) ) then
           allocate(iolevel(0:procnum-1),tptrajindex(0:procnum-1),tplastnkelvin(0:procnum-1)); 
           iolevel=0; tptrajindex=0; tplastnkelvin=0
           do ikelvin=0,procnum-1
              call getfreeunit(iolevel(ikelvin)); write(tpfilename,"(i4)") ikelvin+1
              tpfilename=adjustl(tpfilename)
              tpfilename="procinfo/level_"//trim(tpfilename)//".par"
              open(file=tpfilename, unit=iolevel(ikelvin), action="write")
              tptrajindex(ikelvin) = ikelvin + 1
           end do
           tptotntraj = procnum; tplastnkelvin = tprepnkelvin
        end if

        ifrecord=.false.
        if ( nowtime > lasttptime ) ifrecord=.true.
        do i=0,procnum-1
           ikelvin = tprepnkelvin(i)
           if ( ikelvin /= tplastnkelvin(i) ) then
              tptotntraj = tptotntraj + 1; tptrajindex(i) = tptotntraj
           end if
           if ( ifrecord .eqv. .true. ) then
              write(iolevel(ikelvin),"(i8,i3,f8.2,f11.3,100f8.3)") tptrajindex(i),i,nowtime,tpreppot(i),tprepcv(1:ncv,i)
              flush(iolevel(ikelvin))
           end if
        end do
        tplastnkelvin = tprepnkelvin; lasttptime=nowtime
     end if

1001 continue

  end if
  
  if ( (present(ifonlycv)) .and. (ifonlycv .eqv. .true.) ) return

  if ( (mod(nowstep,savetrajstep) == 0) .or. (nowstep==1) ) then
! save traj on each proc
     call moldata_write_traj()

goto 1002
! save traj on proc0
     if ( procid == 0 ) then
        if ( .not. allocated(ioleveltraj) ) then
           allocate(ioleveltraj(0:procnum-1)); ioleveltraj=0
           allocate(tpallcoor(3,numboxatom,0:procnum-1)); tpallcoor=0d0
           do ikelvin=0,procnum-1
              call getfreeunit(ioleveltraj(ikelvin)); write(tpfilename,"(i4)") ikelvin+1
              tpfilename=adjustl(tpfilename)
              tpfilename="procinfo/level_"//trim(tpfilename)//".mdcrd"
              open(file=tpfilename, unit=ioleveltraj(ikelvin), action="write")
              write(ioleveltraj(ikelvin),"(a,i10)") "amber trajectory ",numres
           end do
        end if
        tpcountarray(0:procnum-1)=3*numboxatom
        do i=0,procnum-1; tpposarray(i)=i*3*numboxatom; end do
     end if
     call mpi_gatherv(tpcoor, 3*numboxatom, mpi_double_precision, &
        tpallcoor, tpcountarray, tpposarray, mpi_double_precision,0,mpi_comm_world, ierr)
     if ( (procid == 0) .and. (ifrecord .eqv. .true.) ) then
        do i=0,procnum-1
           ikelvin = tprepnkelvin(i)
!           write(ioleveltraj(ikelvin),"(10f8.3)") tpallcoor(:,:,i)
           write(ioleveltraj(ikelvin),"(10f8.3)") ((tpallcoor(k,calpha(j),i),k=1,3),j=1,numres)
           flush(ioleveltraj(ikelvin))
        end do
     end if

1002 continue

  end if

  if ( (present(tpmoviecv)) .and. (tpmoviecv/=0) .and. (abs(tpmoviecv) <= ncv) ) then
  if ( iomovie == 0 ) then
     call getfreeunit(iomovie); write(tpfilename,"(i4)") procid; tpfilename=adjustl(tpfilename)
     tpfilename="procinfo/movie_"//trim(tpfilename)//".pdb"
     open(file=tpfilename, unit=iomovie, action="write")
  end if
  if ( ((mod(nowstep,samplestep) == 0)) .or. (nowstep==1) ) then
     ifrecord = .false.
     if ( moviecvvalue == 0.0d0 ) ifrecord = .true.
     if ( (tpmoviecv < 0 ) .and. (tpcvvec(abs(tpmoviecv)) < moviecvvalue) ) ifrecord = .true.
     if ( (tpmoviecv > 0 ) .and. (tpcvvec(abs(tpmoviecv)) > moviecvvalue) ) ifrecord = .true.
     if ( ifrecord .eqv. .true. ) then
        moviecvvalue = tpcvvec(abs(tpmoviecv))
        call gpu_download_crd(tpcoor); call moldata_moveto_comorcenterbox(tpcoor)
        do ires=1,numres
           do iatom=resindex(ires)+1,resindex(ires+1)
              write(iomovie,"(a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3)")   &
                      "ATOM  ",iatom,atomname(iatom)(1:4),resname(ires)(1:3)," ",ires,tpcoor(1:3,iatom)
           end do
           if ( (ires<numres) .and. (resmolindex(ires) /= resmolindex(ires+1)) ) write(iomovie,"(a,i4)") "TER ",resmolindex(ires)
        end do
        if ( ncv <= 10 ) then
           write(iomovie,"(a4,f10.3,10f7.2)") "END ",nowtime,tpcvvec(1:ncv); flush(iomovie)
        else
           write(iomovie,"(a4,f10.3,10f7.2)") "END ",nowtime,tpcvvec(1:10); flush(iomovie)
        end if
     end if
  end if
  end if
end subroutine

subroutine colvar_fsainfo(icv,tpcoor,tpvalue,tpgrad, tpfactor)
  integer,intent(in) :: icv
  real*8,intent(in) :: tpcoor(3,numatom)
  real*8,intent(out) :: tpvalue
  real*8,intent(inout),optional :: tpgrad(3,numatom)
  real*8,intent(in),optional :: tpfactor

  if ( fsacvs(icv)%type == CV_COM_POS ) then
     call colvar_com_pos(icv, tpcoor, tpvalue, tpgrad, tpfactor)
  else if ( fsacvs(icv)%type == CV_COM_DIS ) then
     call colvar_com_dis(icv, tpcoor, tpvalue, tpgrad, tpfactor)
  else if ( fsacvs(icv)%type == CV_COM_POSAXIS ) then
     call colvar_com_posaxis(icv, tpcoor, tpvalue, tpgrad, tpfactor)
  else if ( fsacvs(icv)%type == CV_PCA ) then
     call colvar_pca(icv, tpcoor, tpvalue, tpgrad, tpfactor)
  else if ( fsacvs(icv)%type == CV_RMSD ) then
     call colvar_rmsd(icv, tpcoor, tpvalue, tpgrad, tpfactor)
  else if ( fsacvs(icv)%type == CV_DISSUM ) then
     call colvar_dissum(icv, tpcoor, tpvalue, tpgrad, tpfactor)
  else
     call errormsg("wrong fsacolvar type ! ", int1=fsacvs(icv)%type)
  end if
end subroutine

subroutine colvar_dissum(icv,tpcoor,tpvalue,tpgrad,tpfactor)
  integer,intent(in) :: icv
  real*8,intent(in) :: tpcoor(3,numatom)
  real*8,intent(out) :: tpvalue
  real*8,intent(inout),optional :: tpgrad(3,numatom)
  real*8,intent(in),optional :: tpfactor
  integer :: tpndim,iatom,jatom,i
  real*8 :: tpdis,tpvec(3)

  tpvalue = 0d0
  do i=1,fsacvs(icv)%ni,2
     iatom = fsacvs(icv)%i(i); jatom = fsacvs(icv)%i(i+1)
     tpvec(1:3) = tpcoor(1:3,iatom) - tpcoor(1:3,jatom)
     tpdis = sqrt(dot_product(tpvec,tpvec))
     tpvalue = tpvalue + tpdis
     if ( present(tpgrad) ) then
        if ( present(tpfactor) ) tpvec=tpfactor*tpvec
        tpgrad(1:3,iatom) =  tpgrad(1:3,iatom) + tpvec(1:3)/tpdis
        tpgrad(1:3,jatom) =  tpgrad(1:3,jatom) - tpvec(1:3)/tpdis
     end if 
  end do
end subroutine

subroutine colvar_rmsd_updaterefcoor(icv,tpcoor)
  integer,intent(in) :: icv
  real*8,intent(in) :: tpcoor(3,numatom)
  integer :: i,iatom
  real*8 :: tpvec(3)
  if ( (icv <= 0) .or. (fsacvs(icv)%type /= CV_RMSD) ) then
     call errormsg("error cv index in colvar_rmsd_updaterefcoor ",int1=icv)
  end if
  do i=1,fsacvs(icv)%ni
     iatom=fsacvs(icv)%i(i)
     fsacvs(icv)%r(3*i-2:3*i) = tpcoor(1:3,iatom)
  end do
  tpvec=0d0
  do i=1,fsacvs(icv)%ni
     tpvec(1:3) = tpvec(1:3) + fsacvs(icv)%r(3*i-2:3*i)
  end do
  tpvec(1:3)=tpvec(1:3)/dble(fsacvs(icv)%ni)
  do i=1,fsacvs(icv)%ni
     fsacvs(icv)%r(3*i-2:3*i) = fsacvs(icv)%r(3*i-2:3*i) - tpvec(1:3)
  end do
end subroutine

subroutine colvar_rmsd(icv,tpcoor,tpvalue,tpgrad,tpfactor)
  integer,intent(in) :: icv
  real*8,intent(in) :: tpcoor(3,numatom)
  real*8,intent(out) :: tpvalue
  real*8,intent(inout),optional :: tpgrad(3,numatom)
  real*8,intent(in),optional :: tpfactor
  integer :: tpndim,iatom,i
  real*8,dimension(fsacvs(icv)%ni*3) :: coor1,coor2,tptpgrad

  do i=1,fsacvs(icv)%ni
     iatom = fsacvs(icv)%i(i); coor1(3*i-2:3*i)=tpcoor(1:3,iatom)
  end do
  coor2=fsacvs(icv)%r(1:fsacvs(icv)%nr)
  if ( .not. present(tpgrad) ) then
     call moldata_lRMSD(fsacvs(icv)%nr, coor1, coor2, tpvalue)
  else
     call moldata_lRMSD(fsacvs(icv)%nr, coor1, coor2, tpvalue, tpgrad=tptpgrad)
     if ( present(tpfactor) ) tptpgrad=tpfactor*tptpgrad
     do i=1,fsacvs(icv)%ni
        iatom = fsacvs(icv)%i(i); tpgrad(1:3,iatom)=tpgrad(1:3,iatom)+tptpgrad(3*i-2:3*i)
     end do
  end if
end subroutine

subroutine colvar_com_pos(icv,tpcoor,tpvalue,tpgrad,tpfactor)
  integer,intent(in) :: icv
  real*8,intent(in) :: tpcoor(3,numatom)
  real*8,intent(out) :: tpvalue
  real*8,intent(inout),optional :: tpgrad(3,numatom)
  real*8,intent(in),optional :: tpfactor

  integer :: i,iatom,iaxis
  real*8 :: tppos(3),tpvec(3),tpmass,tpcoeff

  tppos(1:3)=0.0d0; tpmass=0d0
  do i=1,fsacvs(icv)%ni
     iatom = fsacvs(icv)%i(i)
     tppos(1:3)=tppos(1:3)+tpcoor(1:3,iatom)*atommass(iatom)
     tpmass = tpmass + atommass(iatom)
  end do
  tppos(1:3)=tppos(1:3)/tpmass
  tpvec(1:3) = tppos(1:3) - fsacvs(icv)%r(1:3)

  tpvalue = 0d0
  do iaxis=1,3
     if ( fsacvs(icv)%r(iaxis) < -999.0d0 ) cycle
     tpvalue = tpvalue + tpvec(iaxis)**2
  end do
  tpvalue = sqrt(tpvalue)

  if ( .not. present(tpgrad) ) return

  tpcoeff=1.0d0; if ( present(tpfactor) ) tpcoeff=tpfactor
  tpvec(1:3) = tpcoeff*tpvec(1:3)/(tpvalue*tpmass)
  do i=1,fsacvs(icv)%ni
     iatom = fsacvs(icv)%i(i)
     do iaxis=1,3
        if ( fsacvs(icv)%r(iaxis) < -999.0d0 ) cycle
        tpgrad(iaxis,iatom) = tpgrad(iaxis,iatom) + tpvec(iaxis)*atommass(iatom)
     end do
  end do
end subroutine

subroutine colvar_com_posaxis(icv,tpcoor,tpvalue,tpgrad,tpfactor)
  integer,intent(in) :: icv
  real*8,intent(in) :: tpcoor(3,numatom)
  real*8,intent(out) :: tpvalue
  real*8,intent(inout),optional :: tpgrad(3,numatom)
  real*8,intent(in),optional :: tpfactor

  integer :: i,iatom,iaxis
  real*8 :: tppos(3),tpvec(3),tpmass,tpcoeff

  iaxis=0
  do i=1,3
     if ( fsacvs(icv)%r(i) > -999.0d0 ) then
        iaxis=i; exit;
     end if
  end do
  if ( (iaxis < 1) .or. (iaxis > 3) ) call errormsg("wrong axis in posaxis ",int1=iaxis)

  tppos(iaxis)=0.0d0; tpmass=0d0
  do i=1,fsacvs(icv)%ni
     iatom = fsacvs(icv)%i(i)
     tppos(iaxis)=tppos(iaxis)+tpcoor(iaxis,iatom)*atommass(iatom)
     tpmass = tpmass + atommass(iatom)
  end do
  tppos(iaxis)=tppos(iaxis)/tpmass
  tpvalue = tppos(iaxis)

  if ( .not. present(tpgrad) ) return

  tpcoeff=1.0d0; if ( present(tpfactor) ) tpcoeff=tpfactor
  tpcoeff = tpcoeff/tpmass
  do i=1,fsacvs(icv)%ni
     iatom = fsacvs(icv)%i(i)
     tpgrad(iaxis,iatom) = tpgrad(iaxis,iatom) + tpcoeff*atommass(iatom)
  end do
end subroutine

subroutine colvar_com_dis(icv,tpcoor,tpvalue,tpgrad,tpfactor)
  integer,intent(in) :: icv
  real*8,intent(in) :: tpcoor(3,numatom)
  real*8,intent(out) :: tpvalue
  real*8,intent(inout),optional :: tpgrad(3,numatom)
  real*8,intent(in),optional :: tpfactor

  integer :: i,igroup,iatom,iaxis
  real*8 :: tppos1(3),tppos2(3),tpvec(3),tpmass1,tpmass2,tpcoeff

  if ( fsacvs(icv)%ngroup /= 2 ) call errormsg("com_dis error ",int1=fsacvs(icv)%ngroup)
  tppos1(1:3)=0.0d0; tppos2(1:3)=0d0; tpmass1=0d0; tpmass2=0d0; tpcoeff=1.0d0
  if ( present(tpfactor) ) tpcoeff=tpfactor
  do igroup=1,fsacvs(icv)%ngroup
     do i=fsacvs(icv)%groups(1,igroup),fsacvs(icv)%groups(2,igroup)
        iatom = fsacvs(icv)%i(i)
        if ( igroup == 1 ) then
           tppos1(1:3) = tppos1(1:3) + tpcoor(1:3,iatom)*atommass(iatom)
           tpmass1 = tpmass1 + atommass(iatom)
        else if ( igroup == 2 ) then
           tppos2(1:3) = tppos2(1:3) + tpcoor(1:3,iatom)*atommass(iatom)
           tpmass2 = tpmass2 + atommass(iatom)
        end if
     end do
  end do
  tppos1(1:3)=tppos1(1:3)/tpmass1; tppos2(1:3)=tppos2(1:3)/tpmass2
  tpvec(1:3) = tppos1(1:3) - tppos2(1:3)
  tpvalue = 0d0
  do iaxis=1,3
     if ( fsacvs(icv)%r(iaxis) < 0d0 ) cycle
     tpvalue = tpvalue + tpvec(iaxis)**2
  end do
  tpvalue = sqrt(tpvalue)

  if ( .not. present(tpgrad) ) return

  tpvec(1:3) = tpcoeff*tpvec(1:3)/tpvalue
  do igroup=1,fsacvs(icv)%ngroup
     do i=fsacvs(icv)%groups(1,igroup),fsacvs(icv)%groups(2,igroup)
        iatom = fsacvs(icv)%i(i)
        do iaxis=1,3
           if ( fsacvs(icv)%r(iaxis) < 0d0 ) cycle
           if ( igroup == 1 ) then
              tpgrad(iaxis,iatom) = tpgrad(iaxis,iatom) + tpvec(iaxis)*atommass(iatom)/tpmass1
           else if ( igroup == 2 ) then
              tpgrad(iaxis,iatom) = tpgrad(iaxis,iatom) - tpvec(iaxis)*atommass(iatom)/tpmass2
           end if
        end do
     end do
  end do
end subroutine

subroutine colvar_pca(tpcvindex,tptpcoor,tpcv,tpgrad,tpfactor)
  integer,intent(in) :: tpcvindex
  real*8,intent(in) :: tptpcoor(3,numatom)
  real*8,intent(out) :: tpcv
  real*8,intent(inout),optional :: tpgrad(3,numatom)
  real*8,intent(in),optional :: tpfactor
  integer,save :: pcanatom,pcanpc
  integer,allocatable,dimension(:),save :: pcaatoms,pcaindex,pc2cv
  real*8,allocatable,dimension(:),save :: tprefcoor,tppcavec,tpcoor
  real*8,allocatable,dimension(:,:),save :: pcavec,pcavecsum
  real*8,save :: tprotmat(3,3),tpdis
  logical,save :: ifpolar
  integer :: ipc,i,j,icv,iofile,iatom,tppcindex
  real*8 :: tpcoeff,temp
  character*100 :: pcaeigenfile
  logical :: iffind

  if ( .not. allocated(pcavecsum) ) then

     ipc=0; allocate(pcaindex(cvnatom*3),pc2cv(ncv-nfencv)); pcaindex=0; pc2cv=0
     do icv=1,ncv-nfencv
        if ( fsacvs(icv)%type == CV_PCA ) then
           ipc = ipc + 1; pcaindex(ipc)=fsacvs(icv)%cv_index; pc2cv(ipc)=icv
           if ( ipc == 1 ) then
              pcaeigenfile=fsacvs(icv)%cv_file; pcanatom=fsacvs(icv)%ni
              allocate(pcaatoms(pcanatom)); pcaatoms(1:pcanatom)=fsacvs(icv)%i(1:pcanatom)
           end if
        end if
     end do
     if ( pcaindex(1) < 0 ) then
        if ( ipc /= 2 ) call errormsg("wrong number of pca collective variable!")
        pcanpc = abs(pcaindex(1)); ifpolar=.true.
        do i=1,pcanpc; pcaindex(i)=i; end do
     else
        pcanpc = ipc
        do ipc=1,pcanpc
           do i=ipc+1,pcanpc
              if ( pcaindex(ipc) == pcaindex(i) ) call errormsg("pcaindex error ",int1=ipc,int2=i,int3=pcaindex(i))
           end do
        end do
     end if
     allocate(tprefcoor(3*pcanatom),tpcoor(3*pcanatom)); tprefcoor=0d0; tpcoor=0d0
     allocate(pcavec(3*pcanatom,pcanpc),tppcavec(3*pcanatom)); pcavec=0d0; tppcavec=0d0

     call getfreeunit(iofile); open(unit=iofile,file=pcaeigenfile,action="read")
     read(iofile,"(10x,10000f8.3)") tppcavec(1:3*pcanatom)
     tprefcoor(1:3*pcanatom) = tppcavec(1:3*pcanatom)
!     do i=1, pcanatom; tprefcoor(3*i-2:3*i) = nowcoor(1:3, pcaatoms(i)); end do
     read(iofile,*); read(iofile,*)
     i=0; ipc=1
     do
        read(iofile,"(4x,10000f8.3)") tppcavec(1:3*pcanatom)
        i = i + 1
        if ( i == pcaindex(ipc) ) then
           pcavec(1:3*pcanatom,ipc)=tppcavec(1:3*pcanatom)
           tpcoeff = sqrt(dot_product(pcavec(1:3*pcanatom,ipc),pcavec(1:3*pcanatom,ipc)))
           if ( abs(1.0d0-tpcoeff) > 0.01d0 ) call errormsg("pca axis is not normalized ",int1=ipc,real1=tpcoeff)
           if ( ipc == pcanpc ) exit
           ipc = ipc + 1
        end if
     end do
     close(iofile)

     allocate(pcavecsum(3,pcanpc)); pcavecsum=0d0
     do ipc=1,pcanpc
        do i=1,pcanatom
           pcavecsum(1:3,ipc) = pcavecsum(1:3,ipc) + pcavec(3*i-2:3*i,ipc)
        end do
     end do
     pcavecsum = pcavecsum/dble(pcanatom)

     if ( procid == 0 ) then
        write(iosiminfo,*)
        write(iosiminfo,"(a)") "use pca collective variables"
        write(iosiminfo,"(a)") "pca eigen file ",trim(pcaeigenfile)
        write(iosiminfo,"(a,i6,5x,a,i6)") "pca natom ",pcanatom
        write(iosiminfo,"(a)") "pca atoms "
        write(iosiminfo,"(18i6)") pcaatoms(1:pcanatom)
        write(iosiminfo,"(a,i6)") "pca npc ",pcanpc
        if ( ifpolar .eqv. .false. ) then
           write(iosiminfo,"(a,100i4)") "pca indices: ",pcaindex(1:pcanpc)
        else
           write(iosiminfo,"(a)") "use polar coordinates in pca space"
        end if
        flush(iosiminfo)
     end if
  end if

  iffind=.false.; tppcindex=0
  do
     tppcindex = tppcindex + 1
     if ( pc2cv(tppcindex) == tpcvindex ) then; iffind=.true.; exit; end if
     if ( ifpolar .eqv. .false. ) then
        if ( tppcindex == pcanpc ) exit
     else
        if ( tppcindex == 2 ) exit
     end if
  end do
  if ( iffind .eqv. .false. ) call errormsg("wrong tpcvindex in moldata_colvar_pca",int1=tpcvindex+nfencv)
         
  if ( tppcindex == 1 ) then
     do i=1, pcanatom
        tpcoor(3*i-2:3*i) = tptpcoor(1:3, pcaatoms(i))
     end do
     call moldata_lRMSD(pcanatom*3,tpcoor,tprefcoor,tpcoeff,ifrotatecoor1=.true.,rotatemat=tprotmat)
     if ( ifpolar .eqv. .true. ) then
        tppcavec=0d0
        do ipc=1,pcanpc
           tppcavec(ipc)=dot_product(pcavec(1:3*pcanatom,ipc),tpcoor(1:3*pcanatom))
        end do
        tpdis = sqrt(dot_product(tppcavec(1:pcanpc),tppcavec(1:pcanpc)))
     end if
  end if

  if ( ifpolar .eqv. .false. ) then
     tpcv = dot_product(pcavec(1:3*pcanatom,tppcindex),tpcoor(1:3*pcanatom))
  else
     if ( tppcindex == 1 ) then
        tpcv = tpdis
     else if ( tppcindex == 2 ) then
        tpcv = acos(tppcavec(1)/tpdis)
     else
        call errormsg("wrong tppcindex ",int1=tppcindex,int2=tpcvindex)
     end if
  end if

  if ( .not. present(tpgrad) ) return

  tpcoeff=1.0d0; if ( present(tpfactor) ) tpcoeff=tpfactor

  if ( ifpolar .eqv. .false. ) then

     do i=1,pcanatom
        iatom = pcaatoms(i)
        do j=1,3
           tpgrad(1:3,iatom) = tpgrad(1:3,iatom) + tpcoeff*  &
              ( pcavec(3*i-3+j,tppcindex)*tprotmat(j,1:3) - pcavecsum(j,tppcindex)*tprotmat(j,1:3) )
        end do
     end do

  else

     if ( tppcindex == 1 ) then
        tpcoeff = tpcoeff/tpdis
        do ipc=1,pcanpc
           do i=1,pcanatom
              iatom = pcaatoms(i)
              do j=1,3
                 tpgrad(1:3,iatom) = tpgrad(1:3,iatom) + tpcoeff*tppcavec(ipc)*  &
                     ( pcavec(3*i-3+j,ipc)*tprotmat(j,1:3) - pcavecsum(j,ipc)*tprotmat(j,1:3) )
              end do
           end do
        end do
     else if ( tppcindex == 2 ) then
        temp = tppcavec(1)/tpdis; temp = -1d0/sqrt(1-temp**2)
        tpcoeff = tpcoeff*temp
        do ipc=1,pcanpc
           do i=1,pcanatom
              temp = tpcoeff/tpdis
              do j=1,3
                 tpgrad(1:3,iatom) = tpgrad(1:3,iatom) + temp*  &
                    ( pcavec(3*i-3+j,1)*tprotmat(j,1:3) - pcavecsum(j,1)*tprotmat(j,1:3) )
              end do
              temp = tpcoeff*tppcavec(1)*(-1d0)/(tpdis**3)
              do j=1,3
                 tpgrad(1:3,iatom) = tpgrad(1:3,iatom) + temp*tppcavec(ipc)* &
                     ( pcavec(3*i-3+j,ipc)*tprotmat(j,1:3) - pcavecsum(j,ipc)*tprotmat(j,1:3) )
              end do
           end do
        end do
     end if
  end if
end subroutine

end module
