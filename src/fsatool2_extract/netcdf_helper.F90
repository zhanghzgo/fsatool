module NETCDF_FUNC
    use netcdf
    use AmberNetcdf_mod
    implicit none
    integer :: atom_number
    integer :: nframes
    logical :: has_box
contains

    subroutine init_parameter(filename)
        character(*) :: filename
        character(256) :: title
        integer :: ncid
        integer :: coordVID, velocityVID, timeVID, cellLengthVID, cellAngleVID, TempVID
        if (NC_openRead(filename, ncid)) stop "error open"
        if (NC_setupMdcrd(ncid, title, nframes, atom_number, coordVID, &
            velocityVID, timeVID, cellLengthVID, cellAngleVID, TempVID)) stop "error setup  mdcrd"
        if (cellLengthVID == -1) has_box = .false.
        print*, nframes
    end subroutine

    subroutine read_NCfile(filename, ncid, coorDVID)
        implicit none
        character(*),intent(in) :: filename
        integer, intent(out) :: ncid, coorDVID
        character(256) :: title
        integer :: timeVID, cellLengthVID, cellAngleVID, tempVID
        integer :: velocityVID, frameDID
        if (NC_openRead(filename, ncid)) stop "error open"
        if (NC_setupMdcrd(ncid, title, nframes, atom_number, coordVID, &
            velocityVID, timeVID, cellLengthVID, cellAngleVID, TempVID)) stop "error setup  mdcrd"
    end subroutine
        

    subroutine create_NCfile(filename, ncid, coorDVID)
        implicit none
        character(*),intent(in) :: filename
        integer, intent(out) :: ncid, coorDVID
        integer :: timeVID, velocityVID, frcVID, cellLengthVID, cellAngleVID, tempVID
        integer :: frameDID
        if (NC_create(filename, 'O', .false., atom_number, .true.,  &
                       .false., has_box, .false., .true., &
                       .false., "level file", &
                       ncid, timeVID, coordVID, velocityVID, frcVID,&
                       cellLengthVID, cellAngleVID, TempVID, &
                       frameDID)) stop "file create file"
    end subroutine
    
    subroutine write_NCcoor(ncid, coordVID, crd, mdcrd_frame)
        integer, intent(in) :: ncid, coordVID
        real*8, intent(in):: crd(:, :)
        integer, intent(in) :: mdcrd_frame
        call checkNCerror(nf90_put_var(ncid, coordVID, crd(:, :), &
                                       start=(/1,1,mdcrd_frame/), &
                                       count=(/3, atom_number, 1/)), "write atom coords")
    end subroutine
    
    subroutine write_NCcell(ncid, cellLengthVID, cellAngleVID, length, Angle, mdcrd_frame)
        integer :: ncid, cellLengthVID, cellAngleVID 
        real*8 :: length(3), Angle(3)
        integer :: mdcrd_frame
        call checkNCerror(nf90_put_var(ncid, cellLengthVID, length, &
                              start = (/1, mdcrd_frame/), count = (/3,1/)), "write cell length")
        call checkNCerror(nf90_put_var(ncid, cellAngleVID, angle, &
                              start = (/1, mdcrd_frame/), count = (/3,1/)), "write cell length")
    end subroutine

    subroutine get_NCcoor(ncid, coordVID, crd, frame)
        integer, intent(in) :: ncid, coordVID
        integer, intent(in) :: frame
        real*8, dimension(:, :), intent(out) :: crd
        call checkNCerror(nf90_get_var(ncid, coordVID, crd, start=(/1,1,frame/)), &
                    "get coor variable")
    end subroutine

    subroutine finish_NCwrite(ncid)
        integer :: ncid
        call checkNCerror(nf90_sync(ncid))
    end subroutine

    subroutine close_NCfile(ncid)
        integer :: ncid
        call NC_close(ncid)
    end subroutine
end module
