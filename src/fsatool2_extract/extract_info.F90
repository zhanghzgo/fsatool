module extract_mod
    use NETCDF_FUNC
    implicit none 
    integer, PARAMETER :: MAX_File = 100
    integer, PARAMETER :: MAX_CHAR = 256
    character(MAX_CHAR), DIMENSION(MAX_FILE) :: procFile, trajFile
    character(MAX_CHAR) :: trajType
    integer,allocatable, DIMENSION(:) :: r_ncid,  r_coorDVID
    integer,allocatable, DIMENSION(:) :: w_ncid,  w_coorDVID
    integer :: num_file, ncv, frame_cv
    logical :: onlyDealWithProcFile, removeFirstFrame

contains
    subroutine init_nml(inputFile)
        integer :: iofile, ierr
        CHARACTER(MAX_CHAR) :: inputFile
        logical :: file_exist
        integer :: i
        NAMELIST /trajs/ num_file, procFile, ncv, trajFile, trajType, removeFirstFrame

        trajType = "netcdf"
        procFIle = ""
        trajFile = ""
        onlyDealWithProcFile = .false.
        removeFirstFrame = .false.

        inquire(file=inputFile, exist=file_exist)
        if ( .not. file_exist) stop "input file not exist"

        call getFreeUnit(iofile)
        open(iofile, file = inputFile, action="read")
        read(iofile, nml = trajs, iostat=ierr)
        if (ierr < 0 ) STOP "ERROR read namelist trajs!!"
        close(iofile)
        
        do i = 1, num_file
            inquire(file = procFile(i), exist=file_exist)
            if (.not. file_exist) STOP "procsses file not exist"
        end do

        if (trajFile(1) == "") onlyDealWithProcFile = .true.
        if (.not. onlyDealWithProcFile) then
            do i = 1, num_file
            inquire(file = trajFile(i), exist=file_exist)
            if (.not. file_exist) STOP "Trajectory File not exist" 
            enddo
        endif
        if (.not. onlyDealWithProcFile .and. trajType == "netcdf") call init_parameter(trim(trajFile(1)))
    end subroutine

    subroutine extract_run()
        if (onlyDealWithProcFile) then
            call extract_only_procInfo()
        else
            if (trajType ==  "pdb") then
                ! call extract_pdb()
            else if(trajType == "netcdf") then
                allocate(r_ncid(num_file), r_coorDVID(num_file))
                allocate(w_ncid(num_file), w_coorDVID(num_file))
                call extract_netcdf()
            endif
        endif
    end subroutine

    subroutine extract_netcdf()
        implicit none
        integer :: i, j, ierr, index, count
        integer :: r_levelInfo(num_file)
        character(MAX_CHAR), DIMENSION(num_file) :: level_file
        character(MAX_CHAR) :: temp_char, temp_file
        real*8 :: coor(3, atom_number)
        integer :: interval

        call extract_only_procInfo()

        do i = 1, num_file
            write(temp_char, "(i4)") i-1
            temp_file = "levelinfo_"//trim(adjustl(temp_char))//".txt"
            temp_char = "level_"//trim(adjustl(temp_char))//".nc"
            call getFreeUnit(r_levelInfo(i))
            open(file=temp_file, unit = r_levelInfo(i), action="read")
            call create_NCfile(temp_char, w_ncid(i), w_coorDVID(i))
            call read_NCfile(trajFile(i), r_ncid(i), r_coorDVID(i))
        enddo

        interval = nframes / frame_cv
        print*, "interval:", interval
        do i = 1, num_file
            count = 0
            do 
                read (r_levelInfo(i), *, iostat=ierr) j, index
                if (ierr < 0) exit
                count = count + 1
                if (mod(count, interval) == 0) then
                    j = count / interval
                    if (mod(j, 5000) == 0) print*, "finish", j, "steps"
                    if (removeFirstFrame) then
                        call get_NCcoor(r_ncid(index+1), r_coorDVID(index+1), coor, j+1)
                    else
                        call get_NCcoor(r_ncid(index+1), r_coorDVID(index+1), coor, j)
                    endif
                        call write_NCcoor(w_ncid(i), w_coordVID(i), coor, j)
                endif
            enddo
            call finish_NCwrite(w_ncid(i))
        enddo
    end subroutine

    subroutine extract_only_procInfo()
        implicit none
        integer :: i, itraj
        integer :: r_ProcFile(num_file), w_ProcFile(num_file)
        real*8 :: time, potential, cv(ncv), lastTime(num_file)
        integer :: level, lastLevel(num_file), trajIndex(num_file)
        character(MAX_CHAR) :: temp_char
        integer :: ierr

        do i = 1, num_file
            write(temp_char, "(i4)") i-1
            temp_char = "levelinfo_"//trim(adjustl(temp_char))//".txt"
            call getFreeUnit(r_ProcFile(i))
            open(file=trim(procFile(i)), unit = r_ProcFile(i), action="read")
            call getFreeUnit(w_ProcFile(i))
            open(file = temp_char, unit = w_ProcFile(i), action="write")
        enddo

        lastLevel = -1
        lastTime = -1d0
        itraj = 0
        frame_cv = 0

        outer: do 
            do i = 1, num_file
                read(r_ProcFile(i), *, iostat=ierr) time, level, potential, cv
                if (ierr < 0)  then
                    exit outer
                endif
                level = level + 1
                if (level > num_file) then 
                    print*, i, level, time, ncv, potential
                    STOP "Wrong level index"
                else if (level /= lastlevel(i)) then
                    itraj = itraj + 1
                    lastLevel(i) = level
                    trajindex(i) = itraj
                endif
                if (time /= lastTime(i)) then
                    write(w_ProcFile(level), "(i8, i3, 2x, f10.3, 100f8.3)") trajIndex(i), i-1 , potential, cv
                    frame_cv = frame_cv + 1
                endif
                lastTime(i) = time
            enddo
        enddo outer
        frame_cv = frame_cv / num_file

        do i = 1, num_file
            close(r_ProcFile(i)) 
            close(w_ProcFile(i))
        enddo
    end subroutine

    subroutine finish_extract()
        integer :: i
        if (onlyDealWithProcFile) return
        do i = 1, num_file
            call close_NCfile(r_ncid(i))
            call close_NCfile(w_ncid(i))
        enddo
        deallocate(r_ncid, r_coorDVID)
        deallocate(w_ncid, w_coorDVID)
    end subroutine

    subroutine getFreeUnit(subunit)
        implicit none
        integer, intent(out) :: subunit
        logical :: ifUsed
        subunit = 20
        ifUsed = .true.
        do 
            if (ifUsed .eqv. .false.) exit
            subunit = subunit + 1
            INQUIRE(unit=subunit, opened=ifUsed)
        enddo
    end subroutine
end module
