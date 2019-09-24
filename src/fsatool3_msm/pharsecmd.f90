module getcmd
    use mod_global, only: inputfile, resultdir
    implicit none
    integer, parameter :: max_str_len = 256
    character(max_str_len) :: argname, resultfile
    integer, save :: iarg = 0
    public :: get_program
contains
    subroutine get_program(msm)
        logical  :: msm
        msm = .false.
        !if(iargc() == 0) return;
        iarg = iarg + 1;
        call get_command_argument(iarg, argname)
        select case (trim(argname))
            case("cluster")
                call run_cluster()
            case("lumping")
                call run_lumping()
            case("tpt")
                call run_tpt()
            case("check")
                call run_check()
            case("tram")
                call run_tram()
                msm = .true.
            case ("transform")
                call run_transform()
            case("msm")
                call run_msm()
                msm = .true.
            case default
                call run_msm()
                msm = .true.
        end select
    end subroutine get_program

    subroutine get_nextcmd()
        integer :: narg
        narg = iargc()
        do
            if(iarg >= narg) exit
            iarg = iarg + 1
            call get_command_argument(iarg, argname)
            if (trim(argname) == "-i") then
                iarg = iarg + 1
                call get_command_argument(iarg, argname)
                inputfile = argname
            elseif  (trim(argname) == "-o") then
                iarg = iarg + 1
                call get_command_argument(iarg, argname)
                resultfile = argname
            else
                continue
            endif
        enddo
    end subroutine

    subroutine run_cluster()
        use cluster, only: mod_cluster
        inputfile = "cluster.in"
        resultfile = "cluster.out"
        call get_nextcmd()
        call check_file_exist()
        call mod_cluster(inputfile, resultfile)
    end subroutine

    subroutine run_lumping()
        use coarse_grain, only: mod_lumping
        inputfile = "lumping.in"
        resultfile = "lumping.out"
        call get_nextcmd()
        call mod_lumping(inputfile, resultfile)
    end subroutine

    subroutine run_tpt()
        use tpt, only: mod_tpt
        inputfile = "tpt.in"
        resultfile = "tpt.out"
        call get_nextcmd()
        call check_file_exist()
        call mod_tpt(inputfile, resultfile)
    end subroutine

    subroutine run_check()
        use markov, only: markov_check
        inputfile = "check.in"
        resultfile = "check.out"
        call get_nextcmd()
        call check_file_exist()
        call markov_check(inputfile, resultfile)
    end subroutine

    subroutine run_tram()
        use tram, only: mod_tram
        inputfile = "tram.in"
        resultfile = "tram.out"
        call get_nextcmd()
        call check_file_exist()
        call mod_tram(inputfile, resultfile)
    end subroutine

    subroutine run_transform()
        use transform, only: mod_transform
        inputfile = "transform.in"
        resultfile = "transform.out"
        call get_nextcmd()
        call check_file_exist()
        call mod_transform(inputfile, resultfile)
    end subroutine

    subroutine run_msm()
        inputfile = "msm.in"
        resultdir = "info/"
        resultfile = ""
        call get_nextcmd()
        call check_file_exist()
        if (resultfile /= "") resultdir = resultfile
    end subroutine

    subroutine check_file_exist()
        logical :: file_exist
        inquire(file=inputfile, exist=file_exist)
        if( .not. file_exist) stop "input file not exist"
    end subroutine

end module getcmd
