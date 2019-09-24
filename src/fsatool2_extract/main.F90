program Extract_levelInfo
    use extract_mod
    implicit none
    character(256) :: filename = "extract.in"
    call init_nml(filename)
    call extract_run()
    call finish_extract()
end program
