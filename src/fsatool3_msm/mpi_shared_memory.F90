module mpi_shared_memory
    use mpi
    use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
    implicit none
    integer :: shared_num, shared_id
    integer :: shared_comm_world
contains
    subroutine init_shared_memory(mpi_world)
        integer, intent(in) :: mpi_world
        integer :: ierr
        call mpi_comm_split_type(mpi_world, mpi_comm_type_shared, 0, &
                                mpi_info_null, shared_comm_world, ierr)
        call mpi_comm_rank(shared_comm_world, shared_id, ierr)
        call mpi_comm_size(shared_comm_world, shared_num, ierr)
    end subroutine

    subroutine allocate_shared_memory_for_2dim_array(dim_shape, shared_array, shared_win)
        integer, intent(in) :: dim_shape(2)
        integer, intent(out) :: shared_win
        real*8, pointer, intent(out) :: shared_array(:, :)

        type(C_PTR) :: baseptr
        integer :: ierr, disp_unit
        integer(kind=mpi_address_kind) :: windowsize

        if(shared_id == 0) then
            windowsize =  int(dim_shape(1) * dim_shape(2), MPI_ADDRESS_KIND) * 8_MPI_ADDRESS_KIND
        else
            windowsize = 0_MPI_ADDRESS_KIND
        endif
        disp_unit = 1

        call mpi_win_allocate_shared(windowsize, disp_unit, mpi_info_null, shared_comm_world, baseptr, shared_win, ierr)
        if (shared_id /= 0) then
            call mpi_win_shared_query(shared_win, 0, windowsize, disp_unit, baseptr, ierr)
        endif
        call c_f_pointer(baseptr, shared_array, dim_shape)
    end subroutine

    subroutine allocate_shared_memory_for_1dim_int_array(dim_shape, shared_array, shared_win)
        integer, intent(in) :: dim_shape
        integer, intent(out) :: shared_win
        integer, pointer, intent(out) :: shared_array(:)

        type(C_PTR) :: baseptr
        integer :: ierr, disp_unit
        integer(kind=mpi_address_kind) :: windowsize

        if(shared_id == 0) then
            windowsize =  int(dim_shape, MPI_ADDRESS_KIND) * 4_MPI_ADDRESS_KIND
        else
            windowsize = 0_MPI_ADDRESS_KIND
        endif
        disp_unit = 1

        call mpi_win_allocate_shared(windowsize, disp_unit, mpi_info_null, shared_comm_world, baseptr, shared_win, ierr)
        if (shared_id /= 0) then
            call mpi_win_shared_query(shared_win, 0, windowsize, disp_unit, baseptr, ierr)
        endif
        call c_f_pointer(baseptr, shared_array, [dim_shape])
    end subroutine

    subroutine deallocate_shared_memory()
        integer :: ierr
        ! call mpi_win_free(shared_win, ierr)
        call mpi_comm_free(shared_comm_world, ierr)
    end subroutine
end module
