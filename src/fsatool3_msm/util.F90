module util
  implicit none 
  public :: dsort, isort, isSorted

contains 
  recursive subroutine dsort(array, subarray, lo, hi)
    real*8, intent(inout) :: array(:), subarray(:)
    integer :: lo, hi, mid

    if(hi <= lo) return
    mid = (lo + hi) / 2
    call dsort(array, subarray, lo, mid)
    call dsort(array, subarray, mid+1, hi)
    call dmerge(array, subarray, lo, mid, hi)
  end subroutine dsort

  subroutine dmerge(array, subarray, lo, mid, hi)
    real*8, intent(inout) :: array(:), subarray(:)
    integer :: lo, hi, mid, i, j , k
    do i = lo, hi 
          subarray(i) = array(i)
    enddo
    i = lo; j = mid+1;
    do k = lo, hi
          if(i > mid) then
            array(k) = subarray(j); j = j+1
          else if(j > hi) then 
              array(k) = subarray(i); i = i+1
          else if(subarray(j) < subarray(i))  then
            array(k) = subarray(j); j = j+1
          else 
            array(k) = subarray(i); i = i+1
          endif
    enddo
  end subroutine

  recursive subroutine isort(array, subarray, lo, hi)
    integer, intent(inout) :: array(:), subarray(:)
    integer :: lo, hi, mid

    if(hi <= lo) return
    mid = (lo + hi) / 2
    call isort(array, subarray, lo, mid)
    call isort(array, subarray, mid+1, hi)
    call imerge(array, subarray, lo, mid, hi)
  end subroutine isort

  subroutine imerge(array, subarray, lo, mid, hi)
      integer, intent(inout) :: array(:), subarray(:)
      integer :: lo, hi, mid, i, j , k
      do i = lo, hi 
            subarray(i) = array(i)
      enddo
      i = lo; j = mid+1;
      do k = lo, hi
            if(i > mid) then
              array(k) = subarray(j); j = j+1
            else if(j > hi) then 
               array(k) = subarray(i); i = i+1
            else if(subarray(j) < subarray(i))  then
              array(k) = subarray(j); j = j+1
            else 
              array(k) = subarray(i); i = i+1
            endif
      enddo
    end subroutine

  function isSorted(array) result(res)
    real*8 :: array(:)
    integer :: i
    logical :: res
    res = .true. 
    if(size(array) <= 1) return
    do i = 2, size(array)
      if(array(i) < array(i-1)) then
        res = .false.
        return
      endif
    enddo
  end function

subroutine getfreeunit(tpunit)
  implicit none
  integer tpunit
  logical ifused

  tpunit = 10; ifused = .true.
  do
     if ( ifused .eqv. .false. ) exit
     tpunit = tpunit + 1
     if (tpunit > 99) call errormsg ("unit has been used")
     inquire (unit=tpunit,opened=ifused)
  end do
end subroutine getfreeunit

SUBROUTINE init_random_seed()
  INTEGER :: i, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed
  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))
  CALL SYSTEM_CLOCK(COUNT=clock)
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)
  DEALLOCATE(seed)
END SUBROUTINE init_random_seed

subroutine errormsg(msg)
  character(len=*),intent(in) :: msg
  write(*, "(a)")"Error: "//msg
  flush(6)
  stop
end subroutine errormsg

subroutine loginfo(msg)
  character(len=*), intent(in),optional :: msg
  integer :: lenstr
  if (present(msg)) then
     lenstr = len(msg)
     write(*, "(a)") repeat("*", 80)
     write(*, "(a)") repeat(" ", (80 - lenstr)/2) // msg
  else
     write(*, "(a)") repeat("*", 80)
     write(*, *)
  endif
end subroutine loginfo
end module