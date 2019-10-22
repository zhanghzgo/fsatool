module math
  implicit none
contains
  subroutine math_tpm_pi(tpm, nstate, pi)
    integer :: nstate
    real*8, intent(in) :: tpm(nstate, nstate)
    real*8, intent(out) :: pi(nstate)
    real*8 :: eigenvalue(nstate), eigenvector(nstate, nstate)
    call math_solve_eigenproblem_tpm(tpm, eigenvalue, eigenvector, pi)
  end subroutine math_tpm_pi

  subroutine math_solve_eigenproblem_tpm(tpm, eigenvalue, eigenvector, pi)
    real*8, dimension(:, :), intent(in) :: tpm
    real*8, dimension(:, :), intent(out) :: eigenvector
    real*8, dimension(:),intent(out) :: eigenvalue, pi
    integer :: i,j,info, tcmdim
    real*8, allocatable :: eigenvalue_img(:), work(:), tpvec(:), tpmcopy(:,:), eigenleft(:,:)
    real*8 :: temp

    tcmdim = size(tpm, 1)
    allocate(work(4*tcmdim), eigenvalue_img(tcmdim), tpvec(tcmdim), tpmcopy(tcmdim, tcmdim), eigenleft(tcmdim, tcmdim))
    tpmcopy = tpm

    call dgeev('V', 'V', tcmdim, tpmcopy, tcmdim, eigenvalue, eigenvalue_img, eigenleft, &
         tcmdim, eigenvector, tcmdim, work, 4*tcmdim, info )

    do i=1,tcmdim
       do j=i+1,tcmdim
          if ( abs(eigenvalue(j)) > abs(eigenvalue(i)) ) then
             temp = eigenvalue(i); eigenvalue(i) = eigenvalue(j); eigenvalue(j)=temp
             tpvec = eigenvector(:,i); eigenvector(:,i) = eigenvector(:,j); eigenvector(:,j)=tpvec
             tpvec = eigenleft(:, i); eigenleft(:,i) = eigenleft(:, j);eigenleft(:, j) = tpvec
          end if
       end do
    end do
    ! in case there is -1 eigenvalue
    do i = 1, tcmdim
       if(abs(eigenvalue(i)-1)<1e-6) then
          pi(:) = eigenleft(:,1)
          temp = sum(pi); pi = pi/temp
          exit
       endif
    enddo
    deallocate(work, eigenvalue_img, tpvec, tpmcopy,eigenleft)
  end subroutine math_solve_eigenproblem_tpm

  subroutine math_inv(A, Ainv)
    real*8, dimension(:, :), intent(in) :: A
    real*8, dimension(size(A, 1), size(A, 1)) :: Ainv
    real*8, dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n, n, Ainv, n, ipiv, info)

    if (info /= 0) then
       stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n, Ainv, n, ipiv, work, n, info)

    if (info /= 0) then
       stop 'Matrix inversion failed!'
    end if
  end subroutine math_inv

  subroutine math_nelmin ( fn, n, start, xmin, ynewlo, reqmin, step, konvge, kcount, &
       icount, numres, ifault )

    !*****************************************************************************80
    !
    !! NELMIN minimizes a function using the Nelder-Mead algorithm.
    !
    !  Discussion:
    !
    !    This routine seeks the minimum value of a user-specified function.
    !
    !    Simplex function minimisation procedure due to Nelder and Mead (1965),
    !    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
    !    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
    !    25, 97) and Hill(1978, 27, 380-2)
    !
    !    The function to be minimized must be defined by a function of
    !    the form
    !
    !      function fn ( x, f )
    !      real ( kind = 8 ) fn
    !      real ( kind = 8 ) x(*)
    !
    !    and the name of this subroutine must be declared EXTERNAL in the
    !    calling routine and passed as the argument FN.
    !
    !    This routine does not include a termination test using the
    !    fitting of a quadratic surface.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    27 February 2008
    !
    !  Author:
    !
    !    Original FORTRAN77 version by R ONeill.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    John Nelder, Roger Mead,
    !    A simplex method for function minimization,
    !    Computer Journal,
    !    Volume 7, 1965, pages 308-313.
    !
    !    R ONeill,
    !    Algorithm AS 47:
    !    Function Minimization Using a Simplex Procedure,
    !    Applied Statistics,
    !    Volume 20, Number 3, 1971, pages 338-345.
    !
    !  Parameters:
    !
    !    Input, external FN, the name of the function which evaluates
    !    the function to be minimized.
    !
    !    Input, integer ( kind = 4 ) N, the number of variables.
    !    0 < N is required.
    !
    !    Input/output, real ( kind = 8 ) START(N).  On input, a starting point
    !    for the iteration.  On output, this data may have been overwritten.
    !
    !    Output, real ( kind = 8 ) XMIN(N), the coordinates of the point which
    !    is estimated to minimize the function.
    !
    !    Output, real ( kind = 8 ) YNEWLO, the minimum value of the function.
    !
    !    Input, real ( kind = 8 ) REQMIN, the terminating limit for the variance
    !    of the function values.  0 < REQMIN is required.
    !
    !    Input, real ( kind = 8 ) STEP(N), determines the size and shape of the
    !    initial simplex.  The relative magnitudes of its elements should reflect
    !    the units of the variables.
    !
    !    Input, integer ( kind = 4 ) KONVGE, the convergence check is carried out
    !    every KONVGE iterations. 0 < KONVGE is required.
    !
    !    Input, integer ( kind = 4 ) KCOUNT, the maximum number of function
    !    evaluations.
    !
    !    Output, integer ( kind = 4 ) ICOUNT, the number of function evaluations
    !    used.
    !
    !    Output, integer ( kind = 4 ) NUMRES, the number of restarts.
    !
    !    Output, integer ( kind = 4 ) IFAULT, error indicator.
    !    0, no errors detected.
    !    1, REQMIN, N, or KONVGE has an illegal value.
    !    2, iteration terminated because KCOUNT was exceeded without convergence.
    !
    implicit none

    integer ( kind = 4 ) n

    real ( kind = 8 ), parameter :: ccoeff = 0.5D+00
    real ( kind = 8 ) del
    real ( kind = 8 ), parameter :: ecoeff = 2.0D+00
    real ( kind = 8 ), parameter :: eps = 0.001D+00
    real ( kind = 8 ), external :: fn
    integer ( kind = 4 ) i
    integer ( kind = 4 ) icount
    integer ( kind = 4 ) ifault
    integer ( kind = 4 ) ihi
    integer ( kind = 4 ) ilo
    integer ( kind = 4 ) j
    integer ( kind = 4 ) jcount
    integer ( kind = 4 ) kcount
    integer ( kind = 4 ) konvge
    integer ( kind = 4 ) l
    integer ( kind = 4 ) numres
    real ( kind = 8 ) p(n,n+1)
    real ( kind = 8 ) p2star(n)
    real ( kind = 8 ) pbar(n)
    real ( kind = 8 ) pstar(n)
    real ( kind = 8 ), parameter :: rcoeff = 1.0D+00
    real ( kind = 8 ) reqmin
    real ( kind = 8 ) rq
    real ( kind = 8 ) start(n)
    real ( kind = 8 ) step(n)
    real ( kind = 8 ) x
    real ( kind = 8 ) xmin(n)
    real ( kind = 8 ) y(n+1)
    real ( kind = 8 ) y2star
    real ( kind = 8 ) ylo
    real ( kind = 8 ) ynewlo
    real ( kind = 8 ) ystar
    real ( kind = 8 ) z
    !
    !  Check the input parameters.
    !
    if ( reqmin <= 0.0D+00 ) then
       ifault = 1
       return
    end if

    if ( n < 1 ) then
       ifault = 1
       return
    end if

    if ( konvge < 1 ) then
       ifault = 1
       return
    end if
    !
    !  Initialization.
    !
    icount = 0
    numres = 0
    jcount = konvge
    del = 1.0D+00
    rq = reqmin * real ( n, kind = 8 )
    !
    !  Initial or restarted loop.
    !
    do

       p(1:n,n+1) = start(1:n)
       y(n+1) = fn ( start, n )
       icount = icount + 1
       !
       !  Define the initial simplex.
       !
       do j = 1, n
          x = start(j)
          start(j) = start(j) * (step(j) + 1)
          p(1:n,j) = start(1:n)
          y(j) = fn ( start, n )
          icount = icount + 1
          start(j) = x
       end do
       !
       !  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
       !  the vertex of the simplex to be replaced.
       !
       ilo = minloc ( y(1:n+1), 1 )
       ylo = y(ilo)
       !
       !  Inner loop.
       !
       do while ( icount < kcount )
          !
          !  YNEWLO is, of course, the HIGHEST value???
          !
          ihi = maxloc ( y(1:n+1), 1 )
          ynewlo = y(ihi)
          !
          !  Calculate PBAR, the centroid of the simplex vertices
          !  excepting the vertex with Y value YNEWLO.
          !
          do i = 1, n
             pbar(i) = ( sum ( p(i,1:n+1) ) - p(i,ihi) ) / real ( n, kind = 8 )
          end do
          !
          !  Reflection through the centroid.
          !
          pstar(1:n) = pbar(1:n) + rcoeff * ( pbar(1:n) - p(1:n,ihi) )
          ystar = fn ( pstar, n )
          icount = icount + 1
          !
          !  Successful reflection, so extension.
          !
          if ( ystar < ylo ) then

             p2star(1:n) = pbar(1:n) + ecoeff * ( pstar(1:n) - pbar(1:n) )
             y2star = fn ( p2star, n )
             icount = icount + 1
             !
             !  Retain extension or contraction.
             !
             if ( ystar < y2star ) then
                p(1:n,ihi) = pstar(1:n)
                y(ihi) = ystar
             else
                p(1:n,ihi) = p2star(1:n)
                y(ihi) = y2star
             end if
             !
             !  No extension.
             !
          else

             l = 0
             do i = 1, n + 1
                if ( ystar < y(i) ) then
                   l = l + 1
                end if
             end do

             if ( 1 < l ) then

                p(1:n,ihi) = pstar(1:n)
                y(ihi) = ystar
                !
                !  Contraction on the Y(IHI) side of the centroid.
                !
             else if ( l == 0 ) then

                p2star(1:n) = pbar(1:n) + ccoeff * ( p(1:n,ihi) - pbar(1:n) )
                y2star = fn ( p2star, n )
                icount = icount + 1
                !
                !  Contract the whole simplex.
                !
                if ( y(ihi) < y2star ) then

                   do j = 1, n + 1
                      p(1:n,j) = ( p(1:n,j) + p(1:n,ilo) ) * 0.5D+00
                      xmin(1:n) = p(1:n,j)
                      y(j) = fn ( xmin, n )
                      icount = icount + 1
                   end do

                   ilo = minloc ( y(1:n+1), 1 )
                   ylo = y(ilo)

                   cycle
                   !
                   !  Retain contraction.
                   !
                else
                   p(1:n,ihi) = p2star(1:n)
                   y(ihi) = y2star
                end if
                !
                !  Contraction on the reflection side of the centroid.
                !
             else if ( l == 1 ) then

                p2star(1:n) = pbar(1:n) + ccoeff * ( pstar(1:n) - pbar(1:n) )
                y2star = fn ( p2star, n )
                icount = icount + 1
                !
                !  Retain reflection?
                !
                if ( y2star <= ystar ) then
                   p(1:n,ihi) = p2star(1:n)
                   y(ihi) = y2star
                else
                   p(1:n,ihi) = pstar(1:n)
                   y(ihi) = ystar
                end if

             end if

          end if
          !
          !  Check if YLO improved.
          !
          if ( y(ihi) < ylo ) then
             ylo = y(ihi)
             ilo = ihi
          end if

          jcount = jcount - 1

          if ( 0 < jcount ) then
             cycle
          end if
          !
          !  Check to see if minimum reached.
          !
          if ( icount <= kcount ) then

             jcount = konvge

             x = sum ( y(1:n+1) ) / real ( n + 1, kind = 8 )
             z = sum ( ( y(1:n+1) - x )**2 )

             if ( z <= rq ) then
                exit
             end if

          end if

       end do
       !
       !  Factorial tests to check that YNEWLO is a local minimum.
       !
       xmin(1:n) = p(1:n,ilo)
       ynewlo = y(ilo)

       if ( kcount < icount ) then
          ifault = 2
          exit
       end if

       ifault = 0

       do i = 1, n
          del = step(i) * eps
          xmin(i) = xmin(i) + del
          z = fn ( xmin, n )
          icount = icount + 1
          if ( z < ynewlo ) then
             ifault = 2
             exit
          end if
          xmin(i) = xmin(i) - del - del
          z = fn ( xmin, n )
          icount = icount + 1
          if ( z < ynewlo ) then
             ifault = 2
             exit
          end if
          xmin(i) = xmin(i) + del
       end do

       if ( ifault == 0 ) then
          exit
       end if
       !
       !  Restart the procedure.
       !
       start(1:n) = xmin(1:n)
       del = eps
       numres = numres + 1

    end do
    return
  end subroutine math_nelmin

   subroutine math_euclidean_distance(coor1, coor2, dist2, n)
      integer,intent(in) :: n
      real*8, intent(in) :: coor1(n), coor2(n)
      real*8, intent(out) :: dist2
      dist2 = dot_product(coor1-coor2, coor1-coor2)
   end subroutine

   subroutine math_lRMSD(coor1_raw, coor2_raw, rmsd, ndim)
      integer, intent(in) :: ndim
      real*8, intent(in) :: coor1_raw(ndim),coor2_raw(ndim)
      real*8, intent(out) :: rmsd

      integer :: inum, i, j, info, idim
      real*8 :: tpvec(3),umatrix(3,3),vmatrix(3,3),wvector(3),tpmatrix(3,3), work(15), originmatrix(3,3)
      real*8 :: dsign, temp
      real*8 :: tpcoor(ndim), coor1(ndim), coor2(ndim)

      coor1 = coor1_raw
      coor2 = coor2_raw

      if ( mod(ndim, 3) /= 0 ) STOP "wrong ndim!"

      tpvec = 0.0d0

      do i = 1, ndim, 3
         tpvec(1:3) = tpvec(1:3) + coor1(i:i+2)
      end do
      tpvec = tpvec / dble(ndim/3)
      do i =1, ndim, 3
         coor1(i:i+2) = coor1(i:i+2) - tpvec(1:3)
      end do

      tpvec = 0.0d0
      do i = 1, ndim, 3
         tpvec = tpvec + coor2(i:i+2)
      end do
      tpvec = tpvec / dble(ndim/3)
      do i = 1, ndim, 3
         coor2(i:i+2) = coor2(i:i+2) - tpvec
      end do

      originmatrix = 0.0d0

      do i = 1, 3
         do j = 1, 3
            do idim = 1, ndim, 3
               originmatrix(i, j) = originmatrix(i, j) + coor1(idim + i - 1) * coor2(idim + j - 1)
            end do
         end do
      end do
      call dgesvd("A", "A", 3, 3, originmatrix, 3, wvector, umatrix, 3, vmatrix, 3, work, 15 ,info)
      vmatrix = transpose(vmatrix)

      dsign = 1.0d0
      tpvec(1:3) = (/1.0d0,1.0d0,dsign/)
      do i = 1, 3
         do j = 1, 3
            temp = 0.0d0
            do inum = 1, 3
               temp = temp + vmatrix(i, inum) * tpvec(inum) * umatrix(j, inum)
            end do
            tpmatrix(i,j) = temp
         end do
      end do

      rmsd = 0d0; tpcoor = 0d0
      do idim = 1, ndim, 3
         do i = 1, 3
            tpcoor(idim+i-1) = dot_product(tpmatrix(i,1:3),coor1(idim:idim+2))
            rmsd = rmsd + (tpcoor(idim+i-1) - coor2(idim+i-1))**2
         end do
      end do
      rmsd = sqrt( rmsd / dble (ndim / 3) )
   end subroutine

   subroutine math_normalization(traj, normalize_traj, ndim, nsnap)
      real*8, intent(in) :: traj(ndim, nsnap)
      real*8, intent(out) :: normalize_traj(ndim, nsnap)
      integer, intent(in) :: ndim, nsnap

      ! real*8 :: normalize_traj(ndim, nsnap)
      real*8 :: average(ndim), standard_deviation(ndim)
      real*8 :: total
      integer :: i, j

      average = sum(traj, dim=2)/nsnap
      do i = 1, ndim
         total = 0.0d0
         do j = 1, nsnap
            total = total + (traj(i,j) - average(i)) * (traj(i,j)-average(i))
         enddo
         standard_deviation(i) = sqrt(total/dble((nsnap-1)))
      enddo
      do i = 1, ndim
         do j = 1, nsnap
            normalize_traj(i, j) = (traj(i,j) - average(i))/(standard_deviation(i))
         enddo
      enddo
   end subroutine
end module math

