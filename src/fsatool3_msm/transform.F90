module transform
  use util
  implicit none
contains

  subroutine mod_transform(inputfile, resultfile)
    character(256) :: datafile, transformmethod, inputfile, resultfile
    integer :: nsnap, nfeature, ncomponent, iofile, ierr, lagtime, i
    real*8, allocatable :: data(:, :), maparray(:, :)

    namelist /transform/ datafile, transformmethod, nsnap, nfeature, ncomponent, lagtime

    call loginfo("Transformation")

    call getfreeunit(iofile)
    open(unit=iofile, file=trim(inputfile), action="read")
    read(iofile, nml=transform, iostat=ierr)
    if(ierr < 0) call errormsg("error in reading the transform namelist")
    close(iofile)

    allocate(data(nsnap, nfeature), maparray(nsnap, ncomponent))

    call getfreeunit(iofile)
    open(unit=iofile, file=trim(datafile), action="read")
    do i = 1, nsnap
      read(iofile, *) data(i, :)
    enddo
    close(iofile)

    if(trim(transformmethod) == "pca") then
      call pca(data, nsnap, nfeature, ncomponent, maparray)
    else if(trim(transformmethod) == "tica") then
      call tica(data, nsnap, nfeature, ncomponent, maparray, lagtime)
    else
      call errormsg("transformmethod must be pca or tica")
    endif

    open(unit=iofile, file=trim(resultfile), action="write")
    do i = 1, nsnap
      write(iofile, *) maparray(i, :)
    enddo
    close(iofile)
    call loginfo()

  end subroutine

  subroutine pca(tparray, nsnap, nfeature, ncomponent, maparray)
    integer :: i, j, nfeature, info, nsnap, ncomponent
    real*8 :: tparray(nsnap, nfeature), average(nfeature), cov(nfeature, nfeature), temp
    real*8 :: eigenvector(nfeature, nfeature), maparray(nsnap, ncomponent)
    real*8 :: eigenvalue_real(nfeature), eigenvalue_img(nfeature), work(4*nfeature), tempmatrix(nfeature,nfeature)
    integer :: component_index(nfeature), tempint
    do i = 1, nfeature
      average(i) = sum(tparray(:,i))/nsnap
    enddo
    do i = 1, nfeature
      do j = 1, nfeature
          cov(i,j) = sum((tparray(:, i) - average(i)) * (tparray(:, j) - average(j)))/(nsnap - 1)
      enddo
    enddo
    call dgeev("N","V", nfeature, cov, nfeature, eigenvalue_real, eigenvalue_img, &
        tempmatrix, nfeature, eigenvector, nfeature, work, 4*nfeature, info)
    call get_feature_index(eigenvalue_real, nfeature, component_index)
    do i = 1, ncomponent
      tempmatrix(:, i) = eigenvector(:, component_index(i))
    enddo
    maparray = matmul(tparray, tempmatrix(:, 1:ncomponent))
    write(*, "(a,10f8.3)")"PCA eigenvalue", eigenvalue_real
    write(*, "(a,10f8.3)")"PCA eigenvector", eigenvector
  end subroutine pca

  subroutine tica(data, npoint, nfeature, ncomponent, maparray, lagtime)
    use math, only: math_inv
    integer :: lagtime, npoint, nfeature, ncomponent, i, j, offset, info, component_index(nfeature)
    real*8 :: data(npoint, nfeature), cov(nfeature, nfeature), cov_tau(nfeature, nfeature), exx(nfeature, nfeature)
    real*8 :: data_mean(nfeature, 1), cov_inv(nfeature, nfeature)
    real*8 :: eigenvalue_real(nfeature), eigenvalue_img(nfeature), tempmatrix(nfeature, nfeature)
    real*8 :: eigenvector(nfeature,nfeature), work(4*nfeature), maparray(npoint, ncomponent)

    offset = npoint - lagtime
    exx = matmul(transpose(data(:offset, :)), data(:offset, :)) &
        + matmul(transpose(data(lagtime+1:, :)), data(lagtime+1:, :))
    exx = exx / dble(offset * 2)
    do i=1, nfeature
      data_mean(i, 1) = sum(data(:offset, i) + data(lagtime+1:, i)) / dble(offset * 2)
    enddo
    cov  = exx - matmul(data_mean, transpose(data_mean))
    exx = matmul(transpose(data(:offset, :)), data(lagtime+1:, :)) / dble(offset)
    exx = (transpose(exx) + exx) /2.0d0
    cov_tau = exx - matmul(data_mean, transpose(data_mean))

    call dpotrf("L", nfeature, cov, nfeature, info)
    do i = 1, nfeature
      do j = i+1, nfeature
          cov(i,j) = 0.0d0
      enddo
    enddo
    call math_inv(cov, cov_inv)
    !call inv(cov, cov_inv, nfeature)
    !cov_inv = cov
    tempmatrix = matmul(cov_inv, matmul(cov_tau, transpose(cov_inv)))
    call dgeev("N","V", nfeature, tempmatrix, nfeature, eigenvalue_real, eigenvalue_img, &
        tempmatrix, nfeature, eigenvector, nfeature, work, 4*nfeature, info)
    eigenvector = matmul(transpose(cov_inv), eigenvector)

    do i = 1,nfeature
      eigenvector(:,i) = eigenvalue_real(i) * eigenvector(:, i)
    enddo
    call get_feature_index(eigenvector, nfeature, component_index)
    do i = 1, ncomponent
      tempmatrix(:, i) = eigenvector(:, component_index(i))
    enddo
    do i = 1, nfeature
      data(:, i) = data(:, i) - data_mean(i, 1)
    end do
    maparray = matmul(data, tempmatrix(:, 1:ncomponent))
    write(*, "(a,10f8.3)")"tica eigenvalue", eigenvalue_real
    write(*, "(a,10f8.3)")"tica eigenvector", eigenvector
  end subroutine tica

subroutine get_feature_index(tparray, nfeature, sortindex)
  integer :: nfeature
  real*8 :: temp, tparray(nfeature)
  integer :: tempint, sortindex(nfeature)
  integer :: i, j

  do i = 1, nfeature
     sortindex(i) = i
  enddo
  do i = 1, nfeature
     do j = i+1, nfeature
        if (tparray(j) > tparray(i)) then
           temp = tparray(i)
           tparray(i) = tparray(j)
           tparray(j)  = temp

           tempint = sortindex(i)
           sortindex(i) = sortindex(j)
           sortindex(j) = tempint
        endif
     enddo
  enddo
end subroutine

end module