module combine
    use transform
    use cluster
contains

subroutine pca_cluster()
    integer :: nsnap, nfeature, ncomponent
    real*8 :: array(nsnap, nfeature)
    real*8 :: maparray(nsnap, ncomponent)
    call pca(array, nsnap, nfeature, ncomponent, maparray)

end subroutine


end module
