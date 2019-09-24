program main
  use mod_global
  use cluster
  use fileio
  use cluster
  use markov
  use tpt
  use getcmd

  implicit none
  integer :: ierr
  logical :: msm
  include 'mpif.h'

  call mod_global_mpiinit()
  call get_program(msm)

  if (msm) then
    call fileio_init_parameters()
    call mod_global_readtraj()
    if (ifreadcoor .eqv. .false.) then
        call fileio_readclusterinfo()
    else
        call cluster_analysis()
        call fileio_writeclusterinfo()
    endif
    call markov_analysis()
    call fileio_writestateinfo()
    call tpt_msm_analysis()
  endif

  call mpi_finalize(ierr)
end program