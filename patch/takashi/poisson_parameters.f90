module poisson_parameters
  use amr_parameters

  ! Convergence criterion for Poisson solvers
  real(dp)::epsilon=1.0D-4

  ! Type of force computation
  integer ::gravity_type=0

  ! Gravity parameters
  real(dp),dimension(1:10)::gravity_params=0

  ! Maximum level for CIC dark matter interpolation
  integer :: cic_levelmax=0

  ! Min level for CG solver
  ! level < cg_levelmin uses fine multigrid
  ! level >=cg_levelmin uses conjugate gradient
  integer :: cg_levelmin=999

  ! Automatic per-level choice between multigrid and CG (patch/takashi, 2026-11-11).
  ! With cg_auto=.true., cg_levelmin only gives the initial choice; every
  ! cg_auto_every coarse steps each level (> levelmin) solves once with the
  ! other method, the solve times (min over ranks) are compared, and the
  ! faster method is kept. A switch needs the other method to be faster by
  ! the fraction cg_auto_margin (hysteresis). Results change at round-off
  ! level when the method switches.
  logical :: cg_auto=.false.
  integer :: cg_auto_every=20
  real(dp):: cg_auto_margin=0.1d0

  ! Gauss-Seidel smoothing sweeps for fine multigrid
  integer, parameter :: ngs_fine   = 2
  integer, parameter :: ngs_coarse = 2

  ! Number of multigrid cycles for coarse levels *in safe mode*
  !   1 is the fastest,
  !   2 is slower but can give much better convergence in some cases
  integer, parameter :: ncycles_coarse_safe = 1

end module poisson_parameters
