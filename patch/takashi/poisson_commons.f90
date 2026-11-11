module poisson_commons
  use amr_commons
  use poisson_parameters

  real(dp),allocatable,dimension(:)  ::phi,phi_old       ! Potential
  real(dp),allocatable,dimension(:)  ::rho               ! Density
  real(dp),allocatable,dimension(:,:)::f                 ! 3-force

  real(dp),allocatable,dimension(:)  ::rho_top   ! Density at last CIC level

  ! Multigrid lookup table for amr -> mg index mapping
  integer, allocatable, dimension(:) :: lookup_mg   ! Lookup table
  integer, allocatable, dimension(:) :: lookup_mg_backup ! Persistent copy of lookup_mg (upstream PR #345)

  ! Communicator arrays for multigrid levels
#ifdef LIGHT_MPI_COMM
  type(communicator_light), allocatable, dimension(:,:) :: active_mg
  type(communicator_varoct), allocatable, dimension(:) :: emission_mg
#else
  type(communicator), allocatable, dimension(:,:) :: active_mg
  type(communicator), allocatable, dimension(:,:) :: emission_mg
#endif

  ! Send/recv Multigrid temporary communicator (light) used in build_parent_comms_mg subroutine
  type communicator_mg
    integer                       ::ngrid
    integer, dimension(:), pointer::igrid
  end type communicator_mg

  ! Minimum MG level
  integer :: levelmin_mg

  ! Multigrid safety switch
  logical, allocatable, dimension(:) :: safe_mode

  ! Automatic solver choice (cg_auto): per level
  logical, allocatable, dimension(:) :: use_cg        ! current method (.true. = CG)
  logical, allocatable, dimension(:) :: cg_failed     ! CG hit itermax once: never use CG again
  real(dp),allocatable, dimension(:) :: t_solve_mg    ! last local solve time with multigrid (<0: none)
  real(dp),allocatable, dimension(:) :: t_solve_cg    ! last local solve time with CG (<0: none)
  logical, allocatable, dimension(:) :: cg_trial_due  ! next solve at this level tries the other method
  logical, allocatable, dimension(:) :: cg_pending    ! trial done, waiting for the current method's time
  integer :: cg_auto_step=-1                           ! coarse step at which the trial flags were last set

  ! Multipole coefficients
  real(dp),dimension(1:ndim+1)::multipole

end module poisson_commons
