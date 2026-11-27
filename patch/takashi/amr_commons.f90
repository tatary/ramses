module amr_commons
  use amr_parameters

  logical::output_done=.false.                  ! Output just performed
  logical::init=.false.                         ! Set up or run
  logical::balance=.false.                      ! Load balance or run
  logical::shrink=.false.                       ! Shrink mesh or run
  integer::nstep=0                              ! Time step
  integer::nstep_coarse=0                       ! Coarse step
  integer::nstep_coarse_old=0                   ! Old coarse step
  integer::nstep_restart=-1                     ! nstep restored from the restart dump (-1: fresh run); ugvtk files with nstep >= this were written by an abandoned branch (output_analysis)
  integer::nflag,ncreate,nkill                  ! Refinements
  integer::ncoarse                              ! nx.ny.nz
  integer::ngrid_current                        ! Actual number of octs

  real(dp)::emag_tot=0.0D0                      ! Total magnetic energy
  real(dp)::ekin_tot=0.0D0                      ! Total kinetic energy
  real(dp)::eint_tot=0.0D0                      ! Total internal energy
  real(dp)::epot_tot=0.0D0                      ! Total potential energy
  real(dp)::epot_tot_old=0.0D0                  ! Old potential energy
  real(dp)::epot_tot_int=0.0D0                  ! Time integrated potential
  real(dp)::einit=0.0D0                         ! Initial total energy
  real(dp)::aexp_old=1.0D0                      ! Old expansion factor
  real(dp)::rho_tot=0.0D0                       ! Mean density in the box
  real(dp)::t=0.0D0                             ! Time variable

  ! executable identification
  CHARACTER(LEN=80)::builddate,patchdir
  CHARACTER(LEN=80)::gitrepo,gitbranch,githash

  ! Save namelist filename
  CHARACTER(LEN=80)::namelist_file

  ! MPI variables
  integer::ncpu,ndomain,myid,overload=1
#ifdef SOLVERmhd
  logical::mhd_feedback_active=.false.   ! Set by thermal_feedback when a magnetic loop was injected (Mizuki Ono)
#endif

  ! Friedman model variables
  integer::n_frw
  real(dp),allocatable,dimension(:)::aexp_frw,hexp_frw,tau_frw,t_frw

  ! Initial conditions parameters from grafic
  integer                  ::nlevelmax_part
  !real(dp)                 ::aexp_ini=10
  real(dp)                 ::aexp_ini=1.0d0
  real(dp),dimension(1:MAXLEVEL)::dfact=1.0d0,astart
  real(dp),dimension(1:MAXLEVEL)::vfact
  real(dp),dimension(1:MAXLEVEL)::xoff1,xoff2,xoff3,dxini
  integer ,dimension(1:MAXLEVEL)::n1,n2,n3

  ! Level related arrays
  real(dp),dimension(1:MAXLEVEL)::dtold,dtnew ! Time step at each level
  real(dp),dimension(1:MAXLEVEL)::rho_max     ! Maximum density at each level
  integer ,dimension(1:MAXLEVEL)::nsubcycle=2 ! Subcycling at each level

  ! Pointers for each level linked list
  integer,allocatable,dimension(:,:)::headl
  integer,allocatable,dimension(:,:)::taill
  integer,allocatable,dimension(:,:)::numbl
  integer(i8b),allocatable,dimension(:,:)::numbtot

  ! Pointers for each level boundary linked list
  integer,allocatable,dimension(:,:)::headb
  integer,allocatable,dimension(:,:)::tailb
  integer,allocatable,dimension(:,:)::numbb

  ! Pointers for free memory grid linked list
  integer::headf,tailf,numbf,used_mem,used_mem_tot

  ! Tree arrays
  real(dp),allocatable,dimension(:,:)::xg      ! grids position
  integer ,allocatable,dimension(:,:)::nbor    ! neighboring father cells
  integer ,allocatable,dimension(:)  ::father  ! father cell
  integer ,allocatable,dimension(:)  ::next    ! next grid in list
  integer ,allocatable,dimension(:)  ::prev    ! previous grid in list
  integer ,allocatable,dimension(:)  ::son     ! sons grids
  integer ,allocatable,dimension(:)  ::flag1   ! flag for refine
  integer ,allocatable,dimension(:)  ::flag2   ! flag for expansion
  ! FG star formation (sf_fg, star_formation_fg.f90): Pop III mass per cell (allocated in init_amr if sf_fg)
  real(dp),allocatable,dimension(:)  ::flag3
  ! FG star formation: hydro variables at the front of the update (unew + virtual-boundary fluxes),
  ! used by the non-local gas depletion and the SFC model; allocated on first use (my_utils:ensure_ufront)
  real(dp),allocatable,dimension(:,:)::ufront

  ! Global indexing
  integer ,allocatable,dimension(:)  ::cpu_map  ! domain decomposition
  integer ,allocatable,dimension(:)  ::cpu_map2 ! new domain decomposition

  ! Hilbert key
  real(qdp),allocatable,dimension(:)::hilbert_key
  real(qdp),allocatable,dimension(:)::bound_key,bound_key2
  ! two-tier domain decomposition (ported from patch/ks 2026-11-13, docs/merge_ks.md)
  real(qdp),allocatable,dimension(:)::bound_key_f,bound_key2_f ! fine-tier boundaries (grids at level > two_tier_ls), current / new
  integer::two_tier_ls=0,two_tier_ls2=0                        ! last level of the coarse tier, current / new (nlevelmax: single partition)
  ! Fine footprints: the fine tier owns a key only inside the footprint of its domain, i.e. the
  ! key ranges of the level-two_tier_ls cells that contained fine grids at the last load balance,
  ! clipped to the domain interval (sorted, disjoint, grouped by domain: intervals
  ! foot_first(idom-1)+1..foot_first(idom) belong to domain idom). A fine-tier key outside every
  ! footprint (a fine grid born since the last balance in an empty region) belongs to the coarse
  ! owner of the key. So a cpu holds the coarse ancestors of its footprints only, not of its whole
  ! fine key interval.
  integer::two_tier_nfoot=0,two_tier_nfoot2=0
  real(qdp),allocatable,dimension(:)::foot_lo,foot_hi,foot_lo2,foot_hi2
  integer,allocatable,dimension(:)::foot_first,foot_first2
  real(qdp)                         ::order_all_min,order_all_max

  ! Recursive bisection
  real(dp),allocatable,dimension(:)    ::bisec_wall         ! bisection wall positions
  integer ,allocatable,dimension(:,:)  ::bisec_next         ! next 2 child cells in bisection
  integer::bisec_root                                       ! root of bisection tree

  integer,allocatable,dimension(:)     ::bisec_indx         ! map from leaf cell id to cpu id
  real(dp),allocatable,dimension(:,:)  ::bisec_cpubox_min   ! cpu domains boxes
  real(dp),allocatable,dimension(:,:)  ::bisec_cpubox_max
  real(dp),allocatable,dimension(:,:)  ::bisec_cpubox_min2  ! cpu domains boxes for new decomp
  real(dp),allocatable,dimension(:,:)  ::bisec_cpubox_max2

  integer,allocatable,dimension(:)     ::bisec_cpu_load     ! CPU loads (for stats)
  integer,allocatable,dimension(:,:)   ::bisec_hist         ! histograms for load computation
  integer,allocatable,dimension(:)     ::bisec_hist_bounds  ! histogram splitting boundaries
  integer,allocatable,dimension(:)     ::new_hist_bounds
  integer,allocatable,dimension(:)     ::bisec_ind_cell     ! histo swap id -> cell id map (big)
  integer,allocatable,dimension(:)     ::cell_level         ! store the level of the cells (big)

  real(dp)::bisec_res                                       ! resolution parameters
  integer ::bisec_nres

  ! Communication structure
  type communicator
    integer                            ::ngrid
    integer                            ::npart
    integer     ,dimension(:)  ,pointer::igrid
    integer     ,dimension(:,:),pointer::f
    real(kind=8),dimension(:,:),pointer::u
    integer(i8b),dimension(:,:),pointer::fp
    real(kind=8),dimension(:,:),pointer::up
#ifdef ATON
    real(kind=8),dimension(:,:),pointer::u_radiation
#endif
  end type communicator


#ifdef LIGHT_MPI_COMM
  ! ----------- Optimized MPI buffer communicator structures suggested by P. Wautelet on the IDRIS page ------------ !
  ! ----------- http://www.idris.fr/docs/docu/support-avance/ramses.html to reduce the memory footprint ------------ !
  ! ----------- D. Chapon (CEA Saclay - IRFU).                                                          ------------ !
  ! Actual communication structure holding the data allocatable arrays
  type point_comm
    integer          ,dimension(:),  pointer::igrid
    integer          ,dimension(:,:),pointer::f
    integer(kind=i8b),dimension(:,:),pointer::f8
    real(kind=8)     ,dimension(:,:),pointer::u
#ifdef ATON
    real(kind=8)     ,dimension(:,:),pointer::u_radiation
#endif
  end type point_comm

  ! Light communication (intermediate) structure
  type communicator_light
    integer                            ::ngrid
    integer                            ::npart
    type(point_comm), pointer          ::pcomm
  end type communicator_light

  type communicator_varoct
    integer :: nactive    ! Number of processes sharing data with current process
    integer :: ngrids_tot ! Total number of octs to share
    integer,     allocatable,dimension(:)  :: cpuid  ! Number of the active process
    integer,     allocatable,dimension(:)  :: ngrids ! Number of octs to share with each active process
    integer,     allocatable,dimension(:)  :: igrid
    real(kind=8),allocatable,dimension(:,:):: u
    integer,     allocatable,dimension(:,:):: f
#ifdef ATON
    real(kind=8),allocatable,dimension(:) ::u_radiation
#endif

  end type communicator_varoct

  type communicator_varpart
    integer :: nactive    ! Number of processes sharing data with current process
    integer :: nparts_tot ! Total number of particles to share
    integer,     allocatable,dimension(:)       :: cpuid  ! Number of the active process
    integer,     allocatable,dimension(:)       :: nparts ! Number of particles to share with each active process
    real(kind=8),allocatable,dimension(:,:)     :: u
    integer,     allocatable,dimension(:,:)     :: f
    integer(kind=i8b),allocatable,dimension(:,:):: f8
  end type communicator_varpart
#endif

  ! State of the non-blocking reverse exchanges
  ! (see make_virtual_reverse_dp_start/end in virtual_boundaries.f90)
  integer,parameter::nasync=2
  type async_exchange
     logical::active=.false.
     integer::ilevel=0,iv1=0,iv2=0,ivstep=1
     integer::countsend=0,countrecv=0
     integer,allocatable,dimension(:)::reqsend,reqrecv,roff
     real(kind=8),allocatable,dimension(:)::sendbuf,recvbuf
  end type async_exchange
  type(async_exchange),dimension(1:nasync)::async_rev
  type(async_exchange),dimension(1:nasync)::async_fwd   ! owner -> virtual (make_virtual_fine_dp_start/end)
  ! Grids that changed owner to this rank at the last load balance, per level (timer_ranks diagnostics),
  ! and the running total over the run
  integer,allocatable,dimension(:)::nrecv_lb
  integer::nrecv_lb_tot=0
  ! Time-based correction of the load balance cost (cost_time_weighting): factor applied to the cost of
  ! every cell owned by this rank, the compute-timer total at the previous load balance and whether the
  ! interval since then is a clean measurement (a full set of steps on this partition)
  real(dp)::cost_rank_fac=1.0d0
  ! Supernovae injected by this rank since the start of the run (timer_ranks diagnostics)
  integer::nsn_rank=0
  real(kind=8)::cost_time_prev=0d0
  logical::cost_time_valid=.false.

  ! Active grid, emission and reception communicators
  type(communicator),allocatable,dimension(:)  ::active               ! 1D (1:nlevelmax) array => leave as is
  type(communicator),allocatable,dimension(:,:)::boundary             ! 2D (1:MAXBOUND, 1:nlevelmax) array => leave as is
#ifdef LIGHT_MPI_COMM
  type(communicator_light),  allocatable,dimension(:,:)::reception       ! 2D (ncpu,nlevelmax) light reception buffer
  type(communicator_varoct), allocatable,dimension(:)  ::emission        ! 1D (nlevelmax) AMR data emission buffer
  type(communicator_varpart),allocatable,dimension(:)  ::emission_part   ! 1D (nlevelmax) point data emission buffer
! ---------------------------------------------------------------------------------------------------------------------------- !
#else
  type(communicator),allocatable,dimension(:,:)::emission    ! 2D (ncpu,nlevelmax) data emission "heavy" buffer
  type(communicator),allocatable,dimension(:,:)::reception   ! 2D (ncpu, nlevelmax) data reception "heavy" buffer
#endif


  ! Types for physical boundary conditions
  CHARACTER(LEN=20)::type_hydro  ='hydro'
  CHARACTER(LEN=20)::type_accel  ='accel'
  CHARACTER(LEN=20)::type_flag   ='flag'

  ! Units specified by the user in the UNITS_PARAMS namelist for non-cosmo runs.
  ! These values shouldn't be used directly. Instead call units() in amr/units.f90.
  real(dp)::units_density=1.0d0  ! [g/cm^3]
  real(dp)::units_time=1.0d0     ! [seconds]
  real(dp)::units_length=1.0d0   ! [cm]

contains
  ! two-tier domain decomposition helpers (ported from patch/ks).
  ! Fine-tier domain idom belongs to cpu two_tier_fine_cpu(idom); with
  ! two_tier_test_shift=0 this is the usual mod(idom-1,ncpu)+1.
  pure integer function two_tier_fine_cpu(idom)
    integer,intent(in)::idom
    two_tier_fine_cpu=modulo(idom-1+two_tier_test_shift,ncpu)+1
  end function two_tier_fine_cpu
  ! Fine-tier domain index of cpu icpu (overload sub-domain isub)
  pure integer function two_tier_fine_dom(icpu,isub)
    integer,intent(in)::icpu,isub
    two_tier_fine_dom=modulo(icpu-1-two_tier_test_shift,ncpu)+1+(isub-1)*ncpu
  end function two_tier_fine_dom
  ! Is key inside a footprint interval of domain idom?
  pure logical function two_tier_in_foot(key,idom,nfoot,lo,hi,first)
    real(qdp),intent(in)::key
    integer,intent(in)::idom,nfoot
    real(qdp),intent(in),dimension(1:max(nfoot,1))::lo,hi
    integer,intent(in),dimension(0:ndomain)::first
    integer::jlo,jhi,j
    two_tier_in_foot=.false.
    jlo=first(idom-1)+1; jhi=first(idom)
    if(jhi<jlo)return
    do while(jhi>jlo)          ! last interval of the domain with lo<=key
      j=(jlo+jhi+1)/2
      if(lo(j)<=key)then
        jlo=j
      else
        jhi=j-1
      end if
    end do
    if(lo(jlo)<=key.and.key<hi(jlo))two_tier_in_foot=.true.
  end function two_tier_in_foot
  ! Does the key range [kmin,kmax) overlap a footprint interval of domain idom?
  pure logical function two_tier_foot_overlap(kmin,kmax,idom,nfoot,lo,hi,first)
    real(qdp),intent(in)::kmin,kmax
    integer,intent(in)::idom,nfoot
    real(qdp),intent(in),dimension(1:max(nfoot,1))::lo,hi
    integer,intent(in),dimension(0:ndomain)::first
    integer::jlo,jhi,j
    two_tier_foot_overlap=.false.
    jlo=first(idom-1)+1; jhi=first(idom)
    if(jhi<jlo)return
    if(lo(jlo)>=kmax)return
    do while(jhi>jlo)          ! last interval of the domain with lo<kmax
      j=(jlo+jhi+1)/2
      if(lo(j)<kmax)then
        jlo=j
      else
        jhi=j-1
      end if
    end do
    if(hi(jlo)>kmin)two_tier_foot_overlap=.true.
  end function two_tier_foot_overlap
  ! Is [kmin,kmax) inside the footprints (one interval, or touching intervals of
  ! neighbouring domains) of the whole list? False means it contains a gap.
  pure logical function two_tier_in_one_foot(kmin,kmax,nfoot,lo,hi)
    real(qdp),intent(in)::kmin,kmax
    integer,intent(in)::nfoot
    real(qdp),intent(in),dimension(1:max(nfoot,1))::lo,hi
    integer::jlo,jhi,j
    real(qdp)::h
    two_tier_in_one_foot=.false.
    if(nfoot==0)return
    if(lo(1)>kmin)return
    jlo=1; jhi=nfoot
    do while(jhi>jlo)          ! last interval with lo<=kmin
      j=(jlo+jhi+1)/2
      if(lo(j)<=kmin)then
        jlo=j
      else
        jhi=j-1
      end if
    end do
    h=hi(jlo); j=jlo
    do while(j<nfoot)
      if(lo(j+1)/=h)exit
      j=j+1; h=hi(j)
    end do
    two_tier_in_one_foot=(kmax<=h)
  end function two_tier_in_one_foot
end module amr_commons

