module pm_commons

  use amr_parameters
  use pm_parameters
  use random

  implicit none

  ! Sink particle related arrays
  real(dp),allocatable,dimension(:)    ::msink,xmsink
  real(dp),allocatable,dimension(:)    ::msink_new,msink_all
  real(dp),allocatable,dimension(:)    ::msmbh,msmbh_new,msmbh_all
  real(dp),allocatable,dimension(:)    ::dmfsink,dmfsink_new,dmfsink_all !count mass in sink since last stellar object creation
  real(dp),allocatable,dimension(:)    ::oksink_new,oksink_all
  real(dp),allocatable,dimension(:)    ::tsink,tsink_new,tsink_all
  real(dp),allocatable,dimension(:)    ::dMsink_overdt,dMBHoverdt
  real(dp),allocatable,dimension(:)    ::dMsmbh_overdt,dMBHoverdt_smbh
  real(dp),allocatable,dimension(:)    ::rho_gas,volume_gas,eps_sink,c2sink
  real(dp),allocatable,dimension(:,:)  ::vel_gas
  real(dp),allocatable,dimension(:)    ::delta_mass,delta_mass_new,delta_mass_all
  real(dp),allocatable,dimension(:)    ::wden,weth,wvol,wdiv,wden_new,weth_new,wvol_new,wdiv_new
  real(dp),allocatable,dimension(:,:)  ::wmom,wmom_new
  real(dp),allocatable,dimension(:,:)  ::vsink,vsink_new,vsink_all
  real(dp),allocatable,dimension(:,:)  ::fsink,fsink_new,fsink_all
  real(dp),allocatable,dimension(:,:,:)::vsnew,vsold
  real(dp),allocatable,dimension(:,:,:)::fsink_partial,sink_jump
  real(dp),allocatable,dimension(:,:)  ::lsink,lsink_new,lsink_all
  real(dp),allocatable,dimension(:,:)  ::xsink,xsink_new,xsink_all
  real(dp),allocatable,dimension(:)    ::graddescent_over_dt
  real(dp),allocatable,dimension(:,:)  ::xsink_graddescent
  real(dp),allocatable,dimension(:,:)  ::weighted_density,weighted_volume,weighted_ethermal,weighted_divergence
  real(dp),allocatable,dimension(:,:,:)::weighted_momentum
  real(dp),allocatable,dimension(:)    ::rho_sink_tff
  real(dp),allocatable,dimension(:)    ::msum_overlap
  integer,allocatable,dimension(:)     ::idsink,idsink_new,idsink_old,idsink_all
  logical,allocatable,dimension(:)     ::ok_blast_agn,ok_blast_agn_all
  logical,allocatable,dimension(:)     ::direct_force_sink
  logical,allocatable,dimension(:)     ::new_born,new_born_all,new_born_new
  integer,allocatable,dimension(:)     ::idsink_sort
  integer::ncloud_sink,ncloud_sink_massive
  integer::nindsink=0
  integer::sinkint_level=0         ! maximum level currently active is where the global sink variables are updated
  real(dp)::ssoft                  ! sink softening lenght in code units

  ! Particles related arrays
  real(dp),allocatable,dimension(:,:)  ::xp       ! Positions
  real(dp),allocatable,dimension(:,:)  ::vp       ! Velocities
  real(dp),allocatable,dimension(:)    ::mp       ! Masses
  integer,  allocatable, dimension(:)  :: move_flag ! Move flag (for particles), >0 means don't move!
#ifdef OUTPUT_PARTICLE_POTENTIAL
  real(dp),allocatable,dimension(:)    ::ptcl_phi ! Potential of particle added by AP for output purposes
#endif
  real(dp),allocatable,dimension(:)    ::tp       ! Birth epoch
  real(dp),allocatable,dimension(:)    ::zp       ! Birth metallicity
#ifdef INIT_STELLAR_MASS
  real(dp),allocatable,dimension(:)    ::mp0       ! Birth time mass
#endif
#ifdef STELLAR_POPULATION_MASS
  real(dp),allocatable,dimension(:)    ::msp0       ! The mass of stellar population at birth (0 or 100 Msun)
#endif
  integer,  allocatable, dimension(:)  :: itmpp    ! Working array
  integer,  allocatable, dimension(:)  :: partp    ! Particle parent (for tracers only)
  integer ,allocatable,dimension(:)    ::nextp    ! Next particle in list
  integer ,allocatable,dimension(:)    ::prevp    ! Previous particle in list
  integer ,allocatable,dimension(:)    ::levelp   ! Current level of particle
  integer(i8b),allocatable,dimension(:)::idp    ! Identity of particle
  ! Tree related arrays
  integer ,allocatable,dimension(:)    ::headp    ! Head particle in grid
  integer ,allocatable,dimension(:)    ::tailp    ! Tail particle in grid
  integer ,allocatable,dimension(:)    ::numbp    ! Number of particles in grid
  ! Global particle linked lists
  integer::headp_free,tailp_free,numbp_free=0,numbp_free_tot=0
  ! Local and current seed for random number generator
  integer,dimension(IRandNumSize) :: localseed=-1
  integer, dimension(IRandNumSize) :: tracer_seed = -1

  ! Particle types
  ! Two interpretations of the family byte, selected at run time by sf_fg (&SF_PARAMS):
  !   sf_fg=.false. (upstream): 1 DM, 2 star, 3 sink cloud, 4 debris (f_w>0), 5 other (DICE star_as_nbody)
  !   sf_fg=.true.  (patch/ks): 1 DM, 2 Pop II star, 3 Pop III, 4 SN (just exploded), 5 dead Pop II,
  !                             6 BH (Pop III remnant), 7 SFC (active star-forming cloud, zp=-radius),
  !                             8 PSC (SFC past SFC_t_life)
  ! The numbers 3-5 coincide on purpose (ks outputs can be restarted); is_cloud/is_debris and
  ! is_pop3/is_SN/... must therefore only be used in the matching mode (docs/merge_ks.md).
  integer, parameter   :: NFAMILIES=8
  integer(1),parameter :: FAM_DM=1, FAM_STAR=2, FAM_CLOUD=3, FAM_DEBRIS=4, FAM_OTHER=5, FAM_UNDEF=127
  integer(1),parameter :: FAM_POP3=3, FAM_SN=4, FAM_DEAD=5, FAM_BH=6, FAM_SFC=7, FAM_PSC=8   ! sf_fg=.true.
  integer(1),parameter :: FAM_TRACER_GAS=0
  integer(1),parameter :: FAM_TRACER_DM=-1, FAM_TRACER_STAR=-2, FAM_TRACER_CLOUD=-3, FAM_TRACER_DEBRIS=-4, FAM_TRACER_OTHER=-5

  ! SFC (star-forming cloud) particles of the FG star formation (sf_fg): list of the active ones
  ! (collect_active_SFC in star_formation_fg.f90)
  integer,parameter :: SFC_nmax_active = 10000 !maximum number of active SFC
  integer :: SFC_n_active = 0 !number of active SFC
  real(dp),dimension(1:3,1:SFC_nmax_active) :: SFC_pos_active !position of active SFC
  real(dp),dimension(1:SFC_nmax_active) :: SFC_radius_active !radius of active SFC

  ! SN events of the current step on this rank (feedback.f90 -> SN reservoir / radius search; merge_ks)
  integer,parameter :: SNinfo_nSNmax = 10000      !maximum number of SNe per rank and step
  integer :: SNinfo_nSN = 0                       !number recorded on this rank
  integer,dimension(1:SNinfo_nSNmax) :: SNinfo_idp
  real(dp),dimension(1:SNinfo_nSNmax) :: SNinfo_ESN    !pure SN energy [code units] (no winds)
  real(dp),dimension(1:3,1:SNinfo_nSNmax) :: SNinfo_pos

  ! Customize here for particle tags within particle types (e.g. different kind of stars).
  ! Note that the type should be integer(1) (1 byte integers) for memory concerns.
  ! Also don't forget to create a function is_<type>_<tag>. See the wiki for a more complete example.
  ! By default, the tag is always 0.

  ! Particle keys for outputing. They should match the above particle
  ! types, except for 'under' family
  ! (-8..-6 and 6..8 are unused in the upstream mode; particle_family_keys_fg is the sf_fg table)
  character(len=13), dimension(-NFAMILIES:NFAMILIES), parameter :: particle_family_keys = (/ &
    '       unused', '       unused', '       unused', &
    ' other_tracer', 'debris_tracer', ' cloud_tracer', '  star_tracer', ' other_tracer', &
    '   gas_tracer', &
    '           DM', '         star', '        cloud', '       debris', '        other', &
    '       unused', '       unused', '       unused'/)
  character(len=13), dimension(-NFAMILIES:NFAMILIES), parameter :: particle_family_keys_fg = (/ &
    '       unused', '       unused', '       unused', &
    ' other_tracer', 'debris_tracer', ' cloud_tracer', '  star_tracer', ' other_tracer', &
    '   gas_tracer', &
    '           DM', '         star', '         pop3', '    supernova', '    dead_pop2', &
    '      bh_pop3', '          SFC', '          PSC'/)   ! no blanks in a key (Mera/yt parse the header table)

  type(part_t), allocatable, dimension(:) :: typep  ! Particle type array

contains
  function cross(a,b)
    use amr_parameters, only:dp
    real(dp),dimension(1:3)::a,b
    real(dp),dimension(1:3)::cross
    !computes the cross product c= a x b
    cross(1)=a(2)*b(3)-a(3)*b(2)
    cross(2)=a(3)*b(1)-a(1)*b(3)
    cross(3)=a(1)*b(2)-a(2)*b(1)
  end function cross

  elemental logical pure function is_DM(typep)
    type(part_t), intent(in) :: typep
    is_DM = typep%family == FAM_DM
  end function is_DM

  elemental logical pure function is_star(typep)
    type(part_t), intent(in) :: typep
    is_star = typep%family == FAM_STAR
  end function is_star

  ! sf_fg=.true. families (patch/ks); meaningless when sf_fg=.false. (3-5 are cloud/debris/other there)
  elemental logical pure function is_pop3(typep)
    type(part_t), intent(in) :: typep
    is_pop3 = typep%family == FAM_POP3
  end function is_pop3

  elemental logical pure function is_SN(typep)
    type(part_t), intent(in) :: typep
    is_SN = typep%family == FAM_SN
  end function is_SN

  elemental logical pure function is_DEAD(typep)
    type(part_t), intent(in) :: typep
    is_DEAD = typep%family == FAM_DEAD
  end function is_DEAD

  elemental logical pure function is_BH(typep)
    type(part_t), intent(in) :: typep
    is_BH = typep%family == FAM_BH
  end function is_BH

  elemental logical pure function is_SFC(typep)
    type(part_t), intent(in) :: typep
    is_SFC = typep%family == FAM_SFC
  end function is_SFC

  elemental logical pure function is_PSC(typep)
    type(part_t), intent(in) :: typep
    is_PSC = typep%family == FAM_PSC
  end function is_PSC

  elemental logical pure function is_cloud(typep)
    type(part_t), intent(in) :: typep
    is_cloud = typep%family == FAM_CLOUD
  end function is_cloud

  elemental logical pure function is_debris(typep)
    type(part_t), intent(in) :: typep
    is_debris = typep%family == FAM_DEBRIS
  end function is_debris

  elemental logical pure function is_tracer(typep)
    type(part_t), intent(in) :: typep
    is_tracer = typep%family <= 0
  end function is_tracer

  elemental logical pure function is_not_tracer(typep)
    type(part_t), intent(in) :: typep
    is_not_tracer = typep%family > 0
  end function is_not_tracer

  elemental logical pure function is_gas_tracer(typep)
    type(part_t), intent(in) :: typep
    is_gas_tracer = typep%family == FAM_TRACER_GAS
  end function is_gas_tracer

  elemental logical pure function is_not_DM(typep)
    ! Check that the particle is not DM and not a tracer
    type(part_t), intent(in) :: typep
    is_not_DM = typep%family /= FAM_DM .and. is_not_tracer(typep)
  end function is_not_DM

  elemental logical pure function is_star_tracer(typep)
    type(part_t), intent(in) :: typep
    is_star_tracer = typep%family == FAM_TRACER_STAR
  end function is_star_tracer

  elemental logical pure function is_cloud_tracer(typep)
    type(part_t), intent(in) :: typep
    is_cloud_tracer = typep%family == FAM_TRACER_CLOUD
  end function is_cloud_tracer

  elemental function part2int (part)
    ! Convert a particle into an integer
    ! This saves some space e.g. when communicating
    integer :: part2int
    type(part_t), intent(in) :: part

    ! This is the largest value for integer(1)
    integer, parameter :: a = 128, b = 2*a

    part2int = (int(part%family) + a) * b + (int(part%tag) + a)
  end function part2int

  elemental function int2part(index)
    ! Convert from an index to particle type
    type(part_t) :: int2part
    integer, intent(in) :: index

    ! This is the largest value for integer(1)
    integer, parameter :: a = 128, b = 2*a

    int2part%family = int(index / b - a, 1)
    int2part%tag = int(mod(index, b) - a, 1)
  end function int2part

  function props2type(idpii, tpii, mpii)
    use amr_commons
    use pm_parameters, only : part_t

    ! Converts from "old" ramses to "new" ramses
    !
    ! Here's the match, add yours here for backward compatibility purposes
    ! DM     tpii == 0
    ! stars  tpii != 0 and idpii > 0
    ! sinks  tpii != 0 and idpii < 0
    !
    ! This is mostly for support of GRAFFIC I/O.
    ! The reason we use idpii instead of idp is to prevent name clashes
    real(dp), intent(in) :: tpii, mpii
    integer, intent(in)  :: idpii

    type(part_t) :: props2type

    if (tpii == 0) then
      props2type%family = FAM_DM
    else if (idpii > 0) then
      props2type%family = FAM_STAR
    else if (idpii < 0) then
      props2type%family = FAM_CLOUD
    else if (mpii == 0) then
      props2type%family = FAM_TRACER_GAS
    end if
    props2type%tag = 0
  end function props2type
end module pm_commons
