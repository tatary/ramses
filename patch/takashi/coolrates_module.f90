MODULE coolrates_module
  ! Module for returning cooling and thermochemistry interaction rates.
  ! The temperature dependence is tabulated in rate (nonlog) versus the
  ! log of T, because it is expensive to calculate on the fly.
  ! The rates are interpolated using cubic splines, and extrapolated in
  ! log-log space if temperature is above table boundaries.
  ! Joki Rosdahl and Andreas Bleuler, September 2015.
  ! Nickerson, Teyssier, and Rosdahl, March 2018.

  use amr_parameters,only:dp,H2coolfunc,H2_fesc,H2_Hm_cool,H2_chem_cool
  use constants,only:twopi,kB,mH,eV2erg,factG_in_cgs
  implicit none

  private   ! default
  public init_coolrates_tables, update_coolrates_tables                  &
    , inp_coolrates_table, compCoolrate, coolrates_lookup             &
    , comp_Sd, comp_SH2_G09, comp_SH2_DB96, comp_SH2_WGH              &
    , comp_Alpha_H2_Katz, comp_Beta_H2coll, T_H2chem_max              &
    , tbl_alphaZ_H2, tbl_alphaGP_H2, tbl_alphaA_HII, tbl_alphaA_HeII  &
    , tbl_alphaA_HeIII, tbl_alphaB_HII, tbl_alphaB_HeII               &
    , tbl_alphaB_HeIII, tbl_beta_HI, tbl_beta_HeI, tbl_beta_HeII      &
    , tbl_cr_ci_HI, tbl_cr_ci_HeI, tbl_cr_ci_HeII, tbl_cr_ce_HI       &
    , tbl_cr_ce_HeI, tbl_cr_ce_HeII, tbl_cr_r_HII, tbl_cr_r_HeII      &
    , tbl_cr_r_HeIII, tbl_cr_bre, tbl_cr_com, tbl_cr_die              &
    , tbl_beta_H2HI, tbl_beta_H2H2, tbl_cr_H2HI, tbl_cr_H2H2          &
    , tbl_beta_H3B

  ! H2 cooling options (ported from patch/ks 2026-11-12, see docs/merge_ks.md):
  ! the H2/H- fits are evaluated at min(T,T_H2chem_fitmax) and the H2/H-
  ! channels are skipped above T_H2chem_max (H2 is fully dissociated there)
  real(dp),parameter    :: T_H2chem_fitmax = 1d5
  real(dp),parameter    :: T_H2chem_max    = 1d6

  ! Default cooling rates table parameters
  integer,parameter     :: nbinT  = 1001
  real(dp),parameter    :: Tmin   = 1d-2
  real(dp),parameter    :: Tmax   = 1d+9
  real(dp)              :: dlogTinv ! Inverse of the bin space (in K)
  real(dp)              :: hTable, h2Table, h3Table   ! Interpol constants
  real(dp)              :: one_over_lnTen, one_over_hTable, one_over_h2Table
  real(dp)              :: three_over_h2Table, two_over_h3Table

  real(dp),dimension(nbinT) :: T_lookup = 0d0 ! Lookup temperature in log K

  type coolrates_table
    ! Cooling and interaction rates (log):
    real(dp),dimension(nbinT)::rates  = 0d0
    ! Temperature derivatives of those rates (drate/dlog(T)):
    real(dp),dimension(nbinT)::primes = 0d0
  end type coolrates_table

  ! Cache of the last T lookup for inp_coolrates_table (index and spline
  ! weights, one log10 per new T). Declared by the caller and passed as the
  ! optional argument "lookup", so that repeated lookups of different tables
  ! at the same T (e.g. inside compCoolrate) share the work without any
  ! save/threadprivate state (ifx threadprivate access doubled the cooling
  ! time, TODO.md 2d-3).
  type coolrates_lookup
     real(dp)::Tlast=-1d0
     integer ::iT=1
     real(dp)::facT=0d0, yy=0d0, yy2=0d0, yy3=0d0
     logical ::extrap=.false.
  end type coolrates_lookup

  type(coolrates_table),save::tbl_alphaZ_H2  ! H2 dust grain formation
  type(coolrates_table),save::tbl_alphaGP_H2 ! H2 gas-phase formation
  type(coolrates_table),save::tbl_alphaA_HII ! Case A rec. coefficients
  type(coolrates_table),save::tbl_alphaA_HeII
  type(coolrates_table),save::tbl_alphaA_HeIII
  type(coolrates_table),save::tbl_alphaB_HII ! Case B rec. coefficients
  type(coolrates_table),save::tbl_alphaB_HeII
  type(coolrates_table),save::tbl_alphaB_HeIII
  type(coolrates_table),save::tbl_beta_HI ! Collisional ionisation rates
  type(coolrates_table),save::tbl_beta_HeI
  type(coolrates_table),save::tbl_beta_HeII
  type(coolrates_table),save::tbl_cr_ci_HI ! Coll. ionisation cooling
  type(coolrates_table),save::tbl_cr_ci_HeI
  type(coolrates_table),save::tbl_cr_ci_HeII
  type(coolrates_table),save::tbl_cr_ce_HI ! Coll. excitation cooling
  type(coolrates_table),save::tbl_cr_ce_HeI
  type(coolrates_table),save::tbl_cr_ce_HeII
  type(coolrates_table),save::tbl_cr_r_HII ! Recombination cooling
  type(coolrates_table),save::tbl_cr_r_HeII
  type(coolrates_table),save::tbl_cr_r_HeIII
  type(coolrates_table),save::tbl_cr_bre ! Bremsstrahlung cooling rates
  type(coolrates_table),save::tbl_cr_com ! Compton cooling rates
  type(coolrates_table),save::tbl_cr_die ! Dielectronic cooling rates
  type(coolrates_table),save::tbl_beta_H2HI ! Coll. dissociation rates
  type(coolrates_table),save::tbl_beta_H2H2 ! Coll. dissociation rates
  type(coolrates_table),save::tbl_beta_H3B ! Coll. dissociation rates
  type(coolrates_table),save::tbl_cr_H2HI ! Collisional diss. cooling
  type(coolrates_table),save::tbl_cr_H2H2 ! Collisional diss. cooling

CONTAINS

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  SUBROUTINE init_coolrates_tables(aexp)

! Initialise the cooling rates tables.
!-------------------------------------------------------------------------
    use mpi_mod
    implicit none
#ifndef WITHOUTMPI
    integer :: ierr
#endif
    real(dp) :: aexp
    integer :: myid, ncpu, iT
!-------------------------------------------------------------------------
#ifndef WITHOUTMPI
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)
#endif
#ifdef WITHOUTMPI
    myid=0
    ncpu=1
#endif

    ! Initialise the table lookup temperatures -----------------------------
    do iT=1, nbinT
      T_lookup(iT) = log10(Tmin) + (dble(iT)-1d0) / (dble(nbinT)-1d0)     &
        * (log10(Tmax)-log10(Tmin))
    end do
    dlogTinv = dble(nbinT-1)/(T_lookup(nbinT)-T_lookup(1)) ! (space)^-1
    hTable = 1d0/dlogTinv                             !
    h2Table = hTable*hTable                           ! Constants for table
    h3Table = h2Table*hTable                          ! interpolation

    one_over_lnTen = 1d0/log(10d0)
    one_over_hTable = 1d0/hTable
    one_over_h2Table = 1d0/h2Table
    three_over_h2Table = 3d0/h2Table
    two_over_h3Table = 2d0/h3Table

    do iT = myid+1, nbinT, ncpu ! Loop over TK and assign rates
      call comp_table_rates(iT,aexp)
    end do ! end TK loop

    ! Distribute the complete table between cpus ---------------------------
#ifndef WITHOUTMPI
    call mpi_distribute_coolrates_table(tbl_alphaZ_H2)
    call mpi_distribute_coolrates_table(tbl_alphaGP_H2)
    call mpi_distribute_coolrates_table(tbl_alphaA_HII)
    call mpi_distribute_coolrates_table(tbl_alphaA_HeII)
    call mpi_distribute_coolrates_table(tbl_alphaA_HeIII)

    call mpi_distribute_coolrates_table(tbl_alphaB_HII)
    call mpi_distribute_coolrates_table(tbl_alphaB_HeII)
    call mpi_distribute_coolrates_table(tbl_alphaB_HeIII)

    call mpi_distribute_coolrates_table(tbl_beta_HI)
    call mpi_distribute_coolrates_table(tbl_beta_HeI)
    call mpi_distribute_coolrates_table(tbl_beta_HeII)

    call mpi_distribute_coolrates_table(tbl_cr_ci_HI)
    call mpi_distribute_coolrates_table(tbl_cr_ci_HeI)
    call mpi_distribute_coolrates_table(tbl_cr_ci_HeII)

    call mpi_distribute_coolrates_table(tbl_cr_ce_HI)
    call mpi_distribute_coolrates_table(tbl_cr_ce_HeI)
    call mpi_distribute_coolrates_table(tbl_cr_ce_HeII)

    call mpi_distribute_coolrates_table(tbl_cr_r_HII)
    call mpi_distribute_coolrates_table(tbl_cr_r_HeII)
    call mpi_distribute_coolrates_table(tbl_cr_r_HeIII)

    call mpi_distribute_coolrates_table(tbl_cr_bre)
    call mpi_distribute_coolrates_table(tbl_cr_com)
    call mpi_distribute_coolrates_table(tbl_cr_die)

    call mpi_distribute_coolrates_table(tbl_beta_H2HI)
    call mpi_distribute_coolrates_table(tbl_beta_H2H2)
    call mpi_distribute_coolrates_table(tbl_beta_H3B)
    call mpi_distribute_coolrates_table(tbl_cr_H2HI)
    call mpi_distribute_coolrates_table(tbl_cr_H2H2)

#endif

    if(myid==0) print*,'Coolrates tables initialised '

  end subroutine init_coolrates_tables

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  SUBROUTINE update_coolrates_tables(aexp)
! Update cooling rates lookup tables which depend on aexp
!-------------------------------------------------------------------------
    use mpi_mod
    implicit none
#ifndef WITHOUTMPI
    integer :: ierr
#endif
    real(dp) :: aexp
    integer:: myid, ncpu, iT
!-------------------------------------------------------------------------
#ifndef WITHOUTMPI
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)
#endif
#ifdef WITHOUTMPI
    myid=0
    ncpu=1
#endif
    tbl_cr_com%rates  = 0d0 ; tbl_cr_com%primes = 0d0
    do iT = myid+1, nbinT, ncpu ! Loop over TK and assign rates
      call update_table_rates(iT, aexp)
    end do ! end TK loop

    ! Distribute the complete table between cpus ---------------------------
#ifndef WITHOUTMPI
    call mpi_distribute_coolrates_table(tbl_cr_com)
#endif

    if(myid==0) print*,'Coolrates table updated'
  end subroutine update_coolrates_tables

#ifndef WITHOUTMPI
!PRIVATEXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  SUBROUTINE mpi_distribute_coolrates_table(table)
! Distribute table between all cpus, assuming table contains only partial
! entries on each cpu, but the whole table is acquired by summing those
! partial tables
!-------------------------------------------------------------------------
    use mpi_mod
    implicit none
    type(coolrates_table)::table
    real(dp),dimension(:),allocatable :: table_mpi_sum
    integer::ierr
!-------------------------------------------------------------------------
    allocate(table_mpi_sum(nbinT))
    call MPI_ALLREDUCE(table%rates,table_mpi_sum                           &
      ,nbinT,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    table%rates = table_mpi_sum
    call MPI_ALLREDUCE(table%primes,table_mpi_sum                          &
      ,nbinT,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    table%primes = table_mpi_sum
    deallocate(table_mpi_sum)

  end subroutine mpi_distribute_coolrates_table
#endif

!PRIVATEXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  SUBROUTINE comp_table_rates(iT, aexp)
! Fill in index iTK in all rates tables.
!-------------------------------------------------------------------------
    use rt_parameters,only:rt_OTSA
    use constants,only:kB
    implicit none
    integer::iT
    real(dp)::aexp, T, T2, Ta, T5, lambda, f, hf, laHII, laHeII, laHeIII
    real(dp)::lowrleft,lowrright,lowr_hi,lowr_h2
    real(dp)::lowvleft_hi,lowvright_hi,lowv_hi,lowv_h2
    real(dp)::TT
!-------------------------------------------------------------------------
    ! Rates are stored in non-log, while temperature derivatives (primes)
    ! are stored in dRate/dlogT (= dRate/dT * T * ln(10)).

    ! The log-log primes are just the normal primes times T/rate,
    ! i.e. dlogL/dlogT = T/L dL/dT

    T = 10d0**T_lookup(iT)

    ! Creation rate of H2 on dust [cm3 s-1] (functio of T from--------------
    ! Hollenback & McKee 1979, average rate between Habart 2004 dense PDRs
    ! and Jura 1974 diffuse ISM)
    ! When using, must multiply by Z*f_dust
    T2 = T/1d2
    f = (1.0+0.4*sqrt(T2)+0.2*T2+0.08*T2**2)
    tbl_alphaZ_H2%rates(iT)   = 9.0d-17 * sqrt(T2) / f
    tbl_alphaZ_H2%primes(iT)  = (-0.12*T2**2 - 0.1*T2 + 0.5) / f           &
      * log(10d0) * tbl_alphaZ_H2%rates(iT)

    ! Gas phase formation rate rate for H2 on H- [cm3 s-1] assuming---------
    ! equilibrium abundances for H-, as explained in the Appendix of
    ! McKee and Krumholz (2012)
    tbl_alphaGP_H2%rates(iT)  = 8.0d-19*(T/1000.0)**0.88
    tbl_alphaGP_H2%primes(iT) = 0.88                                       &
      * log(10d0) * tbl_alphaGP_H2%rates(iT)
    ! Case A rec. coefficient [cm3 s-1] for HII (Hui&Gnedin'97)-------------
    lambda = 315614./T                                ! 2d0 * 157807d0 / T
    f = 1d0+(lambda/0.522)**0.47
    tbl_alphaA_HII%rates(iT)  =  1.269d-13 * lambda**1.503 / f**1.923
    tbl_alphaA_HII%primes(iT) = ( 0.90381*(f-1.)/f - 1.503 )               &
      * log(10d0) * tbl_alphaA_HII%rates(iT)

    ! Case A rec. coefficient [cm3 s-1] for HeII (Hui&Gnedin'97)------------
    lambda = 570670./T
    tbl_alphaA_HeII%rates(iT)  = 3d-14 * lambda**0.654
    tbl_alphaA_HeII%primes(iT) = -0.654                                    &
      * log(10d0) * tbl_alphaA_HeII%rates(iT)

    ! Case A rec. coefficient [cm3 s-1] for HeIII (Hui&Gnedin'97)-----------
    lambda =  1263030./T
    f= 1d0+(lambda/0.522)**0.47
    tbl_alphaA_HeIII%rates(iT)  =  2.538d-13 * lambda**1.503 / f**1.923
    tbl_alphaA_HeIII%primes(iT) =  ( 0.90381*(f-1.)/f - 1.503 )            &
      * log(10d0) * tbl_alphaA_HeIII%rates(iT)

    ! Case B rec. coefficient [cm3 s-1] for HII (Hui&Gnedin'97)-------------
    lambda = 315614./T
    f= 1d0+(lambda/2.74)**0.407
    tbl_alphaB_HII%rates(iT)  = 2.753d-14 * lambda**1.5 / f**2.242
    tbl_alphaB_HII%primes(iT) = ( 0.912494*(f-1.)/f - 1.5 )                &
      * log(10d0) * tbl_alphaB_HII%rates(iT)

    ! Case B rec. coefficient [cm3 s-1] for HeII (Hui&Gnedin'97)------------
    lambda = 570670./T
    tbl_alphaB_HeII%rates(iT)  = 1.26d-14 * lambda**0.75
    tbl_alphaB_HeII%primes(iT) =  -0.75                                    &
      * log(10d0) * tbl_alphaB_HeII%rates(iT)

    ! Case B rec. coefficient [cm3 s-1] for HeIII (Hui&Gnedin'97)-----------
    lambda = 1263030./T
    f= 1d0+(lambda/2.74)**0.407
    tbl_alphaB_HeIII%rates(iT)  = 5.506d-14 * lambda**1.5 / f**2.242
    tbl_alphaB_HeIII%primes(iT) = ( 0.912494*(f-1.)/f - 1.5 )              &
      * log(10d0) * tbl_alphaB_HeIII%rates(iT)

    ! Collisional ionization rate [cm3 s-1] of HI (Maselli&'03)-------------
    T5 = T/1d5
    f = 1d0+sqrt(T5) ; hf=0.5d0/f
    tbl_beta_HI%rates(iT)  = 5.85d-11 * sqrt(T) / f * exp(-157809.1d0/T)
    tbl_beta_HI%primes(iT) = (hf+157809.1d0/T)                             &
      * log(10d0) * tbl_beta_HI%rates(iT)

    ! Collisional ionization rate [cm3 s-1] of HeI (Maselli&'03)------------
    tbl_beta_HeI%rates(iT)  = 2.38d-11 * sqrt(T) / f * exp(-285335.4d0/T)
    tbl_beta_HeI%primes(iT) = (hf+285335.4d0/T)                            &
      * log(10d0) * tbl_beta_HeI%rates(iT)

    ! Collisional ionization rate [cm3 s-1] of HeII (Maselli&'03)-----------
    tbl_beta_HeII%rates(iT)  = 5.68d-12 * sqrt(T) / f * exp(-631515d0/T)
    tbl_beta_HeII%primes(iT) = (hf+631515d0/T)                            &
      * log(10d0) * tbl_beta_HeII%rates(iT)

    ! BEGIN COOLING RATES-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    T5 = T/1d5
    f = 1d0+sqrt(T5) ; hf=0.5d0/f

    ! Coll. Ionization Cooling from Cen 1992 (via Maselli et al 2003)
    tbl_cr_ci_HI%rates(iT)  = 1.27d-21 * sqrt(T) / f * exp(-157809.1/T)
    tbl_cr_ci_HI%primes(iT) = (hf+157809.1/T)                              &
      * log(10d0) * tbl_cr_ci_HI%rates(iT)

    tbl_cr_ci_HeI%rates(iT)  = 9.38d-22 * sqrt(T) / f * exp(-285335.4/T)
    tbl_cr_ci_HeI%primes(iT) = (hf+285335.4/T)                             &
      * log(10d0) * tbl_cr_ci_HeI%rates(iT)

    tbl_cr_ci_HeII%rates(iT)  = 4.95d-22 * sqrt(T) / f * exp(-631515. /T)
    tbl_cr_ci_HeII%primes(iT) = (hf+631515.0/T)                            &
      * log(10d0) * tbl_cr_ci_HeII%rates(iT)

    ! Collisional excitation cooling from Cen'92
    tbl_cr_ce_HI%rates(iT)  = 7.5d-19 / f * exp(-118348./T)
    tbl_cr_ce_HI%primes(iT) = (118348./T - 0.5d0 * sqrt(T5) / f )          &
      * log(10d0) * tbl_cr_ce_HI%rates(iT)

    tbl_cr_ce_HeI%rates(iT)  = 9.10d-27 * T**(-0.1687) / f * exp(-13179./T)
    tbl_cr_ce_HeI%primes(iT) = (13179./T - 0.1687 - 0.5d0 * sqrt(T5) / f)  &
      * log(10d0) * tbl_cr_ce_HeI%rates(iT)

    tbl_cr_ce_HeII%rates(iT) = 5.54d-17 * T**(-0.397)  / f * exp(-473638./T)
    tbl_cr_ce_HeII%primes(iT) = (473638./T - 0.397 - 0.5d0 * sqrt(T5) / f) &
      * log(10d0) * tbl_cr_ce_HeII%rates(iT)

    ! Recombination Cooling (Hui&Gnedin'97)
    laHII    = 315614./T
    laHeII   = 570670./T
    laHeIII  = 1263030./T
    if(.not. rt_otsa) then ! Case A
      f = 1d0+(laHII/0.541)**0.502
      tbl_cr_r_HII%rates(iT)    = 1.778d-29 * laHII**1.965 / f**2.697 * T
      tbl_cr_r_HII%primes(iT)   = (-0.965 + 1.35389*(f-1.)/f)             &
        * log(10d0) * tbl_cr_r_HII%rates(iT)

      tbl_cr_r_HeII%rates(iT)   = 3d-14 * laHeII**0.654 * kB * T
      tbl_cr_r_HeII%primes(iT)  = 0.346                                   &
        * log(10d0) * tbl_cr_r_HeII%rates(iT)

      f = 1d0+(laHeIII/0.541)**0.502
      tbl_cr_r_HeIII%rates(iT) = 14.224d-29 * laHeIII**1.965 / f**2.697 * T
      tbl_cr_r_HeIII%primes(iT)= (-0.965 + 1.35389*(f-1.)/f)              &
        * log(10d0) * tbl_cr_r_HeIII%rates(iT)
    else ! Case B
      f = 1d0+(laHII/2.25)**0.376
      tbl_cr_r_HII%rates(iT)    = 3.435d-30 * laHII**1.97 / f**3.72 * T
      tbl_cr_r_HII%primes(iT)   = (-0.97 + 1.39827*(f-1.)/f)              &
        * log(10d0) * tbl_cr_r_HII%rates(iT)

      tbl_cr_r_HeII%rates(iT)   = 1.26d-14 * laHeII**0.75 * kB * T
      tbl_cr_r_HeII%primes(iT)  = 0.25                                    &
        * log(10d0) * tbl_cr_r_HeII%rates(iT)

      f = 1d0+(laHeIII/2.25)**0.376
      tbl_cr_r_HeIII%rates(iT)  = 27.48d-30 * laHeIII**1.97 / f**3.72 * T
      tbl_cr_r_HeIII%primes(iT) = (-0.97 + 1.39827*(f-1.)/f)              &
        * log(10d0) * tbl_cr_r_HeIII%rates(iT)
    endif

    ! Bremsstrahlung from Osterbrock & Ferland 2006
    tbl_cr_bre%rates(iT)  = 1.42d-27 * 1.5 * sqrt(T)
    tbl_cr_bre%primes(iT) = 0.5                                            &
      * log(10d0) * tbl_cr_bre%rates(iT)

    ! Compton Cooling from Haimann et al. 96, via Maselli et al.
    ! Need to make sure this is done whenever the redshift changes!
    Ta     = 2.727/aexp
    tbl_cr_com%rates(iT)   = 1.017d-37 * Ta**4 * (T-Ta)
    tbl_cr_com%primes(iT)  = T / (T-Ta)                                    &
      * log(10d0) * tbl_cr_com%rates(iT)

    ! Dielectronic recombination cooling, from Black 1981
    f = 1.24d-13*T**(-1.5d0)*exp(-470000d0/T)
    tbl_cr_die%rates(iT) = f*(1d0+0.3d0*exp(-94000d0/T))
    tbl_cr_die%primes(iT)=0d0
    if(tbl_cr_die%rates(iT) .gt. 0d0) then ! Can simplify w algebra
      tbl_cr_die%primes(iT) = (tbl_cr_die%rates(iT)*(564000.-1.5*T)       &
        - f*94000.) /T**2 * T * log(10d0)
    endif

    !BEGIN H2 STUFF --------------------------------------------------------

    ! Collisional dissociation ground state (Dove&Mandy 1986) [cm3 s-1]
    tbl_Beta_H2HI%rates(iT) =                                              &
      7.073d-19*(T**2.012)*exp(-5.179d4/T)/(1d0+2.130d-5*T)**3.512
    ! Col dissociation ground state (Martin&Keogh&Mandy 1998)  [cm3 s-1]
    tbl_Beta_H2H2%rates(iT) = &
      5.996d-30*(T**4.1881)*exp(-5.466d4/T)/(1d0+6.761d-6*T)**5.6881
    ! Three body H2 formation 3H -> H2 + H (Forrey 2013) [cm6 s-1]
    ! For H2+2H -> 2H2 use Beta_H3B/8 (Palla 1983)
    tbl_Beta_H3B%rates(iT) = &
      6d-32*T**(-0.25)+2d-31*T**(-0.5)
    ! Prefer to use numerical prime, rather than analytic:
    TT=T*1.0001
    tbl_Beta_H2HI%primes(iT) =                                             &
      7.073d-19*(TT**2.012)*exp(-5.179d4/TT)/(1d0+2.130d-5*TT)**3.512
    tbl_Beta_H2HI%primes(iT) = &
      (tbl_Beta_H2HI%primes(iT) - tbl_Beta_H2HI%rates(iT)) &
      / (0.0001*T) * T * log(10d0) ! last two terms for log-log
    tbl_Beta_H2H2%primes(iT) = &
      5.996d-30*(TT**4.1881)*exp(-5.466d4/TT)/(1d0+6.761d-6*TT)**5.6881
    tbl_Beta_H2H2%primes(iT) = &
      (tbl_Beta_H2H2%primes(iT) - tbl_Beta_H2H2%rates(iT))              &
      / (0.0001*T) * T * log(10d0) ! last two terms for log-log
    tbl_Beta_H3B%primes(iT) = &
      6d-32*TT**(-0.25)+2d-31*TT**(-0.5)
    tbl_Beta_H3B%primes(iT) = &
      (tbl_Beta_H3B%primes(iT) - tbl_Beta_H3B%rates(iT))              &
      / (0.0001*T) * T * log(10d0) ! last two terms for log-log

    ! Collisional H2 cooling************************************************
    ! Hallenbach McKee (1979) + Halle Combes (2012) in cgs
    !lowrleft=1.25*exp(-1.70d2/T)*2.35d-14
    !lowrright=1.75*exp(-5.05d2/T)*6.97d-14
    lowrleft=1.25*exp(-1.70d2/T*3d0)*2.35d-14*3d0   !Eqn 6.40 (left term) rotational
    lowrright=1.75*exp(-5.05d2/T*1.683)*6.97d-14*1.683  !Eqn 6.40 (right term) rotational
    lowr_hi=lowrleft*gamma_hi(T,2d0)+lowrright*gamma_hi(T,3d0)
    lowr_h2=lowrleft*gamma_h2(T,2d0)+lowrright*gamma_h2(T,3d0)
    lowvleft_hi=exp(-5860.0/T)*8.09d-13*1.0d-12*sqrt(T)*exp(-1.0d3/T)
    lowvright_hi=exp(-11720.0/T)*1.6d-12*1.6d-12*sqrt(T)*exp(-1.0*(4.0d2/T)**2)
    lowv_hi=lowvleft_hi+lowvright_hi
    lowv_h2=exp(-5860.0/T)*8.09d-13*1.4d-12*sqrt(T)*exp(-1.2d4/(T+1.2d3))

    tbl_cr_H2HI%rates(iT)=lowr_hi+lowv_hi
    tbl_cr_H2H2%rates(iT)=lowr_h2+lowv_h2

    ! Collisional H2 cooling derivatives (calculated numerically)
    TT=T*1.0001
    lowrleft=1.25*exp(-1.70d2/TT)*2.35d-14
    lowrright=1.75*exp(-5.05d2/TT)*6.97d-14
    lowr_hi=lowrleft*gamma_hi(TT,2d0)+lowrright*gamma_hi(TT,3d0)
    lowr_h2=lowrleft*gamma_h2(TT,2d0)+lowrright*gamma_h2(TT,3d0)
    lowvleft_hi=exp(-5860.0/TT)*8.09d-13*1.0d-12*sqrt(TT)*exp(-1.0d3/TT)
    lowvright_hi=exp(-11720.0/TT)*1.6d-12*1.6d-12*sqrt(TT)*exp(-1.0*(4.0d2/TT)**2)
    lowv_hi=lowvleft_hi+lowvright_hi
    lowv_h2=exp(-5860.0/TT)*8.09d-13*1.4d-12*sqrt(TT)*exp(-1.2d4/(TT+1.2d3))

    tbl_cr_H2HI%primes(iT)=lowr_hi+lowv_hi
    tbl_cr_H2HI%primes(iT) = &
      (tbl_cr_H2HI%primes(iT) - tbl_cr_H2HI%rates(iT)) &
      / (0.0001*T) * T * log(10d0) ! last two terms for log-log

    tbl_cr_H2H2%primes(iT)=lowr_h2+lowv_h2
    tbl_cr_H2H2%primes(iT) = &
      (tbl_cr_H2H2%primes(iT) - tbl_cr_H2H2%rates(iT)) &
      / (0.0001*T) * T * log(10d0) ! last two terms for log-log

    !END H2 STUFF ----------------------------------------------------------

  end subroutine comp_table_rates

!PRIVATEXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  SUBROUTINE update_table_rates(iT,aexp)
! Update index iTK in all compton cooling rates tables due to change in
! aexp
!-------------------------------------------------------------------------
    implicit none
    integer::iT
    real(dp)::aexp, T, Ta
!-------------------------------------------------------------------------
    ! Rates are stored in log, while temperature derivatives (primes) are
    ! stored in non-log
    T = 10d0**T_lookup(iT)
    ! Compton Cooling from Haimann et al. 96, via Maselli et al.
    Ta     = 2.727/aexp
    tbl_cr_com%rates(iT)   = 1.017d-37 * Ta**4 * (T-Ta)
    tbl_cr_com%primes(iT)  = T / (T-Ta)                                    &
      * log(10d0) * tbl_cr_com%rates(iT)

  end subroutine update_table_rates

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  FUNCTION inp_coolrates_table(rates_table, T, allow_negative, retPrime, lookup)
! Returns TABULATED rate value from given table
! rates        => Rates (and primes) table to interpolate
! T            => Temperature [K]
! retPrime     <= Temperature derivative of rate at T (optional).
! lookup       <> Caller-owned cache of the last T (optional): the index and
!                 spline weights are reused when T equals lookup%Tlast.
!                 Without it they are recomputed (one log10 per call).
!-------------------------------------------------------------------------
    implicit none
    type(coolrates_table)::rates_table
    !real(dp),dimension(nbinTK)::rates, ratesPrime
    real(dp),intent(in)::T
    logical,intent(in)::allow_negative
    real(dp),optional::retPrime
    type(coolrates_lookup),optional,intent(inout)::lookup
    real(dp)::inp_coolrates_table
    ! No save variables here: the upstream save cache is not thread safe
    ! (a thread could reuse the index of another thread's T), and a
    ! threadprivate cache is slow with ifx (TODO.md 2d-3).
    integer:: iT
    real(dp):: facT, yy, yy2, yy3, fa, fb, fprimea, fprimeb
    real(dp):: alpha, beta, gamma
    logical::extrap, hit
!-------------------------------------------------------------------------
    hit=.false.
    if(present(lookup))then
      if(T .eq. lookup%Tlast)then   ! Reuse index if same T as last call
        hit=.true.
        iT=lookup%iT; facT=lookup%facT; extrap=lookup%extrap
        yy=lookup%yy; yy2=lookup%yy2; yy3=lookup%yy3
      endif
    endif
    if(.not. hit)then
      ! Log of T, snapped to table at the lower boundary, but allowed
      ! to go above upper boundary, in which case we use extrapolation:
      facT = MAX( log10(T), T_lookup(1) )
      extrap=.false.
      if(facT .gt. T_lookup(nbinT)) extrap=.true. ! Above upper limit
      ! Lower closest index in table:
      iT = MIN(MAX(int((facT-T_lookup(1))*dlogTinv)+1, 1), nbinT-1)
      yy=facT-T_lookup(iT)  ! Dist., in log(T), from T to lower table index
      yy2=yy*yy             ! That distance squared
      yy3=yy2*yy            ! ...and cubed
      if(present(lookup))then
        lookup%Tlast=T; lookup%iT=iT; lookup%facT=facT; lookup%extrap=extrap
        lookup%yy=yy; lookup%yy2=yy2; lookup%yy3=yy3
      endif
    endif

    if(extrap) then ! TK above upper table limit, so extrapolate in log-log:
      alpha = (log10(rates_table%rates(nbinT))             &
        -log10(rates_table%rates(nbinT-1)) )        &
        / (T_lookup(nbinT)-T_lookup(nbinT-1))
      inp_coolrates_table = 10d0**(log10(rates_table%rates(nbinT))        &
        + alpha * (facT - T_lookup(nbinT)))
      if( present(retPrime) )                                             &
        retPrime = alpha * inp_coolrates_table / T
      return
    endif

    fa = rates_table%rates(iT)           !
    fb = rates_table%rates(iT+1)         !  Values at neighbouring table
    fprimea=rates_table%primes(iT)       !  indexes
    fprimeb=rates_table%primes(iT+1)     !

    ! Spline interpolation:
    alpha = fprimea
    beta =(fb-fa) * three_over_h2Table - (2d0*fprimea+fprimeb) * one_over_hTable
    gamma = (fprimea+fprimeb) * one_over_h2Table - (fb-fa) * two_over_h3Table
    inp_coolrates_table = fa+alpha*yy+beta*yy2+gamma*yy3
    if(.not. allow_negative) then
      ! Only allow positive rates (spline can go negative):
      inp_coolrates_table = max(0d0,inp_coolrates_table)
    endif
    if( present(retPrime) )                                                &
      retPrime = (alpha+2d0*beta*yy+3d0*gamma*yy2) / T * one_over_lnTen
  end function inp_coolrates_table

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  FUNCTION compCoolrate(T, ne, nH, mu, nN, nI, dcooldT)

! Compute cooling rate in a cell, using interpolation from cooling rate
! tables
! T        => Cell emperature [K]
! ne       => Electron number density [cm-3]
! nH       => Hydrogen number density [cm-3] (H2 options only)
! mu       => Mean molecular weight (H2_fesc only)
! nN       => Neutral and molecular abundances for each ionisation species
! nI       => Ionized and dissociated abundances for each species
! dcooldT  <= Temperature derivative of the rate
! dcooldx  <= Ionized fraction derivative of the rate
! returns:    Cooling rate [erg s-1 cm-3]
!-------------------------------------------------------------------------
    use rt_parameters,only: nIons, isH2, isHe, ixHI, ixHII, ixHeII, ixHeIII
    implicit none
    real(dp)::T, ne, nH, mu
    real(dp),dimension(nIons)::nN,nI
    real(dp)::compCoolrate, dcooldT !---------------------------------------
    real(dp)::beta_esc, xL_H2, xL_HM, dxL_H2, dxL_HM
    logical::H2extra
    real(dp)::ci_HI,ci_HeI,ci_HeII
    real(dp)::ci_HI_prime,ci_HeI_prime,ci_HeII_prime
    real(dp)::cr_H2HI,cr_H2H2,cr_H2HI_prime,cr_H2H2_prime
    real(dp)::ce_HI,ce_HeI,ce_HeII
    real(dp)::ce_HI_prime,ce_HeI_prime,ce_HeII_prime
    real(dp)::r_HII,r_HeII,r_HeIII
    real(dp)::r_HII_prime,r_HeII_prime,r_HeIII_prime
    real(dp)::bre, brefac, bre_prime, com, com_prime, die, die_prime
    type(coolrates_lookup)::lk   ! per-call cache: all lookups below use the same T
!-------------------------------------------------------------------------
    lk%Tlast=-1d0
    ! Coll. Ionization Cooling
    ci_HI   = inp_coolrates_table(tbl_cr_ci_HI, T, .false., ci_HI_prime, lookup=lk)   &
      * ne * nN(ixHII) ! ne * nHI
    ci_HI_prime = ci_HI_prime * ne * nN(ixHII)
    ! Collisional excitation cooling
    ce_HI   = inp_coolrates_table(tbl_cr_ce_HI, T, .false., ce_HI_prime, lookup=lk)   &
      * ne * nN(ixHII) ! ne * nHI
    ce_HI_prime    = ce_HI_prime * ne * nN(ixHII)
    ! Recombination Cooling
    r_HII = inp_coolrates_table(tbl_cr_r_HII, T, .false., r_HII_prime, lookup=lk)     &
      * ne * nI(ixHII) ! ne * nHII
    r_HII_prime   = r_HII_prime * ne * nI(ixHII)

    if(isHe) then
      ! Coll. Ionization Cooling
      ci_HeI      = inp_coolrates_table(tbl_cr_ci_HeI, T, .false.         &
        , ci_HeI_prime, lookup=lk)   &
        * ne * nN(ixHeII)  ! ne * nHeI
      ci_HeI_prime= ci_HeI_prime  * ne * nN(ixHeII)
      ci_HeII     = inp_coolrates_table(tbl_cr_ci_HeII, T, .false.        &
        , ci_HeII_prime, lookup=lk) &
        * ne * nN(ixHeIII) ! ne * nHeII
      ci_HeII_prime= ci_HeII_prime * ne * nN(ixHeIII)
      ! Collisional excitation cooling
      ce_HeI      = inp_coolrates_table(tbl_cr_ce_HeI, T, .false.         &
        , ce_HeI_prime, lookup=lk)   &
        * ne * nN(ixHeII) ! ne * nHeI
      ce_HeI_prime= ce_HeI_prime *ne*nN(ixHeII)
      ce_HeII     = inp_coolrates_table(tbl_cr_ce_HeII, T, .false.        &
        , ce_HeII_prime, lookup=lk) &
        * ne * nN(ixHeIII) ! ne * nHeII
      ce_HeII_prime= ce_HeII_prime * ne * nN(ixHeIII)
      ! Recombination Cooling
      r_HeII      = inp_coolrates_table(tbl_cr_r_HeII, T, .false.         &
        , r_HeII_prime, lookup=lk)   &
        * ne * nI(ixHeII) ! ne * nHeII
      r_HeII_prime= r_HeII_prime * ne * nI(ixHeII)
      r_HeIII     = inp_coolrates_table(tbl_cr_r_HeIII, T, .false.        &
        , r_HeIII_prime, lookup=lk) &
        * ne * nI(ixHeIII) ! ne * nHeIII
      r_HeIII_prime= r_HeIII_prime * ne * nI(ixHeIII)
    endif

    ! Bremsstrahlung
    bre    = inp_coolrates_table(tbl_cr_bre, T, .false., bre_prime, lookup=lk)
    brefac = ne * nI(ixHII)
    if(isHe) then
      brefac = brefac + ne * ( nI(ixHeII) + 4. * nI(ixHeIII) )
    endif
    bre        = bre * brefac
    bre_prime  = bre_prime * brefac

    ! Compton Cooling
    com       = inp_coolrates_table(tbl_cr_com, T, .true., com_prime, lookup=lk) * ne
    com_prime = com_prime * ne

    ! Dielectronic recombination cooling
    if(isHe) then
      die       = inp_coolrates_table(tbl_cr_die, T, .false., die_prime, lookup=lk) &
        * ne * nN(ixHeIII) ! ne * nHeII
      die_prime = die_prime * ne * nN(ixHeIII)
    endif

    H2extra = .false.
    if(isH2) then
      select case (trim(H2coolfunc))
      case ('Galli')
        ! Galli & Palla 1998 (H collisions, HM79 LTE limit); cr_H2H2 unused
        cr_H2HI = 0d0 ; cr_H2HI_prime = 0d0 ; cr_H2H2 = 0d0 ; cr_H2H2_prime = 0d0
        if(T .lt. T_H2chem_max) then
          cr_H2HI       = GalliPalla_rate(nN(ixHII), nN(ixHI), T)
          cr_H2HI_prime = (GalliPalla_rate(nN(ixHII), nN(ixHI), T*1.0001d0)  &
            - cr_H2HI) / (0.0001d0*T)
        endif
      case ('Glover')
        ! Glover & Abel 2008 (H, H2, He, H+, e- collisions) + Glover 2015 LTE
        cr_H2HI = 0d0 ; cr_H2HI_prime = 0d0 ; cr_H2H2 = 0d0 ; cr_H2H2_prime = 0d0
        if(T .lt. T_H2chem_max) then
          cr_H2HI       = Glover_rate(nN(ixHII), nN(ixHI), nN(ixHeII)      &
            , nI(ixHII), ne, T, .false.)
          cr_H2HI_prime = (Glover_rate(nN(ixHII), nN(ixHI), nN(ixHeII)     &
            , nI(ixHII), ne, T, .true.) - cr_H2HI) / (0.0001d0*T)
        endif
      case default
        ! '' or 'Hollenbach': Hollenbach & McKee 1979 + Halle & Combes 2012
        ! low-density rates from the tables (upstream RAMSES-RT)
        cr_H2HI       = inp_coolrates_table(tbl_cr_H2HI, T, .false.       &
          , cr_H2HI_prime, lookup=lk)  &
          * nN(ixHI) * nN(ixHII) ! nH2 * nHI
        cr_H2HI_prime = cr_H2HI_prime * nN(ixHI)*nN(ixHII)

        cr_H2H2       = inp_coolrates_table(tbl_cr_H2H2, T, .false.       &
          , cr_H2H2_prime, lookup=lk)  &
          * nN(ixHI)**2 ! nH2 * nH2
        cr_H2H2_prime = cr_H2H2_prime * nN(ixHI)**2
      end select
      if(H2_fesc) then
        ! H2 line escape fraction (Fukushima, Omukai & Hosokawa 2018)
        beta_esc      = fesc_H2(T, nH, nN(ixHI), mu)
        cr_H2HI       = cr_H2HI * beta_esc
        cr_H2HI_prime = cr_H2HI_prime * beta_esc
        cr_H2H2       = cr_H2H2 * beta_esc
        cr_H2H2_prime = cr_H2H2_prime * beta_esc
      endif
      if((H2_Hm_cool .or. H2_chem_cool) .and. T .lt. T_H2chem_max) then
        ! H- radiative attachment cooling and H2 chemical cooling/heating
        H2extra = .true.
        call chemical_and_Hm_cool(T, nH, nN, nI, xL_H2, xL_HM)
        call chemical_and_Hm_cool(T*1.0001d0, nH, nN, nI, dxL_H2, dxL_HM)
        dxL_H2 = (dxL_H2 - xL_H2) / (0.0001d0*T)
        dxL_HM = (dxL_HM - xL_HM) / (0.0001d0*T)
      endif
    endif

    ! Overall Cooling
    compCoolrate          = ci_HI    + r_HII    + ce_HI    + com + bre
    if(isHe) compCoolrate = compCoolrate                                   &
      + ci_HeI   + r_HeII   + ce_HeI                   &
      + ci_HeII  + r_HeIII  + ce_HeII  + die
    if(isH2) compCoolrate = compCoolrate                                   &
      + cr_H2HI  + cr_H2H2
    if(H2extra) compCoolrate = compCoolrate + xL_H2 + xL_HM

    dCooldT= ci_HI_prime + r_HII_prime + ce_HI_prime + com_prime + bre_prime
    if(isHe) dCooldT = dCooldT                                             &
      + ci_HeI_prime  + r_HeII_prime  + ce_HeI_prime        &
      + ci_HeII_prime + r_HeIII_prime + ce_HeII_prime       &
      + die_prime
    if(isH2) dCooldT = dCooldT                                             &
      + cr_H2HI_prime + cr_H2H2_prime
    if(H2extra) dCooldT = dCooldT + dxL_H2 + dxL_HM

  end function compCoolrate

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ELEMENTAL FUNCTION gamma_hi(T,J)
    !-------------------------------------------------------------------------
    implicit none
    real(dp),intent(in)::T,J
    real(dp)::gamma_hi, T3, jconsts
    !-------------------------------------------------------------------------
    T3 = T/1d3
    jconsts=0.33+0.9*exp(-1.0d0*((J-3.5d0)/0.9d0)**2)
    gamma_hi = &
      jconsts*(1.0d-11*sqrt(T3)/(1.0d0+60.0d0*T3**(-4)) + 1.0d-12*T3)

  end function gamma_hi
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ELEMENTAL FUNCTION gamma_h2(T,J)
    !-------------------------------------------------------------------------
    implicit none
    real(dp),intent(in)::T,J
    real(dp)::gamma_h2, T3, jconsts
    !----------------------------------------------------------
    T3 = T/1d3
    jconsts=(0.276*J**2)*exp(-1.0d0*(J/3.18d0)**1.7)
    gamma_h2 = jconsts*(3.3d-12 +6.6d-12*T3)

  end function gamma_h2

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! H2 cooling options ported from patch/ks (Ramses-FG) 2026-11-12.
! Physical constants come from the constants module (CODATA 2018) instead of
! the hand-written values of the original.
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ELEMENTAL FUNCTION Lrot_LTE(T)
! H2 rotational cooling in LTE, Hollenbach & McKee 1979 [erg s-1 per H2]
    implicit none
    real(dp),intent(in)::T
    real(dp)::Lrot_LTE, T3
    T3=T/1d3
    Lrot_LTE=( (9.5d-22*T3**3.76)/(1d0+0.12*T3**2.1) * exp(-(0.13/T3)**3) + &
      & 3d-24 * exp(-0.51/T3) )
  end function Lrot_LTE

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ELEMENTAL FUNCTION Lvib_LTE(T)
! H2 vibrational cooling in LTE, Hollenbach & McKee 1979 [erg s-1 per H2]
    implicit none
    real(dp),intent(in)::T
    real(dp)::Lvib_LTE, T3
    T3=T/1d3
    Lvib_LTE=( 6.7d-19 * exp(-5.86/T3) + 1.6d-18 * exp(-11.7/T3) )
  end function Lvib_LTE

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ELEMENTAL FUNCTION GalliPalla_rate(nHI, nH2, T)
! H2 cooling rate [erg s-1 cm-3] of Galli & Palla 1998 (Appendix): the
! low-density limit for H collisions (eq. A7, 13 K < T < 1e5 K) joined to
! the HM79 LTE limit as Lambda = L_LTE / (1 + L_LTE/L_0).
    implicit none
    real(dp),intent(in)::nHI, nH2, T
    real(dp)::GalliPalla_rate, logT, Lamb0, LambLTE
    logT=log10(T)
    Lamb0=-103.0+97.59*logT-48.05*logT**2+10.80*logT**3-0.9032*logT**4 ! A7
    Lamb0=10d0**Lamb0                                   ! [erg cm3 s-1]
    LambLTE=Lrot_LTE(T)+Lvib_LTE(T)                     ! [erg s-1]
    ! = nHI*nH2*(LambLTE/nHI)/(1+(LambLTE/nHI)/Lamb0), written without 1/nHI
    GalliPalla_rate = nH2*LambLTE*nHI*Lamb0/(nHI*Lamb0 + LambLTE)
  end function GalliPalla_rate

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  FUNCTION Glover_rate(nHI,nH2,nHe,nHII,ne,T,prime)
! H2 cooling rate [erg s-1 cm-3]: low-density limit of Glover & Abel 2008
! (Table 8; H, H2, He, H+, e- collisions) joined to the LTE limit of
! Glover 2015 (eq. 30) as Lambda = L_LTE / (1 + L_LTE/L_0).
! prime=.true. evaluates at T*1.0001 (for the numerical T derivative).
    implicit none
    real(dp)::Glover_rate
    integer::i
    real(dp),intent(in)::nHI,nH2,nHe,nHII,ne,T
    logical,intent(in)::prime
    real(dp)::logT3,Lamb_nk,Lamb0,LambLTE
    real(dp),dimension(0:5)::a ! table 8 (2008) low density limit
    real(dp),dimension(0:8)::b ! equation 30 (2015) LTE limit
    real(dp)::logT3_LTE,ramp_LTE100
    real(dp),parameter:: dex_ramp_LTE100=0.2d0 ! width (dex) of the ramp to zero below 100 K

    logT3=log10(T/1d3)
    if(prime)logT3=log10(T*1.0001/1d3)
    Lamb0=0d0
    LambLTE = 0d0
    ! low density limit - Glover & Abel 2008 - Table 8
    ! H2 with H
    if(T.gt.10 .and. T.le.6000)then
      Lamb_nk=0d0
      if(T.gt.10 .and. T.le.100)then
        a(0) = -16.818342
        a(1) =  37.383713
        a(2) =  58.145166
        a(3) =  48.656103
        a(4) =  20.159831
        a(5) =  3.8479610
      else if(T.gt.100 .and. T.le.1000)then
        a(0) = -24.311209
        a(1) =  3.5692468
        a(2) = -11.332860
        a(3) = -27.850082
        a(4) = -21.328264
        a(5) = -4.2519023
      else
        a(0) = -24.311209
        a(1) =  4.6450521
        a(2) = -3.7209846
        a(3) =  5.9369081
        a(4) = -5.5108047
        a(5) =  1.5538288
      endif
      do i=0,5
        Lamb_nk = Lamb_nk + a(i) * logT3**i
      enddo
      Lamb0 = Lamb0 + 10d0**Lamb_nk * nHI
    endif
    ! H2 with H2
    if(T.gt.100 .and. T.le.6000)then
      Lamb_nk=0d0
      a(0) = -23.962112
      a(1) =  2.09433740
      a(2) = -0.77151436
      a(3) =  0.43693353
      a(4) = -0.14913216
      a(5) = -0.033638326
      do i=0,5
        Lamb_nk = Lamb_nk + a(i) * logT3**i
      enddo
      Lamb0 = Lamb0 + 10d0**Lamb_nk * nH2
    endif
    ! H2 with He
    if(T.gt.10 .and. T.le.6000)then
      Lamb_nk = 0d0
      a(0) = -23.689237
      a(1) =  2.1892372
      a(2) = -0.81520438
      a(3) =  0.29036281
      a(4) = -0.16596184
      a(5) =  0.19191375
      do i=0,5
        Lamb_nk = Lamb_nk + a(i) * logT3**i
      enddo
      Lamb0 = Lamb0 + 10d0**Lamb_nk * nHe
    endif
    ! H2 with H+
    if(T.gt.10 .and. T.le.10000)then
      Lamb_nk = 0d0
      a(0) = -21.716699
      a(1) =  1.3865783
      a(2) = -0.37915285
      a(3) =  0.11453688
      a(4) = -0.23214154
      a(5) =  0.058538864
      do i=0,5
        Lamb_nk = Lamb_nk + a(i) * logT3**i
      enddo
      Lamb0 = Lamb0 + 10d0**Lamb_nk * nHII
    endif
    ! H2 with e- (10 < T <= 10000 K, two coefficient sets)
    if(T.gt.10 .and. T.le.10000)then
      Lamb_nk = 0d0
      if(T.le.200)then
        a(0) = -34.286155
        a(1) = -48.537163
        a(2) = -77.121176
        a(3) = -51.352459
        a(4) = -15.169160
        a(5) = -0.98120322
      else
        a(0) = -22.190316
        a(1) =  1.5728955
        a(2) = -0.21335100
        a(3) =  0.96149759
        a(4) = -0.91023195
        a(5) =  0.13749749
      endif
      do i=0,5
        Lamb_nk = Lamb_nk + a(i) * logT3**i
      enddo
      Lamb0 = Lamb0 + 10d0**Lamb_nk * ne
    endif
    ! high density (LTE) limit, Glover 2015 eq. 30, valid 100 K < T < 1e4 K.
    ! Below 100 K the fit is evaluated at 100 K and tapered to zero over
    ! dex_ramp_LTE100 in log10(T) (avoids a discontinuity at 100 K).
    if(T.lt.10000)then
      ramp_LTE100 = 1d0
      if(T.lt.100)then
        if(T.lt.100d0*10d0**(-dex_ramp_LTE100)) then
          ramp_LTE100 = 0d0
        else
          ramp_LTE100 = 1d0 - log10(100d0/T)/dex_ramp_LTE100
        endif
      endif
      if(ramp_LTE100.gt.0d0)then
        logT3_LTE=log10(max(T,100d0)/1d3)
        if(prime)logT3_LTE=log10(max(T,100d0)*1.0001/1d3)
        b(0) = -20.584225
        b(1) =  5.0194035
        b(2) = -1.5738805
        b(3) = -4.7155769
        b(4) =  2.4714161
        b(5) =  5.4710750
        b(6) = -3.9467356
        b(7) = -2.2148338
        b(8) =  1.8161874
        do i=0,8
          LambLTE = LambLTE + b(i)*logT3_LTE**i
        enddo
        LambLTE = 10d0**LambLTE * ramp_LTE100          ! [erg s-1]
      endif
    endif
    if(Lamb0.eq.0d0)then
      Glover_rate=0d0
    else
      Glover_rate=LambLTE/(1d0+LambLTE/Lamb0) * nH2
    endif
  end function Glover_rate

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  FUNCTION fesc_H2(T,nH,nH2,mu)
! Escape fraction of H2 line emission, Fukushima, Omukai & Hosokawa 2018
! (eqs. B3, B4, Tables B1, B2), with the H2 column density taken as
! nH2 x (half the Jeans length). Valid for nH < 1e12 cm-3.
    use cooling_module,only:X
    implicit none
    real(dp)::fesc_H2
    real(dp),intent(in)::T,nH,nH2,mu ! nH: total hydrogen number density
    integer::i
    real(dp),dimension(0:5)::a, b
    real(dp)::alpha,Nc,Nc_inv,cs,tff,RJeans,rr,NH2col
    alpha=0d0
    Nc=0d0
    if(T.le.1000)then
      a(0)= 9.78382d-1
      a(1)= 3.99572d-4
      a(2)=-1.73108d-6
      a(3)= 1.15363d-9
      a(4)= 8.24607d-13
      a(5)=-7.65975d-16
    else if(T.gt.1000 .and. T.le.4000)then
      a(0)= 8.27472d-1
      a(1)= 1.07697d-4
      a(2)=-8.25123d-8
      a(3)= 1.92812d-12
      a(4)= 5.7192d-15
      a(5)=-7.869d-19
    else
      a(0)= 1.11569d0
      a(1)=-3.29302d-4
      a(2)= 1.01846d-7
      a(3)=-1.46666d-11
      a(4)= 1.00764d-15
      a(5)=-2.68873d-20
    endif
    b(0)= 24.0561
    b(1)= 1.10043d-3
    b(2)=-2.87224d-7
    b(3)= 6.11525d-11
    b(4)=-6.55034d-15
    b(5)= 2.54997d-19
    do i=0,5
      alpha=alpha+a(i)*T**i
      Nc=Nc+b(i)*T**i
    enddo
    Nc_inv=10d0**(-Nc)                       ! 1/Nc (avoids overflow at large T)
    cs=sqrt( kB*T / ( mu*mH ) )              ! isothermal sound speed
    rr=nH/X*mH                               ! total mass density [g cm-3]
    tff=sqrt( (twopi/2d0) / (factG_in_cgs*rr) ) ! sqrt(pi/(G rho)), as in the original
    RJeans=0.5d0*cs*tff                      ! half the Jeans length
    NH2col=nH2*RJeans                        ! H2 column density [cm-2]
    fesc_H2=(1d0+NH2col*Nc_inv)**(-alpha)
  end function fesc_H2

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ELEMENTAL FUNCTION comp_Beta_H2coll(T,nH,xHI,xH2,xHe,ne,nHI,nH2,nHeI)
! Collisional dissociation rate of H2 [s-1 per H2 ... /nH, i.e. cm3 s-1]
! for H2 + (e-, H, H2, He, H+) from Glover & Abel 2008 (eqs. 13-17 and
! appendix), with the LTE/non-LTE blend set by the critical densities.
! Used by chemical_and_Hm_cool (H2_chem_cool) only.
    use cooling_module,only:X, Y
    implicit none
    real(dp),intent(in)::T,nH,xHI,xH2,xHe,ne,nHI,nH2,nHeI
    real(dp)::comp_Beta_H2coll,ncrH,ncrH2,ncrHe,T4,invncr
    real(dp):: LTEfac, NLTEfac
    real(dp):: k7, k8, k9, k9L, k10, k10L, k11, k11L
    real(dp):: lk7,lk8,lk9,lk10,lk11,lnT,nHII
    T4 = T/1d4
    ! Critical number densities (eqs. 15-17)
    ncrH = 10.0**(3.0 - 0.416*log10(T4) - 0.327*log10(T4)*log10(T4))
    ncrH2 = 10.0**(4.845 - 1.3*log10(T4) + 1.62*log10(T4)*log10(T4))
    ncrHe = 10.0**(5.0792*(1.0 - 1.23d-5*(T-2000.0)))
    ncrHe = max(ncrHe,1d-40)        ! avoids an FPE at T > 1e7 K
    ! 1/ncr (eq. 14); xHe here is per H atom, converted to the He fraction
    invncr = (xHI/ncrH) + (xH2/ncrH2) + (xHe/ncrHe)*(Y/4.*X)
    LTEfac = (nH*invncr)/(1.0 + (nH*invncr))
    NLTEfac = 1.0/(1.0 + (nH*invncr))
    if(ne.gt.nHI)then
      LTEfac=1.d0
      NLTEfac=0.d0
    endif
    nHII = nH-nHI-nH2*2d0
    ! H2 + H+ -> H2+ + H (7th-order polynomial in ln T, evaluated at
    ! min(T,T_H2chem_fitmax); the activation factor uses the true T)
    lnT=log(min(T, T_H2chem_fitmax))
    k7   = -3.3232183d-7         +3.3735382d-7 *lnT    -1.4491368d-7 *lnT**2 &
      & +3.4172805d-8 *lnT**3 -4.7813720d-9 *lnT**4 +3.9731542d-10*lnT**5 &
      & -1.8171411d-11*lnT**6 +3.5311932d-13*lnT**7
    k7   = k7*exp(-21237.15/T)
    k7   = max(k7,1d-40)
    lk7  = log10(k7)
    ! H2 + e- -> H + H + e- (Glover et al. 2010)
    k8 = (3.73d-9)*(T**0.1121)*exp(-99430.0/T)
    ! H2 + H -> H + H + H
    k9 = (6.67d-12)*sqrt(T)*exp(-1.0*(1.0 + (63593.0/T)))
    k9L = (3.52d-9)*exp(-43900.0/T)
    ! H2 + H2 -> H2 + H + H
    k10 = (((5.996d-30)*(T**4.1881))/((1.0 + (6.761d-6)*T)**5.6881))*exp(-54657.4/T)
    k10L = (1.3d-9)*exp(-53300.0/T)
    ! H2 + He -> H + H + He
    k11 = 10.0**(-27.029 + (3.801*log10(T)) - (29487.0/T))
    k11L = 10.0**(-2.729 - (1.75*log10(T)) - (23474.0/T))
    k8   = max(k8,  1d-40)
    k9   = max(k9,  1d-40)
    k10  = max(k10, 1d-40)
    k11  = max(k11, 1d-40)
    k9L  = max(k9L, 1d-40)
    k10L = max(k10L,1d-40)
    k11L = max(k11L,1d-40)
    lk8  = log10(k8)
    lk9  = (LTEfac*log10(k9L)) + (NLTEfac*log10(k9))
    lk10 = (LTEfac*log10(k10L)) + (NLTEfac*log10(k10))
    lk11 = (LTEfac*log10(k11L)) + (NLTEfac*log10(k11))
    comp_Beta_H2coll = (ne*(10.0**lk8)) + (nHI*(10.0**lk9)) + (nH2*(10.0**lk10)) + (nHeI*(10.0**lk11))
    comp_Beta_H2coll = comp_Beta_H2coll + nHII*(10.0**lk7)
    comp_Beta_H2coll = MAX(comp_Beta_H2coll/nH, 1.0d-40)
  end function comp_Beta_H2coll

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  subroutine chemical_and_Hm_cool(T, nH, nN, nI, xL_H2, xL_HM)
! H- radiative attachment cooling (H + e -> H- + hnu) and, with
! H2_chem_cool, the chemical cooling/heating of H2 formation (H- channel,
! three-body) and collisional dissociation (cf. Hollenbach & McKee 1979;
! chemistry module of SFUMATO). Rates: Glover et al. 2010, Forrey 2013,
! Palla, Salpeter & Stahler 1983.
!  xL_H2  <= chemical cooling (>0) / heating (<0) [erg s-1 cm-3]
!  xL_HM  <= H- cooling rate [erg s-1 cm-3]
    use cooling_module,only:X, Y
    use rt_parameters,only: nIons, ixHI, ixHII, ixHeII, ixHeIII, isHe
    implicit none
    real(dp),intent(IN)::T, nH
    real(dp),dimension(nIons),intent(IN)::nN, nI
    real(dp),intent(OUT)::xL_H2, xL_HM
    real(dp):: yHI,yHII,yH2,yHeI,yHeII,yHeIII,yel,&
      logT,lnTe,k1,k2,k5,k13,keff_Hm,kra_Hm,k3b1,k3b2,keff_cdis,&
      rate_raHm,rate_form_Hm,rate_form_3b,rate_cdis,xn_cr,crit
    real(dp)::logT2,logT3,logT4,logT5,logT6,lnTe2,lnTe3,lnTe4,lnTe5,lnTe6,lnTe7,lnTe8
    real(dp)::Tfit ! T clamped to the validity range of the H- fits
    ! abundances per H atom
    yHI = nN(ixHII)/nH
    yHII = nI(ixHII)/nH
    yH2 = nN(ixHI)/nH
    yHeI = 0d0 ; yHeII = 0d0 ; yHeIII = 0d0
    if(isHe) then
      yHeI   = nN(ixHeII)/nH
      yHeII  = nI(ixHeII)/nH
      yHeIII = nI(ixHeIII)/nH
    endif
    yel = yHII + yHeII + 2.d0 * yHeIII
    ! The H- channel fits (k1,k2,k5,k13) blow up when extrapolated to high T,
    ! so they are evaluated at Tfit
    Tfit = min(T, T_H2chem_fitmax)
    logT = log10(Tfit)
    lnTe = log(Tfit*kB/eV2erg)   ! K -> eV
    logT2 = logT*logT
    logT3 = logT2*logT
    logT4 = logT3*logT
    logT5 = logT4*logT
    logT6 = logT5*logT
    ! H + e- -> H- + gamma (Glover et al. 2010)
    if(Tfit<6000)then
      k1 = 10.0**(-17.845 + 0.762*logT + 0.1523*logT2 - 0.03274*logT3)
    else
      k1 = 10.0**(-16.420 + 0.1998*logT2 -5.447d-3*logT4 + 4.0415d-5*logT6)
    endif
    ! H- + H -> H2 + e
    if(Tfit<300)then
      k2 = 1.5d-9
    else
      k2 = 4.0d-9*Tfit**(-0.17)
    endif
    if (H2_chem_cool) then
      ! --- H2 chemical cooling/heating and H- cooling ---
      lnTe2 = lnTe*lnTe
      lnTe3 = lnTe2*lnTe
      lnTe4 = lnTe3*lnTe
      lnTe5 = lnTe4*lnTe
      lnTe6 = lnTe5*lnTe
      lnTe7 = lnTe6*lnTe
      lnTe8 = lnTe7*lnTe
      ! H- + H+ -> H + H
      k5 = 2.4d-6/sqrt(Tfit)*(1.0 + Tfit/20000.)
      ! H- + e -> H + e + e
      k13 = -1.801849334d1 + 2.36085220d0*lnTe -2.82744300d-1*lnTe2 &
        & + 1.62331664d-2*lnTe3 - 3.36501203d-2*lnTe4 + 1.17832978d-2*lnTe5 &
        & - 1.65619470d-3*lnTe6 + 1.06827520d-4*lnTe7 - 2.63128581d-6*lnTe8
      k13 = exp(k13)
      ! effective H2 formation rate via H-: ydot(H2) = keff_Hm*nH*y(HI)*y(el)
      keff_Hm =  k1*k2*yHI/(k2*yHI + k5*yHII + k13*yel)   ! [cm3 s-1]
      kra_Hm = k1
      ! three-body: 3H -> H2 + H (Forrey 2013), 2H + H2 -> 2H2 (Palla+1983)
      k3b1 = 6d-32*T**(-0.25) + 2d-31*T**(-0.5)            ! [cm6 s-1]
      k3b2 = k3b1/8d0
      ! effective collisional dissociation: ydot(H2) = -keff_cdis*nH*y(H2)
      keff_cdis = comp_Beta_H2coll(T,nH,yHI,yH2,yHeI*X/(Y/4.),nH*yel,nN(ixHII),nN(ixHI),nN(ixHeII))
      ! rates [cm-3 s-1]
      rate_raHm    = kra_Hm * nH**2 * yHI * yel
      rate_form_Hm = keff_Hm * nH**2 * yHI * yel
      rate_form_3b = k3b1 * nH**3 * yHI**3 +  k3b2 * nH**3 * yHI**2 * yH2
      rate_cdis    = keff_cdis * nH**2 * yH2
      ! critical density for the rovibrational cascade after H2 formation
      xn_cr =  1.d6/sqrt(T) &
        /(1.6d0*yHI*exp(-(4.d2/T)**2) + 1.4d0*yH2*exp(-1.2d4/(T+1.2d3)))
      crit = 1.d0/(1.d0 + xn_cr/nH)
      ! chemical energy deposited by H2 formation/destruction [erg cm-3 s-1]
      xL_H2 = (-4.48d0*crit*rate_form_3b - 3.73d0*crit*rate_form_Hm  + 4.48d0*rate_cdis) * eV2erg
      ! H- cooling: hnu = 0.76 eV + kinetic energy; the 0.76 eV is the
      ! chemical energy difference already counted in xL_H2
      xL_HM  = kB*T * rate_raHm
    else
      ! --- H- cooling only: all of hnu = 0.76 eV + kT leaves the gas ---
      kra_Hm = k1
      rate_raHm = kra_Hm * nH**2 * yHI * yel
      xL_HM  = (0.76d0*eV2erg + kB*T) * rate_raHm
      xL_H2 = 0d0
    end if
  end subroutine chemical_and_Hm_cool

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! LW self-shielding factors (ported from patch/ks 2026-11-12). The column
! densities are nH2*dx_SS, nHI*dx_SS with dx_SS the shielding length [cm].
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ELEMENTAL FUNCTION comp_Sd(nHI,nH2,dx_SS,Z)
! Dust shielding of the LW band, Gnedin, Tassis & Kravtsov 2009 (sec. 2.2)
    implicit none
    real(dp),intent(in)::nHI,nH2,dx_SS,Z
    real(dp)::comp_Sd,Sdeff,cNHI,cNH2
    Sdeff = 4.0d-21    !dust cross section cm^2
    cNHI = nHI*dx_SS   !HI column density
    cNH2 = nH2*dx_SS   !H2 column density
    comp_Sd = exp(-Sdeff*Z*(cNHI + (2.0*cNH2)))
  end function comp_Sd

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ELEMENTAL FUNCTION comp_SH2_G09(nH2,dx_SS)
! H2 self-shielding in the form used by Gnedin+2009 (DB96 with w=0.2)
    implicit none
    real(dp),intent(in)::nH2,dx_SS
    real(dp)::comp_SH2_G09,xfac,cNH2,wH2,Sa,Sb,Sc
    cNH2 = nH2*dx_SS  !H2 column density
    xfac = cNH2/(5.0d14)
    wH2 = 0.2
    Sa = (1.0-wH2)/((1.0+xfac)*(1.0+xfac))
    Sb = wH2/sqrt(1.0+xfac)
    Sc = exp(-0.00085*sqrt(1.0+xfac))
    comp_SH2_G09 = Sa + (Sb*Sc)
  end function comp_SH2_G09

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ELEMENTAL FUNCTION comp_SH2_DB96(nH2,dx_SS)
! H2 self-shielding, Draine & Bertoldi 1996 eq. 37 (w=0.035, b5=1)
    implicit none
    real(dp),intent(in)::nH2,dx_SS
    real(dp)::comp_SH2_DB96,xfac,cNH2,wH2,Sa,Sb,Sc
    cNH2 = nH2*dx_SS  !H2 column density
    xfac = cNH2/(5.0d14)
    wH2 = 0.035
    Sa = (1.0-wH2)/((1.0+xfac)*(1.0+xfac))
    Sb = wH2/sqrt(1.0+xfac)
    Sc = exp(-0.00085*sqrt(1.0+xfac))
    comp_SH2_DB96 = Sa + (Sb*Sc)
  end function comp_SH2_DB96

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ELEMENTAL FUNCTION comp_SH2_WGH(nH2,dx_SS,n,T)
! H2 self-shielding, Wolcott-Green & Haiman 2019 (T- and n-dependent alpha),
! Doppler parameter b = sqrt(b_turb^2 + b_therm^2) with b_turb = 7.1 km/s
! (Krumholz 2012; Richings+2014). n: total number density, T: temperature.
    implicit none
    real(dp),intent(in)::nH2,dx_SS,n,T
    real(dp)::comp_SH2_WGH,xfac,cNH2,wH2,Sa,Sb,Sc
    real(dp)::alpha,A1,A2,c1,c2,c3,c4,c5,b5,btherm,bturb
    cNH2 = nH2*dx_SS ! H2 column density
    xfac = cNH2/5d14
    wH2 = 0.035
    c1 = 0.2856 ; c2 = 0.8711 ; c3 = 1.928 ; c4 = 0.9639 ; c5 = 3.892
    A1 = c2*log10(T)-c3
    A2 = -c4*log10(T)+c5
    alpha = A1*exp(-c1*log10(n))+A2
    bturb=7.1                              ! km/s
    btherm=sqrt( 2.0*kB*T/2.0/mH )/1d5     ! km/s (H2)
    b5=sqrt( bturb**2 + btherm**2 )
    Sa = (1.0-wH2)/(1.0+xfac/b5)**alpha
    Sb = wH2/sqrt(1.0+xfac)
    Sc = exp(-0.00085*sqrt(1.0+xfac))
    comp_SH2_WGH = Sa + (Sb*Sc)
  end function comp_SH2_WGH

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ELEMENTAL FUNCTION comp_Alpha_H2_Katz(T,Z,xe,xHII,nHI,nH2,xHI,fclump)
! H2 formation rate [cm3 s-1] of the ks/Katz+2017 model (comp_Alpha_H2 with
! isH2Katz in patch/ks): dust surface (3.5e-17 Z, Katz+2017) on neutral H,
! the H- channel with the Glover+2010 rates for H+e -> H-, H-+H -> H2+e,
! H-+H+ -> 2H, H-+e -> H+2e (equilibrium H-), and the three-body channels
! 3H -> H2+H (Forrey 2013) and 2H+H2 -> 2H2 (Palla+1983). The H- fits are
! evaluated at min(T, T_H2chem_fitmax). The creation rate is
! alpha * nH * xHI [s-1 per H].
    implicit none
    real(dp),intent(in)::T,Z,xe,xHII,nHI,nH2,xHI,fclump
    real(dp)::comp_Alpha_H2_Katz
    real(dp)::logT, lnTe, k1, k2, k5, k13, Tfit, R31, R32
    real(dp)::logT2,logT3,logT4,logT6,lnTe2,lnTe3,lnTe4,lnTe5,lnTe6,lnTe7,lnTe8
    Tfit = min(T, T_H2chem_fitmax)
    ! dust surface, no dust in ionised gas
    comp_Alpha_H2_Katz = Z * 3.5d-17 * fclump * (1d0-xHII)
    ! H- channel (Glover+2008 sec 2.1.1; competition of associative
    ! detachment with mutual neutralisation and collisional detachment)
    logT = log10(Tfit)
    lnTe = log(Tfit*kB/eV2erg)                     ! K -> eV
    logT2 = logT*logT ; logT3 = logT2*logT ; logT4 = logT3*logT ; logT6 = logT4*logT2
    lnTe2 = lnTe*lnTe ; lnTe3 = lnTe2*lnTe ; lnTe4 = lnTe3*lnTe ; lnTe5 = lnTe4*lnTe
    lnTe6 = lnTe5*lnTe ; lnTe7 = lnTe6*lnTe ; lnTe8 = lnTe7*lnTe
    if(Tfit<6000)then                              ! H + e- -> H- + gamma
      k1 = 10.0**(-17.845 + 0.762*logT + 0.1523*logT2 - 0.03274*logT3)
    else
      k1 = 10.0**(-16.420 + 0.1998*logT2 -5.447d-3*logT4 + 4.0415d-5*logT6)
    endif
    if(Tfit<300)then                               ! H- + H -> H2 + e
      k2 = 1.5d-9
    else
      k2 = 4.0d-9*Tfit**(-0.17)
    endif
    k5 = 2.4d-6/sqrt(Tfit)*(1.0 + Tfit/20000.)     ! H- + H+ -> H + H
    k13 = -1.801849334d1 + 2.36085220d0*lnTe -2.82744300d-1*lnTe2 &
      & + 1.62331664d-2*lnTe3 - 3.36501203d-2*lnTe4 + 1.17832978d-2*lnTe5 &
      & - 1.65619470d-3*lnTe6 + 1.06827520d-4*lnTe7 - 2.63128581d-6*lnTe8
    k13 = exp(k13)                                 ! H- + e -> H + e + e
    comp_Alpha_H2_Katz = comp_Alpha_H2_Katz + k1*k2*xe*xHI/(k2*xHI + k5*xHII + k13*xe)
    ! three-body formation [cm6 s-1] x density -> [cm3 s-1]
    R31 = 6d-32*T**(-0.25) + 2d-31*T**(-0.5)
    R32 = R31/8d0
    comp_Alpha_H2_Katz = comp_Alpha_H2_Katz + R31*nHI*xHI + R32*nH2*xHI
  end function comp_Alpha_H2_Katz

end module coolrates_module
