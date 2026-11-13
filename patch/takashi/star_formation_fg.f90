!################################################################
! star_formation_fg: FG (first-galaxy) star formation of patch/ks (sf_fg=.true.)
!   Pop II: pop2_model 1-3 local (Jeans threshold, Tpop2, div v, phi, loc_vir, pop2_peak_only),
!           4/5 SFC (star-forming cloud, instantaneous burst, He+2019 efficiency; 5 = star by star)
!   Pop III: pop3_model 1-8 (Zpop3, yH2, candidate separation, mdot-dependent masses)
!   non-local gas depletion (nonlocal_depletion, nH_postSF) through ufront (my_utils)
! Ported from patch/ks star_formation.f90 (merge_ks, docs/merge_ks.md). Differences from ks:
!   - called for every level but acts at nlevelmax only (as ks); the owner->virtual exchange of
!     uold(1:4) started by amr_step is completed here and the virtual cells converted to velocities
!     (ks used the conservative values of the virtual neighbours in div v / the 27-cell dispersion)
!   - particle families: FAM_POP3/SN/DEAD/BH/SFC/PSC of pm_commons (same numbers as ks)
!   - log: per-rank stars_XXXXX.out through sf_log (write_log_star_formation removed); logSFC as ks
!   - ugvtk outputs: output_analysis is a stub; test problems (is_test) and fof removed
!   - mp0/msp0 set under INIT_STELLAR_MASS / STELLAR_POPULATION_MASS (ks: star_by_star.or.fire3_yields)
!   - getnbor / get_27_neighbors are those of star_formation.f90
!################################################################
subroutine star_formation_fg(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X
  use constants, only: rhoc, mH, twopi, M_sun
  use random
  use mpi_mod
  use sf_log
#ifdef RT
  use output_analysis, only: output_around_SFC,output_around_star
  use poisson_commons, only : phi
  use rt_parameters, ONLY: iIons,ixHI,ixHII,ixHeII,ixHeIII
  use sbs_imf
  use my_utils
#endif

  implicit none
  !integer::ilevel
  integer, intent(IN)::ilevel
  integer::ilun_sf=10        ! unit of the per-rank star formation log (sf_log)
  real(dp),allocatable,save :: dep_sumf(:)   ! batched non-local depletion: summed depletion factor per cell (merge_ks)
  integer,allocatable :: pr_star(:), pr_cell(:), pr_lev(:)   ! batched non-local depletion: (star, cell, level) pairs
  integer :: npair=0, npair_max=0

  !----------------------------------------------------------------------
  ! Description: This subroutine spawns star-particle of constant mass
  ! using a Poisson probability law if some gas condition are fulfilled.
  ! It modifies hydrodynamic variables according to mass conservation
  ! and assumes an isothermal transformation...
  ! On exit, the gas velocity and sound speed are unchanged.
  ! New star particles are synchronized with other collisionless particles.
  ! Array flag2 is used as temporary work space.
  ! Yann Rasera  10/2002-01/2003
  !----------------------------------------------------------------------
  ! local constants
  real(dp)::t0,d0,d00,e0,mgas,mcell
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:twotondim,1:3)::xc
  ! other variables
  integer ::ncache,nnew,ivar,irad,ngrid,icpu,index_star,ndebris_tot
  integer ::igrid,ix,iy,iz,ind,i,j,n,iskip,istar,inew,nx_loc
  integer ::ntot,ntot_all,info,nstar_corrected,ideb,ndeb
  logical ::ok_free,ok_all
  real(dp)::d,x,y,z,u,v,w,e,zg,vdisp,dgas
  real(dp)::mstar,dstar,tstar,nISM,nCOM
  real(dp)::T2,nH,T_poly
  real(dp)::velc,uc,vc,wc,mass_load
  real(dp)::vxgauss,vygauss,vzgauss,birth_epoch
  real(kind=8)::mlost,mtot,mlost_all,mtot_all
  real(kind=8)::RandNum,GaussNum,PoissMean
  !real(dp)::vsn,costheta,sintheta,phi,cosphi,sinphi !KS COMMENT OUT
  real(dp),dimension(1:3)::skip_loc
  real(dp)::dx,dx_loc,scale,vol_loc,dx_min,vol_min
  real(dp)::mdebris,vdebris,zdebris,rdebris
  real(dp)::bx1,bx2,by1,by2,bz1,bz2
  integer ,dimension(1:ncpu,1:IRandNumSize)::allseed
  integer ,dimension(1:nvector),save::ind_grid,ind_cell,nstar
  integer ,dimension(1:nvector),save::ind_grid_new,ind_cell_new,ind_part
  integer ,dimension(1:nvector),save::list_debris,ind_debris
  logical ,dimension(1:nvector),save::ok,ok_new=.true.,ok_true=.true.
  logical ,dimension(1:nvector)::ok_sfc ! TO ADDED 2026-11-11: SFC marker candidate cell (pop2_SFC_marker: the SFC-model condition in a local Pop II model)
  logical :: ok_pre_divv, is_peak, vir_ok ! TO ADDED 2026-11-11
  integer ,dimension(1:ncpu)::ntot_star_cpu,ntot_star_all

  !! JP
  logical,dimension(1:nvector)::ok3 ! logical variable for Pop III star
  real(dp)::scale_msun
  logical,dimension(1:nvector)::ok3_new !! distinguish Pop III and Pop I
  integer,dimension(1:nvector)::nstar3 !! number of Pop III star
  ! Note that flag3 is introduced
  real(dp)::overd !! over density
  integer::ipart,jpart
  real(dp)::current_time,birth_time
  real(dp)::factG
  logical::ok_overlap ! check if halo needs to be considered for IMF
  integer::nhalo_overlap ! number of halos to be considered

  ! IMF
  integer,dimension(:),allocatable::halo_overlap ! their indices
  integer::hind ! halo index
  real(dp),dimension(1:3)::pos_prd
  real(dp),dimension(1:nvector)::mpop3imf ! pop3 star mass based on halo mass / detach imf later (test)
  real(dp)::mvir,zcurr ! to use Hirano+ 2014
  logical,dimension(:),allocatable::halo_SF ! check if a halo has another SF site
  logical,dimension(:),allocatable::halo_pop3 ! check if a halo has another Pop III star

  ! Suppress multiple Pop III star formation in one halo
  real(dp),dimension(:,:),allocatable::xp3
  integer::np3
  real(dp)::dist,xx,yy,zz,rsup

  !KS
  integer,save :: ifirst = 0 !KS ADDED (for prcedures at first call)
  real(dp),parameter :: Xhydrogen = 0.76  !(KS MEMO) hydrogen mass fraction (DO NOT CHANGE X (or XH) using namelist file)
  real(dp),parameter :: yHe = 0.25d0*(1d0-Xhydrogen)/Xhydrogen  !(KS MEMO) ratio of helium nuclei and hydrogen nuceli
  real(dp),parameter :: tff0_sec = sqrt(3.*twopi/2./(32.*mH/Xhydrogen*6.67d-8)) !free-fall time [sec] for gas with d*scale_nH = 1cm-3
  real(dp) :: scale_pc
  real(dp) :: Temp, yHI, yHII, yH2, yHeI, yHeII, yHeIII, yel, xmu
  real(dp) :: cI2,njeans_lmax,lamJ_min,dpop2
  real(dp),parameter :: yH2pop3 = 1d-10 !threshold yH2 for pop3 formation (yH2>1e-5 in Ricotti+ 2016)
  !real(dp),parameter :: Tpop2 = 1d4 !threshold T; T<1e4 in Ricotti+ 2016

  !for calculation of divergence
  integer,dimension(1:nvector),save::ind_cell2
  real(dp)::d1,d2,d3,d4,d5,d6
  ! TO for debugging
  real(dp)::ul,ur,divv
  integer :: ncell
  integer ,dimension(1:nvector,0:twondim)::ind_nbor
  integer ,dimension(threetondim)::ind_nborcube

  !for checking potential local minimum
  real(dp)::phi1,phi2,phi3,phi4,phi5,phi6,phimin_nbor

  !for checking local density maximum (TO ADDED)
  real(dp)::dmax_nbor

  !for calculation of the local virial parameter (TO ADDED)
  real(dp) :: vx, vy, vz, vx2, vy2, vz2, vdisp2, m_in_cube, alpha_vir

  !variables for suprresion of simultaneous pop3 formation by distance
  !  integer,parameter :: Npop3_max = 1000 !maximum number of pop3 stars including candidates
  integer,parameter :: Npop3_max = 10000 !maximum number of pop3 stars including candidates
  integer,dimension(1:Npop3_max) :: p3cand_indcell_loc, p3cand_indcell_glob !indcell of pop3 candidates (local/global)
  integer :: p3cand_ncand_loc, p3cand_ncand_tot !number of pop3 candidates in each/all cpu(s)
  integer,dimension(1:ncpu) :: p3cand_ncand_list !list for number of pop3 candidates in each cpu
  real(dp),dimension(1:Npop3_max,1:3) :: p3cand_pos_loc, p3cand_pos_glob !position of pop3 candidates
  integer :: iloc, iglob, jglob
  integer :: p3cand_offset !offset between indices of local and global lists
  real(dp) :: r2, pop3_dist_min
  integer :: p3new_nform_loc, p3new_nform_tot !number of new pop3 formation in each/all cpu(s)
  integer :: p3new_nstar_loc, p3new_nstar_tot !number of new pop3 stars (count individual multiple
  !stars in single star formation event ) in each/all cpu(s)
  integer,dimension(1:ncpu) :: p3new_nform_list !list for number of pop3 formation in each cpu
  integer,dimension(1:ncpu) :: p3new_nstar_list !list for number of pop3 stars in each cpu
  real(dp),dimension(1:3,1:Npop3_max) :: p3new_pos_loc, p3new_pos_glob !position of new pop3
  real(dp),dimension(1:Npop3_max) :: p3new_Mpop3 !mass of new pop3
  !real(dp),parameter :: rsup=2d1 !! physical pc
  !real(dp),parameter :: pop3_dist_min_pc=1d2 !! physical pc (minimum distance) <- now defined in amr_prameters.f90

  !variables for suprresion of simultaneous SFC formation by distance
  integer,parameter :: SFC_n_max = 10000 !maximum number of SFC stars including candidates
  integer,dimension(1:SFC_n_max) :: SFC_indcellcand_loc, SFC_indcellcand_glob !indcell of SFC candidates (local/global)
  integer :: SFC_ncand_loc, SFC_ncand_tot !number of SFC candidates in each/all cpu(s)
  real(dp),dimension(1:3,1:SFC_n_max) :: SFC_poscand_loc, SFC_poscand_glob !position of SFC candidates
  integer :: SFC_n_loc, SFC_n_tot !number of SFCs in each/all cpu(s)
  real(dp),dimension(1:3,1:SFC_n_max) :: SFC_pos_loc, SFC_pos_glob !position of SFCs
  real(dp) ::  SFC_dist_min, SFC_rsize_min
  integer :: iSFC
  integer,dimension(1:ncpu) :: SFC_nstar_list !list for number of SFC stars in each cpu
  integer,dimension(1:ncpu) :: SFC_nSFC_list !list for number of SFCs in each cpu
  integer :: SFC_nSFC_loc !local number of SFCs in each cpu


  !variables for suppression of pop3 formation near pre-existing pop3
  integer :: p3_npop3_loc, p3_npop3_tot
  integer,dimension(1:ncpu) :: p3_npop3_list !list for number of pop3 in each cpu
  real(dp),dimension(1:Npop3_max,1:3) :: p3_pos_loc, p3_pos_glob !position of pop3
  integer :: p3_offset !offset between indices of local and global lists

  !for non-local mass distribution
  logical :: ufront_updated !true if ufront has been updated at the current call of star_formation
  integer :: sweep_level

  !  integer,parameter :: Nnew_max = 10000 !maximum number of newly formed stars each step, changed from 1000 to 10000 -FABG
  !integer,parameter :: Nnew_max = 10000 !maximum number of newly formed stars each step, changed to 10000 for a star-by-star model
  integer,parameter :: Nnew_max = 100000 !maximum number of newly formed stars each step, changed to 100000 for Mstar=10Msun SFC model
  !integer,parameter :: Nnew_max = 300000 !300000 is about the maximum for my Mac environment (for test case)



  !variables for output late SFC
  integer,parameter :: SFCoutput_n_max = 1000 !maximum number of SFC output
  integer :: SFCoutput_n_loc, SFCoutput_n_glob !number of SFC for output
  real(dp),dimension(1:3,1:SFCoutput_n_max) :: SFCoutput_pos_loc, SFCoutput_pos_glob !position of SFC for output
  integer,dimension(1:ncpu) :: SFCoutput_n_list !list for number of SFC for output in each cpu

  !ID of star particle for ugvtk output
  !integer,parameter :: idp_ugvtk = 1 !<- now this parameter is set with namelist


  !for Pop III IMF
  real(dp) ::  Mass, delta_iM, rand
  integer :: iM


  ! !for IMF
  ! real(dp),save,dimension(:),allocatable :: logM_IMF, Pcum_IMF
  ! integer, save:: nmax_IMF




!(KS ADDED) for DM-only simulation
#ifdef RT

  !some of the following lines are commented out for test problems (KS MODIFIED)
  if(numbtot(1,ilevel)==0) return
  !if(.not. hydro)return    !KS COMMENT OUT
  if(ndim.ne.3)return
  !if(static)return         !KS COMMENT OUT

  if(verbose)write(*,*)' Entering star_formation'

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  !! JP
  scale_msun=scale_d*scale_l**3/m_sun
  if(use_proper_time)then
    current_time=texp
  else
    current_time=t
  endif
  zcurr=1d0/aexp-1d0 ! to use Hirano+ 2014
  flag3=0
  !! KS
  scale_pc=scale_l/3.08568d18

  ! Suppress multiple Pop III stars
  rsup=2d1 !! physical pc
  rsup=rsup*3.08568d18/scale_l !! to code unit (merge_ks: no /boxlen)

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)           !(KS COMMENT) skip_loc is always (0,0,0) in our runs
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)

  scale=boxlen/dble(nx_loc) !(KS COMMENT) scale is always 1 in our runsg
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim

  ! Star formation time scale from Gyr to code units
  ! SFR apply here for long lived stars only

  ! ISM density threshold from H/cc to code units
  nISM = n_star
  ! (KS MEMO) 宇宙の平均密度の200倍を最低限要請 (very high-z での星形成を抑制)
  nCOM = del_star*omega_b*rhoc*(h0/100.)**2/aexp**3*XH/mH
  nISM = MAX(nCOM,nISM)
  d0   = nISM/scale_nH
  d00  = n_star/scale_nH

  !**** WARNING WARNING WARNING (KS MODIFIED) WARNING WARNING WARNING ******!
  ! Initial star particle mass
  ! if(m_star < 0d0)then
  !    mstar=n_star/(scale_nH*aexp**3)*vol_min
  ! else
  !    mstar=m_star*mass_sph
  ! endif
  mstar = m_star/scale_msun
  !**** WARNING WARNING WARNING (KS MODIFIED) WARNING WARNING WARNING ******!

  dstar=mstar/vol_loc

  factG = 1d0
  if(cosmo) factG = 3d0/4d0/twopi*omega_m*aexp

  ! Birth epoch as proper time
  if(use_proper_time)then
    birth_epoch=texp
  else
    birth_epoch=t
  endif

  ! Cells center position relative to grid center position
  do ind=1,twotondim
    iz=(ind-1)/4
    iy=(ind-1-4*iz)/2
    ix=(ind-1-2*iy-4*iz)
    xc(ind,1)=(dble(ix)-0.5D0)*dx
    xc(ind,2)=(dble(iy)-0.5D0)*dx
    xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! If necessary, initialize random number generator
  if(localseed(1)==-1)then
    call rans(ncpu,iseed,allseed)
    localseed=allseed(myid,1:IRandNumSize)
  end if

  ! ------------------------- KS ADDED (begin)------------------------!
  !minimum distance of pop3 stars at formation ()
  pop3_dist_min = pop3_dist_min_pc * 3.08568d18/scale_l !! to code unit

  !initialize npop3/ncand count
  p3cand_ncand_loc = 0
  p3_npop3_loc = 0

  !initialize SFC related variables
  SFC_dist_min = SFC_dist_min_pc * 3.08568d18/scale_l !! to code unit
  SFC_rsize_min = SFC_rsize_min_pc * 3.08568d18/scale_l !! to code unit
  SFC_ncand_loc = 0
  SFC_n_loc = 0



  !other initialization (used to judge if star formation occur at the end this subroutine)
  ntot_star_all(:) = 0
  SFC_nstar_list(:) = 0


  !state of ufront
  ufront_updated = .false.
  ! -------------------------- KS ADDED (end) -------------------------!

  !------------------------------------------------
  ! print info at initialization (KS ADDED)
  !------------------------------------------------
  if (ifirst == 0) then
    ! TO ADDED 2026-11-11: the SFC marker is for the local Pop II models only (all ranks must agree)
    if (pop2_SFC_marker .and. (pop2_model == 4 .or. pop2_model == 5)) pop2_SFC_marker = .false.
    if(myid==1)then
      write(*,'(A)') '##################################################################'
      write(*,'(A)') 'Parameters for star formation (first call of star_formation)'
      write(*,'(A)') '##################################################################'

      !check the values of skip_loc and scale
      if (skip_loc(1) /= 0 .or. skip_loc(2) /= 0 .or. skip_loc(1) /= 0 .or. scale /= 1.) then
        ! write (*, '(A)') 'The code assumes skip_loc(:) = 0 and scale = 1, but different value(s) are detected\nstopping...'
        ! stop
        write (*, '(A,/,A)') '*** WARNING *** skip_loc(:) /= 0 or scale /= 1 (boxlen /= nx): the SFC / non-local routines were written for boxlen=1 (merge_ks: cell volumes and the boundary checks are scaled)','keep running...'
      end if

      !nonlocal_depletion
      if (nonlocal_depletion) then
        write(*,'(A,/,A)') 'For both POPII/III',&
          'nonlocal_depletion = true. --> gas is depleted from a spherical region (not just the cell hosting a new star) at star formation'
        !density correction after star fomration
        if (nH_postSF <= 0) then
          write(*,'(A,1PE11.3,A)') 'nH_postSF = ', nH_postSF, ' (< 0) --> do not fix density of post-star-forming cells'
        else
          write(*,'(A,1PE11.3,A)') 'nH_postSF = ', nH_postSF, ' cm-3 (> 0) --> fix density of post-star-forming cells'
        end if

      end if

      !----------------------
      ! Model options
      !----------------------
      write(*,'(A)') '--------------------------------------------------------------------------------------------'

      ! POP II MODEL
      write(*,'(A, I4)') 'pop2_model = ', pop2_model
      if (pop2_model == 1) then
        write(*,'(A)') 'threshold density/temperature is used for Pop II formation condition'
      else if (pop2_model == 2) then
        write(*,'(A)') 'threshold density/temperature and div V is used for Pop II formation condition  *** WARNING NOT WELL TESTED ***'
      else if (pop2_model == 3) then
        write(*,'(A)') 'threshold density/temperature, div V, local pot min is used for Pop II formation condition  *** WARNING NOT WELL TESTED ***'
      else if (pop2_model == 4) then
        write(*,'(A)') 'instantaneous-burst star-forming-cloud model is used for Pop II formation  *** WARNING NOT WELL TESTED ***'
      else if (pop2_model == 5) then
        write(*,'(A)') 'STAR-BY-STAR instantaneous-burst star-forming-cloud model is used for Pop II formation  *** WARNING NOT WELL TESTED ***'
        ! TO ADDED: pop2_model=5 samples the IMF table, which is initialized only if star_by_star (2026-10-30)
        if (.not. star_by_star) then
          write(*,'(A)') 'pop2_model = 5 requires star_by_star = .true. (IMF table is not initialized). stopping...'
          stop
        end if
      else
        write(*,'(A)') 'current code is only compatible with pop2_model = 1 or 2. stopping...'
        stop
      end if

      ! TO ADDED 2026-11-11: SFC markers in the local Pop II models
      if (pop2_SFC_marker) then
        write(*,'(A)') 'pop2_SFC_marker = .true.  -->  a massless FAM_SFC marker is created where the SFC formation condition (density peak, loc_vir, no div v) is met and no active SFC marker lies within SFC_dist_min_pc; formation nest and age outputs as for an SFC'
      else if (pop2_model == 4 .or. pop2_model == 5) then
        write(*,'(A)') 'pop2_SFC_marker = .false. (the SFC model has SFC particles anyway)'
      else
        write(*,'(A)') 'pop2_SFC_marker = .false. -->  no SFC marker particles (no SFC-like ugvtk outputs in this local Pop II model)'
      end if

      ! Feedback model for Pop II stars (star_by_star or not)
      if (star_by_star) then
        write(*,'(A)') 'star_by_star = .true.  -->  radiation/SN is treated assuming each Pop II particle is an individual star'
      else
        write(*,'(A)') 'star_by_star = .false.  -->  radiation/SN is treated assuming each Pop II particle is a group of stars with some population'
      end if

      write(*,'(A)') '--------------------------------------------------------------------------------------------'


      ! POP III MODEL
      write(*,'(A, I4)') 'pop3_model = ', pop3_model
      if (pop3_model == 1) then
        write(*,'(A)') 'always 80+40 Msun binary Pop III stars'
      else if (pop3_model == 2) then
        write(*,'(A)') '1000 Msun single Pop III star only for the first time, and 80+40 Msun binary Pop III stars thereafter'
      else if (pop3_model == 3) then
        write(*,'(A,1PE11.3,A)') 'single Pop III star with mass specified with the namelist parameter mpop3max = ', mpop3max, ' Msun'
        if (mpop3_SNprog .ge. 0)  then
          write(*,'(A,1PE11.3,A)') 'Pop III stellar mass for SN progenitor is set independently to mpop3_SNprog = ', mpop3_SNprog, ' Msun'
        end if
      else if (pop3_model == 4) then
        write(*,'(A,1PE11.3,A)') 'single Pop III star with mass determined by cloud-scale accretion rate'
      else if (pop3_model == 5) then
        write(*,'(A,1PE11.3,A)') 'M1:M2 = 2:1 binary Pop III stars with total mass determined by cloud-scale accretion rate'
      else if (pop3_model == 6) then
        write(*,'(A,1PE11.3,A)') 'M1:M2:..:Mn = 2^n-1:2^n-2:...:1 multiplet Pop III stars with total mass determined by cloud-scale accretion rate (with maximum n that satisfies Mn > 10 Msun)'
      else if (pop3_model == 7) then
        write(*,'(A,1PE11.3,A)') 'Single Pop III with random mass from Salpeter IMF (10-500 Msun)'
      else if (pop3_model == 8) then
        write(*,'(A,1PE11.3,A)') 'Single Pop III with random mass from Flat IMF (10-500 Msun)'
      else
        write(*,'(A)') 'current code is only compatible with pop3_model = 1, 2, 3, or 4 for non-test runs. stopping...'
        stop
      end if

      write(*,'(A)') '--------------------------------------------------------------------------------------------'

      !----------------------------------------------
      ! info about Pop II parameters
      !----------------------------------------------
      if (pop2_njeans > 0) then         !ジーンズ長を用いてthreshold densityを決める場合
        write(*,'(A)') &
          'pop2_njeans > 0  ==> use the ratio of Jeans length and minimum cell size for threshold density of Pop II formation'

        njeans_lmax = jeans_refine(nlevelmax-1)
        if (pop2_njeans > 2*njeans_lmax) then
          write(*,'(2(A,1PE11.3),/,A)') &
            '**WARNING** pop2_njeans = ', pop2_njeans, &
            'is more than two times larger than njeans_lmax = ', njeans_lmax, &
            'This results in Pop II threshold density smaller than the density needed for maximum refinemnt (star formation occurs as soon as the maximum level appears)'
        end if


        cI2 = 1d0 /scale_T2 !T2=1Kに対応する等温音速
        lamJ_min = boxlen / (2.d0)**nlevelmax * pop2_njeans !ジーンズ長がこれ以上短い場合に星形成
        dpop2 = twopi/2.d0*cI2/factG/lamJ_min**2 !lamJ = lamJ_min となる密度 (1K, current aexp))

        !print *, "lamJ@nH=1cm-3, T2=1K [pc]:  ",sqrt((twopi/2.*1d0/scale_T2) / (factG*1d0/scale_nH)) * scale_l/3.08568d18
        write(*,'(3(A, 1PE11.3))') &
          'dx_min [cpc]: ', boxlen / (2.d0)**nlevelmax * scale_l/3.08568d18 / aexp, &
          ', njeans_lmax: ', njeans_lmax, ', pop2_njeans: ', pop2_njeans

        write(*,'(8(A, 1PE11.3))') &
          '** POPII ** dpop2 [cm-3 * ((1+z)/10)^2 * (T2/100K) ]: ',&   !dpop2はTと(1+z)^2に比例
          dpop2*scale_nH * 1d2 * (aexp/0.1)**2,&
          ', Tpop2: ', Tpop2,&         !threshold temperature
          ', eps_star: ', eps_star,&   !t_star = tff/eta_star
          ', tff@nH=1cm-3 [Myr]: ',tff0_sec/(1d6*365.*24.*3600.),&
          ', mstar[Msun]',mstar*scale_msun

      else               !threshold densityを直接決める場合
        if (n_star < 0) then
          write(*,'(2(A,1PE11.3),/,A)') &
            '**ERROR** pop2_njeans = ', pop2_njeans, &
            'and n_star = ', n_star, &
            'n_star should be positive if pop2_njeans is negative. stopping...'
          stop
        end if

        dpop2 = d00 !d0とd00があるが、ひとまずd00を使うようにしておく
        write(*,'(10(A, 1PE11.3))') &
          '** POPII ** d0 [cm-3]: ',d0*scale_nH,&
          ', dpop2 (= d00) [cm-3]: ',d00*scale_nH,&
          ', Tpop2: ', Tpop2,&         !threshold temperature
          ', eps_star: ', eps_star,&   !t_star = tff/eps_star
          ', tff@nH=1cm-3 [Myr]: ',tff0_sec/(1d6*365.*24.*3600.),&
          ', mstar[Msun]',mstar*scale_msun
      end if

      ! TO ADDED 2026-11-11: virial (velocity dispersion) condition and local peak condition for Pop II
      ! (both acted on the Pop III flag before this date, i.e. never)
      if (loc_vir > 0) then
        if (pop2_njeans > 0) then
          write(*,'(A,1PE11.3,A)') 'loc_vir = ', loc_vir, &
            '  ==> Pop II forms only where (cs^2+sigma^2)/(G rho lamJ_min^2) * pi <= loc_vir (Jeans condition incl. velocity dispersion of the 27-cell cube)'
        else
          write(*,'(A,1PE11.3,A)') 'loc_vir = ', loc_vir, &
            '  ==> Pop II forms only where alpha_vir(27-cell cube) <= loc_vir'
        end if
      else
        write(*,'(A)') 'loc_vir <= 0  ==> no virial condition for Pop II formation'
      end if
      if (pop2_peak_only) then
        write(*,'(A)') 'pop2_peak_only = .true.  ==> Pop II forms only in cells that are the density maximum of their 27-cell cube'
      else
        write(*,'(A)') 'pop2_peak_only = .false. ==> no local density peak condition for Pop II formation'
      end if

      !H2-dependent Pop II formation
      if (pop2_H2depSF) then
        write(*,'(A)') 'pop2_H2depSF == .true.  -->  Pop II formation rate is proportional to rho(H2) instead of rho'
      end if


      !SFC PopII formation model
      !if (pop2_model == 4) then
      if (pop2_model == 4 .or. pop2_model == 5) then

        !warning that eps_star not used in SFC model
        write(*,'(A)') 'SFC Pop II star formation  -->  SFC forms as soon as the formation condition is satisfied (eps_star is ignored)'


        ! if ((SFC_nH_min >= 0 .and. SFC_f_nH_min >= 0) .or. (SFC_nH_min <= 0 .and. SFC_f_nH_min <= 0)) then
        !   write(*,'(2(A,1PE11.3),/,A,/,A)') 'Either SFC_nH_min or SFC_f_nH_min should be positive, but SFC_nH_min = ', &
        !     SFC_nH_min, ',  SFC_f_nH_min = ', SFC_f_nH_min,&
        !     'Check namelist file',&
        !     'stopping...'
        !   stop
        ! else if (SFC_nH_min > 0) then
        !   write(*,'(A,1PE11.3,A)') 'SFC_nH_min = ', SFC_nH_min, ' cm-3  -> use fixed threshold density for SFC outer boundary'
        ! else if (SFC_f_nH_min > 0) then
        !   write(*,'(A,1PE11.3,A)') 'SFC_f_nH_min = ', SFC_f_nH_min, ' -> use threshold density proportional to the central value for SFC outer boundary'
        ! end if

        if (SFC_nH_min <= 0 .and. SFC_f_nH_min <= 0) then
          write(*,'(2(A,1PE11.3),/,A,/,A)') 'Either SFC_nH_min or SFC_f_nH_min should be positive, but SFC_nH_min = ', &
            SFC_nH_min, ',  SFC_f_nH_min = ', SFC_f_nH_min,&
            'Check namelist file',&
            'stopping...'
          stop
        else if (SFC_nH_min > 0 .and. SFC_f_nH_min > 0) then
          write(*,'(2(A,1PE11.3))') 'Both SFC_nH_min and SFC_f_nH_min are positive: use the larger threshold density determined by SFC_nH_min = ', &
            SFC_nH_min, ' and  SFC_f_nH_min = ', SFC_f_nH_min
        else if (SFC_nH_min > 0) then
          write(*,'(A,1PE11.3,A)') 'SFC_nH_min = ', SFC_nH_min, ' cm-3  -> use fixed threshold density for SFC outer boundary'
        else if (SFC_f_nH_min > 0) then
          write(*,'(A,1PE11.3,A)') 'SFC_f_nH_min = ', SFC_f_nH_min, ' -> use threshold density proportional to the central value for SFC outer boundary'
        end if


        !recipe for star formation efficiency
        if (SFC_f_star > 0) then
          write(*,'(A,1PE11.3,A)') 'SFC_f_star = ', SFC_f_star, '  -> use fixed star formation efficiency'
        else
          write(*,'(A)') 'SFC_f_star is negative  -> use He, Ricotti, and Geen (2019) formula for star formation efficiency'
          !TO ADDED 2026-11-10: metallicity dependence
          select case (SFC_f_star_Zdep)
          case (0)
            write(*,'(A)') '  with the fixed low-metallicity factor 0.2 (He et al. 2019: 1/5 for Z < 0.1 Zsun; = Sugimura et al. 2024 eq. 2 at Z = 1e-3 Zsun)'
          case (1)
            write(*,'(A,1PE10.2,A)') '  with 0.2 (Z/1e-3 Zsun)^0.25 (Sugimura et al. 2024, eq. 2), Z floored at Zpop3 = ', Zpop3, ' Zsun'
          case (2)
            write(*,'(A)') '  with max(0.2, min(1, Z/Zsun)) (He et al. 2019 sec. 3.4: linear from solar to 0.2 Zsun, saturated below)'
          case default
            write(*,'(A,I3,A)') 'SFC_f_star_Zdep = ', SFC_f_star_Zdep, ' is not 0, 1 or 2\nstopping...'
            call clean_stop
          end select
        end if

        !optically thin assumption for region around  SFC
        if (SFC_flg_thin) then
          write(*,'(A)') &
            'SFC_flg_thin = .true.  ->  assume SFC sphere with r < R_SFC to be optically thin'
        else
          write(*,'(A)') &
            'SFC_flg_thin = .false.  ->  SFC particle is a test particle'
        end if
        write(*,'(A,1PE11.3,A)') 'SFC_t_life = ', SFC_t_life, ' Myr'

      end if




      !----------------------------------------------
      ! info about Pop III parameters
      !----------------------------------------------
      write(*,'(5(A, 1PE11.3))') &
        '** POPIII ** dpop3 [cm-3]: ',dpop3,&
        ', Zpop3 [Zsun]', Zpop3,&
        ', yH2pop3 [Zsun]', yH2pop3,&
        ', pop3_dist_min [pc]', pop3_dist_min_pc
      if (dpop3 < 0d0) then
        write(*,'(A)') 'dpop3 < 0 -> use the same density condition based on the Jeans length as for Pop II'
      end if

      write(*,'(A)') 'Set minimum Pop II particle mass using namelist parameter m_star (KS MODIFIED)'
      write(*,'(A)') 'Pop II formation limited to maximum level (KS MODIFIED)'
      write(*,'(A)') 'Pop III formation independently of host halo (KS DEBUG)'
      write(*,'(A)') 'Pop III formation suppression by pre-existing pop3 modified (KS MODIFIED)'
      write(*,'(A)') 'Pop III formation suppression by other candidates modified (KS MODIFIED)'


      write(*,'(A)') '##################################################################'

    end if !myid == 1

    !initializing IMF (set table for Pcum and so on)
    if (star_by_star) then
      call initialize_IMF
    end if

    !initializing Pop III IMF (set table for Pcum and so on)
    if (pop3_model == 7 .or. pop3_model == 8) then
      call initialize_Pop3_IMF
    end if


    ifirst = 1
  end if

  !output summary if ilevel == levelmin (KS ADDED)
  if (ilevel == levelmin) call star_formation_summary()

  !output SFC at given age (KS ADDED) <- THIS SUBROUTINE IS NOW CALLED IN output_analysis_data
  !call output_around_SFC_at_given_age(ilevel)

  !output star at given age (KS ADDED) <- THIS SUBROUTINE IS NOW CALLED IN output_analysis_data
  !call output_around_star_at_given_age(ilevel)

  !leave star_formation if ilevel != nlevelmax (KS ADDED)
  if (ilevel /= nlevelmax) return

  call sf_log_open(ilun_sf)

  !get list for pre-existing pop3 particles in all cpus (KS ADDED)
  call get_pop3_pos_list()

#if NDIM==3
  !------------------------------------------------
  ! Convert hydro variables to primitive variables
  !------------------------------------------------
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
    ngrid=MIN(nvector,ncache-igrid+1)
    do i=1,ngrid
      ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
    end do
    do ind=1,twotondim
      iskip=ncoarse+(ind-1)*ngridmax
      do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
      end do
      do i=1,ngrid
        d=uold(ind_cell(i),1)
        u=uold(ind_cell(i),2)/d
        v=uold(ind_cell(i),3)/d
        w=uold(ind_cell(i),4)/d
        e=uold(ind_cell(i),5)
#ifdef SOLVERmhd
        bx1=uold(ind_cell(i),6)
        by1=uold(ind_cell(i),7)
        bz1=uold(ind_cell(i),8)
        bx2=uold(ind_cell(i),nvar+1)
        by2=uold(ind_cell(i),nvar+2)
        bz2=uold(ind_cell(i),nvar+3)
        e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)
#endif
        e=e-0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
        do irad=1,nener
          e=e-uold(ind_cell(i),5+irad)
        end do
#endif
        uold(ind_cell(i),1)=d
        uold(ind_cell(i),2)=u
        uold(ind_cell(i),3)=v
        uold(ind_cell(i),4)=w
        uold(ind_cell(i),5)=e/d
      end do
      do ivar=imetal,nvar
        do i=1,ngrid
          d=uold(ind_cell(i),1)
          w=uold(ind_cell(i),ivar)/d
          uold(ind_cell(i),ivar)=w
        end do
      end do
    end do
  end do

  ! Density and velocities of the virtual cells (needed by getnbor / get_27_neighbors at domain
  ! boundaries): complete the owner -> virtual exchange started by amr_step (slot 1) and convert
  ! them like the active cells (ks used the conservative values of the virtual neighbours here).
#ifndef WITHOUTMPI
  if(async_fwd(1)%active)then
#ifdef SOLVERmhd
    call make_virtual_fine_dp_end(uold,nvar+3,1,4,ilevel,1)
#else
    call make_virtual_fine_dp_end(uold,nvar,1,4,ilevel,1)
#endif
  else
#ifdef SOLVERmhd
    call make_virtual_fine_dp_nv(uold,nvar+3,1,4,ilevel)
#else
    call make_virtual_fine_dp_nv(uold,nvar,1,4,ilevel)
#endif
  end if
  do icpu=1,ncpu
    do ind=1,twotondim
      iskip=ncoarse+(ind-1)*ngridmax
      do i=1,reception(icpu,ilevel)%ngrid
#ifdef LIGHT_MPI_COMM
        ncell=reception(icpu,ilevel)%pcomm%igrid(i)+iskip
#else
        ncell=reception(icpu,ilevel)%igrid(i)+iskip
#endif
        d=uold(ncell,1)
        uold(ncell,2)=uold(ncell,2)/d
        uold(ncell,3)=uold(ncell,3)/d
        uold(ncell,4)=uold(ncell,4)/d
      end do
    end do
  end do
#endif

  !**** WARNING WARNING WARNING (KS MODIFIED) WARNING WARNING WARNING ******!
  ! ! Find existing Pop III stars
  ! if(ilevel.eq.nlevelmax)then ! Pop III star forms only at maximum level
  !    np3=0
  !    ! (KS MEMO) 他のcpuが所持するPop IIIの情報は使わなくて大丈夫？
  !    do i=1,npartmax ! Loop over particles in each domain - do not consider ones outside
  !       if( is_pop3(typep(i)) )np3=np3+1
  !    enddo

  !    allocate(xp3(1:np3,1:3))

  !    np3=0
  !    do i=1,npartmax
  !       if( is_pop3(typep(i)) )then
  !          np3=np3+1
  !          xp3(np3,:)=xp(i,:)
  !       endif
  !    enddo

  !    !------------------------------------------------
  !    ! Check if each halo needs to be considered for IMF - JP
  !    !------------------------------------------------
  !    nhalo_overlap=0
  !    do i=1,nhalotot
  !       call check_overlap(halopos(i,1:3),ok_overlap,pos_prd)
  !       if(ok_overlap)nhalo_overlap=nhalo_overlap+1
  !    enddo

  !    allocate(halo_overlap(1:nhalo_overlap))
  !    allocate(halo_SF(1:nhalo_overlap))
  !    allocate(halo_pop3(1:nhalo_overlap))
  !    nhalo_overlap=0
  !    halo_SF=.false.
  !    halo_pop3=.false.
  !    do i=1,nhalotot
  !       call check_overlap(halopos(i,1:3),ok_overlap,pos_prd)
  !       if(ok_overlap)then
  !          nhalo_overlap=nhalo_overlap+1
  !          halo_overlap(nhalo_overlap)=i

  !          ! halos containing Pop III star alive
  !          do j=1,np3
  !             dist=sqrt( (pos_prd(1)-xp3(j,1))**2 + (pos_prd(2)-xp3(j,2))**2 + (pos_prd(3)-xp3(j,3))**2 )
  !             if(dist.le.r200(i))halo_pop3(nhalo_overlap)=.true.
  !          enddo
  !       endif
  !    enddo
  ! endif
  !**** WARNING WARNING WARNING (KS MODIFIED) WARNING WARNING WARNING ******!


  !------------------------------------------------
  ! Compute number of new stars in each cell
  ! set flag2 and flag3
  !------------------------------------------------
  ntot=0
  ndebris_tot=0
  ! Loop over grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
    ngrid=MIN(nvector,ncache-igrid+1)
    do i=1,ngrid
      ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
    end do
    ! Star formation criterion ---> logical array ok(i)
    do ind=1,twotondim
      iskip=ncoarse+(ind-1)*ngridmax
      do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
      end do

      if(ilevel.eq.nlevelmax)mpop3imf=0 ! initialize before using
      do i=1,ngrid

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        !! Pop III. (ok3) Ricotti+ 2016 - JP
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        ! Flag leaf cells - JP
        ok3(i)=son(ind_cell(i))==0

        if(ilevel.ne.nlevelmax)ok3(i)=.false.                     ! Allow only at maximum refinement cell

        !ok3(i)=ilevel.eq.nlevelmax                               ! maximum level! (hence leaf cells)

        if(uold(ind_cell(i),imetal).gt.Zpop3*0.02)ok3(i)=.false. ! Metallicity


        nH=max(uold(ind_cell(i),1),smallr)*scale_nH              ! Density in H/cm3
        if (dpop3 > 0) then  !in the case with independent PopIII threshold density
          if(nH.lt.dpop3)ok3(i)=.false.                            ! Density for clump

        else                 !in the case using the same density threshold as PopII
          d=uold(ind_cell(i),1)
          if (pop2_njeans > 0) then         !ジーンズ長を用いてthreshold densityを決める場合
            cI2 = uold(ind_cell(i),5)*(gamma-1.0) !等温音速 (!上の方ですでにuold(ind_cell(i),1)で割り算されているので注意)
            lamJ_min = boxlen / (2.d0)**nlevelmax * pop2_njeans !ジーンズ長がこれ以上短い場合に星形成
            dpop2 = twopi/2.*cI2/factG/lamJ_min**2 !lamJ = lamJ_min となる密度 (1K, current aexp))
          else                             !固定値を与える場合
            dpop2 = d00 !d0とd00があるが、ひとまずd00を使うようにしておく
          end if
          if(d<=dpop2)ok3(i)=.false.
        end if

        !-------------- KS ADDED (BEGIN) -------------!
        ! Geometrical criterion (Pop IIと同様に、Pop IIIもzoom-in 領域外での星形成は抑制する)
        if(ivar_refine>0)then
          d=uold(ind_cell(i),ivar_refine)
          if(d<=var_cut_refine)ok3(i)=.false.
        endif
        !-------------- KS ADDED (END) -------------!

        !-------------- TO ADDED (BEGIN) -------------!
        if (ok3(i) .and. loc_vir > 0) then

          ! TO: check if the cell is self-gravitating
          ! first we need to compute the local velocity dispersion
          d = uold(ind_cell(i),1)
          cI2 = uold(ind_cell(i),5)*(gamma-1.0) !等温音速 (!上の方ですでにuold(ind_cell(i),1)で割り算されているので注意)
          ncell = 1 ! we just want the neighbors of the current cell
          ind_cell2(1) = ind_cell(i)
          call get_27_neighbors(ind_cell2, ind_nborcube, ilevel)
          vx = 0d0; vy = 0d0; vz = 0d0;
          vx2 = 0d0; vy2 = 0d0; vz2 = 0d0;
          m_in_cube = 0d0
          d = uold(ind_cell(i),1) ! density of the central cell
          dmax_nbor = d
          do j = 1, 27
            d1 = uold(ind_nborcube(j),1)
            dmax_nbor = max(dmax_nbor, d1)
            vx = vx + d1*uold(ind_nborcube(j),2)*vol_loc
            vy = vy + d1*uold(ind_nborcube(j),3)*vol_loc
            vz = vz + d1*uold(ind_nborcube(j),4)*vol_loc
            vx2 = vx2 + d1*vol_loc*(uold(ind_nborcube(j),2))**2
            vy2 = vy2 + d1*vol_loc*(uold(ind_nborcube(j),3))**2
            vz2 = vz2 + d1*vol_loc*(uold(ind_nborcube(j),4))**2
            m_in_cube = m_in_cube + d1*vol_loc
          end do

          !-----------------------------------------------------------------------------------
          if (d < dmax_nbor) then ! check if the cell is a local density peak
            ok3(i) = .false.  ! the cell is not a local density maximum
          else
            vx = vx / m_in_cube
            vy = vy / m_in_cube
            vz = vz / m_in_cube
            vx2 = vx2 / m_in_cube
            vy2 = vy2 / m_in_cube
            vz2 = vz2 / m_in_cube
            vdisp2 = (vx2 - vx**2) + (vy2 - vy**2) + (vz2 - vz**2)
            vdisp2 = max(vdisp2, 0d0)
            ! alpha_vir = 5 * (sigma^2 + c_s^2) * R / (G * M)
            ! here 4 pi / 3 * R^3 = 27 * dx^3  --> R = (3/4 pi * 27 * dV)^(1/3)
            if (dpop3 <= 0 .and. pop2_njeans > 0) then
              lamJ_min = boxlen / (2.d0)**nlevelmax * pop2_njeans !ジーンズ長がこれ以上短い場合に星形成
              dpop2 = twopi/2.*cI2/factG/lamJ_min**2 !lamJ = lamJ_min となる密度 (1K, current aexp))
              alpha_vir = twopi/2.*(cI2 + vdisp2)/(factG*d*(lamJ_min**2))
            else
              alpha_vir = 9.305257d0*(vdisp2 + cI2) * dx_loc / (factG * m_in_cube)
            endif
            ! write(*,*) 'sigma = ', sqrt(vdisp2)*scale_v*1.d-5, ' km/s cs = ', sqrt(cI2)*scale_v*1.d-5, ' km/s a_vir = ', alpha_vir
            if (alpha_vir > loc_vir) ok3(i) = .false.
          end if
        endif
        !-------------- TO ADDED (END) -------------!

        ! --------------------- KS ADDED --------------------- !
        ! H2 fraction check (Ricotti + 2016)
        ! ---------------------------------------------------- !
        yHI = uold(ind_cell(i),iIons-1+ixHI) !上の方ですでにuold(ind_cell(i),1)で割り算されているので注意
        yHII = uold(ind_cell(i),iIons-1+ixHII) !上の方ですでにuold(ind_cell(i),1)で割り算されているので注意
        yH2 = (1.-yHI-yHII)*0.5

        if(yH2 < yH2pop3) ok3(i)=.false. !! yH2>1e-5 in Ricotti+ 2016
        ! ---------------------------------------------------- !
        ! ---------------------------------------------------- !


        ! rejection due to pre-existing Pop III stars
        if(ok3(i))then ! Avoid loop if possible
          do iglob = 1, p3_npop3_tot
            ! Cell position
            xx = (xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
            yy = (xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
            zz = (xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale

            ! Distances to Pop III star particles
            r2 = (xx-p3_pos_glob(iglob,1))**2 + &
              (yy-p3_pos_glob(iglob,2))**2 + &
              (zz-p3_pos_glob(iglob,3))**2
            if(r2 < pop3_dist_min**2) then
              ! !---------------- KS DEBUG ------------------!
              ! print '(A,2(A,3(1PE12.5)),A,1PE12.5)', &
              !      "KS DEBUG (star_formation.f90) rejected due to pre-existing Pop III ", &
              !      ' r_rej(:): ', p3_pos_glob(iglob,:), &
              !      ', r_cell(:): ', xx,yy,zz, &
              !      ', r/rmin: ', sqrt(r2)/pop3_dist_min
              ! !---------------- KS DEBUG ------------------!
              ok3(i)=.false.
              exit
            end if
          end do
          !**** WARNING WARNING WARNING (KS MODIFIED) WARNING WARNING WARNING ******!
        endif

        !mass (KS MODIFIED)
        if(ok3(i))then
          !mpop3imf(i)=100./scale_msun ! (KS MEMO) 100 M_sun決め打ち
          if (pop3_model == 1) then              !POPIII MODEL 1
            mpop3imf(i)=1.2d2/scale_msun ! 80+40=120 Msun (always)
          else if(pop3_model == 2) then          !POPIII MODEL 2
            if (nstar_tot == 0) then
              mpop3imf(i)=1d3/scale_msun ! 1000 Msun (for the first time)
            else
              mpop3imf(i)=1.2d2/scale_msun ! 80+40=120 Msun (otherwise)
            end if
          else if(pop3_model == 3) then          !POPIII MODEL 3
            mpop3imf(i)=mpop3max/scale_msun      !コードに存在していたnamelist変数を流用（maxの意味は無いので注意）
          else if(pop3_model == 4 .or. pop3_model == 5 .or. pop3_model == 6) then          !POPIII MODEL 4, 5 or 6
            mpop3imf(i)=1d-50 !dummy (determine pop3mass later based on cloud-scale Mdot)

          else if (pop3_model == 7 .or. pop3_model == 8 ) then !USE ARTIFICIAL IMF
            call random_value_0to1(rand, initial_seed=nstep)  !get random number between 0 and 1 (pass nstep for initial_seed to avoid the same pattern)
            do iM =1, p3_NM-1
              if (rand <= p3_Pcum(iM+1)) then
                !interpolation
                delta_iM = (rand - p3_Pcum(iM))/(p3_Pcum(iM+1)-p3_Pcum(iM))
                Mass = 10d0**( log10(p3_Mmin_imf) + (iM+delta_iM-1)*(log10(p3_Mmax_imf)-log10(p3_Mmin_imf))/(p3_NM-1) )
                exit
              end if
            end do
            mpop3imf(i)=Mass/scale_msun
            !KS DEBUG
            print '(A,I8, 2(A,1PE12.4))', "(KS DEBUG) i = ", i, ", rand = ", rand, ", M[Msun] = ", Mass

          end if
        end if

        ! !------------------ KS DEBUG ----------------!
        ! if (nH > dpop3 .and. uold(ind_cell(i),imetal) < Zpop3*0.02) then
        !       write(*,'(2(A,L),A,I3, /,3(A, 1P2E10.3), (A, 3(1PE15.8)))') &
        !         'POPIII conditions for d and Z are met. ok3: ', ok3(i), &
        !         ', leaf: ',son(ind_cell(i))==0,&
        !         ', level: ', ilevel,&
        !         'd, dcr [cm-3]:',nH, dpop3,&
        !         ', Z, Zcr [Zsun]:',uold(ind_cell(i),imetal)/0.02, Zpop3,&
        !         ', yH2, yH2cr [K]:',yH2, yH2pop3,&
        !         ', r(:): ', xg(ind_grid(i),1)+xc(ind,1),xg(ind_grid(i),2)+xc(ind,2),xg(ind_grid(i),3)+xc(ind,3)
        ! end if
        ! !------------------ KS DEBUG ----------------!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! --- Pop III end
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Pop II (set ok)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Flag leaf cells - JP
        ok(i)=son(ind_cell(i))==0

        ! Metallicity
        if(uold(ind_cell(i),imetal).le.Zpop3*0.02)ok(i)=.false.

        ! Allow SF only at maximum refinement cell
        if(ilevel.ne.nlevelmax)ok(i)=.false.

        ! Density criterion (KS MODIFIED)
        ! use Jeans length to determine threshold density for Pop II formation
        ! temperature dependent threshold density based on Jeans length (e.g., Ricotti + 2016)
        d=uold(ind_cell(i),1)

        !**** WARNING WARNING WARNING (KS MODIFIED) WARNING WARNING WARNING ******!
        !if(d<=d0)ok(i)=.false.
        !if(d<=d00)ok(i)=.false. ! KS DEBUG
        !**** WARNING WARNING WARNING (KS MODIFIED) WARNING WARNING WARNING ******!

        if (pop2_njeans > 0) then         !ジーンズ長を用いてthreshold densityを決める場合
          cI2 = uold(ind_cell(i),5)*(gamma-1.0) !等温音速 (!上の方ですでにuold(ind_cell(i),1)で割り算されているので注意)
          lamJ_min = boxlen / (2.d0)**nlevelmax * pop2_njeans !ジーンズ長がこれ以上短い場合に星形成
          dpop2 = twopi/2.*cI2/factG/lamJ_min**2 !lamJ = lamJ_min となる密度 (1K, current aexp))
          !if(d<=dpop2)ok(i)=.false.

        else                             !固定値を与える場合
          dpop2 = d00 !d0とd00があるが、ひとまずd00を使うようにしておく
        end if
        if(d<=dpop2)ok(i)=.false.


        ! Temperature criterion (KS MODIFIED)
        ! (KS MEMO) T2の定義が問題無いか要確認 <- 読みにくいが、上の方でuoldを定数倍しているので答えはOK (by Jongwon)
        ! T2=uold(ind_cell(i),5)*scale_T2*(gamma-1.0)
        ! nH=max(uold(ind_cell(i),1),smallr)*scale_nH
        ! T_poly=T2_star*(nH/nISM)**(g_star-1.0)
        ! T2=T2-T_poly
        ! if(T2>1e4)ok(i)=.false. !! T<1e4 in Ricotti+ 2016

        ! Use normal T because it is more intuitive than T2=T/mu
        yHI = uold(ind_cell(i),iIons-1+ixHI) !上の方ですでにuold(ind_cell(i),1)で割り算されているので注意
        yHII = uold(ind_cell(i),iIons-1+ixHII) !上の方ですでにuold(ind_cell(i),1)で割り算されているので注意
        yHeII = uold(ind_cell(i),iIons-1+ixHeII)*yHe !y_i = n(i)/nH
        yHeIII = uold(ind_cell(i),iIons-1+ixHeIII)*yHe !y_i = n(i)/nH
        yH2 = (1d0-yHI-yHII)*0.5d0
        yHeI =  yHe-yHeII-yHeIII
        yel=yHII+yHeII+2*yHeIII
        xmu=(yHI+yHII+2d0*yH2 + 4d0*yHe) / (yHI+yHII+yH2 + yHe + yel) !mean molecular weight
        Temp=uold(ind_cell(i),5)*scale_T2*(gamma-1d0)*xmu
        T2=uold(ind_cell(i),5)*scale_T2*(gamma-1.0)
        ! write(*,'(A, 1P10E9.2)') &
        !      'KS DEBUG (star_formation.f90) nH, Temp, yHI, yHII, yH2, yHeI, yHeII, yHeIII, yel, xmu: ', &
        !      nH, Temp, yHI, yHII, yH2, yHeI, yHeII, yHeIII, yel, xmu

        if(Temp>Tpop2)ok(i)=.false. !! T<1e4 in Ricotti+ 2016


        ! Geometrical criterion
        if(ivar_refine>0)then
          d=uold(ind_cell(i),ivar_refine)
          if(d<=var_cut_refine)ok(i)=.false.
        endif

        ! TO ADDED 2026-11-11: the conditions so far (without div v / potential) are those of the SFC model; kept for the SFC marker
        ok_pre_divv = ok(i)
        ok_sfc(i) = .false.


        !-----------------------------------------------------------------------------------
        ! Compute div v: copy-and-pasted with modification from new star_formation.f90
        !-----------------------------------------------------------------------------------
        ! We need to estimate the norm of the gradient of the velocity field in the cell (tensor of 2nd rank)
        ! i.e. || A ||^2 = trace( A A^T) where A = grad vec(v) is the tensor.
        ! So construct values of velocity field on the 6 faces of the cell using simple linear interpolation
        ! from neighbouring cell values and differentiate.
        ! Get neighbor cells if they exist, otherwise use straight injection from local cell
        if (pop2_model == 2 .or. pop2_model == 3) then
          ncell = 1 ! we just want the neighbors of that cell
          ind_cell2(1) = ind_cell(i)
          call getnbor(ind_cell2,ind_nbor,ncell,ilevel)
          d1           = uold(ind_nbor(1,1),1) ; d2 = uold(ind_nbor(1,2),1) ; d3 = uold(ind_nbor(1,3),1)
          d4           = uold(ind_nbor(1,4),1) ; d5 = uold(ind_nbor(1,5),1) ; d6 = uold(ind_nbor(1,6),1)
          divv = 0d0
          !--------------------
          ! Divergence terms
          !--------------------
          ul        = (d2*uold(ind_nbor(1,2),2) + d*uold(ind_cell(i),2))/(d2+d)
          ur        = (d1*uold(ind_nbor(1,1),2) + d*uold(ind_cell(i),2))/(d1+d)
          divv      = divv + (ur-ul)
          ul        = (d4*uold(ind_nbor(1,4),3) + d*uold(ind_cell(i),3))/(d4+d)
          ur        = (d3*uold(ind_nbor(1,3),3) + d*uold(ind_cell(i),3))/(d3+d)
          divv      = divv + (ur-ul)
          ul        = (d6*uold(ind_nbor(1,6),4) + d*uold(ind_cell(i),4))/(d6+d)
          ur        = (d5*uold(ind_nbor(1,5),4) + d*uold(ind_cell(i),4))/(d5+d)
          divv      = divv + (ur-ul)
          ! Velocity vector divergence
          divv      = divv/dx_loc
          !-----------------------------------------------------------------------------------
          ! END  copy-and-pasted  from new star_formation.f90
          !-----------------------------------------------------------------------------------

          ! Divergence criterion (KS ADDED, cf. new public version of RAMSES as of 2020/03/25)
          if(divv>0d0)ok(i)=.false. !! gas should be converging
        end if

        !-----------------------------------------------------------------------------------
        ! check local potential minimum
        !-----------------------------------------------------------------------------------
        ! Note that getnbor gives the same index as the current cell if there is no next cell in the same level.
        ! The strict local min condition (phi < phimin_nbor) is never satisfied unless the current cell is surrounded
        ! by cells in the same level
        if (pop2_model == 3) then
          phi1           = phi(ind_nbor(1,1)) ; phi2 = phi(ind_nbor(1,2)) ; phi3 = phi(ind_nbor(1,3))
          phi4           = phi(ind_nbor(1,4)) ; phi5 = phi(ind_nbor(1,5)) ; phi6 = phi(ind_nbor(1,6))
          phimin_nbor = minval((/phi1,phi2,phi3,phi4,phi5,phi6/))
          if (phi(ind_cell(i)) >= phimin_nbor) ok(i)=.false. !! potential should be local minimum
        end if

        ! TO MODIFIED 2026-11-11: the peak condition (pop2_peak_only) and the virial condition (loc_vir > 0) are separate;
        ! with pop2_SFC_marker the cube is also evaluated for cells that failed only the div v / potential conditions
        if ((ok(i) .or. (pop2_SFC_marker .and. ok_pre_divv)) .and. (loc_vir > 0 .or. pop2_peak_only .or. pop2_SFC_marker)) then
          ! TO: check if the cell is self-gravitating
          ! first we need to compute the local velocity dispersion
          d = uold(ind_cell(i),1)
          cI2 = uold(ind_cell(i),5)*(gamma-1.0) !等温音速 (!上の方ですでにuold(ind_cell(i),1)で割り算されているので注意)
          ncell = 1 ! we just want the neighbors of the current cell
          ind_cell2(1) = ind_cell(i)
          call get_27_neighbors(ind_cell2, ind_nborcube, ilevel)
          vx = 0d0; vy = 0d0; vz = 0d0;
          vx2 = 0d0; vy2 = 0d0; vz2 = 0d0;
          m_in_cube = 0d0
          d = uold(ind_cell(i),1) ! density of the central cell
          dmax_nbor = d
          do j = 1, 27
            d1 = uold(ind_nborcube(j),1)
            dmax_nbor = max(dmax_nbor, d1)
            vx = vx + d1*uold(ind_nborcube(j),2)*vol_loc
            vy = vy + d1*uold(ind_nborcube(j),3)*vol_loc
            vz = vz + d1*uold(ind_nborcube(j),4)*vol_loc
            vx2 = vx2 + d1*vol_loc*(uold(ind_nborcube(j),2))**2
            vy2 = vy2 + d1*vol_loc*(uold(ind_nborcube(j),3))**2
            vz2 = vz2 + d1*vol_loc*(uold(ind_nborcube(j),4))**2
            m_in_cube = m_in_cube + d1*vol_loc
          end do

          !-----------------------------------------------------------------------------------
          ! TO MODIFIED 2026-11-11: this block (copied from the Pop III one) wrote its result into ok3 instead
          ! of ok, so the Pop II peak/virial condition never acted (ok3 is already .false. here since Z > Zpop3),
          ! and the formula was chosen by dpop3 (a Pop III parameter). Now: ok(i), and the Jeans-length form
          ! (the thermal Jeans threshold with cs^2 -> cs^2 + sigma^2) whenever pop2_njeans > 0.
          ! The peak condition is applied only with pop2_peak_only (a star particle without accretion in every
          ! cell above the threshold is the definition of eps_ff in the local models; the peak makes sense for the
          ! SFC centre).
          is_peak = (d >= dmax_nbor) ! the cell is the density maximum of its 27-cell cube
          vir_ok = .true.
          if (loc_vir > 0) then
            vx = vx / m_in_cube
            vy = vy / m_in_cube
            vz = vz / m_in_cube
            vx2 = vx2 / m_in_cube
            vy2 = vy2 / m_in_cube
            vz2 = vz2 / m_in_cube
            vdisp2 = (vx2 - vx**2) + (vy2 - vy**2) + (vz2 - vz**2)
            vdisp2 = max(vdisp2, 0d0)
            ! alpha_vir = 5 * (sigma^2 + c_s^2) * R / (G * M)
            ! here 4 pi / 3 * R^3 = 27 * dx^3  --> R = (3/4 pi * 27 * dV)^(1/3)
            if (pop2_njeans > 0) then
              lamJ_min = boxlen / (2.d0)**nlevelmax * pop2_njeans !ジーンズ長がこれ以上短い場合に星形成
              alpha_vir = twopi/2.*(cI2 + vdisp2)/(factG*d*(lamJ_min**2))
            else
              alpha_vir = 9.305257d0*(vdisp2 + cI2) * dx_loc / (factG * m_in_cube)
            endif
            ! write(*,*) 'sigma = ', sqrt(vdisp2)*scale_v*1.d-5, ' km/s cs = ', sqrt(cI2)*scale_v*1.d-5, ' km/s a_vir = ', alpha_vir
            vir_ok = (alpha_vir <= loc_vir)
          end if
          if (pop2_peak_only .and. .not. is_peak) ok(i) = .false.
          if (.not. vir_ok) ok(i) = .false.
          ! SFC marker candidate: the SFC-model condition (no div v) + density peak + virial condition
          ok_sfc(i) = pop2_SFC_marker .and. ok_pre_divv .and. is_peak .and. vir_ok
        endif
        !-------------- TO ADDED (END) -------------!


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! --- Pop II end
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! !------------------ KS DEBUG ----------------!
        ! d=uold(ind_cell(i),1)
        ! if (d>dpop2 .and. uold(ind_cell(i),imetal) > Zpop3*0.02) then
        !    if(ivar_refine>0)then
        !       write(*,'(2(A,L),A,I3, /,6(A, 1P2E10.3), (A, 3(1PE15.8)))') &
        !         'POPII conditions for d and Z are met. ok: ', ok(i), &
        !         ', leaf:',son(ind_cell(i))==0,&
        !         ', level:', ilevel,&
        !         'd, dcr [cm-3]:',d*scale_nH, dpop2*scale_nH,&
        !         ', Z, Zcr [Zsun]:',uold(ind_cell(i),imetal)/0.02, Zpop3,&
        !         ', Temp, Tcr [K]:',Temp, Tpop2,&
        !         ', divV, 0.:', divv,0d0,&
        !         ', phi, phimin_nbor:', phi(ind_cell(i)), phimin_nbor,&
        !         ', var_ref, var_cut_ref:', uold(ind_cell(i),ivar_refine), var_cut_refine,&
        !         ', r(:):', xg(ind_grid(i),1)+xc(ind,1),xg(ind_grid(i),2)+xc(ind,2),xg(ind_grid(i),3)+xc(ind,3)

        !    else
        !       write(*,'(2(A,L),A,I3, /,6(A, 1P2E10.3), (A, 3(1PE15.8)))') &
        !         'POPII conditions for d and Z are met. ok: ', ok(i), &
        !         ', leaf: ',son(ind_cell(i))==0,&
        !         ', level: ', ilevel,&
        !         'd, dcr [cm-3]:',d*scale_nH, dpop2*scale_nH,&
        !         ', Z, Zcr [Zsun]:',uold(ind_cell(i),imetal)/0.02, Zpop3,&
        !         ', Temp, Tcr [K]:',Temp, Tpop2,&
        !         ', divV, 0.:', divv,0d0,&
        !         ', phi, phimin_nbor:', phi(ind_cell(i)), phimin_nbor,&
        !         ', r(:): ', xg(ind_grid(i),1)+xc(ind,1),xg(ind_grid(i),2)+xc(ind,2),xg(ind_grid(i),3)+xc(ind,3)
        !    end if


        ! ! write(*,'(A,1P11E16.8)') &
        ! !      'vx,vy,vz,ul1,ur1,ul2,ur2,ul3,ur3,dx_loc,divv: ',&
        ! !      uold(ind_cell(i),2),uold(ind_cell(i),3),uold(ind_cell(i),4),&
        ! !      (d2*uold(ind_nbor(1,2),2) + d*uold(ind_cell(i),2))/(d2+d),(d1*uold(ind_nbor(1,1),2) + d*uold(ind_cell(i),2))/(d1+d),&
        ! !      (d4*uold(ind_nbor(1,4),3) + d*uold(ind_cell(i),3))/(d4+d),(d3*uold(ind_nbor(1,3),3) + d*uold(ind_cell(i),3))/(d3+d),&
        ! !      (d6*uold(ind_nbor(1,6),4) + d*uold(ind_cell(i),4))/(d6+d), (d5*uold(ind_nbor(1,5),4) + d*uold(ind_cell(i),4))/(d5+d),&
        ! !      dx_loc,divv
        ! ! write(*,'(A,1P7E15.8)') &
        ! !      'phi,phi1,phi2,phi3,phi4,phi5,phi6: ',&
        ! !      phi(ind_cell(i)),phi1,phi2,phi3,phi4,phi5,phi6

        ! end if

        ! !------------------ KS DEBUG ----------------!

      enddo

      ! Calculate number of new stars in each cell using Poisson statistics - (Poisson statistics is for Pop II only)
      !--------------------------------------------
      ! POP II
      !--------------------------------------------
      ! - pop2_model = 4 or 5
      !    collect info for SFC candidates
      ! - other models
      !    calculate nstar from ok and increment ntot
      !--------------------------------------------
      ! POP III
      !--------------------------------------------
      ! - pop3_model = 4, 5, or 6
      !    collect info for SFC candidates
      ! - other models
      !    calculate nstar from ok and increment ntot

      do i=1,ngrid
        nstar(i)=0 ; nstar3(i)=0

        !! Pop III - JP
        if(ok3(i))then
          ! Gas mass
          d=uold(ind_cell(i),1)
          mcell=d*vol_loc
          ! A single Pop III star forms.
          nstar3(i)=1

          !create star irrespective of the cell mass if nonlocal_depletion = .true. (KS MODIFIED)
          !or pop3_model = 4, 5, or 6
          !if (.not. nonlocal_depletion) then
          if (.not. (nonlocal_depletion .or. pop3_model == 4 .or. pop3_model == 5 .or. pop3_model == 6)) then

            ! Check if there are remaining gases
            if(mcell.lt.1.2*mpop3imf(i))then ! Guarantee that there are enough gas remaining.
              ! !---------------- KS DEBUG ------------------!
              ! print '(A, 6(1PE12.4))', "KS DEBUG (star_formation.f90) rejected due to cell mass: x,y,z,nH,mcell,mpop3 =",&
              !      xg(ind_grid(i),1)+xc(ind,1),&
              !      xg(ind_grid(i),2)+xc(ind,2),&
              !      xg(ind_grid(i),3)+xc(ind,3),&
              !      d*scale_nH, mcell*scale_msun, mpop3imf(i)*scale_msun
              ! !---------------- KS DEBUG ------------------!
              ok3(i)=.false. ! Star formation does not occur.
              nstar3(i)=0
              mpop3imf(i)=0
            endif

          end if !nonlocal_depletion

          ! Compute new stars local statistics
          mstar_tot=mstar_tot+nstar3(i)*mpop3imf(i)
          if(nstar3(i)>0)then
            ntot=ntot+1
            if(f_w>0)ndebris_tot=ndebris_tot+1
          endif



          ! count/store pop3 candidates (KS ADDED)
          if (nstar3(i)>0) then
            p3cand_ncand_loc = p3cand_ncand_loc + 1 !number of candidates
            ! check the number of candidates
            if (p3cand_ncand_loc > Npop3_max) then
              write(*,*)'Number of local Pop III candidates exceeds the limit ',p3cand_ncand_loc, Npop3_max
              write(*,*)'stopping...'
              stop
            end if
            !store the ind_cell/positions
            iloc = p3cand_ncand_loc
            p3cand_indcell_loc(iloc) = ind_cell(i)
            p3cand_pos_loc(iloc,:) = (/(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale,&
              (xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale, &
              (xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale/)
          end if

          if (nstar3(i)>0) then
            yHI = uold(ind_cell(i),iIons-1+ixHI) !上の方ですでにuold(ind_cell(i),1)で割り算されているので注意
            yHII = uold(ind_cell(i),iIons-1+ixHII) !上の方ですでにuold(ind_cell(i),1)で割り算されているので注意
            yH2 = (1.-yHI-yHII)*0.5
            write(*,'(A, 6(A, 1PE9.2), /, (A, 1PE9.2), (A, 3(1PE15.8)))') &
              '** POP III candidate ** ', &
              'nH [cm-3]:',uold(ind_cell(i),1)*scale_nH,&
              ', T2 [K]:', uold(ind_cell(i),5)*scale_T2*(gamma-1.0), &
              ', yHI:', yHI,&
              ', yHII:', yHII,&
              ', yH2:', yH2,&
              ', Z [Zsun]:',uold(ind_cell(i),imetal)/0.02,&
              'mcell [Msun]:',uold(ind_cell(i),1)*vol_loc*scale_msun,&
              ', r(:):', (xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale,&
              (xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale,&
              (xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
          end if
        endif
        !! --- Pop III end

        if(ok(i))then
          !if (pop2_model /= 4) then !except for the SFC model
          if (pop2_model == 1 .or. pop2_model == 2 .or. pop2_model == 3) then !except for the SFC model

            ! Compute mean number of events
            d=uold(ind_cell(i),1)
            mcell=d*vol_loc

            !**** WARNING WARNING WARNING (KS MODIFIED) WARNING WARNING WARNING ******!
            !tstar=t0*sqrt(d00/d)
            tstar= tff0_sec/scale_t / sqrt(d*scale_nH) / eps_star  !tff0_secはnH=1cm-3のfree-fall time, eps_starは星形成率
            !**** WARNING WARNING WARNING (KS MODIFIED) WARNING WARNING WARNING ******!


            ! Poisson mean of the number of new stars
            PoissMean=dtnew(ilevel)/tstar*mcell/mstar

            !H2-abundance dependent SF rate
            if (pop2_H2depSF) then
              yHI = uold(ind_cell(i),iIons-1+ixHI) !上の方ですでにuold(ind_cell(i),1)で割り算されているので注意
              yHII = uold(ind_cell(i),iIons-1+ixHII) !上の方ですでにuold(ind_cell(i),1)で割り算されているので注意
              yH2 = (1.-yHI-yHII)*0.5
              PoissMean = PoissMean * Xhydrogen * yH2 !SF rate is now proportional to rho(H2) instead of rho
            end if


            ! Compute Poisson realisation
            call poissdev(localseed,PoissMean,nstar(i))
            ! Compute depleted gas mass
            mgas=nstar(i)*mstar

            !create star irrespective of the cell mass if nonlocal_depletion = .true. (KS MODIFIED)
            if (.not. nonlocal_depletion) then

              ! Security to prevent more than 50% of gas depletion
              if (mgas > 0.5*mcell) then
                ! !---------------- KS DEBUG ------------------!
                ! print '(A, 6(1PE12.4),2I4)', "KS DEBUG (star_formation.f90) Pop II corrected/rejected due to cell mass: x,y,z,nH,mcell,mstar,nstar(before),nstar(after) =",&
                !      xg(ind_grid(i),1)+xc(ind,1),&
                !      xg(ind_grid(i),2)+xc(ind,2),&
                !      xg(ind_grid(i),3)+xc(ind,3),&
                !      d*scale_nH, mcell*scale_msun, mstar*scale_msun, nstar(i), int(0.5*mcell/mstar)
                ! !---------------- KS DEBUG ------------------!

                nstar_corrected=int(0.5*mcell/mstar)
                mstar_lost=mstar_lost+(nstar(i)-nstar_corrected)*mstar
                nstar(i)=nstar_corrected
              endif

            end if !nonlocal_depletion

            ! Compute new stars local statistics
            mstar_tot=mstar_tot+nstar(i)*mstar
            if(nstar(i)>0)then
              ntot=ntot+1
              if(f_w>0)ndebris_tot=ndebris_tot+1
            endif

            !------------------ KS ADDED ----------------!
            if (nstar(i)>0) then
              yHI = uold(ind_cell(i),iIons-1+ixHI)   !上の方ですでにuold(ind_cell(i),1)で割り算されているので注意
              yHII = uold(ind_cell(i),iIons-1+ixHII) !上の方ですでにuold(ind_cell(i),1)で割り算されているので注意
              yH2 = (1.-yHI-yHII)*0.5
              write(*,'(A, I3, 4(A, 1PE9.2), /, 6(A, 1PE9.2), (A, 3(1PE15.8)))') &
                '** POP II formation ** nstar:' , nstar(i), &
                ', PoissMean:', PoissMean, &
                ', dt: [Myr]:',dtnew(ilevel)*scale_t/(1d6*365.*24.*3600.),&
                ', tstar: [Myr]:',tstar*scale_t/(1d6*365.*24.*3600.),&
                ', mcell [Msun]:',mcell*scale_msun,&
                'nH [cm-3]:',uold(ind_cell(i),1)*scale_nH,&
              ! ', T [K]:', uold(ind_cell(i),5)*scale_T2*(gamma-1.0),&
                ', T2 [K]:', uold(ind_cell(i),5)*scale_T2*(gamma-1.0),&
                ', yHI:', yHI,&
                ', yHII:', yHII,&
                ', yH2:', yH2,&
                ', Z [Zsun]:',uold(ind_cell(i),imetal)/0.02,&
                ', r(:):', (xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale,&
                (xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale,&
                (xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale

            end if

            ! ! ------------------ KS DEBUG ------------------- !
            ! if (nstar(i)==0) then
            !    write(*,'(A, I3, 4(A, 1PE9.2), /, 6(A, 1PE9.2), (A, 3(1PE15.8)))') &
            !         'POP II unsuccessful candidate. nstar:' , nstar(i), &
            !         ', PoissMean:', PoissMean, &
            !         ', dt: [Myr]:',dtnew(ilevel)*scale_t/(1d6*365.*24.*3600.),&
            !         ', tstar: [Myr]:',tstar*scale_t/(1d6*365.*24.*3600.),&
            !         ', mcell [Msun]:',mcell*scale_msun,&
            !         'nH [cm-3]:',uold(ind_cell(i),1)*scale_nH,&
            !         ', T2 [K]:', uold(ind_cell(i),5)*scale_T2*(gamma-1.0),&
            !         ', yHI:', yHI,&
            !         ', yHII:', yHII,&
            !         ', yH2:', yH2,&
            !         ', Z [Zsun]:',uold(ind_cell(i),imetal)/0.02,&
            !         ', r(:):', xg(ind_grid(i),1)+xc(ind,1),xg(ind_grid(i),2)+xc(ind,2),xg(ind_grid(i),3)+xc(ind,3)
            ! end if
            ! ! ------------------ KS DEBUG ------------------- !
            !------------------ KS ADDED ----------------!


          else !for the SFC model

            !no more SFC formation if test_problem == 'BE_FB_1cluster' .and. nstar_tot > 0

            ! count/store SFC candidates (KS ADDED)
            SFC_ncand_loc = SFC_ncand_loc + 1 !number of candidates
            ! check the number of candidates
            if (SFC_ncand_loc > SFC_n_max) then
              write(*,*)'Number of local SFC candidates exceeds the limit ',SFC_ncand_loc, SFC_n_max
              write(*,*)'stopping...'
              stop
            end if

            !store the ind_cell/positions
            iloc = SFC_ncand_loc
            SFC_indcellcand_loc(iloc) = ind_cell(i)
            SFC_poscand_loc(:,iloc) = (/(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale,&
              (xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale, &
              (xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale/)

            yHI = uold(ind_cell(i),iIons-1+ixHI)   !上の方ですでにuold(ind_cell(i),1)で割り算されているので注意
            yHII = uold(ind_cell(i),iIons-1+ixHII) !上の方ですでにuold(ind_cell(i),1)で割り算されているので注意
            yH2 = (1.-yHI-yHII)*0.5
            write(*,'(A, /, 6(A, 1PE9.2), (A, 3(1PE15.8)))') &
              '** POP II SFC candidate detected **' , &
              'nH [cm-3]:',uold(ind_cell(i),1)*scale_nH,&
              ', T2 [K]:', uold(ind_cell(i),5)*scale_T2*(gamma-1.0),&
              ', yHI:', yHI,&
              ', yHII:', yHII,&
              ', yH2:', yH2,&
              ', Z [Zsun]:',uold(ind_cell(i),imetal)/0.02,&
              ', r(:):', (xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale,&
              (xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale,&
              (xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
            flush(6)


            ! print *,"KS DEBUG (after POP II SFC candidate detection)"
            ! print *,"myid =" ,myid
            ! print *,"ind_cell(i)=",ind_cell(i)

          end if

        endif
        !! --- Pop II end

        ! TO ADDED 2026-11-11: SFC marker candidate (local Pop II models with pop2_SFC_marker): the cell goes on the
        ! SFC candidate list as in the SFC model, independently of the stochastic star formation above
        if (ok_sfc(i)) then
          call store_SFC_candidate(ind_cell(i), (/(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale, &
            (xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale, (xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale/))
        end if
      enddo !i=1,ngrid

      do i=1,ngrid
        ! Store nstar in array flag2
        flag2(ind_cell(i))=nstar(i)
        ! Store pop3 star mass
        if(ilevel.eq.nlevelmax)flag3(ind_cell(i))=mpop3imf(i)
      end do
    end do !ind=1,twotondim
  end do !igrid=1,ncache,nvector

  !**** WARNING WARNING WARNING (KS MODIFIED) WARNING WARNING WARNING ******!
  ! ! JP
  ! if(ilevel.eq.nlevelmax)deallocate(halo_overlap,halo_SF,halo_pop3)
  !**** WARNING WARNING WARNING (KS MODIFIED) WARNING WARNING WARNING ******!

  !--------------------------------------------------------------------
  ! remove pop3 candidates simultaneously forming in close vicinities (KS ADDED)
  !--------------------------------------------------------------------
  call remove_p3cand_by_other_cand()

  !--------------------------------------------------------------------
  ! add Pop III star particles (mdot-dep model)
  !--------------------------------------------------------------------
  if (pop3_model == 4 .or. pop3_model == 5 .or. pop3_model == 6)   then
    call mdot_dependent_pop3_formation
  end if


  !--------------------------------------------------------------------
  ! add Pop II star particles (SFC-based model)
  !--------------------------------------------------------------------
  if (pop2_model == 4 .or. pop2_model == 5)   then
    call SFC_from_candidate !generate star-forming clouds from their candidats

    if (SFC_n_tot > 0)  then
      !KS NOTE: ufront is set to variable corresponding to unew after reflecting the advection effect
      !update ufront if not updated
      if (.not. ufront_updated) then
        do sweep_level=levelmin,nlevelmax
          !call set_ufront(sweep_level)
          call set_ufront(sweep_level,flag_option=3) !use unew to set ufront
        end do
        ufront_updated = .true.
      end if

      ! TO: We will do this later
      !do iSFC = 1, SFC_n_tot call output_around_SFC
      ! call output_around_SFC(iSFC, SFC_pos_glob(1:3,iSFC),use_current_ufront=.true.)      !output info about SFC
      !end do



      ! !KS DEBUG TENTATIVE
      ! do sweep_level = levelmin,nlevelmax
      !    do ind=1,twotondim
      !       iskip=ncoarse+(ind-1)*ngridmax
      !       do i=1,active(sweep_level)%ngrid
      !          if (myid == 5 .and. active(sweep_level)%igrid(i)+iskip == 4821) then
      !             print '(A,2I4,I8,1P3E12.4)', &
      !                  '(KS DEBUG TENTATIVE after output_around_SFC) myid, level, ind_cell, nH(uold), nH(unew), nH(ufront): ',&
      !                  myid, sweep_level, active(sweep_level)%igrid(i)+iskip, uold(active(sweep_level)%igrid(i)+iskip,1)*scale_nH, unew(active(sweep_level)%igrid(i)+iskip,1)*scale_nH,&
      !                  ufront(active(sweep_level)%igrid(i)+iskip,1)*scale_nH
      !          end if
      !       end do
      !    end do
      ! end do


      !------- KS DEBUG --------!
      if (myid == 1) then
        print '(A)', 'call check_cell_health BEFORE star_formation_in_SFC'
      end if
      call check_cell_health !check if all cells are healthy
      !-------------------------!
      call star_formation_in_SFC !form stars in star-forming clouds found above
      !------- KS DEBUG --------!
      if (myid == 1) then
        print '(A)', 'call check_cell_health AFTER star_formation_in_SFC'
      end if
      call check_cell_health !check if all cells are healthy
      !-------------------------!
    end if
  end if

  !--------------------------------------------------------------------
  ! TO ADDED 2026-11-11: SFC marker particles in the local Pop II models (pop2_SFC_marker): the cells that satisfy
  ! the SFC-model condition (candidates collected above; one per SFC_dist_min_pc as in SFC_from_candidate, and none
  ! within SFC_dist_min_pc of an active SFC marker) get a massless FAM_SFC particle; its formation nest
  ! (output_around_SFC) and age outputs (SFC_output_at_given_age) are those of the SFC model. No star formation
  ! and no gas cropping here: the stars come from the local model.
  !--------------------------------------------------------------------
  if (pop2_SFC_marker) then
    call SFC_from_candidate
    if (SFC_n_tot > 0) then
      if (.not. ufront_updated) then
        do sweep_level=levelmin,nlevelmax
          call set_ufront(sweep_level,flag_option=3) !use unew to set ufront
        end do
        ufront_updated = .true.
      end if
      call create_SFC_markers
    end if
  end if



  !---------------------------------
  ! Check for free particle memory
  !---------------------------------
  ok_free=(numbp_free-ntot-ndebris_tot)>=0
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(numbp_free,numbp_free_tot,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  numbp_free_tot=numbp_free
#endif
  if(.not. ok_free)then
    write(*,*)'No more free memory for particles'
    write(*,*)'Increase npartmax'
#ifndef WITHOUTMPI
    call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
    stop
#endif
  end if

  !---------------------------------
  ! Compute global stars statistics
  !---------------------------------
#ifndef WITHOUTMPI
  mlost=mstar_lost; mtot=mstar_tot
  call MPI_ALLREDUCE(ntot,ntot_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mtot,mtot_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mlost,mlost_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  ntot_all=ntot
  mtot_all=mstar_tot
  mlost_all=mstar_lost
#endif
  ntot_star_cpu=0; ntot_star_all=0
  ntot_star_cpu(myid)=ntot
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ntot_star_cpu,ntot_star_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  ntot_star_cpu(1)=ntot_star_all(1)
#endif
  do icpu=2,ncpu
    ntot_star_cpu(icpu)=ntot_star_cpu(icpu-1)+ntot_star_all(icpu)
  end do
  nstar_tot=nstar_tot+ntot_all



  if(myid==1)then
    !------------------------- KS DEBUG ----------------------------!
    ! write(*,'("(KS DEBUG) Level=",I6,", New star=",I6,", Tot=",I10,", Mass=",1PE10.3)')&
    !      & ilevel,ntot_all,nstar_tot,mtot_all*scale_msun
    ! mtot_allは形成時のmassの合計 (SN等によるmass lossの効果は含まない)
    !if (ilevel == levelmin) then
    ! write(*,'("Star formation by this step:  Ntot=",I10,", Mtot=",1PE10.3)')&
    !      nstar_tot,mtot_all*scale_msun
    !end if
    !------------------------- KS DEBUG ----------------------------!
    ! if(ntot_all.gt.0)then
    !      write(*,'(" Level=",I6," New star=",I6," Tot=",I10," Mass=",1PE10.3," Lost=",0PF4.1,"%")')&
    !           & ilevel,ntot_all,nstar_tot,mtot_all,mlost_all/(mlost_all+mtot_all)*100.
    !   endif
  end if

  !------------------------------------------------------------------
  ! Create new star particles
  !------------------------------------------------------------------
  !  POP III in case of pop3_model = 1-3
  !   otherwise Pop III particles are created in mdot_dependent_pop3_formation
  !------------------------------------------------------------------
  !  POP II in case of pop2_model = 1-3
  !   otherwise Pop II particles are created in star_formation_in_SFC
  !------------------------------------------------------------------
  ! Starting identity number
  if(myid==1)then
    index_star=nstar_tot-ntot_all
  else
    index_star=nstar_tot-ntot_all+ntot_star_cpu(myid-1) ! (KS MEMO) 各cpu毎にはじまりのindex_starが違うのでそれに対応
  end if

  ! Loop over grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
    ngrid=MIN(nvector,ncache-igrid+1)
    do i=1,ngrid
      ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
    end do

    ! Loop over cells
    do ind=1,twotondim
      iskip=ncoarse+(ind-1)*ngridmax
      do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
      end do

      ! Gather new star arrays
      nnew=0 ; ok3_new=.false.
      do i=1,ngrid
        if (flag2(ind_cell(i)).gt.0)then
          nnew=nnew+1
          ind_grid_new(nnew)=ind_grid(i)
          ind_cell_new(nnew)=ind_cell(i)
        else if (flag3(ind_cell(i)).gt.0)then !! tell Pop III and Pop II - JP
          nnew=nnew+1
          ind_grid_new(nnew)=ind_grid(i)
          ind_cell_new(nnew)=ind_cell(i)
          ok3_new(nnew)=.true.
        endif
      end do

      ! Update linked list for stars
      call remove_free(ind_part,nnew)
      call add_list(ind_part,ind_grid_new,ok_new,nnew)

      ! Update linked list for debris
      if(f_w>0)then
        call remove_free(ind_debris,nnew)
        call add_list(ind_debris,ind_grid_new,ok_new,nnew)
      endif

      ! Calculate new star particle and modify gas density
      do i=1,nnew
        index_star=index_star+1

        ! Get gas variables
        n=flag2(ind_cell_new(i))
        d=uold(ind_cell_new(i),1)
        u=uold(ind_cell_new(i),2)
        v=uold(ind_cell_new(i),3)
        w=uold(ind_cell_new(i),4)
        x=(xg(ind_grid_new(i),1)+xc(ind,1)-skip_loc(1))*scale
        y=(xg(ind_grid_new(i),2)+xc(ind,2)-skip_loc(2))*scale
        z=(xg(ind_grid_new(i),3)+xc(ind,3)-skip_loc(3))*scale
        if(metal)zg=uold(ind_cell_new(i),imetal)

        ! Set star particle variables
        tp(ind_part(i))=birth_epoch  ! Birth epoch
        mp(ind_part(i))=n*mstar      ! Mass
        levelp(ind_part(i))=ilevel   ! Level
        idp(ind_part(i))=index_star  ! Star identity (negative integers will be allocated to sinks)
        typep(ind_part(i))%family=FAM_STAR  ! Particle type
        if(ok3_new(i))then
          mp(ind_part(i))=flag3(ind_cell_new(i))  ! Pop III
          typep(ind_part(i))%family=FAM_POP3      ! Particle type
        endif
        xp(ind_part(i),1)=x
        xp(ind_part(i),2)=y
        xp(ind_part(i),3)=z
        vp(ind_part(i),1)=u
        vp(ind_part(i),2)=v
        vp(ind_part(i),3)=w
        if(metal)zp(ind_part(i))=zg  ! Initial star metallicity
#ifdef INIT_STELLAR_MASS
        mp0(ind_part(i))=mp(ind_part(i))
#endif
#ifdef STELLAR_POPULATION_MASS
        msp0(ind_part(i))=mp(ind_part(i))
#endif
        if(sf_log_properties) call sf_log_write_star(ilun_sf,ind_part(i),ind_cell_new(i),ilevel,.true.)

        !------------------------- KS DEBUG ----------------------------!
        ! write(*,'(A,1PE9.2,A,I9,A,I3,A,1PE9.2)') &
        !      "** New star added **  mp [Msun]:", mp(ind_part(i))*scale_msun,&
        !      ", idp:",         idp(ind_part(i)),&
        !      ", typep%family:", typep(ind_part(i))%family,&
        !      ", zp [Zsun]:", zp(ind_part(i))/0.02
        write(*,'(A,I6,A,1PE9.2,A,I9,A,I3,A,1PE9.2,A,1P3E20.12)') &
          "** New star added ** myid: ", myid, &
          " mp [Msun]:", mp(ind_part(i))*scale_msun,&
          ", idp:",         idp(ind_part(i)),&
          ", typep%family:", typep(ind_part(i))%family,&
          ", zp [Zsun]:", zp(ind_part(i))/0.02, &
          ", pos:",  xp(ind_part(i),1:3)
        !------------------------- KS DEBUG ----------------------------!


        ! Set GMC particle variables
        if(f_w>0)then
          ! Compute GMC mass without more than 50% of gas depletion
          mdebris=min(f_w*n*mstar,0.5*d*vol_loc-n*mstar)
          ! Add supernova ejecta
          mdebris=mdebris+eta_sn*n*mstar
          ! Remove ejecta from the long lived star mass
          mp(ind_part(i))=n*mstar-eta_sn*n*mstar
          ! Set GMC particle variables
          tp(ind_debris(i))=birth_epoch  ! Birth epoch
          mp(ind_debris(i))=mdebris      ! Mass
          levelp(ind_debris(i))=ilevel   ! Level
          idp(ind_debris(i))=n           ! Number of individual stars
          xp(ind_debris(i),1)=x
          xp(ind_debris(i),2)=y
          xp(ind_debris(i),3)=z
          vp(ind_debris(i),1)=u
          vp(ind_debris(i),2)=v
          vp(ind_debris(i),3)=w
          ! GMC metallicity + yield from ejecta
          if(metal)zp(ind_debris(i))=zg+eta_sn*yield*(1-zg)*n*mstar/mdebris
        endif

      end do
      ! End loop over new star particles



      !depletion will be reflected later if nonlocal_depletion = .true. (KS MODIFIED)
      !depletion is skipped if is_test = .true.
      !if (.not. is_test) then    !KS ADDED
      if (.not. nonlocal_depletion) then

        ! Modify gas density according to mass depletion
        do i=1,ngrid
          n=flag2(ind_cell(i))
          d=uold(ind_cell(i),1)

          if(flag2(ind_cell(i))>0)uold(ind_cell(i),1)=max(d-n*dstar*(1.0+f_w),0.5*d)
          if(flag3(ind_cell(i))>0)uold(ind_cell(i),1)=d-flag3(ind_cell(i))*(1.0+f_w)/vol_loc
          !! At the moment, use f_w=0 - JP
          !! Divde by vol_loc - density !! - JP
          !the above lines implicitly assumes Pop II and III do not simultaneaously form at the same cell (KS NOTE)

          !------------------------- KS DEBUG ----------------------------!
          if(flag2(ind_cell(i))>0 .or. flag3(ind_cell(i))>0) then
            write(*,'(4(A,1PE9.2))') &
              "** gas depletion due to star formation **  Mstar [Msun]:", n*mstar*scale_msun + flag3(ind_cell(i))*scale_msun, &
              "  nHstar [cm-3]: ",   (n*dstar + flag3(ind_cell(i))/vol_loc)*scale_nH, &
              ", nH (before) [cm-3]: ",        d*scale_nH,&
              ", nH (after) [cm-3]: ",         uold(ind_cell(i),1)*scale_nH
          end if
          !------------------------- KS DEBUG ----------------------------!
        end do

      end if !nonlocal_depletion/is_test

    end do
    ! End loop over cells
  end do
  ! End loop over grids

  !**** WARNING WARNING WARNING (KS MODIFIED) WARNING WARNING WARNING ******!
  ! if(ilevel.eq.nlevelmax)deallocate(xp3)
  !**** WARNING WARNING WARNING (KS MODIFIED) WARNING WARNING WARNING ******!

  !---------------------------------------------------------
  ! depletion of gas from spherical regions
  !---------------------------------------------------------
  if (nonlocal_depletion) then

    !if there are new stars (count only Pop III in case of pop3_model = 1-3 and Pop II in case of pop2_model = 1-3)
    if (ntot_all > 0) then
      !update ufront if not updated
      if (.not. ufront_updated) then
        do sweep_level=levelmin,nlevelmax
          !call set_ufront(sweep_level)
          call set_ufront(sweep_level,flag_option=3) !use unew to set ufront
        end do
        ufront_updated = .true.
      end if

      !choose either model of gas depleteion
      if (nH_postSF <= 0) then
        if (nonlocal_depletion_batch) then
          call deplete_gas_from_spheres_batch !all new stars at once (merge_ks)
        else
          call deplete_gas_from_spheres !deplete gas from fintite spheres (ks: one star after the other)
        end if
      else
        call deplete_gas_with_fixed_postSF_density !deplete gas from non-local regions with fixed density of post-star-forimng cells
      end if
    end if
  end if


  !---------------------------------------------------------
  ! Convert hydro variables back to conservative variables
  !---------------------------------------------------------
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
    ngrid=MIN(nvector,ncache-igrid+1)
    do i=1,ngrid
      ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
    end do
    do ind=1,twotondim
      iskip=ncoarse+(ind-1)*ngridmax
      do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
      end do
      do i=1,ngrid
        d=uold(ind_cell(i),1)
        u=uold(ind_cell(i),2)
        v=uold(ind_cell(i),3)
        w=uold(ind_cell(i),4)
        e=uold(ind_cell(i),5)*d
#ifdef SOLVERmhd
        bx1=uold(ind_cell(i),6)
        by1=uold(ind_cell(i),7)
        bz1=uold(ind_cell(i),8)
        bx2=uold(ind_cell(i),nvar+1)
        by2=uold(ind_cell(i),nvar+2)
        bz2=uold(ind_cell(i),nvar+3)
        e=e+0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)
#endif
        e=e+0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
        do irad=1,nener
          e=e+uold(ind_cell(i),5+irad)
        end do
#endif
        uold(ind_cell(i),1)=d
        uold(ind_cell(i),2)=d*u
        uold(ind_cell(i),3)=d*v
        uold(ind_cell(i),4)=d*w
        uold(ind_cell(i),5)=e
      end do
      do ivar=imetal,nvar
        do i=1,ngrid
          d=uold(ind_cell(i),1)
          w=uold(ind_cell(i),ivar)
          uold(ind_cell(i),ivar)=d*w
        end do
      end do
    end do
  end do


  !*** OBSOLETE *** Does Not Work *** OBSOLETE ***
  ! !----------------------------------------------------------
  ! ! update cells at ilevel < nlevelmax (KS ADDED)
  ! !----------------------------------------------------------
  ! ! consistently reflect the update at ilevel < nlevelmax,
  ! ! if star formation occur with non-local gas depletion
  ! !
  ! ! Note: star_formation is called at the end of amr_step at ilevel = nlevelmax,
  ! ! corresoponding to just after calling set_unew(ilevel) at ilevel < nlevelmax
  ! !----------------------------------------------------------
  ! if ((nonlocal_depletion .or. pop2_model == 4) .and. &      !if non-local star fomration occurs
  !      (sum(ntot_star_all) > 0 .or. sum(SFC_nstar_list) > 0)) then
  !    !ntot_star_all: # of normal stars, SFC_nstar_list # of SFC stars
  !    if (myid == 1) then
  !       print '(A, 2(A,I8))', '(star_formation) reflect the update at lower-levels after non-local depletion: ',&
  !            '# of stars = ', sum(ntot_star_all), ', # of SFC stars = ', sum(SFC_nstar_list)
  !    end if
  !    do sweep_level=nlevelmax-1,levelmin,-1
  !       call upload_fine(sweep_level+1)    !sweep_level+1 -> sweep_level
  !       do ivar=1,nvar
  !          call make_virtual_fine_dp(uold(1,ivar),sweep_level) !reset boundary values
  !       enddo
  !       call set_unew(sweep_level)        !reset unew to updated uold
  !    enddo
  ! end if
  ! !----------------------------------------------------------
  ! !----------------------------------------------------------

  !if ((nonlocal_depletion .or. pop2_model == 4) .and. &      !if non-local star fomration occurs
  if ((nonlocal_depletion .or. pop2_model == 4 .or. pop2_model == 5) .and. &  !if non-local star fomration occurs
    (sum(ntot_star_all) > 0 .or. sum(SFC_nstar_list) > 0)) then
    !ntot_star_all: # of normal stars, SFC_nstar_list # of SFC stars
    if (myid == 1) then
      print '(A, 2(A,I8))', '(star_formation) do virtual boundary communication after non-local depletion: ',&
        '# of stars = ', sum(ntot_star_all), ', # of SFC stars = ', sum(SFC_nstar_list)
    end if
    do sweep_level=levelmin, nlevelmax
#ifdef SOLVERmhd
      call make_virtual_fine_dp_nv(uold,nvar+3,1,nvar+3,sweep_level)
#else
      call make_virtual_fine_dp_nv(uold,nvar,1,nvar,sweep_level)
#endif
    enddo

    !KS MEMO: maybe upload_fine is needed but I have no idea

  end if

#endif !NDIM==3

contains
  !------------------------------------------
  ! output summary (KS ADDED)
  !------------------------------------------
  subroutine star_formation_summary()
    ! number and mass of the particles of each FG family (1 Pop II, 2 Pop III, 3 SN, 4 dead, 5 BH, 6 SFC, 7 PSC)
    ! and the accumulated stellar mass; two ALLREDUCEs (ks: 15)
    integer,parameter :: nfam = 7
    integer,dimension(1:nfam) :: nloc, nglob
    real(dp),dimension(1:nfam+1) :: mloc, mglob
    integer :: ifam
    nloc(:) = 0; mloc(:) = 0d0
    do i=1,npartmax ! find particles in each cpu
      if (levelp(i) == 0) cycle   ! free slot
      select case (typep(i)%family)
      case (FAM_STAR); ifam = 1
      case (FAM_POP3); ifam = 2
      case (FAM_SN);   ifam = 3
      case (FAM_DEAD); ifam = 4
      case (FAM_BH);   ifam = 5
      case (FAM_SFC);  ifam = 6
      case (FAM_PSC);  ifam = 7
      case default;    cycle
      end select
      nloc(ifam) = nloc(ifam) + 1
      mloc(ifam) = mloc(ifam) + mp(i)
    enddo
    mloc(nfam+1) = mstar_tot
    call MPI_ALLREDUCE(nloc,nglob,nfam,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(mloc,mglob,nfam+1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)

    !---- output message -----!
    if(myid==1 .and. nstar_tot>0) then
      write(*,'("Star formation by this step:  Ntot=",I10,", Mtot=",1PE12.3)')&
        nstar_tot,mglob(nfam+1)*scale_msun
      write(*,'(7(A, I8))') &
        'Current summary for particle number: Pop II=', nglob(1),&
        ', Pop III=', nglob(2), &
        ', SN=', nglob(3), &
        ', DEAD(PopII)=', nglob(4), &
        ', BH(PopIII)=', nglob(5), &
        ', SFC=', nglob(6), &
        ', PSC=', nglob(7)
      write(*,'(7(A, 1PE12.3))') &
        'Current summary for particle mass[Msun]: Pop II=', mglob(1)*scale_msun, &
        ', Pop III=', mglob(2)*scale_msun, &
        ', SN=', mglob(3)*scale_msun, &
        ', DEAD(PopII)=', mglob(4)*scale_msun, &
        ', BH(PopIII)=', mglob(5)*scale_msun, &
        ', SFC=', mglob(6)*scale_msun, &
        ', PSC=', mglob(7)*scale_msun
    end if
  end subroutine star_formation_summary


  !-----------------------------------------------------------------------
  !get list for pre-existing pop3 particles in all cpus (KS ADDED)
  !-----------------------------------------------------------------------
  subroutine get_pop3_pos_list

    !buffer for MPI send/recv (the size should be larger than both of ncpu and Npop3_max)
    integer,dimension(1:max(ncpu,Npop3_max*3)) :: sbufi,rbufi
    real(dp),dimension(1:max(ncpu,Npop3_max*3)) :: sbufr,rbufr

    real(dp) :: age

    !get local npop3 and position
    do i=1,npartmax
      if( is_pop3(typep(i)) ) then
        !consider only Pop III stars with age < pop3_tMyr_suppress_newstar
        call getAgeGyr(tp(i), age)
        if (age*1d3 > pop3_tMyr_suppress_newstar) then
          cycle !skip
        end if

        p3_npop3_loc=p3_npop3_loc+1
        ! check the number of candidates
        if (p3_npop3_loc > Npop3_max) then
          write(*,*)'Number of local Pop III partilces exceeds the limit ',p3_npop3_loc, Npop3_max
          write(*,*)'stopping...'
          stop
        end if
        !store the positions
        iloc = p3_npop3_loc
        p3_pos_loc(iloc,:) = xp(i,:)
      end if
    enddo

    !create p3_npop3_list
    p3_npop3_list(:) = 0
    p3_npop3_list(myid) = p3_npop3_loc
    sbufi(1:ncpu) = p3_npop3_list(:)
    call MPI_ALLREDUCE(sbufi,rbufi,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
    p3_npop3_list(:) = rbufi(1:ncpu)
    p3_npop3_tot = sum(p3_npop3_list) !total number of candidates

    p3_offset = 0
    do icpu = 1, myid-1  ! "do i = 1,0" loop is just skipped (KS MEMO)
      p3_offset = p3_offset + p3_npop3_list(icpu)
    end do

    !create global pop3 position list
    p3_pos_glob(:,:)=0 !initialization
    do iloc = 1, p3_npop3_loc
      iglob = iloc + p3_offset
      p3_pos_glob(iglob,:) = p3_pos_loc(iloc,:)
    end do
    do iglob = 1, p3_npop3_tot
      sbufr((iglob-1)*3+1) = p3_pos_glob(iglob,1)
      sbufr((iglob-1)*3+2) = p3_pos_glob(iglob,2)
      sbufr((iglob-1)*3+3) = p3_pos_glob(iglob,3)
    end do
    call MPI_ALLREDUCE(sbufr, rbufr, 3*p3_npop3_tot, MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    do iglob = 1, p3_npop3_tot
      p3_pos_glob(iglob,1) = rbufr((iglob-1)*3+1)
      p3_pos_glob(iglob,2) = rbufr((iglob-1)*3+2)
      p3_pos_glob(iglob,3) = rbufr((iglob-1)*3+3)
    end do

  end subroutine get_pop3_pos_list


  !------------------------------------------------------------------------------------------------
  ! subroutine to remove pop3 candidates simultaneously forming in close vicinities (KS ADDED)
  !------------------------------------------------------------------------------------------------
  subroutine remove_p3cand_by_other_cand()

    !buffer for MPI send/recv (the size should be larger than both of ncpu and Npop3_max)
    integer,dimension(1:max(ncpu,Npop3_max*3)) :: sbufi,rbufi
    real(dp),dimension(1:max(ncpu,Npop3_max*3)) :: sbufr,rbufr


    !create p3cand_ncand_list
    p3cand_ncand_list(:) = 0
    p3cand_ncand_list(myid) = p3cand_ncand_loc
    sbufi(1:ncpu) = p3cand_ncand_list(:)
    call MPI_ALLREDUCE(sbufi,rbufi,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
    p3cand_ncand_list(:) = rbufi(1:ncpu)
    p3cand_ncand_tot = sum(p3cand_ncand_list) !total number of candidates

    if(myid==1 .and. p3cand_ncand_tot > 0) then
      write(*,'(A,I5,A)') 'Check simultaneous Pop III formation in close vicinity for',&
        p3cand_ncand_tot, ' candidates'
    end if


    ! !KS DEBUG
    ! write(*,*) '(KS DEBUG) myid, p3cand_ncand_tot: ', myid, p3cand_ncand_tot
    ! write(*,*) '(KS DEBUG) myid, p3cand_ncand_list(:): ',myid,p3cand_ncand_list(:)

    ! check the number of candidates
    if (p3cand_ncand_tot > Npop3_max) then
      write(*,*)'Number of global Pop III candidates exceeds the maximum ',p3cand_ncand_tot, Npop3_max
      write(*,*)'stopping...'
      stop
    end if

    !remove simultaneous close formation, if candiates exist
    if (p3cand_ncand_tot > 0) then
      !determine offset for global cndidate list
      p3cand_offset = 0
      do icpu = 1, myid-1  ! "do i = 1,0" loop is just skipped (KS MEMO)
        p3cand_offset = p3cand_offset + p3cand_ncand_list(icpu)
      end do

      !create global candidate position list
      p3cand_pos_glob(:,:)=0 !initialization
      do iloc = 1, p3cand_ncand_loc
        iglob = iloc + p3cand_offset
        p3cand_pos_glob(iglob,:) = p3cand_pos_loc(iloc,:)
      end do
      do iglob = 1, p3cand_ncand_tot
        sbufr((iglob-1)*3+1) = p3cand_pos_glob(iglob,1)
        sbufr((iglob-1)*3+2) = p3cand_pos_glob(iglob,2)
        sbufr((iglob-1)*3+3) = p3cand_pos_glob(iglob,3)
      end do
      call MPI_ALLREDUCE(sbufr, rbufr, 3*p3cand_ncand_tot, MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
      do iglob = 1, p3cand_ncand_tot
        p3cand_pos_glob(iglob,1) = rbufr((iglob-1)*3+1)
        p3cand_pos_glob(iglob,2) = rbufr((iglob-1)*3+2)
        p3cand_pos_glob(iglob,3) = rbufr((iglob-1)*3+3)
      end do

      ! check if pop3 star (upper in list) forms near current candidate
      do iloc = 1, p3cand_ncand_loc
        iglob = iloc + p3cand_offset
        ! rejected candidates are included in this loop but it should not significantly affect the results
        ! because not all candidates are removed (KS MEMO)
        do jglob = 1,iglob-1
          r2 = (p3cand_pos_glob(iglob,1) - p3cand_pos_glob(jglob,1))**2 + &
            (p3cand_pos_glob(iglob,2) - p3cand_pos_glob(jglob,2))**2 + &
            (p3cand_pos_glob(iglob,3) - p3cand_pos_glob(jglob,3))**2 !distance

          ! !KS DEBUG
          ! print *, 'iloc, iglob, jglob', iloc, iglob, jglob
          ! print '(A,6(1PE13.3))', "p3cand_pos_glob(iglob,:), p3cand_pos_glob(jglob,:)", p3cand_pos_glob(iglob,:), p3cand_pos_glob(jglob,:)
          ! print *, "r2, pop3_dist_min**2: ", r2, pop3_dist_min**2

          ! if too close, undo the pop3 formation
          if (r2 < pop3_dist_min**2) then
            ntot = ntot - 1
            mstar_tot = mstar_tot - flag3(p3cand_indcell_loc(iloc)) !flag3 store the mass of pop3
            flag3(p3cand_indcell_loc(iloc)) = 0

            print '(A,2(A,3(1PE12.5)),A,1PE12.5)', "** Remove Pop III cand by other cand **", &
              ' r_rmv(:): ', p3cand_pos_glob(iglob,:), &
              ', r_frm(:): ', p3cand_pos_glob(jglob,:),&
              ', r/rmin: ', sqrt(r2)/pop3_dist_min

            exit !jump out of the loop
          end if
        end do
      end do
    end if !p3cand_ncand_tot > 0
  end subroutine remove_p3cand_by_other_cand
  !----------------------------------------------------------------------------
  ! mdot_dependent_pop3_formation
  !----------------------------------------------------------------------------
  ! Create Pop III stars with mass depending on cloud-scale Mdot
  !  - calculate mass
  !  - add star particles
  !  - deplete gas
  subroutine mdot_dependent_pop3_formation
    character(len=80) :: prefix, nchar

    !--------------------------------------------------------------
    ! create global pop3 list from pop3 candidate list
    !-------------------------------------------------------------
    call new_pop3_list_from_candidate

    !leave if there is no new pop3 stars
    if (p3new_nform_tot == 0) return

    !--------------------------------------------------------------
    !update ufront if not updated
    !--------------------------------------------------------------
    if (.not. ufront_updated) then
      do sweep_level=levelmin,nlevelmax
        !call set_ufront(sweep_level)
        call set_ufront(sweep_level,flag_option=3) !use unew to set ufront
      end do
      ufront_updated = .true.
    end if

    !--------------------------------------------------------------
    ! output around Pop III formaiton site (TO: we will do this later)
    !-------------------------------------------------------------
    !do iglob = 1, p3new_nform_tot
    ! write(nchar,'(I8)') iglob
    !prefix = 'Pop3_'//TRIM(adjustl(nchar))
    !call output_around_star(prefix, p3new_pos_glob(:,iglob), use_current_ufront=.true.)
    !end do


    !--------------------------------------------------------------
    ! get Pop III mass
    !-------------------------------------------------------------
    call get_mdot_dependent_pop3_mass

    !--------------------------------------------------------------
    ! TO modified
    ! output around Pop III formaiton site
    !-------------------------------------------------------------
    do iglob = 1, p3new_nform_tot
      if (p3_output_minmass < 0 .or. &
        p3new_Mpop3(iglob)*scale_msun .ge. p3_output_minmass) then
        write(nchar,'(I8)') iglob
        prefix = 'Pop3_'//TRIM(adjustl(nchar))
        !call output_around_star_ufront(prefix, p3new_pos_glob(:,iglob))
        call output_around_star(prefix, p3new_pos_glob(:,iglob), use_current_ufront=.true.)
      end if
    end do

    !--------------------------------------------------------------
    ! add Pop III particles
    !-------------------------------------------------------------
    call add_pop3_particles

    !--------------------------------------------------------------
    ! take a log for formed stars
    !--------------------------------------------------------------

    !--------------------------------------------------------------
    ! deplete gas as a back reaction of star formation
    !-------------------------------------------------------------
    call deplete_gas_around_pop3

    !--------------------------------------------------------------
    ! undo pop3 formation in the original method
    !-------------------------------------------------------------
    do iloc = 1, p3cand_ncand_loc
      if (.not. flag3(p3cand_indcell_loc(iloc)) == 0) then     !limit to surviving ones in the candidate screening
        ntot = ntot - 1
        mstar_tot = mstar_tot - flag3(p3cand_indcell_loc(iloc)) !flag3 store the mass of pop3
        flag3(p3cand_indcell_loc(iloc)) = 0
      end if
    end do


  end subroutine mdot_dependent_pop3_formation


  !----------------------------------------------------------------------------
  ! new_pop3_list_from_candidate
  !  - create global pop3 list from pop3 candidate list
  !  - removal has been already made by remove_p3cand_by_other_cand
  !----------------------------------------------------------------------------
  subroutine new_pop3_list_from_candidate

    integer :: offset !offset between indices of local and global lists

    !buffer for MPI send/recv (the size should be larger than both of ncpu and Npop3_max)
    integer,dimension(1:max(ncpu,Npop3_max*3)) :: sbufi,rbufi
    real(dp),dimension(1:max(ncpu,Npop3_max*3)) :: sbufr,rbufr

    !-----------------------------------
    ! create global list of Pop III
    !-----------------------------------
    !Pop III surviving the removal process (local)
    p3new_nform_loc = 0
    do iloc = 1, p3cand_ncand_loc
      if (.not. flag3(p3cand_indcell_loc(iloc)) == 0) then     !limit to surviving ones
        p3new_nform_loc = p3new_nform_loc + 1           !number of Pop III after removal
        p3new_pos_loc(:,p3new_nform_loc) = p3cand_pos_loc(iloc,:) !position of Pop III (KS TODO: index order is opposite between p3new_pos_loc and p3cand_pos_loc)
      end if
    end do

    !list of Pop III formation number (global)
    p3new_nform_list(:) = 0
    p3new_nform_list(myid) = p3new_nform_loc
    sbufi(1:ncpu) = p3new_nform_list(:)
    call MPI_ALLREDUCE(sbufi,rbufi,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
    p3new_nform_list(:) = rbufi(1:ncpu)
    p3new_nform_tot = sum(p3new_nform_list) !total number of candidates


    !leave if there is no new pop3 stars
    if (p3new_nform_tot == 0) return

    !offset between local and global lists
    offset = 0
    do icpu = 1, myid-1  ! "do i = 1,0" loop is just skipped (KS MEMO)
      offset = offset + p3new_nform_list(icpu)
    end do

    !global Pop III position list
    p3new_pos_glob(:,:)=0 !initialization
    do iloc = 1, p3new_nform_loc
      iglob = iloc + offset
      p3new_pos_glob(:,iglob) = p3new_pos_loc(:,iloc)
    end do
    do iglob = 1, p3new_nform_tot
      sbufr((iglob-1)*3+1) = p3new_pos_glob(1,iglob)
      sbufr((iglob-1)*3+2) = p3new_pos_glob(2,iglob)
      sbufr((iglob-1)*3+3) = p3new_pos_glob(3,iglob)
    end do
    call MPI_ALLREDUCE(sbufr, rbufr, 3*p3new_nform_tot, MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    do iglob = 1, p3new_nform_tot
      p3new_pos_glob(1,iglob) = rbufr((iglob-1)*3+1)
      p3new_pos_glob(2,iglob) = rbufr((iglob-1)*3+2)
      p3new_pos_glob(3,iglob) = rbufr((iglob-1)*3+3)
    end do

    !KS DEBUG
    if (myid == 1) then
      print '(A,I5)', 'mdot_dependent_pop3new_formation:  p3new_nform_tot = ', p3new_nform_tot
      do iglob = 1, p3new_nform_tot
        print '(A,I5,1P3E15.7)', 'iglob, p3new_pos = ', iglob, p3new_pos_glob(:,iglob)
      end do
    end if
  end subroutine new_pop3_list_from_candidate

  !----------------------------------------------------------------------------
  ! get_mdot_dependent_pop3_mass
  !----------------------------------------------------------------------------
  ! Get Pop III mass from cloud-scale Mdot
  !  - based on star_formation_in_SFC
  !----------------------------------------------------------------------------
  subroutine get_mdot_dependent_pop3_mass

    integer,parameter :: Nrbin = 60 ! number of radial bins for shells (r_i = r_0 * (r_max/r_0)**(i/(Nrbin-1)))
    real(dp),parameter :: f_rmax = 1d3 !rmax/rmin: the ratio of minmum and maximum radii of radial bins (rmin is set to dx_min)


    ! !FOR DEBUG PURPOSE
    ! integer,parameter :: Nrbin = 20 ! number of radial bins for shells (r_i = r_0 * (r_max/r_0)**(i/(Nrbin-1)))
    ! real(dp),parameter :: f_rmax = 1d1 !rmax/rmin: the ratio of minmum and maximum radii of radial bins (rmin is set to dx_min)

    real(dp),parameter :: pi = 3.1415926d0

    real(dp),dimension(1:Nrbin) :: rbin      !outer radii of shells
    real(dp),dimension(1:Nrbin) :: rcen      !radii of center of shells
    real(dp),allocatable,dimension(:,:) :: Mshell    !mass of shells
    integer,allocatable,dimension(:,:) :: Nshell    !number of cell in shells
    real(dp),allocatable,dimension(:,:) :: Vshell    !volume of shells
    real(dp),allocatable,dimension(:,:) :: RHOshell   !density of shells
    real(dp),allocatable,dimension(:,:) :: MZshell    !metal mass of shells
    real(dp),allocatable,dimension(:,:) :: ETshell    !thermal energy of shells
    real(dp),allocatable,dimension(:,:) :: PRSshell   !pressure of shells
    real(dp),allocatable,dimension(:,:) :: MomRshell  !radial momentum of shells (vbulk subtracted)
    real(dp),allocatable,dimension(:,:) :: RawMomRshell  !radial momentum of shells (raw value)
    real(dp),allocatable,dimension(:,:) :: VelRshell  !radial velocity of shells
    real(dp),allocatable,dimension(:,:) :: MBEshell   !BE mass of shells
    real(dp),dimension(1:Nrbin) :: Menc, Venc, Menc_over_MBE  !enclosed mass, voulume, Menc/MBE

    real(dp) :: Menc_over_MBE_max, R_cloud, M_cloud, rho_cloud, velR_cloud, Mdot_cloud, Menc_over_MBE_cloud
    integer :: ir_cloud, ip3, ir
    real(dp), parameter:: yr = 3.15569d+07

    integer,parameter :: LUN = 11
    character(len=256),parameter :: SubDir  = 'SFdata' !sub directory for SFdata
    character(LEN=80)::filename,cnstep,cip3,command
    integer,save :: ifirst = 0 !for prcedures at first call

    !initialization
    if (ifirst == 0) then
      !if data directory does not exist, create it
      if(myid==1) then
        command = 'mkdir -p '//TRIM(SubDir)
        !call systemcall_wrapper(command)        !I don't know why, but this call sometimes fails (KS NOTE)
        !work around following the way in output_amr.f90
        call execute_command_line(TRIM(command))
      end if
      ifirst = 1
    endif

    ! Conversion factor from user units to cgs units
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    scale_msun=scale_d*scale_l**3/m_sun
    scale_pc=scale_l/3.08568d18
    factG=1
    if(cosmo)factG=3d0/8d0/pi*omega_m*aexp   ! (KS NOTE) = G*(scale_d*scale_t**2) (see amr/units.f90)

    !------------------------------------------
    ! construct 1D profile around each Pop III
    !------------------------------------------
    !define radial bins
    do ir = 1, Nrbin
      !outer radii of bins
      rbin(ir) = dx_min * (1d0-1d-8) * f_rmax**(dble(ir-1)/(Nrbin-1)) !put an extra factor (<1) to avoid ambiguity in defining shell regions
      !central radii of bins
      if (ir == 1) then
        rcen(ir)=0.5d0*rbin(ir)
      else
        rcen(ir)=0.5d0*(rbin(ir-1)+rbin(ir))
      end if
    end do

    !print *, "KS DEBUG (get_mdot_dependent_pop3_mass) myid = ", myid, ", p3new_nform_tot = ", p3new_nform_tot

    allocate(Mshell(Nrbin,p3new_nform_tot),Nshell(Nrbin,p3new_nform_tot),Vshell(Nrbin,p3new_nform_tot),RHOshell(Nrbin,p3new_nform_tot),&
      MZshell(Nrbin,p3new_nform_tot),ETshell(Nrbin,p3new_nform_tot),PRSshell(Nrbin,p3new_nform_tot),MomRshell(Nrbin,p3new_nform_tot),&
      RawMomRshell(Nrbin,p3new_nform_tot),VelRshell(Nrbin,p3new_nform_tot),MBEshell(Nrbin,p3new_nform_tot))
    !get shell quantities around each Pop III (WARNING: overlap of PopIII is not allowed in the current version of code)
    call get_shell_quantities(p3new_nform_tot, p3new_pos_glob(1:3,1:p3new_nform_tot), Nrbin, rbin(1:Nrbin), Mshell(1:Nrbin,1:p3new_nform_tot),&
      Nshell=Nshell(1:Nrbin,1:p3new_nform_tot), Vshell=Vshell(1:Nrbin,1:p3new_nform_tot), &
      MZshell=MZshell(1:Nrbin,1:p3new_nform_tot), ETshell=ETshell(1:Nrbin,1:p3new_nform_tot),&
      MomRshell=MomRshell(1:Nrbin,1:p3new_nform_tot))
    ! call get_shell_quantities(p3new_nform_tot, p3new_pos_glob(1:3,1:p3new_nform_tot), Nrbin, rbin(1:Nrbin), Mshell(1:Nrbin,1:p3new_nform_tot),&
    !      Nshell=Nshell(1:Nrbin,1:p3new_nform_tot), Vshell=Vshell(1:Nrbin,1:p3new_nform_tot), &
    !      MZshell=MZshell(1:Nrbin,1:p3new_nform_tot), ETshell=ETshell(1:Nrbin,1:p3new_nform_tot),&
    !      MomRshell=MomRshell(1:Nrbin,1:p3new_nform_tot),RawMomRshell=RawMomRshell(1:Nrbin,1:p3new_nform_tot))


    !get RHOshell, PRSshell, VelRshell
    do ip3 = 1, p3new_nform_tot
      do ir = 1, Nrbin
        if (Nshell(ir,ip3)>0) then  !if there is at least one cell in the shell
          RHOshell(ir,ip3) = Mshell(ir,ip3)/Vshell(ir,ip3)
          PRSshell(ir,ip3) = (gamma-1.0) * ETshell(ir,ip3)/Vshell(ir,ip3)
          VelRshell(ir,ip3) = MomRshell(ir,ip3)/Mshell(ir,ip3)
          MBEshell(ir,ip3) = &
            1.18d0 * PRSshell(ir,ip3)**1.5 / RHOshell(ir,ip3)**2 / factG**1.5  !BE mass
        else                         !if thre is no cell in the shell
          RHOshell(ir,ip3) = -1d0  !zero for ILL-DEFINED case
          PRSshell(ir,ip3) = -1d0
          VelRshell(ir,ip3) = -1d0
          MBEshell(ir,ip3) = -1d0
        end if
      end do
    end do

    !KS DEBUG
    if (myid == 1) then
      do ip3 = 1, p3new_nform_tot
        if (ip3 == 1) then
          print '(A,I5)', '(get_mdot_dependent_pop3_mass) p3new_nform_tot:', p3new_nform_tot
          print '(A,I5)', 'Nrbin: ', Nrbin
        end if
        print '(A,1P3E20.12)', 'p3new_pos(:): ', p3new_pos_glob(1:3,ip3)
        print '(A,1P100E10.2)', 'rbin(:):     ', rbin(:)*scale_pc
        print '(A,100I9)',     'Nshell(:):   ', Nshell(1:Nrbin,ip3)
        print '(A,1P100E10.2)', 'Mshell(:):   ', Mshell(1:Nrbin,ip3)*scale_msun
        print '(A,1P100E10.2)', 'Vshell(:):   ', Vshell(1:Nrbin,ip3)*scale_pc**3
        print '(A,1P100E10.2)', 'RHOshell(:): ', RHOshell(1:Nrbin,ip3)*scale_nH
        print '(A,1P100E10.2)', 'MZshell(:):  ', MZshell(1:Nrbin,ip3)*scale_msun
        print '(A,1P100E10.2)', 'T2shell(:):  ', PRSshell(1:Nrbin,ip3)/RHOshell(1:Nrbin,ip3)*scale_T2
        print '(A,1P100E10.2)', 'MomRshell(:):', MomRshell(1:Nrbin,ip3)
        print '(A,1P100E10.2)', 'VelRshell(:):', VelRshell(1:Nrbin,ip3)*scale_v/1d5
        print '(A,1P100E10.2)', 'MBEshell(:): ', MBEshell(1:Nrbin,ip3)*scale_msun
      end do
    end if


    !----------------------------------------------------------------------------
    ! get cloud size and cloud-scale accretion rate
    !----------------------------------------------------------------------------
    !loop over PopIII
    do ip3 = 1, p3new_nform_tot
      !calculate Menc
      Menc(:) = 0d0
      do ir = 1, Nrbin
        if (ir == 1) then
          Menc(ir) = Mshell(1,ip3)
        else
          Menc(ir) = Menc(ir-1) + Mshell(ir,ip3)
        end if
      end do

      !Menc_over_MBE (negative value for no-data shell)
      Menc_over_MBE(:) = 0d0
      do ir = 1, Nrbin
        if (Nshell(ir,ip3)>0) then
          Menc_over_MBE(ir) = Menc(ir)/MBEshell(ir,ip3)
        else
          Menc_over_MBE(ir) = -1d0
        end if
      end do

      !determine cloud size from maximum Menc_over_MBE
      Menc_over_MBE_max = 0
      ir_cloud = -1
      do ir = 1, Nrbin
        !if (Menc_over_MBE(ir) < 0) cycle
        if (Nshell(ir,ip3) == 0) cycle
        if (VelRshell(ir, ip3) >= 0) cycle
        if (Menc_over_MBE(ir) > Menc_over_MBE_max) then
          Menc_over_MBE_max = Menc_over_MBE(ir)
          ir_cloud = ir
        end if
      end do

      !check if radius is found
      if (ir_cloud < 0) then
        if (myid == 1) then
          print '(A,I5,A,1PE11.3/,A)', '*** ERROR *** PopIII cloud radius not found in get_mdot_dependent_pop3new_mass. ir = ', ir, ', stopping...'
          stop
        end if
      end if

      !check if SFC touches the boundary (TODO: compatibility with periodic boundary)
      if (VelRshell(ir_cloud, ip3) < 0 .and. (MINVAL(p3new_pos_glob(1:3,ip3))-rbin(ir_cloud) < 0d0 .or. &
      & MAXVAL(p3new_pos_glob(1:3,ip3))+rbin(ir_cloud) > boxlen)) then   ! merge_ks: boxlen (ks: 1d0)
        print '(A,1P4E15.6,/,A)', &
          '*** ERROR *** PopIII cloud may extend beyond the boundary (ToDo: implementation of periodic boundary)', &
          p3new_pos_glob(1:3,ip3),rbin(ir_cloud)*scale_pc, &
          'stopping...'
        stop
      end if

      !cloud-scale variables
      !R_cloud = rbin(ir_cloud)
      R_cloud = rcen(ir_cloud) !values at the center of radial bins
      M_cloud = Menc(ir_cloud)
      Menc_over_MBE_cloud = Menc_over_MBE(ir_cloud)
      rho_cloud = RHOshell(ir_cloud,ip3)
      VelR_cloud = VelRshell(ir_cloud,ip3)
      Mdot_cloud = -4*pi*R_cloud**2*rho_cloud*VelR_cloud


      !calculate Mpop3 from Mdot_cloud
      if (Mdot_cloud > 0) then
        p3new_Mpop3(ip3) = Mpop3_from_Mdot_H15(Mdot_cloud)
      else
        p3new_Mpop3(ip3) = 0
      endif

      !for pop3_model==6, skip Pop III formation entirely (no star, no gas depletion) if the total
      !mass would be below M_multi_min (the minimum mass of a multiple system in add_pop3_particles),
      !instead of forming a single star lighter than that floor
      if (pop3_model == 6 .and. p3new_Mpop3(ip3)*scale_msun < 5d0) then
        p3new_Mpop3(ip3) = 0
      end if

      !print summary
      if (myid == 1) then
        print '(2(A,I5),3(A,1PE15.6),/,4(A,1PE15.6),A)', &
          'Pop III mass from cloud-scale Mdot: ip3 = ', ip3, &
          ', ir_cloud = ', ir_cloud, &
          ', R_cloud = ', R_cloud*scale_pc,  &
          ' [pc], M_cloud = ', M_cloud*scale_msun, &
          ' [Msun], (Menc/MBE)_cloud = ', Menc_over_MBE_cloud,&
          'VelR_cloud = ', VelR_cloud*scale_v/1d5, &
          ' [km/s], nH_cloud = ', rho_cloud*scale_nH, &
          ' [cm-3], Mdot_cloud = ', Mdot_cloud*scale_msun/(scale_t/yr), &
          ' [Msun/yr], Mpop3 = ', p3new_Mpop3(ip3)*scale_msun,' [Msun]'
      end if
      !to make output easier to read
      call flush(6)
      call MPI_BARRIER(MPI_COMM_WORLD,info)


      !------------------------------------------------------------
      ! save summary and 1D profile to a file
      !------------------------------------------------------------
      if (myid == 1 .and. p3new_Mpop3(ip3) > 0) then
        ! get filename (SFdata/pop3prof.[step].[ip3]..dat)
        write(cnstep,'(I8)') nstep
        write(cip3,'(I8)') ip3
        filename = TRIM(SubDir) // '/pop3prof.' // TRIM(ADJUSTL(cnstep)) // &
          '.' // TRIM(ADJUSTL(cip3)) // '.dat'

        print *, "outputting ", filename

        !open data file
        open(LUN, file=filename, status='replace')

        !write metadata
        write (LUN, '(A,1PE14.6,A,1PE14.6,A)') &
          "#Pop III formation at aexp=", aexp
        write (LUN, '(A,1P3E20.12)') &
          "# pos(:)=", p3new_pos_glob(1:3,ip3)
        write (LUN, '(A,I8,3(A,1PE14.6),/,4(A,1PE14.6),A)') &
          '# ir_cloud = ', ir_cloud, &
          ', R_cloud = ', R_cloud*scale_pc,  &
          ' [pc], M_cloud = ', M_cloud*scale_msun, &
          ' [Msun], (Menc/MBE)_cloud = ', Menc_over_MBE_cloud,&
          ' VelR_cloud = ', VelR_cloud*scale_v/1d5, &
          ' [km/s], nH_cloud = ', rho_cloud*scale_nH, &
          ' [cm-3], Mdot_cloud = ', Mdot_cloud*scale_msun/(scale_t/(1d6*365.*24.*3600.)), &
          ' [Msun/yr], Mpop3 = ', p3new_Mpop3(ip3)*scale_msun,' [Msun]'

        !write header
        write (LUN, '(A)') &
          "#i, rcen[pc], rbin[pc], Nshell, Mshell[Msun], Vshell[pc^3], RHOshell[cm-3], MZshell[Msun], T2shell[K], VelRshell[km/s], MBEshell[Msun], Mdotshell[Msun/yr], Menc/MBE"

        !write data
        do ir = 1, Nrbin
          if (Nshell(ir,ip3)>0) then
            write(LUN, '(I5,1P2E16.8,I8,9(1PE16.8))') &
              ir, rcen(ir)*scale_pc, rbin(ir)*scale_pc, Nshell(ir,ip3), Mshell(ir,ip3)*scale_msun, &
              Vshell(ir,ip3)*scale_pc**3, RHOshell(ir,ip3)*scale_nH, &
              MZshell(ir,ip3)*scale_msun, PRSshell(ir,ip3)/RHOshell(ir,ip3)*scale_T2, &
              VelRshell(ir,ip3)*scale_v/1d5, MBEshell(ir,ip3)*scale_msun, &
              -4*pi*rbin(ir)**2*RHOshell(ir,ip3)*VelRshell(ir,ip3)*scale_msun/(scale_t/(1d6*365.*24.*3600.)),&
              Menc(ir)/MBEshell(ir,ip3)
          else
            write(LUN, '(I5,1P2E16.8,I8,9(1PE16.8))') &
              ir, rcen(ir)*scale_pc, rbin(ir)*scale_pc, Nshell(ir,ip3), (/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/) !no-data cells
          end if
        end do
      end if
    end do


  end subroutine get_mdot_dependent_pop3_mass

  !----------------------------------------------------------------------------
  ! add pop3 particles based on pop3_list
  !
  !----------------------------------------------------------------------------
  ! - deplete gas from spherical regions around given poistion  with constant factor f_dpl
  ! - the radius of the sphere is determined so as to contain larger mass than the star
  ! - this subroutine is based partly on star_formation_in_SFC

  subroutine add_pop3_particles

    integer,dimension(1:nvector)::cell_index,cell_levl,cc
    integer::cell_index1,cell_levl1
    real(dp),dimension(1:3)::pos
    integer :: ip3
    !MPI send/recv buffer
    integer,dimension(1:ncpu) :: sbufi,rbufi

    !for multiple Pop III formation
    integer :: i_multi, n_multi
    real(dp),parameter :: M_multi_min = 5d0 !minimum mass in a multiple system
    real(dp) :: M_i, M_n

    p3new_nstar_loc=0
    !loop over PopIII
    do ip3 = 1, p3new_nform_tot
      if (p3new_Mpop3(ip3) > 0) then
        !number of stars in the multile case
        if (pop3_model == 6) then
          !determine n (M_tot = Mn*(2**n-1)   -> n = log((M_tot/Mn+1))/log(2)
          !(KS MEMO) p3new_Mpop3(ip3) is guaranteed to be >= M_multi_min here, since M_tot below
          !M_multi_min is already zeroed out (skipped) in get_mdot_dependent_pop3_mass
          n_multi = int(log(p3new_Mpop3(ip3)*scale_msun/M_multi_min + 1)/log(2d0))
          M_n = p3new_Mpop3(ip3) / (2**n_multi-1)
        end if



        !get cpu responsible for pos
        pos(1:3) = p3new_pos_glob(1:3,ip3)
        !call cmp_cpumap(pos,cc,1)
        call get_active_cell_index_highest_level(pos,cell_levl1,cell_index1)
        !if (cc(1) .eq. myid) then
        if (cell_index1 >= 0) then           !when active grid in the highest level is found on this cpu
          !------------------------------------------------
          !single Pop III star formation case
          !------------------------------------------------
          if (pop3_model == 4) then
            p3new_nstar_loc = p3new_nstar_loc + 1

            !get index of the cell where SFC particle is located
            !call get_cell_index(cell_index,cell_levl,pos,nlevelmax,1)
            !corresponding grid index
            ind_grid(1)=mod(cell_index1-ncoarse-1,ngridmax)+1 !cell -> grid where the cell belongs

            ! Update linked list for stars
            call remove_free(ind_part,1)
            call add_list(ind_part,ind_grid,ok_new,1)  !ok_new are all .true.

            ! Set star particle variables
            index_star=nstar_tot+1         ! particle index
            tp(ind_part(1))=birth_epoch  ! Birth epoch
            mp(ind_part(1))= p3new_Mpop3(ip3)   ! Mass
            levelp(ind_part(1))=cell_levl1   ! Level
            idp(ind_part(1))=index_star  ! Star identity (negative integers will be allocated to sinks)
            typep(ind_part(1))%family=FAM_POP3  ! Particle type
            xp(ind_part(1),1)=pos(1)
            xp(ind_part(1),2)=pos(2)
            xp(ind_part(1),3)=pos(3)

            !set initial star velocity to gas velocity of the host cell
            ! use ufront
            vp(ind_part(1),1)=ufront(cell_index1,2)/ufront(cell_index1,1)
            vp(ind_part(1),2)=ufront(cell_index1,3)/ufront(cell_index1,1)
            vp(ind_part(1),3)=ufront(cell_index1,4)/ufront(cell_index1,1)
            ! TO DEBUG
            !zp(ind_part(1)) = ufront(cell_index1,imetal)
            zp(ind_part(1)) = ufront(cell_index1,imetal)/ufront(cell_index1,1)
#ifdef INIT_STELLAR_MASS
            mp0(ind_part(1))=mp(ind_part(1))
#endif
#ifdef STELLAR_POPULATION_MASS
            msp0(ind_part(1))=mp(ind_part(1))
#endif
            if(sf_log_properties) call sf_log_write_star(ilun_sf,ind_part(1),cell_index1,ilevel,(cell_levl1==nlevelmax))

            write(*,'(A,2I4,I10,A,I8,A,I2,A,1PE8.1,A,1P3E9.1, A,1P3E16.8)') &
              "** New Pop3 particle added **  myid, lev, ind_cell: ", myid, cell_levl1, cell_index1,&
              ", idp:",  idp(ind_part(1)),&
              ", fm:", typep(ind_part(1))%family,&
              ", radius:", -zp(ind_part(1))*scale_pc,&
              ", vp:", vp(ind_part(1),:)*scale_l/scale_t/1d5,&
              ", xp:", xp(ind_part(1),:)

            !------------------------------------------------
            ! binary Pop III star (m1:m2 = 2:1) formation case
            !------------------------------------------------
          else if (pop3_model == 5) then
            p3new_nstar_loc = p3new_nstar_loc + 2 !add two stars (binary)

            !get index of the cell where SFC particle is located
            call get_cell_index(cell_index,cell_levl,pos,nlevelmax,1)
            !corresponding grid index
            ind_grid(1)=mod(cell_index1-ncoarse-1,ngridmax)+1 !cell -> grid where the cell belongs
            ind_grid(2)=mod(cell_index1-ncoarse-1,ngridmax)+1 !cell -> grid where the cell belongs

            ! Update linked list for stars
            call remove_free(ind_part,2)
            call add_list(ind_part,ind_grid,ok_new,2)   !ok_new are all .true.


            !------------------------------------------------------
            ! Set star particle variables for the first binary star
            !------------------------------------------------------
            index_star=nstar_tot+1         ! particle index
            tp(ind_part(1))=birth_epoch  ! Birth epoch
            mp(ind_part(1))= 2d0/3d0 * p3new_Mpop3(ip3)   ! Mass (2/3 of total mass)
            levelp(ind_part(1))=cell_levl1   ! Level
            idp(ind_part(1))=index_star  ! Star identity (negative integers will be allocated to sinks)
            typep(ind_part(1))%family=FAM_POP3  ! Particle type
            xp(ind_part(1),1)=pos(1)
            xp(ind_part(1),2)=pos(2)
            xp(ind_part(1),3)=pos(3)

            !set initial star velocity to gas velocity of the host cell
            ! use ufront
            vp(ind_part(1),1)=ufront(cell_index1,2)/ufront(cell_index1,1)
            vp(ind_part(1),2)=ufront(cell_index1,3)/ufront(cell_index1,1)
            vp(ind_part(1),3)=ufront(cell_index1,4)/ufront(cell_index1,1)
            ! TO DEBUG
            !zp(ind_part(1)) = ufront(cell_index1,imetal)
            zp(ind_part(1)) = ufront(cell_index1,imetal)/ufront(cell_index1,1)
#ifdef INIT_STELLAR_MASS
            mp0(ind_part(1))=mp(ind_part(1))
#endif
#ifdef STELLAR_POPULATION_MASS
            msp0(ind_part(1))=mp(ind_part(1))
#endif
            if(sf_log_properties) call sf_log_write_star(ilun_sf,ind_part(1),cell_index1,ilevel,(cell_levl1==nlevelmax))

            write(*,'(A,2I4,I10,A,I8,A,I2,A,1PE9.2,A,1PE8.1,A,1P3E9.1, A,1P3E16.8)') &
              "** New Pop3 particle added (one of binary with m1 = 2/3 mtot) **  myid, lev, ind_cell: ", myid, cell_levl1, cell_index1,&
              ", idp:",  idp(ind_part(1)),&
              ", family:", typep(ind_part(1))%family,&
              ", mp:", mp(ind_part(1))*scale_msun,&
              ", zp:", zp(ind_part(1))/0.02,&
              ", vp:", vp(ind_part(1),:)*scale_l/scale_t/1d5,&
              ", xp:", xp(ind_part(1),:)

            !------------------------------------------------------
            ! Set star particle variables for the second binary star
            !------------------------------------------------------
            index_star=nstar_tot+2         ! particle index
            tp(ind_part(2))=birth_epoch  ! Birth epoch
            mp(ind_part(2))= 1d0/3d0 * p3new_Mpop3(ip3)   ! Mass (1/3 of total mass)
            levelp(ind_part(2))=cell_levl1   ! Level
            idp(ind_part(2))=index_star  ! Star identity (negative integers will be allocated to sinks)
            typep(ind_part(2))%family=FAM_POP3  ! Particle type
            xp(ind_part(2),1)=pos(1)
            xp(ind_part(2),2)=pos(2)
            xp(ind_part(2),3)=pos(3)

            !set initial star velocity to gas velocity of the host cell
            ! use ufront
            vp(ind_part(2),1)=ufront(cell_index1,2)/ufront(cell_index1,1)
            vp(ind_part(2),2)=ufront(cell_index1,3)/ufront(cell_index1,1)
            vp(ind_part(2),3)=ufront(cell_index1,4)/ufront(cell_index1,1)
            ! TO DEBUG
            ! zp(ind_part(2)) = ufront(cell_index1,imetal)
            zp(ind_part(2)) = ufront(cell_index1,imetal)/ufront(cell_index1,1)
#ifdef INIT_STELLAR_MASS
            mp0(ind_part(2))=mp(ind_part(2))
#endif
#ifdef STELLAR_POPULATION_MASS
            msp0(ind_part(2))=mp(ind_part(2))
#endif
            if(sf_log_properties) call sf_log_write_star(ilun_sf,ind_part(2),cell_index1,ilevel,(cell_levl1==nlevelmax))

            write(*,'(A,2I4,I10,A,I8,A,I2,A,1PE9.2,A,1PE8.1,A,1P3E9.1, A,1P3E16.8)') &
              "** New Pop3 particle added (one of binary with m2 = 1/3 mtot) **  myid, lev, ind_cell: ", myid, cell_levl1, cell_index1,&
              ", idp:",  idp(ind_part(2)),&
              ", family:", typep(ind_part(2))%family,&
              ", mp:", mp(ind_part(2))*scale_msun,&
              ", zp:", zp(ind_part(2))/0.02,&
              ", vp:", vp(ind_part(2),:)*scale_l/scale_t/1d5,&
              ", xp:", xp(ind_part(2),:)

            !------------------------------------------------
            ! case of multiple Pop III star M1:M2:..:Mn = 2^n:2^n-1:...:1,
            ! with maximum n satisfying Mn > 10 Msun
            !------------------------------------------------
          else if (pop3_model == 6) then
            do i_multi = 1, n_multi
              !Mass of i_th star in the multiple system
              !M_i = 2**(n_multi-i_multi) * M_n  !causes error if 2**(n_multi-i_multi)  > 2,147,483,647 (maximum integer(4))
              M_i = 2**(n_multi-i_multi) * M_n

              !add star one-by-one
              p3new_nstar_loc = p3new_nstar_loc + 1

              !get index of the cell where SFC particle is located
              !call get_cell_index(cell_index,cell_levl,pos,nlevelmax,1)
              !corresponding grid index
              ind_grid(1)=mod(cell_index1-ncoarse-1,ngridmax)+1 !cell -> grid where the cell belongs

              ! Update linked list for stars
              call remove_free(ind_part,1)
              call add_list(ind_part,ind_grid,ok_new,1)  !ok_new are all .true.

              ! Set star particle variables
              index_star=nstar_tot+i_multi         ! particle index
              tp(ind_part(1))=birth_epoch  ! Birth epoch
              mp(ind_part(1))= M_i   ! Mass
              levelp(ind_part(1))=cell_levl1   ! Level
              idp(ind_part(1))=index_star  ! Star identity (negative integers will be allocated to sinks)
              typep(ind_part(1))%family=FAM_POP3  ! Particle type
              xp(ind_part(1),1)=pos(1)
              xp(ind_part(1),2)=pos(2)
              xp(ind_part(1),3)=pos(3)

              !set initial star velocity to gas velocity of the host cell
              ! use ufront
              vp(ind_part(1),1)=ufront(cell_index1,2)/ufront(cell_index1,1)
              vp(ind_part(1),2)=ufront(cell_index1,3)/ufront(cell_index1,1)
              vp(ind_part(1),3)=ufront(cell_index1,4)/ufront(cell_index1,1)
              ! TO DEBUG
              ! zp(ind_part(1)) = ufront(cell_index1,imetal)
              zp(ind_part(1)) = ufront(cell_index1,imetal)/ufront(cell_index1,1)
#ifdef INIT_STELLAR_MASS
              mp0(ind_part(1))=mp(ind_part(1))
#endif
#ifdef STELLAR_POPULATION_MASS
              msp0(ind_part(1))=mp(ind_part(1))
#endif
              if(sf_log_properties) call sf_log_write_star(ilun_sf,ind_part(1),cell_index1,ilevel,(cell_levl1==nlevelmax))

              write(*,'(2(A,I4),A,2I4,I10,A,I8,A,I2,A,1PE8.1,/,A,1P3E9.1, A,1P3E16.8)') &
                "** New Pop3 particle added (", i_multi, "th star in the ", &
                n_multi,  "-star system) **  myid, lev, ind_cell: ", &
                myid, cell_levl1, cell_index1,&
                ", idp:",  idp(ind_part(1)),&
                ", fm:", typep(ind_part(1))%family,&
                ", mp:", mp(ind_part(1))*scale_msun,&
                ", vp:", vp(ind_part(1),:)*scale_l/scale_t/1d5,&
                ", xp:", xp(ind_part(1),:)
            end do
          end if !pop3_model

        end if !if (cc(1) .eq. myid) then

        !increase total number of stellar paritlces
        if (pop3_model == 4) then
          nstar_tot =  nstar_tot + 1 !single star
        else if (pop3_model == 5) then
          nstar_tot =  nstar_tot + 2 !binary stars
        else if (pop3_model == 6) then
          nstar_tot =  nstar_tot + n_multi !multiple stars
        end if
      end if
    end do !do ip3 = 1, p3new_nform_tot


    !------------------------------------------------------------------------------
    !list of Pop III star number counting individual multiple stars (global)
    !------------------------------------------------------------------------------
    p3new_nstar_list(:) = 0
    p3new_nstar_list(myid) = p3new_nstar_loc
    sbufi(1:ncpu) = p3new_nstar_list(:)
    call MPI_ALLREDUCE(sbufi,rbufi,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
    p3new_nstar_list(:) = rbufi(1:ncpu)
    p3new_nstar_tot = sum(p3new_nstar_list) !total number of candidates

  end subroutine add_pop3_particles

  !----------------------------------------------------------------------------
  ! deplete gas around pop3 stars
  !----------------------------------------------------------------------------
  subroutine deplete_gas_around_pop3
    integer :: ip3
    do ip3 = 1, p3new_nform_tot
      if (p3new_Mpop3(ip3) > 0) then
        call deplete_gas_around_pos(p3new_Mpop3(ip3),p3new_pos_glob(:,ip3))
      endif
    end do
  end subroutine deplete_gas_around_pop3

  !----------------------------------------------------------------------------
  ! deplete gas around given position
  !----------------------------------------------------------------------------
  ! - deplete gas from spherical regions around given poistion with constant factor f_dpl
  ! - the radius of the sphere is determined so as to contain larger mass than the star
  !----------------------------------------------------------------------------
  subroutine deplete_gas_around_pos(mass,pos)

    !parameters
    real(dp),parameter :: f_dpl_max = 0.2d0 ! depletion radius is determiend as the minimum radius that gives f_dpl smaller than f_dpl_max
    ! integer,parameter :: Nrbin = 40 ! number of radial bin for defining shells (r_i/dx_min = 1,2,4,,,2**(Nrbin-1))
    ! real(dp),parameter :: f_rmax = 1d2 !rmax/rmin: the ratio of minmum and maximum radii of radial bins (rmin is set to dx_min)

    !increase the search radius since the size is sometimes not enough (KS: 2025/07/12)
    integer,parameter :: Nrbin = 60 ! number of radial bin for defining shells (r_i/dx_min = 1,2,4,,,2**(Nrbin-1))
    real(dp),parameter :: f_rmax = 1d3 !rmax/rmin: the ratio of minmum and maximum radii of radial bins (rmin is set to dx_min)



    !variables
    real(dp),dimension(1:Nrbin) :: rbin      !outer radii of shells
    real(dp),dimension(1:Nrbin) :: Mshell    !mass of shells
    integer,dimension(1:Nrbin) :: Nshell    !number of cell in shells
    real(dp),dimension(1:Nrbin) :: Vshell    !volume of shells
    real(dp),dimension(1:Nrbin) :: Menc      !enclosed mass
    real(dp) :: M_sphere, r_sphere
    real(dp) :: f_dpl ! depletion factor = M_dpl/M_gas (0 -> no depletion, 1 -> complete deplition)

    real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v, scale_msun, scale_pc
    integer :: ir

    !for debug
    real(dp),dimension(1:Nrbin) :: Mshell_before, Mshell_after    !mass of shells before/after depletion
    real(dp),dimension(1:Nrbin) :: Menc_before, Menc_after      !enclosed mass before/after depletion

    real(dp),intent(in) :: mass
    real(dp),dimension(1:3),intent(in) :: pos

    ! Conversion factor from user units to cgs units
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    scale_msun=scale_d*scale_l**3/m_sun
    scale_pc=scale_l/3.08568d18

    !define radial bins
    do ir = 1, Nrbin
      !rbin(ir) = dx_min * (1d0-1d-8) * 2**(ir-1) !put an extra factor (<1) to avoid ambiguity in defining shell regions
      rbin(ir) = dx_min * (1d0-1d-8) * f_rmax**(dble(ir-1)/(Nrbin-1)) !put an extra factor (<1) to avoid ambiguity in defining shell regions
    end do


    call get_shell_quantities(1, pos(1:3), Nrbin, rbin, Mshell(1:Nrbin),&
      Nshell=Nshell(1:Nrbin), Vshell=Vshell(1:Nrbin))


    !calculate Menc
    Menc(:) = 0d0
    do ir = 1, Nrbin
      if (ir == 1) then
        Menc(ir) = Mshell(1)
      else
        Menc(ir) = Menc(ir-1) + Mshell(ir)
      end if
    end do

    !determine size of depletion sphere (minimum radius that saisfies mass/M_sphere < f_dpl_max)
    do ir = 1, Nrbin
      if (Menc(ir) > mass/f_dpl_max) then
        r_sphere = rbin(ir)
        M_sphere = Menc(ir)
        !set f_dpl
        f_dpl = mass / Menc(ir)
        exit
      end if
      !最後までMencがmstar/f_dpl_maxを超えられなかった場合はエラーを吐いて止める
      if (ir == Nrbin) then
        if (myid == 1) then
          print *, '*** ERROR *** (deplete_gas_around_pos) mass is too large compared with M_sphere'
          print '(A,1P2E12.4)', 'Menc(ir), mass: ', Menc(ir), mass
          print *, 'ir, Mshell(:) for debug purpose: ', ir, Mshell(:)
          print *, 'Menc(:): ', Menc(:)
          print *, 'something strange may happen, or change Nrbin, f_rmax, or f_dpl_max. stopping...'
          ! stop ! TO MODIFIED: moved out of myid==1 block (2026-10-30)
        end if
        call clean_stop ! TO MODIFIED: condition is global; stop on all ranks (other ranks used to continue with unset values) (2026-10-30)
      end if
    end do

    !KS DEBUG
    if (myid == 1) then
      print '(A, I6)', '(deplete_gas_around_pos) Nrbin: ', Nrbin
      print '(A, 2I6, 1P5E13.5)', 'istar, ir, r_sphere, M_sphere, mass, f_dpl:', &
        istar,  ir, r_sphere*scale_pc, M_sphere*scale_msun, mass*scale_msun, f_dpl
      print '(A, 1P100E9.2)', 'rbin(:) [pc]:   ', rbin(:)*scale_pc
      print '(A, 100I9)',     'Nshell(:):      ', Nshell(:)
      print '(A, 1P100E9.2)', 'Menc(:) [Msun]: ', Menc(:)*scale_msun
    end if


    !deplte gas by f_dpl
    call modify_sphere_by_factor(pos(1:3),r_sphere,1d0-f_dpl)  !argument is f_rho = 1 - f_dpl

  end subroutine deplete_gas_around_pos


  !----------------------------------------------------------------------------
  ! Mpop3_from_Mdot_H15
  !----------------------------------------------------------------------------
  ! calculate Pop III star mass using the formula given by Hirano+15
  !----------------------------------------------------------------------------
  function Mpop3_from_Mdot_H15(Mdot_cloud)
    implicit none
    real(dp), intent(IN):: Mdot_cloud
    real(dp), parameter:: yr = 3.15569d+07
    real(dp) :: Mpop3_from_Mdot_H15
    real(dp) :: Mdot_cloud_Msun_yr

    ! ! Conversion factor from user units to cgs units
    ! call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    ! scale_msun=scale_d*scale_l**3/m_sun

    Mdot_cloud_Msun_yr = Mdot_cloud*scale_msun/(scale_t/yr)

    Mpop3_from_Mdot_H15 = 2.5e2 * (Mdot_cloud_Msun_yr/2.8d-3)**0.7d0 / scale_msun

  end function Mpop3_from_Mdot_H15

  !------------------------------ KS ADDED END ----------------------------------------!







  !----------------------------------------------------------------------------
  ! deplete gas from spheres
  !----------------------------------------------------------------------------
  ! - deplete gas from spherical regions around newly formed stars with constant factor f_dpl
  ! - the radius of the sphere is determined so as to contain larger mass than the star
  !----------------------------------------------------------------------------
  subroutine deplete_gas_from_spheres

    !parameters
    real(dp),parameter :: f_dpl_max = 0.2d0 ! depletion radius is determiend as the minimum radius that gives f_dpl smaller than f_dpl_max
    integer,parameter :: Nrbin = 40 ! number of radial bin for defining shells (r_i/dx_min = 1,2,4,,,2**(Nrbin-1))
    real(dp),parameter :: f_rmax = 1d2 !rmax/rmin: the ratio of minmum and maximum radii of radial bins (rmin is set to dx_min)

    !variables
    real(dp),dimension(1:Nnew_max) :: mass_new     !mass of newly formed stars
    real(dp),dimension(1:3,1:Nnew_max) :: pos_new  !position of newly formed stars
    real(dp),dimension(1:Nrbin) :: rbin      !outer radii of shells
    real(dp),allocatable,dimension(:,:) :: Mshell    !mass of shells
    integer,allocatable,dimension(:,:) :: Nshell    !number of cell in shells
    real(dp),allocatable,dimension(:,:) :: Vshell    !volume of shells
    real(dp),dimension(1:Nrbin) :: Menc      !enclosed mass
    real(dp) :: M_sphere, r_sphere
    real(dp) :: f_dpl ! depletion factor = M_dpl/M_gas (0 -> no depletion, 1 -> complete deplition)

    real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v, scale_msun, scale_pc
    integer :: ir, istar

    !for debug
    real(dp),allocatable,dimension(:,:) :: Mshell_before, Mshell_after    !mass of shells before/after depletion
    real(dp),dimension(1:Nrbin) :: Menc_before, Menc_after      !enclosed mass before/after depletion

    !Note:
    ! ntot_star_all(icpu) is # of new stars in icpu
    ! ntot_all is sum of ntot_star_all(:)

    !skip if there is no new stars
    if (ntot_all == 0) return

    !check number of newly formed stars
    if (ntot_all > Nnew_max) then
      if (myid == 1) then
        print '(A,2I8)', 'number of newly formed stars exceeds the maximum value',ntot_all, Nnew_max
        print '(A)', 'stopping...'
        ! stop ! TO MODIFIED: moved out of myid==1 block (2026-10-30)
      end if
      call clean_stop ! TO MODIFIED: condition is global; stop on all ranks (other ranks used to continue with unset values) (2026-10-30)
    end if

    ! Conversion factor from user units to cgs units
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    scale_msun=scale_d*scale_l**3/m_sun
    scale_pc=scale_l/3.08568d18

    allocate(Mshell(Nrbin,ntot_all),Nshell(Nrbin,ntot_all),Vshell(Nrbin,ntot_all),Mshell_before(Nrbin,ntot_all),Mshell_after(Nrbin,ntot_all))
    !collect info for newly formed stars from all cpus
    call collect_newstar_info(pos_new(1:3,1:ntot_all),mass_new(1:ntot_all))

    !define radial bins
    do ir = 1, Nrbin
      !rbin(ir) = dx_min * (1d0-1d-8) * 2**(ir-1) !put an extra factor (<1) to avoid ambiguity in defining shell regions
      rbin(ir) = dx_min * (1d0-1d-8) * f_rmax**(dble(ir-1)/(Nrbin-1)) !put an extra factor (<1) to avoid ambiguity in defining shell regions
    end do

    ! !get Mshell around each new star
    ! !call get_shell_quantities(ntot_all, pos_new(1:3,1:ntot_all), Nrbin, rbin, Mshell(1:Nrbin,1:ntot_all))
    ! call get_shell_quantities(ntot_all, pos_new(1:3,1:ntot_all), Nrbin, rbin, Mshell(1:Nrbin,1:ntot_all),&
    !      Nshell=Nshell(1:Nrbin,1:ntot_all), Vshell=Vshell(1:Nrbin,1:ntot_all))


    ! ! !KS DEBUG
    ! ! if (myid == 1) then
    ! !    print *, '(KS DEBUG) ntot_all:', ntot_all
    ! !    print *, 'Nrbin, rbin: ', Nrbin, rbin
    ! !    do istar = 1, ntot_all
    ! !       print *, 'pos_new(:): ', pos_new(1:3,istar)
    ! !       print *, 'Mshell(:): ', Mshell(1:Nrbin,istar)*scale_msun
    ! !       print *, 'Nshell(:): ', Nshell(1:Nrbin,istar)
    ! !       print *, 'Vshell(:): ', Vshell(1:Nrbin,istar)*scale_pc**3
    ! !    end do
    ! ! end if

    !depletion around each star
    do istar = 1, ntot_all

      !get Mshell around each new star
      ! note: needs to be called in do-loop to reflect on the former depletion
      ! call get_shell_quantities(ntot_all, pos_new(1:3,1:ntot_all), Nrbin, rbin, Mshell(1:Nrbin,1:ntot_all),&
      !      Nshell=Nshell(1:Nrbin,1:ntot_all), Vshell=Vshell(1:Nrbin,1:ntot_all))
      call get_shell_quantities(1, pos_new(1:3,istar), Nrbin, rbin, Mshell(1:Nrbin,istar),&
        Nshell=Nshell(1:Nrbin,istar), Vshell=Vshell(1:Nrbin,istar))


      ! !KS DEBUG
      ! if (myid == 1) then
      !    print *, '(KS DEBUG) ntot_all:', ntot_all
      !    print *, 'Nrbin, rbin: ', Nrbin, rbin
      !    do istar = 1, ntot_all
      !       print *, 'pos_new(:): ', pos_new(1:3,istar)
      !       print *, 'Mshell(:): ', Mshell(1:Nrbin,istar)*scale_msun
      !       print *, 'Nshell(:): ', Nshell(1:Nrbin,istar)
      !       print *, 'Vshell(:): ', Vshell(1:Nrbin,istar)*scale_pc**3
      !    end do
      ! end if

      !calculate Menc
      Menc(:) = 0d0
      do ir = 1, Nrbin
        if (ir == 1) then
          Menc(ir) = Mshell(1,istar)
        else
          Menc(ir) = Menc(ir-1) + Mshell(ir,istar)
        end if
      end do

      !determine size of depletion sphere (minimum radius that saisfies mass_new/M_sphere < f_dpl_max)
      do ir = 1, Nrbin
        if (Menc(ir) > mass_new(istar)/f_dpl_max) then
          r_sphere = rbin(ir)
          M_sphere = Menc(ir)
          !set f_dpl
          f_dpl = mass_new(istar) / Menc(ir)
          exit
        end if
        !最後までMencがmstar/f_dpl_maxを超えられなかった場合はエラーを吐いて止める
        if (ir == Nrbin) then
          if (myid == 1) then
            print *, '*** ERROR *** (deplete_gas_from_spheres) mass_new is too large compared with M_sphere'
            print '(A,1P2E12.4)', 'M_sphere, mass_new: ', M_sphere, mass_new(istar)
            print *, 'istar, Mshell(:,istar) for debug purpose: ', istar, Mshell(:,istar)
            print *, 'something strange may happen, or change Nrbin, f_rmax, or f_dpl_max. stopping...'
            ! stop ! TO MODIFIED: moved out of myid==1 block (2026-10-30)
          end if
          call clean_stop ! TO MODIFIED: condition is global; stop on all ranks (other ranks used to continue with unset values) (2026-10-30)
        end if
      end do

      !KS DEBUG
      if (myid == 1) then
        print '(A, I6)', '(deplete_gas_from_spheres) Nrbin: ', Nrbin
        print '(A, 2I6, 1P5E13.5)', 'istar, ir, r_sphere, M_sphere, mass_new, f_dpl:', &
          istar,  ir, r_sphere*scale_pc, M_sphere*scale_msun, mass_new(istar)*scale_msun, f_dpl
        print '(A, 1P100E9.2)', 'rbin(:) [pc]:   ', rbin(:)*scale_pc
        print '(A, 100I9)',     'Nshell(:):      ', Nshell(:,istar)
        print '(A, 1P100E9.2)', 'Menc(:) [Msun]: ', Menc(:)*scale_msun
      end if

      !------------------- KS DEBUG ----------------!
      ! ! check mass conservation before/after depletion
      ! !call get_shell_quantities(1, pos_new(1:3,istar), Nrbin, rbin, Mshell_before(1:Nrbin,istar))
      ! call get_shell_quantities(1, pos_new(1:3,1), Nrbin, rbin, Mshell_before(1:Nrbin,istar))
      ! if (myid == 1) then
      !    !calculate Menc
      !    Menc_before(:) = 0d0
      !    do ir = 1, Nrbin
      !       if (ir == 1) then
      !          Menc_before(ir) = Mshell_before(1,istar)
      !       else
      !          Menc_before(ir) = Menc_before(ir-1) + Mshell_before(ir,istar)
      !       end if
      !    end do
      ! end if
      !--------------------------------------------!


      !deplte gas by f_dpl
      !if (myid == 1) print*, '*** WARNING *** SKIP modify_sphere_by_factor'
      call modify_sphere_by_factor(pos_new(1:3,istar),r_sphere,1d0-f_dpl)  !argument is f_rho = 1 - f_dpl

      !------------------- KS DEBUG ----------------!
      ! ! check mass conservation before/after depletion
      ! !call get_shell_quantities(1, pos_new(1:3,istar), Nrbin, rbin, Mshell_after(1:Nrbin,istar))
      ! call get_shell_quantities(1, pos_new(1:3,1), Nrbin, rbin, Mshell_after(1:Nrbin,istar))
      ! if (myid == 1) then
      !    !calculate Menc
      !    Menc_after(:) = 0d0
      !    do ir = 1, Nrbin
      !       if (ir == 1) then
      !          Menc_after(ir) = Mshell_after(1,istar)
      !       else
      !          Menc_after(ir) = Menc_after(ir-1) + Mshell_after(ir,istar)
      !       end if
      !    end do
      !    !print '(A,1P1E10.2,1P3E20.12,/,A,1P4E20.12,/)', '* CHECK MASS CONSERVATION * Rcheck, pos: ',&
      !         !rbin(Nrbin)*scale_pc,pos_new(1:3,istar), &
      !    print '(A,1P1E10.2,1P3E20.12,/,A,1P4E20.12,/)', '* CHECK MASS CONSERVATION * Rcheck, pos1: ',&
      !         rbin(Nrbin)*scale_pc,pos_new(1:3,1), &
      !         'Menc_before, Menc_after, mass_new, err: ',&
      !          Menc_before(Nrbin)*scale_msun, Menc_after(Nrbin)*scale_msun,&
      !          mass_new(istar)*scale_msun, (Menc_after(Nrbin)+mass_new(istar))/Menc_before(Nrbin)-1d0
      ! end if
      !--------------------------------------------!


    end do

  end subroutine deplete_gas_from_spheres

  !----------------------------------------------------------------------------
  ! non-locally deplete gas with fixed post-star-forming cell density
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! - set the density of star-forming region (currently the cell where the star is forming) to a given value
  ! - ensure mass conservation by adjusting density in outer shells by a constant
  !   factor f_rho that satisfies |f_rho - 1| < f_mod_max
  !----------------------------------------------------------------------------
  subroutine deplete_gas_with_fixed_postSF_density

    !parameters
    real(dp),parameter :: f_mod_max = 0.2d0 !maximum change of density in surrounding cells
    integer,parameter :: Nrbin = 20 ! number of radial bins for shells (r_i = r_0 * (r_max/r_0)**(i/(Nrbin-1)))
    real(dp),parameter :: f_rmax = 1d2 !rmax/rmin: the ratio of minmum and maximum radii of radial bins (rmin is set to dx_min)

    !variables
    real(dp),dimension(1:Nnew_max) :: mass_new     !mass of newly formed stars
    real(dp),dimension(1:3,1:Nnew_max) :: pos_new  !position of newly formed stars
    real(dp),dimension(1:Nrbin) :: rbin      !outer radii of shells
    real(dp),allocatable,dimension(:,:) :: Mshell    !mass of shells
    integer,allocatable,dimension(:,:) :: Nshell    !number of cell in shells
    real(dp),allocatable,dimension(:,:) :: Vshell    !volume of shells
    real(dp),dimension(1:Nrbin) :: Menc      !enclosed mass
    real(dp) :: M_sphere, r_sphere
    real(dp) :: f_mod, f_rho ! reduction/amplification factor: f_rho = f_mod + 1
    real(dp) ::  M_mod, M_fix, R_in, rho_after, rho_before

    real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v, scale_msun, scale_pc
    integer :: ir, istar

    !for debug
    real(dp),allocatable,dimension(:,:) :: Mshell_before, Mshell_after    !mass of shells before/after depletion
    real(dp),dimension(1:Nrbin) :: Menc_before, Menc_after      !enclosed mass before/after depletion

    !Note:
    ! ntot_star_all(icpu) is # of new stars in icpu
    ! ntot_all is sum of ntot_star_all(:)

    !skip if there is no new stars
    if (ntot_all == 0) return

    !check number of newly formed stars
    if (ntot_all > Nnew_max) then
      if (myid == 1) then
        print '(A,2I8)', 'number of newly formed stars exceeds the maximum value',ntot_all, Nnew_max
        print '(A)', 'stopping...'
        ! stop ! TO MODIFIED: moved out of myid==1 block (2026-10-30)
      end if
      call clean_stop ! TO MODIFIED: condition is global; stop on all ranks (other ranks used to continue with unset values) (2026-10-30)
    end if

    ! Conversion factor from user units to cgs units
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    scale_msun=scale_d*scale_l**3/m_sun
    scale_pc=scale_l/3.08568d18


    allocate(Mshell(Nrbin,ntot_all),Nshell(Nrbin,ntot_all),Vshell(Nrbin,ntot_all),Mshell_before(Nrbin,ntot_all),Mshell_after(Nrbin,ntot_all))
    !collect info for newly formed stars from all cpus
    call collect_newstar_info(pos_new(1:3,1:ntot_all),mass_new(1:ntot_all))

    !define radial bins
    do ir = 1, Nrbin
      !rbin(ir) = dx_min * (1d0-1d-8) * 2**(ir-1) !put an extra factor (<1) to avoid ambiguity in defining shell regions
      rbin(ir) = dx_min * (1d0-1d-8) * f_rmax**(dble(ir-1)/(Nrbin-1)) !put an extra factor (<1) to avoid ambiguity in defining shell regions
    end do

    !determine radius and f_dpl
    do istar = 1, ntot_all

      !get Mshell around each new star
      ! note: needs to be called in do-loop to reflect on the former depletion
      ! call get_shell_quantities(ntot_all, pos_new(1:3,1:ntot_all), Nrbin, rbin, Mshell(1:Nrbin,1:ntot_all),&
      !      Nshell=Nshell(1:Nrbin,1:ntot_all), Vshell=Vshell(1:Nrbin,1:ntot_all))
      call get_shell_quantities(1, pos_new(1:3,istar), Nrbin, rbin, Mshell(1:Nrbin,istar),&
        Nshell=Nshell(1:Nrbin,istar), Vshell=Vshell(1:Nrbin,istar))

      ! !KS DEBUG
      ! if (myid == 1) then
      !    print *, '(KS DEBUG) ntot_all:', ntot_all
      !    print *, 'Nrbin, rbin: ', Nrbin, rbin
      !    do istar = 1, ntot_all
      !       print *, 'pos_new(:): ', pos_new(1:3,istar)
      !       print *, 'Mshell(:): ', Mshell(1:Nrbin,istar)*scale_msun
      !       print *, 'Nshell(:): ', Nshell(1:Nrbin,istar)
      !       print *, 'Vshell(:): ', Vshell(1:Nrbin,istar)*scale_pc**3
      !    end do
      ! end if

      !calculate Menc
      Menc(:) = 0d0
      do ir = 1, Nrbin
        if (ir == 1) then
          Menc(ir) = 0   !the density of the inner-most shell will be set to a given value
        else
          Menc(ir) = Menc(ir-1) + Mshell(ir,istar) !density of outer shells will be adjusted to ensure mass conservation
        end if
      end do

      !radius of the inner-most shell
      R_in = rbin(1)

      !density before and after modification
      rho_before = Mshell(1,istar)/Vshell(1,istar) !density before modification
      rho_after = nH_postSF/scale_nH           !set to a fixed value

      !mass modificaition corresponding to the
      M_mod = (rho_before - rho_after) * Vshell(1,istar)  - mass_new(istar)

      !mass to be added to the surrounding region (consider new star's mass and density modification)
      M_mod = (rho_before - rho_after) * Vshell(1,istar)  - mass_new(istar)
      M_fix = (rho_before - rho_after) * Vshell(1,istar) !mass difference due to density modification

      !determine size of sphere where mass is added/removed
      ! (minimum radius that saisfies |f_mod| = |M_mod/Menc(ir)| < f_mod_max)
      do ir = 2, Nrbin
        ! f_mod = M_mod/Menc(ir)   !ゼロ割りが発生することがあったので対応
        ! if (abs(f_mod) < f_mod_max) then
        if (abs(M_mod) < f_mod_max*Menc(ir)) then
          r_sphere = rbin(ir)
          M_sphere = Menc(ir)
          f_mod = M_mod/Menc(ir)
          f_rho = 1 + f_mod
          exit
        end if
        !最後まで半径が見つからなかった場合はエラーを吐いて止める
        if (ir == Nrbin) then
          if (myid == 1) then
            print *, '*** ERROR *** (deplete_gas_with_fixed_postSF_density) mass_new is too large compared with M_sphere'
            print '(A,1P4E12.4)', 'Menc(Nrbin), M_mod, M_fix, mass_new: ', Menc(Nrbin), M_mod, M_fix, mass_new(istar)
            print *, 'something strange may happen, or change Nrbin, f_rmax, or f_mod_max. stopping...'
            ! stop ! TO MODIFIED: moved out of myid==1 block (2026-10-30)
          end if
          call clean_stop ! TO MODIFIED: condition is global; stop on all ranks (other ranks used to continue with unset values) (2026-10-30)
        end if
      end do

      !KS DEBUG
      if (myid == 1) then
        print '(A, I4, 1P3E13.5)', '(deplete_gas_with_fixed_postSF_density) Nrbin, rmin, rmax, nH_postSF: ', &
          Nrbin, rbin(1)*scale_pc, rbin(Nrbin)*scale_pc, nH_postSF
        print '(A, I4, 1P5E12.4)', 'istar, dx_min, nH_before, nH_after, M_fix, M_star:', &
          istar,  dx_min*scale_pc, rho_before*scale_nH, rho_after*scale_nH, M_fix*scale_msun, mass_new(istar)*scale_msun
        print '(A, I4, 1P4E12.4)', 'ir, r_sphere, M_sphere, M_mod, f_rho', &
          ir, r_sphere*scale_pc, M_sphere*scale_msun, M_mod*scale_msun, f_rho
        print '(A, 1P100E9.2)', 'rbin(:) [pc]:   ', rbin(:)*scale_pc
        print '(A, 100I9)',     'Nshell(:):      ', Nshell(:,istar)
        print '(A, 1P100E9.2)', 'Menc(:) [Msun]: ', Menc(:)*scale_msun
      end if

      !------------------- KS DEBUG ----------------!
      ! check mass conservation before/after depletion
      !call get_shell_quantities(1, pos_new(1:3,istar), Nrbin, rbin, Mshell_before(1:Nrbin,istar))
      call get_shell_quantities(1, pos_new(1:3,1), Nrbin, rbin, Mshell_before(1:Nrbin,istar))
      if (myid == 1) then
        !calculate Menc
        Menc_before(:) = 0d0
        do ir = 1, Nrbin
          if (ir == 1) then
            Menc_before(ir) = Mshell_before(1,istar)
          else
            Menc_before(ir) = Menc_before(ir-1) + Mshell_before(ir,istar)
          end if
        end do
      end if
      !--------------------------------------------!

      !modify post-star-forming gas density to rho_after (for r < R_in)
      call modify_sphere_to_costant_density(pos_new(1:3,istar), R_in, rho_after)

      !modify surrounding gas density by f_rho (for R_in < r < r_sphere)
      call modify_sphere_by_factor(pos_new(1:3,istar), r_sphere, f_rho, rmin = R_in)

      !------------------- KS DEBUG ----------------!
      ! check mass conservation before/after depletion
      !call get_shell_quantities(1, pos_new(1:3,istar), Nrbin, rbin, Mshell_after(1:Nrbin,istar))
      call get_shell_quantities(1, pos_new(1:3,1), Nrbin, rbin, Mshell_after(1:Nrbin,istar))
      if (myid == 1) then
        !calculate Menc
        Menc_after(:) = 0d0
        do ir = 1, Nrbin
          if (ir == 1) then
            Menc_after(ir) = Mshell_after(1,istar)
          else
            Menc_after(ir) = Menc_after(ir-1) + Mshell_after(ir,istar)
          end if
        end do
        !print '(A,1P1E10.2,1P3E20.12,/,A,1P4E20.12,/)', '* CHECK MASS CONSERVATION * Rcheck, pos: ',&
        !rbin(Nrbin)*scale_pc,pos_new(1:3,istar), &
        print '(A,1P1E10.2,1P3E20.12,/,A,1P4E20.12,/)', '* CHECK MASS CONSERVATION * Rcheck, pos1: ',&
          rbin(Nrbin)*scale_pc,pos_new(1:3,1), &
          'Menc_before, Menc_after, mass_new, err: ',&
          Menc_before(Nrbin)*scale_msun, Menc_after(Nrbin)*scale_msun,&
          mass_new(istar)*scale_msun, (Menc_after(Nrbin)+mass_new(istar))/Menc_before(Nrbin)-1d0
      end if
      !--------------------------------------------!


    end do

  end subroutine deplete_gas_with_fixed_postSF_density



  !----------------------------------------------------------------------------
  ! deplete_gas_from_spheres_batch (merge_ks 2026-11-13): non-local gas depletion of all
  ! new stars of this step at once (nonlocal_depletion_batch=.true., the default).
  !
  ! Same rule as deplete_gas_from_spheres (ks): around every new star the smallest sphere
  ! (radii rbin, geometric from dx_min to f_rmax*dx_min) with M_enc > M_star/f_dpl_max is
  ! depleted uniformly by f = M_star/M_enc. Differences from the sequential ks version:
  !  - the shell masses of all stars are accumulated locally by walking the octree down to
  !    the leaf cells of my active grids inside each sphere and summed with ONE ALLREDUCE
  !    (ks: one sweep over all cells and several ALLREDUCEs per star);
  !  - the depletion is additive, rho -> rho*(1 - sum_i f_i) over the spheres containing
  !    the cell, so the removed mass is exactly sum_i M_i; a cell shared by several spheres
  !    is protected by sum_i f_i <= f_sum_max: when exceeded, the radii of the stars
  !    involved are enlarged (f_i decreases) and the check is repeated (a few iterations,
  !    one ALLREDUCE of n integers each). Without overlap the result is that of the ks
  !    version up to roundoff.
  ! Collective calls: 1 (star list) + 1-3 (shells) + 1-3 (overlap) per step, independent of
  ! the number of stars.
  !----------------------------------------------------------------------------
  subroutine deplete_gas_from_spheres_batch
    real(dp),parameter :: f_dpl_max = 0.2d0 ! depletion radius: minimum radius with M_star/M_enc < f_dpl_max
    real(dp),parameter :: f_sum_max = 0.5d0 ! maximum summed depletion of a cell shared by several spheres
    integer,parameter :: Nrbin = 40         ! radial bins (as deplete_gas_from_spheres)
    real(dp),parameter :: f_rmax = 1d2      ! rmax/rmin (rmin = dx_min)
    integer,parameter :: ir_first = 13      ! shells computed first (rbin(13) ~ 4 dx_min), doubled when insufficient
    integer,parameter :: niter_max = 20     ! overlap iterations

    real(dp),dimension(1:Nrbin) :: rbin
    real(dp),allocatable :: mass_new(:), pos_new(:,:)
    real(dp),allocatable :: Mshell_loc(:,:), Mshell(:,:), Mshell_tmp(:,:)
    real(dp),allocatable :: f_dpl(:), Menc(:)
    integer,allocatable :: ir_need(:), ir_sph(:), flag_loc(:), flag_glob(:)
    integer :: n, k, ir, iter, p, ic, il, nshell_iter, nflag, ncell_tot
    logical :: all_ok
    real(dp) :: fac, m_removed_loc, m_removed, m_target, sumf_max_loc, sumf_max, m_clamped_loc, m_clamped, vol_c
    real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_msun,scale_pc

    if (ntot_all == 0) return
    n = ntot_all

    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    scale_msun=scale_d*scale_l**3/m_sun
    scale_pc=scale_l/3.08568d18

    ! per-cell work array (zero outside the loops below)
    if (.not. allocated(dep_sumf)) then
      allocate(dep_sumf(0:ncoarse+twotondim*ngridmax))
      dep_sumf = 0d0
    end if

    ! positions and masses of all new stars (one cell with n stars counts as one star of mass n*mstar)
    allocate(mass_new(n), pos_new(3,n))
    call collect_newstar_info(pos_new, mass_new)

    ! radial bins
    do ir = 1, Nrbin
      rbin(ir) = dx_min * (1d0-1d-8) * f_rmax**(dble(ir-1)/(Nrbin-1))
    end do

    !------------------------------------------------------------------
    ! shell masses -> depletion radius and factor of every star
    !------------------------------------------------------------------
    allocate(Mshell_loc(Nrbin,n), Mshell(Nrbin,n), Mshell_tmp(Nrbin,n), Menc(Nrbin), f_dpl(n), ir_need(n), ir_sph(n))
    ir_need(:) = ir_first; ir_sph(:) = 0; f_dpl(:) = 0d0
    nshell_iter = 0
    do
      nshell_iter = nshell_iter + 1
      Mshell_loc(:,:) = 0d0
      do k = 1, n
        if (ir_sph(k) > 0) cycle
        call sphere_walk(pos_new(1:3,k), rbin(ir_need(k)), k, 1, Nrbin, rbin, Mshell_loc(1:Nrbin,k))
      end do
      call MPI_ALLREDUCE(Mshell_loc, Mshell, Nrbin*n, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
      all_ok = .true.
      do k = 1, n
        if (ir_sph(k) > 0) cycle
        Menc(1) = Mshell(1,k)
        do ir = 2, ir_need(k)
          Menc(ir) = Menc(ir-1) + Mshell(ir,k)
        end do
        do ir = 1, ir_need(k)
          if (Menc(ir) > mass_new(k)/f_dpl_max) then
            ir_sph(k) = ir
            f_dpl(k) = mass_new(k)/Menc(ir)
            exit
          end if
        end do
        if (ir_sph(k) == 0) then
          if (ir_need(k) == Nrbin) then
            if (myid == 1) then
              print *, '*** ERROR *** (deplete_gas_from_spheres_batch) mass_new is too large compared with M_sphere'
              print '(A,I8,1P2E12.4)', 'istar, Menc(Nrbin), mass_new: ', k, Menc(Nrbin), mass_new(k)
              print *, 'something strange may happen, or change Nrbin, f_rmax, or f_dpl_max. stopping...'
            end if
            call clean_stop
          end if
          ir_need(k) = min(Nrbin, 2*ir_need(k))
          all_ok = .false.
        end if
      end do
      if (all_ok) exit
    end do

    !------------------------------------------------------------------
    ! overlap check: enlarge the spheres whose cells would be depleted by more than f_sum_max
    !------------------------------------------------------------------
    allocate(flag_loc(n), flag_glob(n))
    npair_max = 64*n
    allocate(pr_star(npair_max), pr_cell(npair_max), pr_lev(npair_max))
    do iter = 1, niter_max
      ! (star, cell) pairs of my leaf cells inside the spheres
      npair = 0
      do k = 1, n
        call sphere_walk(pos_new(1:3,k), rbin(ir_sph(k)), k, 2, Nrbin, rbin, Mshell_loc(1:Nrbin,k))
      end do
      do p = 1, npair
        dep_sumf(pr_cell(p)) = dep_sumf(pr_cell(p)) + f_dpl(pr_star(p))
      end do
      flag_loc(:) = 0
      sumf_max_loc = 0d0
      do p = 1, npair
        sumf_max_loc = max(sumf_max_loc, dep_sumf(pr_cell(p)))
        if (dep_sumf(pr_cell(p)) > f_sum_max .and. ir_sph(pr_star(p)) < Nrbin) flag_loc(pr_star(p)) = 1
      end do
      call MPI_ALLREDUCE(flag_loc, flag_glob, n, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, info)
      nflag = sum(flag_glob)
      if (nflag == 0) exit
      ! reset the sums of the touched cells and enlarge the flagged spheres by one bin
      do p = 1, npair
        dep_sumf(pr_cell(p)) = 0d0
      end do
      all_ok = .true.
      flag_loc(:) = 0                         ! reused: 1 = shells must be recomputed to a larger radius
      do k = 1, n
        if (flag_glob(k) == 0) cycle
        ir_sph(k) = ir_sph(k) + 1
        if (ir_sph(k) > ir_need(k)) then      ! shells of this star not yet computed up to the new radius
          ir_need(k) = min(Nrbin, 2*ir_need(k))
          flag_loc(k) = 1
          all_ok = .false.
        end if
      end do
      if (.not. all_ok) then                  ! (all ranks agree: flag_glob is global)
        Mshell_loc(:,:) = 0d0
        do k = 1, n
          if (flag_loc(k) == 0) cycle
          call sphere_walk(pos_new(1:3,k), rbin(ir_need(k)), k, 1, Nrbin, rbin, Mshell_loc(1:Nrbin,k))
        end do
        call MPI_ALLREDUCE(Mshell_loc, Mshell_tmp, Nrbin*n, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
        do k = 1, n
          if (flag_loc(k) == 1) Mshell(1:Nrbin,k) = Mshell_tmp(1:Nrbin,k)
        end do
      end if
      do k = 1, n
        if (flag_glob(k) == 0) cycle
        Menc(1) = Mshell(1,k)
        do ir = 2, ir_sph(k)
          Menc(ir) = Menc(ir-1) + Mshell(ir,k)
        end do
        f_dpl(k) = mass_new(k)/Menc(ir_sph(k))
      end do
    end do

    !------------------------------------------------------------------
    ! apply: rho -> rho*(1 - sum f) once per cell (ufront, uold and unew as modify_sphere_by_factor)
    !------------------------------------------------------------------
    m_removed_loc = 0d0; m_clamped_loc = 0d0; ncell_tot = 0
    do p = 1, npair
      ic = pr_cell(p); il = pr_lev(p)
      if (dep_sumf(ic) <= 0d0) cycle          ! already applied (cell shared by several stars)
      vol_c = (0.5d0**il*scale)**3
      if (dep_sumf(ic) > f_sum_max) then      ! could not be resolved by enlarging (radius limit / iterations)
        m_clamped_loc = m_clamped_loc + (dep_sumf(ic)-f_sum_max)*ufront(ic,1)*vol_c
        dep_sumf(ic) = f_sum_max
      end if
      fac = 1d0 - dep_sumf(ic)
      m_removed_loc = m_removed_loc + dep_sumf(ic)*ufront(ic,1)*vol_c
      ncell_tot = ncell_tot + 1
      do ivar = 1, nvar
        ufront(ic,ivar) = fac*ufront(ic,ivar)
      end do
      if (il == nlevelmax) then               ! uold is primitive here: the density only
        uold(ic,1) = ufront(ic,1)
      else                                    ! conservative uold and unew (difference of the update)
        do ivar = 1, nvar
          uold(ic,ivar) = ufront(ic,ivar)
          unew(ic,ivar) = unew(ic,ivar) - (1d0-fac)/fac*ufront(ic,ivar)
        end do
      end if
      dep_sumf(ic) = 0d0
    end do

    ! report
    call MPI_ALLREDUCE(m_removed_loc, m_removed, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
    call MPI_ALLREDUCE(m_clamped_loc, m_clamped, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
    call MPI_ALLREDUCE(sumf_max_loc, sumf_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, info)
    m_target = sum(mass_new(1:n))
    if (myid == 1) then
      print '(A,I8,A,I3,A,I3,A,1P2E12.4,A,0PF6.3,A,I3)', &
        '(deplete_gas_from_spheres_batch) stars: ', n, ', shell rounds: ', nshell_iter, ', overlap iterations: ', min(iter,niter_max), &
        ', M_star, M_removed [Msun]: ', m_target*scale_msun, m_removed*scale_msun, ', max sum f: ', sumf_max, ', max ir: ', maxval(ir_sph)
      if (m_clamped > 0d0) print '(A,1PE12.4,A)', &
        '*** WARNING *** (deplete_gas_from_spheres_batch) summed depletion clamped at f_sum_max: ', m_clamped*scale_msun, &
        ' Msun of stellar mass was not removed from the gas'
    end if

    deallocate(mass_new, pos_new, Mshell_loc, Mshell, Mshell_tmp, Menc, f_dpl, ir_need, ir_sph, flag_loc, flag_glob, pr_star, pr_cell, pr_lev)

  end subroutine deplete_gas_from_spheres_batch

  !--------------------------------------------------------------------------
  ! octree walk over the leaf cells of my active grids inside the sphere (pos, r):
  ! mode 1: shell masses msh(ir) += vol*rho (ufront), mode 2: append (k, cell, level) pairs
  !--------------------------------------------------------------------------
  subroutine sphere_walk(pos, r, k, mode, Nr, rb, msh)
    real(dp),dimension(1:3),intent(in) :: pos
    real(dp),intent(in) :: r
    integer,intent(in) :: k, mode, Nr
    real(dp),dimension(1:Nr),intent(in) :: rb
    real(dp),dimension(1:Nr),intent(inout) :: msh
    integer :: icoarse_cell, ig
    do icoarse_cell = 1, ncoarse
      ig = son(icoarse_cell)
      if (ig > 0) call sphere_walk_grid(ig, 1, pos, r, k, mode, Nr, rb, msh)
    end do
  end subroutine sphere_walk

  recursive subroutine sphere_walk_grid(ig, il, pos, r, k, mode, Nr, rb, msh)
    integer,intent(in) :: ig, il, k, mode, Nr
    real(dp),dimension(1:3),intent(in) :: pos
    real(dp),intent(in) :: r
    real(dp),dimension(1:Nr),intent(in) :: rb
    real(dp),dimension(1:Nr),intent(inout) :: msh
    integer :: ind, idim, ic, ir, ix, iy, iz
    real(dp) :: h, hs, d2, rc, xc_c, dd
    logical :: mine
    real(dp),dimension(1:3) :: cen
    h = 0.5d0**(il+1)                 ! half cell size at level il (grid units)
    hs = h*scale                      ! in code units
    mine = (cpu_map(father(ig)) == myid)
    do ind = 1, twotondim
      iz = (ind-1)/4; iy = (ind-1-4*iz)/2; ix = (ind-1-2*iy-4*iz)
      cen(1) = (xg(ig,1) + (dble(ix)-0.5d0)*2d0*h - skip_loc(1))*scale
      cen(2) = (xg(ig,2) + (dble(iy)-0.5d0)*2d0*h - skip_loc(2))*scale
      cen(3) = (xg(ig,3) + (dble(iz)-0.5d0)*2d0*h - skip_loc(3))*scale
      ! does the cube intersect the sphere?
      d2 = 0d0
      do idim = 1, 3
        dd = max(abs(pos(idim)-cen(idim)) - hs, 0d0)
        d2 = d2 + dd*dd
      end do
      if (d2 >= r*r) cycle
      ic = ncoarse + (ind-1)*ngridmax + ig
      if (son(ic) > 0) then
        call sphere_walk_grid(son(ic), il+1, pos, r, k, mode, Nr, rb, msh)
      else if (mine) then
        rc = sqrt((pos(1)-cen(1))**2 + (pos(2)-cen(2))**2 + (pos(3)-cen(3))**2)
        if (rc >= r) cycle
        if (mode == 1) then
          do ir = 1, Nr
            if (rc <= rb(ir)) exit
          end do
          if (ir > Nr) cycle
          msh(ir) = msh(ir) + (2d0*hs)**3*ufront(ic,1)
        else
          if (npair == npair_max) call grow_pairs
          npair = npair + 1
          pr_star(npair) = k; pr_cell(npair) = ic; pr_lev(npair) = il
        end if
      end if
    end do
  end subroutine sphere_walk_grid

  subroutine grow_pairs
    integer,allocatable :: t1(:), t2(:), t3(:)
    allocate(t1(2*npair_max), t2(2*npair_max), t3(2*npair_max))
    t1(1:npair) = pr_star(1:npair); t2(1:npair) = pr_cell(1:npair); t3(1:npair) = pr_lev(1:npair)
    call move_alloc(t1, pr_star); call move_alloc(t2, pr_cell); call move_alloc(t3, pr_lev)
    npair_max = 2*npair_max
  end subroutine grow_pairs


  !----------------------------------------------------------------------------
  ! collect info for newly formed stars from all cpus
  !----------------------------------------------------------------------------
  subroutine collect_newstar_info(pos_new,mass_new)
    real(dp), intent(OUT) :: pos_new(:,:), mass_new(:)

    integer :: istar_loc, istar_glob, istar_offset
    real(dp),dimension(1:Nnew_max) :: x_loc, y_loc, z_loc, mass_loc, x_glob, y_glob, z_glob, mass_glob

    !MPI send/recv buffer
    integer,parameter :: N_realdata=4 !number of data in real/integer for MPI communication
    real(dp),dimension(1:Nnew_max*N_realdata) :: sbufr,rbufr

    !----------------------------------------
    ! collect info about newly formed stars
    !----------------------------------------
    !offset in each cpu with respect to global index
    istar_offset = 0
    do icpu = 1, myid-1  ! "do i = 1,0" loop is just skipped (KS MEMO)
      istar_offset = istar_offset + ntot_star_all(icpu)
    end do

    ! collect cell info in each cpu
    istar_loc = 0
    ncache=active(ilevel)%ngrid
    do igrid=1,ncache,nvector
      ngrid=MIN(nvector,ncache-igrid+1)
      do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
      end do
      ! Loop over cells
      do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
          ind_cell(i)=iskip+ind_grid(i)
        end do

        !store information about newly formed stars and their parent cells
        do i=1,ngrid
          if (flag2(ind_cell(i)).gt.0 .or. flag3(ind_cell(i)).gt.0) then !Pop II/Pop III
            istar_loc=istar_loc+1
            x_loc(istar_loc) = (xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale          !x
            y_loc(istar_loc) = (xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale          !y
            z_loc(istar_loc) = (xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale          !z
            !type
            if (flag2(ind_cell(i)).gt.0)then !Pop II
              mass_loc(istar_loc) = flag2(ind_cell(i))*mstar                           !Mass
            else if (flag3(ind_cell(i)).gt.0)then ! Pop III
              mass_loc(istar_loc) = flag3(ind_cell(i))                                 !Mass
            end if
          end if
        end do
      end do
      ! End loop over cells
    end do
    ! End loop over grids


    !----------------------------------------
    ! communiation for creating global data
    !----------------------------------------
    !create MPI sending buffer at each cpu
    sbufr(:)=0d0
    do istar_loc = 1, ntot_star_all(myid)
      istar_glob = istar_loc + istar_offset       !relation between local and global index

      !real quantities (N_realdata)
      sbufr((istar_glob-1)*N_realdata+1) = x_loc(istar_loc)
      sbufr((istar_glob-1)*N_realdata+2) = y_loc(istar_loc)
      sbufr((istar_glob-1)*N_realdata+3) = z_loc(istar_loc)
      sbufr((istar_glob-1)*N_realdata+4) = mass_loc(istar_loc)
    end do

    !MPI communication
    call MPI_ALLREDUCE(sbufr, rbufr, ntot_all*N_realdata, MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)

    !unload received data
    do istar_glob = 1, ntot_all
      !real quantities (N_realdata)
      x_glob(istar_glob) = rbufr((istar_glob-1)*N_realdata+1)
      y_glob(istar_glob) = rbufr((istar_glob-1)*N_realdata+2)
      z_glob(istar_glob) = rbufr((istar_glob-1)*N_realdata+3)
      mass_glob(istar_glob) = rbufr((istar_glob-1)*N_realdata+4)
    end do

    !return data
    mass_new(1:ntot_all) = mass_glob(1:ntot_all)
    pos_new(1,1:ntot_all) = x_glob(1:ntot_all)
    pos_new(2,1:ntot_all) = y_glob(1:ntot_all)
    pos_new(3,1:ntot_all) = z_glob(1:ntot_all)

    ! !KS DEBUG
    ! print *, 'mass_new(1:ntot_all): ', mass_new(1:ntot_all)
    ! print *, 'pos_new(1:3,1:ntot_all): ', pos_new(1:3,1:ntot_all)

  end subroutine collect_newstar_info

  !----------------------------------------------------------------------------
  ! get shell quantities
  !----------------------------------------------------------------------------
  subroutine get_shell_quantities(Npos, pos, Nrbin, rbin, Mshell, Nshell, Vshell, MZshell, ETshell, MomRshell, RawMomRshell)

    !in/out
    integer,intent(in) :: Npos !number of centers
    real(dp),dimension(1:3,1:Npos),intent(in) :: pos  !position of centers
    integer,intent(in) :: Nrbin !number of radial bins
    real(dp),dimension(1:Nrbin),intent(in) :: rbin  !outer radius of shells
    real(dp),dimension(1:Nrbin,1:Npos),intent(out)::Mshell !mass in each shell
    integer,dimension(1:Nrbin,1:Npos),intent(out),optional::Nshell !Number of cell in each shell
    real(dp),dimension(1:Nrbin,1:Npos),intent(out),optional::Vshell !volume in each shell
    real(dp),dimension(1:Nrbin,1:Npos),intent(out),optional::MZshell !metal mass in each shell
    real(dp),dimension(1:Nrbin,1:Npos),intent(out),optional::ETshell !thermal energy in each shell
    real(dp),dimension(1:Nrbin,1:Npos),intent(out),optional::RawMomRshell !radial momentum in each shell (Raw value)
    real(dp),dimension(1:Nrbin,1:Npos),intent(out),optional::MomRshell !radial momentum in each shell (vbulk subtracted)

    !real(dp),intent(out),optional::Mshell(:) !mass of each shell (optional)

    !variables
    integer :: sweep_level !search for all leaf cells through levels (not only ilevel)
    real(dp) :: vol_cell, dx_cell
    real(dp),dimension(1:Nrbin,1:Npos) :: Mshell_loc, Mshell_glob
    integer,dimension(1:Nrbin,1:Npos) :: Nshell_loc, Nshell_glob
    real(dp),dimension(1:Nrbin,1:Npos) :: Vshell_loc, Vshell_glob
    real(dp),dimension(1:Nrbin,1:Npos) :: MZshell_loc, MZshell_glob
    real(dp),dimension(1:Nrbin,1:Npos) :: ETshell_loc, ETshell_glob
    real(dp),dimension(1:Nrbin,1:Npos) :: RawMomRshell_loc, RawMomRshell_glob
    real(dp),dimension(1:Nrbin,1:Npos) :: MomXshell_loc, MomXshell_glob, MomYshell_loc, MomYshell_glob, MomZshell_loc, MomZshell_glob
    real(dp),dimension(1:Nrbin,1:Npos) :: MomRshell_loc, MomRshell_glob
    real(dp) :: x, y, z, r2, r
    integer :: ipos, ir, ix, iy, iz
    real(dp),dimension(1:twotondim,1:3)::xc_level
    real(dp) :: MomR_cell
    real(dp) :: erx,ery,erz !unit radial vector
    real(dp), parameter:: yr = 3.15569d+07
    real(dp) :: vbulk_x,vbulk_y,vbulk_z,rho

    !----------------------------------------
    ! sum variables in each cpu
    !----------------------------------------
    Mshell_loc(:,:) = 0d0
    Nshell_loc(:,:) = 0
    Vshell_loc(:,:) = 0d0
    MZshell_loc(:,:) = 0d0
    ETshell_loc(:,:) = 0d0
    RawMomRshell_loc(:,:) = 0d0
    MomRshell_loc(:,:) = 0d0
    MomXshell_loc(:,:) = 0d0
    MomYshell_loc(:,:) = 0d0
    MomZshell_loc(:,:) = 0d0
    do sweep_level = levelmin,nlevelmax
      dx_cell = 1d0/dble(2**sweep_level)
      vol_cell = (dx_cell*scale)**3   ! merge_ks: cell volume in code units for boxlen /= 1 (ks assumed boxlen=1) !cell volume at sweep_level

      ! Cells center position relative to grid center position
      do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc_level(ind,1)=(dble(ix)-0.5D0)*dx_cell
        xc_level(ind,2)=(dble(iy)-0.5D0)*dx_cell
        xc_level(ind,3)=(dble(iz)-0.5D0)*dx_cell
      end do


      ! sum shell qunatities in each cpu
      ncache=active(sweep_level)%ngrid
      do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
          ind_grid(i)=active(sweep_level)%igrid(igrid+i-1)
        end do
        ! Loop over cells
        do ind=1,twotondim
          iskip=ncoarse+(ind-1)*ngridmax
          do i=1,ngrid
            ind_cell(i)=iskip+ind_grid(i)
          end do

          !store information about newly formed stars and their parent cells
          do i=1,ngrid

            !sum over only leaf cells
            if(son(ind_cell(i))/=0) cycle

            x = (xg(ind_grid(i),1)+xc_level(ind,1)-skip_loc(1))*scale          !x
            y = (xg(ind_grid(i),2)+xc_level(ind,2)-skip_loc(2))*scale          !y
            z = (xg(ind_grid(i),3)+xc_level(ind,3)-skip_loc(3))*scale          !z

            do ipos = 1, Npos
              r2 = (x - pos(1,ipos))**2 +  (y - pos(2,ipos))**2 +  (z - pos(3,ipos))**2
              r = sqrt(r2)


              do ir = 1, Nrbin
                !skip unless the cell is in the shell
                if (ir == 1) then !0 < r <= rbin(1) for initial shell
                  if (r > rbin(1)) cycle
                else                !rbin(i-1) < r <= rbin(i) otherwise
                  if (r <= rbin(ir-1) .or. r > rbin(ir)) cycle
                end if

                ! !KS DEBUG
                ! print *, 'sweep_level, ir, rbin, r, x, y, z: ', sweep_level, ir, rbin(ir), r, x, y, z

                !------------------
                ! sum variables
                !------------------
                Mshell_loc(ir,ipos) = Mshell_loc(ir,ipos) + vol_cell * ufront(ind_cell(i),1)
                !Mshell_loc(ir,ipos) = Mshell_loc(ir,ipos) + vol_cell * uold(ind_cell(i),1)
                ! if (sweep_level == nlevelmax) then !use uold for level = nlevelmax and unew otherwise
                !    Mshell_loc(ir,ipos) = Mshell_loc(ir,ipos) + vol_cell * uold(ind_cell(i),1)
                ! else
                !    Mshell_loc(ir,ipos) = Mshell_loc(ir,ipos) + vol_cell * unew(ind_cell(i),1)
                ! end if

                if (present(Nshell)) Nshell_loc(ir,ipos) = Nshell_loc(ir,ipos) + 1
                if (present(Vshell)) Vshell_loc(ir,ipos) = Vshell_loc(ir,ipos) + vol_cell
                if (present(MZshell)) then
                  !note: metallicity variable of uold, divided by density at ilevel (= nlevelmax), is not used
                  MZshell_loc(ir,ipos) = MZshell_loc(ir,ipos) + vol_cell *  ufront(ind_cell(i),imetal)
                  ! ! uold for metallicity was devided by density in advance only in ilevel (= nlevelmax, in the current code)
                  ! if (sweep_level == nlevelmax) then
                  !    MZshell_loc(ir,ipos) = MZshell_loc(ir,ipos) + vol_cell * uold(ind_cell(i),1) * uold(ind_cell(i),imetal)
                  !    ! print '(A,2I4, 1P2E13.3)', 'sweep_level, ir, MZshell_loc, uold(imetal), unew(imetal)', &
                  !    !      sweep_level, ir, uold(ind_cell(i),imetal), unew(ind_cell(i),imetal)
                  ! else
                  !    MZshell_loc(ir,ipos) = MZshell_loc(ir,ipos) + vol_cell *  uold(ind_cell(i),imetal)
                  !    !MZshell_loc(ir,ipos) = MZshell_loc(ir,ipos) + vol_cell *  unew(ind_cell(i),imetal)
                  !    ! print '(A,2I4, 1P2E13.3)', 'sweep_level, ir, MZshell_loc, uold(imetal), unew(imetal)', &
                  !    !      sweep_level, ir, uold(ind_cell(i),imetal), unew(ind_cell(i),imetal)
                  ! end if
                end if
                if (present(ETshell)) then
                  !E_th = E_tot - E_kin
                  ETshell_loc(ir,ipos) = ETshell_loc(ir,ipos) + vol_cell * (ufront(ind_cell(i),5) &
                    -  0.5d0/ufront(ind_cell(i),1) &
                    * (ufront(ind_cell(i),2)**2 + ufront(ind_cell(i),3)**2 + ufront(ind_cell(i),4)**2))
                end if
                if (present(RawMomRshell)) then
                  if (r > 0) then
                    !unit radial vecotr
                    erx = (x - pos(1,ipos))/r
                    ery = (y - pos(2,ipos))/r
                    erz = (z - pos(3,ipos))/r
                    !radial momentum (positive for outward direction)
                    MomR_cell = vol_cell * &
                      (erx*ufront(ind_cell(i),2) + ery*ufront(ind_cell(i),3) + erz*ufront(ind_cell(i),4))
                    RawMomRshell_loc(ir,ipos) = RawMomRshell_loc(ir,ipos) + MomR_cell
                    ! !KS DEBUG
                    ! if (ir == 2) then
                    !    print *, "KS DEBUG (get_shell_quantities) myid  = ", myid, ", ir = ", ir, ", i = ", i, ", ind_grid(i) = ", ind_grid(i), ", ind_cell(i) = ", ind_cell(i)

                    !    print *, "KS DEBUG (get_shell_quantities) rbin(ir-1) = ", rbin(ir-1), ", r = ", r, ", rbin(ir) = ", rbin(ir)
                    !    print *, "KS DEBUG (get_shell_quantities) erx = ", erx, ", ery = ", ery, ", erz = ", erz, &
                    !         ", M_cell [Msun] = ", vol_cell*ufront(ind_cell(i),1)*scale_msun, &
                    !         ", vel(:) [km/s]= ", ufront(ind_cell(i),2:4)/ufront(ind_cell(i),1)*scale_v/1d5, &
                    !         ", RawMomRshell_loc = ", RawMomRshell_loc(ir,ipos), ", MomR_cell = ", MomR_cell
                    ! end if

                  end if
                end if
                if (present(MomRshell)) then
                  MomXshell_loc(ir,ipos) = MomXshell_loc(ir,ipos) + vol_cell *ufront(ind_cell(i),2)
                  MomYshell_loc(ir,ipos) = MomYshell_loc(ir,ipos) + vol_cell *ufront(ind_cell(i),3)
                  MomZshell_loc(ir,ipos) = MomZshell_loc(ir,ipos) + vol_cell *ufront(ind_cell(i),4)
                end if


              end do
            end do
          end do
        end do
        ! End loop over cells
      end do
      ! End loop over grids
    end do
    ! End loop over levels


    !----------------------------------------
    ! collect global sum
    !----------------------------------------
    call MPI_ALLREDUCE(Mshell_loc(1:Nrbin,1:Npos), Mshell_glob(1:Nrbin,1:Npos), Nrbin*Npos, &
      MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    if (present(Nshell)) call MPI_ALLREDUCE(Nshell_loc(1:Nrbin,1:Npos), Nshell_glob(1:Nrbin,1:Npos), Nrbin*Npos, &
      MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
    if (present(Vshell)) call MPI_ALLREDUCE(Vshell_loc(1:Nrbin,1:Npos), Vshell_glob(1:Nrbin,1:Npos), Nrbin*Npos, &
      MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    if (present(MZshell)) call MPI_ALLREDUCE(MZshell_loc(1:Nrbin,1:Npos), MZshell_glob(1:Nrbin,1:Npos), Nrbin*Npos, &
      MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    if (present(ETshell)) call MPI_ALLREDUCE(ETshell_loc(1:Nrbin,1:Npos), ETshell_glob(1:Nrbin,1:Npos), Nrbin*Npos, &
      MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    if (present(RawMomRshell)) call MPI_ALLREDUCE(RawMomRshell_loc(1:Nrbin,1:Npos), RawMomRshell_glob(1:Nrbin,1:Npos), Nrbin*Npos, &
      MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)

    !---------------------------------------------
    ! obtain vbulk subtracted radial momentum
    !---------------------------------------------
    if (present(MomRshell)) then
      call MPI_ALLREDUCE(MomXshell_loc(1:Nrbin,1:Npos), MomXshell_glob(1:Nrbin,1:Npos), Nrbin*Npos, &
        MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
      call MPI_ALLREDUCE(MomYshell_loc(1:Nrbin,1:Npos), MomYshell_glob(1:Nrbin,1:Npos), Nrbin*Npos, &
        MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
      call MPI_ALLREDUCE(MomZshell_loc(1:Nrbin,1:Npos), MomZshell_glob(1:Nrbin,1:Npos), Nrbin*Npos, &
        MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
      do sweep_level = levelmin,nlevelmax
        dx_cell = 1d0/dble(2**sweep_level)
        vol_cell = (dx_cell*scale)**3   ! merge_ks: cell volume in code units for boxlen /= 1 (ks assumed boxlen=1) !cell volume at sweep_level

        ! Cells center position relative to grid center position
        do ind=1,twotondim
          iz=(ind-1)/4
          iy=(ind-1-4*iz)/2
          ix=(ind-1-2*iy-4*iz)
          xc_level(ind,1)=(dble(ix)-0.5D0)*dx_cell
          xc_level(ind,2)=(dble(iy)-0.5D0)*dx_cell
          xc_level(ind,3)=(dble(iz)-0.5D0)*dx_cell
        end do


        ! sum shell qunatities in each cpu
        ncache=active(sweep_level)%ngrid
        do igrid=1,ncache,nvector
          ngrid=MIN(nvector,ncache-igrid+1)
          do i=1,ngrid
            ind_grid(i)=active(sweep_level)%igrid(igrid+i-1)
          end do
          ! Loop over cells
          do ind=1,twotondim
            iskip=ncoarse+(ind-1)*ngridmax
            do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
            end do

            !store information about newly formed stars and their parent cells
            do i=1,ngrid

              !sum over only leaf cells
              if(son(ind_cell(i))/=0) cycle

              x = (xg(ind_grid(i),1)+xc_level(ind,1)-skip_loc(1))*scale          !x
              y = (xg(ind_grid(i),2)+xc_level(ind,2)-skip_loc(2))*scale          !y
              z = (xg(ind_grid(i),3)+xc_level(ind,3)-skip_loc(3))*scale          !z

              do ipos = 1, Npos
                r2 = (x - pos(1,ipos))**2 +  (y - pos(2,ipos))**2 +  (z - pos(3,ipos))**2
                r = sqrt(r2)


                do ir = 1, Nrbin
                  !skip unless the cell is in the shell
                  if (ir == 1) then !0 < r <= rbin(1) for initial shell
                    if (r > rbin(1)) cycle
                  else                !rbin(i-1) < r <= rbin(i) otherwise
                    if (r <= rbin(ir-1) .or. r > rbin(ir)) cycle
                  end if

                  if (r > 0) then
                    !unit radial vecotr
                    erx = (x - pos(1,ipos))/r
                    ery = (y - pos(2,ipos))/r
                    erz = (z - pos(3,ipos))/r

                    !bulk velocity vector
                    rho = ufront(ind_cell(i),1)
                    vbulk_x = MomXshell_glob(ir,ipos)/Mshell_glob(ir,ipos)
                    vbulk_y = MomYshell_glob(ir,ipos)/Mshell_glob(ir,ipos)
                    vbulk_z = MomZshell_glob(ir,ipos)/Mshell_glob(ir,ipos)

                    !radial momentum (positive for outward direction)
                    MomR_cell = vol_cell * &
                      (erx*(ufront(ind_cell(i),2)-rho*vbulk_x) + ery*(ufront(ind_cell(i),3)-rho*vbulk_y) + erz*(ufront(ind_cell(i),4)-rho*vbulk_z))
                    MomRshell_loc(ir,ipos) = MomRshell_loc(ir,ipos) + MomR_cell
                    ! !KS DEBUG
                    ! if (ir == 2) then
                    !    print *, "KS DEBUG (get_shell_quantities) myid  = ", myid, ", ir = ", ir, ", i = ", i, ", ind_grid(i) = ", ind_grid(i), ", ind_cell(i) = ", ind_cell(i)

                    !    print *, "KS DEBUG (get_shell_quantities) rbin(ir-1) = ", rbin(ir-1), ", r = ", r, ", rbin(ir) = ", rbin(ir)
                    !    print *, "KS DEBUG (get_shell_quantities) vbulk(:) [km/s] = ", (/vbulk_x,vbulk_y,vbulk_z/)*scale_v/1d5
                    !    print *, "KS DEBUG (get_shell_quantities) erx = ", erx, ", ery = ", ery, ", erz = ", erz, &
                    !         ", M_cell [Msun] = ", vol_cell*ufront(ind_cell(i),1)*scale_msun, &
                    !         ", vel(:) [km/s]= ", ufront(ind_cell(i),2:4)/ufront(ind_cell(i),1)*scale_v/1d5, &
                    !         ", MomRshell_loc = ", MomRshell_loc(ir,ipos), ", MomR_cell = ", MomR_cell
                    ! end if

                  end if
                end do
              end do
            end do
          end do
          ! End loop over cells
        end do
        ! End loop over grids
      end do
      ! End loop over levels

      call MPI_ALLREDUCE(MomRshell_loc(1:Nrbin,1:Npos), MomRshell_glob(1:Nrbin,1:Npos), Nrbin*Npos, &
        MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)


    end if !if (present(MomRshell))


    !return data
    Mshell(1:Nrbin,1:Npos) = Mshell_glob(1:Nrbin,1:Npos)
    if (present(Nshell)) Nshell(1:Nrbin,1:Npos) = Nshell_glob(1:Nrbin,1:Npos)
    if (present(Vshell)) Vshell(1:Nrbin,1:Npos) = Vshell_glob(1:Nrbin,1:Npos)
    if (present(MZshell)) MZshell(1:Nrbin,1:Npos) = MZshell_glob(1:Nrbin,1:Npos)
    if (present(ETshell)) ETshell(1:Nrbin,1:Npos) = ETshell_glob(1:Nrbin,1:Npos)
    if (present(RawMomRshell)) RawMomRshell(1:Nrbin,1:Npos) = RawMomRshell_glob(1:Nrbin,1:Npos)
    if (present(MomRshell)) MomRshell(1:Nrbin,1:Npos) = MomRshell_glob(1:Nrbin,1:Npos)

    ! !KS DEBUG
    ! if(myid==1) then
    !    print *, "KS DEBUG (get_shell_quantities) Nshell_glob(1:Nrbin,1) = ", Nshell_glob(1:Nrbin,1)
    !    print *, "KS DEBUG (get_shell_quantities) Vshell_glob(1:Nrbin,1) = ", Vshell_glob(1:Nrbin,1)
    !    print *, "KS DEBUG (get_shell_quantities) Mshell_glob(1:Nrbin,1) = ", Mshell_glob(1:Nrbin,1)
    !    !print *, "KS DEBUG (get_shell_quantities) RawMomRshell_glob(1:Nrbin,1) = ", RawMomRshell_glob(1:Nrbin,1)
    !    print *, "KS DEBUG (get_shell_quantities) RawVelR_shell(1:Nrbin)[km/s] = ", RawMomRshell_glob(1:Nrbin,1)/Mshell_glob(1:Nrbin,1)*scale_v/1d5
    !    print *, "KS DEBUG (get_shell_quantities) VelR_shell(1:Nrbin)[km/s] = ", MomRshell_glob(1:Nrbin,1)/Mshell_glob(1:Nrbin,1)*scale_v/1d5
    !    print *, "KS DEBUG (get_shell_quantities) STOP HERE"
    !    stop
    ! end if

  end subroutine get_shell_quantities

  !----------------------------------------------------------------------------
  ! modify gas quantities in a spherical region by a given factor
  !----------------------------------------------------------------------------
  subroutine modify_sphere_by_factor(pos,radius,f_rho,rmin)

    !in
    real(dp),dimension(1:3),intent(in) :: pos  !center
    real(dp),intent(in) :: radius !radius of modification region
    real(dp),intent(in) :: f_rho  !ratio of new and original density: rho -> f_rho * rho

    real(dp),intent(in),optional :: rmin !minimum radius of modification region (optional)


    !variables
    integer :: sweep_level !search for all leaf cells through levels (not only ilevel)
    real(dp) :: x, y, z, r2, r
    real(dp),dimension(1:twotondim,1:3)::xc_level
    integer :: ix, iy, iz
    real(dp) :: dx_cell
    real(dp) :: radius_min !minimum radius (to be set to rmin, if given)

    !set minimum radius
    if (present(rmin)) then
      radius_min = rmin
    else
      radius_min = -1d0 !dummy value (r is always larger than radius_min)
    end if

    do sweep_level = levelmin,nlevelmax
      dx_cell = 1d0/dble(2**sweep_level)

      ! Cells center position relative to grid center position
      do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc_level(ind,1)=(dble(ix)-0.5D0)*dx_cell
        xc_level(ind,2)=(dble(iy)-0.5D0)*dx_cell
        xc_level(ind,3)=(dble(iz)-0.5D0)*dx_cell
      end do


      ! sum shell qunatities in each cpu
      ncache=active(sweep_level)%ngrid
      do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
          ind_grid(i)=active(sweep_level)%igrid(igrid+i-1)
        end do
        ! Loop over cells
        do ind=1,twotondim
          iskip=ncoarse+(ind-1)*ngridmax
          do i=1,ngrid
            ind_cell(i)=iskip+ind_grid(i)
          end do

          !store information about newly formed stars and their parent cells
          do i=1,ngrid

            !modification only in the leaf cells (OK? KS TODO)
            if(son(ind_cell(i))/=0) cycle

            x = (xg(ind_grid(i),1)+xc_level(ind,1)-skip_loc(1))*scale          !x
            y = (xg(ind_grid(i),2)+xc_level(ind,2)-skip_loc(2))*scale          !y
            z = (xg(ind_grid(i),3)+xc_level(ind,3)-skip_loc(3))*scale          !z

            r2 = (x - pos(1))**2 +  (y - pos(2))**2 +  (z - pos(3))**2
            r = sqrt(r2)

            !modify gas inside the sphere region with (radius_min <) r < radius
            ! Note:
            ! - at sweep_level = ilevel (= nlevelmax), all uold was converted to primitive variables,
            !   and will be converted back to conserved variables later
            ! - for lower levels, both uold and unew need to be updated, because uold will be used to compute flux
            !   while unew contains the difference in the current step
            !if (r < radius) then
            if (radius_min < r .and. r < radius) then
              !update ufront
              do ivar=1,nvar
                ufront(ind_cell(i),ivar) = f_rho * ufront(ind_cell(i),ivar)
              end do

              !update uold/unew
              if (sweep_level == nlevelmax) then                    !at ilevel (= nlevelmax)
                !change only density (others are automatically updated)
                uold(ind_cell(i),1) = ufront(ind_cell(i),1)
              else                                                  !at lower levels
                !change all conservative vars
                do ivar=1,nvar
                  !update uold (to be used for calculating flux by godunov method)
                  uold(ind_cell(i),ivar) = ufront(ind_cell(i),ivar)
                  !updat unew (variation of ufront used to ensure conservation)
                  unew(ind_cell(i),ivar) = unew(ind_cell(i),ivar) - (1d0-f_rho)/f_rho * ufront(ind_cell(i),ivar)
                end do
              end if
            end if
          end do
        end do
        ! End loop over cells
      end do
      ! End loop over grids
    end do
    ! End loop over levels
  end subroutine modify_sphere_by_factor

  !----------------------------------------------------------------------------
  ! modify gas quantities in a spherical region to costant density
  !----------------------------------------------------------------------------
  subroutine modify_sphere_to_costant_density(pos,radius,rho,rmin)

    !in
    real(dp),dimension(1:3),intent(in) :: pos  !center
    real(dp),intent(in) :: radius !radius of modification region
    real(dp),intent(in) :: rho    !target density:  ->  rho

    real(dp),intent(in),optional :: rmin !minimum radius of modification region (optional)


    !variables
    integer :: sweep_level !search for all leaf cells through levels (not only ilevel)
    real(dp) :: x, y, z, r2, r
    real(dp),dimension(1:twotondim,1:3)::xc_level
    integer :: ix, iy, iz
    real(dp) :: dx_cell
    real(dp) :: radius_min !minimum radius (to be set to rmin, if given)
    real(dp) :: f_rho !locally defined ratio of rho_prev and rho

    !set minimum radius
    if (present(rmin)) then
      radius_min = rmin
    else
      radius_min = -1d0 !dummy value (r is always larger than radius_min)
    end if

    do sweep_level = levelmin,nlevelmax
      dx_cell = 1d0/dble(2**sweep_level)

      ! Cells center position relative to grid center position
      do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc_level(ind,1)=(dble(ix)-0.5D0)*dx_cell
        xc_level(ind,2)=(dble(iy)-0.5D0)*dx_cell
        xc_level(ind,3)=(dble(iz)-0.5D0)*dx_cell
      end do


      ! sum shell qunatities in each cpu
      ncache=active(sweep_level)%ngrid
      do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
          ind_grid(i)=active(sweep_level)%igrid(igrid+i-1)
        end do
        ! Loop over cells
        do ind=1,twotondim
          iskip=ncoarse+(ind-1)*ngridmax
          do i=1,ngrid
            ind_cell(i)=iskip+ind_grid(i)
          end do

          !store information about newly formed stars and their parent cells
          do i=1,ngrid

            !modification only in the leaf cells (OK? KS TODO)
            if(son(ind_cell(i))/=0) cycle

            x = (xg(ind_grid(i),1)+xc_level(ind,1)-skip_loc(1))*scale          !x
            y = (xg(ind_grid(i),2)+xc_level(ind,2)-skip_loc(2))*scale          !y
            z = (xg(ind_grid(i),3)+xc_level(ind,3)-skip_loc(3))*scale          !z

            r2 = (x - pos(1))**2 +  (y - pos(2))**2 +  (z - pos(3))**2
            r = sqrt(r2)

            !modify gas inside the sphere region with (radius_min <) r < radius
            ! Note:
            ! - at sweep_level = ilevel (= nlevelmax), all uold was converted to primitive variables,
            !   and will be converted back to conserved variables later
            ! - for lower levels, both uold and unew need to be updated, because uold will be used to compute flux
            !   while unew contains the difference in the current step
            !if (r < radius) then
            if (radius_min < r .and. r < radius) then
              !get f_rho
              f_rho = rho/ufront(ind_cell(i),1)

              !update ufront
              do ivar=1,nvar
                ufront(ind_cell(i),ivar) = f_rho * ufront(ind_cell(i),ivar)
              end do

              !update uold/unew
              if (sweep_level == nlevelmax) then                    !at ilevel (= nlevelmax)
                !change only density (others are automatically updated)
                uold(ind_cell(i),1) = ufront(ind_cell(i),1)
              else                                                  !at lower levels
                !change all conservative vars
                do ivar=1,nvar
                  !update uold (to be used for calculating flux by godunov method)
                  uold(ind_cell(i),ivar) = ufront(ind_cell(i),ivar)
                  !updat unew (variation of ufront used to ensure conservation)
                  unew(ind_cell(i),ivar) = unew(ind_cell(i),ivar) - (1d0-f_rho)/f_rho * ufront(ind_cell(i),ivar)
                end do
              end if
            end if
          end do
        end do
        ! End loop over cells
      end do
      ! End loop over grids
    end do
    ! End loop over levels
  end subroutine modify_sphere_to_costant_density



  !------------------------------------------------------------------------------------------------
  ! subroutine to get SFC list
  !------------------------------------------------------------------------------------------------
  ! - get SFC from their candidates
  ! - candidates simultaneously forming in close vicinities are removed
  !  (copy-and-pasted from remove_p3cand_by_other_cand with modification)
  !------------------------------------------------------------------------------------------------
  subroutine SFC_from_candidate()

    integer :: iloc, iglob, icpu
    integer :: offset !offset between indices of local and global lists
    integer,dimension(1:ncpu) :: SFC_ncand_list, SFC_n_list !list for number of SFC (candidates) in each cpu
    logical,dimension(1:SFC_n_max) :: flg_rmv !remove candidate if true

    !MPI send/recv buffer
    integer,dimension(1:ncpu) :: sbufi,rbufi
    real(dp),dimension(1:SFC_n_max*3) :: sbufr,rbufr

    !-----------------------------------
    ! create global list of SFCs
    !-----------------------------------
    !create SFC_ncand_list (global)
    SFC_ncand_list(:) = 0
    SFC_ncand_list(myid) = SFC_ncand_loc
    sbufi(1:ncpu) = SFC_ncand_list(:)
    call MPI_ALLREDUCE(sbufi,rbufi,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
    SFC_ncand_list(:) = rbufi(1:ncpu)
    SFC_ncand_tot = sum(SFC_ncand_list) !total number of candidates


    if(myid==1 .and. SFC_ncand_tot > 0) then
      write(*,'(A,I5,A)') 'Check simultaneous SFC formation in close vicinity for',&
        SFC_ncand_tot, ' candidates'
    end if


    ! !KS DEBUG
    ! write(*,*) '(KS DEBUG) myid, SFC_ncand_tot: ', myid, SFC_ncand_tot
    ! write(*,*) '(KS DEBUG) myid, SFC_ncand_list(:): ',myid,SFC_ncand_list(:)

    ! check the number of candidates
    if (SFC_ncand_tot > SFC_n_max) then
      write(*,*)'Number of global SFC candidates exceeds the maximum ',SFC_ncand_tot, SFC_n_max
      write(*,*)'stopping...'
      stop
    end if

    !exit from subroutine if there is no candidate
    if (SFC_ncand_tot == 0) then
      SFC_n_tot = 0
      return
    end if

    !determine offset between local and global cndidate lists
    offset = 0
    do icpu = 1, myid-1  ! "do i = 1,0" loop is just skipped (KS MEMO)
      offset = offset + SFC_ncand_list(icpu)
    end do

    !create global candidate position list
    SFC_poscand_glob(:,:)=0 !initialization
    do iloc = 1, SFC_ncand_loc
      iglob = iloc + offset
      SFC_poscand_glob(:,iglob) = SFC_poscand_loc(:,iloc)
    end do
    do iglob = 1, SFC_ncand_tot
      sbufr((iglob-1)*3+1) = SFC_poscand_glob(1,iglob)
      sbufr((iglob-1)*3+2) = SFC_poscand_glob(2,iglob)
      sbufr((iglob-1)*3+3) = SFC_poscand_glob(3,iglob)
    end do
    call MPI_ALLREDUCE(sbufr, rbufr, 3*SFC_ncand_tot, MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    do iglob = 1, SFC_ncand_tot
      SFC_poscand_glob(1,iglob) = rbufr((iglob-1)*3+1)
      SFC_poscand_glob(2,iglob) = rbufr((iglob-1)*3+2)
      SFC_poscand_glob(3,iglob) = rbufr((iglob-1)*3+3)
    end do


    !-----------------------------------
    ! remove candidates in close distance
    !-----------------------------------
    flg_rmv(:) = .false.
    ! TO ADDED 2026-11-11: SFC markers (local Pop II models): no new marker within SFC_dist_min of an active one
    ! (the SFC model needs no such check: the gas of a formed SFC is cropped)
    if (pop2_SFC_marker) then
      call collect_active_SFC   ! SFC_n_active, SFC_pos_active (pm_commons): all FAM_SFC particles
      do iloc = 1, SFC_ncand_loc
        do jglob = 1, SFC_n_active
          r2 = sum((SFC_poscand_loc(1:3,iloc) - SFC_pos_active(1:3,jglob))**2)
          if (r2 < SFC_dist_min**2) then
            flg_rmv(iloc) = .true.
            exit
          end if
        end do
      end do
    end if
    ! check if SFC (upper in list) forms near current candidate
    do iloc = 1, SFC_ncand_loc
      iglob = iloc + offset
      ! rejected candidates are included in this loop but it should not significantly affect the results
      ! because not all candidates are removed (KS MEMO)
      do jglob = 1,iglob-1
        r2 = (SFC_poscand_glob(1,iglob) - SFC_poscand_glob(1,jglob))**2 + &
          (SFC_poscand_glob(2,iglob) - SFC_poscand_glob(2,jglob))**2 + &
          (SFC_poscand_glob(3,iglob) - SFC_poscand_glob(3,jglob))**2 !distance

        ! !KS DEBUG
        ! print *, 'iloc, iglob, jglob', iloc, iglob, jglob
        ! print '(A,6(1PE13.3))', "SFC_poscand_glob(:,iglob), SFC_poscand_glob(:,jglob)", SFC_poscand_glob(:,iglob), SFC_poscand_glob(:,jglob)
        ! print *, "r2, SFC_dist_min**2: ", r2, SFC_dist_min**2

        ! if too close, undo the SFC formation
        if (r2 < SFC_dist_min**2) then
          flg_rmv(iloc) = .true.
          !flag2(SFC_indcellcand_loc(iloc)) = 0
          !ntot = ntot - 1

          print '(A,2(A,3(1PE12.5)),A,1PE12.5)', "** Remove SFC cand by other cand **", &
            ' pos_remove(:): ', SFC_poscand_glob(:,iglob), &
            ', pos_remain(:): ', SFC_poscand_glob(:,jglob),&
            ', r/rmin: ', sqrt(r2)/SFC_dist_min

          exit !jump out of the loop
        end if
      end do
    end do


    !-----------------------------------
    ! save SFCs surviving the removal process
    !-----------------------------------

    !SFCs surviving the removal process (local)
    SFC_n_loc = 0
    do iloc = 1, SFC_ncand_loc
      !if (flag2(SFC_indcellcand_loc(iloc)) > 0) then
      if (.not. flg_rmv(iloc)) then
        SFC_n_loc = SFC_n_loc + 1           !number of SFC after removal
        SFC_pos_loc(:,SFC_n_loc) = SFC_poscand_loc(:,iloc) !position of SFC center
      end if
    end do

    !list of SFC number (global)
    SFC_n_list(:) = 0
    SFC_n_list(myid) = SFC_n_loc
    sbufi(1:ncpu) = SFC_n_list(:)
    call MPI_ALLREDUCE(sbufi,rbufi,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
    SFC_n_list(:) = rbufi(1:ncpu)
    SFC_n_tot = sum(SFC_n_list) !total number of candidates

    !offset between local and global lists
    offset = 0
    do icpu = 1, myid-1  ! "do i = 1,0" loop is just skipped (KS MEMO)
      offset = offset + SFC_n_list(icpu)
    end do

    !global SFC position list
    SFC_pos_glob(:,:)=0 !initialization
    do iloc = 1, SFC_n_loc
      iglob = iloc + offset
      SFC_pos_glob(:,iglob) = SFC_pos_loc(:,iloc)
    end do
    do iglob = 1, SFC_n_tot
      sbufr((iglob-1)*3+1) = SFC_pos_glob(1,iglob)
      sbufr((iglob-1)*3+2) = SFC_pos_glob(2,iglob)
      sbufr((iglob-1)*3+3) = SFC_pos_glob(3,iglob)
    end do
    call MPI_ALLREDUCE(sbufr, rbufr, 3*SFC_n_tot, MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    do iglob = 1, SFC_n_tot
      SFC_pos_glob(1,iglob) = rbufr((iglob-1)*3+1)
      SFC_pos_glob(2,iglob) = rbufr((iglob-1)*3+2)
      SFC_pos_glob(3,iglob) = rbufr((iglob-1)*3+3)
    end do

    !KS DEBUG
    if (myid == 1) then
      print '(A,I5)', 'SFC_from_candidate:  SFC_n_tot = ', SFC_n_tot
      do iglob = 1, SFC_n_tot
        print '(A,I5,1P3E15.7)', 'iglob, SFC_pos = ', iglob, SFC_pos_glob(:,iglob)
      end do
    end if


  end subroutine SFC_from_candidate

  !----------------------------------------------------------------------
  ! TO ADDED 2026-11-11: put a cell on the SFC candidate list (SFC markers in the local Pop II models)
  !----------------------------------------------------------------------
  subroutine store_SFC_candidate(icell, pos)
    integer, intent(IN) :: icell
    real(dp), dimension(1:3), intent(IN) :: pos
    SFC_ncand_loc = SFC_ncand_loc + 1
    if (SFC_ncand_loc > SFC_n_max) then
      write(*,*)'Number of local SFC candidates exceeds the limit ',SFC_ncand_loc, SFC_n_max
      write(*,*)'stopping...'
      stop
    end if
    SFC_indcellcand_loc(SFC_ncand_loc) = icell
    SFC_poscand_loc(1:3,SFC_ncand_loc) = pos(1:3)
    ! (no message per candidate: in a local model the cells near an active marker stay candidates for SFC_t_life;
    !  SFC_from_candidate reports the counts and create_SFC_markers the new markers)
  end subroutine store_SFC_candidate

  !----------------------------------------------------------------------
  ! TO ADDED 2026-11-11: SFC marker particles (local Pop II models): for each SFC of SFC_pos_glob (from
  ! SFC_from_candidate) write the formation nest as for an SFC and add a massless FAM_SFC particle at the
  ! candidate cell (as star_formation_in_SFC does for an SFC; radius 0, no star formation, no cropping).
  ! Afterwards it is treated like an SFC: age outputs (output_around_SFC_at_given_age), FAM_PSC after
  ! SFC_t_life, counted in n_sfc, and while active it blocks new markers within SFC_dist_min_pc.
  !----------------------------------------------------------------------
  subroutine create_SFC_markers()
    integer :: cell_levl1, cell_index1, index_star
    integer, dimension(1:nvector) :: ind_grid1, ind_part1
    logical, dimension(1:nvector) :: ok_new1
    real(dp), dimension(1:3) :: pos1
    ok_new1(:) = .true.
    do iSFC = 1, SFC_n_tot
      call output_around_SFC(iSFC, SFC_pos_glob(1:3,iSFC), use_current_ufront=.true.)  !formation nest, as for an SFC
      nstar_tot = nstar_tot + 1
      pos1(1:3) = SFC_pos_glob(1:3,iSFC)
      call get_active_cell_index_highest_level(pos1(1:3),cell_levl1,cell_index1)
      if (cell_index1 >= 0) then           !the cpu with the active grid containing the position
        ind_grid1(1)=mod(cell_index1-ncoarse-1,ngridmax)+1
        call remove_free(ind_part1,1)
        call add_list(ind_part1,ind_grid1,ok_new1,1)
        index_star=nstar_tot
        tp(ind_part1(1))=birth_epoch
        mp(ind_part1(1))=0d0
        levelp(ind_part1(1))=cell_levl1
        idp(ind_part1(1))=index_star
        typep(ind_part1(1))%family=FAM_SFC
        xp(ind_part1(1),1:3)=pos1(1:3)
        vp(ind_part1(1),1)=ufront(cell_index1,2)/ufront(cell_index1,1)
        vp(ind_part1(1),2)=ufront(cell_index1,3)/ufront(cell_index1,1)
        vp(ind_part1(1),3)=ufront(cell_index1,4)/ufront(cell_index1,1)
        zp(ind_part1(1)) = 0d0   ! -radius for an SFC; a marker has no cloud radius
#ifdef INIT_STELLAR_MASS
        mp0(ind_part1(1))=mp(ind_part1(1))
#endif
#ifdef STELLAR_POPULATION_MASS
        msp0(ind_part1(1))=mp(ind_part1(1))
#endif
        if(sf_log_properties) call sf_log_write_star(ilun_sf,ind_part1(1),cell_index1,ilevel,(cell_levl1==nlevelmax))
        write(*,'(A,2I4,I10,A,I8,A,I2,A,1P3E9.1,A,1P3E16.8)') &
          "** New SFC marker particle added **  myid, lev, ind_cell: ", myid, cell_levl1, cell_index1, &
          ", idp:", idp(ind_part1(1)), ", fm:", typep(ind_part1(1))%family, &
          ", vp:", vp(ind_part1(1),:)*scale_l/scale_t/1d5, ", xp:", xp(ind_part1(1),:)
      end if
      if (myid == 1) then
        print '(A,I5,A,I5,A,I10,A,1P3E20.12)', '** SFC marker **  iSFC/nSFC: ', iSFC, '/', SFC_n_tot, &
          ', idp = ', nstar_tot, ', pos = ', SFC_pos_glob(1:3,iSFC)
      end if
    end do
  end subroutine create_SFC_markers

  !------------------------------------------------------------------------------------------------
  ! subroutine to deal with formation in SFC
  !------------------------------------------------------------------------------------------------
  ! - construct 1D profile around each SFC
  ! - define size and mass of SFC
  ! - convert SFC to stars (+ feedback on SFC)
  !------------------------------------------------------------------------------------------------
  subroutine star_formation_in_SFC

    integer,parameter :: Nrbin = 60 ! number of radial bins for shells (r_i = r_0 * (r_max/r_0)**(i/(Nrbin-1)))
    real(dp),parameter :: f_rmax = 1d3 !rmax/rmin: the ratio of minmum and maximum radii of radial bins (rmin is set to dx_min)
    !integer,parameter :: Nrbin = 20 ! number of radial bins for shells (r_i = r_0 * (r_max/r_0)**(i/(Nrbin-1)))

    !integer,parameter :: Nrbin = 20 ! number of radial bin for defining shells (r_i/dx_min = 1,2,4,,,2**(Nrbin-1))
    !real(dp),parameter :: f_rmax = 1d2 !rmax/rmin: the ratio of minmum and maximum radii of radial bins (rmin is set to dx_min)

    real(dp),dimension(1:Nrbin) :: rbin      !outer radii of shells
    real(dp),allocatable,dimension(:,:) :: Mshell    !mass of shells
    integer,allocatable,dimension(:,:) :: Nshell    !number of cell in shells
    real(dp),allocatable,dimension(:,:) :: Vshell    !volume of shells
    real(dp),allocatable,dimension(:,:) :: RHOshell   !density of shells
    real(dp),allocatable,dimension(:,:) :: MZshell    !metal mass of shells
    real(dp),allocatable,dimension(:,:) :: ETshell    !thermal energy of shells
    real(dp),allocatable,dimension(:,:) :: T2shell    !T2 = T/mu of shells
    real(dp),dimension(1:Nrbin) :: Menc, Venc, RHOenc, MZenc, Zenc  !enclosed mass, voulume, rho(ave), metal mass, metal(ave)
    real(dp) :: radius, mass, rho, Zmetal
    integer :: iSFC, ir
    real(dp) :: rho_SF !threshold density for SF (limit SF to rho > rho_SF)

    integer,parameter :: LUN = 11

    real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v, scale_msun, scale_pc

    real(dp) :: nH_bnd !threshold density for SFC outer boundary

    logical ok_flag ! False if we cannot make stars for a given rho_SF

    real(dp) :: m_star_form
    integer :: n_star_form

    !for adding SFC to particle list
    integer,dimension(1:nvector)::cell_index,cell_levl,cc
    !real(dp),dimension(1:nvector,1:3)::pos
    real(dp),dimension(1:3):: pos1
    integer :: lev, cell_index1, cell_levl1

    !MPI send/recv buffer
    integer,dimension(1:ncpu) :: sbufi,rbufi
    real(dp),dimension(1:SFC_n_max*3) :: sbufr,rbufr

    integer :: sweep_level !KS DEBUG

    ! Conversion factor from user units to cgs units
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    scale_msun=scale_d*scale_l**3/m_sun
    scale_pc=scale_l/3.08568d18

    ! !KS DEBUG TENTATIVE
    ! do sweep_level = levelmin,nlevelmax
    !    do ind=1,twotondim
    !       iskip=ncoarse+(ind-1)*ngridmax
    !       do i=1,active(sweep_level)%ngrid
    !          if (myid == 5 .and. active(sweep_level)%igrid(i)+iskip == 4821) then
    !             print '(A,2I4,I8,1P3E12.4)', &
    !                  '(KS DEBUG TENTATIVE at the beginning of star_formation_in_SFC) myid, level, ind_cell, nH(uold), nH(unew), nH(ufront): ',&
    !                  myid, sweep_level, active(sweep_level)%igrid(i)+iskip, uold(active(sweep_level)%igrid(i)+iskip,1)*scale_nH, unew(active(sweep_level)%igrid(i)+iskip,1)*scale_nH,&
    !                  ufront(active(sweep_level)%igrid(i)+iskip,1)*scale_nH
    !          end if
    !       end do
    !    end do
    ! end do

    !------------------------------------------
    ! construct 1D profile around each SFC
    !------------------------------------------
    !define radial bins
    do ir = 1, Nrbin
      !rbin(ir) = dx_min * (1d0-1d-8) * 2**(ir-1) !put an extra factor (<1) to avoid ambiguity in defining shell regions
      rbin(ir) = dx_min * (1d0-1d-8) * f_rmax**(dble(ir-1)/(Nrbin-1)) !put an extra factor (<1) to avoid ambiguity in defining shell regions
    end do

    allocate(Mshell(Nrbin,SFC_n_tot),Nshell(Nrbin,SFC_n_tot),Vshell(Nrbin,SFC_n_tot),RHOshell(Nrbin,SFC_n_tot),&
      MZshell(Nrbin,SFC_n_tot),ETshell(Nrbin,SFC_n_tot),T2shell(Nrbin,SFC_n_tot))
    !get shell quantities around each SFC (WARNING: overlap of SFCs is not allowed in the current version of code)
    call get_shell_quantities(SFC_n_tot, SFC_pos_glob(1:3,1:SFC_n_tot), Nrbin, rbin(1:Nrbin), Mshell(1:Nrbin,1:SFC_n_tot),&
      Nshell=Nshell(1:Nrbin,1:SFC_n_tot), Vshell=Vshell(1:Nrbin,1:SFC_n_tot), &
      MZshell=MZshell(1:Nrbin,1:SFC_n_tot), ETshell=ETshell(1:Nrbin,1:SFC_n_tot))

    !get RHOshell (and T2shell, which is currently not used)
    do iSFC = 1, SFC_n_tot
      do ir = 1, Nrbin
        if (Nshell(ir,iSFC)>0) then  !if there is at least one cell in the shell
          RHOshell(ir,iSFC) = Mshell(ir,iSFC)/Vshell(ir,iSFC)
          T2shell(ir,iSFC) = (gamma-1.0) * ETshell(ir,iSFC)/Mshell(ir,iSFC)
        else                         !if thre is no cell in the shell
          RHOshell(ir,iSFC) = -1d0/scale_nH  !negative value for ILL-DEFINED case
          T2shell(ir,iSFC) = -1d0/scale_T2
        end if
      end do
    end do

    !KS DEBUG
    if (myid == 1) then
      do iSFC = 1, SFC_n_tot
        if (iSFC == 1) then
          print '(A,I5)', '(star_formation_in_SFC) SFC_n_tot:', SFC_n_tot
          print '(A,I5)', 'Nrbin: ', Nrbin
        end if
        print '(A,1P3E20.12)', 'SFC_pos(:): ', SFC_pos_glob(1:3,iSFC)
        print '(A,1P100E9.2)', 'rbin(:):    ', rbin(:)*scale_pc
        print '(A,100I9)',     'Nshell(:):  ', Nshell(1:Nrbin,iSFC)
        print '(A,1P100E9.2)', 'Mshell(:):  ', Mshell(1:Nrbin,iSFC)*scale_msun
        print '(A,1P100E9.2)', 'Vshell(:):  ', Vshell(1:Nrbin,iSFC)*scale_pc**3
        print '(A,1P100E9.2)', 'RHOshell(:):', RHOshell(1:Nrbin,iSFC)*scale_nH
        print '(A,1P100E9.2)', 'MZshell(:): ', MZshell(1:Nrbin,iSFC)*scale_msun
        print '(A,1P100E9.2)', 'T2shell(:): ', T2shell(1:Nrbin,iSFC)*scale_T2
      end do
    end if

    ! !sanity check
    ! do iSFC = 1, SFC_n_tot
    !    do ir = 1, Nrbin-1
    !       if (RHOshell(ir,iSFC) < 0 .and. RHOshell(ir+1,iSFC) > 0) then
    !          if (myid == 1) then
    !             print '(A,/,A)', 'well-defined shells are sepparated. somehting unexpected may happen.',&
    !                  'stopping...'
    !             stop
    !          end if
    !       end if
    !    end do
    ! end do

    !loop over SFCs
    do iSFC = 1, SFC_n_tot

      !------------------------------------------
      ! define size and mass of SFC
      !------------------------------------------
      !calculate Menc, rho_ave
      Menc(:) = 0d0
      Venc(:) = 0d0
      MZenc(:) = 0d0
      do ir = 1, Nrbin
        if (ir == 1) then
          Menc(ir) = Mshell(1,iSFC)
          Venc(ir) = Vshell(1,iSFC)
          MZenc(ir) = MZshell(1,iSFC)
        else
          Menc(ir) = Menc(ir-1) + Mshell(ir,iSFC)
          Venc(ir) = Venc(ir-1) + Vshell(ir,iSFC)
          MZenc(ir) = MZenc(ir-1) + MZshell(ir,iSFC)
        end if
      end do
      RHOenc(:) = Menc(:)/Venc(:) !average density wtihin the radius
      Zenc(:) = MZenc(:)/Menc(:) !average metallicity wtihin the radius

      ! !**************** WARNING WARNING WARNING WARNING WARNING WARNING *****************!
      ! print *, "*** WARNING *** KS TENTATIVE DEBUG *** WARNING *** "
      ! print *, "RHOshell artificially modified"
      ! print *, "*** WARNING *** KS TENTATIVE DEBUG *** WARNING *** "
      ! RHOshell(1:5,iSFC) = 1d6/scale_nH
      ! RHOshell(6:Nrbin,iSFC) = 1d0/scale_nH
      ! !**************** WARNING WARNING WARNING WARNING WARNING WARNING *****************!

      !threshold density of SFC outer boundary
      ! if (SFC_nH_min > 0) then !fixed density threshold for outer boundary
      !   nH_bnd  = SFC_nH_min
      ! else
      !   nH_bnd  = SFC_f_nH_min * RHOshell(1,iSFC) * scale_nH !density threshold proportional to central density
      ! end if

      if (SFC_nH_min > 0 .and. SFC_f_nH_min > 0) then
        nH_bnd  = max(SFC_nH_min, SFC_f_nH_min * RHOshell(1,iSFC) * scale_nH) !take the maximum of fixed and proportional density thresholds
      else if (SFC_nH_min > 0) then
        nH_bnd  = SFC_nH_min !fixed density threshold for outer boundary
      else
        nH_bnd  = SFC_f_nH_min * RHOshell(1,iSFC) * scale_nH !density threshold proportional to central density
      end if



      !dtermine SFC radius (maximum radius whith which nH > cloud_nHmin)
      radius = -1d0 !dummy
      do ir = 1, Nrbin
        if (RHOshell(ir,iSFC) < 0) cycle !skip undefined shell
        !if (RHOshell(ir,iSFC)*scale_nH < SFC_nH_min) then  !find outer radius
        if (RHOshell(ir,iSFC)*scale_nH < nH_bnd) then  !find outer radius
          if (rbin(ir) < SFC_rsize_min) cycle !impose minimum radius of SFC
          radius = rbin(ir)
          mass = Menc(ir)
          rho = RHOenc(ir)
          Zmetal = Zenc(ir)
          exit
        end if
      end do
      !check if radius is found
      if (radius < 0) then
        if (myid == 1) then
          !------------- ERROR VERSION ------------!
          ! print '(A,I5,A,1PE11.3/,A)', '*** ERROR *** SFC radius cannot be defined in get_SFC_properties. ir = ', ir, &
          !      ', nH_bnd = ', nH_bnd,&
          !      'stopping...'
          ! stop

          !------------- WARNING VERSION ----------!
          print '(A,I5,A,1PE11.3/,A)', '*** WARNING *** SFC radius is set to the largest rbin since nH is always above nH_bnd in get_SFC_properties. ir = ', ir, &
            ', nH_shell(Nrbin) = ', RHOshell(Nrbin,iSFC)*scale_nH, ', nH_bnd = ', nH_bnd
        end if
        !set SFC radius to outermost radial bin
        ir = Nrbin
        radius = rbin(ir)
        mass = Menc(ir)
        rho = RHOenc(ir)
        Zmetal = Zenc(ir)

      end if
      !check if SFC touches the boundary (TODO: compatibility with periodic boundary)
      if (MINVAL(SFC_pos_glob(1:3,iSFC))-rbin(ir) < 0d0 .or. MAXVAL(SFC_pos_glob(1:3,iSFC))+rbin(ir) > boxlen) then   ! merge_ks: boxlen (ks: 1d0)
        print '(A,1P4E15.6,/,A)', &
          '*** ERROR *** SFC may extend beyond the boundary (ToDo: implementation of periodic boundary)', &
          SFC_pos_glob(1:3,iSFC),rbin(ir)*scale_pc, &
          'stopping...'
        stop
      end if

      if (SFC_output_minmass < 0 .or. &
        mass*scale_msun >= SFC_output_minmass) then
        call output_around_SFC(iSFC, SFC_pos_glob(1:3,iSFC),use_current_ufront=.true.)      !output info about SFC
      end if


      !---------------------------------
      ! Add SFC to particle list
      !---------------------------------
      !increase total number of stellar paritlces
      nstar_tot =  nstar_tot + 1
      SFC_nSFC_loc = 0

      !get cpu responsible for pos
      ! pos(1,1:3) = SFC_pos_glob(1:3,iSFC)
      ! call cmp_cpumap(pos,cc,1)
      ! if (cc(1) .eq. myid)then

      !get cpu that has active grid contatining the given position
      pos1(1:3) = SFC_pos_glob(1:3,iSFC)
      !call get_cpu_of_pos(pos(1,1:3),ilevel,icpu)
      call get_active_cell_index_highest_level(pos1(1:3),cell_levl1,cell_index1)

      !print *, "KS DEBUG myid, pos1(1:3),cell_levl1,cell_index1", myid, pos1(1:3),cell_levl1,cell_index1

      ! if (myid==1) then
      !    print *, "KS DEBUG  (in star_formation_in_SFC result of cmp_cpumap)"
      !    print *, "pos(1,:)=", pos(1,:)
      !    print *, "cc(1) from cmp_cpumap=", cc(1)
      !    print *, "icpu from get_cpu_of_pos=", icpu
      ! end if


      !if (icpu .eq. myid)then
      if (cell_index1 >= 0) then           !when active grid in the highest level is found on this cpu
        !local number of SFC particle (0 or 1 for each SF via SFC)
        SFC_nSFC_loc = 1

        !get index of the cell where SFC particle is located
        !call get_cell_index(cell_index,cell_levl,pos,nlevelmax,1)
        !corresponding grid index
        ind_grid(1)=mod(cell_index1-ncoarse-1,ngridmax)+1 !cell -> grid where the cell belongs

        ! Update linked list for stars
        call remove_free(ind_part,1)
        call add_list(ind_part,ind_grid,ok_new,1)


        ! print *,"KS DEBUG (in star_formation_in_SFC)"
        ! print *, "cell_index(1)=",cell_index(1)
        ! print *, "myid=",myid
        ! print *, "ufront(cell_index(1),1)=",ufront(cell_index(1),1)
        ! print *, "uold(cell_index(1),1)=",uold(cell_index(1),1)
        ! print *, "unew(cell_index(1),1)=",unew(cell_index(1),1)
        ! print *, "pos(1,:)=", pos(1,:)


        ! Set star particle variables
        index_star=nstar_tot         ! particle index
        tp(ind_part(1))=birth_epoch  ! Birth epoch
        mp(ind_part(1))= 0      ! Mass
        levelp(ind_part(1))=cell_levl1   ! Level
        idp(ind_part(1))=index_star  ! Star identity (negative integers will be allocated to sinks)
        typep(ind_part(1))%family=FAM_SFC  ! Particle type
        xp(ind_part(1),1)=pos1(1)
        xp(ind_part(1),2)=pos1(2)
        xp(ind_part(1),3)=pos1(3)

        !set initial star velocity to gas velocity of the host cell
        ! use ufront
        vp(ind_part(1),1)=ufront(cell_index1,2)/ufront(cell_index1,1)
        vp(ind_part(1),2)=ufront(cell_index1,3)/ufront(cell_index1,1)
        vp(ind_part(1),3)=ufront(cell_index1,4)/ufront(cell_index1,1)

        !--------------------------------------------------------------------------
        !** WARNING **  WARNING **  WARNING **  WARNING **  WARNING **  WARNING ** !
        !--------------------------------------------------------------------------
        ! use metallicity container zp for SFC radius in order not to increase particle
        ! variables
        ! zp = minus radius (add minus sign to make clear the value is not metallicity)
        zp(ind_part(1)) = - radius ! SFC radius
#ifdef INIT_STELLAR_MASS
        mp0(ind_part(1))=mp(ind_part(1))
#endif
#ifdef STELLAR_POPULATION_MASS
        msp0(ind_part(1))=mp(ind_part(1))
#endif
        if(sf_log_properties) call sf_log_write_star(ilun_sf,ind_part(1),cell_index1,ilevel,(cell_levl1==nlevelmax))
        !--------------------------------------------------------------------------
        !** WARNING **  WARNING **  WARNING **  WARNING **  WARNING **  WARNING ** !
        !--------------------------------------------------------------------------

        write(*,'(A,2I4,I10,A,I8,A,I2,A,1PE8.1,A,1P3E9.1, A,1P3E16.8)') &
          "** New SFC test particle added **  myid, lev, ind_cell: ", myid, cell_levl1, cell_index1,&
          ", idp:",  idp(ind_part(1)),&
          ", fm:", typep(ind_part(1))%family,&
          ", radius:", -zp(ind_part(1))*scale_pc,&
          ", vp:", vp(ind_part(1),:)*scale_l/scale_t/1d5,&
          ", xp:", xp(ind_part(1),:)
      end if


      !---------------------------------
      !create global list for number of SFC particle (0 or 1 for each SF via SFC)
      !---------------------------------
      SFC_nSFC_list(:) = 0
      SFC_nSFC_list(myid) = SFC_nSFC_loc
      sbufi(1:ncpu) = SFC_nSFC_list(:)
      call MPI_ALLREDUCE(sbufi,rbufi,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
      SFC_nSFC_list(:) = rbufi(1:ncpu)


      !------------------------------------------
      ! output SFC info to stdout and a log file
      !------------------------------------------
      if (myid == 1) then
        print '(A,/,A,I5,A,I5,/,A,I5,6(A, 1P1E11.3),/, A,1P3E20.12,/)', '** SFC found **',&
          'iSFC/nSFC: ', iSFC, '/', SFC_n_tot, &
          'ir = ', ir, &
          ', R [pc] =', radius*scale_pc, &
          ', M [Msun]= ', mass*scale_msun,&
          ', nH_cen [cm-3] = ', RHOshell(1,iSFC)*scale_nH,&
          ', nH_ave [cm-3] = ', rho*scale_nH,&
          ', nH_bnd [cm-3] = ', nH_bnd,&
          ', Z [Zsun] = ', Zmetal/0.02,&
          'pos = ', SFC_pos_glob(1:3,iSFC)

        !trik to make the standard output easier to read
        call flush(6) !trick part1

        !also write a log file (logSFC_1Dprof; for quick check by eyes)
        open(LUN, file='logSFC_1Dprof', position='APPEND')
        write (LUN,'(A,2I10,1P1E14.6,I6)') 'coarse_step, fine_step, zred, iSFC:', nstep_coarse, nstep, 1./aexp - 1., iSFC
        write (LUN,'(5(A, 1P1E11.3))') &
          'R [pc] =', radius*scale_pc, &
          ', M [Msun]= ', mass*scale_msun,&
          ', nH_ave [cm-3] = ', rho*scale_nH,&
          ', nH_bnd [cm-3] = ', nH_bnd,&
          ', Z [Zsun] = ', Zmetal/0.02
        write (LUN,'(A,1P3E20.12)') 'SFC_pos(:) = ', SFC_pos_glob(1:3,iSFC)
        write (LUN,'(A,I5,A,I5)') 'Nrbin = ', Nrbin, ', ir = ', ir
        write (LUN,'(A,1P100E10.2)') 'rbin(:):    ', rbin(:)*scale_pc
        write (LUN,'(A,100I10)')     'Nshell(:):  ', Nshell(1:Nrbin,iSFC)
        write (LUN,'(A,1P100E10.2)') 'Mshell(:):  ', Mshell(1:Nrbin,iSFC)*scale_msun
        write (LUN,'(A,1P100E10.2)') 'Vshell(:):  ', Vshell(1:Nrbin,iSFC)*scale_pc**3
        write (LUN,'(A,1P100E10.2)') 'RHOshell(:):', RHOshell(1:Nrbin,iSFC)*scale_nH
        write (LUN,'(A,1P100E10.2)') 'MZshell(:): ', MZshell(1:Nrbin,iSFC)*scale_msun
        write (LUN,'(A,1P100E10.2)') 'T2shell(:): ', T2shell(1:Nrbin,iSFC)*scale_T2
        write (LUN,'(A)') ''
        close(LUN)

      end if

      call MPI_BARRIER(MPI_COMM_WORLD,info)  !trick part2


      !------------------------------------------
      ! convert SFC to stars (+ feedback on SFC)
      !------------------------------------------
      if (SFC_distribution=='mass_weight') then
        rho_SF  = SFC_f_rhoSF * RHOshell(1,iSFC) !threshold density for SF (limit SF to rho > rho_SF)
        call convert_SFC_to_stars(SFC_pos_glob(1:3,iSFC), radius, mass, rho, Zmetal, m_star_form, n_star_form, rho_SF, ok_flag, rbin, Menc, RHOenc, Zenc, Nrbin)
      else
        call convert_SFC_to_stars(SFC_pos_glob(1:3,iSFC), radius, mass, rho, Zmetal, m_star_form, n_star_form, ok_flag=ok_flag,&
          rbin=rbin, Menc=Menc, RHOenc=RHOenc, Zenc=Zenc, Nrbin=Nrbin)
      end if


      !------------------------------------------
      ! write logSFC (for later analysis)
      !------------------------------------------
      call write_log_SFC(iSFC, SFC_pos_glob(1:3,iSFC), radius, mass, rho, Zmetal, m_star_form, n_star_form)

    end do



    ! !KS DEBUG TENTATIVE
    ! do sweep_level = levelmin,nlevelmax
    !    do ind=1,twotondim
    !       iskip=ncoarse+(ind-1)*ngridmax
    !       do i=1,active(sweep_level)%ngrid
    !          if (myid == 5 .and. active(sweep_level)%igrid(i)+iskip == 4821) then
    !             print '(A,2I4,I8,1P3E12.4)', &
    !                  '(KS DEBUG TENTATIVE at the end of star_formation_in_SFC) myid, level, ind_cell, nH(uold), nH(unew), nH(ufront): ',&
    !                  myid, sweep_level, active(sweep_level)%igrid(i)+iskip, uold(active(sweep_level)%igrid(i)+iskip,1)*scale_nH, unew(active(sweep_level)%igrid(i)+iskip,1)*scale_nH,&
    !                  ufront(active(sweep_level)%igrid(i)+iskip,1)*scale_nH
    !          end if
    !       end do
    !    end do
    ! end do

  end subroutine star_formation_in_SFC


  !------------------------------------------------------------------------------------------------
  ! subroutine to instantaneously convert a SFC to cluster stars
  ! * feedback resulting from stellar formationon on the host SFC is also described
  !------------------------------------------------------------------------------------------------
  ! - stars are randomly distributed inside the cloud radius
  ! - use 10 % star formation efficiency for now, but later implement some formula based on
  !   small-scale simulations (KS TODO)
  ! - gas density inside the cloud is depleted and uniformized (tentative)
  !
  !------------------------------------------------------------------------------------------------
  subroutine convert_SFC_to_stars(pos, radius, mass, rho, Zmetal, m_star_form, n_star_form, rho_SF, ok_flag, rbin, Menc, RHOenc, Zenc, Nrbin)
    real(dp),dimension(1:3),intent(in) :: pos  !center
    real(dp),intent(inout) :: mass !SFC mass
    real(dp),intent(inout) :: radius !SFC radius
    real(dp),intent(inout) :: rho !SFC average density
    real(dp),intent(in) :: Zmetal !SFC average metallicity
    real(dp),intent(out) :: m_star_form !mass of cluster stars
    integer,intent(out) :: n_star_form !number of cluster stars
    real(dp),intent(in),optional :: rho_SF !threshold density for SF (limit SF to rho > rho_SF) (optional)

    !integer :: n_star_form, i_star
    integer :: i_star
    real(dp) :: f_star, f_star_form, rho1, runit, theta, phi

    integer,dimension(1:nvector)::cell_index,cell_levl,cc
    real(dp),dimension(1:3)::star_pos
    real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v, scale_msun, scale_pc

    integer,dimension(1:Nnew_max) ::ind_cell_new2,level_new
    real(dp),dimension(1:3,1:Nnew_max) :: pos_new

    !MPI send/recv buffer
    integer,dimension(1:ncpu) :: sbufi,rbufi

    !SFC check (SFC_mass < m_star_form)
    integer :: ntry, ir
    logical, intent(out),optional :: ok_flag ! False if we cannot make stars for a given rho_SF
    integer, intent(in), optional :: Nrbin
    real(dp),dimension(:),intent(in),optional :: rbin, Menc, RHOenc, Zenc

    integer :: nnew, SFC_nstar_tot
    real(dp) :: mnew

    integer :: ivec, nnew_tmp

    !star-by-star
    real(dp),dimension(1:Nnew_max) :: mstar_list, mstar_new
    real(dp) :: mass_cluster !star-cluster mass

    ok_flag = .false.

    ! Conversion factor from user units to cgs units
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    scale_msun=scale_d*scale_l**3/m_sun
    scale_pc=scale_l/3.08568d18

    !star formation efficiency (tentative, KS TODO)
    !f_star = 0.1d0
    !f_star = SFC_f_star
    ntry = 0
    do while(.not. ok_flag .and. ntry < 50)
      !use formula based on small-scale simulations if SFC_f_star is negative
      if (SFC_f_star < 0d0) then
        f_star = f_star_tot_HRG19(rho, mass, Zmetal) !star formation efficiency from RHD MC-scale simulations (He, Ricotti, Geen 2019); TO MODIFIED 2026-11-10: Zmetal for SFC_f_star_Zdep

      else !use given SFC_f_star
        f_star = SFC_f_star
      end if

      ! !**************** WARNING WARNING WARNING WARNING WARNING WARNING *****************!
      ! print *, "*** WARNING *** KS TENTATIVE DEBUG *** WARNING *** "
      ! print *, "f_star is set to 0.99 for test"
      ! print *, "*** WARNING *** KS TENTATIVE DEBUG *** WARNING *** "
      ! f_star = 0.99d0
      ! !**************** WARNING WARNING WARNING WARNING WARNING WARNING *****************!
      !f_star = some_function_of_rho(rho)


      !star-by-star model
      if (pop2_model == 5) then
        mass_cluster = mass*f_star
        call get_cluster_mstar_list(mass_cluster,Zmetal,n_star_form,mstar_list)
        m_star_form = SUM(mstar_list(1:n_star_form)) !stellar mass
        f_star_form = m_star_form / mass !actual efficiency

      else !non-star-by-star model

        ! # of stellar particles to be formed in SFC
        n_star_form = max(int(mass*f_star/mstar), 1) !form at least one star particle

        !actual efficiency
        f_star_form = n_star_form*mstar / mass

        !stellar mass
        m_star_form = n_star_form*mstar

      end if


      !KS DEBUG
      if (myid == 1) then
        print '(A,1P5E10.2,I10)', &
          '(convert_SFC_to_stars) M_SFC, nH_SFC, f_star, M_star_form, f_star_form, n_star_form: ', &
          mass*scale_msun, rho*scale_nH, f_star, m_star_form*scale_msun, f_star_form, n_star_form
      end if

      !check if nnew exceeds the limit
      if (n_star_form > Nnew_max) then
        if (myid == 1) then
          print '(A,2I8)', '(convert_SFC_to_stars) number of newly formed stars exceeds the maximum: ',n_star_form, Nnew_max
          print '(A)', 'stopping...'
          stop
        end if
      end if

      !check if M_star exceed M_SFC
      ! if (m_star_form > mass) then
      !   if (myid == 1) then
      !     print '(A,1PE12.5,I8)', '(convert_SFC_to_stars) mass of newly formed stars exceeds SFC mass: ', m_star_form, mass
      !     print '(A)', 'stopping...'
      !     stop
      !   end if
      ! end if
      ! 星の質量がSFCの質量を超えたらやり直しフラグを立てる
      if (m_star_form > mass) then
        ok_flag = .false.
      else
        ok_flag = .true.
      end if

      if (ok_flag) exit

      if (myid == 1) then
        print *, '---------------------------------------------------'
        print *, '[WARNING] SFC Star Formation Retry Triggered!'
        print *, '  m_star_form :', m_star_form
        print *, '  SFC mass    :', mass
        print *, '  Current rad :', radius
        print *, '  Retry count :', ntry + 1
        print *, '  Expanding SFC radius and recalculating...'
        print *, '---------------------------------------------------'
      end if

      ntry = ntry + 1
      ! --- 半径を拡張して新しいパラメータを取得する処理 ---
      do ir = 1, Nrbin
        if (rbin(ir) > radius * (1.0 + 1d-8)) then
          radius = rbin(ir)
          mass = Menc(ir)   ! ← 直接 mass を上書きする
          rho = RHOenc(ir)  ! ← 直接 rho を上書きする
          exit
        end if
      end do


    end do

    !----------------------
    ! add stars
    !----------------------
    !distribute stars in each cpu according to distribution model
    if (SFC_distribution=='mass_weight') then
      !star-by-star model
      if (pop2_model == 5) then
        call distribute_stars(pos, radius, n_star_form, nnew, ind_cell_new2, level_new, pos_new, rho_SF=rho_SF,mstar_list=mstar_list,mstar_new=mstar_new)

      else  !non-star-by-star model
        call distribute_stars(pos, radius, n_star_form, nnew, ind_cell_new2, level_new, pos_new, rho_SF=rho_SF)
      end if

    else
      call distribute_stars(pos, radius, n_star_form, nnew, ind_cell_new2, level_new, pos_new)
    end if

    !star-by-star model
    if (pop2_model == 5) then
      mnew = SUM(mstar_new(1:nnew)) !sum of new stellar mass in myid

    else !non-star-by-star model
      mnew = nnew * mstar !local stellar mass
    end if




    ! nnew = 0
    ! mnew = 0d0
    ! do i_star = 1, n_star_form
    !    !get random position in a sphere
    !    call random_pos_in_unit_sphere(runit,theta,phi)
    !    star_pos(1) = pos(1) + radius * runit * sin(theta) * cos(phi)
    !    star_pos(2) = pos(2) + radius * runit * sin(theta) * sin(phi)
    !    star_pos(3) = pos(3) + radius * runit * cos(theta)

    !    ! !KS DEBUG
    !    ! print *, 'radius, runit, theta, phi: ', radius, runit, theta, phi
    !    ! print *, "star_pos(:): ", star_pos(:)

    !    !get cpu responsible for pos
    !    call cmp_cpumap(star_pos,cc,1)
    !    if (cc(1) .eq. myid)then
    !       nnew = nnew + 1      !local number of star
    !       mnew = mnew + mstar  !local stellar mass

    !       !check if nnew excceeds the limit
    !       if (nnew > Nnew_max) then
    !          if (myid == 1) then
    !             print '(A,2I8)', 'number of newly formed stars exceeds the maximum in convert_SFC_to_stars: ',nnew, Nnew_max
    !             print '(A)', 'stopping...'
    !             stop
    !          end if
    !       end if

    !       !get index of the cell where sink particle is located
    !       call get_cell_index(cell_index,cell_levl,star_pos,nlevelmax,1)

    !       !corresponding grid index
    !       ind_grid(1)=mod(cell_index1-ncoarse-1,ngridmax)+1 !cell -> grid where the cell belongs

    !       ! !KS DEBUG
    !       ! print '(A,4I8, 1P3E20.12)', &
    !       !      "myid, ind_grid(1), cell_index1, cell_levl1, star_pos(1,:): ", &
    !       !      myid, ind_grid(1), cell_index1, cell_levl1, star_pos(1,:)

    !       ind_cell_new2(nnew) = cell_index1
    !       level_new(nnew) = cell_levl1
    !       pos_new(:,nnew) = star_pos(1,:)
    !    end if
    ! end do

    !---------------------------------
    ! Check for free particle memory
    !---------------------------------
    ok_free=(numbp_free-nnew)>=0
    call MPI_ALLREDUCE(numbp_free,numbp_free_tot,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
    if(.not. ok_free)then
      write(*,*)'No more free memory for particles'
      write(*,*)'Increase npartmax'
      call MPI_ABORT(MPI_COMM_WORLD,1,info)
    end if

    !---------------------------------
    ! Compute global star statistics
    !---------------------------------
    !create SFC_nstar_list (global)
    SFC_nstar_list(:) = 0
    SFC_nstar_list(myid) = nnew
    sbufi(1:ncpu) = SFC_nstar_list(:)
    call MPI_ALLREDUCE(sbufi,rbufi,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
    SFC_nstar_list(:) = rbufi(1:ncpu)
    SFC_nstar_tot = sum(SFC_nstar_list) !total number of stars in SFC

    !starting index of star at each cpu
    index_star=nstar_tot
    do icpu = 1, myid-1  ! "do i = 1,0" loop is just skipped (KS MEMO)
      index_star = index_star + SFC_nstar_list(icpu)
    end do

    !total number/mass of stars including new and old ones
    nstar_tot=nstar_tot+SFC_nstar_tot  !global total number
    mstar_tot = mstar_tot + mnew !I don't know why but this is local quantity

    ! !flag2 will be used only for star-formation log
    ! flag2(:) = 0 !initialization


    !---------------------------------
    ! Add stars to particle list
    !---------------------------------
    ncache=nnew
    do ivec=1,ncache,nvector
      nnew_tmp=MIN(nvector,ncache-ivec+1)
      !get grid index from cell index
      do i=1,nnew_tmp
        inew = ivec+i-1
        ind_grid_new(i)=mod(ind_cell_new2(inew)-ncoarse-1,ngridmax)+1 !cell -> grid where the cell belongs
      end do

      ! Update linked list for stars
      call remove_free(ind_part,nnew_tmp)
      call add_list(ind_part,ind_grid_new,ok_new,nnew_tmp)

      do i=1,nnew_tmp
        inew = ivec+i-1
        ! Set star particle variables
        index_star=index_star+1
        tp(ind_part(i))=birth_epoch  ! Birth epoch

        !star-by-star model
        if (pop2_model == 5) then
          mp(ind_part(i))= mstar_new(inew) ! Mass

        else !non-star-by-star model
          mp(ind_part(i))= mstar      ! Mass
        end if

        levelp(ind_part(i))=level_new(inew)   ! Level
        idp(ind_part(i))=index_star  ! Star identity (negative integers will be allocated to sinks)
        typep(ind_part(i))%family=FAM_STAR  ! Particle type
        xp(ind_part(i),1)=pos_new(1,inew)
        xp(ind_part(i),2)=pos_new(2,inew)
        xp(ind_part(i),3)=pos_new(3,inew)

        !set initial star velocity to gas velocity of the host cell
        ! use ufront
        vp(ind_part(i),1)=ufront(ind_cell_new2(inew),2)/ufront(ind_cell_new2(inew),1)
        vp(ind_part(i),2)=ufront(ind_cell_new2(inew),3)/ufront(ind_cell_new2(inew),1)
        vp(ind_part(i),3)=ufront(ind_cell_new2(inew),4)/ufront(ind_cell_new2(inew),1)

        ! ! uold for velocity was already divided by density at nlevelmax
        ! if (level_new(inew) == nlevelmax) then
        !    vp(ind_part(i),1)=uold(ind_cell_new2(inew),2)
        !    vp(ind_part(i),2)=uold(ind_cell_new2(inew),3)
        !    vp(ind_part(i),3)=uold(ind_cell_new2(inew),4)
        ! else
        !    vp(ind_part(i),1)=uold(ind_cell_new2(inew),2)/uold(ind_cell_new2(inew),1)
        !    vp(ind_part(i),2)=uold(ind_cell_new2(inew),3)/uold(ind_cell_new2(inew),1)
        !    vp(ind_part(i),3)=uold(ind_cell_new2(inew),4)/uold(ind_cell_new2(inew),1)
        ! end if


        ! print '(A,5I6)', "(KS DEBUG) myid, inew, idp, ipart, ind_cell: ", myid, inew, index_star, ind_part(i), ind_cell_new2(inew)

        zp(ind_part(i))= Zmetal !use average metalicity of SFC instead of cell metallicity
#ifdef INIT_STELLAR_MASS
        mp0(ind_part(i))=mp(ind_part(i))
#endif
#ifdef STELLAR_POPULATION_MASS
        msp0(ind_part(i))=mp(ind_part(i))
#endif
        if(sf_log_properties) call sf_log_write_star(ilun_sf,ind_part(i),ind_cell_new2(inew),ilevel,(level_new(inew)==nlevelmax))


        ! !used only for writing star-formation log (the part of No adding star in the original code)
        ! !reset flag2 to 0 again in write_log_star_formation(isSFC = .true.)
        ! flag2(ind_cell_new2(inew)) = flag2(ind_cell_new2(inew)) + 1 !more than two stars can form in a single cell

        !------------------------- KS DEBUG ----------------------------!
        ! write(*,'(A,2I4,I10,A,1PE9.2,A,I9,A,I3,A,1PE9.2, A,1P3E16.8)') &
        !      "** New SFC star added **  myid, lev, ind_cell: ", myid, level_new(inew), ind_cell_new2(inew),&
        !      ", mp [Msun]:", mp(ind_part(i))*scale_msun,&
        !      ", idp:",         idp(ind_part(i)),&
        !      ", typep%family:", typep(ind_part(i))%family,&
        !      ", zp [Zsun]:", zp(ind_part(i))/0.02,&
        !      ", xp:", xp(ind_part(i),:)
        write(*,'(A,2I4,I10,A,1PE10.3,A,I8,A,I2,A,1PE8.1,A,1P3E9.1, A,1P3E16.8)') &
          "** New SFC star added **  myid, lev, ind_cell: ", myid, level_new(inew), ind_cell_new2(inew),&
          ", mp:", mp(ind_part(i))*scale_msun,&
          ", idp:",  idp(ind_part(i)),&
          ", fm:", typep(ind_part(i))%family,&
          ", zp:", zp(ind_part(i))/0.02,&
          ", vp:", vp(ind_part(i),:)*scale_l/scale_t/1d5,&
          ", xp:", xp(ind_part(i),:)
        !------------------------- KS DEBUG ----------------------------!
      end do
    end do

    !----------------------
    ! take a log for formed stars
    !----------------------

    !----------------------
    ! backreaction to the SFC gas by SF
    !----------------------
    if (SFC_fate == 'homog') then !homogenize the SFC by setting density to the value that ensures mass conservation
      !rho after star formation (assume uniform density)
      rho1 = rho * (1-f_star_form)
      call homogenize_SFC(pos,radius,rho1)
    else if (SFC_fate == 'crop') then !crop high-density cells and then adjust density to ensure mass conservation
      call crop_and_adjust_SFC(pos,radius,mass,m_star_form,rho)
    else
      if (myid == 1) then
        print '(A,A,A)', 'SFC_fate = ', trim(SFC_fate), ' but should be either homegenize or crop'
        print '(A)', 'stopping...'
        stop
      end if
    end if

  end subroutine convert_SFC_to_stars

  !-------------------------------------------------------------------------
  ! distribute SFC stars
  ! 3 modes are currently implemented
  !  1. random distribution in a sphere
  !  2. volume-weighted distribution in the region with rho > rho_sf  <-- not implemented
  !  3. mass-weighted random distribution in the region with rho > rho_sf
  !-------------------------------------------------------------------------
  subroutine distribute_stars(pos, radius, n_star_form, nnew, ind_cell_new2, level_new, pos_new, rho_SF, mstar_list,mstar_new)
    !IN/OUT
    real(dp),dimension(1:3),intent(in) :: pos  !center
    real(dp),intent(in) :: radius !SFC radius
    !real(dp),intent(in) :: rho !SFC density
    integer,intent(in) :: n_star_form !number of stars in the SFC
    integer,intent(out) :: nnew  !number of stars in myid
    integer,dimension(1:Nnew_max),intent(out) ::ind_cell_new2,level_new !ind_cell/level of new stars
    real(dp),dimension(1:3,1:Nnew_max),intent(out) :: pos_new !pos of new stars
    real(dp),intent(in),optional :: rho_SF !threshold density for SF (limit SF to rho > rho_SF) (optional)
    real(dp),dimension(1:Nnew_max),intent(in),optional :: mstar_list !global list of stellar mass (for star-by-star model)
    real(dp),dimension(1:Nnew_max),intent(out),optional :: mstar_new  !list of stellar mass in myid (for star-by-star model)

    !variables
    real(dp) :: vol_cell, dx_cell
    integer :: sweep_level !search for all leaf cells through levels (not only ilevel)
    real(dp) :: x, y, z, r2, r
    real(dp),dimension(1:twotondim,1:3)::xc_level
    integer :: ix, iy, iz

    real(dp) :: f_star, f_star_form, m_star_form, rho1, runit, theta, phi

    integer,dimension(1:nvector)::cell_index,cell_levl,cc
    integer::cell_index1,cell_levl1
    real(dp),dimension(1:3)::star_pos
    real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v, scale_msun, scale_pc

    real(dp),dimension(1:ncpu) :: weight_list_tmp, weight_list
    real(dp) :: cell_weight_min, cell_weight_max, weight_loc, my_weight_min, weight_tot
    real(dp),dimension(1:Nnew_max) :: weight_rand
    real(dp) :: rand !, rho_SF
    integer :: i_star
    integer, parameter :: PRIMARY_RANK = 1


    !consistency check for optional arguments (for star-by-star model)
    if (present(mstar_new) .and. (.not. present(mstar_list))) then
      print *, "*** ERROR in distribute_stars *** mstar_list must be provided to return mstar_new"
    end if

    ! Conversion factor from user units to cgs units
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    scale_msun=scale_d*scale_l**3/m_sun
    scale_pc=scale_l/3.08568d18

    nnew = 0
    !-------------------------------------
    ! random distribution in a sphere
    !-------------------------------------
    if (SFC_distribution=='spherical') then
      if (myid == 1) then
        print '(A)', 'place SFC stars randomly in the sphere with SFC radius'
      end if
      do i_star = 1, n_star_form
        !get random position in a sphere
        call random_pos_in_unit_sphere(runit,theta,phi)
        star_pos(1) = pos(1) + radius * runit * sin(theta) * cos(phi)
        star_pos(2) = pos(2) + radius * runit * sin(theta) * sin(phi)
        star_pos(3) = pos(3) + radius * runit * cos(theta)

        ! !KS DEBUG
        ! print *, 'radius, runit, theta, phi: ', radius, runit, theta, phi
        ! print *, "star_pos(:): ", star_pos(:)

        !get cpu responsible for pos
        call get_active_cell_index_highest_level(star_pos,cell_levl1,cell_index1)
        !call cmp_cpumap(star_pos,cc,1)
        if (cell_index1 >= 0) then           !when active grid in the highest level is found on this cpu
          !if (cc(1) .eq. myid)then
          nnew = nnew + 1      !local number of star

          ! !check if nnew excceeds the limit
          ! if (nnew > Nnew_max) then
          !    if (myid == 1) then
          !       print '(A,2I8)', 'number of newly formed stars exceeds the maximum in distribute_stars: ',nnew, Nnew_max
          !       print '(A)', 'stopping...'
          !       stop
          !    end if
          ! end if

          !get index of the cell where sink particle is located
          !call get_cell_index(cell_index,cell_levl,star_pos,nlevelmax,1)

          !corresponding grid index
          ind_grid(1)=mod(cell_index1-ncoarse-1,ngridmax)+1 !cell -> grid where the cell belongs

          ! !KS DEBUG
          ! print '(A,4I8, 1P3E20.12)', &
          !      "myid, ind_grid(1), cell_index1, cell_levl1, star_pos(1,:): ", &
          !      myid, ind_grid(1), cell_index1, cell_levl1, star_pos(1,:)

          ind_cell_new2(nnew) = cell_index1
          level_new(nnew) = cell_levl1
          pos_new(:,nnew) = star_pos(:)
        end if
      end do

      !----------------------------------------------
      ! volume/mass weighted random distribution
      !----------------------------------------------
    else if (SFC_distribution=='mass_weight') then
      !threshold density for SF (limit SF to rho > rho_SF)
      !rho_SF = SFC_f_rhoSF * rho

      !-------- count weight in each cpu -------!
      weight_loc = 0d0
      ! Loop over levels
      do sweep_level = levelmin,nlevelmax
        dx_cell = 1d0/dble(2**sweep_level)
        vol_cell = (dx_cell*scale)**3   ! merge_ks: cell volume in code units for boxlen /= 1 (ks assumed boxlen=1) !cell volume at sweep_level

        ! Cells center position relative to grid center position
        do ind=1,twotondim
          iz=(ind-1)/4
          iy=(ind-1-4*iz)/2
          ix=(ind-1-2*iy-4*iz)
          xc_level(ind,1)=(dble(ix)-0.5D0)*dx_cell
          xc_level(ind,2)=(dble(iy)-0.5D0)*dx_cell
          xc_level(ind,3)=(dble(iz)-0.5D0)*dx_cell
        end do

        ! Loop over grids
        ncache=active(sweep_level)%ngrid
        do igrid=1,ncache,nvector
          ngrid=MIN(nvector,ncache-igrid+1)
          do i=1,ngrid
            ind_grid(i)=active(sweep_level)%igrid(igrid+i-1)
          end do
          ! Loop over cells
          do ind=1,twotondim
            iskip=ncoarse+(ind-1)*ngridmax
            do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
            end do

            ! acuumulate local weight
            do i=1,ngrid
              !only the leaf cells
              if(son(ind_cell(i))/=0) cycle

              !radius from pos
              x = (xg(ind_grid(i),1)+xc_level(ind,1)-skip_loc(1))*scale          !x
              y = (xg(ind_grid(i),2)+xc_level(ind,2)-skip_loc(2))*scale          !y
              z = (xg(ind_grid(i),3)+xc_level(ind,3)-skip_loc(3))*scale          !z
              r2 = (x - pos(1))**2 +  (y - pos(2))**2 +  (z - pos(3))**2
              r = sqrt(r2)

              !limit to a region where rho >  rho_SF and r < R_SFC
              if (ufront(ind_cell(i),1) > rho_SF .and. r < radius) then
                weight_loc = weight_loc + vol_cell * ufront(ind_cell(i),1)
              end if
              ! if (uold(ind_cell(i),1) > rho_SF .and. r < radius) then
              !    weight_loc = weight_loc + vol_cell * uold(ind_cell(i),1)
              ! end if
              ! if (sweep_level == nlevelmax) then
              !    !at maximum level
              !    if (uold(ind_cell(i),1) > rho_SF .and. r < radius) then
              !       weight_loc = weight_loc + vol_cell * uold(ind_cell(i),1)
              !    end if
              ! else
              !    !at lower levels
              !    if (unew(ind_cell(i),1) > rho_SF .and. r < radius) then
              !       weight_loc = weight_loc + vol_cell * unew(ind_cell(i),1)
              !    end if
              ! end if
            end do
          end do
          ! End loop over cells
        end do
        ! End loop over grids
      end do
      ! End loop over levels

      !-------- get global info about weight -------!
      !global list of weight
      weight_list_tmp(1:ncpu) = 0d0
      weight_list_tmp(myid) = weight_loc
      call MPI_ALLREDUCE(weight_list_tmp,weight_list,ncpu,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
      weight_tot = sum(weight_list) !total weight
      !min value of weight in myid
      my_weight_min = 0d0
      do icpu = 1, myid-1  ! "do i = 1,0" loop is just skipped (KS MEMO)
        my_weight_min = my_weight_min + weight_list(icpu)
      end do

      if (myid == 1) then
        ! print '(A,/,A, 1P4E12.4)', 'place SFC stars according to mass-weighted random distribution',&
        !      'SFC_f_rhoSF, nH_ave, nH_SF, M(nH>nH_SF) = ', SFC_f_rhoSF, rho*scale_nH, rho_SF*scale_nH, weight_tot*scale_msun
        print '(A,/,A, 1P3E12.4)', 'place SFC stars according to mass-weighted random distribution',&
          'SFC_f_rhoSF, nH_SF, M(nH>nH_SF) = ', SFC_f_rhoSF, rho_SF*scale_nH, weight_tot*scale_msun


        !-------------------------------------------------------------------------------------
        !I don't think the following error check is not reasonable, so I comment it out
        !-------------------------------------------------------------------------------------
        ! if (weight_tot < n_star_form*mstar) then
        !    print '(A,1P2E12.4)','(M(nH>nH_SF) < m_star_form', weight_tot*scale_msun, n_star_form*mstar*scale_msun
        !    print *, 'stopping...'
        !    stop
        ! end if
        !-------------------------------------------------------------------------------------
        !-------------------------------------------------------------------------------------
      end if

      !check the value of weight_tot
      if (weight_tot .le. 0) then
        if (myid == 1) then
          print '(A,1P1E12.4,A,/,A) ', "*** ERROR *** weight_tot = ", weight_tot, " is non-positive in distribute_stars", &
            'stopping...'
        end if
        stop
      end if


      !----- assign random value in [0, weight_tot] to each star -----!
      if (myid == PRIMARY_RANK) then
        do i_star = 1, n_star_form
          !call random_value_0to1(rand)  !get random number between 0 and 1
          call random_value_0to1(rand, initial_seed=nstep)  !get random number between 0 and 1 (pass nstep for initial_seed to avoid the same pattern)
          weight_rand(i_star) = weight_tot * rand
          !print *, "(KS DEBUG) i_star, weight_rand(i_star) = ",  i_star, weight_rand(i_star)
        end do
      end if



      call mpi_bcast(weight_rand, n_star_form, MPI_DOUBLE_PRECISION, PRIMARY_RANK-1, MPI_COMM_WORLD, info) !one smaller value from id shuould be set in mpi_bcast

      !-------- distribute stars according to the weight values  -------!
      weight_loc = 0d0
      ! Loop over levels
      do sweep_level = levelmin,nlevelmax
        dx_cell = 1d0/dble(2**sweep_level)
        vol_cell = (dx_cell*scale)**3   ! merge_ks: cell volume in code units for boxlen /= 1 (ks assumed boxlen=1) !cell volume at sweep_level

        ! Cells center position relative to grid center position
        do ind=1,twotondim
          iz=(ind-1)/4
          iy=(ind-1-4*iz)/2
          ix=(ind-1-2*iy-4*iz)
          xc_level(ind,1)=(dble(ix)-0.5D0)*dx_cell
          xc_level(ind,2)=(dble(iy)-0.5D0)*dx_cell
          xc_level(ind,3)=(dble(iz)-0.5D0)*dx_cell
        end do

        ! Loop over grids
        ncache=active(sweep_level)%ngrid
        do igrid=1,ncache,nvector
          ngrid=MIN(nvector,ncache-igrid+1)
          do i=1,ngrid
            ind_grid(i)=active(sweep_level)%igrid(igrid+i-1)
          end do
          ! Loop over cells
          do ind=1,twotondim
            iskip=ncoarse+(ind-1)*ngridmax
            do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
            end do

            ! acuumulate local weight
            do i=1,ngrid
              !only the leaf cells
              if(son(ind_cell(i))/=0) cycle

              !radius from pos
              x = (xg(ind_grid(i),1)+xc_level(ind,1)-skip_loc(1))*scale          !x
              y = (xg(ind_grid(i),2)+xc_level(ind,2)-skip_loc(2))*scale          !y
              z = (xg(ind_grid(i),3)+xc_level(ind,3)-skip_loc(3))*scale          !z
              r2 = (x - pos(1))**2 +  (y - pos(2))**2 +  (z - pos(3))**2
              r = sqrt(r2)

              !** acuumulate local weight and get cell min/max of weight **
              !each cell is responsible for the weight values in [cell_weight_min, cell_weight_max]
              cell_weight_min = weight_loc + my_weight_min
              !limit to a region where rho >  rho_SF and r < R_SFC
              if (ufront(ind_cell(i),1) > rho_SF .and. r < radius) then
                weight_loc = weight_loc + vol_cell * ufront(ind_cell(i),1)
              end if
              ! if (uold(ind_cell(i),1) > rho_SF .and. r < radius) then
              !    weight_loc = weight_loc + vol_cell * uold(ind_cell(i),1)
              ! end if
              ! if (sweep_level == nlevelmax) then  !at maximum level
              !    if (uold(ind_cell(i),1) > rho_SF .and. r < radius) then
              !       weight_loc = weight_loc + vol_cell * uold(ind_cell(i),1)
              !    end if
              ! else                                  !at lower levels
              !    if (unew(ind_cell(i),1) > rho_SF .and. r < radius) then
              !       weight_loc = weight_loc + vol_cell * unew(ind_cell(i),1)
              !    end if
              ! end if
              cell_weight_max = weight_loc + my_weight_min

              !loop over stars
              do i_star = 1, n_star_form
                !check if the star is in the cell
                if (cell_weight_min < weight_rand(i_star) .and. weight_rand(i_star) <= cell_weight_max) then
                  nnew = nnew + 1      !count number of new stars (local)

                  !store information about a newly formed star
                  ind_cell_new2(nnew) = ind_cell(i)
                  level_new(nnew) = sweep_level
                  pos_new(:,nnew) = (/x,y,z/)

                  !store stellar mass (for star-by-star model)
                  if (present(mstar_new)) then
                    mstar_new(nnew) = mstar_list(i_star)
                    ! print '(3(A,I8), 1(A,1PE12.4))', "KS DEBUG (distribute_stars) i_star = ", i_star, ", myid = ", myid, ", nnew = ", nnew, &
                    !      ", mstar_new[Msun] = ", mstar_new(nnew)*scale_msun
                  end if

                  ! !KS DEBUG
                  ! print *, '(KS DEBUG) ** cell found ** myid, i_star, ind_cell, weight_rand, cell_weight_min/max', &
                  !      myid, i_star, ind_cell(i), weight_rand(i_star), cell_weight_min, cell_weight_max

                end if
              end do
            end do
          end do
          ! End loop over cells
        end do
        ! End loop over grids
      end do
      ! End loop over levels
    else
      if (myid == 1) then
        print *, 'SFC_distribution = ', SFC_distribution, ' is yet to be implemented..'
        print *, 'stopping...'
        stop
      end if
    end if

  end subroutine distribute_stars



  !-------------------------------------------------------------------------
  ! random position in a unit sphere
  ! assume homogeneous distribution with respect to volume
  ! i.e., each dV = r**2 dr sinth dth d(phi) = 1/3 d(r**3) d(costh) d(phi)
  ! has the equal probability to be picked
  !-------------------------------------------------------------------------
  subroutine random_pos_in_unit_sphere(r,theta,phi)
    real(dp), intent(OUT) :: r, theta, phi
    real(dp) :: rand1,rand2,rand3
    integer i, seedsize
    integer,save,allocatable:: seed(:)
    integer, parameter :: initial_seed=123423424
    real(dp),parameter :: pi = 3.1415926d0
    real(dp) :: costh, r3


    if (.not. allocated(seed)) then
      call random_seed(size=seedsize)  ! シードの格納に必要なサイズを取得する
      allocate(seed(seedsize))         ! シード格納領域を確保
      seed(:) = initial_seed
    end if

    call random_seed(put=seed)         ! シードを渡す

    !random number between 0 and 1
    call random_number(rand1)
    call random_number(rand2)
    call random_number(rand3)


    ! matrix generating rotaion around z-axis
    r3 = rand1
    r = r3**(1.d0/3)
    costh = 2.d0*rand2 - 1.d0
    theta = acos(costh)
    phi = 2.d0*pi*rand3

    ! !KS DEBUG
    ! print *, '########################################################## '
    ! print *, '(KS DEBUG) rand1, 2, 3 = ', rand1, rand2, rand3
    ! print *, '########################################################## '

    call random_seed(get=seed) !現在のseedを取得、次回再利用
    !write (*, *) 'seed (after): ', seed

  end subroutine random_pos_in_unit_sphere

  !-------------------------------------------------------------------------
  ! random number in [0, 1]
  !-------------------------------------------------------------------------
  subroutine random_value_0to1(rand, initial_seed)
    real(dp), intent(OUT) :: rand
    integer i, seedsize
    integer,save,allocatable:: seed(:)
    integer, intent(iN), optional :: initial_seed
    integer, parameter :: initial_seed_default=61234234
    integer, parameter :: n_burn = 10000
    integer,save :: ifirst = 0 !for first call


    !if (.not. allocated(seed)) then
    !---------------------------
    ! first call
    !---------------------------
    if (ifirst == 0) then
      call random_seed(size=seedsize)  ! シードの格納に必要なサイズを取得する
      allocate(seed(seedsize))         ! シード格納領域を確保

      if (present(initial_seed)) then
        seed(:) = initial_seed
      else
        seed(:) = initial_seed_default
      end if

      call random_seed(put=seed)         ! シードを渡す

      !print *, "(random_value_0to1) seed (before burning) =", seed(:)

      !initial burnning phase
      do i = 1, n_burn
        call random_number(rand)
        !print *, "(random_value_0to1) initial burning phase: i, rand =", i, rand
      end do

      call random_seed(get=seed) !現在のseedを取得、次回再利用

      !print *, "(random_value_0to1) seed (after burning) =", seed(:)

      ifirst = 1
    end if

    !print *, "(random_value_0to1) seed (current) =", seed(:)

    call random_seed(put=seed)         ! シードを渡す

    !random number between 0 and 1
    call random_number(rand)
    !print *, "(random_value_0to1) rand =",  rand

    call random_seed(get=seed) !現在のseedを取得、次回再利用
    !write (*, *) 'seed (after): ', seed
  end subroutine random_value_0to1
  !----------------------------------------------------------------------------
  ! homogenize star-forming cloud in consequence of  star formation
  !----------------------------------------------------------------------------
  subroutine homogenize_SFC(pos,radius,rho)
    !in
    real(dp),dimension(1:3),intent(in) :: pos  !center
    real(dp),intent(in) :: radius !radius of depleting region
    real(dp),intent(in) :: rho !new density

    !variables
    integer :: sweep_level !search for all leaf cells through levels (not only ilevel)
    real(dp) :: x, y, z, r2, r
    real(dp),dimension(1:twotondim,1:3)::xc_level
    integer :: ix, iy, iz
    real(dp) :: dx_cell
    real(dp) :: f_rho

    !###################################################################!
    !   OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE  !
    !###################################################################!
    if (myid == 1) then
      print *, '** ERROR ** homogenize_SFC (SFC_fate = homog) is no longer supported'
      print *, 'The code needs some modification to use this method.'
      print *, 'stopping...'
      stop
    end if
    !###################################################################!
    !   OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE  !
    !###################################################################!


    ! !KS DEBUG
    ! if (myid == 1) then
    !    print '(A,1P3E10.2)', &
    !         '(homogenize_SFC)  nH [cm-3] -> ',&
    !         rho*scale_nH
    ! end if

    ! do sweep_level = levelmin,nlevelmax
    !    dx_cell = 1d0/dble(2**sweep_level)

    !    ! Cells center position relative to grid center position
    !    do ind=1,twotondim
    !       iz=(ind-1)/4
    !       iy=(ind-1-4*iz)/2
    !       ix=(ind-1-2*iy-4*iz)
    !       xc_level(ind,1)=(dble(ix)-0.5D0)*dx_cell
    !       xc_level(ind,2)=(dble(iy)-0.5D0)*dx_cell
    !       xc_level(ind,3)=(dble(iz)-0.5D0)*dx_cell
    !    end do


    !    ! sum shell qunatities in each cpu
    !    ncache=active(sweep_level)%ngrid
    !    do igrid=1,ncache,nvector
    !       ngrid=MIN(nvector,ncache-igrid+1)
    !       do i=1,ngrid
    !          ind_grid(i)=active(sweep_level)%igrid(igrid+i-1)
    !       end do
    !       ! Loop over cells
    !       do ind=1,twotondim
    !          iskip=ncoarse+(ind-1)*ngridmax
    !          do i=1,ngrid
    !             ind_cell(i)=iskip+ind_grid(i)
    !          end do

    !          !store information about newly formed stars and their parent cells
    !          do i=1,ngrid

    !             !depleteion only in the leaf cells (OK? KS TODO)
    !             if(son(ind_cell(i))/=0) cycle

    !             x = (xg(ind_grid(i),1)+xc_level(ind,1)-skip_loc(1))*scale          !x
    !             y = (xg(ind_grid(i),2)+xc_level(ind,2)-skip_loc(2))*scale          !y
    !             z = (xg(ind_grid(i),3)+xc_level(ind,3)-skip_loc(3))*scale          !z

    !             r2 = (x - pos(1))**2 +  (y - pos(2))**2 +  (z - pos(3))**2
    !             r = sqrt(r2)

    !             !change gas inside the sphere
    !             ! - at sweep_level = ilevel (= nlevelmax), all uold was converted to primitive variables,
    !             !   and will be converted back to conserved variables later
    !             if (r < radius) then
    !                !at ilevel (= nlevelmax)
    !                if (sweep_level == nlevelmax) then
    !                   uold(ind_cell(i),1) = rho        !change density

    !                else !at lower levels
    !                   !ratio of new and original density
    !                   f_rho = rho/uold(ind_cell(i),1)

    !                   !change all conservative vars
    !                   do ivar=1,nvar
    !                      !updat unew
    !                      unew(ind_cell(i),ivar) = unew(ind_cell(i),ivar) + (f_rho-1d0) * uold(ind_cell(i),ivar)
    !                      !update uold
    !                      uold(ind_cell(i),ivar) = f_rho * uold(ind_cell(i),ivar)
    !                   end do
    !                end if
    !                ! if (sweep_level == nlevelmax) then
    !                !    uold(ind_cell(i),1) = rho !change density
    !                ! else
    !                !    f_rho = rho/uold(ind_cell(i),1) !ratio of new and original density
    !                !    uold(ind_cell(i),1) = f_rho * uold(ind_cell(i),1)        !change density
    !                !    uold(ind_cell(i),2:5) = f_rho * uold(ind_cell(i),2:5)    !change momentum/internal energy
    !                !    do ivar=imetal,nvar
    !                !       uold(ind_cell(i),ivar) = f_rho * uold(ind_cell(i),ivar) !change chemical components and ivar_refine
    !                !    end do
    !                ! end if


    !             end if
    !          end do
    !       end do
    !       ! End loop over cells
    !    end do
    !    ! End loop over grids
    ! end do
    ! ! End loop over levels

  end subroutine homogenize_SFC

  !----------------------------------------------------------------------------
  ! crop and adjust star-forming cloud in consequence of  star formation
  !----------------------------------------------------------------------------
  subroutine crop_and_adjust_SFC(pos,radius,mass_SFC,mass_star,rho_SFC)
    !in
    real(dp),dimension(1:3),intent(in) :: pos  !center
    real(dp),intent(in) :: radius !radius of SFC
    real(dp),intent(in) :: mass_SFC !mass of SFC
    real(dp),intent(in) :: mass_star !mass of newly formed stars
    real(dp),intent(in) :: rho_SFC !average density of SFC

    !variables
    integer :: sweep_level !search for all leaf cells through levels (not only ilevel)
    real(dp) :: x, y, z, r2, r
    real(dp),dimension(1:twotondim,1:3)::xc_level
    integer :: ix, iy, iz
    real(dp) :: dx_cell, vol_cell
    real(dp) :: f_crop, f_adjust
    real(dp) :: mass_crop_loc, mass_crop_glob, mass_crop !mass cropped from high-density cells
    real(dp) :: rho_crop !cells with rho > rho_crop will be cropped to rho_crop

    !parameters
    !real(dp),parameter :: f_crop = 2d0 !rho_crop = f_crop * rho_SFC (MAYBE NEED TRY-AND-ERROR THIS VALUE)


    !----------------------------------
    ! crop high-density regions
    !----------------------------------
    rho_crop = SFC_f_rhocrop * rho_SFC  !threshold density for crop (average density times constant factor)
    mass_crop_loc = 0d0
    do sweep_level = levelmin,nlevelmax
      dx_cell = 1d0/dble(2**sweep_level)
      vol_cell = (dx_cell*scale)**3   ! merge_ks: cell volume in code units for boxlen /= 1 (ks assumed boxlen=1) !cell volume at sweep_level

      ! Cells center position relative to grid center position
      do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc_level(ind,1)=(dble(ix)-0.5D0)*dx_cell
        xc_level(ind,2)=(dble(iy)-0.5D0)*dx_cell
        xc_level(ind,3)=(dble(iz)-0.5D0)*dx_cell
      end do


      ! sum shell qunatities in each cpu
      ncache=active(sweep_level)%ngrid
      do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
          ind_grid(i)=active(sweep_level)%igrid(igrid+i-1)
        end do
        ! Loop over cells
        do ind=1,twotondim
          iskip=ncoarse+(ind-1)*ngridmax
          do i=1,ngrid
            ind_cell(i)=iskip+ind_grid(i)
          end do

          !store information about newly formed stars and their parent cells
          do i=1,ngrid
            !only the leaf cells
            if(son(ind_cell(i))/=0) cycle

            !radius from pos
            x = (xg(ind_grid(i),1)+xc_level(ind,1)-skip_loc(1))*scale          !x
            y = (xg(ind_grid(i),2)+xc_level(ind,2)-skip_loc(2))*scale          !y
            z = (xg(ind_grid(i),3)+xc_level(ind,3)-skip_loc(3))*scale          !z
            r2 = (x - pos(1))**2 +  (y - pos(2))**2 +  (z - pos(3))**2
            r = sqrt(r2)

            !crop if density exceeds rho_crop and r is within radius (R_SFC)
            if (ufront(ind_cell(i),1) > rho_crop .and. r < radius) then
              mass_crop_loc = mass_crop_loc + vol_cell * (ufront(ind_cell(i),1) - rho_crop)

              ! !KS DEBUG TENTATIVE
              ! if (myid == 5 .and. ind_cell(i) == 4821) then
              ! print '(A,2I4,I8,1P4E12.4)', &
              !      '(KS DEBUG TENTATIVE 1) myid, level, ind_cell, nH(uold), nH(unew), nH(ufront), f_crop: ',&
              !      myid, sweep_level, ind_cell(i), uold(ind_cell(i),1)*scale_nH, unew(ind_cell(i),1)*scale_nH,&
              !      ufront(ind_cell(i),1)*scale_nH, f_crop
              ! end if


              !update ufront
              f_crop = rho_crop/ufront(ind_cell(i),1) !ratio of cropped and original density
              do ivar=1,nvar
                ufront(ind_cell(i),ivar) = f_crop * ufront(ind_cell(i),ivar)   !change all conservative vars
              end do

              !update uold/unew
              if (sweep_level == nlevelmax) then                    !at ilevel (= nlevelmax)
                !change only density (others are automatically updated)
                uold(ind_cell(i),1) = ufront(ind_cell(i),1)
              else                                                  !at lower levels
                !change all conservative vars
                do ivar=1,nvar
                  !update uold (to be used for calculating flux by godunov method)
                  uold(ind_cell(i),ivar) = ufront(ind_cell(i),ivar)
                  !updat unew (variation of ufront used to ensure conservation)
                  unew(ind_cell(i),ivar) = unew(ind_cell(i),ivar) - (1d0-f_crop)/f_crop * ufront(ind_cell(i),ivar)
                end do
              end if

              ! !KS DEBUG TENTATIVE
              ! if (myid == 5 .and. ind_cell(i) == 4821) then
              ! print '(A,2I4,I8,1P4E12.4)', &
              !      '(KS DEBUG TENTATIVE 2) myid, level, ind_cell, nH(uold), nH(unew), nH(ufront), f_crop: ',&
              !      myid, sweep_level, ind_cell(i), uold(ind_cell(i),1)*scale_nH, unew(ind_cell(i),1)*scale_nH,&
              !      ufront(ind_cell(i),1)*scale_nH, f_crop
              ! end if
            end if
          end do
        end do
        ! End loop over cells
      end do
      ! End loop over grids
    end do
    ! End loop over levels

    !MPI communication
    call MPI_ALLREDUCE(mass_crop_loc, mass_crop_glob, 1, &
      MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)

    !cropped mass (global)
    mass_crop = mass_crop_glob

    !---------------------------------------------
    ! adjust density to ensure mass conservation
    !---------------------------------------------
    f_adjust = (mass_SFC - mass_star)/(mass_SFC - mass_crop)
    do sweep_level = levelmin,nlevelmax
      dx_cell = 1d0/dble(2**sweep_level)

      ! Cells center position relative to grid center position
      do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc_level(ind,1)=(dble(ix)-0.5D0)*dx_cell
        xc_level(ind,2)=(dble(iy)-0.5D0)*dx_cell
        xc_level(ind,3)=(dble(iz)-0.5D0)*dx_cell
      end do

      ! sum shell qunatities in each cpu
      ncache=active(sweep_level)%ngrid
      do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
          ind_grid(i)=active(sweep_level)%igrid(igrid+i-1)
        end do
        ! Loop over cells
        do ind=1,twotondim
          iskip=ncoarse+(ind-1)*ngridmax
          do i=1,ngrid
            ind_cell(i)=iskip+ind_grid(i)
          end do

          !store information about newly formed stars and their parent cells
          do i=1,ngrid

            !depleteion only in the leaf cells (OK? KS TODO)
            if(son(ind_cell(i))/=0) cycle

            x = (xg(ind_grid(i),1)+xc_level(ind,1)-skip_loc(1))*scale          !x
            y = (xg(ind_grid(i),2)+xc_level(ind,2)-skip_loc(2))*scale          !y
            z = (xg(ind_grid(i),3)+xc_level(ind,3)-skip_loc(3))*scale          !z

            r2 = (x - pos(1))**2 +  (y - pos(2))**2 +  (z - pos(3))**2
            r = sqrt(r2)

            !change density of cells within the SFC according to a given factor
            ! Note:
            ! - at sweep_level = ilevel (= nlevelmax), all uold was converted to primitive variables,
            !   and will be converted back to conserved variables later
            if (r < radius) then
              !update ufront
              do ivar=1,nvar
                ufront(ind_cell(i),ivar) = f_adjust * ufront(ind_cell(i),ivar)
              end do

              !update uold/unew
              if (sweep_level == nlevelmax) then                    !at ilevel (= nlevelmax)
                !change only density (others are automatically updated)
                uold(ind_cell(i),1) = ufront(ind_cell(i),1)
              else                                                  !at lower levels
                !change all conservative vars
                do ivar=1,nvar
                  !update uold (to be used for calculating flux by godunov method)
                  uold(ind_cell(i),ivar) = ufront(ind_cell(i),ivar)
                  !updat unew (variation of ufront used to ensure conservation)
                  unew(ind_cell(i),ivar) = unew(ind_cell(i),ivar) - (1d0-f_adjust)/f_adjust * ufront(ind_cell(i),ivar)
                end do
              end if
            end if
          end do
        end do
        ! End loop over cells
      end do
      ! End loop over grids
    end do
    ! End loop over levels


    !KS DEBUG
    if (myid == 1) then
      print '(A,1P5E10.2,/,A,1P4E10.2)', &
        '(crop_and_adjust_SFC, CROP) SFC_f_rhocrop, nH_SFC, nH_crop, radius, mass_crop: ', &
        SFC_f_rhocrop, rho_SFC*scale_nH, rho_crop*scale_nH, radius*scale_l/3.08568d18, mass_crop*scale_msun, &
        '(crop_and_adjust_SFC, ADJUST) mass_SFC, mass_star, mass_crop, f_adjust: ',&
        mass_SFC*scale_msun, mass_star*scale_msun, mass_crop*scale_msun, f_adjust
    end if


  end subroutine crop_and_adjust_SFC


  !----------------------------------------------------------------------------
  ! logfile for star forming cloud (KS ADDED)
  !----------------------------------------------------------------------------
  ! take a log of information about the cells where stars newly form
  !  - vars saved in the log file
  !   * index_star, type, mass,
  !   * step, cpuid, ind_cell, level, mcell,
  !   * nH, T, yHI, yHII, yH2, yHeI, yHeII, yHeIII, Zsol, pos(:), vel(:)
  !
  !  - input from main code
  !   * nnew: number of new particles in each cpu
  !   * ind_cell_new: cell index corresponding to new particles in each cpu
  !   * ind_part: cell index corresponding to new particles in each cpu
  !
  !----------------------------------------------------------------------------
  subroutine write_log_SFC(iSFC, pos, radius, mass, rho, Zmetal, m_star_form, n_star_form)
    !in
    integer, intent(IN) :: iSFC !index of SFC at each step (usually one)
    real(dp),dimension(1:3),intent(in) :: pos  !center
    real(dp),intent(in) :: radius !radius of SFC
    real(dp),intent(in) :: mass !mass of SFC
    real(dp),intent(in) :: rho !average density of SFC
    real(dp),intent(in) :: Zmetal !average metallicity of SFC
    real(dp),intent(in) :: m_star_form !stellar mass formed from SFC
    integer,intent(in) :: n_star_form !stellar number formed from SFC

    character(LEN=80),parameter::filename='logSFC' !name of log file
    integer,parameter :: LUN = 11, LUN2 = 12

    real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
    real(dp)::zred,scale_m

    real(dp),parameter :: pc_cgs=3.08568d18 !need to be exactly the same as the value in output_analysis.f90
    real(dp),parameter :: box_cen = 0.5d0 ! boxlen=nx_loc=1 is assumed

    real(dp) :: Zmetal_Zsun !metallicity in Zsun
    !real(dp) :: f_star
    !integer :: n_star_form

    !for adjust file
    character(LEN=80)::filename_tmp
    character(1024) :: line
    integer :: nstep_coarse_line, nstep_line, ios
    logical :: is_exist
    integer :: status
    character(200) :: command, message


    integer,save :: ifirst = 0 !for first call

    !-------------------------------
    ! for first call
    !-------------------------------
    ! fix logSFC by discarding (possibly broken) previous data in the future of the restart time
    if (ifirst == 0) then
      if(myid==1) then
        filename_tmp = TRIM(filename)//".tmp"
        inquire(file=filename, exist=is_exist)

        !既にファイルが存在する場合
        if (is_exist) then
          print '(A)', "* discarding (possibly broken) unnecessary data from "//TRIM(filename)

          !*** step 1: copy original file to tmp file except the lines corresponding to the future of restarting time
          !            original file will be deleted after copying its data
          open(LUN,file=filename,status='old')          !元ファイル
          open(LUN2,file=filename_tmp,status='replace') !仮ファイル
          !最初の行はコメントなので単純にコピー
          read(LUN, fmt='(a)', iostat=ios) line
          write(LUN2, fmt='(a)', iostat=ios) TRIM(line)
          !データを順に読み込み
          read(LUN, fmt='(a)', iostat=ios) line
          do while(ios == 0) !ファイル終端まで読み込み
            read(line(1:10), *) nstep_coarse_line      !最初の10桁がnstep_coarse
            read(line(11:20), *) nstep_line            !次の10桁がnstep (fine step)
            !nstep_coarse、または、nstepが現在以上だったらそれ以降のデータは再利用しない
            if ( (nstep_coarse_line .ge. nstep_coarse) .or. (nstep_line .ge. nstep)) exit
            write(LUN2, fmt='(a)', iostat=ios) TRIM(line)
            read(LUN, fmt='(a)', iostat=ios) line
          end do
          close(LUN, status = 'delete') !delte the original file
          close(LUN2)

          !*** step 2: copy tmp file to original file
          !            tmp file will be deteleted after copying its data
          open(LUN2,file=filename_tmp,status='old') !仮ファイル
          open(LUN,file=filename,status='replace') !元ファイル
          !データを順に読み込み
          read(LUN2, fmt='(a)', iostat=ios) line
          do while(ios == 0) !ファイル終端まで読み込み
            write(LUN, fmt='(a)', iostat=ios) TRIM(line)
            read(LUN2, fmt='(a)', iostat=ios) line
          end do
          close(LUN2, status = 'delete') !delte the temporary file
          close(LUN)

        else !初めてファイルを作成する場合

          open(LUN, file=filename)
          write (LUN,'(A, A)')  '#coarse_step, fine_step, zred, iSFC, ',&
            'R_SFC [pc], M_SFC [Msun], N_star, M_star [Msun], nH [cm-3], Zmetal [Zsun], pos(:) [pc]'
          close(LUN)
        end if

      end if
      ifirst = 1
    end if

    !----------------------------------------
    ! write log file
    !----------------------------------------
    !units
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    scale_m=scale_d*scale_l**3

    !save sink properties to log files
    if(myid==1)then
      !redshift
      zred = 1./aexp - 1.
      !write data to file
      open(LUN, file=filename, position='APPEND')

      !ほぼゼロのE-100みたいな値が化けることがあったので対応しておく
      Zmetal_Zsun = Zmetal/0.02d0
      if (Zmetal_Zsun < 1d-50) then
        Zmetal_Zsun = 0d0
      end if

      ! !形成する星の質量をconvert_SFC_to_starsと同様に計算
      ! !use formula based on small-scale simulations if SFC_f_star is negative
      ! if (SFC_f_star < 0d0) then
      !   f_star = f_star_tot_HRG19(rho, mass, Zmetal) !star formation efficiency from RHD MC-scale simulations (He, Ricotti, Geen 2019)
      ! else !use given SFC_f_star
      ! f_star = SFC_f_star
      ! end if
      ! ! f_star = SFC_f_star
      ! n_star_form = max(int(mass*f_star/mstar), 1) !form at least one star particle
      ! !stellar mass
      ! m_star_form = n_star_form*mstar


      !coarse step, fine step, zred, iSFC, R_SFC [pc], M_SFC [Msun], M_star [Msun], nH [cm-3], Zmetal [Zsun], pos(:)
      write (LUN,'(2I10, 1P1E14.6, I4, 1P2E11.3, I8, 1P3E11.3, 1P3E20.12)') &
        nstep_coarse, nstep, zred, &
        iSFC, radius*scale_l/pc_cgs, mass*scale_msun, n_star_form, m_star_form*scale_msun, rho*scale_nH, Zmetal_Zsun, pos(1:3)*scale_l/pc_cgs
      close(LUN)
    end if

  end subroutine write_log_SFC

  !----------------------------------------------------------------------------
  ! check_cell_health (KS ADDED)
  !----------------------------------------------------------------------------
  ! check all cells have healthy values
  ! currently only check density value
  !----------------------------------------------------------------------------
  subroutine check_cell_health

    !variables
    integer :: sweep_level !search for all leaf cells through levels (not only ilevel)
    integer :: ix, iy, iz
    real(dp) :: x, y, z, dx_cell
    real(dp),dimension(1:twotondim,1:3)::xc_level


    real(dp),parameter :: nH_warn_low = 1d-10 !warn if nH < nH_warn_low
    real(dp),parameter :: nH_warn_high = 1d10  !warn if nH > nH_warn_high


    !----------------------------------
    ! crop high-density regions
    !----------------------------------
    do sweep_level = levelmin,nlevelmax
      dx_cell = 1d0/dble(2**sweep_level)

      ! Cells center position relative to grid center position
      do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc_level(ind,1)=(dble(ix)-0.5D0)*dx_cell
        xc_level(ind,2)=(dble(iy)-0.5D0)*dx_cell
        xc_level(ind,3)=(dble(iz)-0.5D0)*dx_cell
      end do

      ! sum shell qunatities in each cpu
      ncache=active(sweep_level)%ngrid
      do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
          ind_grid(i)=active(sweep_level)%igrid(igrid+i-1)
        end do
        ! Loop over cells
        do ind=1,twotondim
          iskip=ncoarse+(ind-1)*ngridmax
          do i=1,ngrid
            ind_cell(i)=iskip+ind_grid(i)
          end do

          !store information about newly pformed stars and their parent cells
          do i=1,ngrid
            !only the leaf cells
            if(son(ind_cell(i))/=0) cycle

            x = (xg(ind_grid(i),1)+xc_level(ind,1)-skip_loc(1))*scale          !x
            y = (xg(ind_grid(i),2)+xc_level(ind,2)-skip_loc(2))*scale          !y
            z = (xg(ind_grid(i),3)+xc_level(ind,3)-skip_loc(3))*scale          !z

            !-------------- warn if nH < nH_warn_low --------------!
            if (uold(ind_cell(i),1)*scale_nH < nH_warn_low .or. &
              unew(ind_cell(i),1)*scale_nH < nH_warn_low .or. &
              ufront(ind_cell(i),1)*scale_nH < nH_warn_low) then
              print '(A,/,A,2I8,I8,1P7E11.3, 1P3E20.12)', '** WARNING ** very small density cell found ** WARNING **',&
                'myid, level, ind_cell, uold(1), unew(1), ufront(1), nHold, nHnew, nHfront, nH_warn_low, x, y, z: ',&
                myid, sweep_level, ind_cell(i), uold(ind_cell(i),1), unew(ind_cell(i),1), ufront(ind_cell(i),1), &
                uold(ind_cell(i),1)*scale_nH, unew(ind_cell(i),1)*scale_nH, ufront(ind_cell(i),1)*scale_nH, &
                nH_warn_low, x, y, z
            end if

            !-------------- warn if nH > nH_warn_high --------------!
            if (uold(ind_cell(i),1)*scale_nH > nH_warn_high .or. &
              unew(ind_cell(i),1)*scale_nH > nH_warn_high .or. &
              ufront(ind_cell(i),1)*scale_nH > nH_warn_high) then
              print '(A,/,A,2I8,I8,1P7E11.3, 1P3E20.12)', '** WARNING ** very small density cell found ** WARNING **',&
                'myid, level, ind_cell, uold(1), unew(1), ufront(1), nHold, nHnew, nHfront, nH_warn_high, x, y, z: ',&
                myid, sweep_level, ind_cell(i), uold(ind_cell(i),1), unew(ind_cell(i),1), ufront(ind_cell(i),1), &
                uold(ind_cell(i),1)*scale_nH, unew(ind_cell(i),1)*scale_nH, ufront(ind_cell(i),1)*scale_nH, &
                nH_warn_high, x, y, z
            end if
          end do
        end do
        ! End loop over cells
      end do
      ! End loop over grids
    end do
    ! End loop over levels
  end subroutine check_cell_health


  ! !----------------------------------------------------------------------------
  ! ! set_ufront
  ! !----------------------------------------------------------------------------
  ! ! set hydro variables ufront (not necessarily synchronized )
  ! ! - ufront: forefront hydro variables
  ! !   * consider all the contribution of advection so far
  ! !   * similar to unew but include the advection effect through virtual (cpu) boundaries
  ! !   * used in the SFC subroutines to detemine the propoerty of SFC and the ammount of gas depletion
  ! !   * used in the non-local gas depletion subtourine
  ! ! - how to get ufront
  ! !   * first, copy unew to ufront
  ! !   * second, consider cpu boundary effect using make_virtual_reverse_dp
  ! !----------------------------------------------------------------------------
  ! subroutine set_ufront(ilevel)
  !   use amr_commons
  !   use hydro_commons
  !   use poisson_commons
  !   implicit none
  !   integer::ilevel
  !   integer::i,ivar,irad,ind,iskip,nx_loc,ind_cell,icpu

  !   !------ copy unew to ufront -----!
  ! do ind=1,twotondim
  !    iskip=ncoarse+(ind-1)*ngridmax
  !     do ivar=1,nvar
  !       !active grids
  !    do i=1,active(ilevel)%ngrid
  !         ufront(active(ilevel)%igrid(i)+iskip,ivar) = unew(active(ilevel)%igrid(i)+iskip,ivar)
  !       end do
  !       !virtual boundary cells
  !       do icpu=1,ncpu
  !         do i=1,reception(icpu,ilevel)%ngrid
  !           ufront(reception(icpu,ilevel)%igrid(i)+iskip,ivar)=unew(reception(icpu,ilevel)%igrid(i)+iskip,ivar)
  !         end do
  !       end do
  !    end do
  ! end do

  !   ! !KS DEBUG TENTATIVE
  !   ! do ind=1,twotondim
  !   !    iskip=ncoarse+(ind-1)*ngridmax
  !   !    do i=1,active(ilevel)%ngrid
  !   !       if (myid == 5 .and. active(ilevel)%igrid(i)+iskip == 4821) then
  !   !          print '(A,2I4,I8,1P3E12.4)', &
  !   !               '(KS DEBUG TENTATIVE in set_ufront before make_virtual_reverse_dp) myid, level, ind_cell, nH(uold), nH(unew), nH(ufront): ',&
  !   !               myid, ilevel, active(ilevel)%igrid(i)+iskip, uold(active(ilevel)%igrid(i)+iskip,1)*scale_nH, unew(active(ilevel)%igrid(i)+iskip,1)*scale_nH,&
  !   !               ufront(active(ilevel)%igrid(i)+iskip,1)*scale_nH
  !   !       end if
  !   !    end do
  !   ! end do

  !   !------ consider virtual boundary advection ------!
  !   do ivar=1,nvar
  !     call make_virtual_reverse_dp(ufront(1,ivar),ilevel)
  ! end do

  !   ! !KS DEBUG TENTATIVE
  !   ! do ind=1,twotondim
  !   !    iskip=ncoarse+(ind-1)*ngridmax
  !   !    do i=1,active(ilevel)%ngrid
  !   !       if (myid == 5 .and. active(ilevel)%igrid(i)+iskip == 4821) then
  !   !          print '(A,2I4,I8,1P3E12.4)', &
  !   !               '(KS DEBUG TENTATIVE in set_ufront after make_virtual_reverse_dp) myid, level, ind_cell, nH(uold), nH(unew), nH(ufront): ',&
  !   !               myid, ilevel, active(ilevel)%igrid(i)+iskip, uold(active(ilevel)%igrid(i)+iskip,1)*scale_nH, unew(active(ilevel)%igrid(i)+iskip,1)*scale_nH,&
  !   !               ufront(active(ilevel)%igrid(i)+iskip,1)*scale_nH
  !   !       end if
  !   !    end do
  !   ! end do



  ! end subroutine set_ufront

  !----------------------------------------------------------------------------
  ! f_star_tot_HRG19
  !----------------------------------------------------------------------------
  ! calculate star formation efficiency of SFC based on Eq. 9 of He, Ricotti, Geen (2019; HRG19)
  !  - assume 1/5 times smaller efficiency considering the metallicity is low (see Fig. 14)
  !  - cap the efficiency at 80%
  !
  ! input:
  !  - rho_SFC: density of SFC
  !  - mass_SFC: mass of SFC
  ! output:
  !  - f_star_tot_HRG19: star formation efficiency (M_star/M_SFC)
  !----------------------------------------------------------------------------
  function f_star_tot_HRG19(rho_SFC, mass_SFC, Z_SFC)
    use amr_commons
    use hydro_commons
    use poisson_commons
    implicit none
    real(dp), intent(IN):: rho_SFC, mass_SFC
    real(dp), intent(IN):: Z_SFC !TO ADDED 2026-11-10: SFC average metallicity (mass fraction; Zsun = 0.02 as elsewhere in this file)

    real(dp) :: nH_SFC, M_SFC, f_star_tot_HRG19, f_metal, Z_Zsun
    !real(dp),parameter :: ncri = 1d3 !critical density in Eq. 9 of HRG19
!*********************************************************
    real(dp),parameter :: ncri = 1d2 !critical density in Eq. 9 of HRG19 with modification (10 times smaller)
!*********************************************************
    real(dp),parameter :: f_lowmetal = 0.2d0 !factor to consider the reduction of efficient for low-metallicity case
    real(dp),parameter :: f_star_tot_max = 0.8d0 !maximum efficiency


    !change units
    nH_SFC = rho_SFC * scale_nH
    M_SFC = mass_SFC * scale_msun

    !low-metallicity factor (TO ADDED 2026-11-10, namelist SFC_f_star_Zdep):
    ! 0: fixed 0.2 = He, Ricotti & Geen 2019 ("f* reduced by ~5 for Z < 0.1 Zsun"); the setting of the Sugimura+2024 runs
    ! 1: 0.2 x (Z/1e-3 Zsun)^0.25 as written in Sugimura, Ricotti, Park, Garcia & Yajima 2024 (ApJ 970, 14) eq. 2
    !    (0.004 = 0.2 x 2e-2 at Z = 1e-3 Zsun), Z floored at the Pop III threshold Zpop3
    ! 2: max(0.2, min(1, Z/Zsun)) = He+2019 sec. 3.4: "f* decreases almost linearly with the metallicity from solar to
    !    Z = 0.2 Zsun" and "saturates when further lowering the metallicity" (1/5 at 0.1 and 0.025 Zsun)
    Z_Zsun = Z_SFC/0.02d0   !SFC metallicity in Zsun (0.02 as elsewhere in this file)
    select case (SFC_f_star_Zdep)
    case (1)
      f_metal = f_lowmetal * (max(Z_Zsun, Zpop3)/1d-3)**0.25d0
    case (2)
      f_metal = max(f_lowmetal, min(1d0, Z_Zsun))
    case default
      f_metal = f_lowmetal
    end select

    !compute the efficiency
    f_star_tot_HRG19 = f_metal * 2d-2 * (M_SFC/1d4)**0.4d0 * (1 + nH_SFC/ncri)**0.91d0

    !impose maximum
    f_star_tot_HRG19 = Min(f_star_tot_HRG19, f_star_tot_max)

  end function f_star_tot_HRG19

  !------------------------------ KS ADDED END ----------------------------------------!


  !----------------------------------------------------------------------------
  ! get_cluster_mstar_list
  !----------------------------------------------------------------------------
  ! get stellar mass list for a star-by-star model
  ! input:
  !  - mass_cluster: mass of SFC
  ! output:
  !  - n_star_form: total number of star formed
  !  - mstar_list: list of stellar mass
  !----------------------------------------------------------------------------

  subroutine get_cluster_mstar_list(mass_cluster,Zmetal,n_star_form,mstar_list)
    use amr_commons
    use hydro_commons
    use poisson_commons
    use sbs_imf
    implicit none
    integer, parameter :: PRIMARY_RANK = 1
    real(dp), intent(IN):: mass_cluster
    real(dp), intent(IN):: Zmetal
    integer, intent(OUT):: n_star_form
    real(dp),dimension(1:Nnew_max),intent(OUT):: mstar_list
    integer :: i_star
    integer :: iM
    real(dp) ::  M_cluster, M_sum, rand, logM, M_resolved, M_unresolved, Mass, delta_iM, Zmetal_Zsun
    real(dp) :: Pcum(sbs_NM), Mfrac_unresolved
    integer :: N_unresolved

    if (myid == PRIMARY_RANK) then !myid = PRIMARY_RANKでlistを作成 -> MPI通信

      !change units
      M_cluster = mass_cluster * scale_msun
      Zmetal_Zsun = Zmetal/0.02d0

      !get Pcum and Mfrac for the given metallicity
      call get_Pcum_at_Z(Zmetal_Zsun, Pcum)
      call get_Mfrac_unresolved_at_Z(Zmetal_Zsun, Mfrac_unresolved)


      !get number of unresolved star particles
      M_unresolved = Mfrac_unresolved * M_cluster
      N_unresolved = nint(M_unresolved/sbs_Mmin) !Nearest Integer
      !M_resolved = M_cluster-M_unresolved
      M_resolved = M_cluster-N_unresolved*sbs_Mmin !subtract actual mass of formed unresolved stars

      print '(2(A,1PE12.4),/,A,1PE12.4,A,I8,/,A,1PE12.4,A)', "star-by-star: get mstar_list according to given IMF for M_cluster[Msun] = ", M_cluster, ", Zmetal[Zsun] = ", Zmetal, &
        "mass fraction of unresolved stars = ", Mfrac_unresolved, " [Msun] ---> number of unresolved stars = ", N_unresolved, &
        "the rest of mass = ", M_resolved, " [Msun] is for resolved stars, which will be drawn soon."

      !check if N_unresolved <= Nnew_max
      if (N_unresolved > Nnew_max) then
        print '(2(A,I8),/,A,/,A)', "N_unresolved = ", N_unresolved, " is larger than  Nnew_max = ", Nnew_max,&
          "increase Nnew_max",&
          "STOP HERE"
        stop
      end if

      !set list for unresolved star particles
      mstar_list(1:N_unresolved) = sbs_Mmin/scale_msun

      print '(A,I8, 3(A,1PE12.4))', "N_unresolved = ", N_unresolved,", sbs_Mmin[Msun] = ", sbs_Mmin, &
        ", M_unresolved[Msun] = ", M_unresolved, ", M_unresolved_form[Msun] = ", N_unresolved * sbs_Mmin



      !set list for resolved star particles
      M_sum = 0d0
      n_star_form = Nnew_max ! TO ADDED: was unset if the loop below ends without exit (e.g. N_unresolved==Nnew_max) (2026-10-30)
      do i_star = N_unresolved+1, Nnew_max
        call random_value_0to1(rand, initial_seed=nstep)  !get random number between 0 and 1 (pass nstep for initial_seed to avoid the same pattern)

        do iM =1, sbs_NM-1
          if (rand <= Pcum(iM+1)) then
            !interpolation
            delta_iM = (rand - Pcum(iM))/(Pcum(iM+1)-Pcum(iM))
            Mass = 10d0**( log10(sbs_Mmin) + (iM+delta_iM-1)*(log10(sbs_Mmax)-log10(sbs_Mmin))/(sbs_NM-1) )
            exit
          end if
        end do

        !KS DEBUG
        ! print *, "(KS DEBUG) Pcum(iM), rand, Pcum(iM+1), iM, delta_iM, Mass: ", Pcum(iM), rand, Pcum(iM+1), iM, delta_iM, Mass

        ! if adding the next star increases the deviation from the target mass, skip this star and end updating the list
        if (abs(M_sum - M_resolved) < (M_sum+Mass - M_resolved)) then
          n_star_form = i_star-1
          exit
        end if

        !update the list
        mstar_list(i_star) = Mass/scale_msun !stellar mass in code_unit
        M_sum = M_sum + Mass

        ! !KS DEBUG
        ! print '(A,I8, 4(A,1PE12.4))', "(KS DEBUG) i_star = ", i_star, ", rand = ", rand, ", M[Msun] = ", Mass, ", M_sum[Msun] = ", M_sum

        !terminal condition
        if (M_sum > M_resolved) then
          n_star_form = i_star
          exit
        end if

        !check if the number of stars is not too large
        if (i_star == Nnew_max .and. M_sum < M_resolved) then
          print '(A,I8,2(A,1PE12.4),/,A,/,A)', "i_star reaches Nnew_max = ", i_star,  ", but M_sum = ", M_sum,&
            " is still smaller than M_cluster = ", M_cluster,&
            "increase Nnew_max",&
            "STOP HERE"
          stop
        end if
      end do
      print '(4(A,1PE12.4),A,I8)', "summary for get_cluster_mstar_list: M_cluster[Msun] = ", M_cluster,&
        ", M_unresolved_form[Msun] = ", N_unresolved*sbs_Mmin,", M_resolved_form[Msun] = ", M_sum,  &

        ", M_form[Msun] = ", M_sum+N_unresolved*sbs_Mmin, ", n_star_form =", n_star_form

    end if
    call mpi_bcast(n_star_form, 1, MPI_INTEGER, PRIMARY_RANK-1, MPI_COMM_WORLD, info) !one smaller value from id shuould be set in mpi_bcast
    call mpi_bcast(mstar_list, n_star_form, MPI_DOUBLE_PRECISION, PRIMARY_RANK-1, MPI_COMM_WORLD, info) !one smaller value from id shuould be set in mpi_bcast
  end subroutine get_cluster_mstar_list


  ! ! -----------------------------------------------------------------
  ! ! read IMF
  ! ! -----------------------------------------------------------------
  ! subroutine readIMF
  !   integer, parameter :: PRIMARY_RANK = 1
  !   !real(dp),dimension(:),pointer,intent(OUT) :: r1d, rho1d
  !   !real(dp),dimension(:),allocatable,intent(OUT) :: r1d, rho1d
  !   integer :: n, nmax
  !   character(len=100):: fn, dir
  !   integer, save :: ifirst = 0

  !   !call only once
  !   if (ifirst == 1) return
  !   if (ifirst == 0) then
  !     ifirst = 1
  !   end if

  !   fn = IMF_file  !path_to_directory/ChabrierIMF.dat

  !   if (myid == PRIMARY_RANK) then !myid = PRIMARY_RANKが読み込み
  !     write(*,*) 'file(IMF)=',fn
  !     open(1, file=fn, status='old')

  !     !nmax読み込み
  !     read(1,'(I17)') nmax_IMF
  !     print *, "nmax_IMF=",nmax_IMF

  !     !data読み込み
  !     allocate( logM_IMF(1:nmax_IMF), Pcum_IMF(1:nmax_IMF) )
  !     do n = 1, nmax_IMF
  !       !read(1,'(2E18.8)') logM_IMF(n), Pcum_IMF(n)
  !       read(1,*) logM_IMF(n), Pcum_IMF(n)
  !       ! print *, "n,logM_IMF(n), Pcum_IMF(n) = ", n,logM_IMF(n), Pcum_IMF(n)
  !     enddo
  !     close(1)
  !   endif

  !   !mpi_bcastで共有
  !   call mpi_bcast(nmax_IMF, 1, MPI_INTEGER, PRIMARY_RANK-1, MPI_COMM_WORLD, info)

  !   if(myid /= PRIMARY_RANK) then
  !     allocate(logM_IMF(1:nmax_IMF), Pcum_IMF(1:nmax_IMF))
  !   end if
  !   call mpi_bcast(logM_IMF, size(logM_IMF), MPI_DOUBLE_PRECISION, PRIMARY_RANK-1, MPI_COMM_WORLD, info) !one smaller value from id shuould be set in mpi_bcast
  !   call mpi_bcast(Pcum_IMF, size(Pcum_IMF), MPI_DOUBLE_PRECISION, PRIMARY_RANK-1, MPI_COMM_WORLD, info) !one smaller value from id shuould be set in mpi_bcast

  !   if(myid == PRIMARY_RANK) &
  !     print *, "(KS DEBUG) myid, logM_IMF(1), logM_IMF(nmax_IMF), Pcum_IMF(1), Pcum_IMF(nmax_IMF): ", &
  !     myid, logM_IMF(1), logM_IMF(nmax_IMF), Pcum_IMF(1), Pcum_IMF(nmax_IMF)

  ! end subroutine readIMF

  ! !-----------------------------------------------------------
  ! ! Get cpu of the grid that contains the given position
  ! ! using cmp_cpumap function
  ! !-----------------------------------------------------------
  ! subroutine get_cpu_of_pos(pos,ilevel,icpu)
  !   real(dp),dimension(1:3),intent(IN)::pos
  !   integer,intent(IN)::ilevel
  !   integer,intent(OUT)::icpu
  !   real(dp),dimension(1:nvector,1:3)::pos2
  !   integer,dimension(1:nvector)::cc
  !   real(dp) :: dx_cell, dx_grid
  !   integer :: ix_cell,iy_cell,iz_cell,ix_grid,iy_grid,iz_grid


  !   dx_cell = 0.5D0**ilevel
  !   dx_grid = 2*dx_cell

  !   !integer position of cell containing the pos
  !   !nint: nearest integer
  !   ! ix_cell = nint(pos(1)/dx_cell+0.5d0)
  !   ! iy_cell = nint(pos(2)/dx_cell+0.5d0)
  !   ! iz_cell = nint(pos(3)/dx_cell+0.5d0)

  !   !integer position of grid containing the pos
  !   ix_grid = int(pos(1)/dx_grid)+1
  !   iy_grid = int(pos(2)/dx_grid)+1
  !   iz_grid = int(pos(3)/dx_grid)+1

  !   !real grid position
  !   pos2(1,1:3) = (/(ix_grid-0.5d0)*dx_grid, (iy_grid-0.5d0)*dx_grid, (iz_grid-0.5d0)*dx_grid/)

  !   ! print '(A,1P3E20.12)', "KS DEBUG (get_cpu_of_pos) pos =         ", pos(1:3)
  !   ! print '(A,1P3E20.12)', "KS DEBUG (get_cpu_of_pos) pos_of_grid = ", pos2(1,1:3)
  !   ! print '(A,3I8)', "KS DEBUG (get_cpu_of_pos) ipos =      ", ix_cell, iy_cell, iz_cell
  !   ! print '(A,3I8)', "KS DEBUG (get_cpu_of_pos) ipos_grid = ", ix_grid, iy_grid, iz_grid
  !   ! print *, "KS DEBUG STOP HERE"
  !   ! stop

  !   !compute cpumap for the grid position
  !   call cmp_cpumap(pos2,cc,1)
  !   icpu = cc(1)
  ! end subroutine get_cpu_of_pos

#endif !RT
end subroutine star_formation_fg
!################################################################
!################################################################
! collect_active_SFC: global list of the active SFC particles (FAM_SFC) in
! SFC_n_active / SFC_pos_active / SFC_radius_active (pm_commons).
! From patch/ks cooling_fine.f90 (merge_ks); used by SFC_from_candidate for the
! SFC markers (pop2_SFC_marker). The cooling-side use (SFC_flg_thin) is not ported.
!################################################################
subroutine collect_active_SFC
  use amr_commons
  use hydro_commons
  use cooling_module
  use pm_commons
  use mpi_mod
  implicit none

  integer :: n_loc, offset, i, info, icpu, iglob, iloc
  real(dp),dimension(1:3,1:SFC_nmax_active) :: pos_loc
  real(dp),dimension(1:SFC_nmax_active) :: radius_loc
  integer,dimension(1:ncpu) :: n_list

  !buffer for MPI send/recv
  integer,dimension(1:ncpu) :: sbufi,rbufi
  real(dp),dimension(1:SFC_nmax_active*4) :: sbufr,rbufr


  !------------------------------------
  ! create global particle number list
  !------------------------------------
  !count local SFC particles
  n_loc = 0   !initialization
  do i=1,npartmax
    if( is_SFC(typep(i)) ) then
      n_loc=n_loc+1
      ! check the number of candidates
      if (n_loc > SFC_nmax_active) then
        write(*,*)'Number of local SFC partilces exceeds the limit ',n_loc, SFC_nmax_active
        write(*,*)'stopping...'
        stop
      end if
      !store the positions/radius
      pos_loc(1:3,n_loc) = xp(i,1:3)  !SFC position
      radius_loc(n_loc) = -zp(i)      !minus SFC radius is stored in zp
    end if
  enddo

  !create n_list
  n_list(:) = 0
  n_list(myid) = n_loc
  sbufi(1:ncpu) = n_list(:)
  call MPI_ALLREDUCE(sbufi,rbufi,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  n_list(:) = rbufi(1:ncpu)
  SFC_n_active = sum(n_list) !number of active SFC

  !offset between local and global indexes
  offset = 0
  do icpu = 1, myid-1  ! "do i = 1,0" loop is just skipped (KS MEMO)
    offset = offset + n_list(icpu)
  end do

  !------------------------------------
  !create global position/radius list
  !------------------------------------
  !initialization
  SFC_pos_active(:,:)=0d0 !position of active SFC
  SFC_radius_active(:)=0d0 !radius of active SFC

  !MPI communication
  do iloc = 1, n_loc
    iglob = iloc + offset
    SFC_pos_active(1:3,iglob) = pos_loc(1:3,iloc)
    SFC_radius_active(iglob) = radius_loc(iloc)
  end do
  do iglob = 1, SFC_n_active
    sbufr((iglob-1)*4+1) = SFC_pos_active(1,iglob)
    sbufr((iglob-1)*4+2) = SFC_pos_active(2,iglob)
    sbufr((iglob-1)*4+3) = SFC_pos_active(3,iglob)
    sbufr((iglob-1)*4+4) = SFC_radius_active(iglob)
  end do
  call MPI_ALLREDUCE(sbufr, rbufr, 4*SFC_n_active, MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  do iglob = 1, SFC_n_active
    SFC_pos_active(1,iglob) = rbufr((iglob-1)*4+1)
    SFC_pos_active(2,iglob) = rbufr((iglob-1)*4+2)
    SFC_pos_active(3,iglob) = rbufr((iglob-1)*4+3)
    SFC_radius_active(iglob) = rbufr((iglob-1)*4+4)
  end do

end subroutine collect_active_SFC
