! merge_ks (2026-11-13): ported from patch/ks (SN reservoir of the FG feedback). Changes: use mpi_mod,
! mH/kB from module constants. Everything else as in ks (TO ADDED/MODIFIED comments are ks history).
! TO ADDED: Read-only sphere statistics for SN feedback design.
module sn_sphere_stats
  use amr_commons, only: dp, ndim, twotondim, ncoarse, ngridmax, levelmin, nlevelmax, &
    active, son, xg, boxlen, icoarse_min, icoarse_max, jcoarse_min, jcoarse_max, &
    kcoarse_min, kcoarse_max, metal, Ncell_SN_max ! TO MODIFIED: Namelist search cap.
  use hydro_parameters, only: imetal, gamma ! TO MODIFIED: Temperature criterion uses the configured EOS.
  use constants, only: mH, kB ! merge_ks: takashi keeps the constants in module constants
#ifdef RT
  ! TO ADDED: Use the same species indices and H/He abundances as getMu.
  use rt_parameters, only: iIons, nIons, isH2, isHe, ixHI, ixHII, ixHeII, ixHeIII
  use cooling_module, only: X, Y
#endif
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mpi_mod ! merge_ks: ks used 'use mpi, only: ...'
  implicit none
  private
  public :: sphere_stats, get_sn_sphere_stats
  public :: get_sn_host_stats_many ! TO ADDED: Owner-only host lookup for globally shared SN events.
  public :: sn_cooling_radius_pc ! TO ADDED: User-specified cooling-radius fit.
  ! TO ADDED: Batched radii and explicit single-cell fallback decisions.
  public :: get_sn_sphere_stats_many, find_sn_feedback_radius, select_sn_feedback_radius
  public :: get_sn_sphere_stats_events !TO ADDED: All event/candidate shells in one collective.
  public :: make_sn_candidate_radii ! TO ADDED: r_max = n_radius * host_dx.
  public :: sn_radius_choice, sn_radius_sphere, sn_radius_unresolved, sn_radius_no_solution
  integer, parameter :: sn_radius_sphere=0, sn_radius_unresolved=1, sn_radius_no_solution=2
  ! TO MODIFIED: The 1e9 K acceptance limit is gone; DeltaT is diagnostic only.
  real(dp), parameter :: pc_in_cm=3.08568e18_dp

  type :: sphere_stats
    real(dp) :: mass = 0.0_dp
    real(dp) :: metal_mass = 0.0_dp
    real(dp) :: volume = 0.0_dp
    real(dp) :: mean_density = 0.0_dp
    real(dp) :: mean_metallicity = 0.0_dp
    ! TO ADDED: Sum(m/mu) has code-mass units; divide physical value by mH
    ! to obtain the total number of gas particles, including electrons.
    real(dp) :: mass_over_mu = 0.0_dp
    real(dp) :: mean_mu = 0.0_dp
    logical :: mu_available = .false.
  end type sphere_stats

  ! TO ADDED: single_cell=true means the HOST CELL ONLY, never a sphere query.
  ! radius is the host half-diagonal in that case (diagnostic code length).
  type :: sn_radius_choice
    logical :: single_cell = .true.
    integer :: reason = sn_radius_no_solution
    integer :: candidate_index = 0
    real(dp) :: radius = 0.0_dp
    real(dp) :: cooling_radius_pc = 0.0_dp
    real(dp) :: delta_temperature = huge(0.0_dp)
    type(sphere_stats) :: stats
  end type sn_radius_choice

  ! TO ADDED: Scalar second argument generates candidates from n*dx;
  ! an array second argument preserves the explicit-candidate API.
  interface find_sn_feedback_radius
    module procedure find_sn_feedback_radius_candidates, find_sn_feedback_radius_ncell
    module procedure find_sn_feedback_radius_configured ! TO ADDED: Read the namelist multiplier.
  end interface find_sn_feedback_radius

contains

  ! TO ADDED: Find the leaf containing each particle position, using half-open
  ! cell bounds. Particle ownership need not equal cell ownership. All ranks
  ! supply identical positions and periodic flags; one reduction gathers all
  ! hosts, including their actual AMR widths. No ghost state is read.
  subroutine get_sn_host_stats_many(positions, state, hosts, host_dx, periodic)
    real(dp), intent(in) :: positions(:,:), state(:,:)
    type(sphere_stats), intent(out) :: hosts(size(positions,2))
    real(dp), intent(out) :: host_dx(size(positions,2))
    logical, intent(in) :: periodic(ndim)
    real(dp) :: local(6,size(positions,2)), global(6,size(positions,2))
    real(dp) :: scale, dx, width, volume, offset(ndim), lower(ndim), point(ndim)
    real(dp) :: origin(3), extent(3)
    integer :: ilevel, ind, idim, ig, grid, cell, event
#ifndef WITHOUTMPI
    integer :: ierr, real_type
#endif
    hosts=sphere_stats()
    host_dx=0.0_dp
    if (size(positions,1) /= ndim) call invalid_input('host position dimension mismatch')
    if (size(positions,2) == 0) return
    if (.not. all(ieee_is_finite(positions))) call invalid_input('host positions must be finite')
    if (size(state,1) < ncoarse+twotondim*ngridmax .or. size(state,2) < 1) &
      call invalid_input('host state has insufficient extent')
    if (metal) then
      if (imetal < 1 .or. imetal > size(state,2)) call invalid_input('invalid host metal column')
    endif
#ifdef RT
    call check_mu_columns(size(state,2))
#else
    call invalid_input('host chemical mu requires an RT build')
#endif
    origin=real([icoarse_min,jcoarse_min,kcoarse_min],dp)
    extent=real([icoarse_max-icoarse_min+1,jcoarse_max-jcoarse_min+1,kcoarse_max-kcoarse_min+1],dp)
    if (boxlen <= 0.0_dp .or. any(extent(1:ndim) <= 0.0_dp)) &
      call invalid_input('invalid host domain dimensions')
    scale=boxlen/extent(1)
    extent=extent*scale
    local=0.0_dp
    do ilevel=levelmin,nlevelmax
      dx=0.5_dp**ilevel
      width=dx*scale
      volume=width**ndim
      do ind=1,twotondim
        do idim=1,ndim
          offset(idim)=(real(ibits(ind-1,idim-1,1),dp)-1.0_dp)*dx
        enddo
        do ig=1,active(ilevel)%ngrid
          grid=active(ilevel)%igrid(ig)
          cell=ncoarse+(ind-1)*ngridmax+grid
          if (son(cell) /= 0) cycle
          lower=(xg(grid,1:ndim)+offset-origin(1:ndim))*scale
          do event=1,size(positions,2)
            point=positions(:,event)
            do idim=1,ndim
              if(periodic(idim))point(idim)=modulo(point(idim),extent(idim))
            enddo
            if (any(point < lower) .or. any(point >= lower+width)) cycle
            local(1,event)=local(1,event)+state(cell,1)*volume
            if(metal)local(2,event)=local(2,event)+state(cell,imetal)*volume
            local(3,event)=local(3,event)+volume
#ifdef RT
            local(4,event)=local(4,event)+cell_mass_over_mu(state,cell)*volume
#endif
            local(5,event)=local(5,event)+width
            local(6,event)=local(6,event)+1.0_dp
          enddo
        enddo
      enddo
    enddo
#ifndef WITHOUTMPI
    call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL,storage_size(0.0_dp)/8,real_type,ierr)
    if(ierr /= MPI_SUCCESS)call invalid_input('no matching MPI host datatype')
    call MPI_ALLREDUCE(local,global,size(local),real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    if(ierr /= MPI_SUCCESS)call invalid_input('host MPI_ALLREDUCE failed')
#else
    global=local
#endif
    do event=1,size(positions,2)
      if(global(6,event) /= 1.0_dp)call invalid_input('SN position must have exactly one owned leaf host')
      hosts(event)%mass=global(1,event)
      hosts(event)%metal_mass=global(2,event)
      hosts(event)%volume=global(3,event)
      hosts(event)%mean_density=global(1,event)/global(3,event)
      hosts(event)%mass_over_mu=global(4,event)
      if(global(1,event)>0.0_dp)hosts(event)%mean_metallicity=global(2,event)/global(1,event)
      if(global(4,event)>0.0_dp)then
        hosts(event)%mean_mu=global(1,event)/global(4,event)
        hosts(event)%mu_available=.true.
      endif
      host_dx(event)=global(5,event)
    enddo
  end subroutine get_sn_host_stats_many

  ! TO ADDED: nH in cm^-3, SN energy in erg, Z as a metal mass fraction.
  ! Zsolar defaults to 0.02, consistent with cooling_fine's Z/Zsolar conversion.
  ! This local calculation performs no MPI communication. Empty/massless
  ! spheres must be excluded by the caller rather than assigning nH=0 here.
  function sn_cooling_radius_pc(nH, energy_erg, metallicity, zsolar) result(radius_pc)
    real(dp), intent(in) :: nH, energy_erg, metallicity
    real(dp), intent(in), optional :: zsolar
    real(dp) :: radius_pc, solar_fraction, fz

    solar_fraction = 0.02_dp
    if (present(zsolar)) solar_fraction = zsolar
    if (.not. ieee_is_finite(nH) .or. nH <= 0.0_dp) &
      call invalid_input('cooling radius requires positive finite nH')
    if (.not. ieee_is_finite(energy_erg) .or. energy_erg <= 0.0_dp) &
      call invalid_input('cooling radius requires positive finite SN energy')
    if (.not. ieee_is_finite(metallicity) .or. metallicity < 0.0_dp) &
      call invalid_input('cooling radius requires nonnegative finite metallicity')
    if (.not. ieee_is_finite(solar_fraction) .or. solar_fraction <= 0.0_dp) &
      call invalid_input('cooling radius requires positive finite Zsolar')
    fz = max(0.01_dp, min(1.0_dp, metallicity/solar_fraction))
    radius_pc = 28.4_dp*nH**(-3.0_dp/7.0_dp) &
      *(energy_erg/1.0e51_dp)**(2.0_dp/7.0_dp)*fz**(-1.0_dp/7.0_dp)
  end function sn_cooling_radius_pc

  ! TO MODIFIED: All ranks must use the same center, radius, periodic, compute_mu.
  ! state is a conservative cell array: column 1 = rho, imetal = rho*Z.
  ! Its owned cells must be current; ghost values are never read. The caller
  ! chooses the snapshot (e.g. uold or ufront); no time synchronization is done.
  ! center and radius use the same code-length coordinates as particle xp.
  ! Include cell centers at distance <= radius, counting the FULL cell volume.
  ! Periodicity is explicit per axis and defaults to false (clipped domain).
  ! A periodic cell is counted once using its nearest image, even for large r.
  ! Empty selection returns all zeros; zero total mass gives metallicity zero.
  ! metal=false gives metal_mass=mean_metallicity=0 without reading imetal.
  ! TO MODIFIED: compute_mu=true reads the current RT chemical composition.
  ! mean_mu = sum(m)/sum(m/mu), not the arithmetic mass-weighted mean of mu.
  ! This describes the whole sphere's mass per particle; a temperature using
  ! this mu is a particle-weighted mean and does not bound every cell's T.
  subroutine get_sn_sphere_stats(center, radius, state, stats, periodic, compute_mu)
    real(dp), intent(in) :: center(ndim), radius
    real(dp), intent(in) :: state(:,:)
    type(sphere_stats), intent(out) :: stats
    logical, intent(in), optional :: periodic(ndim)
    logical, intent(in), optional :: compute_mu ! TO ADDED
    type(sphere_stats) :: batch(1)
    ! TO MODIFIED: Keep the original API as a one-radius batch.
    call get_sn_sphere_stats_many(center, [radius], state, batch, periodic, compute_mu)
    stats = batch(1)
  end subroutine get_sn_sphere_stats

  ! TO ADDED: Increasing candidate radii, one leaf traversal and one collective.
  ! Each cell enters its first enclosing shell. Prefix sums give every sphere.
  subroutine get_sn_sphere_stats_many(center, radii, state, stats, periodic, compute_mu)
    real(dp), intent(in) :: center(ndim), radii(:), state(:,:)
    type(sphere_stats), intent(out) :: stats(size(radii))
    logical, intent(in), optional :: periodic(ndim), compute_mu
    !TO MODIFIED: Preserve the single-event API through the shared batched kernel.
    call get_sn_sphere_stats_events(reshape(center,[ndim,1]),radii,[1,size(radii)+1], &
      state,stats,periodic,compute_mu)
  end subroutine get_sn_sphere_stats_many

  !TO ADDED: Packed, increasing candidates for multiple centers. Event e owns
  ! radii(first(e):first(e+1)-1); equal offsets represent a skipped/empty event.
  ! All ranks pass identical centers, radii, first, periodic and compute_mu.
  ! One owned-leaf traversal and one ALLREDUCE; prefix sums reset at each event.
  subroutine get_sn_sphere_stats_events(centers, radii, first, state, stats, periodic, compute_mu)
    real(dp), intent(in) :: centers(:,:), radii(:), state(:,:)
    integer, intent(in) :: first(:)
    type(sphere_stats), intent(out) :: stats(size(radii))
    logical, intent(in), optional :: periodic(ndim), compute_mu
    logical :: with_mu, have_weights, wrap(ndim)
    integer :: ilevel, igrid, icell, ind, ig, idim, ir, lo, hi, mid, event, begin, last
    real(dp) :: scale, dx, cell_volume, distance(ndim), offset(ndim), radius, separation, position(ndim)
    real(dp) :: skip(3), domain_size(3), weights(4)
    real(dp), allocatable :: local_sum(:,:), global_sum(:,:)
#ifndef WITHOUTMPI
    integer :: ierr, real_type
#endif

    stats = sphere_stats()
    if(size(centers,1)/=ndim.or.size(first)/=size(centers,2)+1) &
      call invalid_input('event center/offset extent mismatch')
    if(first(1)/=1.or.first(size(first))/=size(radii)+1)call invalid_input('invalid event candidate offsets')
    if(any(first<1).or.any(first>size(radii)+1))call invalid_input('event candidate offset out of bounds')
    if(any(first(2:)<first(:size(first)-1)))call invalid_input('event candidate offsets must be nondecreasing')
    do event=1,size(centers,2)
      call check_radii(radii(first(event):first(event+1)-1))
    enddo
    if (size(radii) == 0) return
    with_mu = .false.
    if (present(compute_mu)) with_mu = compute_mu
    if (with_mu) then
#ifdef RT
      call check_mu_columns(size(state,2))
#else
      call invalid_input('compute_mu requires an RT build with chemical species')
#endif
    end if
    if (.not. all(ieee_is_finite(centers))) call invalid_input('center must be finite')
    if (size(state,1) < ncoarse+twotondim*ngridmax .or. size(state,2) < 1) &
      call invalid_input('state has insufficient cell or variable extent')
    if (metal) then
      if (imetal < 1 .or. imetal > size(state,2)) call invalid_input('invalid metal column')
    end if
    if(size(radii)>shiftr(huge(1),2))call invalid_input('too many packed candidate statistics for MPI count')
    allocate(local_sum(4,size(radii)),global_sum(4,size(radii)))

    wrap = .false.
    if (present(periodic)) wrap = periodic
    skip = real((/icoarse_min, jcoarse_min, kcoarse_min/), dp)
    domain_size = real((/icoarse_max-icoarse_min+1, jcoarse_max-jcoarse_min+1, &
      kcoarse_max-kcoarse_min+1/), dp)
    if (boxlen <= 0.0_dp .or. any(domain_size(1:ndim) <= 0.0_dp)) &
      call invalid_input('invalid domain dimensions')
    scale = boxlen/domain_size(1)
    domain_size = domain_size*scale

    local_sum = 0.0_dp
    do ilevel = levelmin, nlevelmax
      dx = 0.5_dp**ilevel
      cell_volume = (dx*scale)**ndim
      do ind = 1, twotondim
        do idim = 1, ndim
          offset(idim) = (real(ibits(ind-1,idim-1,1),dp)-0.5_dp)*dx
        end do
        do ig = 1, active(ilevel)%ngrid
          igrid = active(ilevel)%igrid(ig)
          icell = ncoarse+(ind-1)*ngridmax+igrid
          if (son(icell) /= 0) cycle
          position=(xg(igrid,1:ndim)+offset-skip(1:ndim))*scale
          have_weights=.false.
          do event=1,size(centers,2)
            begin=first(event)
            last=first(event+1)-1
            if(begin>last)cycle
            radius=radii(last)
            distance=position-centers(:,event)
            do idim = 1, ndim
              if (wrap(idim)) distance(idim) = distance(idim) &
                - anint(distance(idim)/domain_size(idim))*domain_size(idim)
            end do
            if (any(abs(distance) > radius)) cycle
            separation = norm2(distance)
            if (separation > radius) cycle
            lo = begin
            hi = last
            do while (lo < hi)
              mid = lo+(hi-lo)/2
              if (separation <= radii(mid)) then
                hi = mid
              else
                lo = mid+1
              end if
            end do
            !TO ADDED: Compute the cell's conserved weights once, even when
            ! several SN spheres overlap. Do not read chemistry outside all spheres.
            if(.not.have_weights)then
              weights=0.0_dp
              weights(1)=state(icell,1)*cell_volume
              if(metal)weights(2)=state(icell,imetal)*cell_volume
              weights(3)=cell_volume
#ifdef RT
              if(with_mu)weights(4)=cell_mass_over_mu(state,icell)*cell_volume
#endif
              have_weights=.true.
            endif
            local_sum(:,lo)=local_sum(:,lo)+weights
          enddo
        end do
      end do
    end do

#ifndef WITHOUTMPI
    call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL, storage_size(0.0_dp)/8, real_type, ierr)
    if (ierr /= MPI_SUCCESS) call invalid_input('no matching MPI real datatype')
    call MPI_ALLREDUCE(local_sum, global_sum, size(local_sum), real_type, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) call invalid_input('MPI_ALLREDUCE failed')
#else
    global_sum = local_sum
#endif
    do event=1,size(centers,2)
      do ir=first(event),first(event+1)-1
        if(ir>first(event))global_sum(:,ir)=global_sum(:,ir)+global_sum(:,ir-1)
        stats(ir)%mass = global_sum(1,ir)
        stats(ir)%metal_mass = global_sum(2,ir)
        stats(ir)%volume = global_sum(3,ir)
        if (stats(ir)%volume > 0.0_dp) stats(ir)%mean_density = stats(ir)%mass/stats(ir)%volume
        if (stats(ir)%mass > 0.0_dp) stats(ir)%mean_metallicity = stats(ir)%metal_mass/stats(ir)%mass
        stats(ir)%mass_over_mu = global_sum(4,ir)
        if (with_mu .and. stats(ir)%mass > 0.0_dp .and. stats(ir)%mass_over_mu > 0.0_dp) then
          stats(ir)%mean_mu = stats(ir)%mass/stats(ir)%mass_over_mu
          stats(ir)%mu_available = .true.
        end if
      enddo
    enddo
  end subroutine get_sn_sphere_stats_events

  ! TO ADDED: All ranks supply the SAME host stats, host dx, candidates, units,
  ! energy and center. Host stats must come from its owner's current leaf state,
  ! not a sphere or a possibly stale ghost copy. Units come from units.f90.
  ! This routine only selects a target; it does not modify any hydro variable.
  subroutine find_sn_feedback_radius_candidates(center, radii, state, host, host_dx, energy_erg, &
      scale_l, scale_d, scale_nH, choice, periodic)
    real(dp), intent(in) :: center(ndim), radii(:), state(:,:), host_dx, energy_erg
    real(dp), intent(in) :: scale_l, scale_d, scale_nH
    type(sphere_stats), intent(in) :: host
    type(sn_radius_choice), intent(out) :: choice
    logical, intent(in), optional :: periodic(ndim)
    type(sphere_stats) :: spheres(size(radii)), empty(0)

    call check_radii(radii)
    ! Decide the unresolved-host shortcut before any sphere traversal/collective.
    call select_sn_feedback_radius(host, host_dx, radii(:0), empty, energy_erg, &
      scale_l, scale_d, scale_nH, choice)
    if (choice%reason == sn_radius_unresolved .or. size(radii) == 0) return
    call get_sn_sphere_stats_many(center, radii, state, spheres, periodic, compute_mu=.true.)
    call select_sn_feedback_radius(host, host_dx, radii, spheres, energy_erg, &
      scale_l, scale_d, scale_nH, choice)
  end subroutine find_sn_feedback_radius_candidates

  ! TO ADDED: Generate increasing candidates from the host half-diagonal to
  ! n_radius*host_dx, using a geometric ratio (default 1.2), including r_max.
  ! n_radius is a positive REAL multiplier, not an AMR level or a cell count.
  function make_sn_candidate_radii(host_dx, n_radius, radius_ratio) result(radii)
    real(dp), intent(in) :: host_dx, n_radius
    real(dp), intent(in), optional :: radius_ratio
    real(dp), allocatable :: radii(:), work(:)
    real(dp) :: ratio, rmin, rmax, steps, next_radius
    integer :: capacity, count
    ratio = 1.2_dp
    if (present(radius_ratio)) ratio = radius_ratio
    if (.not. all(ieee_is_finite([host_dx,n_radius,ratio]))) &
      call invalid_input('nonfinite candidate radius parameters')
    if (host_dx <= 0.0_dp .or. n_radius <= 0.0_dp .or. ratio <= 1.0_dp) &
      call invalid_input('candidate generation requires dx>0, n_radius>0, ratio>1')
    rmin = sqrt(3.0_dp)*host_dx/2.0_dp
    rmax = n_radius*host_dx
    if (.not. all(ieee_is_finite([rmin,rmax]))) call invalid_input('candidate radius overflow')
    if (rmax < rmin) then
      allocate(radii(0))
      return
    end if
    if (rmax == rmin) then
      radii = [rmin]
      return
    end if
    steps = (log(rmax)-log(rmin))/log(ratio)
    if (steps >= real(huge(capacity)-3,dp)) call invalid_input('too many candidate radii')
    capacity = max(3,ceiling(steps)+3)
    allocate(work(capacity))
    count = 1
    work(1) = rmin
    do while (work(count) < rmax/ratio)
      next_radius = work(count)*ratio
      if (next_radius <= work(count)) call invalid_input('radius ratio cannot advance the candidate')
      if (next_radius >= rmax) exit
      if (count >= capacity-1) call invalid_input('candidate radius count exceeds allocation')
      count = count+1
      work(count) = next_radius
    end do
    radii = [work(:count),rmax]
  end function make_sn_candidate_radii

  ! TO ADDED: Convenience overload; all ranks must use the same n_radius/ratio.
  subroutine find_sn_feedback_radius_ncell(center, n_radius, state, host, host_dx, energy_erg, &
      scale_l, scale_d, scale_nH, choice, periodic, radius_ratio)
    real(dp), intent(in) :: center(ndim), n_radius, state(:,:), host_dx, energy_erg
    real(dp), intent(in) :: scale_l, scale_d, scale_nH
    type(sphere_stats), intent(in) :: host
    type(sn_radius_choice), intent(out) :: choice
    logical, intent(in), optional :: periodic(ndim)
    real(dp), intent(in), optional :: radius_ratio
    real(dp), allocatable :: radii(:)
    radii = make_sn_candidate_radii(host_dx,n_radius,radius_ratio)
    call find_sn_feedback_radius_candidates(center,radii,state,host,host_dx,energy_erg, &
      scale_l,scale_d,scale_nH,choice,periodic)
  end subroutine find_sn_feedback_radius_ncell

  ! TO ADDED: Default entry point uses Ncell_SN_max from &PHYSICS_PARAMS.
  subroutine find_sn_feedback_radius_configured(center, state, host, host_dx, energy_erg, &
      scale_l, scale_d, scale_nH, choice, periodic, radius_ratio)
    real(dp), intent(in) :: center(ndim), state(:,:), host_dx, energy_erg
    real(dp), intent(in) :: scale_l, scale_d, scale_nH
    type(sphere_stats), intent(in) :: host
    type(sn_radius_choice), intent(out) :: choice
    logical, intent(in), optional :: periodic(ndim)
    real(dp), intent(in), optional :: radius_ratio
    call find_sn_feedback_radius_ncell(center,Ncell_SN_max,state,host,host_dx,energy_erg, &
      scale_l,scale_d,scale_nH,choice,periodic,radius_ratio)
  end subroutine find_sn_feedback_radius_configured

  ! TO ADDED: Local decision on precomputed statistics (no MPI communication).
  ! Candidates start at or above the host half-diagonal and are increasing.
  ! Host shortcut: rcool_host < sqrt(3)*dx/2. Equality proceeds to the search.
  ! TO MODIFIED: Largest candidate with r < rcool(r)/2; DeltaT no longer gates
  ! acceptance and is recorded for diagnostics only. The largest admissible
  ! sphere minimizes the post-injection sound speed (Courant limit), which is
  ! the SN_reservoir goal. Unresolved hosts and a failed search return
  ! single_cell=true.
  subroutine select_sn_feedback_radius(host, host_dx, radii, spheres, energy_erg, &
      scale_l, scale_d, scale_nH, choice)
    type(sphere_stats), intent(in) :: host, spheres(:)
    real(dp), intent(in) :: host_dx, radii(:), energy_erg, scale_l, scale_d, scale_nH
    type(sn_radius_choice), intent(out) :: choice
    real(dp) :: rcell, rcool_pc, delta_t, length_to_pc, mass_to_g
    integer :: ir

    if (ndim /= 3) call invalid_input('SN radius selection requires NDIM=3')
    call check_radii(radii)
    if (size(spheres) /= size(radii)) call invalid_input('candidate statistics extent mismatch')
    if (.not. all(ieee_is_finite([host_dx,energy_erg,scale_l,scale_d,scale_nH,gamma]))) &
      call invalid_input('nonfinite SN radius selection parameters')
    if (min(host_dx,energy_erg,scale_l,scale_d,scale_nH) <= 0.0_dp .or. gamma <= 1.0_dp) &
      call invalid_input('SN radius selection requires positive units, dx, energy and gamma>1')
    rcell = sqrt(3.0_dp)*host_dx/2.0_dp
    if (any(radii < rcell)) call invalid_input('candidate radius is smaller than the host half-diagonal')
    length_to_pc = scale_l/pc_in_cm
    mass_to_g = scale_d*scale_l**3
    choice = sn_radius_choice()
    choice%radius = rcell
    choice%stats = host
    choice%delta_temperature = temperature_increase(host,energy_erg,mass_to_g)
    if (host%mean_density > 0.0_dp) then
      choice%cooling_radius_pc = sn_cooling_radius_pc(host%mean_density*scale_nH, &
        energy_erg,host%mean_metallicity)
      if (choice%cooling_radius_pc < rcell*length_to_pc) then
        choice%reason = sn_radius_unresolved
        return
      end if
    end if
    ! TO MODIFIED: Scan from largest to smallest; skip failures without assuming
    ! monotonicity. The first success is the largest admissible candidate.
    do ir = size(radii), 1, -1
      if (spheres(ir)%mass <= 0.0_dp .or. spheres(ir)%volume <= 0.0_dp) cycle
      if (.not. spheres(ir)%mu_available) call invalid_input('candidate has no chemical mu')
      rcool_pc = sn_cooling_radius_pc(spheres(ir)%mean_density*scale_nH,energy_erg, &
        spheres(ir)%mean_metallicity)
      if (radii(ir)*length_to_pc >= rcool_pc/2.0_dp) cycle
      delta_t = temperature_increase(spheres(ir),energy_erg,mass_to_g) ! TO MODIFIED: Diagnostic only; no longer gates acceptance.
      choice%single_cell = .false.
      choice%reason = sn_radius_sphere
      choice%candidate_index = ir
      choice%radius = radii(ir)
      choice%cooling_radius_pc = rcool_pc
      choice%delta_temperature = delta_t
      choice%stats = spheres(ir)
      return
    end do
  end subroutine select_sn_feedback_radius

  ! TO ADDED: Particle-weighted mean temperature INCREASE, at fixed chemistry.
  ! A huge value marks undefined temperature for empty/no-chemistry host stats.
  function temperature_increase(stats, energy_erg, mass_to_g) result(delta_t)
    type(sphere_stats), intent(in) :: stats
    real(dp), intent(in) :: energy_erg, mass_to_g
    real(dp) :: delta_t
    delta_t = huge(0.0_dp)
    if (.not. stats%mu_available .or. stats%mass_over_mu <= 0.0_dp) return
    delta_t = ((gamma-1.0_dp)*mH/kB)*(energy_erg/mass_to_g)/stats%mass_over_mu
  end function temperature_increase

  ! TO ADDED: Keep validation identical for scalar and batched sphere queries.
  subroutine check_radii(radii)
    real(dp), intent(in) :: radii(:)
    if (.not. all(ieee_is_finite(radii))) call invalid_input('radius must be finite')
    if (any(radii < 0.0_dp)) call invalid_input('radius must be nonnegative')
    if (size(radii) > 1) then
      if (any(radii(2:) <= radii(:size(radii)-1))) call invalid_input('radii must be strictly increasing')
    end if
  end subroutine check_radii

#ifdef RT
  ! TO ADDED: Validate indices once per sphere, before reading any species.
  subroutine check_mu_columns(ncolumns)
    integer, intent(in) :: ncolumns
    if (iIons < 1) call invalid_input('invalid starting chemical-species column')
    call check_index(ixHII)
    if (isH2) call check_index(ixHI)
    if (isHe) then
      call check_index(ixHeII)
      call check_index(ixHeIII)
    end if
    if (.not. ieee_is_finite(X) .or. .not. ieee_is_finite(Y)) &
      call invalid_input('nonfinite H/He mass fractions')
    if (X <= 0.0_dp .or. Y < 0.0_dp) call invalid_input('invalid H/He mass fractions')
  contains
    subroutine check_index(index)
      integer, intent(in) :: index
      if (index < 1 .or. index > nIons .or. iIons+index-1 > ncolumns) &
        call invalid_input('missing chemical-species column for compute_mu')
    end subroutine check_index
  end subroutine check_mu_columns

  ! TO ADDED: Algebra matches rt_cooling_module:getMu, including H2/electrons.
  ! Metals do not contribute to particle counts in that existing EOS.
  function cell_mass_over_mu(state, icell) result(rho_over_mu)
    real(dp), intent(in) :: state(:,:)
    integer, intent(in) :: icell
    real(dp) :: rho_over_mu, rho, xHI, xHII, xHeII, xHeIII
    real(dp), parameter :: fraction_tolerance = 1.0e-6_dp
    rho = state(icell,1)
    if (.not. ieee_is_finite(rho) .or. rho < 0.0_dp) &
      call invalid_input('compute_mu requires nonnegative finite density')
    rho_over_mu = 0.0_dp
    if (rho == 0.0_dp) return
    xHII = species_fraction(ixHII)
    xHI = 1.0_dp-xHII
    if (isH2) xHI = species_fraction(ixHI)
    xHeII = 0.0_dp
    xHeIII = 0.0_dp
    if (isHe) then
      xHeII = species_fraction(ixHeII)
      xHeIII = species_fraction(ixHeIII)
    end if
    if (xHI+xHII > 1.0_dp+fraction_tolerance .or. xHeII+xHeIII > 1.0_dp+fraction_tolerance) &
      call invalid_input('chemical fractions sum to more than one')
    ! Permit rounding drift, without silently accepting invalid compositions.
    xHI = min(xHI,1.0_dp-xHII)
    xHeIII = min(xHeIII,1.0_dp-xHeII)
    rho_over_mu = rho*(X*(xHI+2.0_dp*xHII+0.5_dp*(1.0_dp-xHI-xHII)) &
      + 0.25_dp*Y*(1.0_dp+xHeII+2.0_dp*xHeIII))
  contains
    function species_fraction(index) result(value)
      integer, intent(in) :: index
      real(dp) :: value
      value = state(icell,iIons+index-1)/rho
      if (.not. ieee_is_finite(value)) call invalid_input('nonfinite chemical fraction')
      if (value < -fraction_tolerance .or. value > 1.0_dp+fraction_tolerance) &
        call invalid_input('chemical fraction outside [0,1]')
      value = max(0.0_dp,min(1.0_dp,value))
    end function species_fraction
  end function cell_mass_over_mu
#endif

  subroutine invalid_input(message)
    use, intrinsic :: iso_fortran_env, only: error_unit ! TO ADDED: Avoid recursive stdout I/O.
    character(len=*), intent(in) :: message
#ifndef WITHOUTMPI
    integer :: ierr
#endif
    write(error_unit,'(a)') 'sn_sphere_stats: '//message ! TO MODIFIED: Shared module diagnostic.
#ifndef WITHOUTMPI
    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#endif
    error stop 1
  end subroutine invalid_input
end module sn_sphere_stats
