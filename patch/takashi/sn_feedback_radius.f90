! merge_ks (2026-11-13): ported from patch/ks unchanged.
! TO ADDED: RAMSES event adapter. Search/diagnostics only; no hydro updates.
module sn_feedback_radius
  use amr_commons, only: dp, ndim, myid, nboundary, simple_boundary, boundary_type, &
    nstep, nstep_coarse, Ncell_SN_max
  use sn_sphere_stats, only: sphere_stats, sn_radius_choice, get_sn_host_stats_many, &
    select_sn_feedback_radius, make_sn_candidate_radii, get_sn_sphere_stats_events, sn_radius_unresolved
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: evaluate_sn_feedback_radii, sn_radius_periodic_axes, report_sn_feedback_radii
  type :: event_candidates !TO ADDED: Candidate counts can differ between events.
    real(dp), allocatable :: radii(:)
  end type event_candidates
contains
  ! RAMSES wraps axes without physical boundaries. Boundary direction codes
  ! 1/2, 3/4, 5/6 denote X/Y/Z; tens encode reflecting/free/inflow conditions.
  function sn_radius_periodic_axes() result(periodic)
    logical :: periodic(ndim)
    integer :: ibound, direction
    periodic=.true.
    if(.not.simple_boundary)return
    do ibound=1,nboundary
      direction=mod(boundary_type(ibound),10)
      if(direction>=1 .and. direction<=2*ndim)periodic((direction+1)/2)=.false.
    enddo
  end function sn_radius_periodic_axes

  ! All event arrays are globally shared, in the same order on every rank.
  ! Energy is ESN alone in erg, excluding simultaneous winds/ejecta motion.
  ! state is the pre-injection uold snapshot at each level's own current time.
  subroutine evaluate_sn_feedback_radii(positions, energy_erg, state, scale_l, scale_d, scale_nH, choices)
    real(dp), intent(in) :: positions(:,:), energy_erg(:), state(:,:)
    real(dp), intent(in) :: scale_l, scale_d, scale_nH
    type(sn_radius_choice), intent(out) :: choices(size(energy_erg))
    type(sphere_stats) :: hosts(size(energy_erg))
    real(dp) :: widths(size(energy_erg))
    logical :: periodic(ndim)
    integer :: event, begin, last
    integer :: first(size(energy_erg)+1) !TO ADDED: Packed event boundaries.
    type(event_candidates) :: candidates(size(energy_erg))
    real(dp), allocatable :: radii(:)
    type(sphere_stats), allocatable :: spheres(:)
    type(sphere_stats) :: empty(0)
    choices=sn_radius_choice()
    if(size(positions,1)/=ndim .or. size(positions,2)/=size(energy_erg))then
      write(*,*)'SN radius event array extent mismatch'
      call clean_stop
    endif
    if(size(energy_erg)==0)return
    if(.not.all(ieee_is_finite(energy_erg)) .or. any(energy_erg<0.0_dp))then
      write(*,*)'SN radius event energy must be nonnegative and finite'
      call clean_stop
    endif
    if(all(energy_erg==0.0_dp))return
    periodic=sn_radius_periodic_axes()
    call get_sn_host_stats_many(positions,state,hosts,widths,periodic)
    !TO MODIFIED: Host lookup is collective 1. Precheck all hosts locally,
    ! pack only candidates that need a sphere search, then collective 2 handles
    ! every event/candidate together. Zero-energy and unresolved events retain order.
    first(1)=1
    do event=1,size(energy_erg)
      first(event+1)=first(event)
      if(energy_erg(event)==0.0_dp)cycle
      call select_sn_feedback_radius(hosts(event),widths(event),[real(dp)::],empty,energy_erg(event), &
        scale_l,scale_d,scale_nH,choices(event))
      if(choices(event)%reason==sn_radius_unresolved)cycle
      candidates(event)%radii=make_sn_candidate_radii(widths(event),Ncell_SN_max)
      if(size(candidates(event)%radii)>huge(1)-first(event))then
        write(*,*)'Too many SN radius candidates'
        call clean_stop
      endif
      first(event+1)=first(event)+size(candidates(event)%radii)
    enddo
    if(first(size(first))==1)return ! No sphere candidates: omit the second collective.
    allocate(radii(first(size(first))-1),spheres(first(size(first))-1))
    do event=1,size(energy_erg)
      begin=first(event); last=first(event+1)-1
      if(begin<=last)radii(begin:last)=candidates(event)%radii
    enddo
    call get_sn_sphere_stats_events(positions,radii,first,state,spheres,periodic,compute_mu=.true.)
    do event=1,size(energy_erg)
      begin=first(event); last=first(event+1)-1
      if(begin>last)cycle
      call select_sn_feedback_radius(hosts(event),widths(event),radii(begin:last),spheres(begin:last), &
        energy_erg(event),scale_l,scale_d,scale_nH,choices(event))
    enddo
  end subroutine evaluate_sn_feedback_radii

  ! One rank emits one record per event to the ordinary RAMSES run log.
  ! reason: 0=sphere, 1=unresolved host, 2=no admissible candidate.
  ! For single_cell=T, radius is a diagnostic half-diagonal, NOT an injection sphere.
  subroutine report_sn_feedback_radii(ilevel, ids, energy_erg, choices, scale_l, scale_nH)
    integer, intent(in) :: ilevel, ids(:)
    real(dp), intent(in) :: energy_erg(:), scale_l, scale_nH
    type(sn_radius_choice), intent(in) :: choices(:)
    integer :: event
    if(myid/=1)return
    do event=1,size(choices)
      if(energy_erg(event)==0.0_dp)cycle
      write(*,'(a,4(i0,1x),a,l1,a,i0,a,7(es16.8,1x))') &
        'SN_RADIUS coarse/fine/level/id=',nstep_coarse,nstep,ilevel,ids(event), &
        'single_cell=',choices(event)%single_cell,' reason=',choices(event)%reason, &
        ' r_pc/rcool_pc/dT_K/nH/Z/mu/E_erg=', &
        choices(event)%radius*scale_l/3.08568e18_dp,choices(event)%cooling_radius_pc, &
        choices(event)%delta_temperature,choices(event)%stats%mean_density*scale_nH, &
        choices(event)%stats%mean_metallicity,choices(event)%stats%mean_mu,energy_erg(event)
    enddo
  end subroutine report_sn_feedback_radii
end module sn_feedback_radius
