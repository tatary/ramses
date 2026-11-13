! merge_ks (2026-11-13): ported from patch/ks (use mpi_mod). Hooks: amr_step (prepare), refine_utils
! (refine/coarsen/ghosts), load_balance (ghosts, defrag reorder), init_hydro/output_hydro (checkpoint
! trailer), feedback (deposit). ivar_coolstop/fcool_stop live in amr_parameters here.
! TO ADDED: Pending SN thermal energy, separate from advected hydro variables.
module sn_energy_reservoir
  use amr_commons
  use hydro_commons
  use sn_sphere_stats, only: sn_radius_choice
  use sn_feedback_radius, only: sn_radius_periodic_axes
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: sn_pending, sn_pending_coolstop, init_sn_reservoir, deposit_sn_reservoir, apply_sn_reservoir
  public :: prepare_sn_reservoir, refine_sn_reservoir, coarsen_sn_reservoir
  public :: reorder_sn_reservoir, write_sn_reservoir, read_sn_reservoir
  ! Physical energy [erg] per CODE cell volume. Only owned leaves enter the budget.
  ! Holding erg fixed avoids changing un-injected energy when aexp/units change.
  ! No pressure, cooling or advection until consumed.
  real(dp), allocatable :: sn_pending(:)
  ! TO ADDED: 0/1 request to set a timer when pending heat is applied; new events mark only their host.
  ! Stored as real for existing AMR/MPI field transfers. It is neither energy nor an advected timer.
  real(dp), allocatable :: sn_pending_coolstop(:)
  ! TO ADDED: Globally identical conservative flag; no collective is needed.
  ! A scheduled levelmin subtree flush consumes every owned pending leaf.
  logical :: pending_possible=.false.
  ! merge_ks: diagnostics, released energy [erg] and cells since the last report (rank local)
  real(dp) :: released_erg=0.0_dp
  integer :: released_cells=0
contains
  subroutine init_sn_reservoir()
    if(.not.SN_reservoir)return
    if(allocated(sn_pending))return
    allocate(sn_pending(ncoarse+twotondim*ngridmax))
    allocate(sn_pending_coolstop(size(sn_pending))) !TO ADDED
    sn_pending=0.0_dp
    sn_pending_coolstop=0.0_dp !TO ADDED
    pending_possible=.true.
  end subroutine init_sn_reservoir

  ! TO MODIFIED: All ranks supply the same event list/decisions and energies in ERG.
  ! Volume weighting gives selected cells the same energy-density increment. Only owners
  ! write the reservoir; no reverse ghost accumulation is needed.
  subroutine deposit_sn_reservoir(positions, energies, choices)
    real(dp), intent(in) :: positions(:,:), energies(:)
    type(sn_radius_choice), intent(in) :: choices(:)
    real(dp) :: origin(ndim), extent(ndim), scale, dx, width, volume
    real(dp) :: offset(ndim), lower(ndim), center(ndim), point(ndim), distance(ndim)
    integer :: lev, ind, idim, ig, cell, event
    logical :: periodic(ndim), selected, host !TO MODIFIED
    call init_sn_reservoir()
    if(.not.allocated(sn_pending))return
    if(size(positions,1)/=ndim.or.size(positions,2)/=size(energies).or.size(choices)/=size(energies)) &
      call fail('event array extent mismatch')
    if(any(energies>0.0_dp))pending_possible=.true. !TO ADDED: Same global events on every rank.
    periodic=sn_radius_periodic_axes()
    origin=real([icoarse_min,jcoarse_min,kcoarse_min],dp)
    extent=real([icoarse_max-icoarse_min+1,jcoarse_max-jcoarse_min+1,kcoarse_max-kcoarse_min+1],dp)
    scale=boxlen/extent(1)
    extent=extent*scale
    do lev=levelmin,nlevelmax
      dx=0.5_dp**lev
      width=dx*scale
      volume=width**ndim
      do ind=1,twotondim
        do idim=1,ndim
          offset(idim)=(real(ibits(ind-1,idim-1,1),dp)-1.0_dp)*dx
        enddo
        do ig=1,active(lev)%ngrid
          cell=ncoarse+(ind-1)*ngridmax+active(lev)%igrid(ig)
          if(son(cell)/=0)cycle
          lower=(xg(active(lev)%igrid(ig),:)+offset-origin)*scale
          center=(xg(active(lev)%igrid(ig),:)+(offset+dx/2.0_dp)-origin)*scale !TO ADDED: Match sphere selection arithmetic.
          do event=1,size(energies)
            if(energies(event)==0.0_dp)cycle
            if(.not.ieee_is_finite(energies(event)).or.energies(event)<0.0_dp)call fail('invalid SN energy')
            point=positions(:,event)
            ! TO ADDED: The host is the leaf containing the event, including periodic wrapping.
            do idim=1,ndim
              if(periodic(idim))point(idim)=modulo(point(idim),extent(idim))
            enddo
            host=all(point>=lower).and.all(point<lower+width)
            if(choices(event)%single_cell)then
              selected=host !TO MODIFIED
              if(selected)sn_pending(cell)=sn_pending(cell)+energies(event)/volume
            else
              distance=center-positions(:,event) !TO MODIFIED: Keep sphere-statistics selection arithmetic.
              do idim=1,ndim
                if(periodic(idim))distance(idim)=distance(idim)-anint(distance(idim)/extent(idim))*extent(idim)
              enddo
              if(norm2(distance)>choices(event)%radius)cycle
              ! TO MODIFIED: E_cell = E_SN * V_cell / sum(V_cell), using selected leaf volumes.
              if(.not.ieee_is_finite(choices(event)%stats%volume).or.choices(event)%stats%volume<=0.0_dp) &
                call fail('sphere has invalid injection volume')
              sn_pending(cell)=sn_pending(cell)+energies(event)/choices(event)%stats%volume
            endif
            ! TO ADDED: Other heated cells retain any old timer but receive no new request.
            if(host.and.fcool_stop>0.0_dp)sn_pending_coolstop(cell)=1.0_dp
          enddo
        enddo
      enddo
    enddo
  end subroutine deposit_sn_reservoir

  ! Consume owned leaves exactly once. Called before this level's dt/set_unew.
  subroutine apply_sn_reservoir(ilevel)
    integer, intent(in) :: ilevel
    integer :: ind, ig, cell
    real(dp) :: rho_cell, eth, cs, dx, sl, st, sd, sv, snh, st2, energy_unit
    if(.not.allocated(sn_pending))return
    dx=0.5_dp**ilevel*boxlen/real(icoarse_max-icoarse_min+1,dp)
    call units(sl,st,sd,sv,snh,st2)
    energy_unit=sd*sl**3*sv**2
    do ind=1,twotondim
      do ig=1,active(ilevel)%ngrid
        cell=ncoarse+(ind-1)*ngridmax+active(ilevel)%igrid(ig)
        if(son(cell)/=0)cycle
        if(sn_pending(cell)==0.0_dp)cycle
        if(.not.ieee_is_finite(sn_pending(cell)).or.sn_pending(cell)<0.0_dp)call fail('invalid pending energy')
        uold(cell,ndim+2)=uold(cell,ndim+2)+sn_pending(cell)/energy_unit
        released_erg=released_erg+sn_pending(cell)*(dx**ndim) !merge_ks: diagnostics (erg per code volume x cell volume)
        released_cells=released_cells+1
        sn_pending(cell)=0.0_dp
        if(fcool_stop>0.0_dp.and.sn_pending_coolstop(cell)>0.0_dp)then !TO MODIFIED: Host only.
          rho_cell=max(uold(cell,1),smallr)
          eth=uold(cell,ndim+2)-0.5_dp*sum(uold(cell,2:ndim+1)**2)/rho_cell
          eth=eth-sum(uold(cell,ndim+3:ndim+2+nener))
          cs=sqrt(max(gamma*(gamma-1.0_dp)*eth/rho_cell,smallc**2))
          uold(cell,ivar_coolstop)=max(uold(cell,ivar_coolstop),fcool_stop*dx/cs*st/3.15569e13_dp)
        endif
        sn_pending_coolstop(cell)=0.0_dp !TO ADDED: Consume the request exactly once.
      enddo
    enddo
  end subroutine apply_sn_reservoir

  ! The scheduled subtree is about to advance. Flush before refinement and dt,
  ! descending so restriction/ghosts include heated finer cells. Ancestors with
  ! an already fixed dt are left untouched (including their pending energy).
  subroutine prepare_sn_reservoir(ilevel)
    integer, intent(in) :: ilevel
    integer :: lev
    if(.not.SN_reservoir.or..not.hydro)return
    call init_sn_reservoir()
    if(.not.pending_possible)return !TO ADDED: No cell traversal or ghost communication between events.
    do lev=nlevelmax,ilevel,-1
      if(lev<nlevelmax)call upload_fine(lev)
      call apply_sn_reservoir(lev)
      call make_virtual_fine_dp(uold(1,ndim+2),lev)
      if(fcool_stop>0.0_dp)call make_virtual_fine_dp(uold(1,ivar_coolstop),lev)
      if(simple_boundary)call make_boundary_hydro(lev)
    enddo
    if(ilevel==levelmin)then
      pending_possible=.false. !TO ADDED: All pending levels have been consumed.
      call report_sn_reservoir !merge_ks
    endif
  end subroutine prepare_sn_reservoir

  ! merge_ks: one line per coarse step on rank 1 (energy released into the gas since the last report)
  subroutine report_sn_reservoir
    use mpi_mod
    real(dp) :: buf(2), tot(2)
    integer :: info
    buf(1)=released_erg; buf(2)=real(released_cells,dp)
    tot=buf
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(buf,tot,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
    if(myid==1.and.tot(2)>0.0_dp)write(*,'(A,1PE12.4,A,I8,A)')' SN reservoir: released ',tot(1),' erg into ',nint(tot(2)),' cells this coarse step'
    released_erg=0.0_dp; released_cells=0
  end subroutine report_sn_reservoir

  ! Piecewise constant prolongation of energy density conserves sum(e*V).
  ! Refined parents are cleared: reservoir is a leaf field, not a second copy.
  subroutine refine_sn_reservoir(parents, grids)
    integer, intent(in) :: parents(:), grids(:)
    integer :: ind
    if(.not.allocated(sn_pending))return
    do ind=1,twotondim
      sn_pending(ncoarse+(ind-1)*ngridmax+grids)=sn_pending(parents)
      ! TO ADDED: If a pending host is refined, its children inherit that host region's request.
      sn_pending_coolstop(ncoarse+(ind-1)*ngridmax+grids)=sn_pending_coolstop(parents)
    enddo
    sn_pending(parents)=0.0_dp
    sn_pending_coolstop(parents)=0.0_dp !TO ADDED
  end subroutine refine_sn_reservoir

  subroutine coarsen_sn_reservoir(parents, grids)
    integer, intent(in) :: parents(:), grids(:)
    integer :: ind
    if(.not.allocated(sn_pending))return
    sn_pending(parents)=0.0_dp
    sn_pending_coolstop(parents)=0.0_dp !TO ADDED
    do ind=1,twotondim
      sn_pending(parents)=sn_pending(parents)+sn_pending(ncoarse+(ind-1)*ngridmax+grids)/real(twotondim,dp)
      sn_pending(ncoarse+(ind-1)*ngridmax+grids)=0.0_dp
      ! TO ADDED: A merged cell retains a request from any of its children.
      sn_pending_coolstop(parents)=max(sn_pending_coolstop(parents), &
        sn_pending_coolstop(ncoarse+(ind-1)*ngridmax+grids))
      sn_pending_coolstop(ncoarse+(ind-1)*ngridmax+grids)=0.0_dp
    enddo
  end subroutine coarsen_sn_reservoir

  ! Same old-grid list/order as defrag; copy before any grid indices are reused.
  subroutine reorder_sn_reservoir(old_grids)
    integer, intent(in) :: old_grids(:)
    real(dp), allocatable :: values(:)
    integer :: ind, offset
    if(.not.allocated(sn_pending))return
    allocate(values(ngridmax))
    do ind=1,twotondim
      offset=ncoarse+(ind-1)*ngridmax
      values=0.0_dp
      values(:size(old_grids))=sn_pending(offset+old_grids)
      sn_pending(offset+1:offset+ngridmax)=values
      values=0.0_dp !TO ADDED: Reorder the request using the same old indices.
      values(:size(old_grids))=sn_pending_coolstop(offset+old_grids)
      sn_pending_coolstop(offset+1:offset+ngridmax)=values
    enddo
  end subroutine reorder_sn_reservoir

  ! Versioned trailer in the hydro checkpoint, after all legacy records.
  ! Grid traversal matches hydro I/O, so restart does not depend on old indices.
  subroutine write_sn_reservoir(unit)
    integer, intent(in) :: unit
    integer :: lev, ib, count, ind
    integer, allocatable :: grids(:)
    ! merge_ks: no trailer at all unless the reservoir is in use, so that the hydro files of
    ! runs without SN_reservoir stay identical to before (read_sn_reservoir treats a missing
    ! trailer as an old checkpoint)
    if(.not.allocated(sn_pending))return
    write(unit)20260907,2,allocated(sn_pending) !TO MODIFIED: v2 adds pending host timer requests.
    do lev=1,nlevelmax
      do ib=1,ncpu+nboundary
        call grid_list(lev,ib,grids)
        count=size(grids)
        write(unit)lev,count
        do ind=1,twotondim
          if(count>0)write(unit)sn_pending(ncoarse+(ind-1)*ngridmax+grids)
          if(count>0)write(unit)sn_pending_coolstop(ncoarse+(ind-1)*ngridmax+grids) !TO ADDED
        enddo
      enddo
    enddo
  end subroutine write_sn_reservoir

  subroutine read_sn_reservoir(unit, levels)
    integer, intent(in) :: unit, levels
    integer :: ios, magic, version, lev, ib, lev2, count, ind
    integer, allocatable :: grids(:)
    real(dp), allocatable :: values(:)
    logical :: present_field
    call init_sn_reservoir()
    pending_possible=.true. !TO ADDED: Restart must recheck all levels even when reusing an allocation.
    if(allocated(sn_pending))sn_pending=0.0_dp
    if(allocated(sn_pending_coolstop))sn_pending_coolstop=0.0_dp !TO ADDED
    read(unit,iostat=ios)magic,version,present_field
    if(is_iostat_end(ios))return ! Old checkpoint: no pending energy field.
    if(ios/=0)call fail('cannot read reservoir checkpoint header')
    if(magic/=20260907.or.(version/=1.and.version/=2))call fail('unknown reservoir checkpoint format') !TO MODIFIED
    if(.not.present_field)return
    do lev=1,levels
      do ib=1,ncpu+nboundary
        call grid_list(lev,ib,grids)
        read(unit,iostat=ios)lev2,count
        if(ios/=0)call fail('truncated reservoir checkpoint')
        if(lev2/=lev.or.count/=size(grids))call fail('reservoir checkpoint grid mismatch')
        allocate(values(count))
        do ind=1,twotondim
          if(count==0)cycle
          read(unit,iostat=ios)values
          if(ios/=0)call fail('truncated reservoir checkpoint cells')
          if(.not.all(ieee_is_finite(values)).or.any(values<0.0_dp))call fail('invalid saved pending energy')
          if(.not.SN_reservoir.and.any(values/=0.0_dp))call fail('enable SN_reservoir to restore pending energy')
          if(allocated(sn_pending))sn_pending(ncoarse+(ind-1)*ngridmax+grids)=values
          if(version==2)then !TO ADDED
            read(unit,iostat=ios)values
            if(ios/=0)call fail('truncated reservoir checkpoint timer requests')
            if(any(values/=0.0_dp.and.values/=1.0_dp))call fail('invalid saved timer request')
            if(.not.SN_reservoir.and.any(values/=0.0_dp))call fail('enable SN_reservoir to restore timer requests')
          else
            ! TO ADDED: v1 lacks host information. Preserve its original all-receiver timer semantics
            ! for already pending heat; events deposited after restart use host-only requests.
            values=merge(1.0_dp,0.0_dp,values>0.0_dp)
          endif
          if(allocated(sn_pending_coolstop))sn_pending_coolstop(ncoarse+(ind-1)*ngridmax+grids)=values
        enddo
        deallocate(values)
      enddo
    enddo
  end subroutine read_sn_reservoir

  subroutine grid_list(lev, ib, grids)
    integer, intent(in) :: lev, ib
    integer, allocatable, intent(out) :: grids(:)
    integer :: count, grid, j
    if(ib<=ncpu)then
      count=numbl(ib,lev); grid=headl(ib,lev)
    else
      count=numbb(ib-ncpu,lev); grid=headb(ib-ncpu,lev)
    endif
    allocate(grids(count))
    do j=1,count
      grids(j)=grid
      grid=next(grid)
    enddo
  end subroutine grid_list

  ! TO MODIFIED: A rank-local failure (e.g. invalid pending energy in apply) must
  ! abort every rank instead of waiting alone in clean_stop's MPI_FINALIZE.
  subroutine fail(message)
    use, intrinsic :: iso_fortran_env, only: error_unit
    use mpi_mod ! merge_ks
    character(len=*), intent(in) :: message
#ifndef WITHOUTMPI
    integer :: ierr
#endif
    write(error_unit,'(a)')'SN reservoir: '//message
    flush(error_unit)
#ifndef WITHOUTMPI
    call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
#endif
    error stop 1
  end subroutine fail
end module sn_energy_reservoir
