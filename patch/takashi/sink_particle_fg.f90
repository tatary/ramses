!################################################################
! sink_particle_fg: FG (patch/ks) sink particles / Pop III remnant black holes (sf_fg=.true.)
!
! Ported from patch/ks sink_particle.f90, init_sink.f90, output_sink.f90 and load_balance.f90
! (check_cloud_count) by scratchpad make_sink_fg.py (merge_ks, docs/merge_ks.md). Every routine
! carries the suffix _fg; the upstream pm/sink_particle.f90 (clump-based creation, FAM_CLOUD) is
! still linked and used when sf_fg=.false. amr_step / init_part / particle_tree / output_amr /
! newdt_fine / force_fine / read_params dispatch on sf_fg.
!
! ks lineage (older than upstream): cloud particles are identified by idp=-isink and, here, also
! carry the family FAM_CLOUD_FG (is_cloud is sf_fg-aware; the shared move_fine / synchro_fine /
! rho_fine kernels need it). Sinks are created from FAM_BH (Pop III remnant, bhpop3 and
! sink_mass_min_pop3BH) or FAM_DEAD (sink_mass_min_pop2BH) particles by make_sink_from_pop3_new_fg,
! which sets mp=0 on the particle. Accretion: BH_model 0 none, 1/3/4 Bondi-Hoyle-Lyttleton (cosmo
! only), 2 Sugimura & Ricotti 2020 (falls back to BHL); acc_sink_boost, eddington_limit; BH_model 3
! deposits 0.5% of the accreted rest mass as heat (sgagn_T2max cap). Luminosity (acc_lum): 0.1 mdot c^2
! (BH_model 1/3) or RIAF/standard/slim (2); the RT sources are in rt_spectra (set_BH_SED,
! BH_RT_feedback). Softening (sink_soft) and the cloud radius are comoving (KS MODIFIED); mergers:
! merge_smbh_sink (smbh) within 2 ssoft, merge_star_sink otherwise. Restart: sink_XXXXX.out<rank>
! (binary, all ranks identical, backup_sink_fg / init_sink_fg). Log: logsink<id> files and the
! stdout table of print_and_save_sink_properties_fg (each coarse step, from create_sink_fg).
!
! Removed from ks: clump-based creation (make_sink_from_clump, true_max), the old accrete_sink
! (weightp / count_clouds / sink_drag / AGN blast), test problems (is_test), sub-grid AGN seeding
! (sgagn_*, needs the halo finder), p3idpsink, nol_accretion/flux_accretion; restart_files -> output_restart.
! Files: every rank writes sink_XXXXX.out<rank> (identical content, as ks) and reads its own on restart; rank 1
! also writes sink_file_descriptor.txt. yt's RAMSES frontend (SinkParticleFileHandler) needs the file of every
! domain to exist (setup_handler raises otherwise) and reads the content of domain 1 only.
! check_sink_vars_among_cpus_fg (consistency ALLREDUCEs) only runs with sink_check_vars=.true.
!################################################################
subroutine create_sink_fg
  use amr_commons
  use pm_commons
  use hydro_commons
  use mpi_mod
  implicit none



  !----------------------------------------------------------------------------
  ! sink creation routine
  ! -remove all cloud particles, keep only global sink arrays
  ! -call clumpfinder for potential relevant peaks
  ! -flag peaks which are eligible for sink formation (flag 2)
  ! -new sink particle are created
  ! -new cloud particles are created
  ! -cloud particles are scattered to grid
  ! -accretion routine is called (with on_creation=.true.)
  !----------------------------------------------------------------------------

  integer::ilevel,ivar
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp):: scale_m !KS DEBUG



  integer :: info !KS DEBUG
  integer,save :: ifirst = 0 !KS ADDED (for prcedures at first call)


  if(verbose)write(*,*)' Entering create_sink_fg'

  !if (myid == 1) print *, "**** KS DEBUG **** create_sink_fg is called"


  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3  !KS DEBUG

  !KS ADDED
  if (ifirst == 0) then
    if (myid == 1) then
      write(*,'(A)') '##############################################################################'
      write(*,'(A)') 'Parameters for sink particle (first call of create_sink_fg) '
      write(*,'(A)') '##############################################################################'
      write(*,'(3(A, 1PE11.3))') &
        'sink_soft: ', sink_soft,&
        ', ssoft [cpc]: ',  sink_soft * 0.5D0**nlevelmax * scale_l/3.08568d18 / aexp,&
        ', r_merge [cpc]: ',  2*sink_soft * 0.5D0**nlevelmax * scale_l/3.08568d18 / aexp
      write(*,'(A,I4,A,1PE11.3,A,I4)') &
        'ir_cloud = ', ir_cloud, &
        ', mass_sink_direct_force = ', mass_sink_direct_force, &
        ', ir_cloud_massive = ', ir_cloud_massive
      write(*,*) 'softening length ssoft is defiened in comving length (KS MODIFIED)'
      write(*,*) 'radius for merger condition is set to 2 * ssoft (KS NOTE)'
      write(*,*) 'escape velocity check is not used as merger criteria (KS TODO)'

      write(*,'(A,1PE11.3,A)') 'sink_mass_min_pop3BH = ',  sink_mass_min_pop3BH, '  [Msun] (no Pop III sink BH if negative)'
      write(*,'(A,1PE11.3,A)') 'sink_mass_min_pop2BH = ',  sink_mass_min_pop2BH, '  [Msun] (no Pop II sink BH if negative)'
      if (no_sinkcloud_refine) write(*,'(A)') 'no_sinkcloud_refine = .true.  -> sink-cloud particles are not used when evaluationg refinement condition'
      write(*,'(A)') '##############################################################################'
    end if
    ifirst = 1
  endif

  ! seed BH formation from DEAD Pop III particles (KS ADDED)
  if(create_sinks .and. bhpop3)then
    do ilevel=levelmin,nlevelmax
      !call make_sink_from_pop3(ilevel)
      call make_sink_from_pop3_new_fg(ilevel)
    enddo
  endif

  ! !--------------- KS DEBUG ------------------!
  ! call sleep(1)
  ! call MPI_BARRIER(MPI_COMM_WORLD,info)
  ! if (myid==1) print *, "KS DEBUG A (create_sink_fg)"
  ! call sleep(1)
  ! call MPI_BARRIER(MPI_COMM_WORLD,info)
  ! do ilevel=1,levelmin
  !    call check_particle_consistency(ilevel)
  ! end do
  ! !--------------- KS DEBUG ------------------!



  ! Merge all particles to level 1
  do ilevel=levelmin-1,1,-1
    call merge_tree_fine(ilevel)
  end do

  ! !--------------- KS DEBUG ------------------!
  ! call sleep(1)
  ! call MPI_BARRIER(MPI_COMM_WORLD,info)
  ! if (myid==1) print *, "KS DEBUG B (create_sink_fg)"
  ! call sleep(1)
  ! call MPI_BARRIER(MPI_COMM_WORLD,info)

  ! do ilevel=1,levelmin
  !    call check_particle_consistency(ilevel)
  ! end do
  ! !--------------- KS DEBUG ------------------!


  ! Remove all particle clouds around old sinks (including the central one)
  call kill_entire_cloud_fg(1)


  ! !--------------- KS DEBUG ------------------!
  ! call sleep(1)
  ! call MPI_BARRIER(MPI_COMM_WORLD,info)
  ! if (myid==1) print *, "KS DEBUG C (create_sink_fg)"
  ! call sleep(1)
  ! call MPI_BARRIER(MPI_COMM_WORLD,info)

  ! do ilevel=1,levelmin
  !    call check_particle_consistency(ilevel)
  ! end do
  ! !--------------- KS DEBUG ------------------!

  ! DO NOT MODIFY FLAG2 BETWEEN CLUMP_FINDER AND MAKE_SINK_FROM_CLUMP
  if (create_sinks)then

!     ! Run the clump finder,(produce no output, keep clump arrays allocated)
!     call clump_finder(.false.,.true.)

!     ! Trim clumps down to R_accretion ball around peaks
!     if(clump_core)call trim_clumps

!     ! Compute simple additive quantities and means (1st moments)
!     call compute_clump_properties(uold(1,1))

!     ! Compute quantities relative to mean (2nd moments)
!     call compute_clump_properties_round2(uold(1,1))

!     ! Apply all checks and flag cells for sink formation
!     call flag_formation_sites

!     ! Create new sink particles if relevant
!     do ilevel=levelmin,nlevelmax
!        call make_sink_from_clump(ilevel)
!     end do

!     ! Deallocate clump finder arrays
!     deallocate(npeaks_per_cpu)
!     deallocate(ipeak_start)
!     if (ntest>0)then
!        deallocate(icellp)
!        deallocate(levp)
    !        deallocate(testp_sort)
!        deallocate(imaxp)
!     endif
!     call deallocate_all

    ! ! seed BH formation from DEAD Pop III particles (KS ADDED)
    ! if(bhpop3)then
    !    do ilevel=levelmin,nlevelmax
    !       !call make_sink_from_pop3(ilevel)
    !       call make_sink_from_pop3_new_fg(ilevel)
    !    enddo
    ! endif

  end if

  ! Merge sink for smbh runs
  if (smbh)then
    call merge_smbh_sink_fg
  else
    call merge_star_sink_fg
  end if

  ! Create new cloud particles
  call create_cloud_from_sink_fg

  !KS DEBUG
  call check_sink_vars_among_cpus_fg("create_sink_fg point A")

  ! Scatter particle to the grid
  ! Compute Bondi parameters
  do ilevel=1,nlevelmax
    call make_tree_fine(ilevel)
    call kill_tree_fine(ilevel)
    call virtual_tree_fine(ilevel)
    call collect_acczone_avg_fg(ilevel)
  end do

  !KS DEBUG
  call check_sink_vars_among_cpus_fg("create_sink_fg point B")

  ! Perform first accretion with seed mass
  ! Gather particles to levelmin
  do ilevel=nlevelmax,levelmin,-1
    call grow_sink_fg(ilevel,.true.)
    call merge_tree_fine(ilevel)
  end do

  !KS DEBUG
  call check_sink_vars_among_cpus_fg("create_sink_fg point C")
  call check_cloud_count_fg('create_sink_fg') !two-tier verification (TO ADDED 2026-11-05)
#ifdef TWO_TIER_CHECK
  call check_particle_tree('create_sink_fg') !TO ADDED
#endif


  ! Update hydro quantities for split cells
  if(hydro)then
    do ilevel=nlevelmax,levelmin,-1
      call upload_fine(ilevel)
      do ivar=1,nvar
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
      end do
      ! Update boundaries
      if(simple_boundary)call make_boundary_hydro(ilevel)
    end do
  end if

  ! Update the cloud particle properties at levelmin
  ! So they match their parent sink properties
  call update_cloud_fg(levelmin)

  ! if (myid==1) print *, "KS DEBUG (create_sink_fg) myid, acc_rate, acc_rate/dtnew, acc_rate/dtnew [Msun/yr]", &
  !      myid, acc_rate(1), acc_rate(1)/dtnew(levelmin), acc_rate(1)/dtnew(levelmin)*scale_m/scale_t/(1.989d33/(365.*24.*3600.))


  ! Effective accretion rate during last coarse step
  acc_rate(1:nsink)=acc_rate(1:nsink)/dtnew(levelmin)

  ! ! ir_eff and 5 are ratio of infalling energy which is radiated and protostellar radius
  ! if(ir_feedback)acc_lum(1:nsink)=ir_eff*acc_rate(1:nsink)*msink(1:nsink)/(5*6.955d10/scale_l)

  ! compute BH luminosity using acc_rate for coarse step (KS ADDED)
  call BH_luminosity_from_acc_rate_fg

  ! Compute and print accretion rates
  call compute_accretion_rate_fg(.true.)

  !KS DEBUG
  call check_sink_vars_among_cpus_fg("create_sink_fg point D")



end subroutine create_sink_fg
!################################################################
subroutine create_cloud_from_sink_fg
  use amr_commons
  use pm_commons
  use hydro_commons
  use mpi_mod
  implicit none

  !----------------------------------------------------------------------
  ! This routine creates the whole cloud of particles for each sink,
  ! Particles are produced in the right MPI domain and inserted in the
  ! linked list at level 1.
  ! The cloud radius is dble(ir_cloud)*dx_min, where dx_min is
  ! the cell size at levelmax. For cosmo runs, the cloud radius is
  ! dx_min/aexp (therefore it is constant in *physical* units).
  !----------------------------------------------------------------------

  real(dp)::scale,dx_min,rr,rmax,rmass
  integer ::i,icpu,isink,indp,ii,jj,kk,nx_loc,idim
  integer ::ntot,ntot_all,info
  logical ::ok_free
  real(dp),dimension(1:ndim)::xrel
  real(dp),dimension(1:nvector,1:ndim)::xtest
  integer ,dimension(1:nvector)::ind_grid,ind_part,cc,ind_cloud
  logical ,dimension(1:nvector)::ok_true
  logical,dimension(1:ndim)::period
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ok_true=.true.

  if(numbtot(1,1)==0) return
  if(verbose)write(*,*)' Entering create_cloud_from_sink_fg'

#if NDIM==3

  ! Level 1 linked list
  do icpu=1,ncpu
    if(numbl(icpu,1)>0)then
      ind_grid(1)=headl(icpu,1)
    endif
  end do

  period(1)=(nx==1)
  if(ndim>1)period(2)=(ny==1)
  if(ndim>2)period(3)=(nz==1)

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  !dx_min=scale*0.5D0**nlevelmax/aexp
  dx_min=scale*0.5D0**nlevelmax !*** KS MODIFEID ***

  rmax=dble(ir_cloud)*dx_min
  rmass=dble(ir_cloud_massive)*dx_min

  do kk=-2*ir_cloud,2*ir_cloud
    xrel(3)=dble(kk)*dx_min/2.0
    do jj=-2*ir_cloud,2*ir_cloud
      xrel(2)=dble(jj)*dx_min/2.0
      do ii=-2*ir_cloud,2*ir_cloud
        xrel(1)=dble(ii)*dx_min/2.0
        rr=sqrt(sum(xrel**2))
        if(rr<=rmax)then
          do isink=1,nsink
            xtest(1,1:3)=xsink(isink,1:3)+xrel(1:3)
            do idim=1,ndim
              if (period(idim) .and. xtest(1,idim)>boxlen)xtest(1,idim)=xtest(1,idim)-boxlen
              if (period(idim) .and. xtest(1,idim)<0.)xtest(1,idim)=xtest(1,idim)+boxlen
            end do
            call cmp_cpumap(xtest,cc,1) !TO MEMO: two-tier: coarse-tier owner; any unique cpu will do, the particle tree moves the cloud particle to its grid
            if(cc(1).eq.myid)then
              call remove_free(ind_cloud,1)
              call add_list(ind_cloud,ind_grid,ok_true,1)

              ! !KS DEBUG
              ! print *, "KS DEBUG (create_cloud_from_sink_fg) kk, jj, ii, ind_cloud(1) = ", kk, jj, ii, ind_cloud(1)

              indp=ind_cloud(1)
              idp(indp)=-isink
              typep(indp)%family=FAM_CLOUD_FG
              typep(indp)%tag=0
              levelp(indp)=levelmin
              !rewrite to avoid floating overflow (KS MODIFIED)
              !if (rr<=rmass .and. msink(isink)<mass_sink_direct_force*1.9891d33/(scale_d*scale_l**3))then
              if (rr<=rmass .and. msink(isink)/(1.9891d33/(scale_d*scale_l**3)) < mass_sink_direct_force)then
                mp(indp)=msink(isink)/dble(ncloud_sink_massive)
              else
                mp(indp)=0.
              end if
              xp(indp,1:3)=xtest(1,1:3)
              vp(indp,1:3)=vsink(isink,1:3)
              tp(indp)=tsink(isink)     ! Birth epoch

              !print *, "(KS DEBUG in create_cloud_from_sink_fg) idp, ii, jj, kk, xp(:), mp = ", idp(indp), ii, jj, kk, xp(indp,:), mp(indp)
            end if
          end do
        end if
      end do
    end do
  end do

  do isink=1,nsink
    !rewrite to avoid floating overflow (KS MODIFIED)
    !direct_force_sink(isink)=(msink(isink) .ge. mass_sink_direct_force*1.9891d33/(scale_d*scale_l**3))
    direct_force_sink(isink)=(msink(isink)/(1.9891d33/(scale_d*scale_l**3)) .ge. mass_sink_direct_force)
  end do

#endif
end subroutine create_cloud_from_sink_fg
!################################################################
subroutine kill_entire_cloud_fg(ilevel)
  use pm_commons
  use amr_commons
  use mpi_mod
  implicit none
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine removes cloud particles (including the central one).
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::ig,ip,npart1,npart2,icpu
  integer,dimension(1:nvector)::ind_grid,ind_part,ind_grid_part
  logical,dimension(1:nvector)::ok=.true.

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  ! Gather sink and cloud particles.
  ! Loop over cpus
  do icpu=1,ncpu
    igrid=headl(icpu,ilevel)
    ig=0
    ip=0
    ! Loop over grids
    do jgrid=1,numbl(icpu,ilevel)
      npart1=numbp(igrid)  ! Number of particles in the grid
      npart2=0
      ! Count sink and cloud particles
      if(npart1>0)then
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
          ! Save next particle   <--- Very important !!!
          next_part=nextp(ipart)
          if(idp(ipart).lt.0)then
            npart2=npart2+1
          endif
          ipart=next_part  ! Go to next particle
        end do
      endif
      !print '(A,4I6)', "KS DEBUG (kill_entire_cloud_fg)  myid, ilevel, npart1, npart2: ", myid, ilevel, npart1, npart2

      ! Gather sink and cloud particles
      if(npart2>0)then
        ig=ig+1
        ind_grid(ig)=igrid
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
          ! Save next particle   <--- Very important !!!
          next_part=nextp(ipart)
          ! Select only sink particles
          if(idp(ipart).lt.0)then
            if(ig==0)then
              ig=1
              ind_grid(ig)=igrid
            end if
            ip=ip+1
            ind_part(ip)=ipart
            ind_grid_part(ip)=ig
          endif
          if(ip==nvector)then
            call remove_list(ind_part,ind_grid_part,ok,ip)
            call add_free_cond(ind_part,ok,ip)
            ip=0
            ig=0
          end if
          ipart=next_part  ! Go to next particle
        end do
        ! End loop over particles
      end if

      igrid=next(igrid)   ! Go to next grid
    end do

    ! End loop over grids
    if(ip>0)then
      call remove_list(ind_part,ind_grid_part,ok,ip)
      call add_free_cond(ind_part,ok,ip)
    end if
  end do
111 format('   Entering kill_cloud for level ',I2)
end subroutine kill_entire_cloud_fg
!################################################################
subroutine collect_acczone_avg_fg(ilevel)
  use pm_commons
  use amr_commons
  use poisson_commons
  !use mpi_mod
  use mpi_mod
  implicit none

#ifndef WITHOUTMPI
  integer::info
#endif
  !-----------------------------------------------------------------------------
  ! This routine is used to collect all relevant information to compute the
  ! the accretion rates. The information is collected level-by-level when
  ! going down in the call tree (leafs at the bottom), while accretion is
  ! performed on the way up.
  ! - first, the volume of each particle is computed. The volume is reduced if
  ! sinks are overlapping.
  ! - then a loop over all particles of ilevel (vectorized) is used to compute
  ! the volume-weighted quantities.
  !-----------------------------------------------------------------------------

  integer::ilevel
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::ig,ip,npart1,npart2,icpu,isink
  integer,dimension(1:nvector)::ind_grid,ind_part,ind_grid_part

  if(ilevel<levelmin)return
  if(verbose)write(*,111)ilevel

  ! print *, "KS DEBUG (collect_acczone_avg_fg) ", ilevel

  level_sink_new(1:nsinkmax,ilevel)=.false. !grid levels that have an overlap with sink radius

  ! Compute (volume weighted) averages over accretion zone
  wden=0d0; wvol=0d0; weth=0d0; wmom=0d0

  ! Loop over cpus
  do icpu=1,ncpu
    igrid=headl(icpu,ilevel)
    ig=0
    ip=0
    ! Loop over grids
    do jgrid=1,numbl(icpu,ilevel)
      npart1=numbp(igrid)  ! Number of particles in the grid
      npart2=0

      ! Count sink and cloud particles
      if(npart1>0)then
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
          ! Save next particle   <--- Very important !!!
          next_part=nextp(ipart)
          !KS MODIFIED
          !if( is_cloud(typep(ipart)) ) then
          if(idp(ipart).lt.0)then
            npart2=npart2+1
          endif
          ipart=next_part  ! Go to next particle
        end do
      endif

      ! Gather sink and cloud particles
      if(npart2>0)then
        ig=ig+1
        ind_grid(ig)=igrid
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
          ! Save next particle   <--- Very important !!!
          next_part=nextp(ipart)
          ! Select only sink particles
          !KS MODIFIED
          !if( is_cloud(typep(ipart)) ) then
          if(idp(ipart).lt.0)then
            if(ig==0)then
              ig=1
              ind_grid(ig)=igrid
            end if
            ip=ip+1
            ind_part(ip)=ipart
            ind_grid_part(ip)=ig
          endif
          if(ip==nvector)then
            call collect_acczone_avg_np_fg(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
            ip=0
            ig=0
          end if
          ipart=next_part  ! Go to next particle
        end do
        ! End loop over particles
      end if
      igrid=next(igrid)   ! Go to next grid
    end do

    ! End loop over grids
    if(ip>0)then
      call collect_acczone_avg_np_fg(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
    end if
  end do
  ! End loop over cpus

  if(nsink>0)then
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(wden,wden_new,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(wvol,wvol_new,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(weth,weth_new,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(wmom,wmom_new,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(level_sink_new,level_sink,nsinkmax*(nlevelmax-levelmin+1),MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,info)
#else
    wden_new=wden
    wvol_new=wvol
    weth_new=weth
    wmom_new=wmom
    level_sink=level_sink_new
#endif
  endif

  do isink=1,nsink
    weighted_density(isink,ilevel)=wden_new(isink)
    weighted_volume(isink,ilevel)=wvol_new(isink)
    weighted_momentum(isink,ilevel,1:ndim)=wmom_new(isink,1:ndim)
    weighted_ethermal(isink,ilevel)=weth_new(isink)
  end do

111 format('   Entering collect_acczone_avg_fg for level ',I2)

end subroutine collect_acczone_avg_fg
!################################################################
subroutine collect_acczone_avg_np_fg(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use mpi_mod
  implicit none
  !-----------------------------------------------------------------------------
  ! inner loop of collect_acczone_avg_fg
  ! weighted gas quantities for each particle are computed
  ! no CIC averaging over quantities anymore as average over whole sink
  ! accretion zone is computed
  !-----------------------------------------------------------------------------

  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid,ind_part,ind_grid_part
  integer::j,nx_loc,isink,idim,ind
#if NENER>0
  integer::irad
#endif
  real(dp)::d,e,v2
  real(dp)::scale,weight,dx_cloud,vol_cloud
  real(dp),dimension(1:ndim)::vv
#ifdef SOLVERmhd
  real(dp)::bx1,bx2,by1,by2,bz1,bz2
#endif
  real(dp),dimension(1:nvector,1:ndim),save::xpart
  integer ,dimension(1:nvector,1:twotondim),save::indp
  real(dp),dimension(1:nvector,1:ndim,1:twotondim)::xx
  real(dp),dimension(1:nvector,1:twotondim)::vol
  logical, dimension(1:nvector,1:twotondim)::ok

  ! Compute volume of each cloud particle
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  !dx_cloud=(0.5D0**nlevelmax)*scale/aexp/2 ! factor of 2 hard-coded
  dx_cloud=(0.5D0**nlevelmax)*scale/2 ! factor of 2 hard-coded *** KS MODIFIED ***
  vol_cloud=dx_cloud**ndim

  ! Copy cloud particle coordinates
  do idim=1,ndim
    do j=1,np
      xpart(j,idim)=xp(ind_part(j),idim)
    end do
  end do

  ! Compute cloud particle CIC weights at the current level
  call cic_get_cells_fg(indp,xx,vol,ok,ind_grid,xpart,ind_grid_part,ng,np,ilevel)

  do ind=1,twotondim
    do j=1,np
      if(ok(j,ind))then

        ! Convert uold to primitive variables
        d=max(uold(indp(j,ind),1),smallr)
        vv(1)=uold(indp(j,ind),2)/d
        vv(2)=uold(indp(j,ind),3)/d
        vv(3)=uold(indp(j,ind),4)/d

        ! Compute fluid total energy
        e=uold(indp(j,ind),5)
#ifdef SOLVERmhd
        bx1=uold(indp(j,ind),6)
        by1=uold(indp(j,ind),7)
        bz1=uold(indp(j,ind),8)
        bx2=uold(indp(j,ind),nvar+1)
        by2=uold(indp(j,ind),nvar+2)
        bz2=uold(indp(j,ind),nvar+3)
        e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)
#endif
#if NENER>0
        do irad=0,nener-1
          e=e-uold(indp(j,ind),inener+irad)
        end do
#endif
        e=e/d ! Specific energy
        v2=sum(vv**2)
        e=e-0.5d0*v2 ! Remove kinetic energy

        ! Get sink index
        isink=-idp(ind_part(j))

        ! Get cloud particle CIC weight
        weight=vol_cloud*vol(j,ind)

        ! Compute sink average quantities
        wvol(isink)=wvol(isink)+weight
        wden(isink)=wden(isink)+weight*d
        wmom(isink,1:ndim)=wmom(isink,1:ndim)+weight*d*vv(1:ndim)
        weth(isink)=weth(isink)+weight*d*e

        if (weight>0.)level_sink_new(isink,ilevel)=.true.
        ! !------------------- KS DEBUG (BEGIN) ----------------!
        ! if (isink == 1) then
        !    print '(A,5I10,/,A,1P8E14.6)', "KS DEBUG (collect_acczone_avg_np_fg) isink, myid, ind_part, ind, indp = ", &
        !         isink, myid, ind_part(j), ind, indp(j,ind),&
        !         "xpart, xx, weight, d = ", &
        !         xpart(j,1:3), xx(j,1:3,ind), weight, d
        ! end if
        ! !------------------- KS DEBUG (END) ----------------!



      endif
    end do
  end do

end subroutine collect_acczone_avg_np_fg
!################################################################
subroutine grow_sink_fg(ilevel,on_creation)
  use pm_commons
  use amr_commons
  use hydro_commons
  use mpi_mod
  implicit none
  integer::ilevel
  logical::on_creation
  !------------------------------------------------------------------------
  ! This routine performs accretion onto the sink. It vectorizes the loop
  ! over all sink cloud particles and calls accrete_sink as soon as nvector
  ! particles are collected
  ! -> replaces grow_bondi and grow_jeans
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,info
  integer::ig,ip,npart1,npart2,icpu,isink,lev,nx_loc
  integer,dimension(1:nvector)::ind_grid,ind_part,ind_grid_part
  real(dp)::scale,dx_min,vol_min
  real(dp),dimension(1:ndim)::old_loc,old_vel
  logical::highest_level

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m !KS DEBUG
  character(LEN=256):: debug_message,nchar !KS DEBUG

  !KS DEBUG
  integer::ilevel_tmp
  real(dp),dimension(1:3)::velocity
  real(dp)::v2,c2,density,volume,ethermal,mgas,mgas_ilevel

  ! print *, "KS DEBUG (grow_sink_fg)", ilevel

!------------------- KS DEBUG (BEGIN) ----------------!
!if (myid == 1) print *, "**** KS DEBUG **** grow_sink_fg is called at level = ", ilevel, " with on_creation = ", on_creation
!if (myid == 1) msink(:)=msink(:)+1d-1
  write(nchar,'(I8)') ilevel
  debug_message = 'BEGINNING OF grow_sink_fg at ilevel = ' // TRIM(ADJUSTL(nchar))
  call check_sink_vars_among_cpus_fg(debug_message)
!------------------- KS DEBUG (END) ----------------!


  ! if (myid == 1) print *, "**** KS DEBUG **** grow_sink_fg is called at level = ", ilevel, " with on_creation = ", on_creation

  ! Conversion factor from user units to cgs units (KS DEBUG)
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3

  if(accretion_scheme=='none')return

  !--------------------------------------------------------------------
  ! following highest level check is deleted in the latest version of Ramses (KS MEMO)
  !--------------------------------------------------------------------
  ! if(.not. on_creation)then
  !    ! Determine if on highest active level...
  !    if (ilevel==nlevelmax)then
  !       highest_level=.true.
  !    else if (numbtot(1,ilevel+1)==0)then
  !       highest_level=.true.
  !    else
  !       highest_level=.false.
  !    end if
  !    if (.not. highest_level)return
  ! end if

  if(verbose)write(*,111)ilevel

  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  !dx_min=scale*(0.5D0**nlevelmax)/aexp
  dx_min=scale*0.5D0**nlevelmax !*** KS MODIFEID ***
  vol_min=dx_min**ndim

#if NDIM==3

  ! Compute sink accretion rates
  call compute_accretion_rate_fg(.false.)

  ! Reset new sink variables
  msink_new=0d0; xsink_new=0.d0; vsink_new=0d0; delta_mass_new=0d0; lsink_new=0d0

  ! Loop over cpus
  do icpu=1,ncpu
    igrid=headl(icpu,ilevel)
    ig=0
    ip=0
    ! Loop over grids
    do jgrid=1,numbl(icpu,ilevel)
      npart1=numbp(igrid)  ! Number of particles in the grid
      npart2=0
      ! Count sink and cloud particles
      if(npart1>0)then
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
          ! Save next particle   <--- Very important !!!
          next_part=nextp(ipart)
          if(idp(ipart).lt.0)then
            npart2=npart2+1
          endif
          ipart=next_part  ! Go to next particle
        end do
      endif
      ! Gather sink and cloud particles
      if(npart2>0)then
        ig=ig+1
        ind_grid(ig)=igrid
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
          ! Save next particle   <--- Very important !!!
          next_part=nextp(ipart)
          ! Select only sink particles
          if(idp(ipart).lt.0)then
            if(ig==0)then
              ig=1
              ind_grid(ig)=igrid
            end if
            ip=ip+1
            ind_part(ip)=ipart
            ind_grid_part(ip)=ig
          endif
          if(ip==nvector)then
            !call accrete_sink(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,on_creation)
            call accrete_sink_new_fg(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,on_creation) !KS MODIFIED
            ip=0
            ig=0
          end if
          ipart=next_part  ! Go to next particle
        end do
        ! End loop over particles
      end if
      igrid=next(igrid)   ! Go to next grid
    end do
    ! End loop over grids
    !if(ip>0)call accrete_sink(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,on_creation)
    if(ip>0)call accrete_sink_new_fg(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,on_creation) !KS MODIFIED
  end do
  ! End loop over cpus
  if(nsink>0)then
#ifndef WITHOUTMPI
    !----------------------------- (KS DEBUG) -------------------------------!
    call MPI_BARRIER(MPI_COMM_WORLD,info)  !KS DEBUG
    !----------------------------- (KS DEBUG) -------------------------------!

    call MPI_ALLREDUCE(msink_new,msink_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(delta_mass_new,delta_mass_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(xsink_new,xsink_all,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(vsink_new,vsink_all,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(lsink_new,lsink_all,nsinkmax*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
    msink_all=msink_new
    delta_mass_all=delta_mass_new
    xsink_all=xsink_new
    vsink_all=vsink_new
    lsink_all=lsink_new
#endif
  endif
  do isink=1,nsink
    ! Reset jump in sink coordinates
    do lev=levelmin,nlevelmax
      sink_jump(isink,1:ndim,lev)=sink_jump(isink,1:ndim,lev)-xsink(isink,1:ndim)
    end do

    !save old velocity and location to compute deltas
    old_loc(1:ndim)=xsink(isink,1:ndim)
    old_vel(1:ndim)=vsink(isink,1:ndim)

    ! Change to conservative quantities
    xsink(isink,1:ndim)=xsink(isink,1:ndim)*msink(isink)
    vsink(isink,1:ndim)=vsink(isink,1:ndim)*msink(isink)


    ! !---------------- KS DEBUG ---------------!
    ! if (myid == 1) then
    !    if (SUM(vsink_all(isink,1:ndim)**2) > 0) then
    !       print '(A,2I4,2(/,1P3E24.16,1PE10.3))', "KS DEBUG (grow_sink_fg, finite v update) isink, ilevel, vsink(before), vsink(after) =", isink, ilevel,&
    !            vsink(isink,1:ndim)/msink(isink), sqrt(SUM(vsink(isink,1:ndim)/msink(isink))**2), &
    !            (vsink(isink,1:ndim)+vsink_all(isink,1:ndim))/msink(isink), sqrt(SUM((vsink(isink,1:ndim)+vsink_all(isink,1:ndim))/msink(isink))**2)
    !    end if
    ! end if
    ! !------------------------------------------!




    ! Accrete to sink variables
    msink(isink)=msink(isink)+msink_all(isink)
    xsink(isink,1:ndim)=xsink(isink,1:ndim)+xsink_all(isink,1:ndim)
    vsink(isink,1:ndim)=vsink(isink,1:ndim)+vsink_all(isink,1:ndim)
    !compute lsink with reference point of old xsink
    lsink(isink,1:3)=lsink(isink,1:3)+lsink_all(isink,1:3)

    ! Change back
    xsink(isink,1:ndim)=xsink(isink,1:ndim)/msink(isink)
    vsink(isink,1:ndim)=vsink(isink,1:ndim)/msink(isink)

    ! !KS DEBUG
    ! if(myid==56) then
    !    print '(/,A,/,A,/, 3I10,1P6E14.6)', "*** grow_sink_fg (myid=56)***",&
    !         "nstep_coarse, nstep, isink, xsink(isink,:), vsink(isink,:): ",&
    !         nstep_coarse, nstep, isink, xsink(isink,:), vsink(isink,:)
    ! end if
    ! if(myid==1) then
    !    print '(/,A,/,A,/, 3I10,1P6E14.6)', "*** grow_sink_fg (myid=1)***",&
    !         "nstep_coarse, nstep, isink, xsink(isink,:), vsink(isink,:): ",&
    !         nstep_coarse, nstep, isink, xsink(isink,:), vsink(isink,:)
    ! end if




    !correct for new center of mass location/velocity
    lsink(isink,1:3)=lsink(isink,1:3)-msink(isink)*cross((xsink(isink,1:ndim)-old_loc(1:ndim)),(vsink(isink,1:ndim)-old_vel(1:ndim)))

    ! Store jump in sink coordinates
    do lev=levelmin,nlevelmax
      sink_jump(isink,1:ndim,lev)=sink_jump(isink,1:ndim,lev)+xsink(isink,1:ndim)
    end do

    ! !KS DEBUG
    ! print *, "KS DEBUG (grow_sink_fg) myid, ilevel, acc_rate, msink_new, msink_all, dtnew(ilevel)", myid, ilevel, acc_rate(isink), msink_new(isink), msink_all(isink), dtnew(ilevel)
    ! print *, "KS DEBUG (grow_sink_fg 2) myid, dMsink_overdt, dM/dt_loc, dM/dt_loc/dMsink_overdt", &
    !      myid,&
    !      dMsink_overdt(isink)*scale_m/scale_t/(1.989d33/(365.*24.*3600.)), &
    !      msink_new(isink)/dtnew(ilevel)*scale_m/scale_t/(1.989d33/(365.*24.*3600.)), &
    !       (msink_new(isink)/dtnew(ilevel))/dMsink_overdt(isink)

    ! print *, "KS DEBUG (grow_sink_fg 3) myid,  mcell/mgas, mgas(code)",  myid,&
    !      msink_new(isink)/(dMsink_overdt(isink)*dtnew(ilevel)),rho_gas(isink)*volume_gas(isink)



    !--------------------------------------------- KS DEBUG ---------------------------------------------------!
    ! if (myid == 1) then
    !    if (msink_all(isink)>0) then
    !       !copy and paste from compute_accretion_rate_fg
    !       density=0.d0; volume=0.d0; velocity=0.d0; ethermal=0d0
    !       do ilevel_tmp=levelmin,nlevelmax
    !          ! print *, "density: ", density
    !          ! print *, "isink,ilevel_tmp, weighted_density(isink,ilevel_tmp): ", isink,ilevel_tmp,weighted_density(isink,ilevel_tmp)
    !          density=density+weighted_density(isink,ilevel_tmp)
    !          ethermal=ethermal+weighted_ethermal(isink,ilevel_tmp)
    !          velocity(1:3)=velocity(1:3)+weighted_momentum(isink,ilevel_tmp,1:3)
    !          volume=volume+weighted_volume(isink,ilevel_tmp)
    !       end do
    !       mgas=density
    !       mgas_ilevel = weighted_density(isink,ilevel)
    !       density=density/(volume+tiny(0.0_dp))
    !       velocity(1:3)=velocity(1:3)/(density*volume+tiny(0.0_dp))
    !       ethermal=ethermal/(density*volume+tiny(0.0_dp))
    !       c2=MAX((gamma-1.0)*ethermal,smallc**2)
    !       v2=SUM((velocity(1:3)-vsink(isink,1:3))**2)


    !       print '(A, 2I5, 1P3E15.6,1PE20.11,1PE12.2,/,8(A,1PE9.2))', &
    !            "KS DEBUG (grow_sink_fg, delta_mass>0) ilevel, isink, acc_rate(old), acc_rate(new), dtnew, t, delta_mass (actual): ", &
    !            ilevel, isink, acc_rate(isink), acc_rate(isink)+msink_all(isink), &
    !            dtnew(ilevel),t, msink_all(isink), &
    !            'msink[Msun]=',msink(isink)*scale_m/1.989d33, &
    !            ', density[cm^3]=', density*scale_nH, &
    !            ', csinf[km/s]=', sqrt(c2)*scale_v/1d5, &
    !            ', vinf[km/s]=', sqrt(v2)*scale_v/1d5, &
    !            ', mdotsink[Msun/yr]=', dMsink_overdt(isink)*scale_m/scale_t/(1.989d33/(365.*24.*3600.)),&
    !            ', mgas[Msun]=', mgas*scale_m/1.989d33,&
    !            ', mgas_ilevel[Msun]=', mgas_ilevel*scale_m/1.989d33,&
    !            ', delta_mass (recalculated) =',mgas_ilevel/mgas * dMsink_overdt(isink)*dtnew(ilevel)
    !    end if
    ! end if

    ! Store accreted mass
    acc_rate(isink)=acc_rate(isink)+msink_all(isink)

    if (agn.and.ok_blast_agn(isink))then
      delta_mass(isink)=delta_mass_all(isink) ! only the most recent accretion
    else
      delta_mass(isink)=delta_mass(isink)+delta_mass_all(isink)
    end if

    ! if (myid==1) print *, "KS DEBUG (grow_sink_fg 4) myid, acc_rate, acc_rate/dtnew [Msun/yr]",  myid, acc_rate(isink), &
    !      acc_rate(isink)/dtnew(ilevel)*scale_m/scale_t/(1.989d33/(365.*24.*3600.))
  end do

!------------------- KS DEBUG (BEGIN) ----------------!
  debug_message = 'END OF grow_sink_fg at ilevel = ' // TRIM(ADJUSTL(nchar))
  call check_sink_vars_among_cpus_fg(debug_message)
!------------------- KS DEBUG (END) ----------------!


#endif
111 format('   Entering grow_sink_fg for level ',I2)

end subroutine grow_sink_fg
!################################################################
subroutine accrete_sink_new_fg(ind_grid,ind_part,ind_grid_part,ng,np,ilevel,on_creation)
  use amr_commons
  use pm_commons
  use hydro_commons
  use mpi_mod
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  logical::on_creation

  !-----------------------------------------------------------------------
  ! This routine is called by subroutine grow_sink_fg. It performs accretion
  ! for nvector particles. Routine is not very efficient. Optimize if taking too long...
  !-----------------------------------------------------------------------

  integer::i,j,nx_loc,isink,ivar,irad,idim,ind
  real(dp)::v2,d,e,d_floor,density,volume
#ifdef SOLVERmhd
  real(dp)::bx1,bx2,by1,by2,bz1,bz2
#endif
  real(dp),dimension(1:nvar)::z
  real(dp)::factG,scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx,dx_loc,dx_min,scale,vol_min,vol_loc,weight,acc_mass,temp,d_jeans
  real(dp)::dx_cloud,vol_cloud !KS ADDED
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim)::xpart
  real(dp),dimension(1:nvector,1:ndim,1:twotondim)::xx
  real(dp),dimension(1:nvector,1:twotondim)::vol
  integer ,dimension(1:nvector)::cell_lev
  ! Particle based arrays
  logical,dimension(1:nvector,1:twotondim)::ok
  integer ,dimension(1:nvector,1:twotondim)::indp
  real(dp),dimension(1:3)::vv

  real(dp),dimension(1:3)::r_rel,x_acc,p_acc,p_rel,p_rel_rad,p_rel_acc,p_rel_tan,delta_x,delta_p,drag
  real(dp)::r_abs,fbk_ener_AGN,fbk_mom_AGN
  logical,dimension(1:ndim)::period
  real(dp)::virt_acc_mass,delta_e_tot,Mred,Macc
  real(dp),dimension(1:nsinkmax)::delta_M

  real(dp)::tan_theta,cone_dist,orth_dist
  real(dp),dimension(1:3)::cone_dir

  real(dp),parameter ::c_light  = 2.9979250d+10 ! [cm/s] !KS ADDED
  real(dp) :: delta_e, delta_e_max, T2gas  !KS ADDED

  !print *, "KS DEBUG (accrete_sink_new_fg)", ilevel

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  period(1)=(nx==1)
#if NDIM>1
  if(ndim>1)period(2)=(ny==1)
#endif
#if NDIM>2
  if(ndim>2)period(3)=(nz==1)
#endif

#if NDIM==3

  ! Gravitational constant
  factG=1d0
  if(cosmo)factG=3d0/8d0/3.1415926*omega_m*aexp

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim
  !dx_min=scale*0.5D0**nlevelmax/aexp
  dx_min=scale*0.5D0**nlevelmax !*** KS MODIFEID ***
  vol_min=dx_min**ndim

  ! Compute volume of each cloud particle
  dx_cloud=dx_min/2 ! factor of 2 hard-coded
  vol_cloud=dx_cloud**ndim

  do idim=1,ndim
    do j=1,np
      xpart(j,idim)=xp(ind_part(j),idim)
    end do
  end do

  call cic_get_cells_fg(indp,xx,vol,ok,ind_grid,xpart,ind_grid_part,ng,np,ilevel)

  do ind=1,twotondim
    do j=1,np
      if(ok(j,ind))then

        ! Convert uold to primitive variables
        d=max(uold(indp(j,ind),1),smallr)
        vv(1)=uold(indp(j,ind),2)/d
        vv(2)=uold(indp(j,ind),3)/d
        vv(3)=uold(indp(j,ind),4)/d
        e=uold(indp(j,ind),5)

        v2=(vv(1)**2+vv(2)**2+vv(3)**2)
        e=e-0.5d0*d*v2
        e=e/d               !thermal energy per mass (KS MEMO)
        do ivar=imetal,nvar
          z(ivar)=uold(indp(j,ind),ivar)/d  !metallicity per mass (KS MEMO)
        end do

        ! Get sink index
        isink=-idp(ind_part(j))

        ! Cloud particle CIC weight (KS MODIFIED)
        weight=vol_cloud*vol(j,ind)

        ! Compute sink average density
        !weight=weightp(ind_part(j),ind)
        density=rho_gas(isink)
        volume=volume_gas(isink)
        if (volume==0. .or. density==0.)print*,'(accrete_sink_new_fg) something might be going wrong here...',weight,volume,density,ilevel

        if (on_creation)then
          if (new_born(isink))then
            !KS MODIFIED
            print *, "** ERROR ** current code can be run only with new_born(isink) = .false."
            print *, "stopping..."
            call clean_stop
          else
            acc_mass = 0. ! on sink creation, sink created from particle (KS MODIFIED)
          end if
        else
          if (bondi_accretion)then
            acc_mass=dMsink_overdt(isink)*dtnew(ilevel)*weight/volume*d/density

            ! print *, "KS DEBUG (accrete_sink_new_fg) dMsink_overdt, dt_new, weight/volume*d/density",  &
            !      dMsink_overdt(isink)*(scale_d*scale_l**3d0/scale_t)/(1.9891d33/(365.*24.*3600.)), &
            !      dtnew(ilevel)*scale_t/(365.*24.*3600.), weight/volume*d/density
          else
            !KS MODIFIED
            print *, "** ERROR ** current code can be run only with bondi_accretion"
            print *, "stopping..."
            call clean_stop
          endif

          ! No neg accretion
          acc_mass=max(acc_mass,0.0_dp)
        end if

        ! reference frame relative to the sink
        r_rel(1:3)=xx(j,1:3,ind)-xsink(isink,1:3)
        do idim=1,ndim
          if (period(idim) .and. r_rel(idim)>boxlen*0.5)xx(j,idim,ind)=xx(j,idim,ind)-boxlen
          if (period(idim) .and. r_rel(idim)<boxlen*(-0.5))xx(j,idim,ind)=xx(j,idim,ind)+boxlen
        end do
        r_rel(1:3)=xx(j,1:3,ind)-xsink(isink,1:3)
        r_abs=sum(r_rel(1:3)**2)**0.5

        ! momentum in relative motion
        p_rel=d*vol_loc*(vv(1:3)-vsink(isink,1:3))

        p_rel_acc=p_rel*acc_mass/(d*vol_loc)

        ! total accreted momentum
        p_acc(1:3)=p_rel_acc+vsink(isink,1:3)*acc_mass

        ! total accreted center of mass
        x_acc(1:3)=acc_mass*xx(j,1:3,ind)

        ! Add accreted properties to sink variables
        msink_new(isink)=msink_new(isink)+acc_mass
        delta_mass_new(isink)=delta_mass_new(isink)+acc_mass
        xsink_new(isink,1:3)=xsink_new(isink,1:3)+x_acc(1:3)
        vsink_new(isink,1:3)=vsink_new(isink,1:3)+p_acc(1:3)
        lsink_new(isink,1:3)=lsink_new(isink,1:3)+cross(r_rel(1:3),p_rel_acc(1:3))

        ! Check for neg density inside sink accretion radius
        if (8*acc_mass/vol_loc>d .and. (.not. on_creation))then
          write(*,*)'====================================================='
          write(*,*)'DANGER of neg density :-( at location'
          write(*,*)xx(j,1:3,ind)
          write(*,*)'due to',isink
          write(*,*)acc_mass/vol_loc,d/d_sink,density/d_sink
          write(*,*)indp(j,ind),myid
          do i=1,nsink
            write(*,*)i,' distance ',sum((xx(j,1:3,ind)-xsink(i,1:3))**2)**0.5/dx_min
          end do
          write(*,*)'try to decrease c_acc in SINK_PARAMS'
          write(*,*)'====================================================='
          !           call clean_stop
        end if

        !KS DEBUG
        ! print *, "KS DEBUG (accrete_sink_new_fg) j, ind_part(j), indp(j,ind), xx(j,1:3,ind) = ", ind, j, ind_part(j), indp(j,ind), xx(j,1:3,ind)


        ! regular accretion from unew
        unew(indp(j,ind),1)=unew(indp(j,ind),1)-acc_mass/vol_loc
        unew(indp(j,ind),2:4)=unew(indp(j,ind),2:4)-vv(1:3)*acc_mass/vol_loc
        ! fix the energy equation to account for the work of the truly accreted material
        ! subtract internal (kinetic+thermal) energy in proportinal to subtracted mass (KS MEMO)
        unew(indp(j,ind),5)=unew(indp(j,ind),5)-uold(indp(j,ind),5)*acc_mass/(d*vol_loc)


        !------------------------------------------------
        !thermal feedback (KS MODIFIED)
        !------------------------------------------------
        if (BH_model == 3) then
          !10 % efficiency for luminosity and 5 % efficiency for thermal energy deposit (cf. Di Matteo 2012)
          delta_e = 0.1 * 0.05 * acc_mass/vol_loc * (c_light/scale_v)**2

          !maximum energy deposit based on sgagn_T2max (require updated T2 < sgagn_T2max)
          T2gas = e*(gamma-1.0)*scale_T2 !e obtained above is energy/mass
          delta_e_max = max((d-acc_mass/vol_loc), 0d0) * max((sgagn_T2max-T2gas)/scale_T2, 0d0) / (gamma-1.0)

          !add energy
          unew(indp(j,ind),5) = unew(indp(j,ind),5) + min(delta_e, delta_e_max)

          ! !KS DEBUG
          ! if (isink == 1) then !sink 1に限定
          !    if (d*scale_nH > 1d2) & !高密度領域に限定
          !    print '(A, A, 2I10, A, 2I10,  A, 1P3E20.12,/ 9(A, 1P1E10.2),/, A, 1P9E10.2, A, 1P9E10.2,/)', &
          !         "KS DEBUG (accrete_sink_new_fg)", &
          !         " j, ind: ", j,ind,&
          !         ", myid, ind_cell: ", myid, indp(j,ind), &
          !         ", xx: ", xx(j,1:3,ind),&
          !         "rho[cm-3]=", d*scale_nH, &
          !         ", T2gas[K]=", T2gas,&
          !         ", weight[pc^3]=", weight*(scale_l/3.08568d18)**3,&
          !         ", mass[Msun]=", weight*d*scale_l**3*scale_d/1.989d33, &
          !         ", mass-frac=", weight/volume*d/density, &
          !         ", drho[cm-3]=", acc_mass/vol_loc*scale_nH, &
          !         ", dM[Msun]=", acc_mass*scale_l**3*scale_d/1.989d33, &
          !         ", dE[erg]=",delta_e*vol_loc*scale_l**3*scale_d*scale_v**2, &
          !         ", dT2[K]=",delta_e/d*(gamma-1.0)*scale_T2,      &
          !         "unew: ", unew(indp(j,ind),1:9),&
          !         ", uold: ", uold(indp(j,ind),1:9)
          ! end if


          ! etot_max = unew(indp(j,ind),1) * sgagn_T2max/scale_T2 / (gamma-1.0) + &
          !      0.5*(unew(indp(j,ind),2)**2+unew(indp(j,ind),3)**2+unew(indp(j,ind),4)**2)/unew(indp(j,ind),1)

          ! !impose maximum energy
          ! unew(indp(j,ind),5) = min(unew(indp(j,ind),5),etot_max)
        end if
        !------------------------------------------------

        do ivar=imetal,nvar
          unew(indp(j,ind),ivar)=unew(indp(j,ind),ivar)-uold(indp(j,ind),ivar)*acc_mass/(d*vol_loc)
        end do

      endif
    end do
  end do
#endif
end subroutine accrete_sink_new_fg
!################################################################
subroutine compute_accretion_rate_fg(write_sinks)
  use pm_commons
  use amr_commons
  use hydro_commons
  use mpi_mod
  implicit none
  logical::write_sinks

  !------------------------------------------------------------------------
  ! This routine computes the accretion rate onto the sink particles based
  ! on the information collected in collect accretion
  ! It also creates output for the sink particle positions
  !------------------------------------------------------------------------

  integer::i,nx_loc,isink
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  real(dp)::factG,d_star,boost,fa_fact
  real(dp)::r2,v2,c2,density,volume,ethermal,dx_min,scale,mgas,rho_inf,v_bondi
  real(dp),dimension(1:3)::velocity
  real(dp),dimension(1:nsinkmax)::dMEDoverdt
  real(dp)::T2_gas,delta_mass_min

  real(dp)::mdot ! KS ADDED

  ! print *, "KS DEBUG (compute_accretion_rate_fg)"

  dt_acc=huge(0._dp)

  ! Gravitational constant
  factG=1d0
  if(cosmo)factG=3d0/8d0/3.1415926*omega_m*aexp

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3d0
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax !*** KS MODIFEID ***
  !dx_min=scale*0.5D0**nlevelmax/aexp
  d_star=n_star/scale_nH   ! (KS MEMO) n_star: 星形成のthreshold density

  ! Compute sink particle accretion rate by averaging contributions from all levels
  do isink=1,nsink
    density=0.d0; volume=0.d0; velocity=0.d0; ethermal=0d0
    do i=levelmin,nlevelmax
      density=density+weighted_density(isink,i)
      ethermal=ethermal+weighted_ethermal(isink,i)
      velocity(1:3)=velocity(1:3)+weighted_momentum(isink,i,1:3)
      volume=volume+weighted_volume(isink,i)
    end do
    mgas=density
    density=density/(volume+tiny(0.0_dp)) !結局、cic-cellのoverlap領域のvolumeで重み付けして平均した密度が得られる (KS MEMO)
    if (density>0)then
      velocity(1:3)=velocity(1:3)/(density*volume+tiny(0.0_dp)) !速度は質量(d*volume)で重み付けした平均 (KS MEMO)
      ethermal=ethermal/(density*volume+tiny(0.0_dp))           !熱エネルギーも同上 (KS MEMO)
      c2=MAX((gamma-1.0)*ethermal,smallc**2)
      c2sink(isink)=c2
      v2=SUM((velocity(1:3)-vsink(isink,1:3))**2)

      ! calculate mdot from gas properties
      if (BH_model == 0) then
        !no accretion
        mdot = 0d0
      else if (BH_model == 1 .or. BH_model == 3 .or. BH_model == 4) then
        !BHL accretion
        call get_BHL_accretion_rate_fg(msink(isink), density, sqrt(c2), sqrt(v2), mdot)
      else if (BH_model == 2) then
        !SR20 formula
        call get_SR20_accretion_rate_fg(msink(isink), density, sqrt(c2), sqrt(v2), mdot)
        if (mdot < 0) then
          call get_BHL_accretion_rate_fg(msink(isink), density, sqrt(c2), sqrt(v2), mdot)
        end if
      end if


      !--------------------------------------------------------------
      ! accretion boost (artificial; should be removed in the future)
      !--------------------------------------------------------------
      mdot = acc_sink_boost * mdot

      ! !---------------------- KS DEBUG ------------- !
      ! msink(isink)  = 1d2 * 2d33 / scale_m ! 100 Msun
      ! density =   1d3 /  scale_nH !1e3 cm-3
      ! c2 = (0.81d5/scale_v)**2 !7km/s
      ! v2 = (60d5/scale_v)**2 !7km/s
      ! call get_BHL_accretion_rate_fg(msink(isink), density, sqrt(c2), sqrt(v2), mdot)
      ! print *, 'BHL: mass,rho,c,v,mdot', msink(isink)*scale_m/2e33, density*scale_nH,&
      !      sqrt(c2)*scale_v/1d5, sqrt(v2)*scale_v/1d5,mdot*scale_m/scale_t/(2e33/3e7)
      ! call get_SR20_accretion_rate_fg(msink(isink), density, sqrt(c2), sqrt(v2), mdot)
      ! print *, 'SR20: mass,rho,c,v,mdot', msink(isink)*scale_m/2e33, density*scale_nH,&
      !      sqrt(c2)*scale_v/1d5, sqrt(v2)*scale_v/1d5,mdot*scale_m/scale_t/(2e33/3e7)
      ! stop
      ! !---------------------- KS DEBUG ------------- !



      ! use other values for lambda depending on the EOS, this is for isothermal
      !dMBHoverdt(isink)=4.*3.1415926*rho_inf*r2*v_bondi
      dMBHoverdt(isink)=mdot
      dMEDoverdt(isink)=4.*3.1415926*6.67d-8*msink(isink)*1.66d-24/(0.1*6.652d-25*3d10)*scale_t


      ! set sink accretion rate to Bondi value
      dMsink_overdt(isink)=dMBHoverdt(isink)

      ! print *, "KS DEBUG (compute_accretion_rate_fg) eddington_limit, dMBHoverdt(isink), dMEDoverdt(isink)", &
      !      eddington_limit, &
      !      dMBHoverdt(isink)*(scale_d*scale_l**3d0/scale_t)/(1.9891d33/(365.*24.*3600.)),&
      !      dMEDoverdt(isink)*(scale_d*scale_l**3d0/scale_t)/(1.9891d33/(365.*24.*3600.))

      ! limit accretion to Eddington rate
      if(eddington_limit)then
        dMsink_overdt(isink)=min(dMBHoverdt(isink),dMEDoverdt(isink))
      endif

      ! store average quantities for accretion and diagnostics
      eps_sink(isink)=ethermal
      rho_gas(isink)=density
      volume_gas(isink)=volume
      vel_gas(isink,1:ndim)=velocity(1:ndim)

      ! make sure, accretion rate is positive
      dMsink_overdt(isink)=max(0.d0,dMsink_overdt(isink))

      ! compute maximum timestep allowed by sink
      if (dMsink_overdt(isink)>0.0)then
        ! make sure that sink doesnt accrete more than the gas mass within one timestep
        dt_acc(isink)=c_acc*mgas/dMsink_overdt(isink)
        ! make sure that sink doesnt accrete more than its own mass within one timestep
        dt_acc(isink)=min(dt_acc(isink),(c_acc*msink(isink)/dMsink_overdt(isink)))
      end if

      !KS DEBUG
      ! if (myid == 1) then
      !    print *, "KS DEBUG (compute_accretion_rate_fg) isink, dMsink_overdt (code_unit), dMsink_overdt, dMBHoverdt, dMEDoverdt", &
      !         isink, dMsink_overdt(isink), &
      !         dMsink_overdt(isink)*(scale_d*scale_l**3d0/scale_t)/(1.9891d33/(365.*24.*3600.)),&
      !         dMBHoverdt(isink)*(scale_d*scale_l**3d0/scale_t)/(1.9891d33/(365.*24.*3600.)),&
      !         dMEDoverdt(isink)*(scale_d*scale_l**3d0/scale_t)/(1.9891d33/(365.*24.*3600.))
      ! end if

      ! KS DEBUG
      ! if (myid == 1) then
      !    !if (isink == 1) then !sink 1に限定
      !    print '(A, 2(A,I3), 9(A,1PE9.2),/,10(A,1PE9.2),/)', 'KS DEBUG (compute_accretion_rate_fg) ',&
      !         'isink = ', isink,&
      !         ', idsink = ', idsink(isink),&
      !         ', mbh[Msun]=',msink(isink)*scale_m/1.989d33, &
      !         ', mgas[Msun]=', mgas*scale_m/1.989d33, &
      !         ', egas[erg]=', ethermal*density*volume, &
      !         ', volume[pc^3]=', volume*(scale_l/3.08568d18)**3, &
      !         ', T2gas[K]=', ethermal*(gamma-1.0)*scale_T2, &
      !         ', rhoinf[cm^3]=', density*scale_nH, &
      !         ', csinf[km/s]=', sqrt(c2)*scale_v/1d5, &
      !         ', vinf[km/s]=', sqrt(v2)*scale_v/1d5, &
      !         ', boost=', acc_sink_boost, &
      !         'mdot[Msun/yr]=', dMsink_overdt(isink)*scale_m/scale_t/(1.989d33/(365.*24.*3600.)),&
      !         ', dt(lmin)[yr]=',dtnew(levelmin) * scale_t / (365.*24.*3600.), &
      !         ', dmass(lmin)[Msun]=', dMsink_overdt(isink)*dtnew(levelmin)*scale_m/1.989d33, &
      !         ', dE(lmin)[erg]=', 0.1 * 0.05 * dMsink_overdt(isink)*dtnew(levelmin)*scale_m * 3d10**2, &
      !         ', dT2(lmin)[K]=', 0.1 * 0.05 * dMsink_overdt(isink)*dtnew(levelmin) * (3d10/scale_v)**2 / (mgas+1d-100) * (gamma-1.0)*scale_T2,&
      !         ', dt(lmax)[yr]=',dtnew(nlevelmax) * scale_t / (365.*24.*3600.), &
      !         ', dmass(lmax)[Msun]=', dMsink_overdt(isink)*dtnew(nlevelmax)*scale_m/1.989d33, &
      !         ', dE(lmax)[erg]=', 0.1 * 0.05 * dMsink_overdt(isink)*dtnew(nlevelmax)*scale_m * 3d10**2, &
      !         ', dT2(lmax)[K]=', 0.1 * 0.05 * dMsink_overdt(isink)*dtnew(nlevelmax) * (3d10/scale_v)**2 /(mgas+1d-100) * (gamma-1.0)*scale_T2,&
      !         ', c_acc =', c_acc
      !    !end if
      ! end if

    end if


  end do

  if (write_sinks)then
    call print_and_save_sink_properties_fg ! KS MODIFIED
  end if

end subroutine compute_accretion_rate_fg
!################################################################
subroutine get_BHL_accretion_rate_fg(mbh, rhoinf, csinf, vinf, mdot_BHL)
  use pm_commons
  use amr_commons
  use hydro_commons
  use mpi_mod
  implicit none

  real(dp), intent(IN) :: mbh, rhoinf, csinf, vinf
  real(dp), intent(OUT) :: mdot_BHL
  real(dp), parameter :: pi=3.1415926
  real(dp)::factG
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3d0

  factG=1d0 ! TO MODIFIED: G=1 in code units when not cosmo (ks stopped here); the convention of the rest of the code
  if(cosmo)factG=3d0/8d0/3.1415926*omega_m*aexp

  !Bondi accretion
  mdot_BHL = 4.*pi*(factG*mbh)**2 * rhoinf / (csinf**2 + vinf**2)**1.5

  ! ! KS DEBUG
  ! if (myid == 1) then
  !    print '(A, 5(A,1PE9.2))', 'KS DEBUG (get_BHL_accretion_rate_fg) ',&
  !         'mbh[Msun]=',mbh*scale_m/1.989d33, &
  !         ', rhoinf[cm^3]=', rhoinf*scale_nH, &
  !         ', csinf[km/s]=', csinf*scale_v/1d5, &
  !         ', vinf[km/s]=', vinf*scale_v/1d5, &
  !         ', mdot_BHL[Msun/yr]=', mdot_BHL*scale_m/scale_t/(1.989d33/(365.*24.*3600.))
  ! end if
end subroutine get_BHL_accretion_rate_fg
!################################################################
subroutine get_SR20_accretion_rate_fg(mbh_code, rhoinf_code, cinf_code, vinf_code, mdot_code)
  use pm_commons
  use amr_commons
  use hydro_commons
  use cooling_module, only : X, Y
  use mpi_mod
  implicit none

  real(dp), intent(IN) :: mbh_code, rhoinf_code, cinf_code, vinf_code
  real(dp), intent(OUT) :: mdot_code

  real(dp) :: mbh, rhoinf, cinf, vinf, mdot ! in cgs
  real(dp) :: cI, cII, mdot_bondi, f_SR20
  real(dp) :: vR_cinf_cI, vR_cinf_cII, Delta_val

  real(dp),parameter :: pi  = 3.1415926535897932d0
  real(dp),parameter :: Gnewton = 6.6726d-8
  real(dp),parameter :: mH = 1.6733d-24

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m

  real(dp), parameter:: TI = 1d4     !neutral gas temperature [K] after X-ray pre-heating
  real(dp), parameter:: TII = 4d4    !ionized gas temperature [K] in HII region

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3

  !input variables in cgs units
  mbh = mbh_code * scale_m
  rhoinf = rhoinf_code * scale_nH * mH/X
  cinf = cinf_code * scale_v
  vinf = vinf_code * scale_v

  !cI, cII
  cI = 8.1d5 * sqrt(TI/1d4)   !sound velocity of neutral gas = 8.1 km/s at TI = 1e4 K
  cII = 23d5 * sqrt(TII/4d4)  !sound velocity of ionized gas = 23 km/s at TII = 4e4 K

  !=====================================================================
  ! Safety checks
  !=====================================================================

  ! Check if cinf > cI: SR20 model is not applicable for hot gas
  ! (vR(cinf, cI) would produce NaN from sqrt of negative number)
  if (cinf >= cI) then
    mdot_code = -1.0d0 !TO: fall back to BHL
    return
  end if

  !=====================================================================
  ! Main calculation
  !=====================================================================

  !Bondi rate
  mdot_bondi = 4.0d0 * pi * Gnewton**2 * mbh**2 * rhoinf / cI**3

  !Compute vR values (cinf < cI < cII is guaranteed here)
  vR_cinf_cI = vR(cinf, cI)
  vR_cinf_cII = vR(cinf, cII)

  !Accretion rate normalized by Bondi rate
  if (vinf < vR_cinf_cI) then
    !----- Sub-sonic regime -----
    f_SR20 = cI**3 / cII**5 * (vinf**2 + cinf**2) * (1.0d0 + (vinf/vR_cinf_cI)**2)**(-2.5d0)

  else if (vinf < vR_cinf_cII) then
    !----- D-type regime -----
    f_SR20 = 2.0d0**(-2.5d0) * cI**3 * (vinf**2 + cinf**2) / cII**5

  else
    !----- R-type regime -----
    Delta_val = Delta_wRIF(vinf, cinf, cII)
    f_SR20 = cI**3 / (cII**2 + vinf**2 * Delta_val**(-2))**(1.5d0) * Delta_val
  end if

  ! f_SR20 should be 0 < f_SR20 < 1
  if (f_SR20 > 1.0d0) then
    mdot_code = -1.0d0
    return
  end if

  !Resulting accretion rate
  f_SR20 = max(f_SR20, 0.0d0)  ! 数値誤差で微小な負になった場合の対策
  mdot = mdot_bondi * f_SR20
  mdot_code = mdot / (scale_m / scale_t)

contains

  !-------------------------------------------------------------
  ! Delta = rho2/rho1 = v1/v2 at weak R-type I-front (minus sign)
  ! Safe version: returns -1 if sqrt argument is negative
  !-------------------------------------------------------------
  function Delta_wRIF(v1, c1, c2)
    use mpi_mod
    implicit none
    real(dp), intent(in) :: v1, c1, c2
    real(dp) :: Delta_wRIF
    real(dp) :: sqrt_arg

    sqrt_arg = (c1**2 + v1**2)**2 - 4.0d0 * c2**2 * v1**2
    sqrt_arg = max(sqrt_arg, 0.0d0)  ! 丸め誤差対策：負にならないようにする

    Delta_wRIF = (c1**2 + v1**2 - sqrt(sqrt_arg)) / (2.0d0 * c2**2)
  end function Delta_wRIF

  !-------------------------------------------------------------
  ! vR = c2 + sqrt(c2^2 - c1^2)
  ! Called only when c1 < c2 is guaranteed by the caller
  !-------------------------------------------------------------
  function vR(c1, c2)
    use mpi_mod
    implicit none
    real(dp), intent(in) :: c1, c2
    real(dp) :: vR
    vR = c2 + sqrt(c2**2 - c1**2)
  end function vR

end subroutine get_SR20_accretion_rate_fg
!################################################################
subroutine print_and_save_sink_properties_fg
  use pm_commons
  use amr_commons
  use hydro_commons
  !use output_analysis, only: systemcall_wrapper
  use mpi_mod
  implicit none

  character(LEN=80)::filename,fileloc,nchar,filename_tmp
  character(LEN=80),parameter::prefix_log='logsink'
  integer,parameter :: LUN = 11,  LUN2 = 12

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m,scale_lum
  real(dp) :: zred, vinf, csinf
  integer :: i, isink

  real(dp),parameter :: pc_cgs=3.08568d18 !need to be exactly the same as the value in output_analysis.f90
  real(dp),parameter :: box_cen = 0.5d0 ! boxlen=nx_loc=1 is assumed
  integer, parameter :: dstep_sinklog = 1 !interval of logsink output in coarse step

  !integer,save :: ifirst = 0 !KS ADDED (for prcedures at first call)
  integer,parameter :: nsinkmax_logsink = 10000
  integer,dimension(1:nsinkmax_logsink),save :: ifirst !KS ADDED (for prcedures at first call for each isink)

  !for adjust file
  character(1024) :: line
  integer :: nstep_coarse_line, nstep_line, ios
  logical :: is_exist
  integer :: status
  character(200) :: command, message

  !-------------------------------
  ! for first call
  !-------------------------------
  ! fix logsink[n] by discarding (possibly broken) previous data in the future of the restart time
  !assure nsink is smaller than nsinkmax_logsink
  if (nsink > nsinkmax_logsink) then
    print *, "*** ERROR *** nsink > nsinkmax_logsink ", nsink, nsinkmax_logsink
    print *, "STOP HERE"
    stop
  end if
  if(myid==1) then
    do isink=1, nsink
      if (ifirst(isink)==0) then
        ! get filename
        write(nchar,'(I5)') idsink(isink)
        filename = TRIM(prefix_log)//TRIM(ADJUSTL(nchar))
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
          ! write (LUN,'(A, A)')  '#coarse_step, fine_step, zred, Msink [Msun], Mdot [Msun/yr], L [erg/s], ',&
          !      'ninf [cm-3], vinf [cm-3], csinf [cm-3], vsink(:) [km/s], xsink(:)'
          write (LUN,'(A, A)')  '#coarse_step, fine_step, zred, Msink [Msun], Mdot [Msun/yr], L [erg/s], ',&
            'ninf [cm-3], vinf [cm-3], csinf [cm-3],     vsink(:) [km/s],                    xsink(:) [pc]'
          close(LUN)
        end if

        ifirst(isink) = 1
      end if !ifirst(isink)==0
    end do !isink
  end if !myid==1

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3
  scale_lum = scale_m * scale_v**2 / scale_t

  !--------------------------------------------------------------------
  !output to standard out
  !--------------------------------------------------------------------
  if(myid==1.and.nsink>0)then
    xmsink(1:nsink)=msink(1:nsink)
    call quick_sort_dp(xmsink(1),idsink_sort(1),nsink)
    write(*,*)'Number of sink = ',nsink
    write(*,'(" =========================================================================================================================================================================")')
    !write(*,'(" Id, Msink[Msun], Mdot[Msun/yr],  L[erg/s],    f_Edd,     ninf[cm-3],  csinf[km/s],  vinf[km/s],  dt[yr],  Mdot_now[Msun/yr],           xsink(:)[code_length]")')
    write(*,'(" Id, Msink[Msun], Mdot[Msun/yr],  L[erg/s],    f_Edd,     ninf[cm-3],  csinf[km/s],  vinf[km/s],  vsink0[code_vel],  dt[yr],  Mdot_now[Msun/yr],           xsink(:)[code_length]")')
    write(*,'(" =========================================================================================================================================================================")')
    !do i=nsink,max(nsink-10,1),-1 !top 10 sinks
    do i=nsink,1,-1 !all sinks
      isink=idsink_sort(i)

      !process some vars
      vinf=sqrt(SUM((vel_gas(isink,1:3)-vsink(isink,1:3))**2))
      csinf = sqrt(c2sink(isink))
      write(*,'(I4,10(3X,1PE10.3),3(3X,1PE13.6))')idsink(isink),msink(isink)*scale_m/1.989d33, &
        acc_rate(isink)*scale_m/scale_t/(1.9891d33/(365.*24.*3600.)), &
        acc_lum(isink) * scale_lum, &
        acc_lum(isink) * scale_lum/(msink(isink)*scale_m/1.989d33*1.26e38), & !Ledd = 1.26e38 erg/s * (M/Msun)
        rho_gas(isink) * scale_nH, csinf * scale_v/1d5,  vinf * scale_v/1d5,&
        sqrt(SUM(vsink(isink,1:3))**2), &
        dtnew(levelmin) * scale_t / (365.*24.*3600.),&
        dMsink_overdt(isink)*scale_m/scale_t/(1.9891d33/(365.*24.*3600.)), &
        xsink(isink,1:3)
      !(xsink(isink,1:3)-box_cen)*scale_l/pc_cgs

    end do
    write(*,'(" =========================================================================================================================================================================")')
  endif

  !--------------------------------------------------------------------
  !save sink properties to log files
  !--------------------------------------------------------------------
  if(myid==1.and.nsink>0.and.MOD(nstep_coarse,dstep_sinklog)==0)then

    do isink=1,nsink
      ! get filename
      write(nchar,'(I5)') idsink(isink)
      filename = TRIM(prefix_log)//TRIM(ADJUSTL(nchar))

      !process some vars
      zred = 1./aexp - 1.
      vinf=sqrt(SUM((vel_gas(isink,1:3)-vsink(isink,1:3))**2))
      csinf = sqrt(c2sink(isink))

      !write data to file
      open(LUN, file=filename, position='APPEND')

      ! coarse_step, fine_step
      ! zred, Msink[Msun], Mdot(coarse)[Msun/yr], L[erg/s],
      ! ninf[cm-3], vinf[cm-3], csinf[cm-3],
      ! vsink(:)[km/s], xsink(:)
      write (LUN,'(2I10, 1P1E14.6, 1P1E16.8, 1P8E11.3, 1P3E20.12)') nstep_coarse, nstep, zred, &
        msink(isink)*scale_m/1.989d33, &
        acc_rate(isink) * scale_m/scale_t/(1.9891d33/(365.*24.*3600.)), acc_lum(isink) * scale_lum, &
        rho_gas(isink) * scale_nH, vinf * scale_v/1d5, csinf * scale_v/1d5, &
        vsink(isink,1:3) * scale_v/1d5, xsink(isink,1:3)*scale_l/pc_cgs
      !, xsink(isink,1:3)
      close(LUN)
    end do

  end if


end subroutine print_and_save_sink_properties_fg
!################################################################
subroutine update_sink_fg(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use mpi_mod
  implicit none
  integer::ilevel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine is called at the leafs of the tree structure (right after
! update time). Here is where the global sink variables vsink and xsink are
! updated by summing the conributions from all levels.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(dp)::dteff
  integer::lev,isink


  character(LEN=256)::nchar !KS DEBUG



  !KS DEBUG
  write(nchar,'(I8)') ilevel
  call check_sink_vars_among_cpus_fg('Beginning of update_sink_fg at ilevel = ' // TRIM(ADJUSTL(nchar)))

#if NDIM==3

  if(verbose)write(*,*)'Entering update_sink_fg for level ',ilevel

  fsink=0.
  call f_sink_sink_fg

  ! !---------------- KS DEBUG ---------------!
  ! if (myid == 1) then
  !    do isink=1,nsink
  !       !if (dot_product(fsink(isink,:), fsink(isink,:)) > 1d-12) then
  !       print '(A,2I4,1P3E20.12,2X,L)', "KS DEBUG (update_sink_fg) soon after f_sink_sink_fg: isink, myid, fsink(isink,:), direct =", isink, myid, fsink(isink,:), direct_force_sink(isink)
  !    end do
  ! end if
  ! !------------------------------------------!

  vsold(1:nsink,1:ndim,ilevel)=vsnew(1:nsink,1:ndim,ilevel)
  vsnew(1:nsink,1:ndim,ilevel)=vsink(1:nsink,1:ndim)

  ! Loop over sinks
  do isink=1,nsink

    !*** WARNING *** WARNING *** WARNING *** WARNING *** WARNING *** WARNING ***!
    !------------------------------------- KS ADDED ----------------------------!
    !*** WARNING *** WARNING *** WARNING *** WARNING *** WARNING *** WARNING ***!
    ! BUG FIXED by KS; fsink_partial sometimes remains at a high level even if the grids in that level are all destroyed
    ! erase zombie fsink_partial
    do lev=levelmin,nlevelmax
      if (numbtot(1,lev)==0) then
        fsink_partial(isink,1:ndim,lev) = 0d0
      end if
    end do
    !*** WARNING *** WARNING *** WARNING *** WARNING *** WARNING *** WARNING ***!
    !---------------------------------------------------------------------------!
    !*** WARNING *** WARNING *** WARNING *** WARNING *** WARNING *** WARNING ***!

    ! !---------------- KS DEBUG ---------------!
    ! if (myid == 1) then
    !    print '(A,2I4,2X,L)', "KS DEBUG (update_sink_fg, sum up contribution to fsink) isink, myid, direct_force_sink = ", isink, myid, direct_force_sink(isink)
    !    print '(A,1P3E24.16)', "fsink from sink-sink =", fsink(isink,:)
    !    do lev=levelmin,nlevelmax
    !       if (SUM(fsink_partial(isink,1:ndim,lev)**2) > 0) then
    !          print '(A,I4,1P3E24.16)', "lev, fsink_partial =", lev, fsink_partial(isink,1:ndim,lev)
    !       end if
    !    end do
    ! end if
    ! !------------------------------------------!



    ! sum force contributions from all levels and gather
    do lev=levelmin,nlevelmax
      fsink(isink,1:ndim)=fsink(isink,1:ndim)+fsink_partial(isink,1:ndim,lev)
    end do


    ! !---------------- KS DEBUG ---------------!
    !    if (myid == 1) then
    !       print '(A,2I4,/,1P3E24.16,/,I8,/,1P3E24.16)', "KS DEBUG (update_sink_fg, dividing by ncloud_sink) isink, myid, fsink(isink,:), ncloud_sink, fsink(isink,:)/dble(ncloud_sink) =", isink, myid,   &
    !            fsink(isink,:), ncloud_sink,  fsink(isink,1:ndim)/dble(ncloud_sink)
    !    end if
    !    !------------------------------------------!



    if (.not. direct_force_sink(isink))then
      fsink(isink,1:ndim)=fsink(isink,1:ndim)/dble(ncloud_sink)
    end if

    ! compute timestep for the synchronization
    if (new_born(isink))then
      ! no sync necessary for newly produced sink
      ! advance only 0.5 dt for the first step after formation of sink  (KS MEMO)
      dteff=0d0
    else
      if(sinkint_level>ilevel)then
        ! level at which sinks are integrated has increased
        ! newdt_fine at coarser level has not been computed
        ! -> use dtnew from coarse level
        ! update_sink_fg is called at maximum existing level (KS MEMO)
        ! if it increases call this line (KS MEMO)
        ! dtnew(sinkint_level) is the dt used in the last call of update_sink_fg (KS MEMO)
        dteff=dtnew(sinkint_level)
      else
        ! normal case: finish timestep using dtold from current level.
        ! if ilevel is larger than or equal to sinkint_level (KS MEMO)
        ! dtold(sinkint_level) is the dt used in the last call of update_sink_fg (KS MEMO)
        dteff=dtold(sinkint_level)
      end if
    end if


    ! !---------------- KS DEBUG ---------------!
    ! if (myid == 1) then
    !    print '(A,2I4,2X,L,2(/,1P3E24.16,1PE10.3),/,1P3E24.16,1PE10.3)', "KS DEBUG (update_sink_fg) isink, ilevel, direct_force_sink(isink), vsink(before), vsink(after), fsink, dt =", isink, ilevel, direct_force_sink(isink), &
    !         vsink(isink,1:ndim), sqrt(SUM(vsink(isink,1:3))**2), &
    !         0.5D0*(dtnew(ilevel)+dteff)*fsink(isink,1:ndim)+vsink(isink,1:ndim), &
    !         sqrt(SUM(0.5D0*(dtnew(ilevel)+dteff)*fsink(isink,1:ndim)+vsink(isink,1:ndim))**2), &
    !         fsink(isink,1:ndim),0.5D0*(dtnew(ilevel)+dteff)
    ! end if
    ! !------------------------------------------!


    ! this is the kick-kick (half old half new timestep)
    vsink(isink,1:ndim)=0.5D0*(dtnew(ilevel)+dteff)*fsink(isink,1:ndim)+vsink(isink,1:ndim)

    ! save the velocity
    vsnew(isink,1:ndim,ilevel)=vsink(isink,1:ndim)


    !-----------------------------------------------------------------------
    !---------- BUG ALREADY FIXED IN THE LATEST VERSION OF RAMSES ----------
    !-----------------------------------------------------------------------
    ! ! this is the kick-kick (half old half new timestep)
    ! vsink(isink,1:ndim)=0.5D0*(dtnew(ilevel)+dteff)*fsink(isink,1:ndim)+vsink(isink,1:ndim)

    ! ! save the velocity
    ! vsnew(isink,1:ndim,ilevel)=vsink(isink,1:ndim)
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------

    ! and this is the drift (only for the global sink variable)
    xsink(isink,1:ndim)=xsink(isink,1:ndim)+vsink(isink,1:ndim)*dtnew(ilevel)
    new_born(isink)=.false.

  end do
  ! End loop over sinks

  ! Store deepest level
  sinkint_level=ilevel

#endif

  !KS DEBUG
  call check_sink_vars_among_cpus_fg('End of update_sink_fg at ilevel = ' // TRIM(ADJUSTL(nchar)))

end subroutine update_sink_fg
!################################################################
subroutine update_cloud_fg(ilevel)
  use amr_commons
  use pm_commons
  use mpi_mod
  implicit none
  integer::ilevel

  !----------------------------------------------------------------------
  ! update sink cloud particle properties
  ! -the particles are moved whenever the level of the grid they sit in is updated
  ! -the amount of drift they get is according to their levelp
  ! -since this is happening on the way down, at level ilevel all particles with
  ! level >= ilevel will be moved. Therefore, the sink_jump for all levels >= ilevel
  ! is set ot zero on exit.
  !----------------------------------------------------------------------

  integer::igrid,jgrid,ipart,jpart,next_part,ig,ip,npart1,isink,nx_loc
  integer,dimension(1:nvector)::ind_grid,ind_part,ind_grid_part
  real(dp)::dx,dx_loc,scale

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Update particles position and velocity
  ig=0
  ip=0
  ! Loop over grids
  igrid=headl(myid,ilevel)
  do jgrid=1,numbl(myid,ilevel)
    npart1=numbp(igrid)  ! Number of particles in the grid
    if(npart1>0)then
      ig=ig+1
      ind_grid(ig)=igrid
      ipart=headp(igrid)
      ! Loop over particles
      do jpart=1,npart1
        ! Save next particle  <---- Very important !!!
        next_part=nextp(ipart) !move only particles which do actually belong to that level
        if(ig==0)then
          ig=1
          ind_grid(ig)=igrid
        end if
        ip=ip+1
        ind_part(ip)=ipart
        ind_grid_part(ip)=ig
        if(ip==nvector)then
          call upd_cloud_fg(ind_part,ip)
          ip=0
          ig=0
        end if
        ipart=next_part  ! Go to next particle
      end do
      ! End loop over particles
    end if
    igrid=next(igrid)   ! Go to next grid
  end do
  ! End loop over grids
  if(ip>0)call upd_cloud_fg(ind_part,ip)

  if (myid==1.and.verbose)then
    write(*,*)'sink drift due to accretion relative to grid size at level ',ilevel
    do isink=1,nsink
      write(*,*)'#sink: ',isink,' drift: ',sink_jump(isink,1:ndim,ilevel)/dx_loc
    end do
  end if

  sink_jump(1:nsink,1:ndim,ilevel:nlevelmax)=0.d0

111 format('   Entering update_cloud_fg for level ',I2)

end subroutine update_cloud_fg
!################################################################
subroutine upd_cloud_fg(ind_part,np)
  use amr_commons
  use pm_commons
  use poisson_commons
  use mpi_mod
  implicit none
  integer::np
  integer,dimension(1:nvector)::ind_part

  !------------------------------------------------------------
  ! Vector loop called by update_cloud_fg
  !------------------------------------------------------------

  integer::j,idim,isink,lev
  real(dp),dimension(1:nvector,1:ndim)::new_xp,new_vp
  integer,dimension(1:nvector)::level_p

  ! Overwrite cloud particle mass with sink mass
  do j=1,np
    isink=-idp(ind_part(j))
    if(isink>0 .and. mp(ind_part(j))>0.)then
      mp(ind_part(j))=msink(isink)/dble(ncloud_sink_massive)
    endif
  end do

  ! store velocity
  do idim=1,ndim
    do j=1,np
      new_vp(j,idim)=vp(ind_part(j),idim)
    end do
  end do

  ! Overwrite cloud particle velocity with sink velocity
  ! is going to be overwritten again before move
  do idim=1,ndim
    do j=1,np
      isink=-idp(ind_part(j))
      if(isink>0)then
        new_vp(j,idim)=vsink(isink,idim)
      end if
    end do
  end do

  ! write back velocity
  do idim=1,ndim
    do j=1,np
      vp(ind_part(j),idim)=new_vp(j,idim)
    end do
  end do

  ! read level
  do j=1,np
    level_p(j)=levelp(ind_part(j))
  end do

  ! Update position
  do idim=1,ndim
    do j=1,np
      new_xp(j,idim)=xp(ind_part(j),idim)
    end do
  end do
  do idim=1,ndim
    do j=1,np
      isink=-idp(ind_part(j))
      if(isink>0)then
        lev=level_p(j)
        new_xp(j,idim)=new_xp(j,idim)+sink_jump(isink,idim,lev)
      endif
    end do
  end do

  ! Write back postion
  do idim=1,ndim
    do j=1,np
      xp(ind_part(j),idim)=new_xp(j,idim)
    end do
  end do

end subroutine upd_cloud_fg
!################################################################
subroutine merge_star_sink_fg
  use pm_commons
  use amr_commons
  use mpi_mod
  implicit none

  !------------------------------------------------------------------------
  ! This routine merges sink particles for the star formation case if
  ! they are too close and one of them is younger than ~1000 years
  !------------------------------------------------------------------------

  integer::j,isink,jsink,i,nx_loc,mergers
  real(dp)::dx_loc,scale,dx_min,rr,rmax2,rmax,mnew,t_larson1
  logical::iyoung,jyoung,merge
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::xcom,vcom,lcom

  if(nsink==0)return

  ! Mesh spacing in that level
  dx_loc=0.5D0**nlevelmax
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax !*** KS MODIFEID ***
  !dx_min=scale*0.5D0**nlevelmax/aexp
  rmax=dble(ir_cloud)*dx_min ! Linking length in physical units
  rmax2=rmax*rmax

  !lifetime of first larson core in code units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  t_larson1=merging_timescale*365.25*24*3600/scale_t

  mergers=0
  !loop over all possible pairs (n-square, problematic when there are zilions of sinks)
  do isink=1,nsink-1
    if (msink(isink)>-1.)then
      do jsink=isink+1,nsink
        iyoung=(t-tsink(isink)<t_larson1)
        jyoung=(t-tsink(jsink)<t_larson1)

        rr=(xsink(isink,1)-xsink(jsink,1))**2&
          +(xsink(isink,2)-xsink(jsink,2))**2&
          +(xsink(isink,3)-xsink(jsink,3))**2

        merge=(iyoung .or. jyoung).and.rr<rmax2
        merge=merge .or. (iyoung .and. jyoung .and. rr<4*rmax2)
        merge=merge .and. msink(jsink)>=0


        if (merge)then
          if (myid==1)write(*,*)'merged ', idsink(jsink),' to ',idsink(isink)
          mergers=mergers+1
          mnew=msink(isink)+msink(jsink)
          ! COM
          xcom(1:3)=(msink(isink)*xsink(isink,1:3)+msink(jsink)*xsink(jsink,1:3))/mnew
          vcom(1:3)=(msink(isink)*vsink(isink,1:3)+msink(jsink)*vsink(jsink,1:3))/mnew
          lcom(1:3)=msink(isink)*cross((xsink(isink,1:3)-xcom(1:3)),vsink(isink,1:3)-vcom(1:3))+ &
          & msink(jsink)*cross((xsink(jsink,1:3)-xcom(1:3)),vsink(jsink,1:3)-vcom(1:3))

          xsink(isink,1:3)=(xsink(isink,1:3)*msink(isink)+xsink(jsink,1:3)*msink(jsink))/mnew
          vsink(isink,1:3)=(vsink(isink,1:3)*msink(isink)+vsink(jsink,1:3)*msink(jsink))/mnew
          lsink(isink,1:3)=lsink(isink,1:3)+lsink(jsink,1:3)+lcom(1:3)
          msink(isink)=mnew

          acc_rate(isink)=acc_rate(isink)+msink(jsink)
          !acc_lum is recomputed anyway before its used
          if(ir_feedback)acc_lum(isink)=ir_eff*acc_rate(isink)/dtnew(levelmin)*msink(isink)/(5*6.955d10/scale_l)

          msink(jsink)=-10.
          tsink(isink)=min(tsink(isink),tsink(jsink))
          idsink(isink)=min(idsink(isink),idsink(jsink))
        endif
      end do
    end if
  end do

  if (myid==1 .and. mergers>0)write(*,*)'merged ',mergers,' sinks'

  i=1
  do while (mergers>0)

    if (msink(i)<-1.)then !if sink has been merged to another one

      mergers=mergers-1
      nsink=nsink-1

      !let them all slide back one index
      do j=i,nsink
        xsink(j,1:3)=xsink(j+1,1:3)
        vsink(j,1:3)=vsink(j+1,1:3)
        lsink(j,1:3)=lsink(j+1,1:3)
        msink(j)=msink(j+1)
        mseed(j)=mseed(j+1)
        new_born(j)=new_born(j+1)
        tsink(j)=tsink(j+1)
        idsink(j)=idsink(j+1)
        acc_rate(j)=acc_rate(j+1)
        acc_lum(j)=acc_lum(j+1)
      end do

      !whipe last position in the sink list
      xsink(nsink+1,1:3)=0.
      vsink(nsink+1,1:3)=0.
      lsink(nsink+1,1:3)=0.
      msink(nsink+1)=0.
      mseed(nsink+1)=0.
      new_born(nsink+1)=.false.
      tsink(nsink+1)=0.
      idsink(nsink+1)=0
      acc_rate(nsink+1)=0.
      acc_lum(nsink+1)=0.

    else
      i=i+1
    end if
  end do

end subroutine merge_star_sink_fg
!################################################################
subroutine merge_smbh_sink_fg
  use pm_commons
  use amr_commons
  use hydro_commons,ONLY:mass_sph
  use mpi_mod
  implicit none

  !------------------------------------------------------------------------
  ! This routine merges sink particles for the smbh formation case if
  ! they are too close
  !------------------------------------------------------------------------

  integer::j,isink,jsink,i,nx_loc,mergers
  real(dp)::dx_loc,scale,dx_min,rr,rmax2,rmax,mnew,v1_v2,factG
  logical::merge
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::xcom,vcom,lcom

  if(nsink==0)return

  ! Mesh spacing in that level
  dx_loc=0.5D0**nlevelmax
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  !dx_min=scale*0.5D0**nlevelmax/aexp
  dx_min=scale*0.5D0**nlevelmax !*** KS MODIFEID ***

  rmax=dble(ir_cloud)*dx_min ! Linking length in physical units
  rmax2=rmax*rmax
  ssoft=sink_soft*dx_min

  ! Gravitational constant
  factG=1d0
  if(cosmo)factG=3d0/8d0/3.1415926*omega_m*aexp

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  mergers=0
  !loop over all possible pairs (n-square, problematic when there are zilions of sinks)
  do isink=1,nsink-1
    if (msink(isink)>-1)then
      do jsink=isink+1,nsink

        ! spacing check
        rr=     (xsink(isink,1)-xsink(jsink,1))**2&
        & +(xsink(isink,2)-xsink(jsink,2))**2&
        & +(xsink(isink,3)-xsink(jsink,3))**2

        merge=rr<4*ssoft**2
        merge=merge .and. msink(jsink)>0

        ! escape velocity check
        ! if (mass_merger_vel_check_AGN>0)then
        !    if((msink(isink)+msink(jsink)).ge.mass_merger_vel_check_AGN*2d33/(scale_d*scale_l**3)) then
        !       v1_v2=(vsink(isink,1)-vsink(jsink,1))**2+(vsink(isink,2)-vsink(jsink,2))**2+(vsink(isink,3)-vsink(jsink,3))**2
        !       merge=merge .and. 2*factG*(msink(isink)+msink(jsink))/sqrt(rr)>v1_v2
        !    end if
        ! end if

        !escape velocity check free from some condition about mass (KS MODIFIED)
        v1_v2=(vsink(isink,1)-vsink(jsink,1))**2+(vsink(isink,2)-vsink(jsink,2))**2+(vsink(isink,3)-vsink(jsink,3))**2
        !No velocity condition imposed (KS TODO)
        !merge=merge .and. 2*factG*(msink(isink)+msink(jsink))/sqrt(rr)>v1_v2

        if (merge)then
          !if (myid==1)write(*,*)'merged ', idsink(jsink),' to ',idsink(isink)
          !KS MODIFIED
          if (myid==1) write(*,'(A, /, 2(A,I6), 2(A,1P3E12.4), /, 4(A,1PE12.4))') &
            '** Sink merger **', 'id=', idsink(jsink),' and id=',idsink(isink), &
            ',xi(:)=', xsink(isink,:), ', xj(:)=', xsink(jsink,:),&
            'r[pc]=', sqrt(rr)*scale_l/3.08568d18, &
            ', r_merge[pc]=', 2*ssoft*scale_l/3.08568d18,&
            ', dv[km/s]=', sqrt(v1_v2)*scale_v/1d5, &
            ', v_esc[km/s]=', sqrt(2*factG*(msink(isink)+msink(jsink))/(sqrt(rr)+1d-100))*scale_v/1d5

          mergers=mergers+1
          mnew=msink(isink)+msink(jsink)
          ! COM
          xcom(1:3)=(msink(isink)*xsink(isink,1:3)+msink(jsink)*xsink(jsink,1:3))/mnew
          vcom(1:3)=(msink(isink)*vsink(isink,1:3)+msink(jsink)*vsink(jsink,1:3))/mnew
          lcom(1:3)=msink(isink)*cross((xsink(isink,1:3)-xcom(1:3)),vsink(isink,1:3)-vcom(1:3))+ &
          & msink(jsink)*cross((xsink(jsink,1:3)-xcom(1:3)),vsink(jsink,1:3)-vcom(1:3))

          xsink(isink,1:3)=(xsink(isink,1:3)*msink(isink)+xsink(jsink,1:3)*msink(jsink))/mnew
          vsink(isink,1:3)=(vsink(isink,1:3)*msink(isink)+vsink(jsink,1:3)*msink(jsink))/mnew
          lsink(isink,1:3)=lsink(isink,1:3)+lsink(jsink,1:3)+lcom(1:3)
          msink(isink)=mnew

          acc_rate(isink)=acc_rate(isink)+acc_rate(jsink)
          delta_mass(isink)=delta_mass(isink)+delta_mass(jsink)


          !acc_lum is recomputed anyway before its used
          !if(ir_feedback)acc_lum(isink)=ir_eff*acc_rate(isink)/dtnew(levelmin)*msink(isink)/(5*6.955d10/scale_l)
          acc_lum(isink)=-1d50 !dummy (KS MODIFIED)

          msink(jsink)=-10.
          tsink(isink)=min(tsink(isink),tsink(jsink))
          idsink(isink)=min(idsink(isink),idsink(jsink))
        endif
      end do
    end if
  end do

  if (myid==1 .and. mergers>0)write(*,*)'merged ',mergers,' sinks'

  i=1
  do while (mergers>0)

    if (msink(i)<-1.)then !if sink has been merged to another one

      mergers=mergers-1
      nsink=nsink-1

      !let them all slide back one index
      do j=i,nsink
        xsink(j,1:3)=xsink(j+1,1:3)
        vsink(j,1:3)=vsink(j+1,1:3)
        lsink(j,1:3)=lsink(j+1,1:3)
        msink(j)=msink(j+1)
        mseed(j)=mseed(j+1)
        new_born(j)=new_born(j+1)
        tsink(j)=tsink(j+1)
        idsink(j)=idsink(j+1)
        acc_rate(j)=acc_rate(j+1)
        delta_mass(j)=delta_mass(j+1)
        acc_lum(j)=acc_lum(j+1)
      end do

      !whipe last position in the sink list
      xsink(nsink+1,1:3)=0.
      vsink(nsink+1,1:3)=0.
      lsink(nsink+1,1:3)=0.
      msink(nsink+1)=0.
      mseed(nsink+1)=0.
      new_born(nsink+1)=.false.
      tsink(nsink+1)=0.
      idsink(nsink+1)=0
      acc_rate(nsink+1)=0.
      delta_mass(nsink+1)=0.
      acc_lum(nsink+1)=0.

    else
      i=i+1
    end if
  end do

end subroutine merge_smbh_sink_fg
!################################################################
subroutine f_gas_sink_fg(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use mpi_mod
  implicit none
  integer::ilevel
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! In this subroutine the sink-gas force contributions are calculated.
  ! A plummer-sphere with radius ssoft is used for softening
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz,isink
  integer::info,nx_loc,idim
  real(dp)::dx,dx_loc,scale,vol_loc,dx_min,factG
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc

  logical ,dimension(1:nvector)::ok
  integer ,dimension(1:nvector)::ind_grid,ind_cell
  real(dp),dimension(1:nvector,1:ndim)::xx,ff
  real(dp),dimension(1:nvector)::d2,mcell,denom
  real(dp)::rho_tff,rho_tff_tot,d_min
  logical,dimension(1:ndim)::period

  ! Gravitational constant
  factG=1d0
  if(cosmo)factG=3d0/8d0/3.1415926*omega_m*aexp

  !  Cell spacing at that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim
  !dx_min=scale*0.5D0**nlevelmax/aexp
  dx_min=scale*0.5D0**nlevelmax !*** KS MODIFEID ***
  ssoft=sink_soft*dx_min

  ! Set position of cell centers relative to grid centre
  do ind=1,twotondim
    iz=(ind-1)/4
    iy=(ind-1-4*iz)/2
    ix=(ind-1-2*iy-4*iz)
    if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
    if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
    if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  rho_tff=0.

  fsink_new=0.

  period(1)=(nx==1)
#if NDIM>1
  if(ndim>1)period(2)=(ny==1)
#endif
#if NDIM>2
  if(ndim>2)period(3)=(nz==1)
#endif

  ! Loop over sinks
  do isink=1,nsink
    if (direct_force_sink(isink))then

      d_min=boxlen

      ! Loop over myid grids by vector sweeps
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

          ! Check if cell is refined
          do i=1,ngrid
            ok(i)=son(ind_cell(i))==0
          end do

          ! Gas and dark matter mass in cell
          do i=1,ngrid
            mcell(i)=rho(ind_cell(i))*vol_loc
          end do

          ! Cell center
          do idim=1,ndim
            do i=1,ngrid
              xx(i,idim)=(xg(ind_grid(i),idim)+xc(ind,idim)-skip_loc(idim))*scale
            end do
          end do

          ! Relative position and distance
          d2=0.d0
          do idim=1,ndim
            if (period(idim))then
              do i=1,ngrid
                ! zero order Ewald sum
                ff(i,idim)=xsink(isink,idim)-xx(i,idim)
                if(ff(i,idim)>0.5*boxlen)ff(i,idim)=ff(i,idim)-boxlen
                if(ff(i,idim)<-0.5*boxlen)ff(i,idim)=ff(i,idim)+boxlen
                d2(i)=d2(i)+ff(i,idim)**2
              end do
            else
              do i=1,ngrid
                ff(i,idim)=xsink(isink,idim)-xx(i,idim)
                d2(i)=d2(i)+ff(i,idim)**2
              end do
            end if
          end do

          ! Store minimum distance of cell in current level to isink
          do i=1,ngrid
            d_min=min(d_min,d2(i))
          end do

          ! Compute sqrt(1/(ssoft**2+d2(i))) to save time
          do i=1,ngrid
            denom(i)=(ssoft**2+d2(i))**(-1.5)
          end do

          ! Compute gas acceleration due to sink
          do i=1,ngrid
            ff(i,1:ndim)=denom(i)*ff(i,1:ndim)
          end do

          ! Add gas acceleration due to sink
          do i=1,ngrid
            f(ind_cell(i),1:ndim)=f(ind_cell(i),1:ndim)+factG*msink(isink)*ff(i,1:ndim)
          end do

          ! Add sink acceleration due to gas
          do i=1,ngrid
            if(ok(i))then
              fsink_new(isink,1:ndim)=fsink_new(isink,1:ndim)-factG*mcell(i)*ff(i,1:ndim)
            end if
          end do
        end do !end loop over cells
      end do !end loop over grids

      d_min=d_min**0.5
      d_min=max(ssoft,d_min)
      rho_tff=max(rho_tff,msink(isink)/(4./3.*4.13145*d_min**3)) !(KS NOTE) 4.13145 is probably a bug of pi (fixed in the newer public version and not affecting result mmuch)

    end if !end if direct force
  end do !end loop over sinks

  !collect sink acceleration from cpus
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(fsink_new,fsink_all,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
  fsink_all=fsink_new
#endif
  do isink=1,nsink


    if (direct_force_sink(isink))then
      fsink_partial(isink,1:ndim,ilevel)=fsink_all(isink,1:ndim)

      ! !------------------ KS DEBUG ----------------!
      !       if (myid == 1) then
      !          print '(A, 3I4, 1P9E13.6)', "KS DEBUG (f_gas_sink_fg)  isink, myid, ilevel, fsink_partial, fsink_all, fsink_new: ", &
      !               isink, myid, ilevel, fsink_partial(isink,1:ndim,ilevel), fsink_all(isink,1:ndim), fsink_new(isink,1:ndim)
      !       end if
      ! !------------------ KS DEBUG ----------------!

    end if
  end do


  do idim=1,ndim
    call make_virtual_fine_dp(f(1,idim),ilevel)
  end do


  !collect rho due to sinks for current level - used for timestep computation
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(rho_tff,rho_tff_tot,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
#else
  rho_tff_tot=rho_tff
#endif
  rho_sink_tff(ilevel)=rho_tff_tot



  if (ilevel==nlevelmax)call make_virtual_fine_dp(phi(1),ilevel)

end subroutine f_gas_sink_fg
!################################################################
subroutine f_sink_sink_fg
  use amr_commons
  use pm_commons
  use mpi_mod
  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! In this subroutine the sink-sink force contribution are calculated by direct
  ! n^2 - summation. A plummer-sphere with radius 4 cells is used for softening
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer::isink,idim,jsink
  real(dp),allocatable,dimension(:)::d2
  real(dp),allocatable,dimension(:,:)::ff
  logical,dimension(1:ndim)::period
  real(dp)::factG

  allocate(d2(1:nsink))
  allocate(ff(1:nsink,1:ndim))

  ! Gravitational constant
  factG=1d0
  if(cosmo)factG=3d0/8d0/3.1415926*omega_m*aexp

!  fsink=0.

  period(1)=(nx==1)
#if NDIM>1
  if(ndim>1)period(2)=(ny==1)
#endif
#if NDIM>2
  if(ndim>2)period(3)=(nz==1)
#endif

  do isink=1,nsink
    if (direct_force_sink(isink))then
      d2=0.d0
      ff=0.d0
      do idim=1,ndim
        !compute relative position and and distances
        if (period(idim))then
          do jsink=1,nsink
            if (direct_force_sink(jsink))then
              ff(jsink,idim)=xsink(jsink,idim)-xsink(isink,idim)
              if(ff(jsink,idim)>0.5*boxlen)ff(jsink,idim)=ff(jsink,idim)-boxlen
              if(ff(jsink,idim)<-0.5*boxlen)ff(jsink,idim)=ff(jsink,idim)+boxlen
              d2(jsink)=d2(jsink)+ff(jsink,idim)**2
            end if
          end do
        else
          do jsink=1,nsink
            if (direct_force_sink(jsink))then
              ff(jsink,idim)=xsink(jsink,idim)-xsink(isink,idim)
              d2(jsink)=d2(jsink)+ff(jsink,idim)**2
            end if
          end do
        end if
      end do
      !compute acceleration
      do jsink=1,nsink
        if (direct_force_sink(jsink))then
          ff(jsink,1:ndim)=factG*msink(jsink)/(ssoft**2+d2(jsink))**1.5*ff(jsink,1:ndim)
        end if
      end do
      do jsink=1,nsink
        fsink(isink,1:ndim)=fsink(isink,1:ndim)+ff(jsink,1:ndim)
        if(jsink<0)then
          print*,'This is just a stupid trick to prevent'
          print*,'the compiler from optimizing this loop!'
        end if
      end do
    end if
  end do

  deallocate(ff,d2) ! JP
end subroutine f_sink_sink_fg
!################################################################
subroutine read_sink_params_fg()
  use pm_commons
  use amr_commons
  use mpi_mod
  implicit none

  !------------------------------------------------------------------------
  ! read sink related parameters and perform some 'sanity chekcs'
  !------------------------------------------------------------------------

  real(dp)::dx_min,scale,cty
  integer::nx_loc
  ! upstream &SINK_PARAMS entries that exist in this pm_parameters, plus the FG (patch/ks) ones
  namelist/sink_params/n_sink,rho_sink,d_sink,accretion_scheme,merging_timescale,&
    ir_cloud_massive,sink_soft,mass_sink_direct_force,ir_cloud,nsinkmax,c_acc,create_sinks,mass_sink_seed,&
    eddington_limit,acc_sink_boost,clump_core,verbose_AGN,T2_AGN,cone_opening,mass_halo_AGN,mass_clump_AGN,&
    smbh,agn,&
    bhpop3,BH_model,sink_mass_min_pop3BH,sink_mass_min_pop2BH,no_sinkcloud_refine,sgagn_T2max,sink_check_vars
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

  if(.not.cosmo) call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  nx_loc=(icoarse_max-icoarse_min+1)
  scale = boxlen/dble(nx_loc)

  ! Read namelist file
  rewind(1)
  read(1,NML=sink_params,END=111)
  goto 112
111 if(myid==1)write(*,*)'You did not set up &SINK_PARAMS in the namelist file'
  if(myid==1)write(*,*)'Using default values '
112 rewind(1)

  if (sink .and. (ndim .ne. 3))then
    if(myid==1)write(*,*)'Sink particles are only implemented for 3d sims.'
    print*, ndim
    call clean_stop
  end if

  if (sink .and. (nlevelmax==levelmin))then
    if(myid==1)write(*,*)'sink particles do currently not work in a single-level cartesian grid'
    if(myid==1)write(*,*)'because they need level 1 to be activated.'
    call clean_stop
  end if

  if (create_sinks .and. accretion_scheme=='none')then
    if(myid==1)write(*,*)'formation of new sinks without subsequent accretion is pointless.'
    if(myid==1)write(*,*)'Choose accretion_scheme=...!'
!     call clean_stop
  end if

  if ((create_sinks .or. accretion_scheme .ne. 'none') .and. (.not. hydro))then
    if(myid==1)write(*,*)'sink creation and accretion require hydro to be turned on'
    call clean_stop
  end if

  !check for accretion scheme
  if (accretion_scheme =='bondi')bondi_accretion=.true.
  if (accretion_scheme =='threshold')threshold_accretion=.true.

  ! for sink formation and accretion a threshold must be given
  if (create_sinks .or. (accretion_scheme .ne. 'none'))then

    ! check for threshold
    if (.not. cosmo)then

      if (rho_sink<0. .and. n_sink<0. .and. d_sink>0.) then
        if(myid==1)write(*,*)'Found d_sink! Assuming code units'
      else if (rho_sink>0. .and. n_sink<0. .and. d_sink<0.)then
        if(myid==1)write(*,*)'Found rho_sink! Assuming g/cc'
        d_sink=rho_sink/scale_d
      else if (rho_sink<0. .and. n_sink>0. .and. d_sink<0.)then
        if(myid==1)write(*,*)'Found n_sink! Assuming H/cc'
        d_sink=n_sink/scale_nH
      else if ((rho_sink>0. .and. n_sink>0.) .or. (rho_sink>0. .and. d_sink>0.) .or. (n_sink>0. .and. d_sink>0.))then
        if (myid==1)write(*,*)'Use n_sink [H/cc] OR rho_sink [g/cc] OR d_sink [code_units]'
        call clean_stop
      else
        if(myid==1)write(*,*)'Trying to setting sink threshold such that jeans length at '
        if(myid==1)write(*,*)'max resolution is resolved by 4 cells, assuming isothermal gas'
        if(T2_star==0.)then
          ! TO MODIFIED: d_sink is only used by the clump-based creation (not ported); bhpop3 sinks do not need it
          if(myid==1)write(*,*)'No value for T2_star given: d_sink left unset (not needed for bhpop3 sinks)'
        else
          dx_min=0.5**nlevelmax*scale
          d_sink=T2_star/scale_T2 *3.14159/16./(dx_min**2)
          if(myid==1)write(*,*)'d_sink = ',d_sink
          if(myid==1)write(*,*)'rho_sink = ',d_sink*scale_d
          if(myid==1)write(*,*)'n_sink = ',d_sink*scale_nH
        end if
      end if

    endif
  end if

  if (.not.smbh)then
    if (merging_timescale<0.)then
      if (myid==1)write(*,*)'You chose sink merging on a timescale but did not provide the timescale'
      if (myid==1)write(*,*)'choosing 1000y as lifetime...'
      merging_timescale=1000.
    end if
    cty=scale_t/(365.25*24.*3600.)
    cont_speed=-1./(merging_timescale/cty)
  end if

  ! nol_accretion requires a somewhat smaller timestep per default
  if(c_acc < 0.)c_acc=0.75

  !check for periodic boundary conditions
  if (nx==1 .or. ny==1 .or. nz==1)then
    if (mass_sink_direct_force .ge. 0.)then
      if(myid==1)print*, 'periodic boundaries in combination with '
      if(myid==1)print*, 'direct force sinks are not treated accurately....'
    end if
  end if

  if(mass_sink_direct_force<0.)then
    mass_sink_direct_force=huge(0._dp)
  end if

end subroutine read_sink_params_fg
!################################################################
subroutine cic_get_cells_fg(indp,xx,vol,ok,ind_grid,xpart,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use mpi_mod
  implicit none

  integer::ng,np,ilevel
  integer ,dimension(1:nvector)::ind_grid,ind_grid_part
  real(dp),dimension(1:nvector,1:ndim)::xpart
  real(dp),dimension(1:nvector,1:twotondim)::vol
  integer ,dimension(1:nvector,1:twotondim)::indp
  logical ,dimension(1:nvector,1:twotondim)::ok
  real(dp),dimension(1:nvector,1:ndim,twotondim)::xx
  !------------------------------------------------------------------
  ! This routine returns the CIC cells and volumes for np particles.
  !------------------------------------------------------------------
  logical::error
  integer::i,j,ind,idim,nx_loc,ix,iy,iz
  real(dp)::dx,dx_loc,scale,vol_loc
  ! Grid-based arrays
  integer ,dimension(1:nvector)::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle-based arrays
  real(dp),dimension(1:nvector),save::mmm
  real(dp),dimension(1:nvector),save::ttt=0d0
  real(dp),dimension(1:nvector),save::vol2
  real(dp),dimension(1:nvector,1:ndim),save::x,dd,dg
  integer ,dimension(1:nvector,1:ndim),save::ig,id,igg,igd,icg,icd
  integer ,dimension(1:nvector,1:twotondim),save::igrid,icell,kg
  real(dp),dimension(1:nvector,1:ndim),save::x0
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:ndim)::xc

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

  do ind=1,twotondim
    iz=(ind-1)/4
    iy=(ind-1-4*iz)/2
    ix=(ind-1-2*iy-4*iz)
    xc(ind,1)=(dble(ix)-0.5D0)*dx
#if NDIM>1
    xc(ind,2)=(dble(iy)-0.5D0)*dx
#endif
#if NDIM>2
    xc(ind,3)=(dble(iz)-0.5D0)*dx
#endif
  end do

  ! Lower left corner of 3x3x3 grid-cube  do idim=1,ndim
  ! 3x3x3 father-cell-cubeと言っても良い (KS MEMO)
  ! 6x6x6のcurrent-level cellのcubeにも対応する (KS MEMO)
  do idim=1,ndim
    do i=1,ng
      x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
    end do
  end do

  do i=1,ng
    ind_cell(i)=father(ind_grid(i))
  end do

  ! Gather neighboring father cells (should be present anytime !)
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Rescale particle position at level ilevel
  do idim=1,ndim
    do j=1,np
      x(j,idim)=xpart(j,idim)/scale+skip_loc(idim)
    end do
  end do
  do idim=1,ndim
    do j=1,np
      x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
    end do
  end do
  do idim=1,ndim
    do j=1,np
      x(j,idim)=x(j,idim)/dx
    end do
  end do

  do idim=1,ndim
    do j=1,np
      ! グリッドの中に粒子がいれば、2 < x < 4 の間にくるハズ (KS MEMO)
      ! 粒子は動くので必ずしもグリッド内にいるとは限らない？？ (KS MEMO)
      if (x(j,idim)<0.5 .or. x(j,idim)>5.5)then   !dx^3のcicセルが3x3x3のcubeからはみ出さない条件 (KS MEMO)
        print*,'particle outside allowed boundary for cic_get_cell'
        print*,x(j,1:ndim)
        if (x(j,idim)<0. .or. x(j,idim)>6.)then
          print*,'particle outside allowed 3by3by3 grid cube'
        end if
        call clean_stop
      end if
    end do
  end do

  ! CIC at level ilevel (dd: right cloud boundary; dg: left cloud boundary)
  ! id,ig -> integer; dd,dg -> double (KS MEMO)
  ! a cube defined with dd (right) and dg (left) has 1.0x1.0x1.0 size (KS MEMO)
  ! this cube is probably the CIC cell which intersects with 8 cells (KS MEMO)
  do idim=1,ndim
    do j=1,np
      dd(j,idim)=x(j,idim)+0.5D0            !0.5 < x < 5.5 -> 1.0 < dd < 6.0 (KS MEMO)
      id(j,idim)=dd(j,idim)                 !1 <= id <= 5 (KS MEMO)
      dd(j,idim)=dd(j,idim)-id(j,idim)      !0.0 < dd < 1.0 (KS MEMO)
      dg(j,idim)=1.0D0-dd(j,idim)           !0.0 < dg < 1.0; dd + dg = 1.0 (KS MEMO)
      ig(j,idim)=id(j,idim)-1               !0 <= ig <= 4 (KS MEMO)
    end do
  end do

  ! Compute cloud volumes
#if NDIM==1
  do j=1,np
    vol(j,1)=dg(j,1)
    vol(j,2)=dd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
    vol(j,1)=dg(j,1)*dg(j,2)
    vol(j,2)=dd(j,1)*dg(j,2)
    vol(j,3)=dg(j,1)*dd(j,2)
    vol(j,4)=dd(j,1)*dd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
    vol(j,1)=dg(j,1)*dg(j,2)*dg(j,3)
    vol(j,2)=dd(j,1)*dg(j,2)*dg(j,3)
    vol(j,3)=dg(j,1)*dd(j,2)*dg(j,3)
    vol(j,4)=dd(j,1)*dd(j,2)*dg(j,3)
    vol(j,5)=dg(j,1)*dg(j,2)*dd(j,3)
    vol(j,6)=dd(j,1)*dg(j,2)*dd(j,3)
    vol(j,7)=dg(j,1)*dd(j,2)*dd(j,3)
    vol(j,8)=dd(j,1)*dd(j,2)*dd(j,3)
  end do
#endif

  ! Compute parent grids
  ! 3x3x3のgridのindexに対応 (元々は6x6x6のセルのindexに対応) (KS MEMO)
  do idim=1,ndim
    do j=1,np
      igg(j,idim)=ig(j,idim)/2   !lower-left corner of cic cell, 0 <= igg <= 2 (KS MEMO)
      igd(j,idim)=id(j,idim)/2   !upper-right corner of cic cell, 0 <= igd <= 2 (KS MEMO)
    end do
  end do
#if NDIM==1
  do j=1,np
    kg(j,1)=1+igg(j,1)
    kg(j,2)=1+igd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
    kg(j,1)=1+igg(j,1)+3*igg(j,2)
    kg(j,2)=1+igd(j,1)+3*igg(j,2)
    kg(j,3)=1+igg(j,1)+3*igd(j,2)
    kg(j,4)=1+igd(j,1)+3*igd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
    kg(j,1)=1+igg(j,1)+3*igg(j,2)+9*igg(j,3) !3x3x3 parent-cell index for lower-left corner of cic cell (KS MEMO)
    kg(j,2)=1+igd(j,1)+3*igg(j,2)+9*igg(j,3)
    kg(j,3)=1+igg(j,1)+3*igd(j,2)+9*igg(j,3)
    kg(j,4)=1+igd(j,1)+3*igd(j,2)+9*igg(j,3)
    kg(j,5)=1+igg(j,1)+3*igg(j,2)+9*igd(j,3)
    kg(j,6)=1+igd(j,1)+3*igg(j,2)+9*igd(j,3)
    kg(j,7)=1+igg(j,1)+3*igd(j,2)+9*igd(j,3)
    kg(j,8)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3) !3x3x3 parent-cell index for upper-right corner of cic cell (KS MEMO)
  end do
#endif
  do ind=1,twotondim
    do j=1,np
      igrid(j,ind)=son(nbors_father_cells(ind_grid_part(j),kg(j,ind))) !grid index of cic-cell corners (KS MEMO)
    end do
  end do

  ! Check if particle has escaped to ilevel-1
  ok(1:np,1:twotondim)=.true.
  do ind=1,twotondim
    do j=1,np
      if (igrid(j,ind)==0)then
        ! print*,'particle has partially escaped to ilevel-1'
        ok(j,ind)=.false.
        !parent-cell index of cic-cell corners (in ilevel-1) if son grid does not exist (KS MEMO)
        indp(j,ind)=nbors_father_cells(ind_grid_part(j),kg(j,ind))
!           xx(j,1:ndim,ind)=(xg(ind_grid(ind_grid_part(j)),1:ndim)+ind3(kg(j,ind),1:ndim)*2.*dx-skip_loc(1:ndim))*scale
      end if
    end do
  end do


  ! Compute parent cell position
  ! 各gird中のcellの相対位置、current-level cellでのindexな気がするが (KS MEMO)
  do idim=1,ndim
    do j=1,np
      icg(j,idim)=ig(j,idim)-2*igg(j,idim)
      icd(j,idim)=id(j,idim)-2*igd(j,idim)
    end do
  end do
#if NDIM==1
  do j=1,np
    icell(j,1)=1+icg(j,1)
    icell(j,2)=1+icd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
    icell(j,1)=1+icg(j,1)+2*icg(j,2)
    icell(j,2)=1+icd(j,1)+2*icg(j,2)
    icell(j,3)=1+icg(j,1)+2*icd(j,2)
    icell(j,4)=1+icd(j,1)+2*icd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
    icell(j,1)=1+icg(j,1)+2*icg(j,2)+4*icg(j,3) !2x2x2 cell index within grid of upper-right corner of cic cell (KS MEMO)
    icell(j,2)=1+icd(j,1)+2*icg(j,2)+4*icg(j,3)
    icell(j,3)=1+icg(j,1)+2*icd(j,2)+4*icg(j,3)
    icell(j,4)=1+icd(j,1)+2*icd(j,2)+4*icg(j,3)
    icell(j,5)=1+icg(j,1)+2*icg(j,2)+4*icd(j,3)
    icell(j,6)=1+icd(j,1)+2*icg(j,2)+4*icd(j,3)
    icell(j,7)=1+icg(j,1)+2*icd(j,2)+4*icd(j,3)
    icell(j,8)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
  end do
#endif


  do ind=1,twotondim
    do j=1,np
      if (ok(j,ind))then
        ! Compute parent cell adress for cells in ilevel or ilevel+1
        ! parent-cell index of cic-cell corners in ilevel (KS MEMO)
        indp(j,ind)=ncoarse+(icell(j,ind)-1)*ngridmax+igrid(j,ind)
        ! Check if particles have leaked into level ilevel+1
        ! if so, set ok to false, but read values from split cell
        ok(j,ind)=(son(indp(j,ind))==0)            ! ok is set to flase if ilevel+1 grid exists (KS MEMO)
      end if
    end do
  end do


  ! get coordinate values for cic-cell corners if in leaf cells (KS MEMO)
  do ind=1,twotondim
    do j=1,np
      if (ok(j,ind))then
        xx(j,1:ndim,ind)=(xg(igrid(j,ind),1:ndim)+xc(icell(j,ind),1:ndim)-skip_loc(1:ndim))*scale
      end if
    end do
  end do

end subroutine cic_get_cells_fg
!################################################################
subroutine make_sink_from_pop3_new_fg(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use mpi_mod
  implicit none
  integer::ilevel

  integer ::icpu,index_sink
  integer ::igrid,jgrid,isink,next_part,npart1
  integer ::ntot,ntot_all,info
  integer ,dimension(1:nvector)::ind_grid,ind_cell
  integer ,dimension(1:nvector)::ind_grid_new,ind_cell_new
  integer ,dimension(1:ncpu)::ntot_sink_cpu,ntot_sink_all
  integer::peak_nr,grid
  logical ::ok_free
  real(dp)::d,u,v,w,e,factG,delta_d,v2
  real(dp)::birth_epoch
  real(dp)::dx,dx_loc,scale,vol_loc
  real(dp)::fourpi,threepi2,tff,tsal
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::d_ball,vol,mclump
  real(dp),dimension(1:nlevelmax)::volume
  real(dp),dimension(1:nvar)::z
  real(dp),dimension(1:3)::skip_loc,x,xcell,xpeak
  real(dp),dimension(1:twotondim,1:3)::xc

  integer::ipart,jpart
  real(dp)::birth_time,t0,dteff,scale_msun
  real(dp)::mprog,t_hne
  real(dp)::cal_tpop3

  integer :: nsink_old
  logical :: idp_exist

  !variables for suppression of pop3 formation near pre-existing pop3
  integer :: sp_nnew_loc, sp_nnew_tot
  integer,dimension(1:ncpu) :: sp_nnew_list !list for number of new in each cpu
  integer :: sp_offset !offset between local count and global isink index
  integer :: sp_offset_ind !offset between  local count and global idsink index

  !buffer for MPI send/recv
  integer,dimension(1:ncpu) :: sbufi,rbufi
  real(dp),dimension(1:ncpu) :: sbufr,rbufr


  if(verbose)write(*,*)'entering make_sink_from_pop3_new_fg for level ',ilevel

  ! Set new sink variables to zero
  msink_new=0d0; mseed_new=0d0; xsink_new=0d0; vsink_new=0d0
  tsink_new=0d0; idsink_new=0; !p3idpsink_new=0

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_msun=scale_d*scale_l**3/1.989d33

  !------------------------------------------------
  ! count the number of new sink particles
  !------------------------------------------------
  sp_nnew_loc=0

  ! （たぶん）receptionグリッドに粒子が存在する可能性があるので、以下の方法でループを回す必要がある
  ! Loop over cpus
  do icpu=1,ncpu
    igrid=headl(icpu,ilevel)
    ! Loop over grids
    do jgrid=1,numbl(icpu,ilevel)
      npart1=numbp(igrid)  ! Number of particles in the grid
      ! Count sink and cloud particles
      if(npart1>0)then
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
          ! Save next particle   <--- Very important !!!
          next_part=nextp(ipart)
          if(idp(ipart) > 0)then !idp > 0 -> either star/Pop III/supernova/dead
            !if (typep(ipart)%family==FAM_BH .and. mp(ipart) > 0.) then

            !both Pop II remnant (FAM_DEAD) and Pop III remnant (FAM_BH) are converted to
            !sink particles if their mass exceeds the threshold
            if ((typep(ipart)%family==FAM_BH .and. mp(ipart)*scale_msun > sink_mass_min_pop3BH .and. sink_mass_min_pop3BH >= 0) .or. &
              (typep(ipart)%family==FAM_DEAD .and. mp(ipart)*scale_msun > sink_mass_min_pop2BH .and. sink_mass_min_pop2BH >= 0)) then

              sp_nnew_loc = sp_nnew_loc + 1
            endif
          endif
          ipart=next_part  ! Go to next particle
        end do
      endif
      igrid=next(igrid)   ! Go to next grid
    end do
    ! End loop over grids
  end do
  ! End loop over cpus

  !create sp_nnew_list
  sp_nnew_list(:) = 0
  sp_nnew_list(myid) = sp_nnew_loc
  sbufi(1:ncpu) = sp_nnew_list(:)
  call MPI_ALLREDUCE(sbufi,rbufi,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  sp_nnew_list(:) = rbufi(1:ncpu)
  sp_nnew_tot = sum(sp_nnew_list) !total number of candidates

  !各cpuで形成するsinkのindexに足すべきoffset (isinkに対応する番号とidsinkに対応する番号)
  sp_offset = nsink
  sp_offset_ind = nindsink
  do icpu = 1, myid-1  ! "do i = 1,0" loop is just skipped (KS MEMO)
    sp_offset = sp_offset + sp_nnew_list(icpu)
    sp_offset_ind = sp_offset_ind + sp_nnew_list(icpu)
  end do

  !total nubmer of sink
  nsink_old = nsink
  nsink = nsink + sum(sp_nnew_list)
  nindsink=nindsink + sum(sp_nnew_list)

  !-------------------------------------------
  ! Check wether max number of sink is reached
  !------------------------------------------
  ok_free=(nsink <=nsinkmax)
  if(.not. ok_free)then
    if(myid==1)write(*,*)'global list of sink particles is too long'
    if(myid==1)write(*,*)'New sink particles',sum(sp_nnew_list)
    if(myid==1)write(*,*)'Increase nsinkmax'
    call clean_stop
  end if

  !------------------------------
  ! Create new sink particles
  !------------------------------
  ! Starting identity number
  isink = sp_offset
  index_sink = sp_offset_ind

  ! （たぶん）receptionグリッドに粒子が存在する可能性があるので、以下の方法でループを回す必要がある
  ! Loop over cpus
  do icpu=1,ncpu
    igrid=headl(icpu,ilevel)
    ! Loop over grids
    do jgrid=1,numbl(icpu,ilevel)
      npart1=numbp(igrid)  ! Number of particles in the grid
      ! Count sink and cloud particles
      if(npart1>0)then
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
          ! Save next particle   <--- Very important !!!
          next_part=nextp(ipart)

          if(idp(ipart) > 0)then !idp > 0 -> either star/Pop III/supernova/dead
            !if (typep(ipart)%family==FAM_BH .and. mp(ipart) > 0.) then

            !both Pop II remnant (FAM_DEAD) and Pop III remnant (FAM_BH) are converted to
            !sink particles if their mass exceeds the threshold
            if ((typep(ipart)%family==FAM_BH .and. mp(ipart)*scale_msun > sink_mass_min_pop3BH .and. sink_mass_min_pop3BH >= 0) .or. &
              (typep(ipart)%family==FAM_DEAD .and. mp(ipart)*scale_msun > sink_mass_min_pop2BH .and. sink_mass_min_pop2BH >= 0)) then

              isink=isink+1
              index_sink=index_sink+1

              ! sink properties
              mseed_new(isink)=mp(ipart)
              msink_new(isink)=mp(ipart)
              ! mseed_new(isink)=1d2/scale_msun !KS DEBUG
              ! msink_new(isink)=1d2/scale_msun !KS DEBUG
              idsink_new(isink)=index_sink
              ! p3idpsink_new(isink)=idp(ipart)
              tsink_new(isink)=tp(ipart)
              xsink_new(isink,1:3)=xp(ipart,:)
              vsink_new(isink,1:3)=vp(ipart,:)

              !change the BH particle into a test particle and follow the sink particle hereafter
              mp(ipart) = 0d0



              if (typep(ipart)%family==FAM_DEAD) then
                print '(A,/, (A,I6), 2(A,1PE9.2), 2(A,1P3E10.2),(A,1PE9.2))', '*** New sink particle generated from Pop2 remnant (FAM_DEAD) ***',&
                  'idsink: ', idsink_new(isink),&
                  ', msink: ', msink_new(isink),&
                  ', msink[Msun]: ', msink_new(isink)*scale_msun,&
                  ', xsink: ', xsink_new(isink,:),&
                  ', vsink: ', vsink_new(isink,:),&
                  ', tsink: ', tsink_new(isink)
              else if (typep(ipart)%family==FAM_BH) then
                print '(A,/, (A,I6), 2(A,1PE9.2), 2(A,1P3E10.2),(A,1PE9.2))', '*** New sink particle generated from Pop3 remnant (FAM_BH) ***',&
                  'idsink: ', idsink_new(isink),&
                  ', msink: ', msink_new(isink),&
                  ', msink[Msun]: ', msink_new(isink)*scale_msun,&
                  ', xsink: ', xsink_new(isink,:),&
                  ', vsink: ', vsink_new(isink,:),&
                  ', tsink: ', tsink_new(isink)
              end if
            endif
          endif
          ipart=next_part  ! Go to next particle
        end do
      endif
      igrid=next(igrid)   ! Go to next grid
    end do
    ! End loop over grids
  end do
  ! End loop over cpus

  ! MPI communication (***_new -> ***_all)
  call MPI_ALLREDUCE(mseed_new ,mseed_all ,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(xsink_new ,xsink_all ,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(vsink_new ,vsink_all ,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(msink_new ,msink_all ,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(tsink_new ,tsink_all ,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(idsink_new,idsink_all,nsinkmax,MPI_INTEGER         ,MPI_SUM,MPI_COMM_WORLD,info)
  !call MPI_ALLREDUCE(p3idpsink_new,p3idpsink_all,nsinkmax,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)

  ! store data to sink particle array
  do isink=nsink_old+1,nsink        !only for new sinks
    msink(isink)=msink_all(isink)
    mseed(isink)=mseed_all(isink)
    xsink(isink,1:ndim)=xsink_all(isink,1:ndim)
    vsink(isink,1:ndim)=vsink_all(isink,1:ndim)
    tsink(isink)=tsink_all(isink)
    idsink(isink)=idsink_all(isink)
    !p3idpsink(isink)=p3idpsink_all(isink)
    new_born(isink)=.false.      !I don't think new_born flag makes significant changes (KS MEMO)

    !define just for consistency with the pre-existing module
    delta_mass(isink)=0d0


    ! if (myid==1) then
    !    print '(A,/, (A,I6), 2(A,1PE9.2), 2(A,1P3E10.2),(A,1PE9.2))', '*** New sink particle generated from Pop3 star ***',&
    !         'idsink: ', idsink(isink),&
    !         ', msink: ', msink(isink),&
    !         ', mseed: ', msink(isink),&
    !         ', xsink: ', xsink(isink,:),&
    !         ', vsink: ', vsink(isink,:),&
    !         ', tsink: ', tsink(isink)
    ! end if
  end do

end subroutine make_sink_from_pop3_new_fg
!################################################################
subroutine check_sink_vars_among_cpus_fg(debug_message)
  use amr_commons
  use pm_commons
  use hydro_commons
  use mpi_mod
  implicit none
  integer :: info

  character(len=*),intent(IN) :: debug_message

  real(dp),dimension(1:3)::xsink_cpu, xsink_ave, vsink_cpu, vsink_ave
  real(dp)::msink_cpu, msink_ave
  integer :: err_flg
  integer isink,idim

  real(dp),parameter :: eps =1d-8 !torrelance

  if(.not.sink_check_vars)return ! &SINK_PARAMS (merge_ks): off by default

  err_flg = 0 !エラーがある場合は値を0以外にする

  do isink=1,nsink

    !msinkのチェック
    msink_cpu=msink(isink) !各cpuが保持している情報
    call MPI_ALLREDUCE(msink_cpu ,msink_ave ,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    msink_ave = msink_ave/ncpu

    if (msink_ave**2+msink_cpu**2 > 0d0) then  !両方ともzeroならok
      if ((msink_cpu - msink_ave)**2 > (msink_ave**2+msink_cpu**2) * eps**2) then  !(ずれの2乗)が(2乗の和) x eps^2より大きければ
        err_flg = 1
      end if
    end if

    !xsinkのチェック
    xsink_cpu(1:3)=xsink(isink,1:3) !各cpuが保持している情報
    call MPI_ALLREDUCE(xsink_cpu ,xsink_ave ,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    xsink_ave(:) = xsink_ave(:)/ncpu
    do idim=1,3
      if (xsink_ave(idim)**2+xsink_cpu(idim)**2 > 0d0) then  !両方ともzeroならok
        if ((xsink_cpu(idim) - xsink_ave(idim))**2 > (xsink_ave(idim)**2+xsink_cpu(idim)**2) * eps**2) then  !(ずれの2乗)が(2乗の和) x eps^2より大きければ
          err_flg = 2
        end if
      end if
    end do


    !vsinkのチェック
    vsink_cpu(1:3)=vsink(isink,1:3) !各cpuが保持している情報
    call MPI_ALLREDUCE(vsink_cpu ,vsink_ave ,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    vsink_ave(:) = vsink_ave(:)/ncpu
    do idim=1,3
      if (vsink_ave(idim)**2+vsink_cpu(idim)**2 > 0d0) then  !両方ともzeroならok
        if ((vsink_cpu(idim) - vsink_ave(idim))**2 > (vsink_ave(idim)**2+vsink_cpu(idim)**2) * eps**2) then  !(ずれの2乗)が(2乗の和) x eps^2より大きければ
          err_flg = 3
        end if
      end if
    end do

    if (err_flg>0) exit
  end do

  if (err_flg>0) then
    if (myid == 1) then
      !デバッグ用メッセージ
      print '(A)', "*** ERROR *** sink properties are not consistent among cpus"
      print '(A)', debug_message
      print '(A,I4,A,I4,A,I4)', "isink/nsink = ", isink, "/", nsink, ", err_flg = ", err_flg
      print '(A)', "stopping..."
    end if

    print '(3I6, 1P7E20.12)', myid, isink, idsink(isink), msink(isink), xsink(isink,1:3), vsink(isink,1:3)

    call clean_stop
  end if

  ! print *, "KS DEBUG (check_sink_vars_among_cpus_fg)", myid

end subroutine check_sink_vars_among_cpus_fg
!################################################################
subroutine BH_luminosity_from_acc_rate_fg()
  use amr_commons
  use pm_commons
  use hydro_commons
  use mpi_mod
  implicit none

  integer :: isink
  real(dp) :: mdot, mass, lum

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m,scale_lum !for test problem

  !loop over sink paritcles
  do isink=1,nsink
    mdot = acc_rate(isink)
    mass = msink(isink)

    if (BH_model == 0 .or. BH_model == 4) then
      lum = 0d0
    else if (BH_model == 1 .or. BH_model == 3) then
      !constant efficiency model
      call simple_BH_luminosity_formula(mass, mdot, lum)
    else if (BH_model == 2) then
      !RIAF/std/slim disk motivated model
      call BH_luminosity_RIAF_Std_Slim(mass, mdot, lum)
    end if

    !store BH luminosity
    acc_lum(isink) = lum

    !test problem with constant luminosity
    ! print *, 'simple: mass, mdot, lum', mass*scale_m/2e33, mdot*scale_m/scale_t/(2e33/3e7), lum * scale_lum
    ! call BH_luminosity_RIAF_Std_Slim(mass, mdot, lum)
    ! print *, 'RIAF_Std_Slim: mass, mdot, lum', mass*scale_m/2e33, mdot*scale_m/scale_t/(2e33/3e7), lum * scale_lum
    ! stop
    ! !---------------------- KS DEBUG ------------- !

  end do


contains
  !---------------------------------
  ! simple luminosity formula
  !---------------------------------
  subroutine simple_BH_luminosity_formula(mass, mdot, lum)
    use amr_commons
    use pm_commons
    use hydro_commons

    use mpi_mod
    implicit none

    real(dp), intent(IN) :: mass, mdot
    real(dp), intent(OUT) :: lum

    real(dp), parameter :: sigma_T = 6.652d-25     ! [cm2]
    real(dp),parameter ::c_light  = 2.9979250d+10 ! [cm/s]

    real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2


    ! Conversion factor from user units to cgs units
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

    !simple 10% efficiency
    lum = 0.1 * mdot * (c_light/scale_v)**2

    !Implement more elaborated model in future (KS MEMO)

  end subroutine simple_BH_luminosity_formula

  !------------------------------------------------------------------
  ! luminosity formula considering RIAF/Std/Slim disk structures
  !
  ! L = eta M c^2
  ! RIAF: eta = 0.1*(Mdot/Medd)                      (Mdot < Medd)
  ! Std:  eta = 0.1                                  (Medd < Mdot < 20 Medd)
  ! Slim: eta = 2 log(Mdot/20 Medd) * (Medd/Mdot)    (20 Medd < Mdot)
  !------------------------------------------------------------------
  subroutine BH_luminosity_RIAF_Std_Slim(mass, mdot, lum)
    use amr_commons
    use pm_commons
    use hydro_commons

    use mpi_mod
    implicit none

    real(dp), intent(IN) :: mass, mdot !in code unit
    real(dp), intent(OUT) :: lum       !in code unit

    real(dp),parameter :: sigma_T = 6.652d-25     ! [cm2]
    real(dp),parameter :: c_light = 2.9979250d+10 ! [cm/s]
    real(dp),parameter :: pi  = 3.1415926
    real(dp),parameter :: Gnewton = 6.6726d-8
    real(dp),parameter :: mH = 1.6733d-24


    real(dp) :: Mdot_edd, L_edd, eta
    real(dp) :: mass_cgs,mdot_cgs

    real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,scale_m


    ! Conversion factor from user units to cgs units
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    scale_m=scale_d*scale_l**3

    !code_unit to cgs
    mass_cgs = mass * scale_m
    mdot_cgs = mdot * scale_m/scale_t

    !Eddington
    L_edd = 4 * pi * Gnewton * mass_cgs  * c_light * mH / sigma_T  ![erg/s]
    Mdot_edd = L_edd/c_light**2                                        ![g/s]

    !luminosity formula
    if (mdot_cgs < Mdot_edd) then
      eta = 0.1d0*(mdot_cgs/Mdot_edd)            !RIAF
    else if (mdot_cgs < 20*Mdot_edd) then
      eta = 0.1d0                            !standard disk
    else !mdot_cgs > 20*Mdot_edd
      eta = 2 * dlog(mdot_cgs/(20*Mdot_edd)) * (Mdot_edd/mdot_cgs)  !slim disk (Watarai+ 2000)
    end if

    !luminosity
    lum = eta * mdot * (c_light/scale_v)**2 !code_unit
  end subroutine BH_luminosity_RIAF_Std_Slim

end subroutine BH_luminosity_from_acc_rate_fg
!################################################################
subroutine init_sink_fg
  use amr_commons
  use pm_commons
  use mpi_mod
  implicit none
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  integer::idim
  integer::i,isink
  integer::ilun,nx_loc
  integer::nsinkold
  real(dp)::xx1,xx2,xx3,vv1,vv2,vv3,mm1,ll1,ll2,ll3
  real(dp),allocatable,dimension(:)::xdp
  integer,allocatable,dimension(:)::isp
  logical,allocatable,dimension(:)::nb
  logical::eof,ic_sink=.false.
  character(LEN=80)::filename
  character(LEN=80)::fileloc
  character(LEN=5)::nchar,ncharcpu

  integer,parameter::tag=1112,tag2=1113
  integer::dummy_io,info2

  !! Weighted divergences have been removed. - JP

  !allocate all sink related quantities...
  allocate(msink(1:nsinkmax))
  allocate(tsink(1:nsinkmax))
  allocate(idsink(1:nsinkmax))
  idsink=0 ! Important: need to set idsink to zero
  allocate(xsink(1:nsinkmax,1:ndim))
  allocate(vsink(1:nsinkmax,1:ndim))
  allocate(vsold(1:nsinkmax,1:ndim,levelmin:nlevelmax))
  allocate(vsnew(1:nsinkmax,1:ndim,levelmin:nlevelmax))
  allocate(fsink_partial(1:nsinkmax,1:ndim,levelmin:nlevelmax))
  fsink_partial=0d0    !KS ADDED
  allocate(fsink(1:nsinkmax,1:ndim))
  allocate(acc_rate(1:nsinkmax))
  acc_rate=0.
  allocate(acc_lum(1:nsinkmax))
  acc_lum=0.
  allocate(dt_acc(1:nsinkmax))
  allocate(rho_sink_tff(levelmin:nlevelmax))
  allocate(lsink(1:nsinkmax,1:3))
  lsink=0.d0
  allocate(level_sink(1:nsinkmax,levelmin:nlevelmax))
  allocate(delta_mass(1:nsinkmax))
  ! Temporary sink variables
  allocate(wden(1:nsinkmax))
  allocate(wmom(1:nsinkmax,1:ndim))
  allocate(weth(1:nsinkmax))
  allocate(wvol(1:nsinkmax))
  allocate(wden_new(1:nsinkmax))
  allocate(wmom_new(1:nsinkmax,1:ndim))
  allocate(weth_new(1:nsinkmax))
  allocate(wvol_new(1:nsinkmax))
  allocate(msink_new(1:nsinkmax))
  allocate(mseed(1:nsinkmax))
  allocate(mseed_new(1:nsinkmax))
  allocate(mseed_all(1:nsinkmax))
  allocate(msink_all(1:nsinkmax))
  allocate(tsink_new(1:nsinkmax))
  allocate(tsink_all(1:nsinkmax))
  allocate(idsink_new(1:nsinkmax))
  allocate(idsink_all(1:nsinkmax))
  allocate(idsink_old(1:nsinkmax))
  allocate(vsink_new(1:nsinkmax,1:ndim))
  allocate(vsink_all(1:nsinkmax,1:ndim))
  allocate(fsink_new(1:nsinkmax,1:ndim))
  allocate(fsink_all(1:nsinkmax,1:ndim))
  allocate(lsink_new(1:nsinkmax,1:3))
  allocate(lsink_all(1:nsinkmax,1:3))
  allocate(xsink_new(1:nsinkmax,1:ndim))
  allocate(xsink_all(1:nsinkmax,1:ndim))
  allocate(sink_jump(1:nsinkmax,1:ndim,levelmin:nlevelmax))
  sink_jump=0.d0
  allocate(level_sink_new(1:nsinkmax,levelmin:nlevelmax))
  allocate(dMsink_overdt(1:nsinkmax))
  allocate(dMBHoverdt(1:nsinkmax))
  allocate(eps_sink(1:nsinkmax))
  eps_sink=0.d0
  allocate(volume_gas(1:nsinkmax))
  volume_gas=0.d0
  allocate(vel_gas(1:nsinkmax,1:ndim))
  vel_gas=0.d0
  allocate(rho_gas(1:nsinkmax))
  rho_gas=0.d0
  allocate(c2sink(1:nsinkmax))
  allocate(weighted_density(1:nsinkmax,1:nlevelmax))
  weighted_density = 0.d0
  allocate(weighted_volume(1:nsinkmax,1:nlevelmax))
  weighted_volume = 0.d0
  allocate(weighted_ethermal(1:nsinkmax,1:nlevelmax))
  weighted_ethermal = 0.d0
  allocate(weighted_momentum(1:nsinkmax,1:nlevelmax,1:ndim))
  weighted_momentum = 0.d0
!  allocate(weighted_divergence(1:nsinkmax,1:nlevelmax))
!  weighted_divergence = 0.d0
  allocate(oksink_new(1:nsinkmax))
  allocate(oksink_all(1:nsinkmax))
  allocate(idsink_sort(1:nsinkmax))
  allocate(xmsink(1:nsinkmax))
  allocate(delta_mass_new(1:nsinkmax),delta_mass_all(1:nsinkmax))
  allocate(ok_blast_agn(1:nsinkmax),ok_blast_agn_all(1:nsinkmax))
  allocate(direct_force_sink(1:nsinkmax))
  allocate(new_born(1:nsinkmax),new_born_all(1:nsinkmax),new_born_new(1:nsinkmax))

  !------------------ KS ADDED (begin) --------------------!
  ! upstream-only arrays touched by shared code (output_sink_csv, newdt_fine with sink_descent)
  allocate(msmbh(1:nsinkmax),dmfsink(1:nsinkmax),graddescent_over_dt(1:nsinkmax))
  msmbh=0d0; dmfsink=0d0; graddescent_over_dt=0d0
  msink=0d0; tsink=0d0; xsink=0d0; vsink=0d0; vsold=0d0; vsnew=0d0; fsink=0d0; dt_acc=huge(0d0)
  delta_mass=0d0; mseed=0d0; dMsink_overdt=0d0; dMBHoverdt=0d0; new_born=.false.; direct_force_sink=.false.
  level_sink=.false.; rho_sink_tff=0d0
  !------------------ KS ADDED (end) ----------------------!

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  ! Compute softening length from minimum cell spacing
  call compute_ncloud_sink

  if(nrestart>0)then
    ilun=4*ncpu+myid+10
    call title(nrestart,nchar)

    if(IOGROUPSIZEREP>0)then
      call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
      fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/sink_'//TRIM(nchar)//'.out'
    else
      fileloc='output_'//TRIM(nchar)//'/sink_'//TRIM(nchar)//'.out'
    endif


    call title(myid,nchar)
    fileloc=TRIM(fileloc)//TRIM(nchar)

    ! Wait for the token
#ifndef WITHOUTMPI
    if(IOGROUPSIZE>0) then
      if (mod(myid-1,IOGROUPSIZE)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
        & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
      end if
    endif
#endif

    open(unit=ilun,file=fileloc,form='unformatted')
    rewind(ilun)
    read(ilun)nsink
    read(ilun)nindsink

    if(nsink>0)then
      allocate(xdp(1:nsink))
      read(ilun)xdp ! Read sink mass
      msink(1:nsink)=xdp
      read(ilun)xdp ! Read sink birth epoch
      tsink(1:nsink)=xdp
      do idim=1,ndim
        read(ilun)xdp ! Read sink position
        xsink(1:nsink,idim)=xdp
      end do
      do idim=1,ndim
        read(ilun)xdp ! Read sink velocity
        vsink(1:nsink,idim)=xdp
      end do
      do idim=1,3
        read(ilun)xdp ! Read sink angular momentum
        lsink(1:nsink,idim)=xdp
      end do
      read(ilun)xdp ! Read sink accumulated rest mass energy
      delta_mass(1:nsink)=xdp
      read(ilun)xdp ! Read sink accretion rate
      acc_rate(1:nsink)=xdp
      deallocate(xdp)
      allocate(isp(1:nsink))
      read(ilun)isp ! Read sink index
      idsink(1:nsink)=isp
      deallocate(isp)
!        read(ilun)ncloud_sink
      allocate(nb(1:nsink))
      read(ilun)nb ! Read newborn boolean
      new_born(1:nsink)=nb
      deallocate(nb)
      read(ilun)sinkint_level
    end if
    close(ilun)
    ! Send the token
#ifndef WITHOUTMPI
    if(IOGROUPSIZE>0) then
      if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
        & MPI_COMM_WORLD,info2)
      end if
    endif
#endif


    call compute_ncloud_sink

    if(ir_feedback)then
      do i=1,nsink
        acc_lum(i)=ir_eff*acc_rate(i)*msink(i)/(5*6.955d10/scale_l)
      end do
    end if
  else  if(nrestart<0)then ! restart from restart files
    ilun=4*ncpu+myid+10
    call title(nrestart,nchar)

    if(IOGROUPSIZEREP>0)then
      call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
      fileloc='output_restart/group_'//TRIM(ncharcpu)//'/sink_restart.out'
    else
      fileloc='output_restart/sink_restart.out'
    endif


    call title(myid,nchar)
    fileloc=TRIM(fileloc)//TRIM(nchar)

    ! Wait for the token
#ifndef WITHOUTMPI
    if(IOGROUPSIZE>0) then
      if (mod(myid-1,IOGROUPSIZE)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
        & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
      end if
    endif
#endif

    open(unit=ilun,file=fileloc,form='unformatted')
    rewind(ilun)
    read(ilun)nsink
    read(ilun)nindsink

    if(nsink>0)then
      allocate(xdp(1:nsink))
      read(ilun)xdp ! Read sink mass
      msink(1:nsink)=xdp
      read(ilun)xdp ! Read sink birth epoch
      tsink(1:nsink)=xdp
      do idim=1,ndim
        read(ilun)xdp ! Read sink position
        xsink(1:nsink,idim)=xdp
      end do
      do idim=1,ndim
        read(ilun)xdp ! Read sink velocity
        vsink(1:nsink,idim)=xdp
      end do
      do idim=1,3
        read(ilun)xdp ! Read sink angular momentum
        lsink(1:nsink,idim)=xdp
      end do
      read(ilun)xdp ! Read sink accumulated rest mass energy
      delta_mass(1:nsink)=xdp
      read(ilun)xdp ! Read sink accretion rate
      acc_rate(1:nsink)=xdp
      deallocate(xdp)
      allocate(isp(1:nsink))
      read(ilun)isp ! Read sink index
      idsink(1:nsink)=isp
      deallocate(isp)
!        read(ilun)ncloud_sink
      allocate(nb(1:nsink))
      read(ilun)nb ! Read newborn boolean
      new_born(1:nsink)=nb
      deallocate(nb)
      read(ilun)sinkint_level
    end if
    close(ilun)
    ! Send the token
#ifndef WITHOUTMPI
    if(IOGROUPSIZE>0) then
      if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
        & MPI_COMM_WORLD,info2)
      end if
    endif
#endif


    call compute_ncloud_sink

    if(ir_feedback)then
      do i=1,nsink
        acc_lum(i)=ir_eff*acc_rate(i)*msink(i)/(5*6.955d10/scale_l)
      end do
    end if
  end if

  if (nrestart>0 .or. nrestart<0)then
    nsinkold=nsink
    if(TRIM(initfile(levelmin)).NE.' ')then
      filename=TRIM(initfile(levelmin))//'/ic_sink_restart'
    else
      filename='ic_sink_restart'
    end if
    INQUIRE(FILE=filename, EXIST=ic_sink)
    if (myid==1)write(*,*)'Looking for file ic_sink_restart: ',filename
    if (.not. ic_sink)then
      filename='ic_sink_restart'
      INQUIRE(FILE=filename, EXIST=ic_sink)
    end if
  else
    nsink=0
    nindsink=0
    nsinkold=0
    if(TRIM(initfile(levelmin)).NE.' ')then
      filename=TRIM(initfile(levelmin))//'/ic_sink'
    else
      filename='ic_sink'
    end if
    INQUIRE(FILE=filename, EXIST=ic_sink)
    if (myid==1)write(*,*)'Looking for file ic_sink: ',filename
    if (.not. ic_sink)then
      filename='ic_sink'
      INQUIRE(FILE=filename, EXIST=ic_sink)
    end if
  end if

  if (ic_sink)then

    ! Wait for the token
#ifndef WITHOUTMPI
    if(IOGROUPSIZE>0) then
      if (mod(myid-1,IOGROUPSIZE)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag2,&
        & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
      end if
    endif
#endif

    open(10,file=filename,form='formatted')
    eof=.false.
    do
      read(10,*,end=102)mm1,xx1,xx2,xx3,vv1,vv2,vv3,ll1,ll2,ll3
      nsink=nsink+1
      nindsink=nindsink+1
      idsink(nsink)=nindsink
      msink(nsink)=mm1
      xsink(nsink,1)=xx1+boxlen/2.0
      xsink(nsink,2)=xx2+boxlen/2.0
      xsink(nsink,3)=xx3+boxlen/2.0
      vsink(nsink,1)=vv1
      vsink(nsink,2)=vv2
      vsink(nsink,3)=vv3
      lsink(nsink,1)=ll1
      lsink(nsink,2)=ll2
      lsink(nsink,3)=ll3
      tsink(nsink)=t
      new_born(nsink)=.true.
    end do
102 continue
    close(10)
    ! Send the token
#ifndef WITHOUTMPI
    if(IOGROUPSIZE>0) then
      if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag2, &
        & MPI_COMM_WORLD,info2)
      end if
    endif
#endif



  end if
  if (myid==1.and.nsink-nsinkold>0)then
    write(*,*)'sinks read from file '//filename
    write(*,'("   Id           M             x             y             z            vx            vy            vz            lx            ly            lz       ")')
    write(*,'("======================================================================================================================================================")')
    do isink=nsinkold+1,nsink
      write(*,'(I8,2X,10(2X,E12.5))')idsink(isink),msink(isink),xsink(isink,1:ndim),&
        vsink(isink,1:ndim),lsink(isink,1:ndim)
    end do
  end if
  do isink=1,nsink
    direct_force_sink(isink)=(msink(isink) .ge. mass_sink_direct_force)
  end do
end subroutine init_sink_fg
!################################################################
!################################################################
!################################################################
!################################################################
!################################################################
subroutine backup_sink_fg(filename)
  use amr_commons
  use pm_commons
  use mpi_mod
  implicit none

  character(LEN=80)::filename

  integer::ilun,idim,i
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  real(dp),allocatable,dimension(:)::xdp
  integer,allocatable,dimension(:)::ii
  logical,allocatable,dimension(:)::nb
  integer,parameter::tag=1135
  integer::dummy_io,info2

  if(.not. sink) return

  if(verbose)write(*,*)'Entering backup_sink_fg'

  ilun=4*ncpu+myid+10

  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)

  ! Wait for the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
    if (mod(myid-1,IOGROUPSIZE)/=0) then
      call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
      & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
    end if
  endif
#endif

  open(unit=ilun,file=TRIM(fileloc),form='unformatted')
  rewind(ilun)

  write(ilun)nsink
  write(ilun)nindsink
  if(nsink>0)then
    allocate(xdp(1:nsink))
    do i=1,nsink
      xdp(i)=msink(i)
    end do
    write(ilun)xdp ! Write sink mass
    do i=1,nsink
      xdp(i)=tsink(i)
    end do
    write(ilun)xdp ! Write sink birth epoch
    do idim=1,ndim
      do i=1,nsink
        xdp(i)=xsink(i,idim)
      end do
      write(ilun)xdp ! Write sink position
    enddo
    do idim=1,ndim
      do i=1,nsink
        xdp(i)=vsink(i,idim)
      end do
      write(ilun)xdp ! Write sink velocity
    enddo
    do idim=1,ndim
      do i=1,nsink
        xdp(i)=lsink(i,idim)
      end do
      write(ilun)xdp ! Write sink angular momentum
    enddo
    do i=1,nsink
      xdp(i)=delta_mass(i)
    end do
    write(ilun)xdp ! Write sink accumulated rest mass energy
    do i=1,nsink
      xdp(i)=acc_rate(i)
    end do
    write(ilun)xdp ! Write sink accretion rate
    deallocate(xdp)
    allocate(ii(1:nsink))
    do i=1,nsink
      ii(i)=idsink(i)
    end do
    write(ilun)ii ! Write sink index
    !        do i=1,nsink
    !           ii(i)=level_sink(i)
    !        end do
    !        write(ilun)ii ! Write sink level
    deallocate(ii)
    !        write(ilun)ncloud_sink   ! Write ncloud
    allocate(nb(1:nsink))
    do i=1,nsink
      nb(i)=new_born(i)
    end do
    write(ilun)nb ! Write level at which sinks where integrated
    deallocate(nb)
    write(ilun)sinkint_level ! Write level at which sinks where integrated
  endif
  close(ilun)

  ! Send the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
    if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
      dummy_io=1
      call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
      & MPI_COMM_WORLD,info2)
    end if
  endif
#endif

end subroutine backup_sink_fg



!################################################################
subroutine output_sink_fg(filename)
  use amr_commons
  use pm_commons
  use mpi_mod
  implicit none
  character(LEN=80)::filename

  integer::i,idim,ipart,isink
  integer::nx_loc,ny_loc,nz_loc,ilun,icpu,idom
  real(dp)::scale,l_abs,rot_period,dx_min
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  character(LEN=80)::fileloc
  character(LEN=5)::nchar

  if(verbose)write(*,*)'Entering output_sink_fg'

  ilun=myid+10

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3d0
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax/aexp

  if(verbose)write(*,*)'Entering output_sink_fg'

  ilun=2*ncpu+myid+10

  fileloc=TRIM(filename)
  open(unit=ilun,file=TRIM(fileloc),form='formatted',status='replace')
  !======================
  ! Write sink properties
  !======================
  write(ilun,*)'Number of sink = ',nsink

  write(ilun,'(" ================================================================================================================================== ")')
  write(ilun,'("        Id       Mass(Msol)             x                y                z               vx               vy               vz      ")')
  write(ilun,'(" ================================================================================================================================== ")')

  do isink=1,nsink
    l_abs=max((lsink(isink,1)**2+lsink(isink,2)**2+lsink(isink,3)**2)**0.5,1.d-50)
    rot_period=32*3.1415*msink(isink)*(dx_min)**2/(5*l_abs+tiny(0.d0))
    write(ilun,'(I10,7(2X,E15.7))')idsink(isink),msink(isink)*scale_m/2d33,xsink(isink,1:ndim),vsink(isink,1:ndim)
  end do
  write(ilun,'(" ================================================================================================================================== ")')
  close(ilun)

end subroutine output_sink_fg
!################################################################
subroutine file_descriptor_sink_fg(filename)
  use amr_commons
  use hydro_commons
  use mpi_mod
  implicit none

  character(LEN=80)::filename
  character(LEN=80)::fileloc
  integer::ivar,ilun

  if(verbose)write(*,*)'Entering file_descriptor_sink_fg'

  ilun=11

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')

  ! Write run parameters
  write(ilun,'("# version: 1")')

  write(ilun,'("")')

  ivar=1
  write(ilun,'(I2,", mass, d")')ivar
  ivar=2
  write(ilun,'(I2,", tsink, d")')ivar
  ivar=3
  write(ilun,'(I2,", position_x, d")')ivar
  ivar=4
  write(ilun,'(I2,", position_y, d")')ivar
  ivar=5
  write(ilun,'(I2,", position_z, d")')ivar
  ivar=6
  write(ilun,'(I2,", velocity_x, d")')ivar
  ivar=7
  write(ilun,'(I2,", velocity_y, d")')ivar
  ivar=8
  write(ilun,'(I2,", velocity_z, d")')ivar
  ivar=9
  write(ilun,'(I2,", BH_spin_x, d")')ivar
  ivar=10
  write(ilun,'(I2,", BH_spin_y, d")')ivar
  ivar=11
  write(ilun,'(I2,", BH_spin_z, d")')ivar
  ivar=12
  write(ilun,'(I2,", delta_mass, d")')ivar
  ivar=13
  write(ilun,'(I2,", real_accretion, d")')ivar
  ivar=14
  write(ilun,'(I2,", identifier, i")')ivar

  close(ilun)

end subroutine file_descriptor_sink_fg

!################################################################
subroutine check_cloud_count_fg(tag)
  use amr_commons
  use pm_commons
  use mpi_mod
  implicit none
  character(len=*)::tag
  !----------------------------------------------------------------------
  ! TO ADDED 2026-11-05 (two-tier verification): after create_sink_fg (lists
  ! merged down to levelmin, see check_particle_tree) every sink must have
  ! exactly ncloud_sink cloud particles (idp=-isink) over all cpus. One
  ! summary line; with -DTWO_TIER_CHECK a mismatch aborts.
  !----------------------------------------------------------------------
  integer::i,igrid,ipart,jpart,next_part,isink,info,nbad,nprint,icpu,ncache
  integer,allocatable,dimension(:)::ncl,ncl_tot

  if(.not.sink.or.nsink==0)return
  allocate(ncl(1:nsink),ncl_tot(1:nsink))
  ncl=0
  do icpu=0,ncpu
    if(icpu==0)then
      ncache=active(levelmin)%ngrid
    else
      if(icpu==myid)cycle
      ncache=reception(icpu,levelmin)%ngrid
    end if
    do i=1,ncache
      if(icpu==0)then
        igrid=active(levelmin)%igrid(i)
      else
#ifdef LIGHT_MPI_COMM
        igrid=reception(icpu,levelmin)%pcomm%igrid(i)
#else
        igrid=reception(icpu,levelmin)%igrid(i)
#endif
      end if
      ipart=headp(igrid)
      do jpart=1,numbp(igrid)
        next_part=nextp(ipart)
        if(idp(ipart)<0)then
          isink=-idp(ipart)
          if(isink>=1.and.isink<=nsink)ncl(isink)=ncl(isink)+1
        end if
        ipart=next_part
      end do
    end do
  end do
  ncl_tot=ncl
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ncl,ncl_tot,nsink,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
  nbad=0; nprint=0
  do isink=1,nsink
    if(ncl_tot(isink)/=ncloud_sink)then
      nbad=nbad+1
      if(myid==1.and.nprint<5)then
        nprint=nprint+1
        write(*,'(" TWO_TIER CHECK sink",I6," (id",I8,") has",I7," cloud particles, expected",I7)')isink,idsink(isink),ncl_tot(isink),ncloud_sink
      end if
    end if
  end do
  if(myid==1)write(*,'(" TWO_TIER CHECK clouds (",A,") step",I7,": nsink",I6,", sinks with wrong cloud count",I6)')trim(tag),nstep_coarse,nsink,nbad
  deallocate(ncl,ncl_tot)
#ifdef TWO_TIER_CHECK
  if(nbad>0)then
    if(myid==1)write(*,*)'TWO_TIER CHECK FAILED (clouds)'
    flush(6)
#ifndef WITHOUTMPI
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
    stop
#endif
  end if
#endif
end subroutine check_cloud_count_fg
