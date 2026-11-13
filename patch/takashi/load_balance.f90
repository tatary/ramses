!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine load_balance
  use amr_commons
  use sn_energy_reservoir, only: sn_pending, sn_pending_coolstop !FG feedback (merge_ks)
  use pm_commons
  use hydro_commons, ONLY: nvar
#ifndef WITHOUTMPI
  use hydro_commons, ONLY: uold, pstarold, rho_eq, p_eq
  use poisson_commons, ONLY: phi, f
#ifdef RT
  use rt_hydro_commons, ONLY: nrtvar, rtuold
#endif
#endif
  use bisection
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  !------------------------------------------------
  ! This routine performs parallel load balancing.
  !------------------------------------------------
  integer::igrid,ncache,ilevel,i,ind,info
  integer::idim,ivar,icpu,jcpu,kcpu
  integer::nxny,ix,iy,iz,iskip
  integer,dimension(nlevelmax,3)::comm_buffin,comm_buffout
#endif

  if(ncpu==1)return

#ifndef WITHOUTMPI
  if(myid==1)write(*,*)'Load balancing AMR grid...'

  ! Put all particle in main tree trunk
  if(pic.and.(.not.init))then
     call make_tree_fine(levelmin)
     do ilevel=levelmin-1,1,-1
        call merge_tree_fine(ilevel)
     end do
  endif

  balance=.true.

  if(.not.allocated(nrecv_lb))allocate(nrecv_lb(1:nlevelmax))
  nrecv_lb=0

  if(verbose)then
     write(*,*)'Input mesh structure'
     do ilevel=1,nlevelmax
        if(numbtot(1,ilevel)>0)write(*,999)ilevel,numbtot(1:4,ilevel)
     end do
  end if

  !-------------------------------------------
  ! Compute new cpu map using chosen ordering
  !-------------------------------------------
  call cmp_new_cpu_map

  !------------------------------------------------------
  ! Expand boundaries to account for new mesh partition
  !------------------------------------------------------
  call flag_coarse
  call refine_coarse
  call build_comm(1)
  call make_virtual_fine_int(cpu_map (1),1)
  call make_virtual_fine_int(cpu_map2(1),1)
  do i=1,nlevelmax-1
     call flag_fine(i,2)
     call refine_fine(i)
     call build_comm(i+1)
     call make_virtual_fine_int(cpu_map (1),i+1)
     call make_virtual_fine_int(cpu_map2(1),i+1)
  end do

  !--------------------------------------
  ! Update physical boundary conditions
  !--------------------------------------
  do ilevel=nlevelmax,1,-1
     if(hydro)then
#ifdef SOLVERmhd
        do ivar=1,nvar+3
#else
        do ivar=1,nvar
#endif
           call make_virtual_fine_dp(uold(1,ivar),ilevel)
#ifdef SOLVERmhd
        end do
#else
        end do
#endif
        if(allocated(sn_pending))call make_virtual_fine_dp(sn_pending(1),ilevel) !FG feedback (merge_ks)
        if(allocated(sn_pending_coolstop))call make_virtual_fine_dp(sn_pending_coolstop(1),ilevel)
        if(momentum_feedback>0)then
           call make_virtual_fine_dp(pstarold(1),ilevel)
        endif
        if(strict_equilibrium>0)then
           call make_virtual_fine_dp(rho_eq(1),ilevel)
           call make_virtual_fine_dp(p_eq(1),ilevel)
        endif
        if(simple_boundary)then
           call make_boundary_hydro(ilevel)
        end if
     end if
#ifdef RT
     if(rt)then
        do ivar=1,nrtvar
           call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
        end do
        if(simple_boundary)then
           call rt_make_boundary_hydro(ilevel)
        end if
     endif
#endif
     if(poisson)then
        call make_virtual_fine_dp(phi(1),ilevel)
        do idim=1,ndim
           call make_virtual_fine_dp(f(1,idim),ilevel)
        end do
     end if
  end do

  !--------------------------------------
  ! Rearrange octs between cpus
  !--------------------------------------
  do ilevel=1,nlevelmax
     do icpu=1,ncpu
        if(icpu==myid)then
           ncache=active(ilevel)%ngrid
        else
           ncache=reception(icpu,ilevel)%ngrid
        end if
        ! Disconnect from old linked list
        do i=1,ncache
           if(icpu==myid)then
              igrid=active(ilevel)%igrid(i)
           else
#ifdef LIGHT_MPI_COMM
              igrid=reception(icpu,ilevel)%pcomm%igrid(i)
#else
              igrid=reception(icpu,ilevel)%igrid(i)
#endif
           end if
           kcpu=cpu_map (father(igrid))
           jcpu=cpu_map2(father(igrid))
           if(kcpu.ne.jcpu)then
              if(prev(igrid).ne.0) then
                 if(next(igrid).ne.0)then
                    next(prev(igrid))=next(igrid)
                    prev(next(igrid))=prev(igrid)
                 else
                    next(prev(igrid))=0
                    taill(kcpu,ilevel)=prev(igrid)
                 end if
              else
                 if(next(igrid).ne.0)then
                    prev(next(igrid))=0
                    headl(kcpu,ilevel)=next(igrid)
                 else
                    headl(kcpu,ilevel)=0
                    taill(kcpu,ilevel)=0
                 end if
              end if
              numbl(kcpu,ilevel)=numbl(kcpu,ilevel)-1
           end if
        end do
        ! Connect to new linked list
        do i=1,ncache
           if(icpu==myid)then
              igrid=active(ilevel)%igrid(i)
           else
#ifdef LIGHT_MPI_COMM
              igrid=reception(icpu,ilevel)%pcomm%igrid(i)
#else
              igrid=reception(icpu,ilevel)%igrid(i)
#endif
           end if
           kcpu=cpu_map (father(igrid))
           jcpu=cpu_map2(father(igrid))
           if(kcpu.ne.jcpu.and.jcpu==myid)then
              nrecv_lb(ilevel)=nrecv_lb(ilevel)+1
              nrecv_lb_tot=nrecv_lb_tot+1
           end if
           if(kcpu.ne.jcpu)then
              if(numbl(jcpu,ilevel)>0)then
                 next(igrid)=0
                 prev(igrid)=taill(jcpu,ilevel)
                 next(taill(jcpu,ilevel))=igrid
                 taill(jcpu,ilevel)=igrid
                 numbl(jcpu,ilevel)=numbl(jcpu,ilevel)+1
              else
                 next(igrid)=0
                 prev(igrid)=0
                 headl(jcpu,ilevel)=igrid
                 taill(jcpu,ilevel)=igrid
                 numbl(jcpu,ilevel)=1
              end if
           end if
        end do
     end do
  end do
  !--------------------------------------
  ! Compute new grid number statistics
  !--------------------------------------
  do ilevel=1,nlevelmax
     comm_buffin(ilevel,1)=numbl(myid,ilevel)
     comm_buffin(ilevel,2)=numbl(myid,ilevel)
     comm_buffin(ilevel,3)=numbl(myid,ilevel)
  end do
  call MPI_ALLREDUCE(comm_buffin(1,1),comm_buffout(1,1),nlevelmax,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(comm_buffin(1,2),comm_buffout(1,2),nlevelmax,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(comm_buffin(1,3),comm_buffout(1,3),nlevelmax,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(used_mem        ,used_mem_tot     ,1        ,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
  do ilevel=1,nlevelmax
     numbtot(1,ilevel)=comm_buffout(ilevel,1)
     numbtot(2,ilevel)=comm_buffout(ilevel,2)
     numbtot(3,ilevel)=comm_buffout(ilevel,3)
     numbtot(4,ilevel)=numbtot(1,ilevel)/ncpu
  end do

  !--------------------------------------
  ! Set old cpu map to new cpu map
  !--------------------------------------
  if(ordering/='bisection') then
     bound_key=bound_key2
     if(two_tier_balance)then   ! two-tier: fine-tier boundaries, Ls and footprints
       bound_key_f=bound_key2_f
       two_tier_ls=two_tier_ls2
       call two_tier_set_foot(1,two_tier_nfoot2,foot_lo2,foot_hi2,foot_first2)
     end if
  else
     bisec_cpubox_min=bisec_cpubox_min2
     bisec_cpubox_max=bisec_cpubox_max2
  end if

  nxny=nx*ny
  do iz=kcoarse_min,kcoarse_max
  do iy=jcoarse_min,jcoarse_max
  do ix=icoarse_min,icoarse_max
     ind=1+ix+iy*nx+iz*nxny
     cpu_map(ind)=cpu_map2(ind)
  end do
  end do
  end do
  do ilevel=1,nlevelmax
     ! Build new communicators
     call build_comm(ilevel)
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           cpu_map(active(ilevel)%igrid(i)+iskip)=cpu_map2(active(ilevel)%igrid(i)+iskip)
        end do
     end do
     call make_virtual_fine_int(cpu_map(1),ilevel)
  end do

  if(pic.and.(.not.init))then
     ! Sort particles down to nlevelmax
     do ilevel=1,nlevelmax-1
        call kill_tree_fine(ilevel)
        call virtual_tree_fine(ilevel)
     end do
     call virtual_tree_fine(nlevelmax)
     do ilevel=nlevelmax-1,levelmin,-1
        call merge_tree_fine(ilevel)
     end do
  end if

  !--------------------------------------------
  ! Shrink boundaries around new mesh partition
  !--------------------------------------------
  shrink=.true.
  do i=nlevelmax-1,1,-1
     call flag_fine(i,2)
     call refine_fine(i)
     call build_comm(i+1)
  end do
  call flag_coarse
  call refine_coarse
  call build_comm(1)
  shrink=.false.
  if(two_tier_balance)call two_tier_report      ! grid/ghost statistics per level
#ifdef TWO_TIER_CHECK
  if(pic.and.(.not.init))call check_particle_tree('load_balance')
#endif
  if(two_tier_balance.and.myid==1)write(*,'(" TWO_TIER STATUS step",I7," load balance ok")')nstep_coarse

  balance=.false.

  if(verbose)then
     write(*,*)'Output mesh structure'
     do ilevel=1,nlevelmax
        if(numbtot(1,ilevel)>0)write(*,999)ilevel,numbtot(1:4,ilevel)
     end do
  end if

999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')
#endif

end subroutine load_balance
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_new_cpu_map
  use amr_commons
  use pm_commons
  use bisection
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  !---------------------------------------------------
  ! This routine computes the new cpu map using
  ! the choosen ordering to balance load across cpus.
  !---------------------------------------------------
  integer::igrid,ncell,ncell_loc,ncache,ngrid
  integer::ilevel,i,ind,idim
  integer::nx_loc
  integer::icpu,isub,idom
  integer::nxny,ix,iy,iz,iskip
  integer::ind_long
  integer,dimension(1:nvector),save::ind_grid,ind_cell

  real(dp)::dx,scale
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(kind=8)::incost_tot,local_cost,cell_cost
  real(kind=8),dimension(0:ndomain)::incost_new,incost_old
  integer(kind=8),dimension(1:overload)::npart_sub
  integer(kind=8)::wflag
  integer,dimension(1:overload)::ncell_sub
  real(kind=8),dimension(1:ndomain)::cost_loc,cost_old,cost_new
  real(qdp),dimension(0:ndomain)::bound_key_loc
  integer,dimension(1:nvector),save::dom
  real(qdp),dimension(1:nvector),save::order_min,order_max
  integer,dimension(1:MAXLEVEL),save::niter_cost
  real(kind=8)::tcost,wcost,rate_loc
  real(kind=8),dimension(1:2)::rate_in,rate_out
  real(kind=8),dimension(1:2)::fac_in,fac_out

  real(dp),dimension(1:1,1:ndim),save :: xx_tmp
  integer,dimension(1:1),save :: c_tmp
#ifdef QUADHILBERT
  real(kind=8),dimension(0:ndomain)::bigdbl,bigtmp
#endif

  ! Local constants
  nxny=nx*ny
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)

  ! Compute time step related cost
  if(cost_weighting)then
     niter_cost(levelmin)=1
     if (nlevelmax - levelmin - 1 > 31) write(*,*) 'Warning load_balance: niter_cost may need to become a kind=8 integer'
     do ilevel=levelmin+1,nlevelmax
        niter_cost(ilevel)=nsubcycle(ilevel-1)*niter_cost(ilevel-1)
     end do
     ! Optional power on the level weights (cost_weight_power=1: unchanged):
     ! >1 balances the finest levels more strictly at the price of a less
     ! even grid count (memory) per cpu.
     if(cost_weight_power/=1.0d0)then
        do ilevel=levelmin,nlevelmax
           niter_cost(ilevel)=max(1,nint(dble(niter_cost(ilevel))**cost_weight_power))
        end do
     end if
  else
     niter_cost(levelmin:nlevelmax)=1
  endif

  if(verbose) print *,"Entering cmp_new_cpu_map"

  ! two-tier decomposition has its own routine (below)
  if(two_tier_balance.and.ordering/='bisection')then
    call cmp_new_cpu_map_two_tier
    return
  end if

  if(ordering/='bisection') then      ! begin if not bisection

  !----------------------------------------
  ! Compute cell ordering and cost
  ! for leaf cells with cpu map = myid.
  ! Store cost in flag1 and MAXIMUM
  ! ordering key in hilbert_key of kind=16
  !----------------------------------------
  ncell=0
  npart_sub=0
  ncell_sub=0
  ncell_loc=1
  dx=1.0d0*scale
  do iz=0,nz-1
  do iy=0,ny-1
  do ix=0,nx-1
     ind=1+ix+iy*nx+iz*nxny
     if(cpu_map(ind)==myid.and.son(ind)==0)then
        xx(1,1)=(dble(ix)+0.5d0-dble(icoarse_min))*scale
#if NDIM>1
        xx(1,2)=(dble(iy)+0.5d0-dble(jcoarse_min))*scale
#endif
#if NDIM>2
        xx(1,3)=(dble(iz)+0.5d0-dble(kcoarse_min))*scale
#endif
        call cmp_minmaxorder(xx,order_min,order_max,dx,ncell_loc)
        call cmp_dommap(xx,dom,ncell_loc)
        ncell=ncell+1
        isub=(dom(1)-1)/ncpu+1
        ncell_sub(isub)=ncell_sub(isub)+1
        flag1(ncell)=0
        hilbert_key(ncell)=order_max(1)
     end if
  end do
  end do
  end do
  ! Loop over levels
  do ilevel=1,nlevelmax
     ! Cell size and cell center offset
     dx=0.5d0**ilevel
     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5d0)*dx-dble(icoarse_min)
#if NDIM>1
        xc(ind,2)=(dble(iy)-0.5d0)*dx-dble(jcoarse_min)
#endif
#if NDIM>2
        xc(ind,3)=(dble(iz)-0.5d0)*dx-dble(kcoarse_min)
#endif
     end do
     ! Loop over cpus
     do icpu=1,ncpu
        if(icpu==myid)then
           ncache=active(ilevel)%ngrid
        else
           ncache=reception(icpu,ilevel)%ngrid
        end if
        ! Loop over grids by vector sweeps
        do igrid=1,ncache,nvector
           ! Gather nvector grids
           ngrid=MIN(nvector,ncache-igrid+1)
           if(icpu==myid)then
              do i=1,ngrid
                 ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
              end do
           else
              do i=1,ngrid
#ifdef LIGHT_MPI_COMM
                 ind_grid(i)=reception(icpu,ilevel)%pcomm%igrid(igrid+i-1)
#else
                 ind_grid(i)=reception(icpu,ilevel)%igrid(igrid+i-1)
#endif
              end do
           end if
           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=ind_grid(i)+iskip
              end do
              do idim=1,ndim
              ncell_loc=0
              do i=1,ngrid
              if(cpu_map(ind_cell(i))==myid.and.son(ind_cell(i))==0)then
                 ncell_loc=ncell_loc+1
                 xx(ncell_loc,idim)=(xg(ind_grid(i),idim)+xc(ind,idim))*scale
              end if
              end do
              end do
              if(ncell_loc>0)then
                 call cmp_minmaxorder(xx,order_min,order_max,dx*scale,ncell_loc)
                 call cmp_dommap(xx,dom,ncell_loc)
              end if
              ncell_loc=0
              do i=1,ngrid
                 if(cpu_map(ind_cell(i))==myid.and.son(ind_cell(i))==0)then
                    ncell    =ncell    +1
                    ncell_loc=ncell_loc+1
                    isub=(dom(ncell_loc)-1)/ncpu+1
                    ncell_sub(isub)=ncell_sub(isub)+1
                    flag1(ncell)=8*10 ! Magic number
                    if(pic)then
                       ! Add more load for tracer particles
                       if (tracer .and. ilevel >= tracer_first_balance_levelmin) then
                          flag1(ncell) = flag1(ncell) + numbp(ind_grid(i)) + &
                               tracer_first_balance_part_per_cell
                       else
                          flag1(ncell)=flag1(ncell)+numbp(ind_grid(i))
                       endif
                    end if
                    wflag = flag1(ncell)*niter_cost(ilevel)
                    if (wflag > 2147483647) then
                       write(*,*) ' wrong type for flag1 --> change to integer kind=8: ',wflag
                       stop
                    endif
                    flag1(ncell)=flag1(ncell)*niter_cost(ilevel)
                    npart_sub(isub)=npart_sub(isub)+flag1(ncell)
                    hilbert_key(ncell)=order_max(ncell_loc)
                 end if
              end do
           end do
           ! End loop over cells
        end do
        ! End loop over grids
     end do
     ! End loop over cpus
  end do
  ! End loop over levels

  !------------------------------------------------------------------
  ! Time-based correction (cost_time_weighting): scale the cost of every
  ! cell this rank owns by (measured compute time / cost) relative to the
  ! average over ranks. The measured time is the sum of the compute timers
  ! in cost_timer_total since the previous load balance; it is a clean
  ! measurement only when the previous load balance was followed by the
  ! steps of a full coarse step on this partition (cost_time_valid), so
  ! the first load balance (initial conditions, restart) keeps factor 1.
  ! Clipped to [0.5,2] and relaxed with cost_time_relax so that one bad
  ! measurement cannot swing the partition.
  !------------------------------------------------------------------
  call cost_timer_total(tcost)
  if(cost_weighting.and.cost_time_weighting)then
     wcost=0d0
     do isub=1,overload
        wcost=wcost+dble(npart_sub(isub))
     end do
     rate_in=0d0
     if(cost_time_valid.and.wcost>0d0.and.tcost>cost_time_prev)then
        rate_loc=(tcost-cost_time_prev)/wcost
        rate_in(1)=rate_loc; rate_in(2)=1d0
     else
        rate_loc=0d0
     end if
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(rate_in,rate_out,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
     rate_out=rate_in
#endif
     if(rate_out(2)>0d0.and.rate_loc>0d0)then
        ! Normalise by the average rate over the ranks that have a measurement
        rate_loc=rate_loc/(rate_out(1)/rate_out(2))
        rate_loc=min(max(rate_loc,0.5d0),2.0d0)
        cost_rank_fac=(1d0-cost_time_relax)*cost_rank_fac+cost_time_relax*rate_loc
     else
        cost_rank_fac=1d0
     end if
     fac_in(1)=cost_rank_fac; fac_in(2)=-cost_rank_fac
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(fac_in,fac_out,2,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
#else
     fac_out=fac_in
#endif
     if(myid==1.and.rate_out(2)>0d0)write(*,'(a,f6.3,a,f6.3,a,i5,a)') &
          & ' Load balance: time cost factor min=',-fac_out(2),' max=',fac_out(1), &
          & ' (',nint(rate_out(2)),' ranks measured)'
  else
     cost_rank_fac=1d0
  end if
  cost_time_prev=tcost
  cost_time_valid=.true.

  !------------------------------------------------
  ! Sort ordering key and store new index in flag2
  !------------------------------------------------
  if (ncell>0) call quick_sort(hilbert_key(1),flag2(1),ncell)

  !-----------------------------
  ! Balance cost across cpus
  !-----------------------------
  cost_loc = 0 ! Compute local and global cost
  do isub=1,overload
     cost_loc(myid+(isub-1)*ncpu) = dble(npart_sub(isub))*cost_rank_fac
  end do
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(cost_loc,cost_old,ndomain,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
  incost_tot = 0D0
  incost_old(0) = 0D0
  do idom = 1,ndomain
     incost_tot = incost_tot + cost_old(idom)
     incost_old(idom) = incost_tot
  end do
  incost_new(0) = 0D0
  do idom = 1,ndomain
     cost_new(idom) = incost_tot/dble(ndomain) ! Exact load balancing
     incost_new(idom) = incost_new(idom-1) + cost_new(idom)
  end do

  !-----------------------------
  ! Compute new cpu boundaries
  !-----------------------------
  bound_key_loc=0.0d0; bound_key2=0.0d0
  ncell_loc=0
  do isub=1,overload
     if(ncell_sub(isub)>0)then
        ! First cpu on local domain
        idom=0
        do while(incost_new(idom)<incost_old(myid-1+(isub-1)*ncpu))
           idom=idom+1
           if (idom > ndomain) exit
        end do
        ! Compute Hilbert key at boundaries
        i=idom
        local_cost=incost_old(myid-1+(isub-1)*ncpu)
        do ind_long=1,ncell_sub(isub)
           cell_cost=dble(flag1(flag2(ind_long+ncell_loc)))*cost_rank_fac
           local_cost=local_cost+cell_cost
           if (i > ndomain) exit
           if(incost_new(i)<local_cost)then
              bound_key_loc(i)=hilbert_key(ind_long+ncell_loc)
              i=i+1
           endif
        end do
     end if
     ncell_loc=ncell_loc+ncell_sub(isub)
  end do
#ifndef WITHOUTMPI
#ifdef QUADHILBERT
  bigdbl= real(bound_key_loc,kind=8)
  bigtmp= 0.0d0
  call MPI_ALLREDUCE(bigdbl,bigtmp,ndomain+1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  ! if call to mpi_sum with mpi_type=mpi_real16 is supported by mpi_allreduce we can do:
  !call MPI_ALLREDUCE(bound_key_loc,bound_key2,ndomain+1,MPI_REAL16,MPI_SUM,MPI_COMM_WORLD,info)
  bound_key2         = real(bigtmp,kind=qdp)
#else
  call MPI_ALLREDUCE(bound_key_loc,bound_key2,ndomain+1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#endif
  bound_key2(0)      =order_all_min
  bound_key2(ndomain)=order_all_max

  else     ! doing bisection
     ! update the bisection
     call build_bisection(update=.true.)
  end if   ! end if not bisection

  !----------------------------------------
  ! Compute new cpu map
  !----------------------------------------
  cpu_map2=0
  ncell_loc=1
  do iz=0,nz-1
  do iy=0,ny-1
  do ix=0,nx-1
     ind=1+ix+iy*nx+iz*nxny
     xx(1,1)=(dble(ix)+0.5d0-dble(icoarse_min))*scale
#if NDIM>1
     xx(1,2)=(dble(iy)+0.5d0-dble(jcoarse_min))*scale
#endif
#if NDIM>2
     xx(1,3)=(dble(iz)+0.5d0-dble(kcoarse_min))*scale
#endif
     cpu_map2(ind)=ncpu ! default value

     if(ordering/='bisection') then
        call cmp_ordering(xx,order_max,ncell_loc)
        cpu_map2(ind)=ncpu ! default value
        do idom=1,ndomain
           if( order_max(1).ge.bound_key2(idom-1).and. &
                & order_max(1).lt.bound_key2(idom))then
              cpu_map2(ind)=mod(idom-1,ncpu)+1
           endif
        end do
     else
        xx_tmp(1,:) = xx(1,:)
        call cmp_bisection_cpumap(xx_tmp,c_tmp,1)
        cpu_map2(ind) = c_tmp(1)
     end if
  end do
  end do
  end do
  ! Loop over levels
  do ilevel=1,nlevelmax
     ! Cell size and cell center offset
     dx=0.5d0**ilevel
     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5d0)*dx-dble(icoarse_min)
#if NDIM>1
        xc(ind,2)=(dble(iy)-0.5d0)*dx-dble(jcoarse_min)
#endif
#if NDIM>2
        xc(ind,3)=(dble(iz)-0.5d0)*dx-dble(kcoarse_min)
#endif
     end do
     ncache=active(ilevel)%ngrid
     ! Loop over grids by vector sweeps
     do igrid=1,ncache,nvector
        ! Gather nvector grids
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=ind_grid(i)+iskip
           end do
           do idim=1,ndim
              do i=1,ngrid
                 xx(i,idim)=(xg(ind_grid(i),idim)+xc(ind,idim))*scale
              end do
           end do
           if(ordering/='bisection') then
              if(ngrid>0)call cmp_ordering(xx,order_max,ngrid)
              do i=1,ngrid
                 cpu_map2(ind_cell(i))=ncpu ! default value
                 do idom=1,ndomain
                    if( order_max(i).ge.bound_key2(idom-1).and. &
                         & order_max(i).lt.bound_key2(idom))then
                       cpu_map2(ind_cell(i))=mod(idom-1,ncpu)+1
                    endif
                 end do
              end do
           else
              do i=1,ngrid
                 ! compute cpu_map2 using bisection
                 xx_tmp(1,:) = xx(i,:)
                 call cmp_bisection_cpumap(xx_tmp,c_tmp,1)
                 cpu_map2(ind_cell(i)) = c_tmp(1)
              end do
           endif
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels

  ! Update virtual boundaries for new cpu map
  call make_virtual_coarse_int(cpu_map2(1))
  do ilevel=1,nlevelmax
     call make_virtual_fine_int(cpu_map2(1),ilevel)
  end do

end subroutine cmp_new_cpu_map
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_cpumap(x,c,nn)
  use amr_parameters
  use amr_commons
  use bisection
  implicit none
  integer, intent(in) ::nn
  real(dp),dimension(1:nvector,1:ndim), intent(in)::x
  integer ,dimension(1:nvector), intent(out)::c

  integer::i,idom
  real(qdp),dimension(1:nvector),save::order

  if(ordering /= 'bisection') then
     call cmp_ordering(x,order,nn)
     do i=1,nn
        c(i)=ndomain ! default value
        do idom=1,ndomain
           if(    order(i).ge.bound_key(idom-1).and. &
                & order(i).lt.bound_key(idom  ))then
              c(i)=idom
           endif
        end do
     end do
     do i=1,nn
        c(i)=MOD(c(i)-1,ncpu)+1
!        c(i)=c(i)-((c(i)-1)/ncpu)*ncpu
     end do
  else
     call cmp_bisection_cpumap(x,c,nn)
  end if

end subroutine cmp_cpumap
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_dommap(x,c,nn)
  use amr_parameters
  use amr_commons
  use bisection
  implicit none
  integer ::nn
  integer ,dimension(1:nvector)::c
  real(dp),dimension(1:nvector,1:ndim)::x

  integer::i,idom
  real(qdp),dimension(1:nvector),save::order

  call cmp_ordering(x,order,nn)
  do i=1,nn
     c(i)=ndomain ! default value
     do idom=1,ndomain
        if(    order(i).ge.bound_key(idom-1).and. &
             & order(i).lt.bound_key(idom  ))then
           c(i)=idom
        endif
     end do
  end do

end subroutine cmp_dommap
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_ordering(x,order,nn)
  use amr_parameters
  use amr_commons
  use mpi_mod
  implicit none
  integer ::nn
#ifndef WITHOUTMPI
  integer::info
#endif
  real(dp),dimension(1:nvector,1:ndim)::x
  real(qdp),dimension(1:nvector)::order
  !--------------------------------------------------------
  ! This routine computes the index key of the input cell
  ! according to its position in space and for the chosen
  ! ordering. Position x are in user units.
  !-----------------------------------------------------
  integer,dimension(1:nvector),save::ix,iy,iz
  integer::i,ncode,bit_length,nx_loc
  integer::temp
  real(kind=8)::scale,bscale
#if NDIM>1
  real(kind=8)::xx,yy,xc,yc
#endif
#if NDIM>2
  real(kind=8)::zz,zc
#endif

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)

  if(ordering=='planar')then
     ! Planar domain decomposition
     do i=1,nn
        order(i)=x(i,1)
     end do
  end if

#if NDIM>1
  if(ordering=='angular')then
     ! Angular domain decomposition
     xc=boxlen/2
     yc=boxlen/2
#if NDIM>2
     zc=boxlen/2
#endif
     do i=1,nn
        xx=x(i,1)-xc+1d-10
        yy=x(i,2)-yc
#if NDIM>2
        zz=x(i,3)
#endif
        if(xx>0.)then
           order(i)=atan(yy/xx)+acos(-1d0)/2d0
        else
           order(i)=atan(yy/xx)+acos(-1d0)*3d0/2d0
        endif
#if NDIM>2
        if(zz.gt.zc)order(i)=order(i)+2d0*acos(-1d0)
#endif
     end do
  end if
#endif

  if(ordering=='hilbert')then
     ! Hilbert curve domain decomposition
     bscale=2**(nlevelmax+1)
     ncode=nx_loc*int(bscale)
     bscale=bscale/scale

     temp=ncode
     do bit_length=1,32
        ncode=ncode/2
        if(ncode<=1) exit
     end do
     if(bit_length==32) then
        write(*,*)'Error in cmp_minmaxorder'
#ifndef WITHOUTMPI
        call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
        stop
#endif
     end if

     do i=1,nn
        ix(i)=int(x(i,1)*bscale)
#if NDIM>1
        iy(i)=int(x(i,2)*bscale)
#endif
#if NDIM>2
        iz(i)=int(x(i,3)*bscale)
#endif
     end do

     if(ndim==1)then
        call hilbert1d(ix,order,nn)
     else if(ndim==2)then
        call hilbert2d(ix,iy,order,bit_length,nn)
     else if (ndim==3)then
        call hilbert3d(ix,iy,iz,order,bit_length,nn)
     end if

  end if

end subroutine cmp_ordering
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_minmaxorder(x,order_min,order_max,dx,nn)
  use amr_parameters
  use amr_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  integer ::nn
  integer ::temp
  real(dp)::dx
  real(dp),dimension(1:nvector,1:ndim)::x
  real(qdp),dimension(1:nvector)::order_min,order_max
  !-----------------------------------------------------
  ! This routine computes the minimum and maximum index
  ! key contained in the input cell and for the chosen
  ! ordering.
  !-----------------------------------------------------
  integer,dimension(1:nvector),save::ix,iy,iz
  integer::i,ncode,bit_length,nx_loc

  real(dp)::dxx,dxmin
#if NDIM>1
  real(dp)::theta1,theta2,theta3,theta4
  real(kind=8)::xx,yy,xc,yc
#endif
#if NDIM>2
  real(kind=8)::zz,zc
#endif
  real(kind=8)::scale,bscaleloc,bscale
  real(qdp)::dkey,oneqdp=1

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dxmin=scale/dble(2**nlevelmax)

  if(ordering=='planar')then
     ! Planar domain decomposition
     dxx=0.5d0*dx
     do i=1,nn
        order_min(i)=x(i,1)-dxx
        order_max(i)=x(i,1)+dxx
     end do
  end if

#if NDIM>1
  if(ordering=='angular')then
     ! Angular domain decomposition
     dxx=0.5d0*dx
     xc=boxlen/2
     yc=boxlen/2
#if NDIM>2
     zc=boxlen/2
#endif
     do i=1,nn
        if(dx==boxlen)then
           order_min(i)=0
           order_max(i)=4d0*acos(-1d0)
        else
           ! x- y-
           yy=x(i,2)-yc-dxx
           xx=x(i,1)-xc-dxx
           if(xx.ge.0.)then
              xx=xx+1d-10
              theta1=atan(yy/xx)+acos(-1d0)/2
           else
              xx=xx-1d-10
              theta1=atan(yy/xx)+acos(-1d0)*3d0/2d0
           endif
           ! x+ y-
           xx=x(i,1)-xc+dxx
           if(xx.gt.0.)then
              xx=xx+1d-10
              theta2=atan(yy/xx)+acos(-1d0)/2
           else
              xx=xx-1d-10
              theta2=atan(yy/xx)+acos(-1d0)*3d0/2d0
           endif

           ! x+ y+
           yy=x(i,2)-yc+dxx
           if(xx.gt.0.)then
              xx=xx+1d-10
              theta3=atan(yy/xx)+acos(-1d0)/2d0
           else
              xx=xx-1d-10
              theta3=atan(yy/xx)+acos(-1d0)*3d0/2d0
           endif
           ! x- y+
           xx=x(i,1)-xc-dxx
           if(xx.ge.0.)then
              xx=xx+1d-10
              theta4=atan(yy/xx)+acos(-1d0)/2d0
           else
              xx=xx-1d-10
              theta4=atan(yy/xx)+acos(-1d0)*3d0/2d0
           endif
           order_min(i)=min(theta1,theta2,theta3,theta4)
           order_max(i)=max(theta1,theta2,theta3,theta4)
#if NDIM>2
           zz=x(i,3)
           if(zz.gt.zc)then
              order_min(i)=order_min(i)+2d0*acos(-1d0)
              order_max(i)=order_max(i)+2d0*acos(-1d0)
           endif
#endif
        endif
     end do
  end if
#endif

  if(ordering=='hilbert')then
     ! Hilbert curve domain decomposition
     bscale=2**(nlevelmax+1)
     bscaleloc=2**nlevelmax*dxmin/dx
     ncode=nx_loc*int(bscaleloc)
     bscaleloc=bscaleloc/scale
     bscale   =bscale   /scale

     temp=ncode
     do bit_length=1,32
        ncode=ncode/2
        if(ncode<=1) exit
     end do
     if(bit_length==32) then
        write(*,*)'Error in cmp_minmaxorder'
#ifndef WITHOUTMPI
        call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
        stop
#endif
     end if

     do i=1,nn
        ix(i)=int(x(i,1)*bscaleloc)
#if NDIM>1
        iy(i)=int(x(i,2)*bscaleloc)
#endif
#if NDIM>2
        iz(i)=int(x(i,3)*bscaleloc)
#endif
     end do

     if(ndim==1)then
        call hilbert1d(ix,order_min,nn)
     else if(ndim==2)then
        call hilbert2d(ix,iy,order_min,bit_length,nn)
     else if (ndim==3)then
        call hilbert3d(ix,iy,iz,order_min,bit_length,nn)
     end if

     dkey=(real(bscale,kind=qdp)/real(bscaleloc,kind=qdp))**ndim
     do i=1,nn
        order_max(i)=(order_min(i)+oneqdp)*dkey
        order_min(i)=(order_min(i))*dkey
     end do

  end if

end subroutine cmp_minmaxorder
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine defrag
  use amr_commons
  use pm_commons
  use poisson_commons
  use sn_energy_reservoir, only: sn_pending, reorder_sn_reservoir !FG feedback (merge_ks)
  use hydro_commons
#ifdef RT
  use rt_hydro_commons
#endif
  implicit none
  integer,allocatable::reservoir_order(:) !FG feedback (merge_ks)

  integer::ncache,ngrid2,igridmax,i,igrid,ibound,ilevel
  integer::iskip1,iskip2,igrid1,igrid2,ind1,icell1,icell2
  integer::ind,idim,ivar,istart

  if(verbose)write(*,*)'Defragmenting main memory...'

  if(allocated(sn_pending))allocate(reservoir_order(ngridmax)) !FG feedback (merge_ks)
  ngrid2=0
  igridmax=0
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              if(allocated(sn_pending))reservoir_order(ngrid2+i)=igrid !FG feedback (merge_ks)
              cpu_map2(igrid)=ngrid2+i
              igridmax=max(igridmax,igrid)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do

  if(allocated(sn_pending))call reorder_sn_reservoir(reservoir_order(:ngrid2)) !FG feedback (merge_ks)

  ngrid2=0
  do igrid=1,igridmax
     flag2(igrid)=0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              icell1=father(igrid)
              if(icell1>ncoarse)then
                 ind1=(icell1-ncoarse-1)/ngridmax+1
                 iskip1=ncoarse+(ind1-1)*ngridmax
                 igrid1=(icell1-iskip1)
                 igrid2=cpu_map2(igrid1)
                 iskip2=ncoarse+(ind1-1)*ngridmax
                 icell2=iskip2+igrid2
              else
                 icell2=icell1
              end if
              flag2(ngrid2+i)=icell2
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     father(igrid)=flag2(igrid)
  end do

  do ind=1,twondim
  ngrid2=0
  do igrid=1,igridmax
     flag2(igrid)=0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              icell1=nbor(igrid,ind)
              if(icell1>ncoarse)then
                 ind1=(icell1-ncoarse-1)/ngridmax+1
                 iskip1=ncoarse+(ind1-1)*ngridmax
                 igrid1=(icell1-iskip1)
                 igrid2=cpu_map2(igrid1)
                 iskip2=ncoarse+(ind1-1)*ngridmax
                 icell2=iskip2+igrid2
              else
                 icell2=icell1
              end if
              flag2(ngrid2+i)=icell2
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     nbor(igrid,ind)=flag2(igrid)
  end do
  end do

  do idim=1,ndim
  ngrid2=0
  do igrid=1,igridmax
     hilbert_key(igrid)=0.0D0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              hilbert_key(ngrid2+i)=real(xg(igrid,idim),kind=qdp)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     xg(igrid,idim)=real(hilbert_key(igrid),kind=8)
  end do
  end do

  if(pic)then

  ngrid2=0
  do igrid=1,igridmax
     flag2(igrid)=0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              flag2(ngrid2+i)=headp(igrid)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     headp(igrid)=flag2(igrid)
  end do

  ngrid2=0
  do igrid=1,igridmax
     flag2(igrid)=0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              flag2(ngrid2+i)=tailp(igrid)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     tailp(igrid)=flag2(igrid)
  end do

  ngrid2=0
  do igrid=1,igridmax
     flag2(igrid)=0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              flag2(ngrid2+i)=numbp(igrid)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     numbp(igrid)=flag2(igrid)
  end do

  endif

  do ind=1,twotondim
  iskip2=ncoarse+(ind-1)*ngridmax
  ngrid2=0
  do igrid=1,igridmax
     flag2(igrid)=0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              igrid1=son(iskip2+igrid)
              if(igrid1>0)then
                 igrid2=cpu_map2(igrid1)
              else
                 igrid2=0
              end if
              flag2(ngrid2+i)=igrid2
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     son(iskip2+igrid)=flag2(igrid)
  end do
  end do

  do ind=1,twotondim
  iskip2=ncoarse+(ind-1)*ngridmax
  ngrid2=0
  do igrid=1,igridmax
     flag2(igrid)=0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              flag2(ngrid2+i)=cpu_map(iskip2+igrid)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     cpu_map(iskip2+igrid)=flag2(igrid)
  end do
  end do

  do ind=1,twotondim
  iskip2=ncoarse+(ind-1)*ngridmax
  ngrid2=0
  do igrid=1,igridmax
     flag2(igrid)=0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              flag2(ngrid2+i)=flag1(iskip2+igrid)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     flag1(iskip2+igrid)=flag2(igrid)
  end do
  end do

  if(hydro)then

#ifdef SOLVERmhd
  do ivar=1,nvar+3
#else
  do ivar=1,nvar
#endif
  do ind=1,twotondim
  iskip2=ncoarse+(ind-1)*ngridmax
  ngrid2=0
  do igrid=1,igridmax
     hilbert_key(igrid)=0.0D0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              hilbert_key(ngrid2+i)=real(uold(iskip2+igrid,ivar),kind=qdp)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     uold(iskip2+igrid,ivar)=real(hilbert_key(igrid),kind=8)
  end do
  end do
  end do

  end if

  if(momentum_feedback>0)then

  do ind=1,twotondim
  iskip2=ncoarse+(ind-1)*ngridmax
  ngrid2=0
  do igrid=1,igridmax
     hilbert_key(igrid)=0.0D0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              hilbert_key(ngrid2+i)=real(pstarold(iskip2+igrid),kind=qdp)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     pstarold(iskip2+igrid)=real(hilbert_key(igrid),kind=8)
  end do
  end do

  endif

  if(strict_equilibrium>0)then

  do ind=1,twotondim
  iskip2=ncoarse+(ind-1)*ngridmax
  ngrid2=0
  do igrid=1,igridmax
     hilbert_key(igrid)=0.0D0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              hilbert_key(ngrid2+i)=real(rho_eq(iskip2+igrid),kind=qdp)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     rho_eq(iskip2+igrid)=real(hilbert_key(igrid),kind=8)
  end do
  end do

  do ind=1,twotondim
  iskip2=ncoarse+(ind-1)*ngridmax
  ngrid2=0
  do igrid=1,igridmax
     hilbert_key(igrid)=0.0D0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              hilbert_key(ngrid2+i)=real(p_eq(iskip2+igrid),kind=qdp)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     p_eq(iskip2+igrid)=real(hilbert_key(igrid),kind=8)
  end do
  end do

  endif

#ifdef RT
  if(rt)then

  do ivar=1,nrtvar
  do ind=1,twotondim
  iskip2=ncoarse+(ind-1)*ngridmax
  ngrid2=0
  do igrid=1,igridmax
     hilbert_key(igrid)=0.0D0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              hilbert_key(ngrid2+i)=real(rtuold(iskip2+igrid,ivar),kind=qdp)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     rtuold(iskip2+igrid,ivar)=real(hilbert_key(igrid),kind=8)
  end do
  end do
  end do

  end if
#endif

  if(poisson)then

  do ind=1,twotondim
  iskip2=ncoarse+(ind-1)*ngridmax
  ngrid2=0
  do igrid=1,igridmax
     hilbert_key(igrid)=0.0D0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              hilbert_key(ngrid2+i)=real(phi(iskip2+igrid),kind=qdp)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     phi(iskip2+igrid)=real(hilbert_key(igrid),kind=8)
  end do
  end do

  do idim=1,ndim
  do ind=1,twotondim
  iskip2=ncoarse+(ind-1)*ngridmax
  ngrid2=0
  do igrid=1,igridmax
     hilbert_key(igrid)=0.0D0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              hilbert_key(ngrid2+i)=real(f(iskip2+igrid,idim),kind=qdp)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     f(iskip2+igrid,idim)=real(hilbert_key(igrid),kind=8)
  end do
  end do
  end do

  end if

  ngrid2=0
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           if(ibound<=ncpu)then
              headl(ibound,ilevel)=ngrid2+1
              taill(ibound,ilevel)=ngrid2+ncache
           else
              headb(ibound-ncpu,ilevel)=ngrid2+1
              tailb(ibound-ncpu,ilevel)=ngrid2+ncache
           end if
           prev(ngrid2+1)=0
           do i=2,ncache
              prev(ngrid2+i)=ngrid2+i-1
           end do
           do i=1,ncache-1
              next(ngrid2+i)=ngrid2+i+1
           end do
           next(ngrid2+ncache)=0
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  headf=ngrid2+1
  tailf=ngridmax
  numbf=ngridmax-ngrid2
  prev(headf)=0
  next(tailf)=0
  do i=ngrid2+2,ngridmax
     prev(i)=i-1
  end do
  do i=ngrid2+1,ngridmax-1
     next(i)=i+1
  end do

  do i=1,nlevelmax
     call build_comm(i)
  end do

  ngrid_current=ngrid2

end subroutine defrag
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################

!#########################################################################
!#########################################################################
! two-tier domain decomposition (ported from patch/ks 2026-11-13; design in
! ks NOTES_handover_two_tier_2026-11-05.md, summary in docs/merge_ks.md).
! MPI calls in these routines use mpi_mod like the rest of the patch (ks used
! include 'mpif.h'; gfortran rejects mixed-type implicit calls in one file).
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
! TO ADDED 2026-11-05: two-tier domain decomposition
! (NOTES_handover_two_tier_2026-11-05.md). Grids at level <= two_tier_ls are
! owned according to bound_key (coarse tier, balanced by memory: octs and
! particles), grids at level > two_tier_ls according to bound_key_f (fine
! tier, balanced by cost = octs x fine steps per coarse step). A cell at
! level l stores the owner of its child grid (level l+1), so a cell uses the
! tier of level l+1. The owner of a grid is decided by the Hilbert key of its
! centre (= centre of its father cell), as in cmp_cpumap.
!#########################################################################
subroutine two_tier_find_dom(key,bkey,idom)
  use amr_commons
  use mpi_mod
  implicit none
  real(qdp)::key
  real(qdp),dimension(0:ndomain)::bkey
  integer::idom
  !----------------------------------------------------------------------
  ! Domain idom with bkey(idom-1) <= key < bkey(idom); ndomain if none
  ! (same convention as cmp_cpumap, binary search instead of a scan).
  !----------------------------------------------------------------------
  integer::jlo,jhi,j
  idom=ndomain
  if(key<bkey(0).or.key>=bkey(ndomain))return
  jlo=1; jhi=ndomain
  do while(jhi>jlo)
    j=(jlo+jhi)/2
    if(key<bkey(j))then
      jhi=j
    else
      jlo=j+1
    end if
  end do
  idom=jlo
end subroutine two_tier_find_dom
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine two_tier_cpu_of_key(order,c,nn,lgrid,bkey_c,bkey_f,ls,nfoot,lo,hi,first)
  use amr_commons
  use mpi_mod
  implicit none
  integer::nn,lgrid,ls,nfoot
  real(qdp),dimension(1:nvector)::order
  integer,dimension(1:nvector)::c
  real(qdp),dimension(0:ndomain)::bkey_c,bkey_f
  real(qdp),dimension(1:max(nfoot,1))::lo,hi
  integer,dimension(0:ndomain)::first
  !----------------------------------------------------------------------
  ! Owner cpu of grids at level lgrid whose centre key is order(:), for the
  ! decomposition (bkey_c, bkey_f, ls, footprints). TO MODIFIED 2026-11-06:
  ! a fine-tier key outside the footprint of its fine domain belongs to the
  ! coarse owner of the key (see amr_commons).
  !----------------------------------------------------------------------
  integer::i,idom
  if(lgrid<=ls)then
    do i=1,nn
      call two_tier_find_dom(order(i),bkey_c,idom)
      c(i)=mod(idom-1,ncpu)+1
    end do
  else
    do i=1,nn
      call two_tier_find_dom(order(i),bkey_f,idom)
      if(two_tier_in_foot(order(i),idom,nfoot,lo,hi,first))then
        c(i)=two_tier_fine_cpu(idom)
      else
        call two_tier_find_dom(order(i),bkey_c,idom)
        c(i)=mod(idom-1,ncpu)+1
      end if
    end do
  end if
end subroutine two_tier_cpu_of_key
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_cpumap_lev(x,c,nn,ilevel)
  use amr_parameters
  use amr_commons
  use mpi_mod
  implicit none
  integer ::nn,ilevel
  integer ,dimension(1:nvector)::c
  real(dp),dimension(1:nvector,1:ndim)::x
  !----------------------------------------------------------------------
  ! Owner cpu of the grids at level ilevel centred at x (user units, as in
  ! cmp_cpumap; x is also the centre of the father cell whose cpu_map is
  ! being set). Identical to cmp_cpumap unless two_tier_balance puts level
  ! ilevel in the fine tier (ilevel > two_tier_ls).
  !----------------------------------------------------------------------
  real(qdp),dimension(1:nvector),save::order
  if(.not.two_tier_balance.or.ordering=='bisection'.or.ilevel<=two_tier_ls)then
    call cmp_cpumap(x,c,nn)
    return
  end if
  call cmp_ordering(x,order,nn)
  call two_tier_cpu_of_key(order,c,nn,ilevel,bound_key,bound_key_f,two_tier_ls, &
  & two_tier_nfoot,foot_lo,foot_hi,foot_first)
end subroutine cmp_cpumap_lev
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_cpumap_pos(x,c,nn)
  use amr_parameters
  use amr_commons
  use mpi_mod
  implicit none
  integer ::nn
  integer ,dimension(1:nvector)::c
  real(dp),dimension(1:nvector,1:ndim)::x
  !----------------------------------------------------------------------
  ! Owner of the leaf cell that contains the point x (user units). Used as
  ! "if(c==myid) then act on the local cell" (sinks, BH photons, debug).
  ! Single partition: cmp_cpumap, as before. Two-tier: a Hilbert key does
  ! not decide this (the tier depends on the leaf level), so the local tree
  ! is walked down: the finest local grid containing x is the true leaf
  ! grid if and only if it is active (all children of active grids are
  ! held locally), so c==myid exactly for the owner; other cpus get the
  ! owner of the deepest grid they hold (never myid).
  !----------------------------------------------------------------------
  integer::i,idim,nx_loc,nxny,ind,igrid,jind,jcell,jgrid
  real(dp)::scale
  real(dp),dimension(1:3)::xx,skip_loc
  if(.not.two_tier_balance.or.ordering=='bisection')then
    call cmp_cpumap(x,c,nn)
    return
  end if
  nxny=nx*ny
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  do i=1,nn
    xx=0d0
    do idim=1,ndim
      xx(idim)=x(i,idim)/scale+skip_loc(idim)
    end do
    if(xx(1)<0d0)xx(1)=xx(1)+dble(nx)
    if(xx(1)>=dble(nx))xx(1)=xx(1)-dble(nx)
    if(xx(2)<0d0)xx(2)=xx(2)+dble(ny)
    if(xx(2)>=dble(ny))xx(2)=xx(2)-dble(ny)
    if(xx(3)<0d0)xx(3)=xx(3)+dble(nz)
    if(xx(3)>=dble(nz))xx(3)=xx(3)-dble(nz)
    ind=1+int(xx(1))+int(xx(2))*nx+int(xx(3))*nxny
    igrid=son(ind)
    if(igrid==0)then
      c(i)=cpu_map(ind)
      cycle
    end if
    do
      jind=1
      do idim=1,ndim
        if(xx(idim)>xg(igrid,idim))jind=jind+2**(idim-1)
      end do
      jcell=ncoarse+(jind-1)*ngridmax+igrid
      jgrid=son(jcell)
      if(jgrid==0)exit
      igrid=jgrid
    end do
    c(i)=cpu_map(father(igrid))
  end do
end subroutine cmp_cpumap_pos
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine two_tier_cut(key,cost,n,bkey)
  use amr_commons
  use mpi_mod
  implicit none
  integer::n
  real(qdp),dimension(1:n)::key
  real(kind=8),dimension(1:n)::cost
  real(qdp),dimension(0:ndomain)::bkey
  !----------------------------------------------------------------------
  ! Cut the Hilbert curve into ndomain intervals of equal total cost.
  ! Items (key, cost) are distributed over the cpus in any way (no
  ! contiguity assumed, unlike cmp_new_cpu_map): boundary i is found by
  ! bisection in key space on the global cumulative cost C(b) = sum of
  ! cost with key < b, one MPI_ALLREDUCE of ndomain-1 doubles per round
  ! (about 3*(nlevelmax+1) rounds). On return key is sorted.
  !----------------------------------------------------------------------
  integer,allocatable,dimension(:)::perm
  real(kind=8),allocatable,dimension(:)::csum,cloc,cglob
  real(qdp),allocatable,dimension(:)::lo,hi,mid
  real(kind=8)::tot_loc,tot,target
  integer::i,idom,info,iround,j,jlo,jhi
  logical::done

  allocate(perm(1:max(n,1)),csum(0:max(n,1)))
  if(n>1)call quick_sort(key,perm,n)
  if(n==1)perm(1)=1
  csum(0)=0d0
  do i=1,n
    csum(i)=csum(i-1)+cost(perm(i))
  end do
  tot_loc=csum(n)
  tot=tot_loc
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(tot_loc,tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif

  allocate(lo(1:max(ndomain-1,1)),hi(1:max(ndomain-1,1)),mid(1:max(ndomain-1,1)))
  allocate(cloc(1:max(ndomain-1,1)),cglob(1:max(ndomain-1,1)))
  lo=order_all_min; hi=order_all_max
  do iround=1,400
    done=.true.
    do idom=1,ndomain-1
      if(hi(idom)>lo(idom))then
        done=.false.
        mid(idom)=lo(idom)+aint((hi(idom)-lo(idom))/2)
      else
        mid(idom)=lo(idom)
      end if
    end do
    if(done)exit
    do idom=1,ndomain-1
      ! local cost of the items with key < mid
      jlo=0; jhi=n
      do while(jhi>jlo)
        j=(jlo+jhi+1)/2
        if(key(j)<mid(idom))then
          jlo=j
        else
          jhi=j-1
        end if
      end do
      cloc(idom)=csum(jlo)
    end do
    cglob=cloc
#ifndef WITHOUTMPI
    if(ndomain>1)call MPI_ALLREDUCE(cloc,cglob,ndomain-1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
    do idom=1,ndomain-1
      if(hi(idom)>lo(idom))then
        target=tot*dble(idom)/dble(ndomain)
        if(cglob(idom)>=target)then
          hi(idom)=mid(idom)
        else
          lo(idom)=mid(idom)+1
        end if
      end if
    end do
  end do
  if(.not.done.and.myid==1)write(*,*)'two_tier_cut: bisection did not converge'
  bkey(0)=order_all_min
  do idom=1,ndomain-1
    bkey(idom)=lo(idom)
  end do
  bkey(ndomain)=order_all_max
  deallocate(perm,csum,lo,hi,mid,cloc,cglob)
end subroutine two_tier_cut
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_new_cpu_map_two_tier
  use amr_commons
  use mpi_mod
  use pm_commons
  implicit none
  !----------------------------------------------------------------------
  ! Two-tier version of cmp_new_cpu_map: computes two_tier_ls2, bound_key2
  ! (coarse tier) and bound_key2_f (fine tier), then cpu_map2 for every
  ! cell (coarse cells and cells of active grids, spread to virtual cells).
  !----------------------------------------------------------------------
  integer::ilevel,igrid,i,ind,iskip,idim,info,ls,nx_loc,nxny,ix,iy,iz,ngrid,ncache
  integer::n_c,n_f,nmax
  integer(kind=8)::ntot,nfine
  integer(kind=8),dimension(1:MAXLEVEL)::nlev_loc,nlev_tot
  real(dp)::scale,dx,wpart
  real(kind=8),dimension(1:MAXLEVEL)::niter_cost
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(qdp),dimension(1:nvector),save::order
  integer,dimension(1:nvector),save::ind_grid,ind_cell,cc
  real(qdp),allocatable,dimension(:)::key_c,key_f
  real(kind=8),allocatable,dimension(:)::cost_c,cost_f

  nxny=nx*ny
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)

  ! Fine steps per coarse step (the cost_weighting=.true. weights, whatever
  ! cost_weighting says: the fine tier is balanced by time, the coarse by memory)
  niter_cost=1d0
  do ilevel=levelmin+1,nlevelmax
    niter_cost(ilevel)=dble(nsubcycle(ilevel-1))*niter_cost(ilevel-1)
  end do

  ! Global number of grids per level
  nlev_loc=0; nlev_tot=0
  do ilevel=1,nlevelmax
    nlev_loc(ilevel)=active(ilevel)%ngrid
  end do
  nlev_tot=nlev_loc
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nlev_loc,nlev_tot,nlevelmax,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#endif
  ntot=sum(nlev_tot(1:nlevelmax))

  ! Last level of the coarse tier
  if(two_tier_level>0)then
    ls=min(max(two_tier_level,1),nlevelmax)
  else
    ls=max(levelmin,1)
    do while(ls<nlevelmax)
      nfine=sum(nlev_tot(ls+1:nlevelmax))
      if(dble(nfine)<=two_tier_fine_frac*dble(ntot))exit
      ls=ls+1
    end do
  end if
  nfine=0
  if(ls<nlevelmax)nfine=sum(nlev_tot(ls+1:nlevelmax))
  two_tier_ls2=ls

  ! Items: my active grids, key = centre of the father cell
  nmax=int(sum(nlev_loc(1:nlevelmax)))
  allocate(key_c(1:max(nmax,1)),cost_c(1:max(nmax,1)))
  allocate(key_f(1:max(nmax,1)),cost_f(1:max(nmax,1)))
  n_c=0; n_f=0
  wpart=0d0
  if(pic.and.npartmax>0)wpart=dble(ngridmax)/dble(npartmax)
  do ilevel=1,nlevelmax
    ncache=active(ilevel)%ngrid
    do igrid=1,ncache,nvector
      ngrid=MIN(nvector,ncache-igrid+1)
      do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
      end do
      do idim=1,ndim
        do i=1,ngrid
          xx(i,idim)=(xg(ind_grid(i),idim)-skip_loc(idim))*scale
        end do
      end do
      call cmp_ordering(xx,order,ngrid)
      if(ilevel<=ls)then
        ! memory: one oct (1/ngridmax) plus its particles (1/npartmax), in oct units
        do i=1,ngrid
          n_c=n_c+1
          key_c(n_c)=order(i)
          cost_c(n_c)=1d0
          if(pic)cost_c(n_c)=cost_c(n_c)+dble(numbp(ind_grid(i)))*wpart
        end do
      else
        ! time: (cells x 80 + particles) x fine steps per coarse step (cf. cmp_new_cpu_map)
        do i=1,ngrid
          n_f=n_f+1
          key_f(n_f)=order(i)
          cost_f(n_f)=dble(twotondim*80)
          if(pic)cost_f(n_f)=cost_f(n_f)+dble(numbp(ind_grid(i)))
          cost_f(n_f)=cost_f(n_f)*niter_cost(ilevel)
        end do
      end if
    end do
  end do

  ! New boundaries
  call two_tier_cut(key_c,cost_c,n_c,bound_key2)
  if(nfine>0)then
    call two_tier_cut(key_f,cost_f,n_f,bound_key2_f)
  else
    bound_key2_f=bound_key2
  end if
  deallocate(key_c,cost_c,key_f,cost_f)
  ! Footprints of the fine domains (TO ADDED 2026-11-06)
  call two_tier_build_foot(ls,bound_key2_f)
  if(myid==1)write(*,'(" Two-tier balance: Ls=",I3,"  octs <= Ls:",I12,"  octs > Ls:",I12,"  (",F6.2,"%)  footprints:",I9)') &
  & ls,ntot-nfine,nfine,100d0*dble(nfine)/dble(max(ntot,1_8)),two_tier_nfoot2

  ! New cpu map: coarse cells own level-1 grids
  cpu_map2=0
  do iz=0,nz-1
    do iy=0,ny-1
      do ix=0,nx-1
        ind=1+ix+iy*nx+iz*nxny
        xx(1,1)=(dble(ix)+0.5d0-dble(icoarse_min))*scale
#if NDIM>1
        xx(1,2)=(dble(iy)+0.5d0-dble(jcoarse_min))*scale
#endif
#if NDIM>2
        xx(1,3)=(dble(iz)+0.5d0-dble(kcoarse_min))*scale
#endif
        call cmp_ordering(xx,order,1)
        call two_tier_cpu_of_key(order,cc,1,1,bound_key2,bound_key2_f,ls,two_tier_nfoot2,foot_lo2,foot_hi2,foot_first2)
        cpu_map2(ind)=cc(1)
      end do
    end do
  end do
  ! Cells of level ilevel grids own level ilevel+1 grids
  do ilevel=1,nlevelmax
    dx=0.5d0**ilevel
    do ind=1,twotondim
      iz=(ind-1)/4
      iy=(ind-1-4*iz)/2
      ix=(ind-1-2*iy-4*iz)
      xc(ind,1)=(dble(ix)-0.5d0)*dx-dble(icoarse_min)
#if NDIM>1
      xc(ind,2)=(dble(iy)-0.5d0)*dx-dble(jcoarse_min)
#endif
#if NDIM>2
      xc(ind,3)=(dble(iz)-0.5d0)*dx-dble(kcoarse_min)
#endif
    end do
    ncache=active(ilevel)%ngrid
    do igrid=1,ncache,nvector
      ngrid=MIN(nvector,ncache-igrid+1)
      do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
      end do
      do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
          ind_cell(i)=ind_grid(i)+iskip
        end do
        do idim=1,ndim
          do i=1,ngrid
            xx(i,idim)=(xg(ind_grid(i),idim)+xc(ind,idim))*scale
          end do
        end do
        call cmp_ordering(xx,order,ngrid)
        call two_tier_cpu_of_key(order,cc,ngrid,ilevel+1,bound_key2,bound_key2_f,ls, &
        & two_tier_nfoot2,foot_lo2,foot_hi2,foot_first2)
        do i=1,ngrid
          cpu_map2(ind_cell(i))=cc(i)
        end do
      end do
    end do
  end do

  ! Update virtual boundaries for new cpu map
  call make_virtual_coarse_int(cpu_map2(1))
  do ilevel=1,nlevelmax
    call make_virtual_fine_int(cpu_map2(1),ilevel)
  end do

end subroutine cmp_new_cpu_map_two_tier
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine two_tier_build_foot(ls,bkey_f)
  use amr_commons
  use mpi_mod
  implicit none
  integer::ls
  real(qdp),dimension(0:ndomain)::bkey_f
  !----------------------------------------------------------------------
  ! TO ADDED 2026-11-06: footprints of the fine domains of the new
  ! partition (two_tier_nfoot2, foot_lo2, foot_hi2, foot_first2): the key
  ! ranges of the level-ls cells that contain my grids of level > ls,
  ! clipped to the fine domain of each grid, gathered from all cpus,
  ! sorted and merged (see amr_commons).
  !----------------------------------------------------------------------
  integer::ilevel,igrid,i,idim,ngrid,ncache,nloc,n,info,idom,nx_loc,nb
  integer,dimension(1:nvector),save::ind_grid
  real(dp),dimension(1:nvector,1:ndim),save::xx,xa
  real(qdp),dimension(1:nvector),save::order,omin,omax
  real(dp)::dxl,scale
  real(dp),dimension(1:3)::skip_loc
  real(qdp),allocatable,dimension(:)::lo,hi,glo,ghi
  integer,allocatable,dimension(:)::dom,gdom,cnt,disp,cntb,dispb,first
  real(qdp)::oneqdp=1.0

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  dxl=0.5d0**ls

  n=0
  do ilevel=ls+1,nlevelmax
    n=n+active(ilevel)%ngrid
  end do
  allocate(lo(1:max(n,1)),hi(1:max(n,1)),dom(1:max(n,1)))
  nloc=0
  do ilevel=ls+1,nlevelmax
    ncache=active(ilevel)%ngrid
    do igrid=1,ncache,nvector
      ngrid=MIN(nvector,ncache-igrid+1)
      do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
      end do
      do idim=1,ndim
        do i=1,ngrid
          xx(i,idim)=(xg(ind_grid(i),idim)-skip_loc(idim))*scale
          ! centre of the level-ls ancestor cell of the grid
          xa(i,idim)=((dble(int(xg(ind_grid(i),idim)/dxl))+0.5d0)*dxl-skip_loc(idim))*scale
        end do
      end do
      call cmp_ordering(xx,order,ngrid)
      call cmp_minmaxorder(xa,omin,omax,dxl*scale,ngrid)
      do i=1,ngrid
        call two_tier_find_dom(order(i),bkey_f,idom)
        if(min(omax(i),bkey_f(idom))>max(omin(i),bkey_f(idom-1)))then
          nloc=nloc+1
          lo(nloc)=max(omin(i),bkey_f(idom-1))
          hi(nloc)=min(omax(i),bkey_f(idom))
          dom(nloc)=idom
        end if
      end do
    end do
  end do
  call two_tier_merge_foot(nloc,lo,hi,dom)

  ! Gather from all cpus
  allocate(cnt(1:ncpu),disp(1:ncpu),cntb(1:ncpu),dispb(1:ncpu))
  cnt=nloc
#ifndef WITHOUTMPI
  call MPI_ALLGATHER(nloc,1,MPI_INTEGER,cnt,1,MPI_INTEGER,MPI_COMM_WORLD,info)
#endif
  n=sum(cnt)
  disp(1)=0
  do i=2,ncpu
    disp(i)=disp(i-1)+cnt(i-1)
  end do
  nb=storage_size(oneqdp)/8
  cntb=cnt*nb; dispb=disp*nb
  allocate(glo(1:max(n,1)),ghi(1:max(n,1)),gdom(1:max(n,1)))
#ifndef WITHOUTMPI
  call MPI_ALLGATHERV(lo,nloc*nb,MPI_BYTE,glo,cntb,dispb,MPI_BYTE,MPI_COMM_WORLD,info)
  call MPI_ALLGATHERV(hi,nloc*nb,MPI_BYTE,ghi,cntb,dispb,MPI_BYTE,MPI_COMM_WORLD,info)
  call MPI_ALLGATHERV(dom,nloc,MPI_INTEGER,gdom,cnt,disp,MPI_INTEGER,MPI_COMM_WORLD,info)
#else
  if(n>0)then
    glo(1:n)=lo(1:n); ghi(1:n)=hi(1:n); gdom(1:n)=dom(1:n)
  end if
#endif
  call two_tier_merge_foot(n,glo,ghi,gdom)
  allocate(first(0:ndomain))
  first=0
  do i=1,n
    first(gdom(i))=first(gdom(i))+1
  end do
  do idom=1,ndomain
    first(idom)=first(idom)+first(idom-1)
  end do
  call two_tier_set_foot(2,n,glo,ghi,first)
  deallocate(lo,hi,dom,glo,ghi,gdom,cnt,disp,cntb,dispb,first)
end subroutine two_tier_build_foot
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine two_tier_merge_foot(n,lo,hi,dom)
  use amr_commons
  use mpi_mod
  implicit none
  integer::n
  real(qdp),dimension(1:max(n,1))::lo,hi
  integer,dimension(1:max(n,1))::dom
  !----------------------------------------------------------------------
  ! TO ADDED 2026-11-06: sort the intervals by lo and merge overlapping or
  ! touching intervals of the same domain; n is updated. Intervals of
  ! different domains never overlap (each lies inside its domain interval).
  !----------------------------------------------------------------------
  integer::i,k
  integer,allocatable,dimension(:)::perm,tdom
  real(qdp),allocatable,dimension(:)::thi
  if(n<=1)return
  allocate(perm(1:n),thi(1:n),tdom(1:n))
  call quick_sort(lo,perm,n)
  do i=1,n
    thi(i)=hi(perm(i)); tdom(i)=dom(perm(i))
  end do
  hi(1:n)=thi(1:n); dom(1:n)=tdom(1:n)
  k=1
  do i=2,n
    if(dom(i)==dom(k).and.lo(i)<=hi(k))then
      hi(k)=max(hi(k),hi(i))
    else
      k=k+1
      lo(k)=lo(i); hi(k)=hi(i); dom(k)=dom(i)
    end if
  end do
  n=k
  deallocate(perm,thi,tdom)
end subroutine two_tier_merge_foot
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine two_tier_set_foot(which,n,lo,hi,first)
  use amr_commons
  use mpi_mod
  implicit none
  integer::which,n
  real(qdp),dimension(1:max(n,1))::lo,hi
  integer,dimension(0:ndomain)::first
  !----------------------------------------------------------------------
  ! TO ADDED 2026-11-06: store a footprint list as the current (which=1)
  ! or the new (which=2) one.
  !----------------------------------------------------------------------
  if(which==1)then
    if(allocated(foot_lo))deallocate(foot_lo,foot_hi,foot_first)
    allocate(foot_lo(1:max(n,1)),foot_hi(1:max(n,1)),foot_first(0:ndomain))
    two_tier_nfoot=n
    foot_lo(1:max(n,1))=lo(1:max(n,1)); foot_hi(1:max(n,1))=hi(1:max(n,1)); foot_first=first
  else
    if(allocated(foot_lo2))deallocate(foot_lo2,foot_hi2,foot_first2)
    allocate(foot_lo2(1:max(n,1)),foot_hi2(1:max(n,1)),foot_first2(0:ndomain))
    two_tier_nfoot2=n
    foot_lo2(1:max(n,1))=lo(1:max(n,1)); foot_hi2(1:max(n,1))=hi(1:max(n,1)); foot_first2=first
  end if
end subroutine two_tier_set_foot
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine two_tier_bcast_foot
  use amr_commons
  use mpi_mod
  implicit none
  !----------------------------------------------------------------------
  ! TO ADDED 2026-11-06: after a restart, cpu 1 has read the footprint
  ! list from its amr file (read_two_tier_record); give it to everybody
  ! (outside the token-passing read loop) and make it the new list too.
  !----------------------------------------------------------------------
  integer::n,info,nb
  real(qdp),allocatable,dimension(:)::lo,hi
  integer,allocatable,dimension(:)::first
  real(qdp)::oneqdp=1.0
  n=two_tier_nfoot
#ifndef WITHOUTMPI
  call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
#endif
  allocate(lo(1:max(n,1)),hi(1:max(n,1)),first(0:ndomain))
  if(myid==1)then
    lo(1:max(n,1))=foot_lo(1:max(n,1)); hi(1:max(n,1))=foot_hi(1:max(n,1)); first=foot_first
  end if
#ifndef WITHOUTMPI
  nb=storage_size(oneqdp)/8
  call MPI_BCAST(lo,max(n,1)*nb,MPI_BYTE,0,MPI_COMM_WORLD,info)
  call MPI_BCAST(hi,max(n,1)*nb,MPI_BYTE,0,MPI_COMM_WORLD,info)
  call MPI_BCAST(first,ndomain+1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
#endif
  call two_tier_set_foot(1,n,lo,hi,first)
  call two_tier_set_foot(2,n,lo,hi,first)
  deallocate(lo,hi,first)
  if(myid==1.and.n>0)write(*,*)'Two-tier fine footprints read from the restart files:',n
end subroutine two_tier_bcast_foot
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine two_tier_cpu_hull(icpu,key_lo,key_hi)
  use amr_commons
  use mpi_mod
  implicit none
  integer::icpu
  real(qdp)::key_lo,key_hi
  !----------------------------------------------------------------------
  ! Hull of all key intervals (both tiers, all overload sub-domains) of
  ! cpu icpu, for the DOMAIN table of the info file.
  !----------------------------------------------------------------------
  integer::isub,idom,jdom
  key_lo=order_all_max; key_hi=order_all_min
  do isub=1,overload
    idom=icpu+(isub-1)*ncpu
    key_lo=min(key_lo,bound_key(idom-1)); key_hi=max(key_hi,bound_key(idom))
    if(two_tier_balance)then
      jdom=two_tier_fine_dom(icpu,isub)
      if(foot_first(jdom)>foot_first(jdom-1))then ! TO MODIFIED 2026-11-06: footprints, not the whole interval
        key_lo=min(key_lo,foot_lo(foot_first(jdom-1)+1)); key_hi=max(key_hi,foot_hi(foot_first(jdom)))
      end if
    end if
  end do
end subroutine two_tier_cpu_hull
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine two_tier_report
  use amr_commons
  use mpi_mod
  implicit none
  !----------------------------------------------------------------------
  ! Grid statistics after a load balance: active and ghost (virtual) octs
  ! per level, total and maximum over cpus (compare with
  ! tools/lb_two_tier_estimate.py).
  !----------------------------------------------------------------------
  integer::ilevel,icpu,info,nfmax
  integer(kind=8),dimension(1:MAXLEVEL,1:2)::nloc,nsum,nmax
  integer(kind=8)::mloc,msum,mmax
  nloc=0
  do ilevel=1,nlevelmax
    nloc(ilevel,1)=active(ilevel)%ngrid
    do icpu=1,ncpu
      nloc(ilevel,2)=nloc(ilevel,2)+reception(icpu,ilevel)%ngrid
    end do
  end do
  mloc=sum(nloc(1:nlevelmax,1))+sum(nloc(1:nlevelmax,2))
  nsum=nloc; nmax=nloc; msum=mloc; mmax=mloc
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nloc,nsum,2*MAXLEVEL,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(nloc,nmax,2*MAXLEVEL,MPI_INTEGER8,MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mloc,msum,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mloc,mmax,1,MPI_INTEGER8,MPI_MAX,MPI_COMM_WORLD,info)
#endif
  if(myid==1)then
    write(*,'(" Two-tier balance (Ls=",I3,"): level   active_tot   active_max    ghost_tot    ghost_max  ghost/active")')two_tier_ls
    do ilevel=1,nlevelmax
      if(nsum(ilevel,1)==0)cycle
      write(*,'(31X,I3,4I13,F12.3)')ilevel,nsum(ilevel,1),nmax(ilevel,1),nsum(ilevel,2),nmax(ilevel,2), &
      & dble(nsum(ilevel,2))/dble(nsum(ilevel,1))
    end do
    write(*,'("   octs per cpu (active+ghost): max",I12," (",F6.2,"% of ngridmax) avg",I12)') &
    & mmax,100d0*dble(mmax)/dble(ngridmax),msum/ncpu
    nfmax=0
    do icpu=1,ndomain
      nfmax=max(nfmax,foot_first(icpu)-foot_first(icpu-1))
    end do
    write(*,'("   fine footprints:",I10," intervals, max per domain",I8)')two_tier_nfoot,nfmax
  end if
end subroutine two_tier_report
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine check_particle_tree(tag)
  use amr_commons
  use mpi_mod
  use pm_commons
  implicit none
  character(len=*)::tag
  !----------------------------------------------------------------------
  ! TO ADDED 2026-11-05 (two-tier verification, -DTWO_TIER_CHECK): called
  ! after the particle sort that ends with merge_tree_fine down to levelmin
  ! (end of load_balance, end of create_sink). At that point the lists of
  ! the active levelmin grids are the authoritative ones (merge_tree_fine
  ! leaves stale copies in the finer grids, so finer levels must not be
  ! counted). Checks: every live particle (npart) is in exactly one
  ! levelmin list (active or virtual), and inside that grid. Collective; aborts.
  !----------------------------------------------------------------------
  integer::i,igrid,ipart,jpart,next_part,idim,nbad,info,nx_loc,nprint,icpu,ncache
  integer(kind=8)::npl,npt,nlist,nlist_tot,nbad_tot
  real(dp)::dx,scale,xxx
  real(dp),dimension(1:3)::skip_loc
  logical::ok

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  ! Lists of levelmin grids, active and virtual (after merge_tree_fine a particle whose finest
  ! grid is owned by me can sit in a virtual levelmin ancestor; the next step's kill_tree_fine /
  ! virtual_tree_fine sort it out, as in the single partition).
  nbad=0; nprint=0; nlist=0
  dx=0.5d0**levelmin
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
        nlist=nlist+1
        ok=.true.
        do idim=1,ndim
          xxx=(xp(ipart,idim)/scale+skip_loc(idim)-(xg(igrid,idim)-dx))/dx
          if(xxx<0d0.or.xxx>=2d0)ok=.false.
        end do
        if(.not.ok)then
          nbad=nbad+1
          if(nprint<5)then
            nprint=nprint+1
            write(*,'(" TWO_TIER CHECK particle outside its grid: cpu",I6," list",I6," grid",I9," idp",I12," xp",3(1X,ES14.7)," xg",3(1X,ES14.7))') &
            & myid,icpu,igrid,idp(ipart),xp(ipart,1:ndim),xg(igrid,1:ndim)
          end if
        end if
        ipart=next_part
      end do
    end do
  end do
  npl=npart
  nbad_tot=nbad; npt=npl; nlist_tot=nlist
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(npl,npt,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(nlist,nlist_tot,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
  npl=nbad
  call MPI_ALLREDUCE(npl,nbad_tot,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#endif
  if(myid==1)write(*,'(" TWO_TIER CHECK particles (",A,") step",I7,":",I13," live,",I13," in levelmin lists,",I9," outside their grid")') &
  & trim(tag),nstep_coarse,npt,nlist_tot,nbad_tot
  if(nbad_tot>0.or.nlist_tot/=npt)then
    if(myid==1)write(*,*)'TWO_TIER CHECK FAILED (particles)'
    flush(6)
#ifndef WITHOUTMPI
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
    stop
#endif
  end if
end subroutine check_particle_tree
