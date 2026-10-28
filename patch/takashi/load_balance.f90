!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine load_balance
  use amr_commons
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
        call make_virtual_fine_dp_nv(uold,nvar+3,1,nvar+3,ilevel)
#else
        call make_virtual_fine_dp_nv(uold,nvar,1,nvar,ilevel)
#endif
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
        call make_virtual_fine_dp_nv(rtuold,nrtvar,1,nrtvar,ilevel)
        if(simple_boundary)then
           call rt_make_boundary_hydro(ilevel)
        end if
     endif
#endif
     if(poisson)then
        call make_virtual_fine_dp(phi(1),ilevel)
        call make_virtual_fine_dp_nv(f,ndim,1,ndim,ilevel)
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

  balance=.false.

  if(verbose)then
     write(*,*)'Output mesh structure'
     do ilevel=1,nlevelmax
        if(numbtot(1,ilevel)>0)write(*,999)ilevel,numbtot(1:4,ilevel)
     end do
  end if

  !--------------------------------------------------------------
  ! Report the size of the virtual boundaries (ghost octs), which
  ! grows when the decomposition differs between levels.
  !--------------------------------------------------------------
  do ilevel=1,nlevelmax
     ncache=0
     do icpu=1,ncpu
        if(icpu/=myid)ncache=ncache+numbl(icpu,ilevel)
     end do
     comm_buffin(ilevel,1)=ncache
     comm_buffin(ilevel,2)=ncache
     comm_buffin(ilevel,3)=ncache
  end do
  call MPI_ALLREDUCE(comm_buffin(1,1),comm_buffout(1,1),nlevelmax,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(comm_buffin(1,2),comm_buffout(1,2),nlevelmax,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(comm_buffin(1,3),comm_buffout(1,3),nlevelmax,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
  if(myid==1)then
     do ilevel=levelmin,nlevelmax
        if(numbtot(1,ilevel)>0)write(*,998)ilevel,comm_buffout(ilevel,1),comm_buffout(ilevel,2),comm_buffout(ilevel,3), &
             & dble(comm_buffout(ilevel,1))/dble(max(numbtot(1,ilevel),1_i8b))
     end do
  end if

998 format(' Level ',I2,' has ',I10,' ghost octs (min=',I8,', max=',I8,', ghost/active=',F6.3,')')
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
  !---------------------------------------------------------------------
  ! This routine computes the new cpu map using the chosen ordering
  ! to balance load across cpus.
  ! If balance_per_level is true, an independent set of domain boundaries
  ! bound_key2(:,ilevel) is computed for every grid level (ilevel is the
  ! level of the octs; bound_key2(:,1) also serves the coarse cells).
  ! Otherwise a single set is computed from the leaf cells (original
  ! RAMSES algorithm) and copied to all levels.
  !---------------------------------------------------------------------
  integer::igrid,ncache,ngrid,jlevel
  integer::ilevel,i,ind,idim
  integer::nx_loc
  integer::idom
  integer::nxny,ix,iy,iz,iskip
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::dx,scale
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(qdp),dimension(1:nvector),save::order_max
  real(dp),dimension(1:1,1:ndim),save :: xx_tmp
  integer,dimension(1:1),save :: c_tmp

  ! Local constants
  nxny=nx*ny
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)

  if(verbose) print *,"Entering cmp_new_cpu_map"

  if(ordering/='bisection') then
     if(balance_per_level)then
        call cmp_new_bound_key_per_level
     else
        call cmp_new_bound_key_single
     end if
  else
     ! update the bisection
     call build_bisection(update=.true.)
  end if

  !----------------------------------------
  ! Compute new cpu map
  !----------------------------------------
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
     cpu_map2(ind)=ncpu ! default value

     if(ordering/='bisection') then
        call cmp_ordering(xx,order_max,1)
        cpu_map2(ind)=ncpu ! default value
        do idom=1,ndomain
           if( order_max(1).ge.bound_key2(idom-1,1).and. &
                & order_max(1).lt.bound_key2(idom,1))then
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
     ! Level of the octs that the cells of level ilevel can host
     jlevel=min(ilevel+1,nlevelmax)
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
                    if( order_max(i).ge.bound_key2(idom-1,jlevel).and. &
                         & order_max(i).lt.bound_key2(idom,jlevel))then
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
subroutine cmp_new_bound_key_single
  use amr_commons
  use pm_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  !---------------------------------------------------------------------
  ! Original RAMSES load balancing: the cost of all leaf cells (weighted
  ! by the number of sub-cycles of their level) is distributed along the
  ! ordering curve. The result is stored in bound_key2(:,1) and copied
  ! to all levels. This requires the current decomposition to be the
  ! same on all levels.
  !---------------------------------------------------------------------
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
  logical::uniform
#ifdef QUADHILBERT
  real(kind=8),dimension(0:ndomain)::bigdbl,bigtmp
#endif

  call check_bound_key_uniform(uniform)
  if(.not.uniform)then
     if(myid==1)then
        write(*,*)'cmp_new_bound_key_single: the current decomposition differs between levels.'
        write(*,*)'Use balance_per_level=.true.'
     end if
     call clean_stop
  end if

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
  else
     niter_cost(levelmin:nlevelmax)=1
  endif

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
        call cmp_dommap(xx,dom,ncell_loc,1)
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
                 call cmp_dommap(xx,dom,ncell_loc,min(ilevel+1,nlevelmax))
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

  !------------------------------------------------
  ! Sort ordering key and store new index in flag2
  !------------------------------------------------
  if (ncell>0) call quick_sort(hilbert_key(1),flag2(1),ncell)

  !-----------------------------
  ! Balance cost across cpus
  !-----------------------------
  cost_loc = 0 ! Compute local and global cost
  do isub=1,overload
     cost_loc(myid+(isub-1)*ncpu) = dble(npart_sub(isub))
  end do
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(cost_loc,cost_old,ndomain,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
  cost_old=cost_loc
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
           cell_cost=dble(flag1(flag2(ind_long+ncell_loc)))
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
  bound_key2(:,1)    = real(bigtmp,kind=qdp)
#else
  call MPI_ALLREDUCE(bound_key_loc,bound_key2(0,1),ndomain+1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#else
  bound_key2(:,1)=bound_key_loc
#endif
  bound_key2(0,1)      =order_all_min
  bound_key2(ndomain,1)=order_all_max

  ! Same decomposition on all levels
  do ilevel=2,nlevelmax
     bound_key2(0:ndomain,ilevel)=bound_key2(0:ndomain,1)
  end do

end subroutine cmp_new_bound_key_single
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_new_bound_key_per_level
  use amr_commons
  use pm_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  !---------------------------------------------------------------------
  ! Per-level load balancing. For every grid level ilevel the cost of
  ! the local octs (twotondim*80 + twotondim*numbp per oct, weighted by
  ! the number of sub-cycles) is distributed along the ordering curve,
  ! independently of the other levels, giving bound_key2(:,ilevel).
  ! Boundaries are placed at the maximum key of the father cell of an
  ! oct, so that the ownership of an oct (decided from the key of its
  ! father cell centre) is never ambiguous.
  ! With balance_tolerance>0 the boundaries of a level are moved, within
  ! the allowed imbalance, towards those of the next finer level (nested
  ! decompositions, fewer ghost octs); levels below balance_level_min
  ! copy the decomposition of that level.
  ! Work space: hilbert_key (keys), flag1 (cost), flag2 (sort index).
  !---------------------------------------------------------------------
  integer::igrid,ncache,ngrid,ngrid_loc,ngrid_lev
  integer::ilevel,i,ioff
  integer::nx_loc
  integer::isub,idom
  integer::ind_long
  integer,dimension(1:nvector),save::ind_grid
  integer,dimension(1:nvector),save::dom
  real(dp)::dx,scale
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(qdp),dimension(1:nvector),save::order_min,order_max
  real(kind=8)::incost_tot,local_cost,cell_cost
  real(kind=8),dimension(0:ndomain)::incost_new,incost_old
  real(kind=8),dimension(1:ndomain,1:nlevelmax)::cost_loc,cost_old
  real(kind=8),dimension(1:overload,1:nlevelmax)::cost_sub
  integer,dimension(1:overload,1:nlevelmax)::ngrid_sub
  integer,dimension(1:nlevelmax)::ngrid_off
  integer(kind=8)::wflag
  real(qdp),dimension(0:ndomain)::bkey_loc
  real(kind=8),dimension(0:ndomain)::bex,bex_tot,klo,khi,klo_tot,khi_tot
  real(kind=8)::dom_cost,wtol,tol,bnew
  integer,dimension(1:MAXLEVEL),save::niter_cost
  logical,dimension(1:nlevelmax)::has_cost
  integer::j,lmax_pop,lmin_bal,lref

  ! Local constants
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)

  ! Compute time step related cost
  if(cost_weighting)then
     niter_cost(1:levelmin)=1
     if (nlevelmax - levelmin - 1 > 31) write(*,*) 'Warning load_balance: niter_cost may need to become a kind=8 integer'
     do ilevel=levelmin+1,nlevelmax
        niter_cost(ilevel)=nsubcycle(ilevel-1)*niter_cost(ilevel-1)
     end do
  else
     niter_cost(1:nlevelmax)=1
  endif

  !----------------------------------------------------------------
  ! Gather, for every level, the key (maximum key of the father
  ! cell) and the cost of each local oct. Octs of level ilevel are
  ! stored in hilbert_key/flag1 from ngrid_off(ilevel)+1 on.
  !----------------------------------------------------------------
  cost_sub=0d0
  ngrid_sub=0
  ngrid_loc=0
  do ilevel=1,nlevelmax
     ngrid_off(ilevel)=ngrid_loc
     ! Size of the father cells
     dx=0.5d0**(ilevel-1)*scale
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        ! Oct centre = father cell centre, in user units
        do i=1,ngrid
           xx(i,1)=(xg(ind_grid(i),1)-dble(icoarse_min))*scale
#if NDIM>1
           xx(i,2)=(xg(ind_grid(i),2)-dble(jcoarse_min))*scale
#endif
#if NDIM>2
           xx(i,3)=(xg(ind_grid(i),3)-dble(kcoarse_min))*scale
#endif
        end do
        call cmp_minmaxorder(xx,order_min,order_max,dx,ngrid)
        call cmp_dommap(xx,dom,ngrid,ilevel)
        do i=1,ngrid
           ngrid_loc=ngrid_loc+1
           isub=(dom(i)-1)/ncpu+1
           ngrid_sub(isub,ilevel)=ngrid_sub(isub,ilevel)+1
           wflag=8*10*twotondim
           if(pic)then
              if (tracer .and. ilevel >= tracer_first_balance_levelmin) then
                 wflag=wflag+twotondim*(numbp(ind_grid(i))+tracer_first_balance_part_per_cell)
              else
                 wflag=wflag+twotondim*numbp(ind_grid(i))
              endif
           end if
           wflag=wflag*niter_cost(ilevel)
           if (wflag > 2147483647_8) then
              write(*,*) ' wrong type for flag1 --> change to integer kind=8: ',wflag
              stop
           endif
           flag1(ngrid_loc)=int(wflag)
           cost_sub(isub,ilevel)=cost_sub(isub,ilevel)+dble(wflag)
           hilbert_key(ngrid_loc)=order_max(i)
        end do
     end do
  end do

  !------------------------------------------------
  ! Sort ordering keys level by level (index in flag2)
  !------------------------------------------------
  do ilevel=1,nlevelmax
     ngrid_lev=ngrid_loc_of_level(ilevel)
     if(ngrid_lev>0) call quick_sort(hilbert_key(ngrid_off(ilevel)+1),flag2(ngrid_off(ilevel)+1),ngrid_lev)
  end do

  !-----------------------------
  ! Balance cost across cpus
  !-----------------------------
  cost_loc=0d0
  do ilevel=1,nlevelmax
     do isub=1,overload
        cost_loc(myid+(isub-1)*ncpu,ilevel)=cost_sub(isub,ilevel)
     end do
  end do
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(cost_loc,cost_old,ndomain*nlevelmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
  cost_old=cost_loc
#endif

  ! Finest populated level, lowest level balanced on its own, tolerance
  has_cost=.false.
  lmax_pop=0
  do ilevel=1,nlevelmax
     if(sum(cost_old(1:ndomain,ilevel))>0d0)then
        has_cost(ilevel)=.true.
        lmax_pop=ilevel
     end if
  end do
  lmin_bal=min(max(balance_level_min,levelmin),nlevelmax)
  tol=min(max(balance_tolerance,0d0),0.9d0)   ! max relative deviation of a domain cost

  !----------------------------------------------------------------
  ! Levels are processed from the finest populated one downwards.
  ! The finest level is balanced exactly. For a coarser level each
  ! boundary may be placed anywhere in the window of octs whose
  ! cumulative cost is within +-tol/2*(cost/ndomain) of the exact
  ! quantile (so a domain cost deviates by at most +-tol from the
  ! mean); inside that window the position closest to the
  ! boundary of the reference (next finer balanced) level is taken,
  ! which nests the decompositions and reduces the ghost octs.
  ! Levels below lmin_bal copy the decomposition of lmin_bal, empty
  ! levels copy the finer level, levels above the finest populated
  ! level copy it.
  !----------------------------------------------------------------
  bound_key2=0.0d0
  lref=0
  do ilevel=lmax_pop,1,-1
     if(.not.has_cost(ilevel))then
        bound_key2(0:ndomain,ilevel)=bound_key2(0:ndomain,ilevel+1)
        cycle
     end if
     if(ilevel<lmin_bal.and.lref>0)then
        bound_key2(0:ndomain,ilevel)=bound_key2(0:ndomain,lref)
        cycle
     end if

     incost_tot = 0D0
     incost_old(0) = 0D0
     do idom = 1,ndomain
        incost_tot = incost_tot + cost_old(idom,ilevel)
        incost_old(idom) = incost_tot
     end do
     dom_cost=incost_tot/dble(ndomain)
     do idom = 0,ndomain
        incost_new(idom) = dble(idom)*dom_cost ! Exact load balancing
     end do
     wtol=0d0
     if(lref>0)wtol=0.5d0*tol*dom_cost   ! half per boundary: windows never overlap

     ! Scan the local octs (sorted by key): exact boundaries and windows
     bkey_loc=0.0d0
     klo=huge(1d0); khi=-1d0
     ioff=ngrid_off(ilevel)
     do isub=1,overload
        if(ngrid_sub(isub,ilevel)>0)then
           ! First cpu on local domain
           idom=0
           do while(incost_new(idom)<incost_old(myid-1+(isub-1)*ncpu))
              idom=idom+1
              if (idom > ndomain) exit
           end do
           ! Compute Hilbert key at boundaries
           i=idom
           local_cost=incost_old(myid-1+(isub-1)*ncpu)
           do ind_long=1,ngrid_sub(isub,ilevel)
              cell_cost=dble(flag1(ioff+flag2(ioff+ind_long)))
              ! This oct is a candidate first oct of domain j+1 if the
              ! cost accumulated before it lies in the window of boundary j
              if(wtol>0d0)then
                 j=nint(local_cost/dom_cost)
                 if(j>=1.and.j<=ndomain-1)then
                    if(abs(local_cost-incost_new(j))<=wtol)then
                       klo(j)=min(klo(j),real(hilbert_key(ioff+ind_long),kind=8))
                       khi(j)=max(khi(j),real(hilbert_key(ioff+ind_long),kind=8))
                    end if
                 end if
              end if
              local_cost=local_cost+cell_cost
              if (i > ndomain) exit
              if(incost_new(i)<local_cost)then
                 bkey_loc(i)=hilbert_key(ioff+ind_long)
                 i=i+1
              endif
           end do
        end if
        ioff=ioff+ngrid_sub(isub,ilevel)
     end do

     bex=real(bkey_loc,kind=8)
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(bex,bex_tot,ndomain+1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     if(wtol>0d0)then
        call MPI_ALLREDUCE(klo,klo_tot,ndomain+1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
        call MPI_ALLREDUCE(khi,khi_tot,ndomain+1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
     end if
#else
     bex_tot=bex; klo_tot=klo; khi_tot=khi
#endif

     do idom=1,ndomain-1
        bnew=bex_tot(idom)
        if(wtol>0d0)then
           if(klo_tot(idom)<=khi_tot(idom))then
              if(real(bound_key2(idom,lref),kind=8)<klo_tot(idom))then
                 bnew=klo_tot(idom)
              else if(real(bound_key2(idom,lref),kind=8)>khi_tot(idom))then
                 bnew=khi_tot(idom)
              else
                 bound_key2(idom,ilevel)=bound_key2(idom,lref) ! nested boundary
                 cycle
              end if
           end if
        end if
        bound_key2(idom,ilevel)=real(bnew,kind=qdp)
     end do
     bound_key2(0,ilevel)      =order_all_min
     bound_key2(ndomain,ilevel)=order_all_max
     lref=ilevel
  end do
  do ilevel=lmax_pop+1,nlevelmax
     bound_key2(0:ndomain,ilevel)=bound_key2(0:ndomain,lmax_pop)
  end do

  if(myid==1.and.verbose)then
     do ilevel=1,nlevelmax
        write(*,'(" Level ",I2," new bound_key:",100(1X,E12.5))')ilevel,(real(bound_key2(idom,ilevel),kind=8),idom=0,min(ndomain,20))
     end do
  end if

contains

  integer function ngrid_loc_of_level(il)
    integer,intent(in)::il
    if(il<nlevelmax)then
       ngrid_loc_of_level=ngrid_off(il+1)-ngrid_off(il)
    else
       ngrid_loc_of_level=ngrid_loc-ngrid_off(il)
    end if
  end function ngrid_loc_of_level

end subroutine cmp_new_bound_key_per_level
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine check_bound_key_uniform(uniform)
  use amr_commons
  implicit none
  logical,intent(out)::uniform
  !-----------------------------------------------------------------
  ! True if the domain decomposition is the same on all levels.
  !-----------------------------------------------------------------
  integer::ilevel,idom
  uniform=.true.
  if(ordering=='bisection')return
  do ilevel=2,nlevelmax
     do idom=0,ndomain
        if(bound_key(idom,ilevel)/=bound_key(idom,1))then
           uniform=.false.
           return
        end if
     end do
  end do
end subroutine check_bound_key_uniform
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_cpumap(x,c,nn,ilevel)
  use amr_parameters
  use amr_commons
  use bisection
  implicit none
  integer, intent(in) ::nn,ilevel
  real(dp),dimension(1:nvector,1:ndim), intent(in)::x
  integer ,dimension(1:nvector), intent(out)::c
  !-----------------------------------------------------------------
  ! Cpu owning the oct of level ilevel whose centre (= centre of its
  ! father cell) is at position x (user units).
  !-----------------------------------------------------------------
  integer::i,idom,jlevel
  real(qdp),dimension(1:nvector),save::order

  jlevel=min(max(ilevel,1),nlevelmax)
  if(ordering /= 'bisection') then
     call cmp_ordering(x,order,nn)
     do i=1,nn
        c(i)=ndomain ! default value
        do idom=1,ndomain
           if(    order(i).ge.bound_key(idom-1,jlevel).and. &
                & order(i).lt.bound_key(idom  ,jlevel))then
              c(i)=idom
           endif
        end do
     end do
     do i=1,nn
        c(i)=MOD(c(i)-1,ncpu)+1
     end do
  else
     call cmp_bisection_cpumap(x,c,nn)
  end if

end subroutine cmp_cpumap
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_dommap(x,c,nn,ilevel)
  use amr_parameters
  use amr_commons
  use bisection
  implicit none
  integer ::nn,ilevel
  integer ,dimension(1:nvector)::c
  real(dp),dimension(1:nvector,1:ndim)::x
  !-----------------------------------------------------------------
  ! Domain index (1:ndomain) of the oct of level ilevel centred at x.
  !-----------------------------------------------------------------
  integer::i,idom,jlevel
  real(qdp),dimension(1:nvector),save::order

  jlevel=min(max(ilevel,1),nlevelmax)
  call cmp_ordering(x,order,nn)
  do i=1,nn
     c(i)=ndomain ! default value
     do idom=1,ndomain
        if(    order(i).ge.bound_key(idom-1,jlevel).and. &
             & order(i).lt.bound_key(idom  ,jlevel))then
           c(i)=idom
        endif
     end do
  end do

end subroutine cmp_dommap
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_cpumap_leaf(x,c,nn)
  use amr_parameters
  use amr_commons
  implicit none
  integer, intent(in) ::nn
  real(dp),dimension(1:nvector,1:ndim), intent(in)::x
  integer ,dimension(1:nvector), intent(out)::c
  !-----------------------------------------------------------------
  ! Cpu "owning" the leaf cell that contains position x (user units),
  ! i.e. cpu_map of that leaf cell. The local tree is walked down as
  ! far as it exists: c(i)==myid if and only if the leaf cell is a
  ! myid cell (a cell with cpu_map==myid whose son grid exists
  ! globally is always present locally), so exactly one cpu answers
  ! myid for a given position. This replaces cmp_cpumap(x,c,nn) for
  ! point queries when the decomposition depends on the level.
  !-----------------------------------------------------------------
  integer::i,idim,ix,iy,iz,ind,icell,igrid,nx_loc,nxny
  real(dp)::scale
  real(dp),dimension(1:3)::skip_loc,xx

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  nxny=nx*ny
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)

  do i=1,nn
     xx=0d0
     do idim=1,ndim
        xx(idim)=x(i,idim)/scale+skip_loc(idim)
     end do
     ix=0; iy=0; iz=0
     if(ndim>0)ix=min(max(int(xx(1)),0),nx-1)
     if(ndim>1)iy=min(max(int(xx(2)),0),ny-1)
     if(ndim>2)iz=min(max(int(xx(3)),0),nz-1)
     icell=1+ix+iy*nx+iz*nxny
     c(i)=cpu_map(icell)
     igrid=son(icell)
     do while(igrid>0)
        ind=1
        if(xx(1)>xg(igrid,1))ind=ind+1
#if NDIM>1
        if(xx(2)>xg(igrid,2))ind=ind+2
#endif
#if NDIM>2
        if(xx(3)>xg(igrid,3))ind=ind+4
#endif
        icell=ncoarse+(ind-1)*ngridmax+igrid
        c(i)=cpu_map(icell)
        igrid=son(icell)
     end do
  end do

end subroutine cmp_cpumap_leaf
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
  use hydro_commons
#ifdef RT
  use rt_hydro_commons
#endif
  implicit none

  integer::ncache,ngrid2,igridmax,i,igrid,ibound,ilevel
  integer::iskip1,iskip2,igrid1,igrid2,ind1,icell1,icell2
  integer::ind,idim,ivar,istart

  if(verbose)write(*,*)'Defragmenting main memory...'

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
              cpu_map2(igrid)=ngrid2+i
              igridmax=max(igridmax,igrid)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
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
