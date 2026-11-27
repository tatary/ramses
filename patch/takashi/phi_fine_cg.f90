!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine phi_fine_cg(ilevel,icount)
  use amr_commons
!$ use omp_lib
  use pm_commons
  use poisson_commons
  use mpi_mod
  use constants, only: twopi
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  integer::ilevel,icount
  !=========================================================
  ! Iterative Poisson solver with Conjugate Gradient method
  ! to solve A x = b
  ! r  : stored in f(i,1)
  ! p  : stored in f(i,2)
  ! A p: stored in f(i,3)
  ! x  : stored in phi(i)
  ! b  : stored in rho(i)
  !=========================================================
  integer::i,ind,iter,iskip,itermax,nx_loc
  integer::idx,ngrid_act,tid
  logical::warm,thr,prec
  real(kind=8)::rz,rz_loc,rz_old,red2(2),red2_all(2)
  real(dp)::wall_time,tw0,tw1,t_comm,t_red,t_setup
  real(kind=8)::r2_loc,psum
  real(kind=8),dimension(0:255)::part   ! per-thread partial sums, added in thread order (reproducible)
  real(dp)::error,error_ini
  real(dp)::dx2,fourpi,scale,oneoversix,fact,fact2
  real(dp)::r2_old=0,alpha_cg,beta_cg
  real(kind=8)::r2,pAp,rhs_norm
#ifndef WITHOUTMPI
  real(kind=8) :: rhs_norm_all, pAp_all, r2_all
#endif

  if(gravity_type>0)return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Time breakdown (cg_stat): ghost exchange and MPI_ALLREDUCE of the iterations; the rest of the
  ! solve is local work (plus the set-up). The ALLREDUCE time includes the wait for the slowest rank.
  tw0=wall_time(); t_comm=0d0; t_red=0d0

  ! Set constants
  dx2=(0.5D0**ilevel)**2
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  fourpi=2*twopi*scale
  if(cosmo)fourpi=1.5D0*omega_m*aexp*scale
  oneoversix=1.0D0/dble(twondim)
  fact=oneoversix*fourpi*dx2
  fact2 = fact*fact

  !===============================
  ! Compute initial phi
  !===============================
   warm=cg_warm_start.and.ilevel>levelmin
   if(warm)warm=cg_warm_ok(ilevel)
   if(warm)then
      call make_warm_phi(ilevel,icount)          ! cg_warm_start: previous solution (phi_old)
   else if(ilevel>levelmin)then
      call make_initial_phi(ilevel,icount)              ! Interpolate phi down
   else
      call make_multipole_phi(ilevel)            ! Fill up with simple initial guess
   endif
   call make_virtual_fine_dp(phi(1),ilevel)      ! Update boundaries
   call make_boundary_phi(ilevel)                ! Update physical boundaries

  !===============================
  ! Compute right-hand side norm
  !===============================
  rhs_norm=0
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        idx=active(ilevel)%igrid(i)+iskip
        rhs_norm=rhs_norm+fact2*(rho(idx)-rho_tot)*(rho(idx)-rho_tot)
     end do
  end do
  ! Compute global norms
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(rhs_norm,rhs_norm_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       & MPI_COMM_WORLD,info)
  rhs_norm=rhs_norm_all
#endif
  rhs_norm=DSQRT(rhs_norm/dble(twotondim*numbtot(1,ilevel)))

  !==============================================
  ! Compute r = b - Ax and store it into f(i,1)
  ! Also set p = r and store it into f(i,2)
  !==============================================
  call cmp_residual_cg(ilevel,icount)

  !====================================
  ! Main iteration loop
  !====================================
  ! OpenMP (omp_poisson, TODO.md 2d-4): the loops of one iteration are threaded as three parallel regions
  ! (p update; Ap and p.Ap; x and r updates with the next r.r). The sums keep the serial order (ind outer,
  ! grid inner), so one thread is bit-identical to the serial code; with several threads the reductions
  ! change the rounding. The partial sums of the threads are added in thread order (not with a
  ! reduction clause, whose combination order is unspecified), so a run is reproducible for a given
  ! number of threads. The regions are short (a fraction of a ms), so the workers must not go to
  ! sleep between them: run with KMP_BLOCKTIME=1 (Intel) / GOMP_SPINCOUNT (GNU) rather than 0.
  ! Two forms of the (ind,i) loops, chosen per solve (thr):
  ! - not threaded (default: omp_cg=.false., or few grids): plain loops without any OpenMP construct.
  !   With ifx a collapse(2) loop is ~28% slower than the plain one even on one thread (it recovers
  !   (ind,i) from a flat index), and at ~2400 grids per rank the threaded iteration is slower than
  !   the plain one in either form (cluster 2026-11-17/18, level 15: plain 0.0237 s, threaded
  !   collapse 0.0281 s, threaded without collapse 0.0411 s).
  ! - threaded (omp_poisson .and. omp_cg): collapse(2), the faster of the two threaded forms.
  ngrid_act=active(ilevel)%ngrid
  thr=.false.
!$ thr=omp_poisson .and. omp_cg .and. ngrid_act>=omp_mg_mingrid
  ! cg_precond=1: preconditioned CG, z = M^-1 r with M one symmetric Gauss-Seidel (SSOR) sweep over
  ! the cells of this rank only (block Jacobi over the ranks: no communication). Plain loops only.
  prec=cg_precond>0
  if(prec)thr=.false.
  rz=0.0d0; rz_loc=0.0d0; rz_old=0.0d0
  if(thr)then
     part=0.0d0
!$omp parallel private(ind,i,idx,tid,psum)
     tid=0
!$   tid=omp_get_thread_num()
     psum=0.0d0
!$omp do collapse(2) schedule(static)
     do ind=1,twotondim
        do i=1,ngrid_act
           idx=active(ilevel)%igrid(i)+ncoarse+(ind-1)*ngridmax
           psum=psum+f(idx,1)*f(idx,1)
        end do
     end do
!$omp end do nowait
     part(tid)=psum
!$omp end parallel
     r2_loc=0.0d0
     do tid=0,255
        r2_loc=r2_loc+part(tid)
     end do
  else
     r2_loc=0.0d0
     do ind=1,twotondim
        do i=1,ngrid_act
           idx=active(ilevel)%igrid(i)+ncoarse+(ind-1)*ngridmax
           r2_loc=r2_loc+f(idx,1)*f(idx,1)
        end do
     end do
  end if

  if(prec)then
     if(.not.allocated(cg_zw))then
        allocate(cg_zw(0:ncoarse+twotondim*ngridmax),cg_yw(0:ncoarse+twotondim*ngridmax),cg_own(1:ngridmax))
        cg_zw=0.0d0; cg_yw=0.0d0; cg_own=0; cg_stamp=0
     end if
     cg_stamp=cg_stamp+1
     do i=1,ngrid_act
        cg_own(active(ilevel)%igrid(i))=cg_stamp
     end do
     call cg_precond_table(ilevel)
     call cg_precond_sgs(ilevel)
     do ind=1,twotondim
        do i=1,ngrid_act
           idx=active(ilevel)%igrid(i)+ncoarse+(ind-1)*ngridmax
           rz_loc=rz_loc+f(idx,1)*cg_zw(idx)
        end do
     end do
  end if

  t_setup=wall_time()-tw0   ! initial guess, right-hand side norm, first residual

  iter=0; itermax=10000
  error=1.0D0; error_ini=1.0D0
  do while(iter<itermax)
     ! stock test: relative to the initial residual; warm start: relative to the right-hand side
     if(warm)then
        if(iter>0.and.error<=epsilon*rhs_norm)exit
     else
        if(error<=epsilon*error_ini)exit
     end if

     iter=iter+1

     !====================================
     ! Residual norm (local sum computed with the r update of the previous iteration)
     !====================================
     r2=r2_loc
     ! Compute global norm
#ifndef WITHOUTMPI
     tw1=wall_time()
     if(prec)then
        red2(1)=r2; red2(2)=rz_loc
        call MPI_ALLREDUCE(red2,red2_all,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
        r2=red2_all(1); rz=red2_all(2)
     else
        call MPI_ALLREDUCE(r2,r2_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
             & MPI_COMM_WORLD,info)
        r2=r2_all
     end if
     t_red=t_red+(wall_time()-tw1)
#else
     rz=rz_loc
#endif

     !====================================
     ! Compute beta factor
     !====================================
     if(iter==1)then
        beta_cg=0
     else if(prec)then
        beta_cg=rz/rz_old
     else
        beta_cg=r2/r2_old
     end if
     r2_old=r2
     rz_old=rz

     !====================================
     ! Recurrence on p
     !====================================
     if(prec)then
        do ind=1,twotondim
           do i=1,ngrid_act
              idx=active(ilevel)%igrid(i)+ncoarse+(ind-1)*ngridmax
              f(idx,2)=cg_zw(idx)+beta_cg*f(idx,2)
           end do
        end do
     else if(thr)then
!$omp parallel do collapse(2) private(ind,i,idx) schedule(static)
        do ind=1,twotondim
           do i=1,ngrid_act
              idx=active(ilevel)%igrid(i)+ncoarse+(ind-1)*ngridmax
              f(idx,2)=f(idx,1)+beta_cg*f(idx,2)
           end do
        end do
!$omp end parallel do
     else
        do ind=1,twotondim
           do i=1,ngrid_act
              idx=active(ilevel)%igrid(i)+ncoarse+(ind-1)*ngridmax
              f(idx,2)=f(idx,1)+beta_cg*f(idx,2)
           end do
        end do
     end if
     ! Update boundaries
     tw1=wall_time()
     call make_virtual_fine_dp(f(1,2),ilevel)
     t_comm=t_comm+(wall_time()-tw1)

     !==============================================
     ! Compute z = Ap and store it into f(i,3)
     !==============================================
     if(thr)then
        part=0.0d0
!$omp parallel private(ind,i,idx,tid,psum)
        call cmp_Ap_cg(ilevel)      ! orphaned worksharing loop (implicit barrier at its end)
        ! p.Ap scalar product
        tid=0
!$      tid=omp_get_thread_num()
        psum=0.0d0
!$omp do collapse(2) schedule(static)
        do ind=1,twotondim
           do i=1,ngrid_act
              idx=active(ilevel)%igrid(i)+ncoarse+(ind-1)*ngridmax
              psum=psum+f(idx,2)*f(idx,3)
           end do
        end do
!$omp end do nowait
        part(tid)=psum
!$omp end parallel
        pAp=0.0d0
        do tid=0,255
           pAp=pAp+part(tid)
        end do
     else
        call cmp_Ap_cg(ilevel)      ! outside a parallel region: its worksharing loop runs serially
        pAp=0.0d0
        do ind=1,twotondim
           do i=1,ngrid_act
              idx=active(ilevel)%igrid(i)+ncoarse+(ind-1)*ngridmax
              pAp=pAp+f(idx,2)*f(idx,3)
           end do
        end do
     end if
     ! Compute global sum
#ifndef WITHOUTMPI
     tw1=wall_time()
     call MPI_ALLREDUCE(pAp,pAp_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
          & MPI_COMM_WORLD,info)
     pAp=pAp_all
     t_red=t_red+(wall_time()-tw1)
#endif

     !====================================
     ! Compute alpha factor
     !====================================
     alpha_cg = r2/pAp
     if(prec)alpha_cg = rz/pAp

     !====================================
     ! Recurrence on x
     !====================================
     ! Recurrence on r, and r.r for the next iteration (same loop)
     !====================================
     if(thr)then
        part=0.0d0
!$omp parallel private(ind,i,idx,tid,psum)
        tid=0
!$      tid=omp_get_thread_num()
        psum=0.0d0
!$omp do collapse(2) schedule(static)
        do ind=1,twotondim
           do i=1,ngrid_act
              idx=active(ilevel)%igrid(i)+ncoarse+(ind-1)*ngridmax
              phi(idx)=phi(idx)+alpha_cg*f(idx,2)
              f(idx,1)=f(idx,1)-alpha_cg*f(idx,3)
              psum=psum+f(idx,1)*f(idx,1)
           end do
        end do
!$omp end do nowait
        part(tid)=psum
!$omp end parallel
        r2_loc=0.0d0
        do tid=0,255
           r2_loc=r2_loc+part(tid)
        end do
     else
        r2_loc=0.0d0
        do ind=1,twotondim
           do i=1,ngrid_act
              idx=active(ilevel)%igrid(i)+ncoarse+(ind-1)*ngridmax
              phi(idx)=phi(idx)+alpha_cg*f(idx,2)
              f(idx,1)=f(idx,1)-alpha_cg*f(idx,3)
              r2_loc=r2_loc+f(idx,1)*f(idx,1)
           end do
        end do
     end if

     ! Preconditioner: z = M^-1 r, and r.z for the next iteration
     if(prec)then
        call cg_precond_sgs(ilevel)
        rz_loc=0.0d0
        do ind=1,twotondim
           do i=1,ngrid_act
              idx=active(ilevel)%igrid(i)+ncoarse+(ind-1)*ngridmax
              rz_loc=rz_loc+f(idx,1)*cg_zw(idx)
           end do
        end do
     end if

     ! Compute error
     error=DSQRT(r2/dble(twotondim*numbtot(1,ilevel)))
     if(iter==1)error_ini=error
     if(verbose)write(*,112)iter,error/rhs_norm,error/error_ini

  end do
  ! End main iteration loop

  if(myid==1)write(*,115)ilevel,iter,error/rhs_norm,error/error_ini
  if(iter >= itermax)then
     if(myid==1)write(*,*)'Poisson failed to converge...'
     cg_failed(ilevel)=.true.   ! cg_auto: fall back to multigrid for this level
  end if

  if(cg_warm_start)cg_warm_ok(ilevel)=iter<itermax

  ! Update boundaries
  call make_virtual_fine_dp(phi(1),ilevel)

  cg_stat(1,ilevel)=cg_stat(1,ilevel)+t_comm
  cg_stat(2,ilevel)=cg_stat(2,ilevel)+t_red
  cg_stat(3,ilevel)=cg_stat(3,ilevel)+(wall_time()-tw0)
  cg_stat(4,ilevel)=cg_stat(4,ilevel)+1d0
  cg_stat(5,ilevel)=cg_stat(5,ilevel)+dble(iter)
  cg_stat(6,ilevel)=cg_stat(6,ilevel)+t_setup

111 format('   Entering phi_fine_cg for level ',I2)
112 format('   ==> Step=',i5,' Error=',2(1pe10.3,1x))
115 format('   ==> Level=',i5,' Step=',i5,' Error=',2(1pe10.3,1x))

end subroutine phi_fine_cg
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cg_precond_table(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------
  ! Neighbour table of the preconditioner, built once per solve: for the k-th cell of the sweep
  ! (grids of this rank in the order of the active list, the 8 cells of a grid in order),
  ! cg_cell(k) is its cell index and cg_nb(1:6,k) those of its neighbours. A neighbour that is
  ! not a cell of this rank on this level (other rank, coarser cell, boundary) points to element
  ! 0 of the work arrays, which stays zero.
  !------------------------------------------------------------------
  integer::i,j,g,ig,idim,ind,n,k
  integer,dimension(0:twondim)::nbg
  integer,dimension(1:3,1:2,1:8)::iii,jjj

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  n=active(ilevel)%ngrid
  if(allocated(cg_cell))then
     if(size(cg_cell)<twotondim*n)deallocate(cg_cell,cg_nb)
  end if
  if(.not.allocated(cg_cell))allocate(cg_cell(1:twotondim*n+twotondim*n/4+64),cg_nb(1:twondim,1:twotondim*n+twotondim*n/4+64))

  k=0
  do i=1,n
     g=active(ilevel)%igrid(i)
     nbg(0)=g
     do j=1,twondim
        ig=son(nbor(g,j))
        if(ig>0)then
           if(cg_own(ig)/=cg_stamp)ig=0
        end if
        nbg(j)=ig
     end do
     do ind=1,twotondim
        k=k+1
        cg_cell(k)=g+ncoarse+(ind-1)*ngridmax
        do idim=1,ndim
           ig=nbg(iii(idim,1,ind))
           cg_nb(2*idim-1,k)=0
           if(ig>0)cg_nb(2*idim-1,k)=ig+ncoarse+(jjj(idim,1,ind)-1)*ngridmax
           ig=nbg(iii(idim,2,ind))
           cg_nb(2*idim,k)=0
           if(ig>0)cg_nb(2*idim,k)=ig+ncoarse+(jjj(idim,2,ind)-1)*ngridmax
        end do
     end do
  end do

end subroutine cg_precond_table
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cg_precond_sgs(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------
  ! Preconditioner of phi_fine_cg (cg_precond=1): cg_zw = M^-1 r (r in f(:,1)), with M one SSOR sweep
  ! (relaxation factor cg_precond_omega; 1 = symmetric Gauss-Seidel) of B = I - (1/6) sum over
  ! the neighbours, restricted to the cells of this rank (table of cg_precond_table).
  ! Neighbours of other ranks and outside the level count as zero, so nothing is communicated and
  ! M is symmetric positive definite. The sign differs from A = -B, which CG does not see
  ! (z -> -z flips p and alpha together). Forward sweep into cg_yw, backward sweep into cg_zw:
  ! cells not yet swept hold zero, so each sweep adds all neighbours without testing the order.
  !------------------------------------------------------------------
  integer::k,nc,idx
  real(dp)::oneoversix,om,om6,om2

  oneoversix=1.0D0/dble(twondim)
  om=cg_precond_omega
  om6=om*oneoversix
  om2=2.0d0-om
  nc=twotondim*active(ilevel)%ngrid

  do k=1,nc
     idx=cg_cell(k)
     cg_yw(idx)=0.0d0
     cg_zw(idx)=0.0d0
  end do

  ! Forward sweep: (I/om - L) y = r
  do k=1,nc
     idx=cg_cell(k)
     cg_yw(idx)=om*f(idx,1)+om6*(cg_yw(cg_nb(1,k))+cg_yw(cg_nb(2,k))+cg_yw(cg_nb(3,k)) &
          &                     +cg_yw(cg_nb(4,k))+cg_yw(cg_nb(5,k))+cg_yw(cg_nb(6,k)))
  end do

  ! Backward sweep: (I/om - U) z = ((2-om)/om) y
  do k=nc,1,-1
     idx=cg_cell(k)
     cg_zw(idx)=om2*cg_yw(idx)+om6*(cg_zw(cg_nb(1,k))+cg_zw(cg_nb(2,k))+cg_zw(cg_nb(3,k)) &
          &                        +cg_zw(cg_nb(4,k))+cg_zw(cg_nb(5,k))+cg_zw(cg_nb(6,k)))
  end do

end subroutine cg_precond_sgs
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_residual_cg(ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use constants, only: twopi
  implicit none
  integer::ilevel,icount
  !------------------------------------------------------------------
  ! This routine computes the residual for the Conjugate Gradient
  ! Poisson solver. The residual is stored in f(i,1).
  !------------------------------------------------------------------
  integer::i,idim,igrid,ngrid,ncache,ind,iskip,nx_loc
  integer::id1,id2,ig1,ig2,ih1,ih2,nn
  real(dp)::dx2,fourpi,scale,oneoversix,fact
  integer,dimension(1:3,1:2,1:8)::iii,jjj

  ! Work arrays are ordinary locals (not save): the sweep loop is threaded (omp_poisson). Every
  ! sweep writes only f(:,1:2) of the cells of its own grids, so the result does not depend on threads.
  integer ,dimension(1:nvector)::ind_grid,ind_cell,ind_tmp,ind_map
  integer ,dimension(1:nvector,0:twondim)::igridn
  integer ,dimension(1:nvector,1:ndim)::ind_left,ind_right
  real(dp),dimension(1:nvector,1:ndim)::phig,phid
  real(dp),dimension(1:nvector,1:twotondim,1:ndim)::phi_left,phi_right
  real(dp),dimension(1:nvector,1:twotondim)::phi_tmp
  real(dp),dimension(1:nvector)::residu

  ! Set constants
  dx2=(0.5D0**ilevel)**2
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  fourpi=2*twopi*scale
  if(cosmo)fourpi=1.5D0*omega_m*aexp*scale
  oneoversix=1.0D0/dble(twondim)
  fact=oneoversix*fourpi*dx2

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
!$omp parallel do schedule(static) if(omp_poisson .and. ncache>=omp_mg_mingrid) &
!$omp & private(igrid,ngrid,i,idim,ind,iskip,id1,id2,ig1,ig2,ih1,ih2,nn,ind_grid,ind_cell,ind_tmp,ind_map) &
!$omp & private(igridn,ind_left,ind_right,phig,phid,phi_left,phi_right,phi_tmp,residu)
  do igrid=1,ncache,nvector

     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Gather neighboring grids
     do i=1,ngrid
        igridn(i,0)=ind_grid(i)
     end do
     do idim=1,ndim
        do i=1,ngrid
           ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
           ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
           igridn(i,2*idim-1)=son(ind_left (i,idim))
           igridn(i,2*idim  )=son(ind_right(i,idim))
        end do
     end do

     ! Interpolate potential from upper level, only for the grids without a neighbouring grid in
     ! that direction (the others take phi of the neighbour below; the stock code interpolated for
     ! every grid and discarded most of it). interpol_phi works element by element: same values.
     do idim=1,ndim
        nn=0
        do i=1,ngrid
           if(igridn(i,2*idim-1)==0)then
              nn=nn+1; ind_map(nn)=i; ind_tmp(nn)=ind_left(i,idim)
           end if
        end do
        if(nn>0)then
           call interpol_phi(ind_tmp,phi_tmp,nn,ilevel,icount)
           do ind=1,twotondim
              do i=1,nn
                 phi_left(ind_map(i),ind,idim)=phi_tmp(i,ind)
              end do
           end do
        end if
        nn=0
        do i=1,ngrid
           if(igridn(i,2*idim)==0)then
              nn=nn+1; ind_map(nn)=i; ind_tmp(nn)=ind_right(i,idim)
           end if
        end do
        if(nn>0)then
           call interpol_phi(ind_tmp,phi_tmp,nn,ilevel,icount)
           do ind=1,twotondim
              do i=1,nn
                 phi_right(ind_map(i),ind,idim)=phi_tmp(i,ind)
              end do
           end do
        end if
     end do

     ! Loop over cells
     do ind=1,twotondim
        ! Gather neighboring potential
        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           ih1=ncoarse+(id1-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig1)>0)then
                 phig(i,idim)=phi(igridn(i,ig1)+ih1)
              else
                 phig(i,idim)=phi_left(i,id1,idim)
              end if
           end do
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           ih2=ncoarse+(id2-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig2)>0)then
                 phid(i,idim)=phi(igridn(i,ig2)+ih2)
              else
                 phid(i,idim)=phi_right(i,id2,idim)
              end if
           end do
        end do

        ! Compute central cell index
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Compute residual using 6 neighbors potential
        do i=1,ngrid
           residu(i)=phi(ind_cell(i))
        end do
        do idim=1,ndim
           do i=1,ngrid
              residu(i)=residu(i)-oneoversix*(phig(i,idim)+phid(i,idim))
           end do
        end do
        do i=1,ngrid
           residu(i)=residu(i)+fact*(rho(ind_cell(i))-rho_tot)
        end do

        ! Store results in f(i,1)
        do i=1,ngrid
           f(ind_cell(i),1)=residu(i)
        end do

        ! Store results in f(i,2)
        do i=1,ngrid
           f(ind_cell(i),2)=residu(i)
        end do

     end do
     ! End loop over cells

  end do
!$omp end parallel do
  ! End loop over grids

end subroutine cmp_residual_cg
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_Ap_cg(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------
  ! This routine computes Ap for the Conjugate Gradient
  ! Poisson Solver and store the result into f(i,3).
  !------------------------------------------------------------------
  integer::i,idim,igrid,ngrid,ncache,ind,iskip
  integer::id1,id2,ig1,ig2,ih1,ih2
  real(dp)::oneoversix
  integer,dimension(1:3,1:2,1:8)::iii,jjj

  ! Work arrays are ordinary locals (not save): the routine is called inside a parallel region
  ! of phi_fine_cg and its grid loop is an orphaned worksharing loop.
  integer,dimension(1:nvector)::ind_grid,ind_cell
  integer,dimension(1:nvector,0:twondim)::igridn
  real(dp),dimension(1:nvector,1:ndim)::phig,phid
  real(dp),dimension(1:nvector)::residu

  ! Set constants
  oneoversix=1.0D0/dble(twondim)

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
!$omp do private(igrid,ngrid,i,idim,ind,iskip,id1,id2,ig1,ig2,ih1,ih2) schedule(static)
  do igrid=1,ncache,nvector

     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Gather neighboring grids
     do i=1,ngrid
        igridn(i,0)=ind_grid(i)
     end do
     do idim=1,ndim
        do i=1,ngrid
           igridn(i,2*idim-1)=son(nbor(ind_grid(i),2*idim-1))
           igridn(i,2*idim  )=son(nbor(ind_grid(i),2*idim  ))
        end do
     end do

     ! Loop over cells
     do ind=1,twotondim

        ! Gather neighboring potential
        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           ih1=ncoarse+(id1-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig1)>0)then
                 phig(i,idim)=f(igridn(i,ig1)+ih1,2)
              else
                 phig(i,idim)=0
              end if
           end do
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           ih2=ncoarse+(id2-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig2)>0)then
                 phid(i,idim)=f(igridn(i,ig2)+ih2,2)
              else
                 phid(i,idim)=0
              end if
           end do
        end do

        ! Compute central cell index
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Compute Ap using neighbors potential
        do i=1,ngrid
           residu(i)=-f(ind_cell(i),2)
        end do
        do idim=1,ndim
           do i=1,ngrid
              residu(i)=residu(i)+oneoversix*(phig(i,idim)+phid(i,idim))
           end do
        end do
        ! Store results in f(i,3)
        do i=1,ngrid
           f(ind_cell(i),3)=residu(i)
        end do

     end do
     ! End loop over cells

  end do
!$omp end do
  ! End loop over grids

end subroutine cmp_Ap_cg
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine make_warm_phi(ilevel,icount)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
  integer::ilevel,icount
  !------------------------------------------------------------------
  ! cg_warm_start: initial guess of phi_fine_cg from the previous solution.
  ! phi itself was used as work space by rho_fine; phi_old is the copy made
  ! just before (save_phi_old in amr_step). A grid created since the last
  ! solve has no solution: make_grid_fine gave its eight cells the value of
  ! the father cell. Such grids (eight identical values) get the
  ! interpolation from level-1, as in make_initial_phi.
  !------------------------------------------------------------------
  integer::igrid,ncache,i,ngrid,ind,iskip,nnew
  integer ,dimension(1:nvector)::ind_grid,ind_cell_father
  logical ,dimension(1:nvector)::isnew
  real(dp),dimension(1:nvector,1:twotondim)::phi_int
  real(dp)::p1

  ncache=active(ilevel)%ngrid
!$omp parallel do schedule(static) if(omp_poisson .and. ncache>=omp_mg_mingrid) &
!$omp & private(igrid,ngrid,i,ind,iskip,nnew,ind_grid,ind_cell_father,isnew,phi_int,p1)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     isnew(1:ngrid)=.true.
     do ind=2,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           p1=phi_old(ncoarse+ind_grid(i))
           if(phi_old(iskip+ind_grid(i))/=p1)isnew(i)=.false.
        end do
     end do
     nnew=count(isnew(1:ngrid))
     if(nnew>0)then
        do i=1,ngrid
           ind_cell_father(i)=father(ind_grid(i))
        end do
        call interpol_phi(ind_cell_father,phi_int,ngrid,ilevel,icount)
     end if
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           if(isnew(i))then
              phi(iskip+ind_grid(i))=phi_int(i,ind)
           else
              phi(iskip+ind_grid(i))=phi_old(iskip+ind_grid(i))
           end if
        end do
     end do
  end do
!$omp end parallel do

end subroutine make_warm_phi
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine make_initial_phi(ilevel,icount)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
  integer::ilevel,icount
  !
  !
  !
  integer::igrid,ncache,i,ngrid,ind,iskip,idim
  ! not save: the sweep loop is threaded (omp_poisson); every sweep writes only its own grids
  integer ,dimension(1:nvector)::ind_grid,ind_cell,ind_cell_father
  real(dp),dimension(1:nvector,1:twotondim)::phi_int

  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
!$omp parallel do schedule(static) if(omp_poisson .and. ncache>=omp_mg_mingrid) &
!$omp & private(igrid,ngrid,i,ind,iskip,idim,ind_grid,ind_cell,ind_cell_father,phi_int)
  do igrid=1,ncache,nvector
     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     if(ilevel==1)then
        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           do i=1,ngrid
              phi(ind_cell(i))=0
           end do
           do idim=1,ndim
              do i=1,ngrid
                 f(ind_cell(i),idim)=0
              end do
           end do
        end do
        ! End loop over cells
     else
        ! Compute father cell index
        do i=1,ngrid
           ind_cell_father(i)=father(ind_grid(i))
        end do

        ! Interpolate
        call interpol_phi(ind_cell_father,phi_int,ngrid,ilevel,icount)

        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           do i=1,ngrid
              phi(ind_cell(i))=phi_int(i,ind)
           end do
           do idim=1,ndim
              do i=1,ngrid
                 f(ind_cell(i),idim)=0
              end do
           end do
        end do
        ! End loop over cells
     end if

  end do
!$omp end parallel do
  ! End loop over grids

end subroutine make_initial_phi
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine make_multipole_phi(ilevel)
  use amr_commons
  use pm_commons
  use poisson_commons
  use constants, only: twopi
  implicit none
  integer::ilevel

  integer::idim
  integer::i,ncache,igrid,ngrid,ind
  integer::iskip,nx_loc,ix,iy,iz
  integer,dimension(1:nvector),save::ind_grid,ind_cell

  real(dp)::dx,dx_loc,scale,boxlen2,eps
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector),save::rr,pp
  real(dp),dimension(1:nvector,1:ndim),save::xx

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  boxlen2=boxlen**2
  eps=dx_loc

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
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
           ind_cell(i)=iskip+ind_grid(i)
        end do

        if(simple_boundary)then
           ! Compute cell center in code units
           do idim=1,ndim
               do i=1,ngrid
                  xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
               end do
           end do

           ! Rescale position from code units to user units
           rr(1:ngrid)=0.0d0
           do idim=1,ndim
               do i=1,ngrid
                  xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
                  rr(i)=rr(i)+(xx(i,idim)-multipole(idim+1)/multipole(1))**2
               end do
           end do

           do i=1,ngrid
               rr(i)=max(eps,sqrt(rr(i)))       ! Cutoff
           end do

           if(ngrid>0) call phi_ana(rr,pp,ngrid)

           ! Scatter variables
           do i=1,ngrid
               phi(ind_cell(i))=pp(i)/scale
           end do

        else
           do i=1,ngrid
               phi(ind_cell(i))=0
           end do
        endif

        ! End loop over cells
     end do

  end do
  ! End loop over grids

end subroutine make_multipole_phi
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cg_timer_output
  use amr_commons
  use poisson_commons
  use mpi_mod
  implicit none
  !------------------------------------------------------------------
  ! Print the time breakdown of phi_fine_cg accumulated over the run
  ! (cg_stat), per level: number of solves, iterations per solve, time
  ! of the solves (average over ranks), the set-up per solve (initial
  ! guess, first residual) and, per iteration, the ghost exchange, the
  ! two MPI_ALLREDUCE and the local work, as average / maximum over ranks.
  !------------------------------------------------------------------
  integer::ilevel
#ifndef WITHOUTMPI
  integer::info
#endif
  real(dp),dimension(1:6,1:nlevelmax)::sav,smx
  real(dp),dimension(1:3)::av,mx
  real(dp)::nit,nsol

  if(.not.poisson)return
  if(.not.allocated(cg_stat))return
  ! local work of the iterations = whole solve - set-up - exchange - reduce
  cg_stat(3,:)=cg_stat(3,:)-cg_stat(1,:)-cg_stat(2,:)-cg_stat(6,:)
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(cg_stat,sav,6*nlevelmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(cg_stat,smx,6*nlevelmax,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
  sav=sav/dble(ncpu)
#else
  sav=cg_stat; smx=cg_stat
#endif
  cg_stat(3,:)=cg_stat(3,:)+cg_stat(1,:)+cg_stat(2,:)+cg_stat(6,:)
  if(myid/=1)return
  if(sum(smx(4,:))==0d0)return
  write(*,'(" CG solver breakdown (ms, average / maximum over ranks; set-up per solve, the others per iteration)")')
  write(*,'(" level  solves  iter/solve  total(s)      set-up          exchange          allreduce         local")')
  do ilevel=1,nlevelmax
     if(smx(4,ilevel)==0d0)cycle
     nit=max(smx(5,ilevel),1d0); nsol=smx(4,ilevel)
     av=1d3*sav(1:3,ilevel)/nit; mx=1d3*smx(1:3,ilevel)/nit
     write(*,'(I6,I8,F12.1,F10.3,4(F9.4," /",F7.4))')ilevel,nint(nsol), &
          & smx(5,ilevel)/nsol,sum(sav(1:3,ilevel))+sav(6,ilevel),1d3*sav(6,ilevel)/nsol,1d3*smx(6,ilevel)/nsol, &
          & av(1),mx(1),av(2),mx(2),av(3),mx(3)
  end do

end subroutine cg_timer_output
