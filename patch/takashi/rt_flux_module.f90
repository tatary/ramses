MODULE rt_flux_module

  use rt_parameters
  implicit none

  private   ! default

  public read_hll_eigenvalues, cmp_rt_faces

CONTAINS

!************************************************************************
SUBROUTINE read_hll_eigenvalues()

! Read M1 eigenvalues for hll method from file 'hll_evals.list'
! into              in the 3*4 flux tensor
!------------------------------------------------------------------------
  use amr_commons,only:myid
  use rt_parameters,only:hll_evals_file
  use mpi_mod
#ifndef WITHOUTMPI
  use amr_commons,only:IOGROUPSIZE,ncpu
  integer::dummy_io,info2
#endif
  integer::i,j,ii,jj
  real(dp)::dummy
  logical::ok
  integer,parameter::tag=1128

!------------------------------------------------------------------------
  if(myid==1) write(*,*) 'Reading HLL eigenvalues from file'
  if(hll_evals_file.eq.'') &
          call get_environment_variable('RAMSES_HLLFILE', hll_evals_file)
  inquire(file=TRIM(hll_evals_file), exist=ok)
  if(.not. ok)then
     if(myid.eq.1) then
        write(*,*)'Cannot read hll eigenvalues file...'
        write(*,*)'File '//TRIM(hll_evals_file)//' not found'
        write(*,*)'You need to set the RAMSES_HLLFILE envvar' // &
                  ' to the correct path, or use the namelist.'
     endif
     call clean_stop
  end if

  ! Wait for the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if (mod(myid-1,IOGROUPSIZE)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
     end if
  endif
#endif


  open(unit=10,file=TRIM(hll_evals_file),status='old',form='formatted')
  read(10,*)i
  !write(*,*)i
  allocate(lambda1(0:100,0:100)) ; allocate(lambda4(0:100,0:100))
  do i=0,100
     do j=0,100
        read(10,*)ii,jj,lambda1(ii,jj),dummy,dummy,lambda4(ii,jj)
     end do
  end do
  close(10)

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


END SUBROUTINE read_hll_eigenvalues

!************************************************************************
SUBROUTINE cmp_eigenvals(uin, iP0, ngrid, lmin, lmax)

!  Compute Jacobian eigenvalues for given vector of sub-grids.
!
!  inputs/outputs
!  uin         => input cell states
!  ngrid       => number of sub-grids of 3^ndim cells
!  lmin       <=  return minimum cell eigenvalues
!  lmax       <=  return maximum cell eigenvalues
!
!  other vars
!  iu1,iu2     |first and last index of input array,
!  ju1,ju2     |cell centered,
!  ku1,ku2     |including buffer cells.
!------------------------------------------------------------------------
  real(dp), dimension(nvector, iu1:iu2, ju1:ju2, ku1:ku2, nrtvar),  &
                                                          intent(in)::uin
  real(dp),dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,ndim)::  lmin, lmax
  integer,intent(in)::iP0, ngrid!----------------------------------------
  real(dp), dimension(ndim)::f, costheta
  real(dp)::ff, np
  integer::i, j, k, n, id, nedge
!------------------------------------------------------------------------
  do k=kfrt1,kf2                                !
  do j=jfrt1,jf2                                !  Loop each cell in grid
  do i=ifrt1,if2                                !

     nedge=0                ! Check if we're at a corner and if so, cycle
     if(mod(i,3).eq.0) nedge=nedge+1
     if(ndim.gt.1 .and. mod(j,3).eq.0) nedge=nedge+1
     if(ndim.gt.2 .and. mod(k,3).eq.0) nedge=nedge+1
     if(nedge.ge.2) cycle

  do n=1,ngrid                                     ! Loop buffer of grids
     ! calculate eigenvals
     np = uin(n, i, j, k, iP0)              !              Photon density
     f =  uin(n, i, j, k, iP0+1:iP0+nDim)   !   Ph fluxes, all directions
     ff=sqrt(sum(f(:)**2))                  !              Flux magnitude

     if(ff.gt.0.)then
        costheta=f/ff
     else
        costheta=0d0
     endif
     if(np>0d0)then
        ff=ff/rt_c/np
     else
        ff=0.0d0
     endif
     ff=max(min(ff,1d0),0d0)
     do id=1, ndim
        costheta(id)=max(min(costheta(id),1d0),-1d0)
     end do
     do id=1, ndim
        call inp_eigenvals(ff,costheta(id), lmin(n, i, j, k, id),&
             lmax(n, i, j, k, id))
     end do
  end do
  end do
  end do
  end do

END SUBROUTINE cmp_eigenvals

!************************************************************************
SUBROUTINE inp_eigenvals(ff, omega, lmin, lmax)

! Compute eigenvalues by interpolation from table
! input/output:
! ff        => photon flux magnitude in cell
! omega     => angle of flux with cell face in question
! lmin     <=  Smallest eigenvalue in cell for given direction
! lmax      => Biggest eigenvalue in cell for given direction
!------------------------------------------------------------------------
 use amr_commons
 real(dp), intent(in):: ff,omega
 real(dp):: lmin, lmax
 real(dp)::theta,pi,dd1,dd2,de1,de2,lff,ltt
 integer::ii,jj
!------------------------------------------------------------------------
 theta=ACOS(omega)
 pi=ACOS(-1.0d0)
 lff=ff*100.
 ltt=theta/pi*100.

 ii  = MIN(int(lff),99)
 jj  = MIN(int(ltt),99)
 dd1 = lff-float(ii)
 dd2 = ltt-float(jj)
 de1 = 1.0-dd1
 de2 = 1.0-dd2

 lmin=0.0
 lmin=lmin+de1*de2*lambda1(ii,jj)
 lmin=lmin+dd1*de2*lambda1(ii+1,jj)
 lmin=lmin+de1*dd2*lambda1(ii,jj+1)
 lmin=lmin+dd1*dd2*lambda1(ii+1,jj+1)

 lmax=0.0
 lmax=lmax+de1*de2*lambda4(ii,jj)
 lmax=lmax+dd1*de2*lambda4(ii+1,jj)
 lmax=lmax+de1*dd2*lambda4(ii,jj+1)
 lmax=lmax+dd1*dd2*lambda4(ii+1,jj+1)

END SUBROUTINE inp_eigenvals

!************************************************************************
subroutine cmp_flux_tensors(uin, iP0, nGrid, F)

! Compute central fluxes for a photon group, for each cell in a vector
! of grids.
! The flux tensor is a three by four tensor (2*3 and 1*2 in 1D and 2D,
! respectively) where the first column is photon flux (x,y,z) and
! the other three columns compose the Eddington tensor (see Aubert &
! Teyssier '08), times c^2.
! input/output:
! uin       => RT variables of all cells in a vector of grids
!              (photon energy densities and photon fluxes).
! iP0       => Starting index of photon group among the RT variables.
! ngrid     => Number of 'valid' grids in uin.
! F        <=  Group flux tensors for all the cells.
!------------------------------------------------------------------------
  real(dp),dimension(1:nvector, iu1:iu2, ju1:ju2, ku1:ku2, 1:nrtvar)::    uin
  real(dp),dimension(1:nvector, iu1:iu2, ju1:ju2, ku1:ku2, 1:nDim+1, 1:ndim)::F
  integer::iP0, nGrid!---------------------------------------------------
  real(dp)::Np, Np_c_sq, pflux_sq, chi, iterm, oterm, den, rn, fx, ux
#if NDIM>1
  real(dp)::fy, uy
#endif
#if NDIM>2
  real(dp)::fz, uz
#endif
  logical ::nonzero, posflux
  integer::i, j, k, n, nedge
!------------------------------------------------------------------------
  ! Loop 4X4X4 cells in grid. All go from 0 to 3, out of 6X6X6 cells.
  ! We only need to calculate tensors for those cells which have faces to
  ! the 2X2X2 center cells, so by skipping the 'corners' we are reduced
  ! to half of the 4X4X4 cells.
  !
  ! The arithmetic below is the original (array/loop) formulation written
  ! out per component, in the same order of operations, so that the
  ! innermost loop over grids has no branches and vectorizes. The zero
  ! density and zero flux cases are handled with merge instead of cycle.

  do k = kfrt1, kf2
  do j = jfrt1, jf2
  do i = ifrt1, if2

     nedge=0                 ! Check if we're at a corner and if so, cycle
     if(mod(i,3).eq.0) nedge=nedge+1
     if(ndim.gt.1 .and. mod(j,3).eq.0) nedge=nedge+1
     if(ndim.gt.2 .and. mod(k,3).eq.0) nedge=nedge+1
     if(nedge.ge.2) cycle

     if(any(uin(1:ngrid, i, j, k, iP0) .lt. 0d0)) then
        write(*,*)'negative photon density in cmp_eddington. -EXITING-'
        call clean_stop
     endif

     do n=1,ngrid                        !        Also loop over all grids
        Np = uin(n, i, j, k, iP0)        !          Photon density in cell
        fx = uin(n, i, j, k, iP0+1)      !                    Photon flux
        F(n,i,j,k,1,1) = fx              !   First col is photon flux
        pflux_sq = fx*fx                 !  Sq. photon flux magnitude
#if NDIM>1
        fy = uin(n, i, j, k, iP0+2)
        F(n,i,j,k,1,2) = fy
        pflux_sq = pflux_sq + fy*fy
#endif
#if NDIM>2
        fz = uin(n, i, j, k, iP0+3)
        F(n,i,j,k,1,3) = fz
        pflux_sq = pflux_sq + fz*fz
#endif
        Np_c_sq = Np*rt_c2*Np
        nonzero = Np_c_sq .ne. 0d0           !Zero density => no pressure
        den = merge(Np_c_sq, 1d0, nonzero)
        posflux = pflux_sq .gt. 0d0
        rn = merge(sqrt(pflux_sq), 1d0, posflux)
        ux = merge(fx/rn, 0d0, posflux)      !           Flux unit vector
#if NDIM>1
        uy = merge(fy/rn, 0d0, posflux)
#endif
#if NDIM>2
        uz = merge(fz/rn, 0d0, posflux)
#endif
        pflux_sq = pflux_sq/den               !      Reduced flux, squared
        chi = max(4d0-3d0*pflux_sq, 0d0)   !           Eddington factor
        chi = (3d0+ 4d0*pflux_sq)/(5d0 + 2d0*sqrt(chi))

        iterm = (1d0-chi)/2d0               !    Identity term in tensor
        oterm = (3d0*chi-1d0)/2d0          !         Outer product term
        ! Eddington tensor: (oterm*u(p)*u(q) [+ iterm on the diagonal])
        ! times c^2*Np, zero where Np_c_sq is zero
        F(n,i,j,k,2,1) = merge((oterm*ux*ux + iterm)*rt_c2*Np, 0d0, nonzero)
#if NDIM>1
        F(n,i,j,k,2,2) = merge((oterm*ux*uy)*rt_c2*Np, 0d0, nonzero)
        F(n,i,j,k,3,1) = merge((oterm*uy*ux)*rt_c2*Np, 0d0, nonzero)
        F(n,i,j,k,3,2) = merge((oterm*uy*uy + iterm)*rt_c2*Np, 0d0, nonzero)
#endif
#if NDIM>2
        F(n,i,j,k,2,3) = merge((oterm*ux*uz)*rt_c2*Np, 0d0, nonzero)
        F(n,i,j,k,3,3) = merge((oterm*uy*uz)*rt_c2*Np, 0d0, nonzero)
        F(n,i,j,k,4,1) = merge((oterm*uz*ux)*rt_c2*Np, 0d0, nonzero)
        F(n,i,j,k,4,2) = merge((oterm*uz*uy)*rt_c2*Np, 0d0, nonzero)
        F(n,i,j,k,4,3) = merge((oterm*uz*uz + iterm)*rt_c2*Np, 0d0, nonzero)
#endif
     enddo
  enddo
  enddo
  enddo

end subroutine cmp_flux_tensors

!************************************************************************
FUNCTION cmp_face(fdn, fup, udn, uup, lminus, lplus)

! Compute intercell fluxes for all (four) RT variables, using the
! Harten-Lax-van Leer method (see eq. 30 in Aubert&Teyssier(2008).
! fdn    => flux function in the cell downwards from the border
! fup    => flux function in the cell upwards from the border
! udn    => state (photon density and flux downwards from the border
! uup    => state (photon density and flux upwards from the border
! lminus => minimum intercell eigenvalue
! lplus  => maximum intercell eigenvalue
! returns      flux vector for the given state variables, i.e. line nr dim
!              in the 3*4 flux function tensor
!------------------------------------------------------------------------
  real(dp),dimension(nDim+1)::fdn, fup, udn, uup, cmp_face
  real(dp)::lminus, lplus
  real(dp)::div
!------------------------------------------------------------------------
  if(rt_use_hll) then
     div=lplus-lminus
     !if(div .le. 1d-20) then
     !   cmp_face=0d0
     !   print *,'eigenvalue problem: ',lminus,lplus
     !   return
     !endif
     cmp_face =                                                          &
          ( lplus*fdn - lminus*fup + lplus*lminus*rt_c*( uup-udn )) / div
  else
     cmp_face = ( fdn + fup - rt_c*( uup-udn )) / 2d0
  endif
  return
END FUNCTION cmp_face


!************************************************************************
SUBROUTINE cmp_rt_faces(uin, iFlx, dx, dy, dz, dt, iP0, ngrid)

!  Compute intercell fluxes for one photon group in all dimensions,
!  using the Eddington tensor with the M1 closure relation.
!  The intercell fluxes are the right-hand sides of the equations:
!      dN/dt = - nabla(F),
!      dF/dt = - nabla(c^2*P),
!  where N is photon density, F is photon flux, c the light speed and P
!  the Eddington pressure tensor. A flux at index i,j,k represents
!  flux across the lower faces of that cell, i.e. at i-1/2 etc.
!
!  inputs/outputs
!  uin         => input states
!  iFlx       <=  return fluxes in the 3 coord directions.
!  dx,dy,dz    => (dx,dy,dz)
!  dt          => time step
!  iP0         => Starting index, among the RT variables, of the group.
!  ngrid       => number of sub-grids
!
!  other vars
!  iu1,iu2     |First and last index of input array,
!  ju1,ju2     |cell centered,
!  ku1,ku2     |including buffer cells (6x6x6).
!  if1,if2     |First and last index of output array,
!  jf1,jf2     |edge centered, for active
!  kf1,kf2     |cells only (3x3x3).
!------------------------------------------------------------------------
  use amr_parameters
  use amr_commons
  use const
  real(dp),dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2, 1:nrtvar)::     uin
  real(dp),dimension(nvector,if1:if2,jf1:jf2,kf1:kf2, 1:nrtvar, 1:ndim) &
       ::iFlx
  real(dp)::dx, dy, dz, dt
  integer ::iP0, iP1, nGrid!---------------------------------------------
  real(dp),save, &                                     !   Central fluxes
           dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2, ndim+1, ndim)::cFlx
  real(dp),save, &                                     ! Cell eigenvalues
           dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2, ndim)::  lmin, lmax
  ! Upwards and downwards fluxes and states of the group
  real(dp),dimension(nDim+1),save:: fdn, fup, udn, uup
  real(dp):: lminus, lplus                        ! Intercell eigenvalues
  real(dp)::dtdx
  integer ::i, j, k, n, m
  integer ::if2c, jf2c, kf2c   ! last face index across the flux direction
!------------------------------------------------------------------------
  iP1=iP0+nDim                                ! end index of photon group

  ! Only the faces of the 2x2x2 grid cells are used by rt_godfine1:
  ! index 1..3 along the flux direction and 1..2 across it (1..1 in the
  ! directions that do not exist). The faces at index 3 across the
  ! direction belong to cells outside the grid and were never read.
  if2c=if2-(if2-if1)/2
  jf2c=jf2-(jf2-jf1)/2
  kf2c=kf2-(kf2-kf1)/2

  ! compute flux tensors for all the cells
  call cmp_flux_tensors(uin, iP0, ngrid, cFlx)    ! flux tensors of cells

  if(rt_use_hll) &
       call cmp_eigenvals(uin, iP0, ngrid, lmin, lmax) ! eigenvs of cells

  ! Solve for 1D flux in X direction
  !----------------------------------------------------------------------
  dtdx=dt/dx
  if(rt_use_hll)then
  do i=if1,if2                                 !
  do j=jf1,jf2c                                !        each cell in grid
  do k=kf1,kf2c                                !
     do n=1,ngrid                              ! <- buffer of grids
           fdn = cFlx(n,  i-1, j, k, :, 1    )    !
           fup = cFlx(n,  i,   j, k, :, 1    )   !  upwards and downwards
           udn = uin( n,  i-1, j, k, iP0:iP1 )   !  conditions
           uup = uin( n,  i,   j, k, iP0:iP1 )    !
           lminus=MIN(lmin(n,i-1,j,k,1), lmin(n,i,j,k,1), 0d0)
           lplus =MAX(lmax(n,i-1,j,k,1), lmax(n,i,j,k,1), 0d0)
           iFlx(n, i, j, k, iP0:iP1, 1)=&
                        cmp_face( fdn, fup, udn, uup, lminus, lplus)*dtdx
      end do
  end do
  end do
  end do
  else
  ! GLF flux, same arithmetic as cmp_face but written out per component
  ! so that the innermost loop over grids vectorizes.
  do m=1,ndim+1
  do k=kf1,kf2c
  do j=jf1,jf2c
  do i=if1,if2
     do n=1,ngrid
           iFlx(n,i,j,k,iP0+m-1,1)=( cFlx(n,i-1,j,k,m,1) + cFlx(n,i,j,k,m,1) &
                - rt_c*( uin(n,i,j,k,iP0+m-1)-uin(n,i-1,j,k,iP0+m-1) )) / 2d0 * dtdx
     end do
  end do
  end do
  end do
  end do
  end if

  ! Solve for 1D flux in Y direction
  !----------------------------------------------------------------------
#if NDIM>1
  dtdx=dt/dy
  if(rt_use_hll)then
  do i=if1,if2c
  do j=jf1,jf2
  do k=kf1,kf2c
     do n=1,ngrid
           fdn = cFlx(n,  i, j-1, k, :, 2    )
           fup = cFlx(n,  i, j,   k, :, 2    )
           udn = uin( n,  i, j-1, k, iP0:iP1 )
           uup = uin( n,  i, j,   k, iP0:iP1 )
           lminus=MIN(lmin(n,i,j-1,k,2), lmin(n,i,j,k,2), 0d0)
           lplus =MAX(lmax(n,i,j-1,k,2), lmax(n,i,j,k,2), 0d0)
           iFlx(n,i,j,k,iP0:iP1,2)=&
                       cmp_face( fdn, fup, udn, uup, lminus, lplus )*dtdx
     end do
  end do
  end do
  end do
  else
  do m=1,ndim+1
  do k=kf1,kf2c
  do j=jf1,jf2
  do i=if1,if2c
     do n=1,ngrid
           iFlx(n,i,j,k,iP0+m-1,2)=( cFlx(n,i,j-1,k,m,2) + cFlx(n,i,j,k,m,2) &
                - rt_c*( uin(n,i,j,k,iP0+m-1)-uin(n,i,j-1,k,iP0+m-1) )) / 2d0 * dtdx
     end do
  end do
  end do
  end do
  end do
  end if
#endif

  ! Solve for 1D flux in Z direction
  !----------------------------------------------------------------------
#if NDIM>2
  dtdx=dt/dz
  if(rt_use_hll)then
  do i=if1,if2c
  do j=jf1,jf2c
  do k=kf1,kf2
     do n=1,ngrid
           fdn = cFlx(n,  i, j, k-1, :, 3    )
           fup = cFlx(n,  i, j, k,   :, 3    )
           udn = uin( n,  i, j, k-1, iP0:iP1 )
           uup = uin( n,  i, j, k,   iP0:iP1 )
           lminus=MIN(lmin(n,i,j,k-1,3), lmin(n,i,j,k,3), 0d0)
           lplus =MAX(lmax(n,i,j,k-1,3), lmax(n,i,j,k,3), 0d0)
           iFlx(n,i,j,k,iP0:iP1,3)=&
                       cmp_face( fdn, fup, udn, uup, lminus, lplus )*dtdx
      end do
  end do
  end do
  end do
  else
  do m=1,ndim+1
  do k=kf1,kf2
  do j=jf1,jf2c
  do i=if1,if2c
     do n=1,ngrid
           iFlx(n,i,j,k,iP0+m-1,3)=( cFlx(n,i,j,k-1,m,3) + cFlx(n,i,j,k,m,3) &
                - rt_c*( uin(n,i,j,k,iP0+m-1)-uin(n,i,j,k-1,iP0+m-1) )) / 2d0 * dtdx
     end do
  end do
  end do
  end do
  end do
  end if
#endif

end subroutine cmp_rt_faces


END MODULE rt_flux_module
