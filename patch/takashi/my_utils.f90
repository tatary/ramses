!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
! my_utils: helpers of the FG (patch/ks) star formation: ufront (hydro state at the front of the
! update, for the non-local gas depletion and the SFC model) and position -> cpu / cell lookups.
! Ported from patch/ks (merge_ks): include 'mpif.h' -> use mpi_mod, ufront allocated on first use.
#ifdef LIGHT_MPI_COMM
#define RECV_IGRID(i) reception(icpu,ilevel)%pcomm%igrid(i)
#else
#define RECV_IGRID(i) reception(icpu,ilevel)%igrid(i)
#endif
module my_utils
  implicit none
  public ! default

CONTAINS

!----------------------------------------------------------------------------
! ensure_ufront: allocate ufront (same size as uold) the first time it is needed
!----------------------------------------------------------------------------
  subroutine ensure_ufront
    use amr_commons
    use hydro_commons
    implicit none
    integer::ncell
    if(.not.allocated(ufront))then
      ncell=ncoarse+twotondim*ngridmax
      allocate(ufront(0:ncell,1:nvar))
      ufront=0d0
      if(myid==1)write(*,'(A,F8.3,A)')' star_formation_fg: ufront allocated (',8d0*dble(ncell+1)*nvar/1024d0**3,' GB per rank)'
    end if
  end subroutine ensure_ufront


!----------------------------------------------------------------------------
! set_ufront
!----------------------------------------------------------------------------
! set hydro variables ufront (not necessarily synchronized between levels)
! - ufront: forefront hydro variables
!   * used in star_fomration.f90 and feedback.f90 to realize non-local gas depletion
!
! how to get ufront depend on the flag
!
! flag_option = 1 -> use uold
! flag_option = 2 -> use uold after uploading higher level cells to the
!                    lower level cells by calling upload_fine (averaging)
! flag_option = 3 -> use unew after reflecting the advection effect by calling
!                    virtual (cpu) boundary communications with make_virtual_reverse_dp
!----------------------------------------------------------------------------
  subroutine set_ufront(ilevel,flag_option)
    use amr_commons
    use hydro_commons
    use poisson_commons
    implicit none
    integer::ilevel
    integer::flag_option
    integer::i,ivar,ind,iskip,icpu

    call ensure_ufront

    !--------------------------------------------------
    ! get uold
    !--------------------------------------------------
    if (flag_option == 1) then
      !------ copy uold to ufront -----!
      do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do ivar=1,nvar
          !active grids
          do i=1,active(ilevel)%ngrid
            ufront(active(ilevel)%igrid(i)+iskip,ivar) = uold(active(ilevel)%igrid(i)+iskip,ivar)
          end do
          !virtual boundary cells
          do icpu=1,ncpu
            do i=1,reception(icpu,ilevel)%ngrid
              ufront(RECV_IGRID(i)+iskip,ivar)=uold(RECV_IGRID(i)+iskip,ivar)
            end do
          end do
        end do
      end do

      !--------------------------------------------------
      ! get uold and upload from higher to lower levels
      !--------------------------------------------------
    else if (flag_option == 2) then
      !------ copy uold to ufront -----!
      do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do ivar=1,nvar
          !active grids
          do i=1,active(ilevel)%ngrid
            ufront(active(ilevel)%igrid(i)+iskip,ivar) = uold(active(ilevel)%igrid(i)+iskip,ivar)
          end do
          !virtual boundary cells
          do icpu=1,ncpu
            do i=1,reception(icpu,ilevel)%ngrid
              ufront(RECV_IGRID(i)+iskip,ivar)=uold(RECV_IGRID(i)+iskip,ivar)
            end do
          end do
        end do
      end do
      !------ upload from ilevel+1 to ilevel -----!
      call upload_fine_ufront(ilevel)


      !--------------------------------------------------
      ! get unew and do virtual_boundary communication
      !--------------------------------------------------
    else if (flag_option == 3) then

      !------ copy unew to ufront -----!
      do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do ivar=1,nvar
          !active grids
          do i=1,active(ilevel)%ngrid
            ufront(active(ilevel)%igrid(i)+iskip,ivar) = unew(active(ilevel)%igrid(i)+iskip,ivar)
          end do
          !virtual boundary cells
          do icpu=1,ncpu
            do i=1,reception(icpu,ilevel)%ngrid
              ufront(RECV_IGRID(i)+iskip,ivar)=unew(RECV_IGRID(i)+iskip,ivar)
            end do
          end do
        end do
      end do
      !------ consider virtual boundary advection ------!
      do ivar=1,nvar
        call make_virtual_reverse_dp(ufront(1,ivar),ilevel)
      end do

      !--------------------------------------------------
      !flag_option must be either 1, 2, or 3
      !--------------------------------------------------
    else
      print *, "***ERROR *** flag_option must be either 1, 2, or 3, but is set to ", flag_option
      print *, "stopping..."
      stop
    end if


  end subroutine set_ufront

!##########################################################################
! copy the subroutine for uold from interpol_hydro.f90 and modify it for ufront
!##########################################################################
  subroutine upload_fine_ufront(ilevel)
    use amr_commons
    use hydro_commons
    implicit none
    integer::ilevel
    !----------------------------------------------------------------------
    ! This routine performs a restriction operation (averaging down)
    ! for the hydro variables.
    !----------------------------------------------------------------------
    integer::i,ncache,igrid,ngrid,ind,iskip,nsplit,icell
    integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_split
    logical,dimension(1:nvector),save::ok

    if(ilevel==nlevelmax)return
    if(numbtot(1,ilevel)==0)return
    if(verbose)write(*,111)ilevel

    ! Loop over active grids by vector sweeps
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

        ! Gather split cells
        do i=1,ngrid
          ok(i)=son(ind_cell(i))>0
        end do

        ! Count split cells
        nsplit=0
        do i=1,ngrid
          if(ok(i))nsplit=nsplit+1
        end do

        ! Upload for selected cells
        if(nsplit>0)then
          icell=0
          do i=1,ngrid
            if(ok(i))then
              icell=icell+1
              ind_split(icell)=ind_cell(i)
            end if
          end do
          call upl_ufront(ind_split,nsplit)
        end if

      end do
      ! End loop over cells

    end do
    ! End loop over grids

111 format('   Entering upload_fine for level',i2)

  end subroutine upload_fine_ufront
!##########################################################################
! copy the subroutine for uold from interpol_hydro.f90 and modify it for ufront
!##########################################################################
  subroutine upl_ufront(ind_cell,ncell)
    use amr_commons
    use hydro_commons
    implicit none
    integer::ncell
    integer,dimension(1:nvector)::ind_cell
    !---------------------------------------------------------------------
    ! This routine performs a restriction operation (averaging down)
    ! for the following variables:
    ! interpol_var=0: use rho, rho u and E
    ! interpol_tar=1: use rho, rho u and rho epsilon
    !---------------------------------------------------------------------
    integer ::ivar,irad,i,idim,ind_son,iskip_son
    integer ,dimension(1:nvector),save::igrid_son,ind_cell_son
    real(dp),dimension(1:nvector),save::getx,ekin,erad

    ! Get child oct index
    do i=1,ncell
      igrid_son(i)=son(ind_cell(i))
    end do

    !-------------------------------
    ! Average conservative variables
    !-------------------------------
    ! Loop over variables
    do ivar=1,nvar

      getx(1:ncell)=0.0d0
      do ind_son=1,twotondim
        iskip_son=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncell
          ind_cell_son(i)=iskip_son+igrid_son(i)
        end do
        ! Update average
        do i=1,ncell
          getx(i)=getx(i)+ufront(ind_cell_son(i),ivar)
        end do
      end do

      ! Scatter result to cells
      do i=1,ncell
        ufront(ind_cell(i),ivar)=getx(i)/dble(twotondim)
      end do

    end do
    ! End loop over variables

    !------------------------------------------------
    ! Average internal energy instead of total energy
    !------------------------------------------------
    if(interpol_var==1 .or. interpol_var==2)then

      getx(1:ncell)=0.0d0
      do ind_son=1,twotondim
        iskip_son=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncell
          ind_cell_son(i)=iskip_son+igrid_son(i)
        end do
        ! Compute child kinetic energy
        ekin(1:ncell)=0.0d0
        do idim=1,ndim
          do i=1,ncell
            ekin(i)=ekin(i)+0.5d0*ufront(ind_cell_son(i),1+idim)**2 &
            &               /max(ufront(ind_cell_son(i),1),smallr)
          end do
        end do
        ! Compute child radiative energy
        erad(1:ncell)=0.0d0
#if NENER>0
        do irad=1,nener
          do i=1,ncell
            erad(i)=erad(i)+ufront(ind_cell_son(i),ndim+2+irad)
          end do
        end do
#endif
        ! Update average
        do i=1,ncell
          getx(i)=getx(i)+ufront(ind_cell_son(i),ndim+2)-ekin(i)-erad(i)
        end do
      end do

      ! Compute new kinetic energy
      ekin(1:ncell)=0.0d0
      do idim=1,ndim
        do i=1,ncell
          ekin(i)=ekin(i)+0.5d0*ufront(ind_cell(i),1+idim)**2 &
          &               /max(ufront(ind_cell(i),1),smallr)
        end do
      end do
      ! Compute new radiative energy
      erad(1:ncell)=0.0d0
#if NENER>0
      do irad=1,nener
        do i=1,ncell
          erad(i)=erad(i)+ufront(ind_cell(i),ndim+2+irad)
        end do
      end do
#endif

      ! Scatter result to cells
      do i=1,ncell
        ufront(ind_cell(i),ndim+2)=getx(i)/dble(twotondim)+ekin(i)+erad(i)
      end do

    end if

  end subroutine upl_ufront


  !-----------------------------------------------------------
  ! Get cpu of the grid that contains the given position
  ! using cmp_cpumap function (copied from star_formation.f90, ADDED BY KS)
  !-----------------------------------------------------------
  subroutine get_cpu_of_pos(pos,ilevel,icpu)
    use amr_commons
    real(dp),dimension(1:3),intent(IN)::pos !given position
    integer,intent(IN)::ilevel              !given level
    integer,intent(OUT)::icpu               !cpu responsible for the grid that contatins point at pos at ilevel (no meaning if such a grid does not exist in the simulaiton)
    real(dp),dimension(1:nvector,1:3)::pos2
    integer,dimension(1:nvector)::cc
    real(dp) :: dx_cell, dx_grid
    integer :: ix_cell,iy_cell,iz_cell,ix_grid,iy_grid,iz_grid

    ! TO MODIFIED (merge_ks, 2026-11-14): pos and cmp_cpumap_lev are in user units ([0,boxlen)), so the grid
    ! size must be scaled by boxlen/nx_loc (ks has boxlen=1; with boxlen=100 the old snapping moved the point
    ! by up to one code-unit grid, i.e. into a neighbouring real grid when pos was within 1/boxlen of its edge)
    dx_cell = 0.5D0**ilevel*boxlen/dble(icoarse_max-icoarse_min+1)
    dx_grid = 2*dx_cell

    !integer position of cell containing the pos
    !nint: nearest integer
    ! ix_cell = nint(pos(1)/dx_cell+0.5d0)
    ! iy_cell = nint(pos(2)/dx_cell+0.5d0)
    ! iz_cell = nint(pos(3)/dx_cell+0.5d0)

    !integer position of grid containing the pos
    ix_grid = int(pos(1)/dx_grid)+1
    iy_grid = int(pos(2)/dx_grid)+1
    iz_grid = int(pos(3)/dx_grid)+1

    !real grid position
    pos2(1,1:3) = (/(ix_grid-0.5d0)*dx_grid, (iy_grid-0.5d0)*dx_grid, (iz_grid-0.5d0)*dx_grid/)

    ! print '(A,1P3E20.12)', "KS DEBUG (get_cpu_of_pos) pos =         ", pos(1:3)
    ! print '(A,1P3E20.12)', "KS DEBUG (get_cpu_of_pos) pos_of_grid = ", pos2(1,1:3)
    ! print '(A,3I8)', "KS DEBUG (get_cpu_of_pos) ipos =      ", ix_cell, iy_cell, iz_cell
    ! print '(A,3I8)', "KS DEBUG (get_cpu_of_pos) ipos_grid = ", ix_grid, iy_grid, iz_grid
    ! print *, "KS DEBUG STOP HERE"
    ! stop

    !compute cpumap for the grid position
    call cmp_cpumap_lev(pos2,cc,1,ilevel) !TO MODIFIED: two-tier: owner of the level-ilevel grid centred at pos2
    icpu = cc(1)
  end subroutine get_cpu_of_pos

  !-----------------------------------------------------------
  ! Get cell_index for given pos at lev
  ! Return cell_index = -1 if no active grid is found at this level on myid
  !-----------------------------------------------------------
  subroutine get_active_cell_index_this_level(pos,lev,cell_index)
    use amr_commons
    real(dp),dimension(1:3),intent(IN)::pos !given position
    integer,intent(IN)::lev              !given level
    integer,intent(OUT)::cell_index               !cpu responsible for the grid that contatins point at pos at ilevel (no meaning if such a grid does not exist in the simulaiton)
    integer :: cc               !cpu responsible for the grid that contatins point at pos at ilevel (no meaning if such a grid does not exist in the simulaiton)
    real(dp),dimension(1:nvector,1:3)::pos2
    integer,dimension(1:nvector):: cell_index2, cell_levl2
    integer :: nok

    !set vector variable
    nok = 1
    pos2(1,1:3) = pos(:)

    !take care of grid poisition when searching for cpu
    call get_cpu_of_pos(pos,lev,cc)
    call get_cell_index(cell_index2,cell_levl2,pos2,lev,nok) !get cell_index/cell_levl

    !check if active cell exists at lev on this cpu
    if (cc == myid .and. cell_levl2(1) == lev) then
      cell_index = cell_index2(1)
    else
      !set invalid value when active grid not found on this cpu at this level
      cell_index=-1
    end if
  end subroutine get_active_cell_index_this_level

  !-----------------------------------------------------------
  ! Get cell_index for given pos at highest level (corresponding to a leaf cell)
  ! Return cell_index = -1 if active grid is owned by another cpu
  !-----------------------------------------------------------
  subroutine get_active_cell_index_highest_level(pos,lev,cell_index)
    use amr_commons
    use mpi_mod
    real(dp),dimension(1:3),intent(IN)::pos !given position
    integer,intent(OUT)::lev              !highest level of an active cell that contains pos
    integer,intent(OUT)::cell_index               !cpu responsible for the grid that contatins point at pos at ilevel (no meaning if such a grid does not exist in the simulaiton)
    integer :: cc               !cpu responsible for the grid that contatins point at pos at ilevel (no meaning if such a grid does not exist in the simulaiton)
    real(dp),dimension(1:nvector,1:3)::pos2
    integer,dimension(1:nvector):: cell_index2, cell_levl2
    integer :: nok
    integer :: ifound
    integer :: info


    !set vector variable
    nok = 1
    pos2(1,1:3) = pos(:)

    !initialization
    ifound = 0
    cell_index=-1

    do lev = nlevelmax, levelmin, -1
      !take care of grid poisition when searching for cpu
      call get_cpu_of_pos(pos,lev,cc)
      call get_cell_index(cell_index2,cell_levl2,pos2,lev,nok) !get cell_index/cell_levl

      !check if active cell exists at lev on this cpu
      if (cc == myid .and. cell_levl2(1) == lev) then
        ifound = 1
        cell_index = cell_index2(1)
      end if

      !tell all cpu if the cell is found at lev
      call MPI_ALLREDUCE(MPI_IN_PLACE,ifound,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)

      !exit if the cell is found
      if (ifound == 1) then
        exit
      end if
    end do !lev
    ! print *, "(KS DEBUG) myid, lev, cell_index = ", myid, lev, cell_index
  end subroutine get_active_cell_index_highest_level



!#########################################################################
!#########################################################################


end module my_utils
