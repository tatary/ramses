!=======================================================================
real(kind=8) function wallclock()
  use mpi_mod
  implicit none
#ifdef WITHOUTMPI
  integer,      save :: tstart
  integer            :: tcur
  integer            :: count_rate
#else
  real(kind=8), save :: tstart
  real(kind=8)       :: tcur
#endif
  logical,      save :: first_call=.true.
  real(kind=8), save :: norm, offset=0
  !---------------------------------------------------------------------
  if (first_call) then
#ifdef WITHOUTMPI
    call system_clock(count=tstart, count_rate=count_rate)
    norm=1d0/count_rate
#else
    norm = 1d0
    tstart = MPI_Wtime()
#endif
    first_call=.false.
  end if
#ifdef WITHOUTMPI
  call system_clock(count=tcur)
#else
  tcur = MPI_Wtime()
#endif
  wallclock = (tcur-tstart)*norm + offset
  if (wallclock < 0.) then
    offset = offset + 24d0*3600d0
    wallclock = wallclock + 24d0*3600d0
  end if
end function wallclock
!=======================================================================
module timer_m
  implicit none
  integer,            parameter         :: mtimer=200                    ! max nr of timers
  real(kind=8),       dimension(mtimer) :: start, time
  integer                               :: ntimer=0, itimer
  character(len=72), dimension(mtimer)  :: labels
contains
!-----------------------------------------------------------------------
  subroutine findit (label)
    implicit none
    character(len=*) label
    do itimer=1,ntimer
      if (trim(label) == trim(labels(itimer))) return
    end do
    ntimer = ntimer+1
    itimer = ntimer
    labels(itimer) = label
    time(itimer) = 0
  end subroutine findit
end module timer_m
!=======================================================================
subroutine timer (label, cmd)
  use timer_m
  implicit none
  character(len=*)::label,cmd
  real(kind=8)::wallclock,current
!-----------------------------------------------------------------------
  current = wallclock()                                                 ! current time
  if (itimer > 0) then                                                  ! if timer is active ..
    time(itimer) = time(itimer) + current - start(itimer)              ! add to it
  end if
  call findit (label)                                                   ! locate timer slot
  if (cmd == 'start') then                                              ! start command
    start(itimer) = current                                            ! register start time
  else if (cmd == 'stop') then                                          ! stop command
    itimer = 0                                                         ! turn off timer
  end if
end subroutine timer
!=======================================================================
subroutine gather_hostnames (names)
  ! Host (processor) name of every rank, gathered on rank 1 (names is only
  ! meaningful there). Used for the startup log and the timer_ranks header
  ! so that a slow rank can be traced to the node it ran on.
  use amr_commons
  use mpi_mod
  implicit none
  character(len=64),dimension(1:ncpu)::names
  character(len=64)::myname
#ifndef WITHOUTMPI
  character(len=MPI_MAX_PROCESSOR_NAME)::pname
  integer::plen,mpi_err
  call MPI_GET_PROCESSOR_NAME(pname,plen,mpi_err)
  myname=pname(1:min(plen,64))
  call MPI_GATHER(myname,64,MPI_CHARACTER,names,64,MPI_CHARACTER,0,MPI_COMM_WORLD,mpi_err)
#else
  myname='localhost'
  names(1)=myname
#endif
end subroutine gather_hostnames
!=======================================================================
subroutine output_hosts
  ! Startup log: which node each range of consecutive ranks runs on
  use amr_commons
  implicit none
  character(len=64),dimension(1:ncpu)::names
  integer::i,i0
  call gather_hostnames(names)
  if(myid/=1)return
  i0=1
  do i=2,ncpu+1
    if(i>ncpu)then
      write(*,'(a,i6,a,i6,a,a)')' Ranks ',i0,' -',i-1,' on ',trim(names(i0))
    else if(names(i)/=names(i0))then
      write(*,'(a,i6,a,i6,a,a)')' Ranks ',i0,' -',i-1,' on ',trim(names(i0))
      i0=i
    end if
  end do
end subroutine output_hosts
!=======================================================================
subroutine cost_timer_total (tsum)
  ! Sum of the compute timers used by cost_time_weighting (load_balance.f90):
  ! the parts whose per-rank time was found to vary at fixed cell count
  ! (patch/takashi/TODO.md 2b). Completed intervals only.
  use timer_m
  implicit none
  real(kind=8)::tsum
  integer::i,j
  character(len=24),dimension(3),parameter::cost_labels=(/'hydro - godunov','cooling        ','rt - godunov   '/)
  tsum=0d0
  do i=1,ntimer
    do j=1,size(cost_labels)
      if(trim(labels(i))==trim(cost_labels(j)))tsum=tsum+time(i)
    end do
  end do
end subroutine cost_timer_total
!=======================================================================
subroutine output_timer(write_file, filename)
  use amr_parameters
  use amr_commons
  use timer_m
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  real(kind=8) :: gtotal, avtime, rmstime
  real(kind=8), dimension(ncpu) :: vtime
  integer,      dimension(ncpu) :: all_ntimer
  logical,      dimension(ncpu) :: gprint_timer
  integer      :: imn, imx, mpi_err, icpu
  logical      :: print_timer
#endif
  real(kind=8) :: total
  integer      :: i
  logical      :: id_is_one, write_file
  integer      :: ilun=11
  character(LEN=80)::filename, fileloc !Optional for writing timing info
  character(LEN=80)::fileloc2          ! per-rank dump of all timers and grid counts
!-----------------------------------------------------------------------
  id_is_one = myid == 1
  total = 1e-9
  if (.not. write_file) ilun=6 ! 6 = std output
  if (id_is_one .and. write_file) then
    fileloc=TRIM(filename) ! Open file for timing info
    open(unit=ilun,file=fileloc,form='formatted')
  endif

  if (id_is_one .and. ncpu==1) write (ilun,'(/a,i7,a)') '     seconds         %    STEP (rank=',myid,')'
  do i = 1,ntimer
    total = total + time(i)
  end do
  if (ncpu==1) then
    do i = 1,ntimer
      if (id_is_one .and. time(i)/total > 0.001) write (ilun,'(f12.3,4x,f6.1,4x,a24)') &
        time(i), 100.*time(i)/total,labels(i)
    end do
    if (id_is_one) write (ilun,'(f12.3,4x,f6.1,4x,a)') total, 100., 'TOTAL'
  end if
#ifndef WITHOUTMPI
  if (ncpu > 1) then
    ! Check that timers are consistent across ranks
    call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
    call MPI_GATHER(ntimer,1,MPI_INTEGER,all_ntimer,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
    if (id_is_one) then
      if (maxval(all_ntimer) .ne. minval(all_ntimer)) then
        write (ilun,*)
        write (ilun,*) '--------------------------------------------------------------------'
        write (ilun,*) 'Error: Inconsistent number of timers on each rank. Min, max nr:', minval(all_ntimer), maxval(all_ntimer)
        write (ilun,*) 'Timing summary below can be misleading'
        write (ilun,*) 'Labels of timer on rank==1 :'
        write (ilun,*) '--------------------------------------------------------------------'
        do i=1,ntimer
          write(ilun,'(i3,1x,a)') i, labels(i)
        enddo
      endif
      ! Find first occurence of a rank with a different number of timers -- if it exists
      gprint_timer=.false.
      do icpu=1,ncpu
        if (all_ntimer(icpu) .ne. ntimer) then
          gprint_timer(icpu) = .true.
          exit
        endif
      enddo
      if (any(gprint_timer)) call sleep(1) ! Make sure that master rank finished, before we print from other rank.
    endif
    call MPI_SCATTER(gprint_timer,1,MPI_LOGICAL,print_timer,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_err)
    if (print_timer) then
      write (ilun,*)
      write (ilun,*) 'Labels of timer on rank==',myid
      write (ilun,*) '--------------------------------------------------------------------'
      do i=1,ntimer
        write(ilun,'(i3,1x,a)') i, labels(i)
      enddo
      write (ilun,*)
    endif

    call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
    call MPI_ALLREDUCE(total,gtotal,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
    gtotal = gtotal / ncpu

    if (id_is_one) write (ilun,*) '--------------------------------------------------------------------'
    if (id_is_one) write (ilun,'(/a)') '     minimum       average       maximum' // &
      '  standard dev        std/av       %   rmn   rmx  TIMER'
    ! Per-rank dump of all timers and grid counts: <timing file>.ranks when writing
    ! a timing file (outputs), else timer_ranks_<job id>.dat (end of run)
    if (write_file) then
      fileloc2=TRIM(filename)//'.ranks'
    else
      call timer_ranks_filename(0,fileloc2)
    end if
    do i = 1,ntimer
      call MPI_GATHER(real(time(i),kind=8),1,MPI_REAL8,vtime,1,MPI_REAL8,0,MPI_COMM_WORLD,mpi_err)
      if (id_is_one) then
        if (maxval(vtime)/gtotal > 0.001) then
          avtime  = sum(vtime) / ncpu ! average time used
          imn     = minloc(vtime,1)
          imx     = maxloc(vtime,1)
          rmstime = sqrt(sum((vtime - avtime)**2)/ncpu)
          write (ilun,'(5(f12.3,2x),f6.1,2x,2i4,4x,a24)') &
            vtime(imn), avtime, vtime(imx), rmstime, rmstime/avtime, 100.*avtime/gtotal, imn, imx, labels(i)
        endif
      endif
    end do
    if (id_is_one) write (ilun,'(f12.3,4x,f6.1,4x,a)') total, 100., 'TOTAL'
    call output_timer_ranks(fileloc2)
  endif
#endif
  if (id_is_one) close(ilun)
end subroutine output_timer
!=======================================================================
subroutine timer_ranks_filename(istep,fileloc2)
  ! timer_ranks_<job id>[_n<istep>].dat; the job id is PBS_JOBID (digits before
  ! the first dot) or SLURM_JOB_ID, else timer_ranks_n<nstep_coarse>.dat
  use amr_commons, only: nstep_coarse
  implicit none
  integer::istep
  character(LEN=80)::fileloc2
  character(LEN=32)::jobid
  integer::idot,lenv
  call get_environment_variable('PBS_JOBID',jobid,lenv)
  if(lenv<=0)call get_environment_variable('SLURM_JOB_ID',jobid,lenv)
  if(lenv>0)then
     idot=index(jobid,'.')
     if(idot>1)jobid=jobid(1:idot-1)
     if(istep>0)then
        write(fileloc2,'(a,a,a,i0,a)')'timer_ranks_',TRIM(jobid),'_n',istep,'.dat'
     else
        fileloc2='timer_ranks_'//TRIM(jobid)//'.dat'
     end if
  else
     write(fileloc2,'(a,i0,a)')'timer_ranks_n',max(istep,nstep_coarse),'.dat'
  end if
end subroutine timer_ranks_filename
!=======================================================================
subroutine output_timer_ranks(fileloc2)
  ! Per-rank dump (rank 1 writes): one line per timer (cumulative seconds on
  ! rank 1..ncpu), then per level the active, ghost and received-at-last-load-
  ! balance grids, the particles, the cost_time_weighting factor and the stars
  ! on each rank. Lines starting with # are the header (ncpu, host of each rank).
  ! Called at outputs and at the end of the run (output_timer), and every
  ! timer_ranks_every coarse steps (adaptive_loop) so that the difference of
  ! consecutive files gives the per-step time of each rank.
  use amr_parameters
  use amr_commons
  use pm_parameters, only: npartmax
  use pm_commons, only: npart, nextp, typep, is_star, headp_free, numbp_free
  use timer_m
  use mpi_mod
  implicit none
  character(LEN=80)::fileloc2
#ifndef WITHOUTMPI
  integer      :: ipart, jpart, nstar_loc, nghost_loc
  logical, allocatable, dimension(:) :: inuse
  character(len=64), dimension(ncpu) :: hostnames
  real(kind=8), dimension(ncpu) :: vtime
  integer      :: mpi_err, icpu
  integer      :: i, ilevel
  logical      :: id_is_one
  integer      :: ilun2=12
  integer,      dimension(ncpu) :: all_numbl
  id_is_one = myid == 1
  if (id_is_one) then
    open(unit=ilun2,file=fileloc2,form='formatted')
    write (ilun2,'(a,i6,a)') '# ncpu=',ncpu,'  columns: label, value on rank 1..ncpu'
  end if
  call gather_hostnames(hostnames)
  if (id_is_one) write (ilun2,'(a,100000(1x,a))') '# host', (trim(hostnames(i)),i=1,ncpu)
  do i = 1,ntimer
    call MPI_GATHER(real(time(i),kind=8),1,MPI_REAL8,vtime,1,MPI_REAL8,0,MPI_COMM_WORLD,mpi_err)
    if (id_is_one) write (ilun2,'(a24,100000(1x,ES12.5))') labels(i), vtime
  end do
  do ilevel=levelmin,nlevelmax
    call MPI_GATHER(numbl(myid,ilevel),1,MPI_INTEGER,all_numbl,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
    if (id_is_one) write (ilun2,'(a,i2,a,100000(1x,i12))') 'grids level ',ilevel,'          ', all_numbl
  end do
  ! Ghost (virtual) grids per level
  do ilevel=levelmin,nlevelmax
    nghost_loc=0
    do icpu=1,ncpu
      nghost_loc=nghost_loc+reception(icpu,ilevel)%ngrid
    end do
    call MPI_GATHER(nghost_loc,1,MPI_INTEGER,all_numbl,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
    if (id_is_one) write (ilun2,'(a,i2,a,100000(1x,i12))') 'ghost level ',ilevel,'          ', all_numbl
  end do
  ! Grids received (owner changed to this rank) at the last load balance, per level, and over the run
  do ilevel=levelmin,nlevelmax
    nghost_loc=0
    if(allocated(nrecv_lb))nghost_loc=nrecv_lb(ilevel)
    call MPI_GATHER(nghost_loc,1,MPI_INTEGER,all_numbl,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
    if (id_is_one) write (ilun2,'(a,i2,a,100000(1x,i12))') 'recv  level ',ilevel,'          ', all_numbl
  end do
  call MPI_GATHER(nrecv_lb_tot,1,MPI_INTEGER,all_numbl,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
  if (id_is_one) write (ilun2,'(a24,100000(1x,i12))') 'recv total', all_numbl
  call MPI_GATHER(npart,1,MPI_INTEGER,all_numbl,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
  if (id_is_one) write (ilun2,'(a24,100000(1x,i12))') 'particles', all_numbl
  call MPI_GATHER(cost_rank_fac,1,MPI_DOUBLE_PRECISION,vtime,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpi_err)
  if (id_is_one) write (ilun2,'(a24,100000(1x,es12.5))') 'cost factor', vtime
  ! Stars: count over the particle slots in use (all slots minus the free
  ! list). The grid linked lists cannot be used here: between merge_tree_fine
  ! and the next kill_tree_fine the lists of the children grids are stale.
  nstar_loc=0
  if(star.and.npartmax>0)then
    allocate(inuse(1:npartmax))
    inuse=.true.
    ipart=headp_free
    do jpart=1,numbp_free
      inuse(ipart)=.false.
      ipart=nextp(ipart)
    end do
    do ipart=1,npartmax
      if(inuse(ipart))then
        if(is_star(typep(ipart)))nstar_loc=nstar_loc+1
      end if
    end do
    deallocate(inuse)
  end if
  call MPI_GATHER(nstar_loc,1,MPI_INTEGER,all_numbl,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
  if (id_is_one) write (ilun2,'(a24,100000(1x,i12))') 'stars', all_numbl
  call MPI_GATHER(nsn_rank,1,MPI_INTEGER,all_numbl,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
  if (id_is_one) write (ilun2,'(a24,100000(1x,i12))') 'SN events', all_numbl
  if (id_is_one) close(ilun2)
#endif
end subroutine output_timer_ranks
!=======================================================================
subroutine reset_timer
  use timer_m
  use mpi_mod
  implicit none

!-----------------------------------------------------------------------
  do itimer = 1,ntimer
    time(itimer)=0
  end do
end subroutine reset_timer
!=======================================================================
subroutine update_time(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module
#if USE_TURB==1
  use turb_commons
#endif
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  real(kind=8)::ttend
  real(kind=8),save::ttstart=0
#endif
  integer::ilevel

  real(dp)::dt,econs,mcons
#ifdef SOLVERmhd
  real(dp)::sqrt_aexp_prev
#endif
#if USE_TURB==1
  real(kind=dp) :: cur_turb_rms
#endif
  integer::i,itest
  real(dp)::rho_max_loc,nH_max
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  integer::ilev,ind,igrid,icell,iskip
#ifndef WITHOUTMPI
  real(kind=8)::rho_max_all
  integer::info
#endif

  ! Local constants
  dt=dtnew(ilevel)
  itest=0

#ifndef WITHOUTMPI
  if(myid==1)then
    if(ttstart.eq.0.0)ttstart=MPI_WTIME()
  endif
#endif

  !-------------------------------------------------------------
  ! At this point, IF nstep_coarse has JUST changed, all levels
  ! are synchronised, and all new refinements have been done.
  !-------------------------------------------------------------
  if(nstep_coarse .ne. nstep_coarse_old)then

    !--------------------------
    ! Check mass conservation
    !--------------------------
    if(mass_tot_0==0.0D0)then
      mass_tot_0=mass_tot
      mcons=0.0D0
    else
      mcons=(mass_tot-mass_tot_0)/mass_tot_0
    end if

    !----------------------------
    ! Check energy conservation
    !----------------------------
    if(epot_tot_old.ne.0)then
      epot_tot_int=epot_tot_int + &
      & 0.5D0*(epot_tot_old+epot_tot)*log(aexp/aexp_old)
    end if
    epot_tot_old=epot_tot
    aexp_old=aexp
    if(einit==0.0D0)then
      einit=epot_tot+ekin_tot  ! initial total energy
      econs=0.0D0
    else
      econs=(ekin_tot+epot_tot-epot_tot_int-einit) / &
      &(-(epot_tot-epot_tot_int-einit)+ekin_tot)
    end if

    if(mod(nstep_coarse,ncontrol)==0.or.output_done)then

      !-------------------------------------------------------
      ! Compute maximum gas density over all leaf cells
      ! and convert it to hydrogen number density [cm^-3]
      !-------------------------------------------------------
      rho_max_loc=0.0D0
      if(hydro)then
        do ilev=levelmin,nlevelmax
          do ind=1,twotondim
            iskip=ncoarse+(ind-1)*ngridmax
            do igrid=1,active(ilev)%ngrid
              icell=active(ilev)%igrid(igrid)+iskip
              if(son(icell)==0)then
                rho_max_loc=max(rho_max_loc,uold(icell,1))
              end if
            end do
          end do
        end do
      end if
#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(rho_max_loc,rho_max_all,1,MPI_DOUBLE_PRECISION, &
      & MPI_MAX,MPI_COMM_WORLD,info)
      rho_max_loc=rho_max_all
#endif
      call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
      nH_max=rho_max_loc*scale_nH

      if(myid==1)then

        !-------------------------------
        ! Output AMR structure to screen
        !-------------------------------
        write(*,*)'Mesh structure'
        do i=1,nlevelmax
          if(numbtot(1,i)>0)write(*,999)i,numbtot(1:4,i)
        end do

        !----------------------------------------------
        ! Output mass and energy conservation to screen
        !----------------------------------------------
        if(cooling.or.pressure_fix)then
          write(*,778)nstep_coarse,mcons,econs,epot_tot,ekin_tot,eint_tot,nH_max
        else
          write(*,777)nstep_coarse,mcons,econs,epot_tot,ekin_tot,nH_max
        end if
#ifdef SOLVERmhd
        write(*,'(" emag=",ES9.2)') emag_tot
#endif
        if(pic)then
          write(*,888)nstep,t,dt,aexp,&
          & real(100.0D0*dble(used_mem_tot)/dble(ngridmax+1)),&
          & real(100.0D0*dble(npartmax-numbp_free_tot)/dble(npartmax+1))
        else
          write(*,888)nstep,t,dt,aexp,&
          & real(100.0D0*dble(used_mem_tot)/dble(ngridmax+1))
        endif
        itest=1
      end if
      output_done=.false.
    end if

    !---------------
    ! Exit program
    !---------------
    if(t>=tout(noutput).or.aexp>=aout(noutput).or. &
    & nstep_coarse>=nstepmax)then
      if(myid==1)then
        if(.not.finish_run)write(*,*)'Run completed'
#ifndef WITHOUTMPI
        ttend=MPI_WTIME()
        write(*,*)'Total elapsed time:',ttend-ttstart
#endif
      endif
      call clean_end
    end if

  end if
  nstep_coarse_old=nstep_coarse

  !----------------------------
  ! Output controls to screen
  !----------------------------
  if(mod(nstep,ncontrol)==0)then
    if(myid==1.and.itest==0)then
      if(pic)then
        write(*,888)nstep,t,dt,aexp,&
        & real(100.0D0*dble(used_mem_tot)/dble(ngridmax+1)),&
        & real(100.0D0*dble(npartmax-numbp_free_tot)/dble(npartmax+1))
      else
        write(*,888)nstep,t,dt,aexp,&
        & real(100.0D0*dble(used_mem_tot)/dble(ngridmax+1))
      endif
    end if
  end if

  !------------------------
  ! Update time variables
  !------------------------
  t=t+dt
  nstep=nstep+1
  if(cosmo)then
#ifdef SOLVERmhd
    ! Keep for magnetic field expansion
    sqrt_aexp_prev = SQRT(aexp)
#endif
    ! Find neighboring times
    i=1
    do while(tau_frw(i)>t.and.i<n_frw)
      i=i+1
    end do
    ! Interpolate expansion factor
    aexp = aexp_frw(i  )*(t-tau_frw(i-1))/(tau_frw(i  )-tau_frw(i-1))+ &
    & aexp_frw(i-1)*(t-tau_frw(i  ))/(tau_frw(i-1)-tau_frw(i  ))
    hexp = hexp_frw(i  )*(t-tau_frw(i-1))/(tau_frw(i  )-tau_frw(i-1))+ &
    & hexp_frw(i-1)*(t-tau_frw(i  ))/(tau_frw(i-1)-tau_frw(i  ))
    texp =    t_frw(i  )*(t-tau_frw(i-1))/(tau_frw(i  )-tau_frw(i-1))+ &
    &    t_frw(i-1)*(t-tau_frw(i  ))/(tau_frw(i-1)-tau_frw(i  ))

#ifdef SOLVERmhd
    do i=1,ilevel
      call update_cosmomag(i,SQRT(aexp)/sqrt_aexp_prev)
    end do
#endif
  else
    aexp = 1
    hexp = 0
    texp = t
  end if

#if USE_TURB==1
  if (turb) then
    call turb_check_time
    if (myid==1) then
      call current_turb_rms(cur_turb_rms)
      write (6,*) ' Current turbulent rms: ', cur_turb_rms
    end if
  end if
#endif

777 format(' Main step=',i7,' mcons=',1pe9.2,' econs=',1pe9.2, &
  & ' epot=',1pe9.2,' ekin=',1pe9.2,' nHmax=',1pe9.2)
778 format(' Main step=',i7,' mcons=',1pe9.2,' econs=',1pe9.2, &
  & ' epot=',1pe9.2,' ekin=',1pe9.2,' eint=',1pe9.2,' nHmax=',1pe9.2)
888 format(' Fine step=',i7,' t=',1pe12.5,' dt=',1pe10.3, &
  & ' a=',1pe10.3,' mem=',0pF4.1,'% ',0pF4.1,'%')
999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

end subroutine update_time

!------------------------------------------------------------------------
SUBROUTINE getProperTime(tau,tproper)
! Calculate proper time tproper corresponding to conformal time tau (both
! in code units).
!------------------------------------------------------------------------
  use amr_commons
  implicit none
  real(dp)::tau, tproper
  integer::i
  if(.not. cosmo .or. tau .eq. 0d0) then ! this might happen quite often
    tproper = tau
    return
  endif
  i = 1
  do while( tau_frw(i) > tau .and. i < n_frw )
    i = i+1
  end do
  tproper = t_frw(i  )*(tau-tau_frw(i-1))/(tau_frw(i  )-tau_frw(i-1))+ &
  & t_frw(i-1)*(tau-tau_frw(i  ))/(tau_frw(i-1)-tau_frw(i  ))
end subroutine getProperTime
!------------------------------------------------------------------------
SUBROUTINE getAgeGyr(t_birth_proper, age)
! Calculate proper time passed, in Gyrs, since proper time t_birth_proper
! (given in code units) until the current time.
!------------------------------------------------------------------------
  use amr_commons
  use pm_commons
  use constants,only: Gyr2sec
  implicit none
  real(dp):: t_birth_proper, age
  real(dp),save:: scale_t_Gyr
  logical,save::scale_init=.false.
  real(dp):: scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
  if( .not. scale_init) then
    ! The timescale has not been initialized
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    scale_t_Gyr = (scale_t/aexp**2)/Gyr2sec
    scale_init=.true.
  endif
  age = (texp - t_birth_proper) * scale_t_Gyr
end subroutine getAgeGyr
!------------------------------------------------------------------------
SUBROUTINE getAgeSec(t_birth_proper, age)
! Calculate proper time passed, in sec, since proper time t_birth_proper
! (given in code units) until the current time.
!------------------------------------------------------------------------
  use amr_commons
  use pm_commons
  implicit none
  real(dp):: t_birth_proper, age
  real(dp),save:: scale_t_sec
  logical::scale_init=.false.
  real(dp):: scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
  if( .not. scale_init) then
    ! The timescale has not been initialized
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    scale_t_sec = (scale_t/aexp**2)
    scale_init=.true.
  endif
  age = (texp - t_birth_proper) * scale_t_sec
end subroutine getAgeSec
!------------------------------------------------------------------------
