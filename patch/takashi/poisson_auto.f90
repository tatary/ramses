!################################################################
! Choice between multigrid and conjugate gradient for the Poisson
! solve at level ilevel (patch/takashi, 2026-11-11).
!
! cg_auto=.false.: level >= cg_levelmin uses CG, the others multigrid
!   (upstream behaviour, bit-identical to the old amr_step branch).
! cg_auto=.true. : cg_levelmin only gives the initial choice. Each
!   level > levelmin carries a flag cg_trial_due. It is set for all
!   levels at the start of the run and every cg_auto_every coarse
!   steps (not for a level whose previous trial is still waiting for
!   the selected method's time); a level that does not exist yet (no
!   grids: amr_step never calls us) simply keeps the flag until its
!   first solve, so a level that appears later is evaluated at once. When a level is solved
!   with the flag set, this solve uses the method that is NOT selected
!   ("trial") and the stored time of the selected method is discarded,
!   so the comparison uses the trial solve and the next solve with the
!   selected method (normally the next fine step of the same coarse
!   step). Both times are reduced with MPI_MIN (the last-arriving rank
!   measures the pure solve time, the others include waiting) and the
!   faster method is selected; a switch requires the other method to
!   be faster by the fraction cg_auto_margin (hysteresis). All ranks
!   decide on the same reduced numbers.
!   A CG solve that hits itermax disables CG for that level for good.
!################################################################
subroutine poisson_solve_level(ilevel,icount)
  use amr_commons
  use poisson_commons
  use mpi_mod
  implicit none
  integer,intent(in)::ilevel,icount
  integer::i
  logical::docg
  real(dp)::t0,t1
  real(dp)::wall_time

  ! Fixed choice (upstream behaviour)
  if(ilevel==levelmin.or.(.not.cg_auto))then
     if(ilevel>levelmin.and.ilevel>=cg_levelmin)then
        call phi_fine_cg(ilevel,icount)
     else
        call multigrid_fine(ilevel,icount)
     end if
     return
  end if

  ! Every cg_auto_every coarse steps: all levels are due for a trial,
  ! except those still waiting for the selected method's time after the
  ! last trial (a level solved once per coarse step, e.g. levelmin+1 with
  ! nsubcycle(levelmin)=1, gets that time only in the next coarse step).
  if(nstep_coarse/=cg_auto_step)then
     if(cg_auto_step>=0.and.mod(nstep_coarse,max(cg_auto_every,1))==0)then
        do i=levelmin+1,nlevelmax
           if(.not.cg_pending(i))cg_trial_due(i)=.true.
        end do
     end if
     cg_auto_step=nstep_coarse
  end if

  ! Which method for this solve
  docg=use_cg(ilevel)
  if(cg_trial_due(ilevel))then
     docg=.not.use_cg(ilevel)          ! try the other method once
     if(use_cg(ilevel))then            ! and measure the selected one afresh
        t_solve_cg(ilevel)=-1d0
     else
        t_solve_mg(ilevel)=-1d0
     end if
     cg_trial_due(ilevel)=.false.
     cg_pending(ilevel)=.true.
  end if
  if(cg_failed(ilevel))docg=.false.

  ! Solve and time it
  t0=wall_time()
  if(docg)then
     call phi_fine_cg(ilevel,icount)
  else
     call multigrid_fine(ilevel,icount)
  end if
  t1=wall_time()-t0
  if(docg)then
     t_solve_cg(ilevel)=t1
  else
     t_solve_mg(ilevel)=t1
  end if

  ! CG did not converge in this solve: multigrid from now on
  if(docg.and.cg_failed(ilevel))then
     use_cg(ilevel)=.false.
     cg_pending(ilevel)=.false.
     if(myid==1)write(*,'(" cg_auto: level",I3," CG did not converge, multigrid from now on")')ilevel
     return
  end if
  if(cg_failed(ilevel))cg_pending(ilevel)=.false.

  ! Decide as soon as both methods have been timed since the trial
  if(cg_pending(ilevel).and.t_solve_mg(ilevel)>=0d0.and.t_solve_cg(ilevel)>=0d0)then
     call cg_auto_decide(ilevel)
     cg_pending(ilevel)=.false.
  end if

end subroutine poisson_solve_level
!################################################################
subroutine cg_auto_decide(ilevel)
  use amr_commons
  use poisson_commons
  use mpi_mod
  implicit none
  integer,intent(in)::ilevel
  real(dp),dimension(2)::tsol_in,tsol_out
  real(dp)::tmg,tcg
  character(LEN=2)::choice
#ifndef WITHOUTMPI
  integer::info
#endif

  tsol_in(1)=t_solve_mg(ilevel); tsol_in(2)=t_solve_cg(ilevel)
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(tsol_in,tsol_out,2,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
#else
  tsol_out=tsol_in
#endif
  tmg=tsol_out(1); tcg=tsol_out(2)
  if(use_cg(ilevel))then
     if(tmg<(1d0-cg_auto_margin)*tcg)use_cg(ilevel)=.false.
  else
     if(tcg<(1d0-cg_auto_margin)*tmg)use_cg(ilevel)=.true.
  end if
  choice='MG'; if(use_cg(ilevel))choice='CG'
  if(myid==1)write(*,'(" cg_auto: level",I3,"  MG",F9.4," s  CG",F9.4," s  -> ",A2)')ilevel,tmg,tcg,choice

end subroutine cg_auto_decide
!################################################################
function wall_time()
  use amr_parameters, only: dp
  use mpi_mod
  implicit none
  real(dp)::wall_time
#ifndef WITHOUTMPI
  wall_time=MPI_WTIME()
#else
  integer(kind=8)::c,r
  call system_clock(c,r)
  wall_time=dble(c)/dble(r)
#endif
end function wall_time
