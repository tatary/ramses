!################################################################
! Choice between multigrid and conjugate gradient for the Poisson
! solve at level ilevel (patch/takashi, 2026-11-11).
!
! cg_auto=.false.: level >= cg_levelmin uses CG, the others multigrid
!   (upstream behaviour, bit-identical to the old amr_step branch).
! cg_auto=.true. : cg_levelmin only gives the initial choice. Each
!   level > levelmin carries a flag cg_trial_due. It is set for all
!   levels at the start of the run and every cg_auto_every coarse
!   steps; a level that does not exist yet (no grids: amr_step never
!   calls us) simply keeps the flag until its first solve, so a level
!   that appears later is evaluated at once. When a level is solved
!   with the flag set, it is solved twice: first with the method that
!   is NOT selected ("trial"), then again with the selected method.
!   The second solve starts from the same state (the initial guess is
!   made from phi_old or from level-1, rho is not touched), so the
!   solution that is kept is the one a run without trials would have.
!   (Keeping the trial solution instead put a solver-tolerance jump
!   between phi_old and phi of the level; finer levels interpolate
!   their boundary between the two in time, and with cg_warm_start
!   their iteration counts stayed high until the level was solved
!   again: TODO.md 2d-4e.)
!   Both times are reduced with MPI_MIN (the last-arriving rank
!   measures the pure solve time, the others include waiting) and the
!   faster method is selected for the following solves; a switch
!   requires the other method to be faster by the fraction
!   cg_auto_margin (hysteresis). All ranks decide on the same reduced
!   numbers.
!   A CG solve that hits itermax disables CG for that level for good.
!################################################################
subroutine poisson_solve_level(ilevel,icount)
  use amr_commons
  use poisson_commons
  use mpi_mod
  implicit none
  integer,intent(in)::ilevel,icount
  logical::docg,trial
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

  ! Every cg_auto_every coarse steps: all levels are due for a trial
  if(nstep_coarse/=cg_auto_step)then
     if(cg_auto_step>=0.and.mod(nstep_coarse,max(cg_auto_every,1))==0)then
        cg_trial_due(levelmin+1:nlevelmax)=.true.
     end if
     cg_auto_step=nstep_coarse
  end if

  ! Trial: solve once with the method that is not selected, only to time it
  trial=cg_trial_due(ilevel).and.(.not.cg_failed(ilevel))
  cg_trial_due(ilevel)=.false.
  if(trial)then
     t0=wall_time()
     if(use_cg(ilevel))then
        call multigrid_fine(ilevel,icount)
        t_solve_mg(ilevel)=wall_time()-t0
     else
        call phi_fine_cg(ilevel,icount)
        t_solve_cg(ilevel)=wall_time()-t0
        ! CG did not converge: no comparison, multigrid stays
        if(cg_failed(ilevel))then
           trial=.false.
           if(myid==1)write(*,'(" cg_auto: level",I3," CG did not converge, multigrid from now on")')ilevel
        end if
     end if
  end if

  ! Solve with the selected method (this solution is kept) and time it
  docg=use_cg(ilevel).and.(.not.cg_failed(ilevel))
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
     if(myid==1)write(*,'(" cg_auto: level",I3," CG did not converge, multigrid from now on")')ilevel
     return
  end if

  ! Both methods have just been timed on the same problem: choose
  if(trial)call cg_auto_decide(ilevel)

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
