!################################################################
! output_analysis: STUB (merge_ks).
!
! patch/ks has a 3600-line module of the same name that writes nested
! uniform-grid VTK files ("ugvtk") around new SFCs / Pop III stars / sinks,
! and at given ages. Whether to port it is still open (docs/merge_ks.md,
! TODO 3b); the FG star formation (star_formation_fg.f90) only needs the
! entry points below to exist. They do nothing except warn once.
!################################################################
module output_analysis
  implicit none
  logical,save,private :: warned=.false.
contains
  subroutine warn_once(what)
    use amr_commons, only: myid
    implicit none
    character(len=*),intent(in)::what
    if(.not.warned)then
      if(myid==1)write(*,'(A)')' output_analysis: '//trim(what)//' requested but the ugvtk output of patch/ks is not ported (stub); skipped'
      warned=.true.
    end if
  end subroutine warn_once

  subroutine output_around_SFC(iSFC, pos, suffix, use_current_ufront)
    use amr_parameters, only: dp
    implicit none
    integer, intent(IN) :: iSFC
    real(dp), dimension(1:3), intent(IN) :: pos
    character(len=*), intent(IN), optional :: suffix
    logical, intent(IN), optional :: use_current_ufront
    call warn_once('output_around_SFC')
  end subroutine output_around_SFC

  subroutine output_around_star(prefix, pos, use_current_ufront)
    use amr_parameters, only: dp
    implicit none
    character(len=*), intent(IN) :: prefix
    real(dp), dimension(1:3), intent(IN) :: pos
    logical, intent(IN), optional :: use_current_ufront
    call warn_once('output_around_star')
  end subroutine output_around_star

  subroutine output_around_SFC_at_given_age(dt)
    use amr_parameters, only: dp
    implicit none
    real(dp), intent(IN) :: dt
    call warn_once('output_around_SFC_at_given_age')
  end subroutine output_around_SFC_at_given_age

  subroutine output_around_star_at_given_age(dt)
    use amr_parameters, only: dp
    implicit none
    real(dp), intent(IN) :: dt
    call warn_once('output_around_star_at_given_age')
  end subroutine output_around_star_at_given_age

  subroutine systemcall_wrapper(command, status)
    implicit none
    character(len=*),intent(IN) :: command
    integer,intent(OUT),optional :: status
    integer :: cmdstat
    call execute_command_line(TRIM(command), exitstat=cmdstat)
    if(present(status)) status=cmdstat
  end subroutine systemcall_wrapper
end module output_analysis
