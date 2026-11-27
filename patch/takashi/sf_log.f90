!################################################################
! sf_log: per-rank star formation / feedback event log (sf_log_properties)
!
! One file per rank and per output interval, output_XXXXX/stars_XXXXX.outNNNNN,
! written directly by the rank that creates the particle (no collective
! communication). The columns are those of the upstream sf_log_properties
! (event, id, ilevel, mp, xp, vp, u1..u_nvar, tag) plus the particle family
! (merge_ks: the FG star formation needs to tell Pop II / Pop III / SFC apart).
! A restart rewinds the file from the output_restart snapshot (init_amr).
!
! Used by star_formation (sf_fg=.false.), star_formation_fg (sf_fg=.true.)
! and feedback (SN events, event=1, which it still writes itself).
!################################################################
module sf_log
  implicit none
contains
  !----------------------------------------------------------------------
  ! Open (create with a header, or append to) this rank's log file.
  ! Mirrors the former inline code of star_formation/feedback, including the
  ! IOGROUPSIZE token receive (there has never been a matching send; the
  ! production runs use IOGROUPSIZE=0).
  !----------------------------------------------------------------------
  subroutine sf_log_open(ilun)
    use amr_commons
    use hydro_commons, only: nvar
    use mpi_mod
    implicit none
    integer,intent(inout)::ilun
    character(LEN=80)::filename,filedir,fileloc,filedirini
    character(LEN=5)::nchar,ncharcpu
    logical::file_exist
    integer::idim,ivar
#ifndef WITHOUTMPI
    integer::info2,dummy_io
    integer,parameter::tag=1120
#endif

    if(.not.(sf_log_properties.and.ifout.gt.1))return
    call title(ifout-1,nchar)
    if(IOGROUPSIZEREP>0) then
      call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
      filedirini='output_'//TRIM(nchar)//'/'
      filedir='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/'
    else
      filedir='output_'//TRIM(nchar)//'/'
    endif
    filename=TRIM(filedir)//'stars_'//TRIM(nchar)//'.out'
    ilun=myid+103
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

    inquire(file=fileloc,exist=file_exist)
    if((.not.file_exist).or.(ifout-1.ne.ifout_stars_log)) then
      open(ilun, file=fileloc, status='replace', form='formatted')
      write(ilun,'(A24)',advance='no') '# event id  ilevel  mp  '
      do idim=1,ndim
        write(ilun,'(A2,I1,A2)',advance='no') 'xp',idim,'  '
      enddo
      do idim=1,ndim
        write(ilun,'(A2,I1,A2)',advance='no') 'vp',idim,'  '
      enddo
      do ivar=1,nvar
        if(ivar.ge.10) then
          write(ilun,'(A1,I2,A2)',advance='no') 'u',ivar,'  '
        else
          write(ilun,'(A1,I1,A2)',advance='no') 'u',ivar,'  '
        endif
      enddo
      write(ilun,'(A5)',advance='no') 'tag  '
      write(ilun,'(A8)',advance='no') 'family  '
      write(ilun,'(A1)') ' '
      ifout_stars_log = ifout-1
    else
      open(ilun, file=fileloc, status="old", position="append", action="write", form='formatted')
    endif
  end subroutine sf_log_open

  !----------------------------------------------------------------------
  subroutine sf_log_close(ilun)
    use amr_parameters, only: sf_log_properties
    implicit none
    integer,intent(in)::ilun
    if(sf_log_properties) close(ilun)
  end subroutine sf_log_close

  !----------------------------------------------------------------------
  ! One line (event 0 = star formation) for the new particle ipart born in
  ! cell icell. primitive=.true. when uold(icell,:) currently holds the
  ! primitive variables of star_formation (d, u, v, w, e_th/d, passive/d);
  ! .false. when it holds the conservative variables (the line is then
  ! written with the same meaning after converting here).
  !----------------------------------------------------------------------
  subroutine sf_log_write_star(ilun,ipart,icell,ilevel,primitive)
    use amr_commons
    use pm_commons
    use hydro_commons
    implicit none
    integer,intent(in)::ilun,ipart,icell,ilevel
    logical,intent(in)::primitive
    integer::idim,ivar
    real(dp)::uvar,d,ekin
    real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
#if NENER>0
    integer::irad
#endif

    if(.not.sf_log_properties)return
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    write(ilun,'(I10)',advance='no') 0
    write(ilun,'(2I10,E24.12)',advance='no') idp(ipart),ilevel,mp(ipart)
    do idim=1,ndim
      write(ilun,'(E24.12)',advance='no') xp(ipart,idim)
    enddo
    do idim=1,ndim
      write(ilun,'(E24.12)',advance='no') vp(ipart,idim)
    enddo
    d=uold(icell,1)
    write(ilun,'(E24.12)',advance='no') d
    if(primitive)then
      do ivar=2,nvar
        if(ivar.eq.ndim+2)then
          ! Temperature
          uvar=(gamma-1.0d0)*(uold(icell,ndim+2))*scale_T2
#ifdef SOLVERmhd
        else if(ivar.ge.ndim+3.and.ivar.le.ndim+2+ndim)then
          ! Average of left and right magnetic fields (Mizuki Ono, 2026/06/13)
          uvar=0.5d0*(uold(icell,ivar)+ &
          & uold(icell,nvar+ivar-(ndim+2)))
#endif
        else
          uvar=uold(icell,ivar)
        endif
        write(ilun,'(E24.12)',advance='no') uvar
      enddo
    else
      ! conservative -> the same quantities as above
      ekin=0d0
      do idim=1,ndim
        ekin=ekin+0.5d0*uold(icell,idim+1)**2/max(d,smallr)
      enddo
#if NENER>0
      do irad=0,nener-1
        ekin=ekin+uold(icell,inener+irad)
      enddo
#endif
#ifdef SOLVERmhd
      do idim=1,ndim
        ekin=ekin+0.125d0*(uold(icell,idim+ndim+2)+uold(icell,idim+nvar))**2
      enddo
#endif
      do ivar=2,nvar
        if(ivar.le.ndim+1)then
          uvar=uold(icell,ivar)/max(d,smallr)
        else if(ivar.eq.ndim+2)then
          uvar=(gamma-1.0d0)*(uold(icell,ndim+2)-ekin)/max(d,smallr)*scale_T2
#ifdef SOLVERmhd
        else if(ivar.ge.ndim+3.and.ivar.le.ndim+2+ndim)then
          uvar=0.5d0*(uold(icell,ivar)+uold(icell,nvar+ivar-(ndim+2)))
#endif
        else
          uvar=uold(icell,ivar)/max(d,smallr)
        endif
        write(ilun,'(E24.12)',advance='no') uvar
      enddo
    endif
    write(ilun,'(I10)',advance='no') typep(ipart)%tag
    write(ilun,'(I10)',advance='no') typep(ipart)%family
    write(ilun,'(A1)') ' '
  end subroutine sf_log_write_star
end module sf_log
