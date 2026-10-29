
subroutine clean_end
  !---------------------------
  ! Properly end the run. 
  !---------------------------
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  character(LEN=80)::str

  call output_timer(.false., str)

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD, info)
  ! call MPI_FINALIZE(info)
#endif

  stop
end subroutine clean_end

subroutine clean_stop
  !-----------------------------------------------------
  ! This subroutine brings the program to a halt after
  ! an error.
  !-----------------------------------------------------
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif

#ifndef WITHOUTMPI
  call MPI_ABORT(MPI_COMM_WORLD, 2, info)
#endif

  stop
end subroutine clean_stop

subroutine deallocate_amr
  implicit none
end subroutine deallocate_amr

subroutine deallocate_pm
  implicit none
end subroutine deallocate_pm

subroutine deallocate_poisson
  implicit none
end subroutine deallocate_poisson 
