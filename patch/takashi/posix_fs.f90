!=======================================================================
! posix_fs: minimal directory operations through libc (ISO_C_BINDING).
!
! Used by dump_restart to write restart files into a temporary directory
! and publish them with an atomic rename once every rank has finished.
! No shell is spawned (works on Cray compute nodes where system() does
! not) and no compiler specific PXF* routines are needed, so the same
! code compiles with Intel (ifort/ifx) and gfortran.
!
! fs_mkdir(path)        mkdir path, mode 0755      (0 on success)
! fs_rename(old,new)    rename(2)                  (0 on success)
! fs_rmtree(path)       rm -rf path via nftw+remove (0 on success,
!                       -1 if path does not exist)
!=======================================================================
module posix_fs
  use iso_c_binding
  implicit none
  private
  public :: fs_mkdir, fs_rename, fs_rmtree

  ! flags of nftw(3) as defined in glibc <ftw.h>
  integer(c_int), parameter :: FTW_PHYS  = 1   ! do not follow symbolic links
  integer(c_int), parameter :: FTW_DEPTH = 8   ! post-order: children before their directory

  interface
    function c_remove(path) bind(C, name='remove') result(r)
      import :: c_char, c_int
      character(kind=c_char), dimension(*) :: path
      integer(c_int) :: r
    end function c_remove
    function c_rename(old, new) bind(C, name='rename') result(r)
      import :: c_char, c_int
      character(kind=c_char), dimension(*) :: old, new
      integer(c_int) :: r
    end function c_rename
    function c_mkdir(path, mode) bind(C, name='mkdir') result(r)
      import :: c_char, c_int
      character(kind=c_char), dimension(*) :: path
      integer(c_int), value :: mode
      integer(c_int) :: r
    end function c_mkdir
    function c_nftw(path, fn, nopenfd, flags) bind(C, name='nftw') result(r)
      import :: c_char, c_int, c_funptr
      character(kind=c_char), dimension(*) :: path
      type(c_funptr), value :: fn
      integer(c_int), value :: nopenfd, flags
      integer(c_int) :: r
    end function c_nftw
  end interface

contains

  function fs_mkdir(path) result(ierr)
    character(len=*), intent(in) :: path
    integer :: ierr
    ierr = int(c_mkdir(trim(path)//c_null_char, int(o'755', c_int)))
  end function fs_mkdir

  function fs_rename(old, new) result(ierr)
    character(len=*), intent(in) :: old, new
    integer :: ierr
    ierr = int(c_rename(trim(old)//c_null_char, trim(new)//c_null_char))
  end function fs_rename

  function fs_rmtree(path) result(ierr)
    character(len=*), intent(in) :: path
    integer :: ierr
    ierr = int(c_nftw(trim(path)//c_null_char, c_funloc(rm_entry), 64_c_int, ior(FTW_PHYS, FTW_DEPTH)))
  end function fs_rmtree

  ! callback of nftw: with FTW_DEPTH a directory is visited after its
  ! content, so remove() (unlink for files, rmdir for directories) suffices
  function rm_entry(fpath, sb, typeflag, ftwbuf) bind(C) result(r)
    character(kind=c_char), dimension(*) :: fpath
    type(c_ptr), value :: sb, ftwbuf
    integer(c_int), value :: typeflag
    integer(c_int) :: r
    r = c_remove(fpath)
  end function rm_entry

end module posix_fs
