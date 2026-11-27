! sbs_imf: IMF tables for the star-by-star mode of the FG star formation (pop2_model=5, star_by_star)
! and the Pop III IMF of pop3_model=7/8. Ported from patch/ks (merge_ks); include 'mpif.h' removed.
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
module sbs_imf
  use amr_parameters,only:dp,sbs_Mmin,sbs_imf_name
  use amr_parameters,only:pop3_model
  implicit none
  private
  public ChabrierIMF, ChonIMF, initialize_IMF, get_Mfrac_unresolved_at_Z, get_Pcum_at_Z, get_massbin_unresolved,&
    initialize_Pop3_IMF
  !integer,parameter :: dp=kind(1.0D0)

  real(dp),allocatable,dimension(:,:),public :: sbs_Pcum ! cumulative probability for resolved stars
  real(dp),allocatable,dimension(:),public :: sbs_Mfrac_unresolved ! mass fraction of unresolved stars

  !character(len=16),parameter,public :: sbs_imf_name='Chabrier' ! IMF name (defined in amr_parameters.f90)
  !character(len=16),parameter,public :: sbs_imf_name='Chon' ! IMF name

  ! Parameters for resolved IMF
  !real(dp),parameter,public :: sbs_Mmin = 1d0   ! minimum resolved stellar mass [Msun] (defined in amr_parameters.f90)
  real(dp),parameter,public :: sbs_Mmax = 1d2   ! maximum stellar mass [Msun]
  integer,parameter,public :: sbs_NM   = 100    ! number of mass bins

  ! Parameters for metallicity grid
  real(dp),parameter,public :: sbs_Zmin = 1d-4  ! minimum metallicity [Zsun]
  real(dp),parameter,public :: sbs_Zmax = 1d0   ! maximum metallicity [Zsun]
  integer,parameter,public :: sbs_NZ   = 25     ! number of metallicity bins


  ! Parameters for IMF integration
  real(dp),parameter,public :: Mmin_imf = 0.1d0  ! lower bound of IMF [Msun]
  real(dp),parameter,public :: Mmax_imf = 100d0  ! upper bound of IMF [Msun]
  integer,parameter,public :: Ninteg    = 10000  ! number of bins for integration

  ! artificial IMF for Pop III stars (used only when pop3_model = 7 or 8)
  real(dp),allocatable,dimension(:),public :: p3_Pcum ! cumulative probability for Pop III stars
  real(dp),parameter,public :: p3_Mmin_imf = 1d1  ! lower bound of Pop III IMF [Msun]
  real(dp),parameter,public :: p3_Mmax_imf = 5d2  ! upper bound of Pop III IMF [Msun]
  integer,parameter,public :: p3_NM = 100  ! number of mass bins
  integer,parameter,public :: p3_Ninteg    = 10000  ! number of bins for integration

contains

  !====================================================================
  ! Chabrier IMF (see Hirai et al. 2021)
  !====================================================================
  subroutine ChabrierIMF(M,dndM)
    real(dp), intent(in) :: M
    real(dp), intent(out):: dndM
    ! if (0.1d0 < M .and. M <= 1d0) then
    !    dndM = exp(-log10(M/0.079d0)**2/(2*0.69d0**2)) / &
    !           exp(-log10(1d0/0.079d0)**2/(2*0.69d0**2))
    ! else if (1d0 < M .and. M <= 100d0) then
    !    dndM = M**(-2.3d0)
    ! else
    !    dndM = 0d0
    ! end if

    !upper and lower limits are determined by Mmin_imf and Mmax_imf, respectively
    if (M <= 1d0) then
      dndM = exp(-log10(M/0.079d0)**2/(2*0.69d0**2)) / &
        exp(-log10(1d0/0.079d0)**2/(2*0.69d0**2))
    else
      dndM = M**(-2.3d0)
    end if
  end subroutine ChabrierIMF

  !====================================================================
  ! Salpeter IMF (1Msun < M < 100 Msun)
  !====================================================================
  subroutine SalpeterIMF1to100Msun(M,dndM)
    real(dp), intent(in) :: M
    real(dp), intent(out):: dndM
    if (M >= 1d0 .and. M <= 1d2) then
      dndM = M**(-2.35d0)
    else
      dndM = 0d0
    end if
  end subroutine SalpeterIMF1to100Msun


  !====================================================================
  ! Salpeter IMF (10Msun < M < 500 Msun)
  !====================================================================
  subroutine SalpeterIMF10to500Msun(M,dndM)
    real(dp), intent(in) :: M
    real(dp), intent(out):: dndM
    if (M >= 1d1 .and. M <= 5d2) then
      dndM = M**(-2.35d0)
    else
      dndM = 0d0
    end if
  end subroutine SalpeterIMF10to500Msun

  !====================================================================
  ! Flat IMF (10Msun < M < 500 Msun)
  !====================================================================
  subroutine FlatIMF10to500Msun(M,dndM)
    real(dp), intent(in) :: M
    real(dp), intent(out):: dndM
    if (M >= 1d1 .and. M <= 5d2) then
      dndM = M**(0d0) !flat IMF
    else
      dndM = 0d0
    end if
  end subroutine FlatIMF10to500Msun


  !====================================================================
  ! Metallicity-dependent IMF (Chon et al. 2024)
  !====================================================================
  subroutine ChonIMF(M,Z,dndM)
    real(dp), intent(in) :: M,Z !M[Msun] and Z[Zsun]
    real(dp), intent(out):: dndM
    real(dp) :: m_min, m_max, c, logZ, alpha, logm0, m0
    m_min = 0.04d0
    m_max = 150d0
    c = 1.6d0
    logZ = log10(Z)
    alpha = 2.3d0 + 0.33d0*logZ
    logm0 = 0.2d0 + 0.45d0*logZ
    m0 = 10d0**logm0
    dndM = M**(-alpha) * (1d0 - exp(-(M/m0)**c)) * exp(-(m_min/M) - (M/m_max))
  end subroutine ChonIMF

  !====================================================================
  ! Initialize IMF tables
  !====================================================================
  subroutine initialize_IMF

    use amr_commons,only:myid
    integer :: iM, iZ
    real(dp), allocatable :: M(:), dndM(:), cum_num(:), cum_mass(:)
    real(dp) :: Z, dM, total_num, total_mass, Mval
    integer, save :: ifirst = 0

    if (ifirst /= 0) return
    ifirst = 1


    !-----------------------------------------------------------
    ! Print parameter settings with formatted output
    !-----------------------------------------------------------
    if (myid == 1) then
      print *, "============================================="
      print *, "Initialize IMF table for star-by-star treatment"
      write(*,'(A,1X,A)')    " IMF name             = ", trim(sbs_imf_name)
      write(*,'(A,1PE12.4,A,1PE12.4,A)') " Mass range (resolved)= ", sbs_Mmin, " - ", sbs_Mmax, " Msun"
      write(*,'(A,I6)')       " Number of mass bins  = ", sbs_NM
      write(*,'(A,1PE12.4,A,1PE12.4,A)') " Metallicity range    = ", sbs_Zmin, " - ", sbs_Zmax, " Zsun"
      write(*,'(A,I6)')       " Number of Z bins     = ", sbs_NZ
      write(*,'(A,1PE12.4,A,1PE12.4,A)') " Integration range    = ", Mmin_imf, " - ", Mmax_imf, " Msun"
      write(*,'(A,I8)')       " Integration bins     = ", Ninteg
      print *, "============================================="
    end if

    allocate(sbs_Pcum(sbs_NM, sbs_NZ))
    allocate(sbs_Mfrac_unresolved(sbs_NZ))

    ! Loop over metallicity bins
    do iZ = 1, sbs_NZ
      Z = 10d0**( log10(sbs_Zmin) + (iZ-1)*(log10(sbs_Zmax)-log10(sbs_Zmin))/(sbs_NZ-1) )

      allocate(M(Ninteg), dndM(Ninteg), cum_num(Ninteg), cum_mass(Ninteg))

      ! Log-spaced mass bins
      do iM = 1, Ninteg
        M(iM) = 10d0**( log10(Mmin_imf) + (iM-1)*(log10(Mmax_imf)-log10(Mmin_imf))/(Ninteg-1) )
        if (trim(sbs_imf_name) == 'Chabrier') then
          call ChabrierIMF(M(iM), dndM(iM))
        else
          call ChonIMF(M(iM), Z, dndM(iM))
        end if
      end do

      !-----------------------------------------------------------
      ! (1) Cumulative number distribution
      cum_num(1) = 0d0
      do iM=2,Ninteg
        dM = M(iM) - M(iM-1)
        cum_num(iM) = cum_num(iM-1) + 0.5d0*(dndM(iM) + dndM(iM-1))*dM
      end do
      total_num = cum_num(Ninteg)

      ! Resample to sbs_NM points between sbs_Mmin and sbs_Mmax
      do iM=1,sbs_NM
        Mval = 10d0**( log10(sbs_Mmin) + (iM-1)*(log10(sbs_Mmax)-log10(sbs_Mmin))/(sbs_NM-1) )
        call interp_log(M, cum_num/total_num, Ninteg, Mval, sbs_Pcum(iM,iZ))
      end do

      ! cumulative probability for resolved stars (from 0 to 1 between sbs_Mmin and sbs_Mmax)
      sbs_Pcum(:,iZ) = (sbs_Pcum(:,iZ)-sbs_Pcum(1,iZ))/(sbs_Pcum(sbs_NM,iZ)-sbs_Pcum(1,iZ))

      !-----------------------------------------------------------
      ! (2) Cumulative mass distribution (for unresolved fraction)
      cum_mass(1) = 0d0
      do iM=2,Ninteg
        dM = M(iM) - M(iM-1)
        cum_mass(iM) = cum_mass(iM-1) + 0.5d0*((M(iM)*dndM(iM))+(M(iM-1)*dndM(iM-1)))*dM
      end do
      total_mass = cum_mass(Ninteg)

      ! Fraction of stellar mass below the resolved minimum mass
      call interp_log(M, cum_mass/total_mass, Ninteg, sbs_Mmin, sbs_Mfrac_unresolved(iZ))

      deallocate(M,dndM,cum_num,cum_mass)
    end do

  end subroutine initialize_IMF

  !====================================================================
  ! Log-linear interpolation: x in log space, y in linear space
  !====================================================================
  subroutine interp_log(xarr,yarr,N,xval,yval)
    real(dp), intent(in) :: xarr(:), yarr(:)
    integer, intent(in)  :: N
    real(dp), intent(in) :: xval
    real(dp), intent(out):: yval
    integer :: i

    if (xval <= xarr(1)) then
      yval = yarr(1)
    else if (xval >= xarr(N)) then
      yval = yarr(N)
    else
      do i=1,N-1
        if (xarr(i) <= xval .and. xval < xarr(i+1)) then
          yval = yarr(i) + (yarr(i+1)-yarr(i))*(log10(xval)-log10(xarr(i)))/ &
            (log10(xarr(i+1))-log10(xarr(i)))
          return
        end if
      end do
    end if
  end subroutine interp_log

  !====================================================================
  ! Interpolate cumulative probability array for given Z
  !====================================================================
  subroutine get_Pcum_at_Z(Zval, Pcum)
    real(dp), intent(in)  :: Zval
    real(dp), intent(out) :: Pcum(sbs_NM)
    integer :: i, iz
    real(dp) :: logZval, logZmin, logZmax, dlogZ, w

    ! Guard input range
    if (Zval <= sbs_Zmin) then
      Pcum(:) = sbs_Pcum(:,1)
      return
    else if (Zval >= sbs_Zmax) then
      Pcum(:) = sbs_Pcum(:,sbs_NZ)
      return
    end if

    ! Grid in logZ
    logZmin = log10(sbs_Zmin)
    logZmax = log10(sbs_Zmax)
    dlogZ   = (logZmax - logZmin) / real(sbs_NZ-1,dp)
    logZval = log10(Zval)

    ! Find interval index
    iz = int( (logZval - logZmin)/dlogZ ) + 1

    ! Linear weight
    w = (logZval - (logZmin + (iz-1)*dlogZ)) / dlogZ

    ! Interpolation
    Pcum(:) = (1.0d0-w)*sbs_Pcum(:,iz) + w*sbs_Pcum(:,iz+1)
  end subroutine get_Pcum_at_Z

  !====================================================================
  ! Interpolate unresolved mass fraction for given Z
  !====================================================================
  subroutine get_Mfrac_unresolved_at_Z(Zval, Mfrac_unres)
    real(dp), intent(in)  :: Zval
    real(dp), intent(out) :: Mfrac_unres
    integer :: iz
    real(dp) :: logZval, logZmin, logZmax, dlogZ, w

    ! Guard input range
    if (Zval <= sbs_Zmin) then
      Mfrac_unres = sbs_Mfrac_unresolved(1)
      return
    else if (Zval >= sbs_Zmax) then
      Mfrac_unres = sbs_Mfrac_unresolved(sbs_NZ)
      return
    end if

    ! Grid in logZ
    logZmin = log10(sbs_Zmin)
    logZmax = log10(sbs_Zmax)
    dlogZ   = (logZmax - logZmin) / real(sbs_NZ-1,dp)
    logZval = log10(Zval)

    ! Find interval index
    iz = int( (logZval - logZmin)/dlogZ ) + 1

    ! Linear weight
    w = (logZval - (logZmin + (iz-1)*dlogZ)) / dlogZ

    ! Interpolation
    Mfrac_unres = (1.0d0-w)*sbs_Mfrac_unresolved(iz) + w*sbs_Mfrac_unresolved(iz+1)
  end subroutine get_Mfrac_unresolved_at_Z

  !====================================================================
  ! Get expected number of unresolved stars per logarithmic mass bin
  ! normalized so that the total unresolved stellar mass = sbs_Mmin
  !====================================================================
  subroutine get_massbin_unresolved(Zval, Nbin_unres, Mbin, dN)
    implicit none
    integer, intent(in) :: Nbin_unres
    real(dp), intent(in) :: Zval
    real(dp), intent(out) :: Mbin(Nbin_unres), dN(Nbin_unres)
    !
    integer :: i, j, Ntmp
    real(dp) :: logMmin, logMmax, dlogM, logM1, logM2
    real(dp) :: M1, M2, dndM1, dndM2, dM
    real(dp) :: integralN, Mtot_unres, scale
    real(dp), allocatable :: M(:), dndM(:)
    !
    !-----------------------------------------------------------
    ! Step 0: setup integration arrays
    !-----------------------------------------------------------
    Ntmp = Ninteg
    allocate(M(Ntmp), dndM(Ntmp))

    do i=1,Ntmp
      M(i) = 10d0**( log10(Mmin_imf) + (i-1)*(log10(Mmax_imf)-log10(Mmin_imf))/(Ntmp-1) )
      if (trim(sbs_imf_name) == 'Chabrier') then
        call ChabrierIMF(M(i), dndM(i))
      else
        call ChonIMF(M(i), Zval, dndM(i))
      end if
    end do

    !-----------------------------------------------------------
    ! Step 1: compute the total unresolved stellar mass
    !         in the range [Mmin_imf, sbs_Mmin]
    !-----------------------------------------------------------
    Mtot_unres = 0d0
    do i=2,Ntmp
      if (M(i) <= sbs_Mmin) then
        dM = M(i) - M(i-1)
        Mtot_unres = Mtot_unres + 0.5d0*((M(i)*dndM(i))+(M(i-1)*dndM(i-1)))*dM
        !print * , "i, M(i), dM, dndM(i) = " , i, M(i), dM, dndM(i)
      else if (M(i-1) < sbs_Mmin .and. M(i) > sbs_Mmin) then
        ! last partial interval up to sbs_Mmin
        dM = sbs_Mmin - M(i-1)
        dndM1 = dndM(i-1)
        dndM2 = dndM(i-1) + (dndM(i)-dndM(i-1))*(sbs_Mmin-M(i-1))/(M(i)-M(i-1))
        Mtot_unres = Mtot_unres + 0.5d0*((sbs_Mmin*dndM2)+(M(i-1)*dndM1))*dM
        !print * , "last i, M(i), dM, dndM1, dndM2, Mtot_unres = " , i, M(i), dM, dndM1,dndM2, Mtot_unres
      end if

    end do

    !-----------------------------------------------------------
    ! Step 2: determine scaling factor so that
    !         total unresolved stellar mass = sbs_Mmin
    !-----------------------------------------------------------
    scale = sbs_Mmin / Mtot_unres

    !-----------------------------------------------------------
    ! Step 3: define logarithmic bins and compute expected counts
    !-----------------------------------------------------------
    logMmin = log10(Mmin_imf)
    logMmax = log10(sbs_Mmin)
    dlogM   = (logMmax - logMmin)/real(Nbin_unres,dp)

    do j=1,Nbin_unres
      logM1 = logMmin + (j-1)*dlogM
      logM2 = logMmin + j*dlogM
      M1 = 10d0**logM1
      M2 = 10d0**logM2

      ! arithmetic mean of the bin edges
      Mbin(j) = (M1+M2)/2d0

      ! integrate number distribution in this bin
      integralN = 0d0
      do i=2,Ntmp
        if (M(i) >= M1 .and. M(i-1) <= M2) then
          dM = min(M(i),M2) - max(M(i-1),M1)
          if (dM > 0d0) then
            dndM1 = dndM(i-1)
            dndM2 = dndM(i)
            integralN = integralN + 0.5d0*(dndM1+dndM2)*dM
            !print *, "i, M1, M2, dM, dndM1, dndM2, integralN = ", i, M1, M2, dM, dndM1, dndM2, integralN
          end if
        end if
      end do

      dN(j) = scale * integralN

      !print *, "j, dN(j) = ", j, dN(j)
    end do

    !print *, "scale = ", scale

    deallocate(M,dndM)
  end subroutine get_massbin_unresolved


  !====================================================================
  ! Initialize Pop III IMF tables
  !====================================================================
  subroutine initialize_Pop3_IMF

    use amr_commons,only:myid
    integer :: iM
    real(dp), allocatable :: M(:), dndM(:), cum_num(:), cum_mass(:)
    real(dp) ::  dM, total_num, total_mass, Mval
    integer, save :: ifirst = 0

    if (ifirst /= 0) return
    ifirst = 1


    !-----------------------------------------------------------
    ! Print parameter settings with formatted output
    !-----------------------------------------------------------
    if (myid == 1) then
      print *, "============================================="
      print '(A, I7)', "Initialize Pop III IMF table for pop3_model = ", pop3_model
      if (pop3_model == 7) then
        write(*,'(A)')    "Salpeter IMF (10-500 Msun)"
      else if (pop3_model == 8) then
        write(*,'(A)')    "Flat IMF (10-500 Msun)"
      end if
      write(*,'(A,1PE12.4,A,1PE12.4,A)') " Mass range = ", p3_Mmin_imf, " - ", p3_Mmax_imf, " Msun"
      write(*,'(A,I6)')       " Number of mass bins  = ", p3_NM
      write(*,'(A,I8)')       " Integration bins     = ", p3_Ninteg
      print *, "============================================="
    end if

    allocate(p3_Pcum(p3_NM))

    allocate(M(Ninteg), dndM(Ninteg), cum_num(Ninteg), cum_mass(Ninteg))

    ! Log-spaced mass bins
    do iM = 1, Ninteg
      M(iM) = 10d0**( log10(p3_Mmin_imf) + (iM-1)*(log10(p3_Mmax_imf)-log10(p3_Mmin_imf))/(Ninteg-1) )
      if (pop3_model == 7) then
        call SalpeterIMF10to500Msun(M(iM), dndM(iM))
        !print *, "iM, M(iM), dndM(iM): ", iM, M(iM), dndM(iM)
      else if (pop3_model == 8) then
        call FlatIMF10to500Msun(M(iM), dndM(iM))
      end if
    end do

    !-----------------------------------------------------------
    ! (1) Cumulative number distribution
    cum_num(1) = 0d0
    do iM=2,Ninteg
      dM = M(iM) - M(iM-1)
      cum_num(iM) = cum_num(iM-1) + 0.5d0*(dndM(iM) + dndM(iM-1))*dM
    end do
    total_num = cum_num(Ninteg)

    ! Resample to p3_NM points
    do iM=1,p3_NM
      Mval = 10d0**( log10(p3_Mmin_imf) + (iM-1)*(log10(p3_Mmax_imf)-log10(p3_Mmin_imf))/(p3_NM-1) )
      call interp_log(M, cum_num/total_num, Ninteg, Mval, p3_Pcum(iM))
      if (myid == 1) &
        print *, "KS DEBUG iM, Mval, p3_Pcum(iM): ",  iM, Mval, p3_Pcum(iM)
    end do

    !print *, "KS DEBUG p3_Pcum(:): ",  p3_Pcum(:)
    ! print *, "KS DEBUG STOP HERE"
    ! stop

    deallocate(M,dndM,cum_num)


  end subroutine initialize_Pop3_IMF



end module sbs_imf

