module besselfilter_design
  !! Final Fortran 90 module for automatic design of Bessel IIR filters
  !! - LPF / HPF follow the coefficient formulas in the paper’s Appendix
  !! - BPF follows the standard “analog-LP poles → analog-BP poles → bilinear” flow
  !!   (mathematically equivalent to the paper’s A-variable derivation)
  !! - Includes: Bessel polynomial generator, robust root solver,
  !!   frequency response / phase unwrapping / group delay utilities
  !! reference: Katsumata, A. (1993) Automatic designing of bessel digital filters, Q. J. Seismo., 56, 17-34
  use nrtype, only : dp
  use constants
  implicit none


contains
  !======================================================================
  ! Utility: log-factorial with real(dp) accumulation (stable for n<=50)
  !======================================================================
  pure real(dp) function log_factorial(m) result(v)
    integer, intent(in) :: m
    integer :: i
    v = 0.0_dp
    if (m <= 1) return
    do i=2,m
      v = v + log(real(i,dp))
    end do
  end function log_factorial

  !======================================================================
  ! GETBPL
  !   Build Bessel polynomial T_n(s) (ascending coeffs), normalized: a(0)=1
  !   COEF(k) = (2n-k)! / ( 2^(n-k) (n-k)! k! ) , then divided by COEF(0)
  !======================================================================
  subroutine GETBPL(N, COEF, MAXDBS, IERR)
    integer, intent(in)  :: N
    real(dp), intent(out) :: COEF(0:50)
    integer, intent(out) :: MAXDBS, IERR
    integer :: k
    real(dp) :: c0

    IERR = 0
    if (N <= 0 .or. N > 20) then
      IERR = 1
      return
    end if

    COEF(:) = 0.0_dp
    MAXDBS  = N

    do k = 0, N
      COEF(k) = exp( log_factorial(2*N-k) - (N-k)*log(2.0_dp) &
                   - log_factorial(N-k) - log_factorial(k) )
    end do
    c0 = COEF(0)
    if (c0 <= 0.0_dp) then
      IERR = 2
      return
    end if
    COEF(0:N) = COEF(0:N) / c0   ! normalize constant term to 1.0
  end subroutine GETBPL

  !======================================================================
  ! VALBPL
  !   Evaluate |T_n(i*x)| for a given x>=0 using Horner’s rule
  !======================================================================
  subroutine VALBPL(X, COEF, MAXDBS, VAL)
    real(dp), intent(in) :: X
    real(dp), intent(in) :: COEF(0:50)
    integer, intent(in)  :: MAXDBS
    real(dp), intent(out) :: VAL
    complex(dp) :: s, p
    integer :: k

    s = cmplx(0.0_dp, X, kind=dp)    ! s = i*x
    p = cmplx(0.0_dp, 0.0_dp, kind=dp)
    do k = MAXDBS, 0, -1
      p = p*s + cmplx(COEF(k), 0.0_dp, kind=dp)
    end do
    VAL = abs(p)
  end subroutine VALBPL

  !======================================================================
  ! DETBTX
  !   Find TX (=x_p) such that |T_n(i*TX)| = sqrt(1 + AP^2)
  !   Robust bracketing + bisection (monotone in x>=0)
  !======================================================================
  subroutine DETBTX(AP, COEF, MAXDBS, TX, IERR)
    real(dp), intent(in) :: AP
    real(dp), intent(in) :: COEF(0:50)
    integer, intent(in)  :: MAXDBS
    real(dp), intent(out) :: TX
    integer, intent(out) :: IERR

    real(dp) :: BP, lo, hi, vlo, vhi, mid, vmid
    integer :: it

    IERR = 0
    if (AP < 0.0_dp) then
      IERR = 1
      return
    end if
    BP = sqrt(1.0_dp + AP*AP)

    lo = 0.0_dp
    call VALBPL(lo, COEF, MAXDBS, vlo)      ! vlo ~ 1.0
    if (vlo > BP) then
      IERR = 2
      return
    end if

    hi = 1.0_dp
    do it = 1, 200
      call VALBPL(hi, COEF, MAXDBS, vhi)
      if (vhi >= BP) exit
      hi = hi * 2.0_dp
      if (hi > 1.0e6_dp) then
        IERR = 3
        return
      end if
    end do

    do it = 1, 200
      mid = 0.5_dp*(lo+hi)
      call VALBPL(mid, COEF, MAXDBS, vmid)
      if (abs(vmid - BP) < 1.0e-12_dp*BP) exit
      if (vmid >= BP) then
        hi = mid
      else
        lo = mid
      end if
    end do
    TX = 0.5_dp*(lo+hi)
  end subroutine DETBTX

  !======================================================================
  ! cauchy_bound
  !   Cauchy root bound for ascending coefficients:
  !   |z| <= 1 + max_{k<N} |a_k/a_N|
  !======================================================================
  pure real(dp) function cauchy_bound(COEF, N) result(rad)
    real(dp), intent(in) :: COEF(0:50)
    integer, intent(in)  :: N
    integer :: k
    real(dp) :: an, mx
    an = abs(COEF(N))
    if (an <= 0.0_dp) then
      rad = 1.0_dp
      return
    end if
    mx = 0.0_dp
    do k=0,N-1
      mx = max(mx, abs(COEF(k))/an)
    end do
    rad = 1.0_dp + mx
  end function cauchy_bound

  !======================================================================
  ! poly_eval_asc
  !   Evaluate ascending polynomial: P(z) = sum_{k=0..N} c(k) z^k
  !======================================================================
  subroutine poly_eval_asc(c, N, z, pz)
    real(dp), intent(in)    :: c(0:50)
    integer , intent(in)    :: N
    complex(dp), intent(in) :: z
    complex(dp), intent(out):: pz
    integer :: k
    pz = cmplx(0.0_dp, 0.0_dp, kind=dp)
    do k = N, 0, -1
      pz = pz*z + cmplx(c(k), 0.0_dp, kind=dp)
    end do
  end subroutine poly_eval_asc

  !======================================================================
  ! BAIRSTOW (replacement)
  !   Robust root solver using scaled Durand–Kerner (Weierstrass) iteration
  !   Interface kept identical:
  !     COEF ascending (0..N), ROOT(1..N), IERR=0 ok / 2 not converged
  !======================================================================
  subroutine BAIRSTOW(COEF, N, ROOT, IERR)
    real(dp), intent(in)     :: COEF(0:50)
    integer , intent(in)     :: N
    complex(dp), intent(out) :: ROOT(50)
    real(dp) :: b(0:50)
    integer , intent(out)    :: IERR

    complex(dp), allocatable :: z(:), znew(:)
    real(dp) :: rad, alpha, tol, maxchg
    complex(dp) :: pz, denom, dz
    integer :: k, j, it, maxit

    IERR = 0
    ROOT(:) = cmplx(0.0_dp, 0.0_dp, kind=dp)
    if (N <= 0 .or. N > 20) then
      IERR = 1
      return
    end if

    ! Scale variable: s = alpha * w  (use Cauchy bound for alpha ≈ rad/2)
    rad   = cauchy_bound(COEF, N)
    if (rad <= 0.0_dp) rad = 1.0_dp
    alpha = max(1.0_dp, rad/2.0_dp)

    ! Build monic Q(w) with ascending coeffs: w^N + sum_k b_k w^k
    !   b_k = (a_k/a_N)*alpha^(k-N)
    do k=0,N-1
      b(k) = (COEF(k)/COEF(N)) * alpha**(real(k-N,dp))
    end do
    b(N) = 1.0_dp

    ! Durand–Kerner on Q(w) (unit circle initial guesses)
    allocate(z(N), znew(N))
    do k=1,N
      z(k) = cmplx( cos(2.0_dp*PI*real(k-1,dp)/real(N,dp)), &
                    sin(2.0_dp*PI*real(k-1,dp)/real(N,dp)), kind=dp )
    end do
    tol   = 1.0e-14_dp
    maxit = 1000

    do it=1,maxit
      maxchg = 0.0_dp
      do k=1,N
        call poly_eval_asc(b, N, z(k), pz)
        denom = cmplx(1.0_dp, 0.0_dp, kind=dp)
        do j=1,N
          if (j /= k) denom = denom * (z(k) - z(j))
        end do
        if (abs(denom) < 1.0e-30_dp) denom = denom + cmplx(1.0e-16_dp, 1.0e-16_dp, kind=dp)
        dz      = pz / denom
        znew(k) = z(k) - dz
        maxchg  = max(maxchg, abs(dz))
      end do
      z = znew
      if (maxchg < tol) exit
    end do

    if (maxchg >= tol) then
      IERR = 2
    end if

    do k=1,N
      ROOT(k) = alpha * z(k)     ! map back to s-plane
    end do

    call sort_roots_by_imag_desc(ROOT(1:N))  ! stable ordering
    deallocate(z, znew)
  end subroutine BAIRSTOW

  !======================================================================
  ! sort_roots_by_imag_desc (assumed-shape, safe for any length)
  !======================================================================
  subroutine sort_roots_by_imag_desc(z)
    complex(dp), intent(inout) :: z(:)
    integer :: i, j, n
    complex(dp) :: tmp
    n = size(z)
    do i=1, n-1
      do j=i+1, n
        if (aimag(z(j)) > aimag(z(i))) then
          tmp = z(i); z(i) = z(j); z(j) = tmp
        end if
      end do
    end do
  end subroutine sort_roots_by_imag_desc

  !======================================================================
  ! FREQRESP_POINT
  !   Single-point frequency response (linear AMP, principal PHASE)
  !   Transfer per section: (1 + a1 z^-1 + a2 z^-2)/(1 + b1 z^-1 + b2 z^-2)
  !======================================================================
  subroutine FREQRESP_POINT(H, M, GN, DT, FREQ, AMP, PHASE)
    real(dp), intent(in) :: H(*), GN, DT, FREQ
    integer , intent(in) :: M
    real(dp), intent(out):: AMP, PHASE

    real(dp) :: omega
    complex(dp) :: z1, z2, num, den, Hc
    integer :: i
    real(dp) :: a1, a2, b1, b2

    omega = 2.0_dp*PI*FREQ
    z1 = exp(cmplx(0.0_dp, -omega*DT, kind=dp))   ! z^-1
    z2 = z1*z1

    Hc = cmplx(GN,0.0_dp,kind=dp)
    do i=1,M
      a1 = H(4*i-3); a2 = H(4*i-2)
      b1 = H(4*i-1); b2 = H(4*i)
      num = (1.0_dp + a1*z1 + a2*z2)
      den = (1.0_dp + b1*z1 + b2*z2)
      Hc  = Hc * (num/den)
    end do
    AMP   = abs(Hc)
    PHASE = atan2(aimag(Hc), real(Hc))
  end subroutine FREQRESP_POINT

  !======================================================================
  ! FREQRESP_SWEEP (with phase unwrapping)
  !======================================================================
  subroutine FREQRESP_SWEEP(H, M, GN, DT, FREQ, NF, AMP, PHASEU)
    real(dp), intent(in) :: H(*), GN, DT
    integer , intent(in) :: M, NF
    real(dp), intent(in) :: FREQ(:)
    real(dp), intent(out):: AMP(:), PHASEU(:)

    integer :: k
    real(dp) :: ph, phw_prev, dph

    do k=1,NF
      call FREQRESP_POINT(H, M, GN, DT, FREQ(k), AMP(k), ph)
      if (k == 1) then
        PHASEU(k) = ph
      else
        ! unwrap by nearest to previous wrapped value
        phw_prev = atan2(sin(PHASEU(k-1)), cos(PHASEU(k-1)))
        dph = ph - phw_prev
        do while (dph >  PI); dph = dph - 2.0_dp*PI; end do
        do while (dph < -PI); dph = dph + 2.0_dp*PI; end do
        PHASEU(k) = PHASEU(k-1) + dph
      end if
    end do
  end subroutine FREQRESP_SWEEP

  !======================================================================
  ! GROUPDELAY_SWEEP
  !   tau_g(f) = - dphi/domega   (central difference)
  !======================================================================
  subroutine GROUPDELAY_SWEEP(FREQ, NF, PHASEU, TAU)
    real(dp), intent(in)  :: FREQ(:), PHASEU(:)
    integer , intent(in)  :: NF
    real(dp), intent(out) :: TAU(:)

    integer :: k
    real(dp) :: domega

    TAU(:) = 0.0_dp
    if (NF < 3) return
    do k=2, NF-1
      domega = 2.0_dp*PI * (FREQ(k+1) - FREQ(k-1))
      if (abs(domega) > 0.0_dp) then
        TAU(k) = -(PHASEU(k+1) - PHASEU(k-1)) / domega
      else
        TAU(k) = 0.0_dp
      end if
    end do
    TAU(1)  = TAU(2)
    TAU(NF) = TAU(NF-1)
  end subroutine GROUPDELAY_SWEEP

  !======================================================================
  ! BESLPN (Bessel LPF design)
  !   - Follows Appendix coefficient form for each quadratic factor
  !   - Final 1-point normalization at dp to meet |H(dp)|=1/sqrt(1+AP^2)
  !======================================================================
  subroutine BESLPN(H, M, GN, N, FPK, AP, DT, IERR)
    real(dp), intent(out), allocatable :: H(:) !!output: :H(4 * M)
    real(dp), intent(out) :: gn
    integer , intent(out) :: M, IERR  !!M: filter order
    integer , intent(in)  :: N        !!N: order of polynominal
    real(dp), intent(in)  :: FPK, AP, DT   !!FPK: not normalized

    real(dp) :: COEF(0:50)
    complex(dp) :: ZROOT(50)
    integer :: MAXDBS, i, sec, idx_real
    real(dp) :: wp2, tp, TX, C, sr, sa2, alpha, b1, b2, gsec
    real(dp) :: amp_dp, ph_dp, target
    logical :: has_real

    IERR=0; GN=1.0_dp; M=0
    if (N<=0 .or. N>20) then; IERR=1; return; end if
    if (DT<=0.0_dp) then; IERR=2; return; end if
    if (FPK<=0.0_dp .or. FPK>=0.5_dp/DT) then; IERR=3; return; end if
    if (AP<0.0_dp) then; IERR=4; return; end if

    wp2 = PI*FPK*DT                ! = ωp*T/2
    if (wp2<=0.0_dp .or. wp2>=PI/2.0_dp) then; IERR=5; return; end if
    tp  = tan(wp2)                 ! tan(ωp*T/2)

    call GETBPL(N, COEF, MAXDBS, IERR); if (IERR/=0) return
    call DETBTX(AP, COEF, MAXDBS, TX, IERR); if (IERR/=0) return

    C = TX / tp                    ! LPF: C = xp / tan(ωp*T/2)

    call BAIRSTOW(COEF, MAXDBS, ZROOT, IERR); if (IERR/=0) return

    M = N/2
    if (mod(N,2) == 1) M = M + 1
    allocate(h(4*M))

    has_real = .false.; idx_real = -1
    do i=1,N
      if (abs(aimag(ZROOT(i))) < 1.0e-10_dp) then
        has_real = .true.; idx_real=i; exit
      end if
    end do

    sec = 0
    do i=1,N
      if (i == idx_real) cycle
      if (aimag(ZROOT(i)) <= 1.0e-10_dp) cycle
      sr  = real(ZROOT(i))
      sa2 = real(ZROOT(i)*conjg(ZROOT(i)))
      alpha = C*C - 2.0_dp*sr*C + sa2    ! α

      b1   = (-2.0_dp*C*C + 2.0_dp*sa2)/alpha
      b2   = ( C*C + 2.0_dp*sr*C + sa2)/alpha
      gsec = 1.0_dp/alpha

      sec = sec + 1
      H(4*sec-3) =  2.0_dp   ! alpha1
      H(4*sec-2) =  1.0_dp   ! alpha2
      H(4*sec-1) =  b1
      H(4*sec  ) =  b2
      GN = GN * gsec
    end do

    if (mod(N,2) == 1) then
      if (.not.has_real) then; IERR=6; return; end if
      sr = real(ZROOT(idx_real))
      if (abs(C-sr) < 1.0e-14_dp) then; IERR=7; return; end if

      b1   = (C + sr)/(C - sr)
      gsec = 1.0_dp/(C - sr)

      sec = sec + 1
      H(4*sec-3) = 1.0_dp
      H(4*sec-2) = 0.0_dp
      H(4*sec-1) = b1
      H(4*sec  ) = 0.0_dp
      GN = GN * gsec
    end if

    if (sec /= M) then
      IERR = 8
      return
    end if

    ! Final 1-point normalization at dp
    target = 1.0_dp / sqrt(1.0_dp + AP*AP)
    call FREQRESP_POINT(H, M, GN, DT, FPK, amp_dp, ph_dp)
    if (amp_dp <= 0.0_dp) then
      IERR = 9
      return
    end if
    GN = GN * (target / amp_dp)
  end subroutine BESLPN

  !======================================================================
  ! BESHPN (Bessel HPF design)
  !   - Uses HPF transform x = C * cot(ωT/2) -> C = xp * tan(ωp*T/2)
  !   - Final 1-point normalization at dp to meet |H(dp)|=1/sqrt(1+AP^2)
  !======================================================================
  subroutine BESHPN(H, M, GN, N, FPK, AP, DT, IERR)
    real(dp), intent(out), allocatable :: H(:) !!output: H(4 * M)
    real(dp), intent(out) :: GN
    integer , intent(out) :: M, IERR           !!M: filter order
    integer , intent(in)  :: N                 !!N: order of polynominal
    real(dp), intent(in)  :: FPK, AP, DT

    real(dp) :: COEF(0:50)
    complex(dp) :: ZROOT(50)
    integer :: MAXDBS, i, sec, idx_real
    real(dp) :: wp2, tp, TX, C, sr, sa2, alpha, b1, b2, gsec
    real(dp) :: amp_dp, ph_dp, target
    logical :: has_real

    IERR=0; GN=1.0_dp; M=0
    if (N<=0 .or. N>20) then; IERR=1; return; end if
    if (DT<=0.0_dp) then; IERR=2; return; end if
    if (FPK<=0.0_dp .or. FPK>=0.5_dp/DT) then; IERR=3; return; end if
    if (AP<0.0_dp) then; IERR=4; return; end if

    wp2 = PI*FPK*DT
    if (wp2<=0.0_dp .or. wp2>=PI/2.0_dp) then; IERR=5; return; end if
    tp  = tan(wp2)

    call GETBPL(N, COEF, MAXDBS, IERR); if (IERR/=0) return
    call DETBTX(AP, COEF, MAXDBS, TX, IERR); if (IERR/=0) return

    C = TX * tp                    ! HPF: C = xp * tan(ωp*T/2)  (important!)

    call BAIRSTOW(COEF, MAXDBS, ZROOT, IERR); if (IERR/=0) return

    M = N/2
    if (mod(N,2) == 1) M = M + 1
    allocate(H(4*M))

    has_real = .false.; idx_real = -1
    do i=1,N
      if (abs(aimag(ZROOT(i))) < 1.0e-10_dp) then
        has_real = .true.; idx_real=i; exit
      end if
    end do

    sec = 0
    do i=1,N
      if (i == idx_real) cycle
      if (aimag(ZROOT(i)) <= 1.0e-10_dp) cycle
      sr  = real(ZROOT(i))
      sa2 = real(ZROOT(i)*conjg(ZROOT(i)))
      alpha = C*C - 2.0_dp*sr*C + sa2

      b1   = ( 2.0_dp*C*C - 2.0_dp*sa2)/alpha
      b2   = ( C*C + 2.0_dp*sr*C + sa2)/alpha
      gsec = (-sa2)/alpha

      sec = sec + 1
      H(4*sec-3) = -2.0_dp
      H(4*sec-2) =  1.0_dp
      H(4*sec-1) =  b1
      H(4*sec  ) =  b2
      GN = GN * gsec
    end do

    if (mod(N,2) == 1) then
      if (.not.has_real) then; IERR=6; return; end if
      sr = real(ZROOT(idx_real))
      if (abs(C-sr) < 1.0e-14_dp) then; IERR=7; return; end if

      b1   = (C + sr)/(C - sr)
      gsec = (-sr)/(C - sr)

      sec = sec + 1
      H(4*sec-3) = -1.0_dp
      H(4*sec-2) =  0.0_dp
      H(4*sec-1) =  b1
      H(4*sec  ) =  0.0_dp
      GN = GN * gsec
    end if

    if (sec /= M) then
      IERR = 8
      return
    end if

    ! Final 1-point normalization at dp
    target = 1.0_dp / sqrt(1.0_dp + AP*AP)
    call FREQRESP_POINT(H, M, GN, DT, FPK, amp_dp, ph_dp)
    if (amp_dp <= 0.0_dp) then
      IERR = 9
      return
    end if
    GN = GN * (target / amp_dp)
  end subroutine BESHPN

  !======================================================================
  ! BESPSN (Bessel BPF design)
  !   - Analog prewarp: Ωa = (2/T) tan(ωT/2)
  !   - LP poles s_i -> analog BP poles p: s^2 - B s_i s + Ω0^2 = 0
  !   - Bilinear: z = (1 + pT/2) / (1 - pT/2)
  !   - Numerator per biquad: (1 - z^-2) => alpha1=0, alpha2=-1
  !   - Denominator from conjugate (or real-real) pairs
  !   - 1-point normalization at digital center (A-variable geometric mean)
  !======================================================================
  subroutine BESPSN(H, M, GN, N, FLK, FHK, AP, DT, IERR)
    real(dp), intent(out), allocatable :: H(:)  !!output: H(4 * M)
    real(dp), intent(out) :: GN
    integer , intent(out) :: M, IERR            !!M: filter order
    integer , intent(in)  :: N                  !!N: order of polynominal
    real(dp), intent(in)  :: FLK, FHK, AP, DT

    real(dp) :: COEF(0:50)
    complex(dp) :: SROOT(50)
    complex(dp), allocatable :: PBP(:)    ! analog BP poles
    complex(dp), allocatable :: ZPOLE(:)  ! digital poles
    integer :: MAXDBS, i, k, sec
    real(dp) :: wL, wH, TL, TH, Om_aL, Om_aH, B, Om0
    real(dp) :: amp0, ph0, t0, omega0, f0
    complex(dp) :: disc, p1, p2, z1, z2
    real(dp) :: b1, b2
    real(dp) :: tol
    integer :: nZ, ic
    complex(dp), allocatable :: pos(:), reals(:)

    IERR=0; GN=1.0_dp; M=0
    if (N<=0 .or. N>20) then; IERR=1; return; end if
    if (DT<=0.0_dp) then; IERR=2; return; end if
    if (FLK<=0.0_dp .or. FHK<=0.0_dp .or. FLK>=FHK) then; IERR=3; return; end if
    if (FHK >= 0.5_dp/DT) then; IERR=4; return; end if
    if (AP < 0.0_dp) then; IERR=5; return; end if

    ! Digital angular freqs and prewarp to analog (bilinear)
    wL = 2.0_dp*PI*FLK
    wH = 2.0_dp*PI*FHK
    TL    = tan(0.5_dp*wL*DT)
    TH    = tan(0.5_dp*wH*DT)
    Om_aL = (2.0_dp/DT) * TL
    Om_aH = (2.0_dp/DT) * TH
    if (Om_aL<=0.0_dp .or. Om_aH<=Om_aL) then; IERR=6; return; end if

    Om0 = sqrt(Om_aL*Om_aH)
    B   = Om_aH - Om_aL

    ! LP prototype (ascending coeffs) and roots
    call GETBPL(N, COEF, MAXDBS, IERR); if (IERR/=0) return
    call BAIRSTOW(COEF, MAXDBS, SROOT, IERR); if (IERR/=0) return

    ! Build 2N analog BP poles per prototype pole
    allocate(PBP(2*N))
    k = 0
    do i=1, N
      disc = sqrt( (B*SROOT(i))*(B*SROOT(i)) - cmplx(4.0_dp*Om0*Om0, 0.0_dp, kind=dp) )
      p1   = 0.5_dp*( B*SROOT(i) + disc )
      p2   = 0.5_dp*( B*SROOT(i) - disc )
      k = k+1; PBP(k)=p1
      k = k+1; PBP(k)=p2
    end do

    ! Map analog poles -> digital poles (bilinear)
    allocate(ZPOLE(2*N))
    do i=1,2*N
      if (abs(1.0_dp - 0.5_dp*PBP(i)*DT) < 1.0e-14_dp) then
        IERR = 7; return
      end if
      ZPOLE(i) = ( 1.0_dp + 0.5_dp*PBP(i)*DT ) / ( 1.0_dp - 0.5_dp*PBP(i)*DT )
    end do

    ! Classification by imaginary part to build stable biquads
    call sort_roots_by_imag_desc(ZPOLE)

    M   = N
    sec = 0
    tol = 1.0e-12_dp
    nZ  = 2*N
    allocate(h(4 * M))

    ! Extract positive-imag poles
    ic = 0
    do i=1, nZ
      if (aimag(ZPOLE(i)) >  tol) ic = ic + 1
    end do
    allocate(pos(ic))
    ic = 0
    do i=1, nZ
      if (aimag(ZPOLE(i)) >  tol) then
        ic = ic + 1
        pos(ic) = ZPOLE(i)
      end if
    end do

    ! Extract (near-)real poles
    ic = 0
    do i=1, nZ
      if (abs(aimag(ZPOLE(i))) <= tol) ic = ic + 1
    end do
    allocate(reals(ic))
    ic = 0
    do i=1, nZ
      if (abs(aimag(ZPOLE(i))) <= tol) then
        ic = ic + 1
        reals(ic) = ZPOLE(i)
      end if
    end do

    ! Sanity: number of sections must be N
    if ( size(pos) + (size(reals)/2) /= N ) then
      IERR = 91
      return
    end if
    if (mod(size(reals),2) /= 0) then
      IERR = 92
      return
    end if

    ! Biquads from positive-imag poles (pair with conjugate)
    do i=1, size(pos)
      z1 = pos(i)
      z2 = conjg(z1)
      b1 = - real(z1 + z2, kind=dp)   ! = -2*Re(z1)
      b2 =   real(z1 * z2, kind=dp)   ! = |z1|^2
      sec = sec + 1
      H(4*sec-3) =  0.0_dp
      H(4*sec-2) = -1.0_dp           ! (1 - z^-2)
      H(4*sec-1) =  b1
      H(4*sec  ) =  b2
    end do

    ! Biquads from real-real pairs
    do i=1, size(reals), 2
      z1 = reals(i)
      z2 = reals(i+1)
      b1 = - real(z1 + z2, kind=dp)
      b2 =   real(z1 * z2, kind=dp)
      sec = sec + 1
      H(4*sec-3) =  0.0_dp
      H(4*sec-2) = -1.0_dp
      H(4*sec-1) =  b1
      H(4*sec  ) =  b2
    end do

    if (sec /= M) then
      IERR = 93
      return
    end if

    ! 1-point normalization at the digital center:
    ! tan(ω0*T/2) = sqrt( tan(ωL*T/2) * tan(ωH*T/2) )
    t0     = sqrt( TL * TH )
    omega0 = (2.0_dp/DT) * atan( t0 )
    f0     = omega0 / (2.0_dp*PI)

    call FREQRESP_POINT(H, M, 1.0_dp, DT, f0, amp0, ph0)
    if (amp0 <= 0.0_dp) then
      IERR = 94
      return
    end if
    GN = 1.0_dp / amp0
  end subroutine BESPSN

end module besselfilter_design
