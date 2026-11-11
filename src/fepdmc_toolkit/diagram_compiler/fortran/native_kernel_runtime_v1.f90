! native_kernel_v1 primitive runtime. Versioned; do not change semantics in place.
! Formulas copied from Task-4 native_interpret.py. Closed-form 2x2 real-symmetric
! Hermitian exponential via H=mI+K, K^2=r^2 I (no LAPACK, no eigenvectors).
! Matrix layout: mathematical (i,j) stored as a(i,j), Fortran column-major.
module native_kernel_runtime_v1
  use, intrinsic :: iso_fortran_env, only: real64, int64
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer(int64), parameter :: NK_OK = 0_int64
  integer(int64), parameter :: NK_INVALID_EXACT_TARGET = 1_int64
  integer(int64), parameter :: NK_INVALID_CHEAP_WEIGHT = 2_int64
  integer(int64), parameter :: NK_INVALID_PROPOSAL_RATIO = 3_int64
  integer(int64), parameter :: NK_INVALID_RANDOM_UNIFORM = 4_int64
  integer(int64), parameter :: NK_CACHE_KEY_MISMATCH = 5_int64
  integer(int64), parameter :: NK_INVALID_IR = 6_int64
  character(len=*), parameter :: NK_RUNTIME_NAME = "native_kernel_v1"
  integer(int64) :: nk_n_electron_twoband = 0_int64
  integer(int64) :: nk_n_vertex_sigmaz = 0_int64
  integer(int64) :: nk_n_matmul = 0_int64
  integer(int64) :: nk_n_trace = 0_int64
  integer(int64) :: nk_n_exact_graph = 0_int64
  integer(int64) :: nk_n_cheap_graph = 0_int64
contains

  subroutine nk_reset_counters()
    nk_n_electron_twoband = 0_int64
    nk_n_vertex_sigmaz = 0_int64
    nk_n_matmul = 0_int64
    nk_n_trace = 0_int64
    nk_n_exact_graph = 0_int64
    nk_n_cheap_graph = 0_int64
  end subroutine nk_reset_counters

  subroutine nk_bump_exact_graph()
    nk_n_exact_graph = nk_n_exact_graph + 1_int64
  end subroutine nk_bump_exact_graph

  subroutine nk_bump_cheap_graph()
    nk_n_cheap_graph = nk_n_cheap_graph + 1_int64
  end subroutine nk_bump_cheap_graph

  pure function nk_xi(k_eff, t) result(xi)
    real(real64), intent(in) :: k_eff, t
    real(real64) :: xi
    xi = 2.0_real64 * t * (1.0_real64 - cos(k_eff))
  end function nk_xi

  pure function nk_electron_scalar(k_eff, t, dtau) result(z)
    real(real64), intent(in) :: k_eff, t, dtau
    complex(real64) :: z
    z = cmplx(exp(-nk_xi(k_eff, t) * dtau), 0.0_real64, kind=real64)
  end function nk_electron_scalar

  subroutine nk_electron_twoband(k_eff, t, dtau, delta, gap, u)
    ! exp(-dtau H) for real-symmetric H=[[xi,delta],[delta,xi+gap]].
    ! Uses H = m I + K, K^2 = r^2 I, so no eigenvectors (stable for tiny delta).
    real(real64), intent(in) :: k_eff, t, dtau, delta, gap
    complex(real64), intent(out) :: u(2, 2)
    real(real64) :: xi, a, b, c, m, d, r, scale, ch, sh_over_r
    xi = nk_xi(k_eff, t)
    a = xi
    b = xi + gap
    c = delta
    m = 0.5_real64 * (a + b)
    d = 0.5_real64 * (a - b)
    r = sqrt(d * d + c * c)
    scale = exp(-dtau * m)
    if (r == 0.0_real64) then
      u = cmplx(0.0_real64, 0.0_real64, kind=real64)
      u(1, 1) = cmplx(scale, 0.0_real64, kind=real64)
      u(2, 2) = cmplx(scale, 0.0_real64, kind=real64)
    else
      ch = cosh(dtau * r)
      sh_over_r = sinh(dtau * r) / r
      u(1, 1) = cmplx(scale * (ch - sh_over_r * d), 0.0_real64, kind=real64)
      u(1, 2) = cmplx(scale * (-sh_over_r * c), 0.0_real64, kind=real64)
      u(2, 1) = cmplx(scale * (-sh_over_r * c), 0.0_real64, kind=real64)
      u(2, 2) = cmplx(scale * (ch + sh_over_r * d), 0.0_real64, kind=real64)
    end if
    nk_n_electron_twoband = nk_n_electron_twoband + 1_int64
  end subroutine nk_electron_twoband

  pure function nk_phonon(omega, dtau) result(z)
    real(real64), intent(in) :: omega, dtau
    complex(real64) :: z
    z = cmplx(exp(-omega * dtau), 0.0_real64, kind=real64)
  end function nk_phonon

  pure function nk_prefactor(g, L, n) result(z)
    real(real64), intent(in) :: g, L
    integer(int64), intent(in) :: n
    complex(real64) :: z
    z = cmplx((g / sqrt(L)) ** (2_int64 * n), 0.0_real64, kind=real64)
  end function nk_prefactor

  subroutine nk_vertex_sigmaz(g, L, m)
    real(real64), intent(in) :: g, L
    complex(real64), intent(out) :: m(2, 2)
    real(real64) :: s
    s = g / sqrt(L)
    m = cmplx(0.0_real64, 0.0_real64, kind=real64)
    m(1, 1) = cmplx(s, 0.0_real64, kind=real64)
    m(2, 2) = cmplx(-s, 0.0_real64, kind=real64)
    nk_n_vertex_sigmaz = nk_n_vertex_sigmaz + 1_int64
  end subroutine nk_vertex_sigmaz

  subroutine nk_matmul(a, b, c)
    complex(real64), intent(in) :: a(2, 2), b(2, 2)
    complex(real64), intent(out) :: c(2, 2)
    c = matmul(a, b)
    nk_n_matmul = nk_n_matmul + 1_int64
  end subroutine nk_matmul

  subroutine nk_scale_m2(a, s, c)
    complex(real64), intent(in) :: a(2, 2)
    complex(real64), intent(in) :: s
    complex(real64), intent(out) :: c(2, 2)
    c = s * a
  end subroutine nk_scale_m2

  subroutine nk_add_m2(a, b, c)
    complex(real64), intent(in) :: a(2, 2), b(2, 2)
    complex(real64), intent(out) :: c(2, 2)
    c = a + b
  end subroutine nk_add_m2

  pure function nk_trace(a) result(z)
    complex(real64), intent(in) :: a(2, 2)
    complex(real64) :: z
    z = a(1, 1) + a(2, 2)
  end function nk_trace

  subroutine nk_trace_counted(a, z)
    complex(real64), intent(in) :: a(2, 2)
    complex(real64), intent(out) :: z
    z = nk_trace(a)
    nk_n_trace = nk_n_trace + 1_int64
  end subroutine nk_trace_counted

  subroutine nk_log_positive_real(z, imag_atol, logw, status, fail_code)
    complex(real64), intent(in) :: z
    real(real64), intent(in) :: imag_atol
    real(real64), intent(out) :: logw
    integer(int64), intent(inout) :: status
    integer(int64), intent(in) :: fail_code
    real(real64) :: re, im
    re = real(z, real64)
    im = aimag(z)
    if (.not. ieee_is_finite(re) .or. .not. ieee_is_finite(im)) then
      status = fail_code
      logw = 0.0_real64
      return
    end if
    if (abs(im) > imag_atol .or. re <= 0.0_real64) then
      status = fail_code
      logw = 0.0_real64
      return
    end if
    logw = log(re)
  end subroutine nk_log_positive_real

  subroutine nk_log_uniform(u, floor, logu, status)
    real(real64), intent(in) :: u, floor
    real(real64), intent(out) :: logu
    integer(int64), intent(inout) :: status
    if (.not. ieee_is_finite(u) .or. u < 0.0_real64 .or. u > 1.0_real64) then
      status = NK_INVALID_RANDOM_UNIFORM
      logu = 0.0_real64
      return
    end if
    logu = log(max(u, floor))
  end subroutine nk_log_uniform

  subroutine nk_check_log_q(value, kind_sym, logq, status)
    real(real64), intent(in) :: value
    integer(int64), intent(in) :: kind_sym
    real(real64), intent(out) :: logq
    integer(int64), intent(inout) :: status
    if (.not. ieee_is_finite(value)) then
      status = NK_INVALID_PROPOSAL_RATIO
      logq = 0.0_real64
      return
    end if
    if (kind_sym == 1_int64 .and. value /= 0.0_real64) then
      status = NK_INVALID_PROPOSAL_RATIO
      logq = 0.0_real64
      return
    end if
    logq = value
  end subroutine nk_check_log_q

end module native_kernel_runtime_v1
