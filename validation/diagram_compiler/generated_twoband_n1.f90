! backend=fortran_v1 primitive_set=native_kernel_v1
! schema_version=1 native_ir_digest=4932e8b066554e966b71efb9b54c69974e17dcca82868d8df511927d2f03be31
! ir_kind=diagram design=B cheap_policy=propagator_only_v1
! cache_contract=caller_guarantees_exact_x_valid_belongs_to_x
module generated_da_kernel
  use, intrinsic :: iso_fortran_env, only: real64, int64
  use native_kernel_runtime_v1
  implicit none
  character(len=*), parameter :: NK_BACKEND = 'fortran_v1'
  character(len=*), parameter :: NK_PRIMITIVE_SET = 'native_kernel_v1'
  character(len=*), parameter :: NK_IR_DIGEST = '4932e8b066554e966b71efb9b54c69974e17dcca82868d8df511927d2f03be31'
  integer, parameter :: NK_SCHEMA_VERSION = 1
contains
  subroutine eval_graph_cheap(&
      b_k, b_q, b_tau, b_t, b_omega, b_g, b_L, b_delta, b_gap, graph_value)
    real(real64), intent(in) :: b_k
    real(real64), intent(in) :: b_q(1)
    real(real64), intent(in) :: b_tau(2)
    real(real64), intent(in) :: b_t
    real(real64), intent(in) :: b_omega
    real(real64), intent(in) :: b_g
    integer(int64), intent(in) :: b_L
    real(real64), intent(in) :: b_delta
    real(real64), intent(in) :: b_gap
    complex(real64), intent(out) :: graph_value
    complex(real64) :: c1
    real(real64) :: c10
    real(real64) :: c11
    real(real64) :: c12
    real(real64) :: c13
    real(real64) :: c14
    real(real64) :: c15
    complex(real64) :: c16
    complex(real64) :: c17
    real(real64) :: c18
    real(real64) :: c19
    real(real64) :: c2
    real(real64) :: c20
    real(real64) :: c21
    complex(real64) :: c22
    complex(real64) :: c23
    complex(real64) :: c24
    real(real64) :: c3
    complex(real64) :: c4
    real(real64) :: c5
    real(real64) :: c6
    real(real64) :: c7
    real(real64) :: c8
    real(real64) :: c9
    call nk_bump_cheap_graph()
    c1 = cmplx(1.00000000000000000e+00_real64, 0.00000000000000000e+00_real64, kind=real64)
    c2 = b_g
    c3 = real(b_L, real64)
    c4 = nk_prefactor(c2, c3, 1_int64)
    c5 = b_k
    c6 = 1.00000000000000000e+00_real64
    c7 = c6 * c5
    c8 = b_q(1)
    c9 = -1.00000000000000000e+00_real64
    c10 = c9 * c8
    c11 = c7 + c10
    c12 = b_t
    c13 = b_tau(2)
    c14 = b_tau(1)
    c15 = c13 - c14
    c16 = nk_electron_scalar(c11, c12, c15)
    c17 = c1 * c16
    c18 = b_omega
    c19 = b_tau(2)
    c20 = b_tau(1)
    c21 = c19 - c20
    c22 = nk_phonon(c18, c21)
    c23 = c17 * c22
    c24 = c4 * c23
    graph_value = c24
  end subroutine eval_graph_cheap
  subroutine eval_graph_exact(&
      b_k, b_q, b_tau, b_t, b_omega, b_g, b_L, b_delta, b_gap, graph_value)
    real(real64), intent(in) :: b_k
    real(real64), intent(in) :: b_q(1)
    real(real64), intent(in) :: b_tau(2)
    real(real64), intent(in) :: b_t
    real(real64), intent(in) :: b_omega
    real(real64), intent(in) :: b_g
    integer(int64), intent(in) :: b_L
    real(real64), intent(in) :: b_delta
    real(real64), intent(in) :: b_gap
    complex(real64), intent(out) :: graph_value
    complex(real64) :: e1(2, 2)
    real(real64) :: e10
    real(real64) :: e11
    real(real64) :: e12
    real(real64) :: e13
    real(real64) :: e14
    real(real64) :: e15
    real(real64) :: e16
    real(real64) :: e17
    real(real64) :: e18
    real(real64) :: e19
    complex(real64) :: e2
    complex(real64) :: e20(2, 2)
    complex(real64) :: e21(2, 2)
    real(real64) :: e22
    real(real64) :: e23
    complex(real64) :: e24(2, 2)
    complex(real64) :: e25(2, 2)
    real(real64) :: e26
    real(real64) :: e27
    real(real64) :: e28
    real(real64) :: e29
    real(real64) :: e3
    complex(real64) :: e30
    complex(real64) :: e31
    complex(real64) :: e32(2, 2)
    complex(real64) :: e33
    real(real64) :: e4
    complex(real64) :: e5(2, 2)
    complex(real64) :: e6(2, 2)
    real(real64) :: e7
    real(real64) :: e8
    real(real64) :: e9
    call nk_bump_exact_graph()
    e1(1, 1) = cmplx(1.00000000000000000e+00_real64, 0.00000000000000000e+00_real64, kind=real64)
    e1(1, 2) = cmplx(0.00000000000000000e+00_real64, 0.00000000000000000e+00_real64, kind=real64)
    e1(2, 1) = cmplx(0.00000000000000000e+00_real64, 0.00000000000000000e+00_real64, kind=real64)
    e1(2, 2) = cmplx(1.00000000000000000e+00_real64, 0.00000000000000000e+00_real64, kind=real64)
    e2 = cmplx(1.00000000000000000e+00_real64, 0.00000000000000000e+00_real64, kind=real64)
    e3 = b_g
    e4 = real(b_L, real64)
    call nk_vertex_sigmaz(e3, e4, e5)
    call nk_matmul(e5, e1, e6)
    e7 = b_k
    e8 = 1.00000000000000000e+00_real64
    e9 = e8 * e7
    e10 = b_q(1)
    e11 = -1.00000000000000000e+00_real64
    e12 = e11 * e10
    e13 = e9 + e12
    e14 = b_t
    e15 = b_tau(2)
    e16 = b_tau(1)
    e17 = e15 - e16
    e18 = b_delta
    e19 = b_gap
    call nk_electron_twoband(e13, e14, e17, e18, e19, e20)
    call nk_matmul(e20, e6, e21)
    e22 = b_g
    e23 = real(b_L, real64)
    call nk_vertex_sigmaz(e22, e23, e24)
    call nk_matmul(e24, e21, e25)
    e26 = b_omega
    e27 = b_tau(2)
    e28 = b_tau(1)
    e29 = e27 - e28
    e30 = nk_phonon(e26, e29)
    e31 = e2 * e30
    call nk_scale_m2(e25, e31, e32)
    call nk_trace_counted(e32, e33)
    graph_value = e33
  end subroutine eval_graph_exact
  subroutine evaluate_da_kernel(&
      x_k, x_q, x_tau, x_t, x_omega, x_g, x_L, x_delta, x_gap, y_k, y_q, y_tau, y_t, &
        y_omega, y_g, y_L, y_delta, y_gap, log_q_reverse_minus_forward, u1, u2, &
        exact_x_valid, exact_log_weight_x, accepted, reject_stage, ell_hat, &
        ell_R_valid, ell_R, exact_y_evaluated, status)
    real(real64), intent(in) :: x_k
    real(real64), intent(in) :: x_q(1)
    real(real64), intent(in) :: x_tau(2)
    real(real64), intent(in) :: x_t
    real(real64), intent(in) :: x_omega
    real(real64), intent(in) :: x_g
    integer(int64), intent(in) :: x_L
    real(real64), intent(in) :: x_delta
    real(real64), intent(in) :: x_gap
    real(real64), intent(in) :: y_k
    real(real64), intent(in) :: y_q(1)
    real(real64), intent(in) :: y_tau(2)
    real(real64), intent(in) :: y_t
    real(real64), intent(in) :: y_omega
    real(real64), intent(in) :: y_g
    integer(int64), intent(in) :: y_L
    real(real64), intent(in) :: y_delta
    real(real64), intent(in) :: y_gap
    real(real64), intent(in) :: log_q_reverse_minus_forward, u1, u2
    logical, intent(in) :: exact_x_valid
    real(real64), intent(in) :: exact_log_weight_x
    logical, intent(out) :: accepted, ell_R_valid, exact_y_evaluated
    integer(int64), intent(out) :: reject_stage, status
    real(real64), intent(out) :: ell_hat, ell_R
    real(real64) :: cached_px
    complex(real64) :: cheap_x
    complex(real64) :: cheap_y
    real(real64) :: computed_px
    real(real64) :: delta
    real(real64) :: diff_pi
    complex(real64) :: exact_x
    complex(real64) :: exact_y
    real(real64) :: log_px
    real(real64) :: log_py
    real(real64) :: log_q
    real(real64) :: log_u1
    real(real64) :: log_u2
    real(real64) :: log_wx
    real(real64) :: log_wy
    real(real64) :: min0_delta
    real(real64) :: min0_hat
    logical :: stage1_pass
    logical :: stage2_pass
    real(real64) :: zero
    call nk_reset_counters()
    status = NK_OK
    accepted = .false.
    reject_stage = 0_int64
    ell_R_valid = .false.
    exact_y_evaluated = .false.
    ell_hat = 0.0_real64
    ell_R = 0.0_real64
    goto 100
    ! block entry
100 continue
    call eval_graph_cheap(x_k, x_q, x_tau, x_t, x_omega, x_g, x_L, x_delta, x_gap, cheap_x)
    call eval_graph_cheap(y_k, y_q, y_tau, y_t, y_omega, y_g, y_L, y_delta, y_gap, cheap_y)
    call nk_log_positive_real(cheap_x, 9.99999999999999980e-13_real64, log_wx, status, NK_INVALID_CHEAP_WEIGHT)
    if (status /= NK_OK) goto 900
    call nk_log_positive_real(cheap_y, 9.99999999999999980e-13_real64, log_wy, status, NK_INVALID_CHEAP_WEIGHT)
    if (status /= NK_OK) goto 900
    ell_hat = log_wy - log_wx
    call nk_check_log_q(log_q_reverse_minus_forward, 0_int64, log_q, status)
    if (status /= NK_OK) goto 900
    call nk_log_uniform(u1, 1.00000000000000003e-300_real64, log_u1, status)
    if (status /= NK_OK) goto 900
    zero = 0.00000000000000000e+00_real64
    min0_hat = min(zero, ell_hat)
    stage1_pass = (log_u1 < min0_hat)
    if (stage1_pass) then
      goto 300
    else
      goto 200
    end if
    ! block reject_stage1
200 continue
    accepted = .false.
    reject_stage = 1_int64
    exact_y_evaluated = .false.
    ell_R_valid = .false.
    goto 900
    ! block exact_candidate
300 continue
    if (exact_x_valid) then
      goto 310
    else
      goto 320
    end if
    ! block use_cache
310 continue
    cached_px = exact_log_weight_x
    log_px = cached_px
    goto 330
    ! block compute_exact_x
320 continue
    call eval_graph_exact(x_k, x_q, x_tau, x_t, x_omega, x_g, x_L, x_delta, x_gap, exact_x)
    call nk_log_positive_real(exact_x, 9.99999999999999980e-13_real64, computed_px, status, NK_INVALID_EXACT_TARGET)
    if (status /= NK_OK) goto 900
    log_px = computed_px
    goto 330
    ! block after_exact_x
330 continue
    call eval_graph_exact(y_k, y_q, y_tau, y_t, y_omega, y_g, y_L, y_delta, y_gap, exact_y)
    call nk_log_positive_real(exact_y, 9.99999999999999980e-13_real64, log_py, status, NK_INVALID_EXACT_TARGET)
    if (status /= NK_OK) goto 900
    diff_pi = log_py - log_px
    ell_R = diff_pi + log_q
    delta = ell_R - ell_hat
    call nk_log_uniform(u2, 1.00000000000000003e-300_real64, log_u2, status)
    if (status /= NK_OK) goto 900
    min0_delta = min(zero, delta)
    stage2_pass = (log_u2 < min0_delta)
    if (stage2_pass) then
      goto 500
    else
      goto 400
    end if
    ! block reject_stage2
400 continue
    accepted = .false.
    reject_stage = 2_int64
    exact_y_evaluated = .true.
    ell_R_valid = .true.
    goto 900
    ! block accept
500 continue
    accepted = .true.
    reject_stage = 0_int64
    exact_y_evaluated = .true.
    ell_R_valid = .true.
    goto 900
900 continue
    if (status /= NK_OK) then
      accepted = .false.
      reject_stage = 0_int64
      ell_R_valid = .false.
      exact_y_evaluated = .false.
    end if
  end subroutine evaluate_da_kernel
end module generated_da_kernel
