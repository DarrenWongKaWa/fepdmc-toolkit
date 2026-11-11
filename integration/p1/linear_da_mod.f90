! Public delayed-acceptance module for pinned FEP-DMC JJ update_swap.
! Architecture C: ell_R = log(P_accept) from native P_accept.
! LINEAR_DA=off preserves the vanilla swap kernel and native seed sequence.
module linear_da_mod
  use, intrinsic :: ieee_arithmetic
  use DiagMC, only: fynman, dmc_band, cal_ek_int, Gel, Dph
  use pert_const, only: dp
  use pert_param, only: DMC_Method, zeroTMC, sample_gt
  implicit none
  character(len=16) :: linear_da = 'off'
  character(len=16) :: score_kind = 'prop'
  logical :: da_cfg_init = .false.
  integer(kind=8) :: aux_state = 1_8
  integer(kind=8) :: aux_seed_used = 0_8
  integer :: n_seen = 0, n_eligible = 0
  integer :: n_s1_reject = 0, n_expensive = 0, n_s2_reject = 0, n_accept = 0
  real(dp), parameter :: clip_log10 = log(10.0_dp)
contains

  subroutine linear_da_ensure()
    character(len=32) :: env
    integer :: st
    if (da_cfg_init) return
    linear_da = 'off'
    score_kind = 'prop'
    call get_environment_variable('LINEAR_DA', env, status=st)
    if (st == 0 .and. len_trim(env) > 0) linear_da = trim(adjustl(env))
    if (trim(linear_da) /= 'off' .and. trim(linear_da) /= 'on') then
       error stop 'linear_da: LINEAR_DA must be off or on'
    endif
    call get_environment_variable('LINEAR_DA_SCORE', env, status=st)
    if (st == 0 .and. len_trim(env) > 0) score_kind = trim(adjustl(env))
    if (trim(score_kind) /= 'prop') then
       error stop 'linear_da: LINEAR_DA_SCORE must be prop'
    endif
    call get_environment_variable('LINEAR_DA_AUX_SEED', env, status=st)
    if (st == 0 .and. len_trim(env) > 0) then
       read(env, *) aux_seed_used
    else
       aux_seed_used = 9142871_8
    endif
    if (aux_seed_used == 0_8) aux_seed_used = 1_8
    aux_state = aux_seed_used
    if (trim(linear_da) == 'on') then
       if (DMC_Method /= 0) error stop 'linear_da: DMC_Method must be 0'
       if (dmc_band /= 1) error stop 'linear_da: dmc_band must be 1'
       if (sample_gt) error stop 'linear_da: sample_gt must be false'
       if (.not. zeroTMC) error stop 'linear_da: zeroTMC must be true'
    endif
    da_cfg_init = .true.
  end subroutine linear_da_ensure

  logical function linear_da_is_on()
    call linear_da_ensure()
    linear_da_is_on = (trim(linear_da) == 'on')
  end function linear_da_is_on

  subroutine linear_da_aux_uniform(u)
    real(dp), intent(out) :: u
    integer(kind=8) :: x
    x = aux_state
    if (x == 0_8) x = 1_8
    x = ieor(x, ishft(x, 13))
    x = ieor(x, ishft(x, -7))
    x = ieor(x, ishft(x, 17))
    aux_state = x
    u = 0.5_dp + 0.5_dp * real(x, dp) / real(huge(x), dp)
    if (u <= 0.0_dp) u = 1.0e-16_dp
    if (u >= 1.0_dp) u = 1.0_dp - 1.0e-16_dp
  end subroutine linear_da_aux_uniform

  subroutine linear_da_eval(diagram, iv1, iv2, eligible, ell_hat)
    type(fynman), intent(in) :: diagram
    integer, intent(in) :: iv1, iv2
    logical, intent(out) :: eligible
    real(dp), intent(out) :: ell_hat
    integer :: kout(3)
    real(dp) :: ekv(dmc_band), tau12, w1, w2, tL, tR, tPL, tPR, pk, gnew, gold, d1n, d1d, d2n, d2d
    eligible = .false.
    ell_hat = 0.0_dp
    call linear_da_ensure()
    !$omp atomic update
    n_seen = n_seen + 1
    kout = diagram%vertexList(iv1)%i_kin + diagram%vertexList(iv2)%i_q
    call cal_ek_int(kout, ekv)
    if (.not. ieee_is_finite(minval(ekv))) then
       error stop 'linear_da: cheap energy is nonfinite'
    endif
    tL = diagram%vertexList(iv1)%tau
    tR = diagram%vertexList(iv2)%tau
    tau12 = tR - tL
    tPL = diagram%vertexList(diagram%vertexList(iv1)%link(2))%tau
    tPR = diagram%vertexList(diagram%vertexList(iv2)%link(2))%tau
    w1 = diagram%vertexList(iv1)%wq(diagram%vertexList(iv1)%nu)
    w2 = diagram%vertexList(iv2)%wq(diagram%vertexList(iv2)%nu)
    gnew = Gel(minval(ekv), tau12)
    gold = Gel(minval(diagram%vertexList(iv1)%ekout), tau12)
    d1n = Dph(w1, abs(tR - tPL))
    d1d = Dph(w1, abs(tL - tPL))
    d2n = Dph(w2, abs(tL - tPR))
    d2d = Dph(w2, abs(tR - tPR))
    if (gold <= 0.0_dp .or. d1d <= 0.0_dp .or. d2d <= 0.0_dp) then
       ell_hat = -clip_log10
    else
       pk = (gnew / gold) * (d1n / d1d) * (d2n / d2d)
       if (.not. ieee_is_finite(pk) .or. pk <= 0.0_dp) then
          ell_hat = -clip_log10
       else
          ell_hat = log(pk)
          if (ell_hat > clip_log10) ell_hat = clip_log10
          if (ell_hat < -clip_log10) ell_hat = -clip_log10
       endif
    endif
    eligible = .true.
    !$omp atomic update
    n_eligible = n_eligible + 1
  end subroutine linear_da_eval

  subroutine linear_da_stage2(ran, P_accept, ell_hat, da_acc)
    real(dp), intent(in) :: ran, P_accept, ell_hat
    logical, intent(out) :: da_acc
    real(dp) :: ell_R, u
    da_acc = .false.
    if (.not. ieee_is_finite(P_accept)) then
       error stop 'linear_da: P_accept is nonfinite'
    endif
    if (P_accept <= 0.0_dp) return
    ell_R = log(P_accept)
    u = ran
    if (.not. ieee_is_finite(u) .or. u <= 0.0_dp) u = 1.0e-300_dp
    da_acc = (log(max(u, 1.0e-300_dp)) < min(0.0_dp, ell_R - ell_hat))
  end subroutine linear_da_stage2

  subroutine linear_da_note_s1_reject()
    !$omp atomic update
    n_s1_reject = n_s1_reject + 1
  end subroutine linear_da_note_s1_reject

  subroutine linear_da_note_expensive()
    !$omp atomic update
    n_expensive = n_expensive + 1
  end subroutine linear_da_note_expensive

  subroutine linear_da_note_accept()
    !$omp atomic update
    n_accept = n_accept + 1
  end subroutine linear_da_note_accept

  subroutine linear_da_note_stage2_reject()
    !$omp atomic update
    n_s2_reject = n_s2_reject + 1
  end subroutine linear_da_note_stage2_reject

  subroutine linear_da_report()
    call linear_da_ensure()
    write(*,'(A,A,1X,A,6I12,1X,I18)') 'LINEAR_DA_COUNTS ', trim(linear_da), trim(score_kind), &
         n_seen, n_eligible, n_s1_reject, n_expensive, n_s2_reject, n_accept, aux_seed_used
  end subroutine linear_da_report
end module linear_da_mod
