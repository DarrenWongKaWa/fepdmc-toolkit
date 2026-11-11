! Delayed-acceptance helpers for update_swap. Default linear_da='off'
! reproduces the original single-stage path. Do not call Python per move.
module linear_da_mod
  use DiagMC, only: fynman, vertex, maxN, config, cal_ek_int, dmc_band
  use pert_const, only: dp
  use pert_param, only: phfreq_cutoff, DMC_Method
  implicit none
  character(len=16) :: linear_da = 'off'
  character(len=16) :: score_kind = 'prop'
  logical :: da_cfg_init = .false.
  logical :: linear_model_loaded = .false.
  integer(kind=8) :: aux_state = 1_8
  integer(kind=8) :: aux_seed_used = 0_8
  integer :: n_invoked = 0, n_eligible = 0, n_score = 0
  integer :: n_stage1_rej = 0, n_stage2 = 0, n_accept = 0
  real(dp), parameter :: log_dph_cut = log(1.0e-15_dp)
  real(dp), parameter :: clip_log10 = log(10.0_dp)
  integer, parameter :: ngrid_da = 20
contains

  subroutine linear_da_ensure()
    character(len=32) :: env
    integer :: st
    if (da_cfg_init) return
    call get_environment_variable('LINEAR_DA', env, status=st)
    if (st == 0 .and. len_trim(env) > 0) linear_da = trim(adjustl(env))
    call get_environment_variable('LINEAR_DA_SCORE', env, status=st)
    if (st == 0 .and. len_trim(env) > 0) score_kind = trim(adjustl(env))
    call get_environment_variable('LINEAR_DA_AUX_SEED', env, status=st)
    if (st == 0 .and. len_trim(env) > 0) then
       read(env, *) aux_seed_used
    else
       aux_seed_used = 9142871_8
    endif
    if (aux_seed_used == 0_8) aux_seed_used = 1_8
    aux_state = aux_seed_used
    if ((trim(score_kind) == 'linear' .or. trim(score_kind) == 'lr') .and. &
         .not. linear_model_loaded) then
       write(*,*) 'MODEL_NOT_READY score_kind=', trim(score_kind)
       error stop 'MODEL_NOT_READY'
    endif
    da_cfg_init = .true.
  end subroutine linear_da_ensure

  logical function linear_da_is_on()
    call linear_da_ensure()
    linear_da_is_on = (trim(linear_da) == 'on')
  end function linear_da_is_on

  logical function linear_da_is_shadow()
    call linear_da_ensure()
    linear_da_is_shadow = (trim(linear_da) == 'shadow')
  end function linear_da_is_shadow

  logical function p1_fixture_on()
    character(len=8) :: env
    integer :: st
    p1_fixture_on = .false.
    call get_environment_variable('P1_FIXTURE', env, status=st)
    if (st==0) then
      if (trim(adjustl(env))=='1' .or. trim(adjustl(env))=='on') p1_fixture_on = .true.
    endif
  end function p1_fixture_on

  logical function linear_da_hook()
    call linear_da_ensure()
    linear_da_hook = (trim(linear_da) == 'on' .or. trim(linear_da) == 'shadow')
  end function linear_da_hook

  real(dp) function linear_ld(w, d)
    real(dp), intent(in) :: w, d
    real(dp) :: cut
    cut = phfreq_cutoff
    ! vertex wq already eV; native Dph compares w < phfreq_cutoff*ryd2ev
    ! after pert_param converted the namelist meV cutoff into Ry.
    ! Python offline uses 1e-3 eV. Here keep native Dph threshold via
    ! the same comparison as diagMC Dph: w < phfreq_cutoff*ryd2ev is
    ! applied in Dph itself. For log form use w in eV vs 1e-3.
    if (w < 1.0e-3_dp) then
       linear_ld = log_dph_cut
    else
       linear_ld = -w * abs(d)
    endif
  end function linear_ld

  subroutine linear_da_aux_uniform(u)
    real(dp), intent(out) :: u
    integer(kind=8) :: x
    ! Independent xorshift64; does not consume native diagram%seed.
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

  subroutine extract_prefix4_ids(diagram, ids, nfound)
    type(fynman), intent(in) :: diagram
    integer, intent(out) :: ids(4), nfound
    integer :: iv
    ids = 0
    nfound = 0
    iv = diagram%vertexList(1)%link(3)
    do while (iv /= maxN .and. nfound < 4 .and. iv > 0)
       nfound = nfound + 1
       ids(nfound) = iv
       iv = diagram%vertexList(iv)%link(3)
    end do
  end subroutine extract_prefix4_ids

  logical function ids_closed(diagram, ids, nfound)
    type(fynman), intent(in) :: diagram
    integer, intent(in) :: ids(4), nfound
    integer :: i, p, j
    logical :: inside
    ids_closed = .false.
    if (nfound < 4) return
    if (diagram%nph_ext /= 0) return
    do i = 1, 4
       p = diagram%vertexList(ids(i))%link(2)
       if (p == 1 .or. p == maxN) return
       inside = .false.
       do j = 1, 4
          if (ids(j) == p) inside = .true.
       end do
       if (.not. inside) return
    end do
    ids_closed = .true.
  end function ids_closed

  logical function prefix_gamma(diagram, ids, nfound)
    type(fynman), intent(in) :: diagram
    integer, intent(in) :: ids(4), nfound
    integer :: k(3)
    prefix_gamma = .false.
    if (nfound < 4) return
    k = modulo(diagram%vertexList(ids(1))%i_kin, ngrid_da)
    if (maxval(abs(k)) /= 0) return
    k = modulo(diagram%vertexList(ids(4))%i_kout, ngrid_da)
    if (maxval(abs(k)) /= 0) return
    prefix_gamma = .true.
  end function prefix_gamma

  subroutine overlay_ids_walk(diagram, iv1, iv2, ids, nfound)
    type(fynman), intent(in) :: diagram
    integer, intent(in) :: iv1, iv2
    integer, intent(out) :: ids(4), nfound
    integer :: iv, vLL, vRR, nxt
    vLL = diagram%vertexList(iv1)%link(1)
    vRR = diagram%vertexList(iv2)%link(3)
    ids = 0
    nfound = 0
    iv = 1
    do while (nfound < 4)
       if (iv == vLL) then
          nxt = iv2
       else if (iv == iv2) then
          nxt = iv1
       else if (iv == iv1) then
          nxt = vRR
       else
          nxt = diagram%vertexList(iv)%link(3)
       endif
       if (nxt == maxN .or. nxt <= 0) exit
       nfound = nfound + 1
       ids(nfound) = nxt
       iv = nxt
    end do
  end subroutine overlay_ids_walk

  logical function pair_in_ids(ids, nfound, iv1, iv2)
    integer, intent(in) :: ids(4), nfound, iv1, iv2
    integer :: i
    logical :: a, b
    a = .false.
    b = .false.
    do i = 1, nfound
       if (ids(i) == iv1) a = .true.
       if (ids(i) == iv2) b = .true.
    end do
    pair_in_ids = a .and. b
  end function pair_in_ids

  subroutine linear_da_eval(diagram, iv1, iv2, eligible, ell_hat)
    type(fynman), intent(in) :: diagram
    integer, intent(in) :: iv1, iv2
    logical, intent(out) :: eligible
    real(dp), intent(out) :: ell_hat
    integer :: idsX(4), idsY(4), nX, nY
    real(dp) :: lx, ly, d
    eligible = .false.
    ell_hat = 0.0_dp
    call linear_da_ensure()
    n_invoked = n_invoked + 1
    if (diagram%order <= 3) return
    if (DMC_Method /= 0) return
    if (dmc_band /= 1) return
    if (trim(score_kind) == 'prop') then
       eligible = .true.
       n_eligible = n_eligible + 1
       n_score = n_score + 1
       call p1_prop_ell(diagram, iv1, iv2, ell_hat)
       return
    endif
    if (diagram%nph_ext /= 0) return
    call extract_prefix4_ids(diagram, idsX, nX)
    call overlay_ids_walk(diagram, iv1, iv2, idsY, nY)
    if (.not. ids_closed(diagram, idsX, nX)) return
    if (.not. ids_closed(diagram, idsY, nY)) return
    if (.not. prefix_gamma(diagram, idsX, nX)) return
    if (.not. pair_in_ids(idsX, nX, iv1, iv2)) return
    if (.not. pair_in_ids(idsY, nY, iv1, iv2)) return
    eligible = .true.
    n_eligible = n_eligible + 1
    n_score = n_score + 1
    if (trim(score_kind) /= 'prop') then
       ! Must have aborted at init if linear/lr requested without a model.
       if (.not. linear_weights_ready()) then
          write(*,*) 'MODEL_NOT_READY during eval'
          error stop 'MODEL_NOT_READY'
       endif
    endif
    lx = prefix_logp(diagram, idsX, nX, .false., iv1, iv2)
    ly = prefix_logp(diagram, idsY, nY, .true., iv1, iv2)
    d = ly - lx
    if (d > clip_log10) d = clip_log10
    if (d < -clip_log10) d = -clip_log10
    ell_hat = d
  end subroutine linear_da_eval

  logical function linear_weights_ready()
    linear_weights_ready = linear_model_loaded
  end function linear_weights_ready

  subroutine linear_da_report()
    write(*,'(A,A,1X,A,6I12,1X,I18)') 'LINEAR_DA_COUNTS ', trim(linear_da), trim(score_kind), &
         n_invoked, n_eligible, n_score, n_stage1_rej, n_stage2, n_accept, aux_seed_used
  end subroutine linear_da_report

  real(dp) function prefix_logp(diagram, ids, nfound, as_y, iv1, iv2)
    type(fynman), intent(in) :: diagram
    integer, intent(in) :: ids(4), nfound, iv1, iv2
    logical, intent(in) :: as_y
    integer :: i, a, p, kout(3)
    real(dp) :: tau(4), dtau, w, t0, t1, ekv(dmc_band)
    prefix_logp = 0.0_dp
    if (nfound < 4) then
       prefix_logp = -1.0e30_dp
       return
    endif
    do i = 1, 4
       tau(i) = vertex_tau(diagram, ids(i), as_y, iv1, iv2)
    end do
    do i = 1, 3
       call vertex_kout(diagram, ids(i), as_y, iv1, iv2, kout)
       call cal_ek_int(kout, ekv)
       dtau = tau(i + 1) - tau(i)
       prefix_logp = prefix_logp - minval(ekv) * dtau
    end do
    do i = 1, 4
       a = ids(i)
       p = diagram%vertexList(a)%link(2)
       if (p < a) cycle
       w = diagram%vertexList(a)%wq(diagram%vertexList(a)%nu)
       t0 = vertex_tau(diagram, a, as_y, iv1, iv2)
       t1 = vertex_tau(diagram, p, as_y, iv1, iv2)
       prefix_logp = prefix_logp + linear_ld(w, t1 - t0)
    end do
  end function prefix_logp

  real(dp) function vertex_tau(diagram, iv, as_y, iv1, iv2)
    type(fynman), intent(in) :: diagram
    integer, intent(in) :: iv, iv1, iv2
    logical, intent(in) :: as_y
    if (as_y .and. iv == iv1) then
       vertex_tau = diagram%vertexList(iv2)%tau
    else if (as_y .and. iv == iv2) then
       vertex_tau = diagram%vertexList(iv1)%tau
    else
       vertex_tau = diagram%vertexList(iv)%tau
    endif
  end function vertex_tau

  subroutine vertex_kout(diagram, iv, as_y, iv1, iv2, kout)
    type(fynman), intent(in) :: diagram
    integer, intent(in) :: iv, iv1, iv2
    logical, intent(in) :: as_y
    integer, intent(out) :: kout(3)
    integer :: kin(3)
    if (.not. as_y) then
       kout = diagram%vertexList(iv)%i_kout
       return
    endif
    if (iv == iv2) then
       kin = diagram%vertexList(iv1)%i_kin
       kout = kin + diagram%vertexList(iv2)%i_q
    else if (iv == iv1) then
       kout = diagram%vertexList(iv2)%i_kout
    else
       kout = diagram%vertexList(iv)%i_kout
    endif
  end subroutine vertex_kout


  subroutine p1_prop_ell(diagram, iv1, iv2, ell)
    type(fynman), intent(in) :: diagram
    integer, intent(in) :: iv1, iv2
    real(dp), intent(out) :: ell
    integer :: kout(3), st
    real(dp) :: ekv(dmc_band), tau12, w1, w2, tL, tR, tPL, tPR
    character(len=8) :: env
    kout = diagram%vertexList(iv1)%i_kin + diagram%vertexList(iv2)%i_q
    call cal_ek_int(kout, ekv)
    tL = diagram%vertexList(iv1)%tau; tR = diagram%vertexList(iv2)%tau
    tau12 = tR - tL
    tPL = diagram%vertexList(diagram%vertexList(iv1)%link(2))%tau
    tPR = diagram%vertexList(diagram%vertexList(iv2)%link(2))%tau
    w1 = diagram%vertexList(iv1)%wq(diagram%vertexList(iv1)%nu)
    w2 = diagram%vertexList(iv2)%wq(diagram%vertexList(iv2)%nu)
    ell = -(minval(ekv) - minval(diagram%vertexList(iv1)%ekout)) * tau12
    ell = ell + linear_ld(w1, tR - tPL) - linear_ld(w1, tL - tPL)
    ell = ell + linear_ld(w2, tL - tPR) - linear_ld(w2, tR - tPR)
    if (ell > clip_log10) ell = clip_log10
    if (ell < -clip_log10) ell = -clip_log10
    ! C2_FORCE_A1 is fixture-only. Production (P1_FIXTURE off) ignores it
    ! so a stale env var cannot break occupancy reverse / detailed balance.
    if (p1_fixture_on()) then
      call get_environment_variable('C2_FORCE_A1', env, status=st)
      if (st==0) then
         if (trim(adjustl(env))=='1') ell = -clip_log10
      endif
    endif
  end subroutine

  integer*8 function c2_live_fp(diagram)
    type(fynman), intent(in) :: diagram
    integer :: i
    integer*8 :: s
    s = int(diagram%order,8) + 10007_8*int(diagram%nph_ext,8)
    do i = 1, diagram%order
      s = s + int(diagram%vertexList(i)%i_q(1)+20*diagram%vertexList(i)%i_q(2)+400*diagram%vertexList(i)%i_q(3),8)
      s = s + int(diagram%vertexList(i)%nu,8)
      s = s + nint(diagram%vertexList(i)%tau * 1.0e8_dp)
      s = s + int(diagram%vertexList(i)%i_kin(1)+20*diagram%vertexList(i)%i_kout(1),8)
    enddo
    i = maxN
    s = s + nint(diagram%vertexList(i)%tau * 1.0e8_dp)
    c2_live_fp = s
  end function

  subroutine c2_dump_swap(stage, iv1, iv2, diagram, ell_hat, ellR, ng, ne, ekold, eknew, re_old, re_new)
    use closure_observer
    integer, intent(in) :: stage, iv1, iv2
    type(fynman), intent(in) :: diagram
    real(dp), intent(in) :: ell_hat, ellR, ekold, eknew, re_old, re_new
    integer*8, intent(in) :: ng, ne
    integer, save :: u
    integer :: ios, st
    integer*8 :: fp
    logical, save :: opened = .false., dump_cfg = .false., dump_on = .true.
    character(len=8) :: env
    if (.not. dump_cfg) then
      dump_cfg = .true.
      call get_environment_variable('C2_DUMP', env, status=st)
      if (st==0) then
        if (trim(adjustl(env))=='0' .or. trim(adjustl(env))=='off') dump_on = .false.
      endif
    endif
    if (.not. dump_on) return
    if (.not. opened) then
      open(newunit=u, file='c2_swap.tsv', status='replace', action='write', recl=4096, iostat=ios)
      if (ios /= 0) error stop 'C2 swap dump open failed'
      opened = .true.
      write(u,'(a)') 'stage attempt iv1 iv2 kinx kiny kinz qrx qry qrz nu1 nu2 tauL tauR tpl tpr w1 w2 ell_hat ellR ekold eknew re_old re_new ng_delta nenv_delta order nph_ext fp'
    endif
    ! stage-1 count is in update_swap, independent of dump
    fp = c2_live_fp(diagram)
    write(u,'(12(i0,1x),12(es26.17e3,1x),4(i0,1x),i0)') &
      stage, c1_event, iv1, iv2, diagram%vertexList(iv1)%i_kin, diagram%vertexList(iv2)%i_q, &
      diagram%vertexList(iv1)%nu, diagram%vertexList(iv2)%nu, &
      diagram%vertexList(iv1)%tau, diagram%vertexList(iv2)%tau, &
      diagram%vertexList(diagram%vertexList(iv1)%link(2))%tau, &
      diagram%vertexList(diagram%vertexList(iv2)%link(2))%tau, &
      diagram%vertexList(iv1)%wq(diagram%vertexList(iv1)%nu), &
      diagram%vertexList(iv2)%wq(diagram%vertexList(iv2)%nu), ell_hat, ellR, ekold, eknew, re_old, re_new, ng, ne, &
      diagram%order, diagram%nph_ext, fp
  end subroutine
end module linear_da_mod
