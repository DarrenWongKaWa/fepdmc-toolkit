! Future B research: grouped-measure chain B with phonon modes summed over a
! non-crossing (laminar) set of lines, run on top of the native FEP-DMC kernel.
!
! Opt-in with FUTUREB_BCHAIN=1; otherwise every entry point returns at once.
! DMC_Method 0, any dmc_band, electron or hole. OMP_NUM_THREADS=1.
!
! Group. S(C) = greedy non-crossing set of internal lines by opening order
! (topology only; FUTUREB_BCHAIN_SRULE=max: maximum non-crossing set).
! S_eff(C) = lines of S whose current mode is live
! (omega >= phfreq cutoff). G(C) = all live-mode assignments on S_eff, other
! lines fixed. S_eff is invariant on G(C), so the groups partition the space.
!
! Target. w_B(C) = |Re F_G| / |G| with F_G = sum_{C' in G} D(C').
! F_G is exact: lines of S_eff are folded inside out, each fold replacing the
! vertices a..b by  S = sum_nu Dph_nu g_b(nu) M g_a(nu)  where M already holds any
! folded inner lines. Modes do not change momenta, so no energy is recomputed.
!
! Kernel. The native update K_A (reversible for pi_A ~ |Re D|) is the proposal:
! when it changes the state, a second stage accepts with min(1, r'/r),
! r = w_B/w_A; otherwise the snapshot is restored. After every step the modes
! of S_eff are redrawn uniformly among live modes every FUTUREB_BCHAIN_REFRESH
! steps (default 10; Gibbs, since w_B is constant on G; r is then recomputed).
! FUTUREB_BCHAIN_GROUP=0 disables grouping (S empty, r = 1): the chain is then
! bit-identical to the native one (control).
!
! Estimator. Blocks carry (V, H, W): value, energy insertion (G -> HG on one
! segment), phonon-energy insertion (omega_nu dtau on one grouped line). Then
!   f_O = Re[tr F_H + tr F_W + (W_fixed - (order-1)) tr F_V] / |Re tr F_V|,
!   f_1 = sign(Re tr F_V),   E = sum f_O / (tau_max sum f_1),
! the grouped version of the native measure_EZ_wfn increment.
!
! Self-checks (bchain_summary.dat): the current-mode fold reproduces the native
! Etrue increment and sign; for |S_eff| <= 2 the folded sums equal explicit
! enumeration of all members.
!
! Change-q move (FUTUREB_CHQ=1, probability FUTUREB_CHQ_P per step, default 0.05).
! The native kernel refreshes a line's momentum only by removing it, and only
! span-1 lines are removable. change_q is an exact MH move for pi_A: pick an
! internal line uniformly (as remove_ph), draw q' ~ Pq for its earlier vertex and
! nu' ~ Pnu(q') (as add_ph), shift every electron segment across the line by
! q' - q, recompute energies, eigenvectors and vertex matrices on the span, and
! accept with min(1, w'/w * Pq(q) Pnu(nu) / (Pq(q') Pnu'(nu'))), w = |Re D| from
! eval_state. Without FUTUREB_BCHAIN it runs on the native chain and writes
! chq_trace.dat (native increment, sign, order); with it, change_q is one more
! pi_A-reversible piece of the composite proposal and the second stage corrects it.
! FUTUREB_CHQ_EXT=1 also moves external (boundary-wrapping) pairs: the momentum of
! every segment outside the pair shifts (change_q_ext).
! FUTUREB_ANY=1 (probability FUTUREB_ANY_P per step) adds a general-span add/remove
! pair for internal lines (any_move): remove any internal line, or add one between
! tau1 ~ U(0, tau_max) and tau2 ~ truncated Exp(omega_nu) on (tau1, tau_max), with
! q ~ Pq and nu ~ Pnu; the span's momenta shift by -q / +q. Native add/remove only
! touch span-1 lines, which is what limits order mixing.
! FUTUREB_EXTAR=1 (probability FUTUREB_EXTAR_P) adds a general external add/remove
! (ext_move): a boundary-wrapping pair inserted or removed at any position, not only
! outermost as in native (whose LIFO rule makes the external-pair count the slowest
! variable, tau_int ~ 1800 measurements at LiF-hole 500 K).
! FUTUREB_MV_UNTIL=N limits the Future B moves (change-q, any_move) to the first
! N steps (e.g. the burn-in); afterwards the chain is purely native.
! Self-checks (chq_summary.dat): momentum conservation and continuity along the
! chain, and the native measure equals evaluate() on the current state.
module bchain_mod
   use DiagMC
   use pert_param, only : nsvd, bc_method => DMC_Method, bc_cut => phfreq_cutoff
   use pert_const, only : bc_ryd2ev => ryd2ev
   implicit none
   private
   public :: bchain_pre_update, bchain_post_update, bchain_pre_measure, bchain_post_measure

   logical, save :: bc_checked = .false., bc_on = .false., bc_group = .true., bc_init = .false.
   integer, save :: bc_unit = 0, nbb = 0, snap_order = 0, snap_nph_ext = 0, refresh_every = 10
   integer(8), save :: acc0 = 0, n_prop = 0, n_acc2 = 0, n_meas = 0, n_bad_o = 0, n_bad_sign = 0
   integer(8), save :: accv0(7) = 0, prop_t(7) = 0, acc2_t(7) = 0
   logical, save :: bc_addrem = .false.
   integer, save :: srule = 0          ! 0 greedy by opening order, 1 maximum non-crossing set
   integer, save :: blk = 1            ! native steps per second-stage correction (K_A^blk)
   logical, save :: blk_changed = .false.
   integer(8), save :: n_blocks = 0, n_blocks_changed = 0, n_blocks_acc = 0
   integer, allocatable, save :: fdp(:,:)
   real(dp), save :: p_addrem = 0.2d0
   integer(8), save :: n_addB = 0, n_addB_acc = 0, n_remB = 0, n_remB_acc = 0
   integer(8), save :: n_bad_enum = 0, n_enum = 0, n_steps = 0
   real(dp), save :: r_snap = 1.d0, e_before = 0.d0, g_before = 0.d0, max_o_rel = 0.d0, max_enum_rel = 0.d0
   real(dp), save :: sum_seff = 0.d0, sum_loggroup = 0.d0
   integer(8), save :: rng_state = 88172645463325252_8
   type(vertex), allocatable, save :: snap(:)
   complex(dp), allocatable, save :: UV(:,:,:), UH(:,:,:), UW(:,:,:)
   real(dp), allocatable, save :: GdB(:,:), HGdB(:,:)
   real(dp), parameter :: tol = 1.d-8
   ! wall-clock profile (seconds)
   real(dp), save :: t_native = 0, t_ratio = 0, t_snap = 0, t_bmove = 0, t_refresh = 0, t_meas = 0, t_natmeas = 0
   real(dp), save :: t_mark = 0, t_meas_mark = 0
   ! cached absolute logs of the current state: W_B (group sum), w_A (|Re D|), |G|
   real(dp), save :: lwb_cur = 0, lwa_cur = 0, lgrp_cur = 0
   integer(8), save :: n_cache_checks = 0, n_cache_bad = 0
   real(dp), save :: max_cache_err = 0
   ! change-q move
   logical, save :: chq_on = .false., chq_ext_on = .false.
   integer(8), save :: n_chq_ext_try = 0, n_chq_ext_acc = 0, n_rt = 0
   logical, save :: chq_rt = .false.
   real(dp), save :: rt_la = 0, rt_pnu = 0, rt_g = 0, rt_u = 0
   type(vertex), allocatable, save :: csnap2(:)
   ! general-span add/remove
   logical, save :: any_on = .false., mv_on = .false.
   real(dp), save :: p_any = 0.05d0
   integer(8), save :: n_anyA = 0, n_anyA_acc = 0, n_anyR = 0, n_anyR_acc = 0, n_gchecks = 0, n_gbad = 0
   real(dp), save :: sum_spanA = 0, sum_spanR = 0, t_any = 0
   integer, save :: fs_order = 0, fs_nph_ext = 0
   integer(8), save :: mv_until = huge(1_8), n_mvsteps = 0
   ! general external add/remove
   logical, save :: ext_on = .false., ext_outer = .false.   ! FUTUREB_EXTAR_OUTER=1: native support only
   logical, save :: ext_mimic = .false.   ! FUTUREB_EXTAR_MIMIC=1: native proposal densities, eval_state weight
   logical, save :: ext_rpos = .false., ext_rtau = .false.   ! outermost-only / tau1 < tmax/2 < tau2 only
   ! native's diagram space: every head-attached external vertex precedes every
   ! tail-attached one (update_swap forbids crossing them; add is outermost). The
   ! general add must respect it (FUTUREB_EXTAR_UNORDERED=1 drops the constraint).
   logical, save :: ext_ordered = .true.
   real(dp), save :: p_ext = 0.05d0, t_ext = 0
   integer(8), save :: n_extA = 0, n_extA_acc = 0, n_extR = 0, n_extR_acc = 0, n_rt_ext = 0
   real(dp), save :: rt_ext_la = 0, rt_ext_pa = 0, rt_ext_pr = 0
   logical, save :: pnu_fro = .false.          ! FUTUREB_PNU_FRO=1: gauge-invariant Pnu ~ sum |g_ij|^2
   real(dp), save :: gauge_dpnu = 0            ! max |Pnu(stored) - Pnu(fresh)| seen in any_remove
   type(vertex), save :: vtmp
   logical, save :: vtmp_init = .false.
   real(dp), save :: p_chq = 0.05d0
   integer(8), save :: n_chq = 0, n_chq_acc = 0, n_chq_ext = 0, n_kchecks = 0, n_kbad = 0
   real(dp), save :: sum_span_try = 0, sum_span_acc = 0, t_chq = 0
   type(vertex), allocatable, save :: csnap(:)

contains

   real(dp) function now()
      integer(8) :: c, r
      call system_clock(c, r)
      now = real(c, dp) / real(r, dp)
   end function now

   subroutine bc_setup(diagram)
      type(fynman), intent(inout) :: diagram
      character(len=32) :: env
      integer :: st, iv, sd
      bc_checked = .true.
      call get_environment_variable('FUTUREB_BCHAIN', env, status=st)
      bc_on = (st == 0 .and. trim(adjustl(env)) == '1')
      call get_environment_variable('FUTUREB_CHQ', env, status=st)
      chq_on = (st == 0 .and. trim(adjustl(env)) == '1')
      call get_environment_variable('FUTUREB_ANY', env, status=st)
      any_on = (st == 0 .and. trim(adjustl(env)) == '1')
      call get_environment_variable('FUTUREB_EXTAR', env, status=st)
      ext_on = (st == 0 .and. trim(adjustl(env)) == '1')
      mv_on = chq_on .or. any_on .or. ext_on
      if (.not. (bc_on .or. mv_on)) return
      if (bc_method /= 0) stop 'FUTUREB_BCHAIN/CHQ: DMC_Method 0 only'
      call get_environment_variable('FUTUREB_CHQ_P', env, status=st)
      if (st == 0 .and. len_trim(env) > 0) read(env, *) p_chq
      call get_environment_variable('FUTUREB_HERM_TEST', env, status=st)
      if (st == 0 .and. trim(adjustl(env)) == '1') call herm_test(diagram)
      call get_environment_variable('FUTUREB_EXTAR_OUTER', env, status=st)
      ext_outer = (st == 0 .and. trim(adjustl(env)) == '1')
      call get_environment_variable('FUTUREB_EXTAR_MIMIC', env, status=st)
      ext_mimic = (st == 0 .and. trim(adjustl(env)) == '1')
      if (ext_mimic) ext_outer = .true.
      call get_environment_variable('FUTUREB_EXTAR_UNORDERED', env, status=st)
      if (st == 0 .and. trim(adjustl(env)) == '1') ext_ordered = .false.
      call get_environment_variable('FUTUREB_EXTAR_RPOS', env, status=st)
      ext_rpos = ext_outer .or. (st == 0 .and. trim(adjustl(env)) == '1')
      call get_environment_variable('FUTUREB_EXTAR_RTAU', env, status=st)
      ext_rtau = ext_outer .or. (st == 0 .and. trim(adjustl(env)) == '1')
      call get_environment_variable('FUTUREB_EXTAR_P', env, status=st)
      if (st == 0 .and. len_trim(env) > 0) read(env, *) p_ext
      call get_environment_variable('FUTUREB_MV_UNTIL', env, status=st)
      if (st == 0 .and. len_trim(env) > 0) read(env, *) mv_until
      call get_environment_variable('FUTUREB_PNU_FRO', env, status=st)
      pnu_fro = (st == 0 .and. trim(adjustl(env)) == '1')
      call get_environment_variable('FUTUREB_ANY_P', env, status=st)
      if (st == 0 .and. len_trim(env) > 0) read(env, *) p_any
      call get_environment_variable('FUTUREB_CHQ_EXT', env, status=st)
      chq_ext_on = (st == 0 .and. trim(adjustl(env)) == '1')
      call get_environment_variable('FUTUREB_CHQ_RT', env, status=st)
      chq_rt = (st == 0 .and. trim(adjustl(env)) == '1')
      call get_environment_variable('FUTUREB_BCHAIN_GROUP', env, status=st)
      if (st == 0 .and. trim(adjustl(env)) == '0') bc_group = .false.
      if (.not. bc_on) bc_group = .false.
      call get_environment_variable('FUTUREB_BCHAIN_BLOCK', env, status=st)
      if (st == 0 .and. len_trim(env) > 0) read(env, *) blk
      if (mod(100, blk) /= 0) stop 'FUTUREB_BCHAIN_BLOCK must divide 100 (measurement stride)'
      call get_environment_variable('FUTUREB_BCHAIN_SRULE', env, status=st)
      if (st == 0 .and. trim(adjustl(env)) == 'max') srule = 1
      if (srule == 1) allocate(fdp(0:maxN + 1, 0:maxN + 1))
      call get_environment_variable('FUTUREB_BCHAIN_ADDREM', env, status=st)
      bc_addrem = bc_group .and. (st == 0 .and. trim(adjustl(env)) == '1')
      call get_environment_variable('FUTUREB_BCHAIN_PAR', env, status=st)
      if (st == 0 .and. len_trim(env) > 0) read(env, *) p_addrem
      call get_environment_variable('FUTUREB_BCHAIN_REFRESH', env, status=st)
      if (st == 0 .and. len_trim(env) > 0) read(env, *) refresh_every
      call get_environment_variable('FUTUREB_SEED', env, status=st)
      if (st == 0 .and. len_trim(env) > 0) then
         read(env, *) sd
         rng_state = 88172645463325252_8 + int(sd, 8) * 2654435761_8
      end if
      nbb = dmc_band
      allocate(snap(maxN))
      do iv = 1, maxN
         call CreateVertex(Nph, dmc_band, nbnd, nsvd, snap(iv))
      end do
      allocate(UV(nbb,nbb,maxN), UH(nbb,nbb,maxN), UW(nbb,nbb,maxN), GdB(nbb,maxN), HGdB(nbb,maxN))
      if (mv_on) then
         allocate(csnap(maxN), csnap2(maxN))
         do iv = 1, maxN
            call CreateVertex(Nph, dmc_band, nbnd, nsvd, csnap(iv))
            call CreateVertex(Nph, dmc_band, nbnd, nsvd, csnap2(iv))
         end do
      end if
      if (bc_on) then
         open(newunit=bc_unit, file='bchain_trace.dat', status='replace', action='write')
         write(bc_unit, '(A)') '# f_O f_1 n_Seff log_group native_inc native_sign order'
      else
         open(newunit=bc_unit, file='chq_trace.dat', status='replace', action='write')
         write(bc_unit, '(A)') '# native_inc native_sign order nph_ext khead2 n_int_lines'
      end if
      if (diagram%order < 0) stop 'bchain: bad diagram'
   end subroutine bc_setup

   real(dp) function urand()
      integer(8) :: x
      x = rng_state
      x = ieor(x, ishft(x, 13))
      x = ieor(x, ishft(x, -7))
      x = ieor(x, ishft(x, 17))
      rng_state = x
      urand = 0.5d0 + 0.5d0 * real(x, dp) / real(huge(x), dp)
      if (urand <= 0.d0) urand = 1.d-16
      if (urand >= 1.d0) urand = 1.d0 - 1.d-16
   end function urand

   complex(dp) function ctrb(A)
      complex(dp), intent(in) :: A(:,:)
      integer :: i
      ctrb = (0.d0, 0.d0)
      do i = 1, size(A, 1)
         ctrb = ctrb + A(i, i)
      end do
   end function ctrb

   subroutine eyeb(A)
      complex(dp), intent(out) :: A(:,:)
      integer :: i
      A = (0.d0, 0.d0)
      do i = 1, size(A, 1)
         A(i, i) = (1.d0, 0.d0)
      end do
   end subroutine eyeb

   logical function live(w)
      real(dp), intent(in) :: w
      live = .not. (w < bc_cut * bc_ryd2ev)
   end function live


   ! The line set S (topology only, never modes). srule 0: greedy by opening order.
   ! srule 1: a maximum-cardinality non-crossing subset, by the O(n^2) chord DP
   !   f(i,j) = max( f(i+1,j), 1 + f(i+1,p-1) + f(p+1,j) )  if i < p = partner(i) <= j,
   ! ties resolved towards including the chord (deterministic, mode-independent).
   subroutine select_S(nl, la, lb, inS, ns)
      integer, intent(in) :: nl, la(maxN), lb(maxN)
      integer, intent(out) :: inS(maxN), ns
      integer :: l, l2, j, npos, pl(0:maxN + 1), lineat(0:maxN + 1), i, p, top, stk_i(maxN), stk_j(maxN), inc
      logical :: ok
      ns = 0
      if (srule == 0) then
         do l = 1, nl
            ok = .true.
            do l2 = 1, ns
               j = inS(l2)
               if ((la(j) < la(l) .and. la(l) < lb(j) .and. lb(j) < lb(l)) .or. &
                   (la(l) < la(j) .and. la(j) < lb(l) .and. lb(l) < lb(j))) then
                  ok = .false.
                  exit
               end if
            end do
            if (ok) then
               ns = ns + 1
               inS(ns) = l
            end if
         end do
         return
      end if
      npos = 0
      pl = 0
      lineat = 0
      do l = 1, nl
         pl(la(l)) = lb(l)
         lineat(la(l)) = l
         npos = max(npos, lb(l))
      end do
      do i = npos + 1, 1, -1
         fdp(i, i - 1) = 0
      end do
      do i = npos, 1, -1
         do j = i, npos
            fdp(i, j) = fdp(i + 1, j)
            p = pl(i)
            if (p > i .and. p <= j) then
               inc = 1 + fdp(i + 1, p - 1)
               if (p + 1 <= j) inc = inc + fdp(p + 1, j)
               if (inc >= fdp(i, j)) fdp(i, j) = inc
            end if
         end do
      end do
      top = 1
      stk_i(1) = 1
      stk_j(1) = npos
      do while (top > 0)
         i = stk_i(top)
         j = stk_j(top)
         top = top - 1
         if (i > j) cycle
         p = pl(i)
         inc = -1
         if (p > i .and. p <= j) then
            inc = 1 + fdp(i + 1, p - 1)
            if (p + 1 <= j) inc = inc + fdp(p + 1, j)
         end if
         if (inc >= 0 .and. inc >= fdp(i + 1, j)) then
            ns = ns + 1
            inS(ns) = lineat(i)
            top = top + 1; stk_i(top) = i + 1; stk_j(top) = p - 1
            top = top + 1; stk_i(top) = p + 1; stk_j(top) = j
         else
            top = top + 1; stk_i(top) = i + 1; stk_j(top) = j
         end if
      end do
   end subroutine select_S

   ! ---- snapshot of the active diagram (slots 1..order and maxN) --------------
   subroutine take_snapshot(diagram)
      type(fynman), intent(inout) :: diagram
      integer :: iv
      do iv = 1, diagram%order
         call copy_vtex(Nph, dmc_band, nbnd, nsvd, diagram%vertexList(iv), snap(iv))
      end do
      call copy_vtex(Nph, dmc_band, nbnd, nsvd, diagram%vertexList(maxN), snap(maxN))
      snap_order = diagram%order
      snap_nph_ext = diagram%nph_ext
   end subroutine take_snapshot

   subroutine restore_snapshot(diagram)
      type(fynman), intent(inout) :: diagram
      integer :: iv
      do iv = 1, snap_order
         call copy_vtex(Nph, dmc_band, nbnd, nsvd, snap(iv), diagram%vertexList(iv))
      end do
      call copy_vtex(Nph, dmc_band, nbnd, nsvd, snap(maxN), diagram%vertexList(maxN))
      diagram%order = snap_order
      diagram%nph_ext = snap_nph_ext
   end subroutine restore_snapshot

   ! ---- the grouped evaluation --------------------------------------------------
   ! mode_all = .true.: live-mode sums on S_eff (F_G); .false.: current modes only (D).
   ! Returns tr of V, H, W, the fixed phonon-energy term, |S_eff|, log|G|.
   subroutine evaluate(diagram, mode_all, trV, trH, trW, wfixed, nseff, loggroup, forced_nl, forced, logpc)
      type(fynman), intent(inout), target :: diagram
      logical, intent(in) :: mode_all
      complex(dp), intent(out) :: trV, trH, trW
      real(dp), intent(out) :: wfixed, loggroup
      integer, intent(out) :: nseff
      integer, intent(in), optional :: forced_nl, forced(:)
      real(dp), intent(out), optional :: logpc
      real(dp) :: lpc
      integer :: ids(maxN), posof(maxN), la(maxN), lb(maxN), inS(maxN), nxt(maxN), endp(maxN)
      integer :: nv, iv, j, l, l2, nl, ns, partner, s, u, a, b, nu, nlive, ord(maxN), tmp, i
      real(dp) :: emin, el, dt, cnorm(maxN), dtl, d
      logical :: ok, grouped(maxN)
      complex(dp) :: MV(nbb,nbb), MH(nbb,nbb), MW(nbb,nbb), SV(nbb,nbb), SH(nbb,nbb), SW(nbb,nbb)
      complex(dp) :: TV(nbb,nbb), TH(nbb,nbb), TW(nbb,nbb), ga(nbb,nbb), gb(nbb,nbb), X1(nbb,nbb)
      type(vertex), pointer :: vi, vp, va, vb

      nv = 0
      iv = diagram%vertexList(1)%link(3)
      do while (iv /= maxN)
         nv = nv + 1
         ids(nv) = iv
         iv = diagram%vertexList(iv)%link(3)
      end do
      posof = 0
      do j = 1, nv
         posof(ids(j)) = j
      end do
      ! segments 1..nv+1 (native propagator / Hpropagator)
      lpc = 0.d0
      do s = 1, nv + 1
         if (s == 1) then
            vi => diagram%vertexList(1)
         else
            vi => diagram%vertexList(ids(s - 1))
         end if
         if (s <= nv) then
            vp => diagram%vertexList(ids(s))
         else
            vp => diagram%vertexList(maxN)
         end if
         dt = vp%tau - vi%tau
         emin = minval(vi%ekout(1:nbb))
         lpc = lpc - emin * dt
         do i = 1, nbb
            el = vi%ekout(i) - emin
            GdB(i, s) = exp(-min(dt * el, 100.0_dp))
            HGdB(i, s) = exp(-dt * el) * (el + emin) * dt
         end do
      end do
      ! lines and the laminar set
      nl = 0
      wfixed = 0.d0
      do j = 1, nv
         vi => diagram%vertexList(ids(j))
         partner = vi%link(2)
         if (partner == 1 .or. partner == maxN) then
            wfixed = wfixed + vi%wq(vi%nu) * abs(vi%tau - diagram%vertexList(partner)%tau)
            lpc = lpc + log(Dph(vi%wq(vi%nu), abs(vi%tau - diagram%vertexList(partner)%tau)))
            cycle
         end if
         if (posof(partner) > j) then
            nl = nl + 1
            la(nl) = j
            lb(nl) = posof(partner)
         end if
      end do
      ns = 0
      grouped = .false.
      if (present(forced)) then
         do l = 1, forced_nl
            grouped(forced(l)) = .true.
         end do
      else if (bc_group) then
         call select_S(nl, la, lb, inS, ns)
         do l2 = 1, ns
            l = inS(l2)
            va => diagram%vertexList(ids(la(l)))
            if (live(va%wq(va%nu))) grouped(l) = .true.
         end do
      end if
      nseff = 0
      loggroup = 0.d0
      do l = 1, nl
         va => diagram%vertexList(ids(la(l)))
         dtl = diagram%vertexList(ids(lb(l)))%tau - va%tau
         if (grouped(l)) then
            nseff = nseff + 1
            nlive = 0
            do nu = 1, Nph
               if (live(va%wq(nu))) nlive = nlive + 1
            end do
            loggroup = loggroup + log(real(nlive, dp))
         else
            wfixed = wfixed + va%wq(va%nu) * dtl
            lpc = lpc + log(Dph(va%wq(va%nu), dtl))
         end if
      end do
      ! units: vertex matrices normalized by the current-mode max (fixed per state)
      do j = 1, nv
         vi => diagram%vertexList(ids(j))
         cnorm(j) = maxval(abs(vi%gkq(:, :, vi%nu)))
         lpc = lpc + log(cnorm(j))
         UV(:, :, j) = vi%gkq(:, :, vi%nu) / cnorm(j)
         UH(:, :, j) = (0.d0, 0.d0)
         UW(:, :, j) = (0.d0, 0.d0)
         nxt(j) = j + 1
         endp(j) = j
      end do
      ! fold grouped lines, innermost (shortest) first
      ns = 0
      do l = 1, nl
         if (grouped(l)) then
            ns = ns + 1
            ord(ns) = l
         end if
      end do
      do i = 2, ns
         do j = i, 2, -1
            if (lb(ord(j - 1)) - la(ord(j - 1)) > lb(ord(j)) - la(ord(j))) then
               tmp = ord(j); ord(j) = ord(j - 1); ord(j - 1) = tmp
            end if
         end do
      end do
      do i = 1, ns
         l = ord(i)
         a = la(l)
         b = lb(l)
         va => diagram%vertexList(ids(a))
         vb => diagram%vertexList(ids(b))
         dtl = vb%tau - va%tau
         ! interior product M (segments and units strictly between a and b)
         call eyeb(MV)
         MH = (0.d0, 0.d0)
         MW = (0.d0, 0.d0)
         s = a + 1
         call seg_mult(s, MV, MH, MW)
         u = nxt(a)
         do while (u /= b)
            call unit_mult(u, MV, MH, MW)
            s = endp(u) + 1
            call seg_mult(s, MV, MH, MW)
            u = nxt(u)
         end do
         SV = (0.d0, 0.d0)
         SH = (0.d0, 0.d0)
         SW = (0.d0, 0.d0)
         do nu = 1, Nph
            if (mode_all) then
               if (.not. live(va%wq(nu))) cycle
            else
               if (nu /= va%nu) cycle
            end if
            d = Dph(va%wq(nu), dtl)
            ga = va%gkq(:, :, nu) / cnorm(a)
            gb = vb%gkq(:, :, nu) / cnorm(b)
            X1 = matmul(gb, matmul(MV, ga))
            SV = SV + d * X1
            SH = SH + d * matmul(gb, matmul(MH, ga))
            SW = SW + d * (matmul(gb, matmul(MW, ga)) + va%wq(nu) * dtl * X1)
         end do
         UV(:, :, a) = SV
         UH(:, :, a) = SH
         UW(:, :, a) = SW
         nxt(a) = nxt(b)
         endp(a) = b
      end do
      ! full chain
      call eyeb(TV)
      TH = (0.d0, 0.d0)
      TW = (0.d0, 0.d0)
      call seg_mult(1, TV, TH, TW)
      if (nv > 0) then
         u = 1
         do while (u <= nv)
            call unit_mult(u, TV, TH, TW)
            s = endp(u) + 1
            call seg_mult(s, TV, TH, TW)
            u = nxt(u)
         end do
      end if
      trV = ctrb(TV)
      trH = ctrb(TH)
      trW = ctrb(TW)
      if (present(logpc)) logpc = lpc
      ! ungrouped lines carry their fixed Dph in the value only through native
      ! weights; for ratios every configuration of the group shares them, but the
      ! current-mode evaluation must include the grouped lines' Dph (done in folds).
   contains
      subroutine seg_mult(sidx, AV, AH, AW)
         integer, intent(in) :: sidx
         complex(dp), intent(inout) :: AV(nbb,nbb), AH(nbb,nbb), AW(nbb,nbb)
         integer :: q
         complex(dp) :: OV(nbb,nbb)
         OV = AV
         do q = 1, nbb
            AV(q, :) = GdB(q, sidx) * OV(q, :)
            AH(q, :) = GdB(q, sidx) * AH(q, :) + HGdB(q, sidx) * OV(q, :)
            AW(q, :) = GdB(q, sidx) * AW(q, :)
         end do
      end subroutine seg_mult
      subroutine unit_mult(uidx, AV, AH, AW)
         integer, intent(in) :: uidx
         complex(dp), intent(inout) :: AV(nbb,nbb), AH(nbb,nbb), AW(nbb,nbb)
         complex(dp) :: OV(nbb,nbb)
         OV = AV
         AV = matmul(UV(:, :, uidx), OV)
         AH = matmul(UV(:, :, uidx), AH) + matmul(UH(:, :, uidx), OV)
         AW = matmul(UV(:, :, uidx), AW) + matmul(UW(:, :, uidx), OV)
      end subroutine unit_mult
   end subroutine evaluate

   ! r(C) = w_B / w_A  (common positive factors cancel)
   real(dp) function ratio_r(diagram)
      type(fynman), intent(inout), target :: diagram
      complex(dp) :: fv, fh, fw, dv, dh, dw
      real(dp) :: wf, lg, wf2, lg2
      integer :: ns, ns2
      call evaluate(diagram, .true., fv, fh, fw, wf, ns, lg)
      call evaluate(diagram, .false., dv, dh, dw, wf2, ns2, lg2)
      if (abs(real(dv, dp)) <= 0.d0) then
         ratio_r = huge(1.d0)
      else
         ratio_r = abs(real(fv, dp)) / abs(real(dv, dp)) * exp(-lg)
      end if
   end function ratio_r

   ! Uniform redraw of the live modes of S_eff (group-invariant Gibbs move).
   subroutine refresh_modes(diagram)
      type(fynman), intent(inout), target :: diagram
      integer :: ids(maxN), posof(maxN), la(maxN), lb(maxN), inS(maxN), nv, iv, j, l, l2, nl, ns, partner
      integer :: livelist(Nph), nlive, nu, k
      logical :: ok
      type(vertex), pointer :: va, vb
      if (.not. bc_group) return
      nv = 0
      iv = diagram%vertexList(1)%link(3)
      do while (iv /= maxN)
         nv = nv + 1
         ids(nv) = iv
         iv = diagram%vertexList(iv)%link(3)
      end do
      posof = 0
      do j = 1, nv
         posof(ids(j)) = j
      end do
      nl = 0
      do j = 1, nv
         partner = diagram%vertexList(ids(j))%link(2)
         if (partner == 1 .or. partner == maxN) cycle
         if (posof(partner) > j) then
            nl = nl + 1
            la(nl) = j
            lb(nl) = posof(partner)
         end if
      end do
      ns = 0
      call select_S(nl, la, lb, inS, ns)
      do l2 = 1, ns
         l = inS(l2)
         va => diagram%vertexList(ids(la(l)))
         vb => diagram%vertexList(ids(lb(l)))
         if (.not. live(va%wq(va%nu))) cycle
         nlive = 0
         do nu = 1, Nph
            if (live(va%wq(nu))) then
               nlive = nlive + 1
               livelist(nlive) = nu
            end if
         end do
         k = 1 + int(urand() * nlive)
         if (k > nlive) k = nlive
         va%nu = livelist(k)
         vb%nu = livelist(k)
      end do
   end subroutine refresh_modes



   ! Fast state evaluation (value component only). Absolute logs, up to one global
   ! constant shared by every state:
   !   lwa  = log|Re D(C)|          (current modes; grouped lines' Dph included)
   !   lwb  = log|Re F_G(C)|        (live-mode sums on S_eff, exact fold)  [if do_fold]
   !   lgrp = log|G(C)|
   subroutine eval_state(diagram, do_fold, lwb, lwa, lgrp)
      type(fynman), intent(inout), target :: diagram
      logical, intent(in) :: do_fold
      real(dp), intent(out) :: lwb, lwa, lgrp
      integer :: ids(maxN), posof(maxN), la(maxN), lb(maxN), inS(maxN), nxt(maxN), endp(maxN), ord(maxN)
      integer :: nv, iv, j, l, l2, nl, ns, partner, s, u, a, b, nu, nlive, tmp, i, q
      real(dp) :: emin, el, dt, cn, dtl, d, lpc, lgd
      logical :: ok, grouped(maxN)
      complex(dp) :: TV(nbb,nbb), MV(nbb,nbb), SV(nbb,nbb), X1(nbb,nbb), OV(nbb,nbb)
      type(vertex), pointer :: vi, vp, va, vb
      nv = 0
      iv = diagram%vertexList(1)%link(3)
      do while (iv /= maxN)
         nv = nv + 1
         ids(nv) = iv
         iv = diagram%vertexList(iv)%link(3)
      end do
      posof = 0
      do j = 1, nv
         posof(ids(j)) = j
      end do
      lpc = 0.d0
      do s = 1, nv + 1
         if (s == 1) then
            vi => diagram%vertexList(1)
         else
            vi => diagram%vertexList(ids(s - 1))
         end if
         if (s <= nv) then
            vp => diagram%vertexList(ids(s))
         else
            vp => diagram%vertexList(maxN)
         end if
         dt = vp%tau - vi%tau
         emin = minval(vi%ekout(1:nbb))
         lpc = lpc - emin * dt
         do i = 1, nbb
            el = vi%ekout(i) - emin
            GdB(i, s) = exp(-min(dt * el, 100.0_dp))
         end do
      end do
      nl = 0
      do j = 1, nv
         vi => diagram%vertexList(ids(j))
         partner = vi%link(2)
         if (partner == 1 .or. partner == maxN) then
            lpc = lpc + log(Dph(vi%wq(vi%nu), abs(vi%tau - diagram%vertexList(partner)%tau)))
            cycle
         end if
         if (posof(partner) > j) then
            nl = nl + 1
            la(nl) = j
            lb(nl) = posof(partner)
         end if
      end do
      grouped = .false.
      ns = 0
      if (bc_group) then
         call select_S(nl, la, lb, inS, ns)
         do l2 = 1, ns
            l = inS(l2)
            va => diagram%vertexList(ids(la(l)))
            if (live(va%wq(va%nu))) grouped(l) = .true.
         end do
      end if
      lgrp = 0.d0
      lgd = 0.d0
      do l = 1, nl
         va => diagram%vertexList(ids(la(l)))
         dtl = diagram%vertexList(ids(lb(l)))%tau - va%tau
         if (grouped(l)) then
            nlive = 0
            do nu = 1, Nph
               if (live(va%wq(nu))) nlive = nlive + 1
            end do
            lgrp = lgrp + log(real(nlive, dp))
            lgd = lgd + log(Dph(va%wq(va%nu), dtl))
         else
            lpc = lpc + log(Dph(va%wq(va%nu), dtl))
         end if
      end do
      do j = 1, nv
         vi => diagram%vertexList(ids(j))
         cn = maxval(abs(vi%gkq(:, :, vi%nu)))
         lpc = lpc + log(cn)
         UV(:, :, j) = vi%gkq(:, :, vi%nu) / cn
         nxt(j) = j + 1
         endp(j) = j
      end do
      ! plain chain (current modes)
      call eyeb(TV)
      do q = 1, nbb
         TV(q, :) = GdB(q, 1) * TV(q, :)
      end do
      do j = 1, nv
         TV = matmul(UV(:, :, j), TV)
         do q = 1, nbb
            TV(q, :) = GdB(q, j + 1) * TV(q, :)
         end do
      end do
      d = abs(real(ctrb(TV), dp))
      if (d > 0.d0) then
         lwa = log(d) + lpc + lgd
      else
         lwa = -huge(1.d0)
      end if
      lwb = lwa
      if (.not. do_fold) return
      ! exact live-mode fold of the grouped lines, innermost first
      ns = 0
      do l = 1, nl
         if (grouped(l)) then
            ns = ns + 1
            ord(ns) = l
         end if
      end do
      do i = 2, ns
         do j = i, 2, -1
            if (lb(ord(j - 1)) - la(ord(j - 1)) > lb(ord(j)) - la(ord(j))) then
               tmp = ord(j); ord(j) = ord(j - 1); ord(j - 1) = tmp
            end if
         end do
      end do
      do i = 1, ns
         l = ord(i)
         a = la(l)
         b = lb(l)
         va => diagram%vertexList(ids(a))
         vb => diagram%vertexList(ids(b))
         dtl = vb%tau - va%tau
         call eyeb(MV)
         do q = 1, nbb
            MV(q, :) = GdB(q, a + 1) * MV(q, :)
         end do
         u = nxt(a)
         do while (u /= b)
            MV = matmul(UV(:, :, u), MV)
            s = endp(u) + 1
            do q = 1, nbb
               MV(q, :) = GdB(q, s) * MV(q, :)
            end do
            u = nxt(u)
         end do
         SV = (0.d0, 0.d0)
         cn = maxval(abs(va%gkq(:, :, va%nu)))
         d = maxval(abs(vb%gkq(:, :, vb%nu)))
         do nu = 1, Nph
            if (.not. live(va%wq(nu))) cycle
            X1 = matmul(vb%gkq(:, :, nu) / d, matmul(MV, va%gkq(:, :, nu) / cn))
            SV = SV + Dph(va%wq(nu), dtl) * X1
         end do
         UV(:, :, a) = SV
         nxt(a) = nxt(b)
         endp(a) = b
      end do
      call eyeb(TV)
      do q = 1, nbb
         TV(q, :) = GdB(q, 1) * TV(q, :)
      end do
      u = 1
      do while (u <= nv)
         OV = TV
         TV = matmul(UV(:, :, u), OV)
         s = endp(u) + 1
         do q = 1, nbb
            TV(q, :) = GdB(q, s) * TV(q, :)
         end do
         u = nxt(u)
      end do
      d = abs(real(ctrb(TV), dp))
      if (d > 0.d0) then
         lwb = log(d) + lpc
      else
         lwb = -huge(1.d0)
      end if
   end subroutine eval_state

   subroutine set_current(diagram)
      type(fynman), intent(inout), target :: diagram
      call eval_state(diagram, .true., lwb_cur, lwa_cur, lgrp_cur)
      r_snap = exp(max(min(lwb_cur - lgrp_cur - lwa_cur, 700.d0), -700.d0))
   end subroutine set_current

   ! Periodic consistency check of the cache against the slow full evaluation.
   subroutine check_cache(diagram)
      type(fynman), intent(inout), target :: diagram
      real(dp) :: lb1, la1, lg1, rr, err
      call eval_state(diagram, .true., lb1, la1, lg1)
      rr = ratio_r(diagram)
      err = max(abs(lb1 - lwb_cur), abs(la1 - lwa_cur), abs(lg1 - lgrp_cur), &
                abs((lb1 - lg1 - la1) - log(rr)))
      n_cache_checks = n_cache_checks + 1
      max_cache_err = max(max_cache_err, err)
      if (err > 1.d-8) n_cache_bad = n_cache_bad + 1
   end subroutine check_cache

   ! log of the absolute group weight W_B = |Re F_G| (sum over members), up to a
   ! global constant: log|Re tr F_norm| + log of the positive factors dropped by
   ! the normalized fold (segment exp(-Emin dtau), vertex max|g|, Dph of the
   ! ungrouped and external lines).
   real(dp) function log_wb(diagram)
      type(fynman), intent(inout), target :: diagram
      complex(dp) :: fv, fh, fw
      real(dp) :: wf, lg, lpc
      integer :: ns
      call evaluate(diagram, .true., fv, fh, fw, wf, ns, lg, logpc=lpc)
      if (abs(real(fv, dp)) <= 0.d0) then
         log_wb = -huge(1.d0)
      else
         log_wb = log(abs(real(fv, dp))) + lpc
      end if
   end function log_wb

   subroutine live_modes(wq, nlive, wbar)
      real(dp), intent(in) :: wq(Nph)
      integer, intent(out) :: nlive
      real(dp), intent(out) :: wbar
      integer :: nu
      nlive = 0
      wbar = 0.d0
      do nu = 1, Nph
         if (live(wq(nu))) then
            nlive = nlive + 1
            wbar = wbar + wq(nu)
         end if
      end do
      if (nlive > 0) wbar = wbar / nlive
   end subroutine live_modes

   ! Add a span-1 line on a random electron segment. As native add_ph, except that
   ! the mode is summed: tau2 is drawn with the mode-independent decay
   ! wbar + Emin_out - Emin_in, the representative mode is uniform over live modes,
   ! and acceptance uses W_B. The new line crosses nothing, so it joins S_eff and
   ! |G'| = |G| n_live, which cancels the 1/n_live of the representative choice.
   subroutine add_B(diagram)
      type(fynman), intent(inout), target :: diagram
      integer :: ivL, ivR, ivn1, ivn2, nlive, k, nu, order0
      real(dp) :: Pq, tauL, tauR, tau1, tau2, temp, Ptau12, decay, wbar, lw0, lw1, p_topo, p_rem, acc
      real(dp) :: la1, lg1
      integer :: livelist(Nph)
      type(vertex), pointer :: vL, vR, vn1, vn2
      if (diagram%order + 2 > maxOrder) return
      n_addB = n_addB + 1
      order0 = diagram%order
      ivn1 = order0 + 1
      ivn2 = order0 + 2
      vn1 => diagram%vertexList(ivn1)
      vn2 => diagram%vertexList(ivn2)
      call uniform_int_omp(diagram%seed, 1, order0, ivL)
      vL => diagram%vertexList(ivL)
      ivR = vL%link(3)
      vR => diagram%vertexList(ivR)
      call sample_q_omp_int(diagram%seed, vn1%i_q, vn1%Pq, vn1)
      Pq = vn1%Pq
      vn1%i_kin = vL%i_kout; vn1%ekin = vL%ekout; vn1%ukin = vL%ukout; vn1%vkin = vL%vkout
      vn1%i_kout = vn1%i_kin + vn1%i_q
      call cal_ek_int(vn1%i_kout, vn1%ekout, vn1%ukout, vn1%vkout)
      call cal_wq_int(vn1%i_q, vn1%wq)
      call cal_gkq_vtex_int(vn1, vn1%gkq)
      call cal_gkq_full_vtex_int(vn1, vn1%gkq_full)
      call live_modes(vn1%wq, nlive, wbar)
      if (nlive == 0) return
      k = 0
      do nu = 1, Nph
         if (live(vn1%wq(nu))) then
            k = k + 1
            livelist(k) = nu
         end if
      end do
      k = 1 + int(urand() * nlive)
      if (k > nlive) k = nlive
      vn1%nu = livelist(k)
      decay = wbar + minval(vn1%ekout) - minval(vn1%ekin)
      tauL = vL%tau
      tauR = vR%tau
      call uniform_real_omp(diagram%seed, tauL, tauR, tau1, temp, 1)
      call exp_sample_omp(diagram%seed, tau1, tauR, decay, tau2, Ptau12, 1)
      Ptau12 = Ptau12 / (tauR - tauL)
      if (tau2 <= tau1 .or. tau2 >= tauR) return
      p_topo = Pq * (1.d0 / order0) * Ptau12
      lw0 = lwb_cur
      ! tentative insertion (as the native acceptance branch)
      vn1%tau = tau1
      vn2%tau = tau2
      vn2%i_kin = vn1%i_kout
      vn2%i_kout = vn1%i_kin
      vn2%i_q = nk_svd - vn1%i_q
      vn2%Pq = vn1%Pq
      vn2%nu = vn1%nu
      vn2%wq = vn1%wq
      vn2%ekin = vn1%ekout; vn2%Einmin = minval(vn2%ekin)
      vn2%vkin = vn1%vkout
      vn2%ekout = vn1%ekin; vn2%Eoutmin = minval(vn2%ekout)
      vn2%vkout = vn1%vkin
      vn2%ukin = vn1%ukout
      vn2%ukout = vn1%ukin
      do nu = 1, Nph
         vn2%gkq(:, :, nu) = conjg(transpose(vn1%gkq(:, :, nu)))
         vn2%gkq_full(:, :, nu) = conjg(transpose(vn1%gkq_full(:, :, nu)))
      end do
      vn1%link = (/ivL, ivn2, ivn2/)
      vn2%link = (/ivn1, ivn1, ivR/)
      diagram%vertexList(ivL)%link(3) = ivn1
      diagram%vertexList(ivR)%link(1) = ivn2
      diagram%order = order0 + 2
      call eval_state(diagram, .true., lw1, la1, lg1)
      p_rem = 1.d0 / ((diagram%order - 1) / 2)
      ! pi_B = W_B/|G|; forward proposal carries 1/n_live for the representative mode
      acc = exp(min((lw1 - lg1) - (lw0 - lgrp_cur) + log(real(nlive, dp)), 700.d0)) * p_rem / p_topo
      if (urand() < acc) then
         n_addB_acc = n_addB_acc + 1
         lwb_cur = lw1; lwa_cur = la1; lgrp_cur = lg1
         r_snap = exp(max(min(lwb_cur - lgrp_cur - lwa_cur, 700.d0), -700.d0))
         call take_snapshot(diagram)
      else
         diagram%vertexList(ivL)%link(3) = ivR
         diagram%vertexList(ivR)%link(1) = ivL
         diagram%order = order0
      end if
   end subroutine add_B

   ! Exact reverse of add_B on a span-1 line with a live mode.
   subroutine remove_B(diagram)
      use multiphonon_update_matrix, only : delete_vertex
      type(fynman), intent(inout), target :: diagram
      integer :: iv1, iv2, ivL, ivR, nlive
      real(dp) :: Pq, tauL, tauR, tau1, tau2, Ptau12, decay, wbar, lw0, lw1, p_topo, p_rem, acc, errG
      real(dp) :: la1, lg1
      type(vertex), pointer :: v1, v2, vL, vR
      if (diagram%order <= 1) return
      call uniform_int_omp(diagram%seed, 2, diagram%order, iv1)
      v1 => diagram%vertexList(iv1)
      iv2 = v1%link(2)
      if (iv2 == 1 .or. iv2 == maxN) return
      v2 => diagram%vertexList(iv2)
      if (v2%tau < v1%tau) then
         iv1 = iv2
         v1 => diagram%vertexList(iv1)
         iv2 = v1%link(2)
         v2 => diagram%vertexList(iv2)
      end if
      if (v1%link(3) /= v1%link(2)) return
      if (.not. live(v1%wq(v1%nu))) return
      n_remB = n_remB + 1
      ivL = v1%link(1)
      ivR = v2%link(3)
      vL => diagram%vertexList(ivL)
      vR => diagram%vertexList(ivR)
      tauL = vL%tau
      tauR = vR%tau
      tau1 = v1%tau
      tau2 = v2%tau
      call live_modes(v1%wq, nlive, wbar)
      decay = wbar + minval(v1%ekout) - minval(v1%ekin)
      call exp_sample_omp(diagram%seed, tau1, tauR, decay, tau2, Ptau12)
      Ptau12 = Ptau12 / (tauR - tauL)
      Pq = v1%Pq
      p_topo = Pq * (1.d0 / (diagram%order - 2)) * Ptau12
      p_rem = 1.d0 / ((diagram%order - 1) / 2)
      lw0 = lwb_cur
      errG = maxval(abs(vR%ukin - vL%ukout)) / maxval(abs(vR%ukin))
      if (errG > 1e-10) then
         if (ivL /= 1) then
            vL%i_kout = vR%i_kin; vL%ekout = vR%ekin; vL%vkout = vR%vkin; vL%ukout = vR%ukin
         else if (ivR /= maxN) then
            vR%i_kin = vL%i_kout; vR%ekin = vL%ekout; vR%vkin = vL%vkout; vR%ukin = vL%ukout
         end if
      end if
      call delete_vertex(diagram, iv1, iv2)
      call eval_state(diagram, .true., lw1, la1, lg1)
      acc = exp(min((lw1 - lg1) - (lw0 - lgrp_cur) - log(real(nlive, dp)), 700.d0)) * p_topo / p_rem
      if (urand() < acc) then
         n_remB_acc = n_remB_acc + 1
         lwb_cur = lw1; lwa_cur = la1; lgrp_cur = lg1
         r_snap = exp(max(min(lwb_cur - lgrp_cur - lwa_cur, 700.d0), -700.d0))
         call take_snapshot(diagram)
      else
         call restore_snapshot(diagram)
      end if
   end subroutine remove_B

   ! ---- change-q move (pi_A-reversible) -------------------------------------------
   ! Returns accepted = .true. if the state changed.
   subroutine change_q(diagram, accepted)
      type(fynman), intent(inout), target :: diagram
      logical, intent(out) :: accepted
      integer :: iv1, iv2, ie, il, iu, nspan, k, nu, nu_old, ivs(maxN), qn(3)
      real(dp) :: Pq_old, Pq_new, pnu_old(Nph), pnu_new(Nph), la0, la1, lb1, lg1, acc, t0, x, c
      type(vertex), pointer :: ve, vl, vu, vprev
      accepted = .false.
      if (diagram%order < 3) return
      t0 = now()
      iv1 = 2 + int(urand() * (diagram%order - 1))
      if (iv1 > diagram%order) iv1 = diagram%order
      iv2 = diagram%vertexList(iv1)%link(2)
      if (iv2 == 1 .or. iv2 == maxN) then
         n_chq_ext = n_chq_ext + 1
         if (chq_ext_on) call change_q_ext(diagram, iv1, accepted)
         t_chq = t_chq + (now() - t0)
         return
      end if
      if (diagram%vertexList(iv2)%tau < diagram%vertexList(iv1)%tau) then
         ie = iv2; il = iv1
      else
         ie = iv1; il = iv2
      end if
      ve => diagram%vertexList(ie)
      vl => diagram%vertexList(il)
      nspan = 1
      ivs(1) = ie
      iu = ve%link(3)
      do while (iu /= il)
         if (iu == maxN) stop 'change_q: partner not reached along the chain'
         nspan = nspan + 1
         ivs(nspan) = iu
         iu = diagram%vertexList(iu)%link(3)
      end do
      nspan = nspan + 1
      ivs(nspan) = il
      n_chq = n_chq + 1
      sum_span_try = sum_span_try + (nspan - 1)
      call eval_state(diagram, .false., lb1, la0, lg1)
      do k = 1, nspan
         call copy_vtex(Nph, dmc_band, nbnd, nsvd, diagram%vertexList(ivs(k)), csnap(k))
      end do
      Pq_old = ve%Pq
      nu_old = ve%nu
      do nu = 1, Nph
         pnu_old(nu) = maxval(abs(ve%gkq(:, :, nu)))**2
      end do
      pnu_old = pnu_old / sum(pnu_old)
      ! proposal: new momentum on the earlier vertex, then the mode
      call sample_q_omp_int(diagram%seed, qn, Pq_new, ve)
      ve%i_q = qn
      ve%Pq = Pq_new
      ve%i_kout = modulo(ve%i_kin + ve%i_q, nk_svd)
      call cal_ek_int(ve%i_kout, ve%ekout, ve%ukout, ve%vkout)
      ve%Eoutmin = minval(ve%ekout)
      call cal_wq_int(ve%i_q, ve%wq)
      call cal_gkq_vtex_int(ve, ve%gkq)
      call cal_gkq_full_vtex_int(ve, ve%gkq_full)
      do nu = 1, Nph
         pnu_new(nu) = maxval(abs(ve%gkq(:, :, nu)))**2
      end do
      if (sum(pnu_new) <= 0.d0) then
         call chq_restore(diagram, ivs, nspan)
         t_chq = t_chq + (now() - t0)
         return
      end if
      pnu_new = pnu_new / sum(pnu_new)
      x = urand()
      c = 0.d0
      nu = Nph
      do k = 1, Nph
         c = c + pnu_new(k)
         if (x < c) then
            nu = k
            exit
         end if
      end do
      ve%nu = nu
      ! shift the span: interior vertices keep their own q, the partner takes -q'
      do k = 2, nspan
         vu => diagram%vertexList(ivs(k))
         vprev => diagram%vertexList(ivs(k - 1))
         vu%i_kin = vprev%i_kout; vu%ekin = vprev%ekout; vu%ukin = vprev%ukout; vu%vkin = vprev%vkout
         vu%Einmin = minval(vu%ekin)
         if (k < nspan) then
            vu%i_kout = modulo(vu%i_kin + vu%i_q, nk_svd)
            call cal_ek_int(vu%i_kout, vu%ekout, vu%ukout, vu%vkout)
            vu%Eoutmin = minval(vu%ekout)
         else
            vu%i_q = nk_svd - ve%i_q
            vu%Pq = ve%Pq
            vu%nu = ve%nu
            vu%wq = ve%wq
         end if
         call cal_gkq_vtex_int(vu, vu%gkq)
         call cal_gkq_full_vtex_int(vu, vu%gkq_full)
      end do
      call eval_state(diagram, .false., lb1, la1, lg1)
      acc = exp(max(min(la1 - la0, 700.d0), -700.d0)) * (Pq_old * pnu_old(nu_old)) / (Pq_new * pnu_new(nu))
      if (urand() < acc) then
         n_chq_acc = n_chq_acc + 1
         sum_span_acc = sum_span_acc + (nspan - 1)
         accepted = .true.
      else
         call chq_restore(diagram, ivs, nspan)
      end if
      t_chq = t_chq + (now() - t0)
   end subroutine change_q

   ! External (boundary-wrapping) pair: vertex eh (link(2) = 1) and its plink et
   ! (link(2) = maxN) carry q and -q; the phonon runs et -> tau_max = 0 -> eh, so
   ! every segment outside [eh, et] (head .. eh, et .. tail) carries -q. New q' ~ Pq
   ! on eh, nu' ~ Pnu(eh) (as add_external_ph); those segments shift by -(q' - q).
   ! Head, tail and the last vertex share one eigenvector set (the periodic segment),
   ! exactly as the native code keeps it.
   subroutine change_q_ext(diagram, ivx, accepted)
      type(fynman), intent(inout), target :: diagram
      integer, intent(in) :: ivx
      logical, intent(out) :: accepted
      integer :: ieh, iet, iu, k, n, nu, nu_old, ivs(maxN), qn(3), q_old(3)
      real(dp) :: Pq_old, Pq_new, pnu_old(Nph), pnu_new(Nph), pnu_rt(Nph), la0, la1, la_rt, lb1, lg1, acc, x, c
      logical :: ok
      type(vertex), pointer :: veh, vet, vh, vt
      accepted = .false.
      if (diagram%vertexList(ivx)%link(2) == 1) then
         ieh = ivx
         iet = diagram%vertexList(ivx)%plink
      else
         iet = ivx
         ieh = diagram%vertexList(ivx)%plink
      end if
      if (ieh < 2 .or. ieh > diagram%order .or. iet < 2 .or. iet > diagram%order) return
      veh => diagram%vertexList(ieh)
      vet => diagram%vertexList(iet)
      if (veh%link(2) /= 1 .or. vet%link(2) /= maxN .or. vet%plink /= ieh) return
      vh => diagram%vertexList(1)
      vt => diagram%vertexList(maxN)
      ! the vertices of the outer region, in chain order: head side then tail side
      n = 0
      iu = vh%link(3)
      do while (iu /= ieh)
         if (iu == iet .or. iu == maxN) return      ! et before eh: leave the state alone
         n = n + 1
         ivs(n) = iu
         iu = diagram%vertexList(iu)%link(3)
      end do
      n = n + 1
      ivs(n) = ieh
      iu = veh%link(3)
      do while (iu /= iet)
         if (iu == maxN) return
         iu = diagram%vertexList(iu)%link(3)
      end do
      do while (iu /= maxN)
         n = n + 1
         ivs(n) = iu
         iu = diagram%vertexList(iu)%link(3)
      end do
      n_chq_ext_try = n_chq_ext_try + 1
      call eval_state(diagram, .false., lb1, la0, lg1)
      do k = 1, n
         call copy_vtex(Nph, dmc_band, nbnd, nsvd, diagram%vertexList(ivs(k)), csnap(k))
      end do
      call copy_vtex(Nph, dmc_band, nbnd, nsvd, vh, csnap(n + 1))
      call copy_vtex(Nph, dmc_band, nbnd, nsvd, vt, csnap(n + 2))
      Pq_old = veh%Pq
      nu_old = veh%nu
      q_old = veh%i_q
      do nu = 1, Nph
         pnu_old(nu) = maxval(abs(veh%gkq(:, :, nu)))**2
      end do
      pnu_old = pnu_old / sum(pnu_old)
      call sample_q_omp_int(diagram%seed, qn, Pq_new, veh)
      call ext_shift(diagram, ivs, n, ieh, iet, qn, Pq_new, pnu_new, ok)
      if (.not. ok) then
         call chq_restore_ext(diagram, ivs, n)
         return
      end if
      x = urand()
      c = 0.d0
      nu = Nph
      do k = 1, Nph
         c = c + pnu_new(k)
         if (x < c) then
            nu = k
            exit
         end if
      end do
      veh%nu = nu
      vet%nu = nu
      call eval_state(diagram, .false., lb1, la1, lg1)
      ! round trip C -> C' -> C'': the reverse proposal must see exactly pnu_old and w(C)
      if (chq_rt .and. mod(n_chq_ext_try, 20_8) == 0) then
         do k = 1, n
            call copy_vtex(Nph, dmc_band, nbnd, nsvd, diagram%vertexList(ivs(k)), csnap2(k))
         end do
         call copy_vtex(Nph, dmc_band, nbnd, nsvd, vh, csnap2(n + 1))
         call copy_vtex(Nph, dmc_band, nbnd, nsvd, vt, csnap2(n + 2))
         call ext_shift(diagram, ivs, n, ieh, iet, q_old, Pq_old, pnu_rt, ok)
         veh%nu = nu_old
         vet%nu = nu_old
         call eval_state(diagram, .false., lb1, la_rt, lg1)
         n_rt = n_rt + 1
         rt_la = max(rt_la, abs(la_rt - la0))
         rt_pnu = max(rt_pnu, maxval(abs(pnu_rt - pnu_old)))
         do k = 1, n
            rt_g = max(rt_g, maxval(abs(diagram%vertexList(ivs(k))%gkq - csnap(k)%gkq)))
            rt_u = max(rt_u, maxval(abs(diagram%vertexList(ivs(k))%ukout - csnap(k)%ukout)))
         end do
         do k = 1, n
            call copy_vtex(Nph, dmc_band, nbnd, nsvd, csnap2(k), diagram%vertexList(ivs(k)))
         end do
         call copy_vtex(Nph, dmc_band, nbnd, nsvd, csnap2(n + 1), vh)
         call copy_vtex(Nph, dmc_band, nbnd, nsvd, csnap2(n + 2), vt)
      end if
      acc = exp(max(min(la1 - la0, 700.d0), -700.d0)) * (Pq_old * pnu_old(nu_old)) / (Pq_new * pnu_new(nu))
      if (urand() < acc) then
         n_chq_ext_acc = n_chq_ext_acc + 1
         accepted = .true.
      else
         call chq_restore_ext(diagram, ivs, n)
      end if
   end subroutine change_q_ext

   ! Put momentum qn on the head-side vertex of the pair and shift the outer region;
   ! returns the normalized Pnu of the new head-side vertex (ok = .false. if all zero).
   ! Modes are not touched.
   subroutine ext_shift(diagram, ivs, n, ieh, iet, qn, Pqn, pnu, ok)
      type(fynman), intent(inout), target :: diagram
      integer, intent(in) :: ivs(:), n, ieh, iet, qn(3)
      real(dp), intent(in) :: Pqn
      real(dp), intent(out) :: pnu(Nph)
      logical, intent(out) :: ok
      integer :: k, nu, dq(3), kh(3)
      real(dp) :: ekh(dmc_band), vkh(3, dmc_band)
      complex(dp) :: ukh(nbnd, nbnd)
      type(vertex), pointer :: veh, vu, vprev, vh, vt
      veh => diagram%vertexList(ieh)
      vh => diagram%vertexList(1)
      vt => diagram%vertexList(maxN)
      ok = .true.
      dq = qn - veh%i_q
      ! new periodic-segment momentum and its eigen-system, shared by head, tail, last vertex
      kh = modulo(vh%i_kout - dq, nk_svd)
      call cal_ek_int(kh, ekh, ukh, vkh)
      vh%i_kout = kh; vh%ekout = ekh; vh%ukout = ukh; vh%vkout = vkh
      vh%Eoutmin = minval(ekh)
      do k = 1, n
         vu => diagram%vertexList(ivs(k))
         vprev => diagram%vertexList(vu%link(1))
         if (vprev%link(3) /= ivs(k)) stop 'change_q_ext: link(1) is not the chain predecessor'
         if (ivs(k) /= iet) then
            vu%i_kin = vprev%i_kout; vu%ekin = vprev%ekout; vu%ukin = vprev%ukout; vu%vkin = vprev%vkout
            vu%Einmin = minval(vu%ekin)
         end if
         if (ivs(k) == ieh) then
            vu%i_q = qn
            vu%Pq = Pqn
            call cal_wq_int(vu%i_q, vu%wq)
            call cal_gkq_vtex_int(vu, vu%gkq)
            call cal_gkq_full_vtex_int(vu, vu%gkq_full)
            do nu = 1, Nph
               pnu(nu) = maxval(abs(vu%gkq(:, :, nu)))**2
            end do
            if (sum(pnu) <= 0.d0) then
               ok = .false.
               return
            end if
            pnu = pnu / sum(pnu)
            cycle
         end if
         if (ivs(k) == iet) then
            vu%i_q = -veh%i_q
            vu%Pq = veh%Pq; vu%wq = veh%wq
         end if
         vu%i_kout = modulo(vu%i_kout - dq, nk_svd)
         if (vu%link(3) == maxN) then
            vu%i_kout = kh; vu%ekout = ekh; vu%ukout = ukh; vu%vkout = vkh
         else
            call cal_ek_int(vu%i_kout, vu%ekout, vu%ukout, vu%vkout)
         end if
         vu%Eoutmin = minval(vu%ekout)
         call cal_gkq_vtex_int(vu, vu%gkq)
         call cal_gkq_full_vtex_int(vu, vu%gkq_full)
      end do
      vt%i_kin = kh; vt%ekin = ekh; vt%ukin = ukh; vt%vkin = vkh
      vt%Einmin = minval(ekh)
   end subroutine ext_shift

   subroutine chq_restore_ext(diagram, ivs, n)
      type(fynman), intent(inout), target :: diagram
      integer, intent(in) :: ivs(:), n
      integer :: k
      do k = 1, n
         call copy_vtex(Nph, dmc_band, nbnd, nsvd, csnap(k), diagram%vertexList(ivs(k)))
      end do
      call copy_vtex(Nph, dmc_band, nbnd, nsvd, csnap(n + 1), diagram%vertexList(1))
      call copy_vtex(Nph, dmc_band, nbnd, nsvd, csnap(n + 2), diagram%vertexList(maxN))
   end subroutine chq_restore_ext

   subroutine chq_restore(diagram, ivs, nspan)
      type(fynman), intent(inout), target :: diagram
      integer, intent(in) :: ivs(:), nspan
      integer :: k
      do k = 1, nspan
         call copy_vtex(Nph, dmc_band, nbnd, nsvd, csnap(k), diagram%vertexList(ivs(k)))
      end do
   end subroutine chq_restore

   ! Hermiticity of the tabulated vertex matrices: native sets a partner vertex's matrix
   ! to conjg(transpose(g)) of the first vertex, while direct evaluation computes
   ! g(k+q, -q) from the tables. Samples random (k, q) with q ~ Pq and reports the
   ! relative deviation (FUTUREB_HERM_TEST=1, written to herm_test.dat).
   subroutine herm_test(diagram)
      type(fynman), intent(inout), target :: diagram
      type(vertex) :: va, vb
      integer :: it, nu, u, k(3), q(3), nsamp
      real(dp) :: x, Pq, d, gmax, dmax, dsum, wsum, dw, pa(Nph), pb(Nph), dp_max
      complex(dp) :: gt(dmc_band, dmc_band)
      call CreateVertex(Nph, dmc_band, nbnd, nsvd, va)
      call CreateVertex(Nph, dmc_band, nbnd, nsvd, vb)
      nsamp = 4000
      dmax = 0; dsum = 0; wsum = 0; dp_max = 0
      do it = 1, nsamp
         do u = 1, 3
            k(u) = int(urand() * nk_svd)
         end do
         call sample_q_omp_int(diagram%seed, q, Pq, va)
         va%i_kin = k
         va%i_q = q
         va%i_kout = modulo(k + q, nk_svd)
         call cal_ek_int(va%i_kin, va%ekin, va%ukin, va%vkin)
         call cal_ek_int(va%i_kout, va%ekout, va%ukout, va%vkout)
         call cal_gkq_vtex_int(va, va%gkq)
         vb%i_kin = va%i_kout
         vb%i_q = nk_svd - q
         vb%i_kout = va%i_kin
         vb%ukin = va%ukout; vb%ukout = va%ukin; vb%ekin = va%ekout; vb%ekout = va%ekin
         call cal_gkq_vtex_int(vb, vb%gkq)
         do nu = 1, Nph
            gt = conjg(transpose(va%gkq(:, :, nu)))
            gmax = maxval(abs(va%gkq(:, :, nu)))
            if (gmax <= 0.d0) cycle
            d = maxval(abs(vb%gkq(:, :, nu) - gt)) / gmax
            dmax = max(dmax, d)
            dw = gmax**2
            dsum = dsum + d * dw
            wsum = wsum + dw
         end do
         do nu = 1, Nph
            pa(nu) = maxval(abs(va%gkq(:, :, nu)))**2
            pb(nu) = maxval(abs(vb%gkq(:, :, nu)))**2
         end do
         if (sum(pa) > 0 .and. sum(pb) > 0) dp_max = max(dp_max, maxval(abs(pa / sum(pa) - pb / sum(pb))))
      end do
      open(newunit=u, file='herm_test.dat', status='replace', action='write')
      write(u, '(A,I8)') 'samples ', nsamp
      write(u, '(A,ES12.4)') 'max_rel_dev_g ', dmax
      write(u, '(A,ES12.4)') 'weighted_mean_rel_dev_g ', dsum / max(wsum, 1.d-300)
      write(u, '(A,ES12.4)') 'max_dev_Pnu ', dp_max
      close(u)
   end subroutine herm_test

   ! Momentum conservation at every internal vertex and continuity along the chain.
   subroutine k_check(diagram)
      type(fynman), intent(inout), target :: diagram
      integer :: iv, nx
      logical :: bad
      type(vertex), pointer :: v
      bad = .false.
      iv = 1
      do while (iv /= maxN)
         v => diagram%vertexList(iv)
         nx = v%link(3)
         if (iv /= 1) bad = bad .or. any(modulo(v%i_kin + v%i_q - v%i_kout, nk_svd) /= 0)
         if (nx /= maxN) bad = bad .or. any(modulo(diagram%vertexList(nx)%i_kin - v%i_kout, nk_svd) /= 0)
         if (nx == maxN) bad = bad .or. any(modulo(diagram%vertexList(maxN)%i_kin - v%i_kout, nk_svd) /= 0)
         iv = nx
      end do
      ! periodic segment: head k_out = tail k_in
      bad = bad .or. any(modulo(diagram%vertexList(1)%i_kout - diagram%vertexList(maxN)%i_kin, nk_svd) /= 0)
      n_kchecks = n_kchecks + 1
      if (bad) n_kbad = n_kbad + 1
   end subroutine k_check

   ! Gauge consistency: every segment's eigenvectors are the same on both ends,
   ! and the periodic segment (last -> tail = head -> first) uses the head's.
   subroutine g_check(diagram)
      type(fynman), intent(inout), target :: diagram
      integer :: iv, nx
      logical :: bad
      type(vertex), pointer :: v
      bad = .false.
      iv = 1
      do while (iv /= maxN)
         v => diagram%vertexList(iv)
         nx = v%link(3)
         if (nx /= maxN) then
            bad = bad .or. maxval(abs(diagram%vertexList(nx)%ukin - v%ukout)) > 1.d-12
         else
            bad = bad .or. maxval(abs(v%ukout - diagram%vertexList(1)%ukout)) > 1.d-12
         end if
         iv = nx
      end do
      bad = bad .or. maxval(abs(diagram%vertexList(maxN)%ukin - diagram%vertexList(1)%ukout)) > 1.d-12
      n_gchecks = n_gchecks + 1
      if (bad) n_gbad = n_gbad + 1
   end subroutine g_check

   ! ---- general-span add/remove of internal lines (pi_A-reversible pair) ----------
   subroutine full_snapshot(diagram)
      type(fynman), intent(inout), target :: diagram
      integer :: iv
      do iv = 1, diagram%order
         call copy_vtex(Nph, dmc_band, nbnd, nsvd, diagram%vertexList(iv), csnap(iv))
      end do
      call copy_vtex(Nph, dmc_band, nbnd, nsvd, diagram%vertexList(maxN), csnap(maxN))
      fs_order = diagram%order
      fs_nph_ext = diagram%nph_ext
   end subroutine full_snapshot

   subroutine full_restore(diagram)
      type(fynman), intent(inout), target :: diagram
      integer :: iv
      do iv = 1, fs_order
         call copy_vtex(Nph, dmc_band, nbnd, nsvd, csnap(iv), diagram%vertexList(iv))
      end do
      call copy_vtex(Nph, dmc_band, nbnd, nsvd, csnap(maxN), diagram%vertexList(maxN))
      diagram%order = fs_order
      diagram%nph_ext = fs_nph_ext
   end subroutine full_restore

   ! mode proposal from a vertex's matrices: native max|g_ij|^2, or the
   ! gauge-invariant Frobenius norm sum |g_ij|^2 (FUTUREB_PNU_FRO=1)
   subroutine mode_probs(g, pnu)
      complex(dp), intent(in) :: g(:, :, :)
      real(dp), intent(out) :: pnu(Nph)
      integer :: nu
      do nu = 1, Nph
         if (pnu_fro) then
            pnu(nu) = sum(abs(g(:, :, nu))**2)
         else
            pnu(nu) = maxval(abs(g(:, :, nu)))**2
         end if
      end do
      if (sum(pnu) > 0.d0) pnu = pnu / sum(pnu)
   end subroutine mode_probs

   ! Native exp_sample_omp semantics: on [a, b], density ~ exp(-|beta| xdis) with
   ! xdis = x - a (beta > 0) or b - x (beta < 0). sample = .true. draws x.
   subroutine nexp(a, b, beta, x, px, sample)
      real(dp), intent(in) :: a, b, beta
      real(dp), intent(inout) :: x
      real(dp), intent(out) :: px
      logical, intent(in) :: sample
      real(dp) :: L, yb, xdis, ab, xt
      L = b - a
      ab = abs(beta)
      if (ab * L < 1.d-10) then
         if (sample) x = a + urand() * L
         px = 1.d0 / L
         return
      end if
      yb = exp(-ab * L)
      if (sample) then
         xt = urand()
         xdis = abs(log(xt + yb * (1.d0 - xt)) / ab)
         xdis = min(xdis, L)
         if (beta > 0) then
            x = a + xdis
         else
            x = b - xdis
         end if
      end if
      if (beta > 0) then
         xdis = x - a
      else
         xdis = b - x
      end if
      px = ab * exp(-ab * xdis) / (1.d0 - yb)
   end subroutine nexp

   ! density of tau2 on (tau1, tmax) for rate lam (truncated exponential)
   real(dp) function ftau2(tau1, tau2, tmax, lam)
      real(dp), intent(in) :: tau1, tau2, tmax, lam
      real(dp) :: L
      L = tmax - tau1
      if (lam * L < 1.d-8) then
         ftau2 = 1.d0 / L
      else
         ftau2 = lam * exp(-lam * (tau2 - tau1)) / (1.d0 - exp(-lam * L))
      end if
   end function ftau2

   subroutine copy_in_from_out(vu, vsrc)
      type(vertex), intent(inout) :: vu
      type(vertex), intent(in) :: vsrc
      vu%i_kin = vsrc%i_kout; vu%ekin = vsrc%ekout; vu%ukin = vsrc%ukout; vu%vkin = vsrc%vkout
      vu%Einmin = minval(vu%ekin)
   end subroutine copy_in_from_out

   subroutine copy_out_from_in(vu, vsrc)
      type(vertex), intent(inout) :: vu
      type(vertex), intent(in) :: vsrc
      vu%i_kout = vsrc%i_kin; vu%ekout = vsrc%ekin; vu%ukout = vsrc%ukin; vu%vkout = vsrc%vkin
      vu%Eoutmin = minval(vu%ekout)
   end subroutine copy_out_from_in

   ! Shift the momentum of every vertex marked in shiftv by dk (its out segment), keep
   ! each segment's eigen-system identical on both ends, give the periodic segment
   ! (head out = last out = tail in) the head's eigen-system, and recompute vertex
   ! matrices wherever an in or out segment changed.
   subroutine resync_shift(diagram, dk, shiftv)
      type(fynman), intent(inout), target :: diagram
      integer, intent(in) :: dk(3)
      logical, intent(in) :: shiftv(:)
      integer :: iv, ip, kh(3)
      logical :: prev_changed, need
      type(vertex), pointer :: v, vp, vh, vt
      vh => diagram%vertexList(1)
      vt => diagram%vertexList(maxN)
      kh = modulo(vh%i_kout + dk, nk_svd)
      call cal_ek_int(kh, vh%ekout, vh%ukout, vh%vkout)
      vh%i_kout = kh
      vh%Eoutmin = minval(vh%ekout)
      prev_changed = .true.
      ip = 1
      iv = vh%link(3)
      do while (iv /= maxN)
         v => diagram%vertexList(iv)
         vp => diagram%vertexList(ip)
         need = prev_changed
         call copy_in_from_out(v, vp)
         prev_changed = .false.
         if (shiftv(iv)) then
            v%i_kout = modulo(v%i_kout + dk, nk_svd)
            if (v%link(3) == maxN) then
               v%i_kout = vh%i_kout; v%ekout = vh%ekout; v%ukout = vh%ukout; v%vkout = vh%vkout
            else
               call cal_ek_int(v%i_kout, v%ekout, v%ukout, v%vkout)
            end if
            v%Eoutmin = minval(v%ekout)
            need = .true.
            prev_changed = .true.
         else if (v%link(3) == maxN) then
            v%i_kout = vh%i_kout; v%ekout = vh%ekout; v%ukout = vh%ukout; v%vkout = vh%vkout
            v%Eoutmin = minval(v%ekout)
            need = .true.
         end if
         if (need) then
            call cal_gkq_vtex_int(v, v%gkq)
            call cal_gkq_full_vtex_int(v, v%gkq_full)
         end if
         ip = iv
         iv = v%link(3)
      end do
      vt%i_kin = vh%i_kout; vt%ekin = vh%ekout; vt%ukin = vh%ukout; vt%vkin = vh%vkout
      vt%Einmin = minval(vt%ekin)
   end subroutine resync_shift

   ! ---- general external add/remove (boundary-wrapping pairs at any position) -----
   subroutine ext_move(diagram, accepted)
      type(fynman), intent(inout), target :: diagram
      logical, intent(out) :: accepted
      real(dp) :: t0
      t0 = now()
      if (urand() < 0.5d0) then
         call ext_add(diagram, accepted)
      else
         call ext_remove(diagram, accepted)
      end if
      t_ext = t_ext + (now() - t0)
   end subroutine ext_move

   ! head-side vertex eh at tau1 ~ U(0, tmax), q ~ Pq, nu ~ Pnu(eh); tail-side vertex
   ! et at tau2 with tmax - tau2 ~ truncated Exp(omega_nu) on (0, tmax - tau1).
   ! Every segment outside [tau1, tau2] (head side and tail side) carries -q.
   subroutine ext_add(diagram, accepted)
      type(fynman), intent(inout), target :: diagram
      logical, intent(out) :: accepted
      integer :: ipa, ina, ipb, inb, ieh, iet, iv, k, nu, qn(3)
      real(dp) :: tmax, tau1, tau2, Pqn, pnu(Nph), lam, L, x, c, ft, la0, la1, lb1, lg1, acc, p_r, p_a, u
      real(dp) :: tauc, pt1, pt2
      integer :: ilast
      logical :: shiftv(maxN), before
      type(vertex), pointer :: veh, vet, vpa
      accepted = .false.
      if (diagram%order + 2 > maxOrder) return
      n_extA = n_extA + 1
      tmax = diagram%vertexList(maxN)%tau
      if (ext_mimic) then
         ipa = 1                     ! native inserts right after the head; tau1 drawn later
         tau1 = 0.d0
      else
         tau1 = urand() * tmax
         ipa = 1
         do while (diagram%vertexList(ipa)%link(3) /= maxN)
            if (diagram%vertexList(diagram%vertexList(ipa)%link(3))%tau >= tau1) exit
            ipa = diagram%vertexList(ipa)%link(3)
         end do
      end if
      ina = diagram%vertexList(ipa)%link(3)
      if (.not. ext_mimic) then
         if (ext_rpos .and. ipa /= 1) return
         if (ext_rtau .and. tau1 >= 0.5d0 * tmax) return
      end if
      call eval_state(diagram, .false., lb1, la0, lg1)
      call full_snapshot(diagram)
      ieh = diagram%order + 1
      iet = diagram%order + 2
      veh => diagram%vertexList(ieh)
      vet => diagram%vertexList(iet)
      vpa => diagram%vertexList(ipa)
      veh%tau = tau1
      veh%i_kout = vpa%i_kout; veh%ekout = vpa%ekout; veh%ukout = vpa%ukout; veh%vkout = vpa%vkout
      veh%Eoutmin = minval(veh%ekout)
      call sample_q_omp_int(diagram%seed, qn, Pqn, veh)
      veh%i_q = qn
      veh%Pq = Pqn
      veh%i_kin = modulo(veh%i_kout - qn, nk_svd)
      call cal_ek_int(veh%i_kin, veh%ekin, veh%ukin, veh%vkin)
      veh%Einmin = minval(veh%ekin)
      call cal_wq_int(veh%i_q, veh%wq)
      call cal_gkq_vtex_int(veh, veh%gkq)
      call cal_gkq_full_vtex_int(veh, veh%gkq_full)
      call mode_probs(veh%gkq, pnu)
      if (sum(pnu) <= 0.d0) then
         call full_restore(diagram)
         return
      end if
      x = urand()
      c = 0.d0
      nu = Nph
      do k = 1, Nph
         c = c + pnu(k)
         if (x < c) then
            nu = k
            exit
         end if
      end do
      veh%nu = nu
      if (ext_mimic) then
         ! native add_external_ph proposal: tau1 on [0, min(tau_first, tmax/2)], tau2 on
         ! [max(tau_last, tmax/2), tmax], signed decays from the pair's energies and omega
         if (ina == maxN) then
            tauc = 0.5d0 * tmax
         else
            tauc = min(diagram%vertexList(ina)%tau, 0.5d0 * tmax)
         end if
         call nexp(0.d0, tauc, minval(veh%ekin) + veh%wq(nu) - minval(veh%ekout), tau1, pt1, .true.)
         veh%tau = tau1
         ilast = diagram%vertexList(maxN)%link(1)
         if (ilast == 1) then
            tauc = 0.5d0 * tmax
         else
            tauc = max(diagram%vertexList(ilast)%tau, 0.5d0 * tmax)
         end if
         ! et: in = last segment (momentum k, = head out before the shift), out = k - q
         call nexp(tauc, tmax, minval(diagram%vertexList(ilast)%ekout) - veh%wq(nu) - minval(veh%ekin), &
                   tau2, pt2, .true.)
         ft = pt1 * pt2 * tmax          ! p_a below carries 1/tmax
      else
         lam = max(veh%wq(nu), 0.d0)
         L = tmax - tau1
         x = urand()
         if (lam * L < 1.d-8) then
            u = x * L
         else
            u = -log(1.d0 - x * (1.d0 - exp(-lam * L))) / lam
         end if
         tau2 = tmax - u
         ft = ftau2(0.d0, u, L, lam)
      end if
      if (.not. (tau2 > tau1 .and. tau2 < tmax)) then
         call full_restore(diagram)
         return
      end if
      ipb = ipa
      do while (diagram%vertexList(ipb)%link(3) /= maxN)
         if (diagram%vertexList(diagram%vertexList(ipb)%link(3))%tau >= tau2) exit
         ipb = diagram%vertexList(ipb)%link(3)
      end do
      inb = diagram%vertexList(ipb)%link(3)
      if ((ext_rpos .and. inb /= maxN) .or. (ext_rtau .and. tau2 <= 0.5d0 * tmax)) then
         call full_restore(diagram)
         return
      end if
      if (ext_ordered) then
         ! no tail-attached vertex before tau1, no head-attached vertex after tau2
         iv = diagram%vertexList(1)%link(3)
         do while (iv /= maxN)
            if (diagram%vertexList(iv)%tau < tau1 .and. diagram%vertexList(iv)%link(2) == maxN .or. &
                diagram%vertexList(iv)%tau > tau2 .and. diagram%vertexList(iv)%link(2) == 1) then
               call full_restore(diagram)
               return
            end if
            iv = diagram%vertexList(iv)%link(3)
         end do
      end if
      vet%tau = tau2
      if (ipb == ipa) then
         vet%i_kin = veh%i_kout; vet%ekin = veh%ekout; vet%ukin = veh%ukout; vet%vkin = veh%vkout
      else
         call copy_in_from_out(vet, diagram%vertexList(ipb))
      end if
      vet%Einmin = minval(vet%ekin)
      vet%i_q = -qn
      vet%Pq = Pqn
      vet%nu = nu
      vet%wq = veh%wq
      vet%i_kout = vet%i_kin
      vet%ekout = vet%ekin; vet%ukout = vet%ukin; vet%vkout = vet%vkin
      ! link: eh after ipa, et after ipb (eh -> et directly if they share a segment)
      if (ipb == ipa) then
         veh%link = (/ipa, 1, iet/)
         vet%link = (/ieh, maxN, ina/)
         diagram%vertexList(ipa)%link(3) = ieh
         diagram%vertexList(ina)%link(1) = iet
      else
         veh%link = (/ipa, 1, ina/)
         vet%link = (/ipb, maxN, inb/)
         diagram%vertexList(ipa)%link(3) = ieh
         diagram%vertexList(ina)%link(1) = ieh
         diagram%vertexList(ipb)%link(3) = iet
         diagram%vertexList(inb)%link(1) = iet
      end if
      veh%plink = iet
      vet%plink = ieh
      diagram%order = diagram%order + 2
      diagram%nph_ext = diagram%nph_ext + 1
      ! outside region: vertices before eh, et itself (its out segment) and after et
      shiftv = .false.
      before = .true.
      iv = diagram%vertexList(1)%link(3)
      do while (iv /= maxN)
         if (iv == ieh) before = .false.
         if (before) shiftv(iv) = .true.
         if (iv == iet) before = .true.
         iv = diagram%vertexList(iv)%link(3)
      end do
      shiftv(ieh) = .false.
      shiftv(iet) = .true.
      call resync_shift(diagram, -qn, shiftv)
      call eval_state(diagram, .false., lb1, la1, lg1)
      p_r = 0.5d0 * 2.d0 / real(diagram%order - 1, dp)
      p_a = 0.5d0 * (1.d0 / tmax) * Pqn * pnu(nu) * ft
      if (chq_rt .and. mod(n_extA, 20_8) == 0) call ext_roundtrip(diagram, ieh, iet, la0, p_a, p_r)
      acc = exp(max(min(la1 - la0, 700.d0), -700.d0)) * p_r / p_a
      if (urand() < acc) then
         n_extA_acc = n_extA_acc + 1
         accepted = .true.
      else
         call full_restore(diagram)
      end if
   end subroutine ext_add

   subroutine ext_remove(diagram, accepted)
      type(fynman), intent(inout), target :: diagram
      logical, intent(out) :: accepted
      integer :: iv1, ieh, iet, iv, order_c
      real(dp) :: la0, la1, lb1, lg1, acc, p_r, p_a
      logical :: seen_eh
      accepted = .false.
      if (diagram%order < 3) return
      iv1 = 2 + int(urand() * (diagram%order - 1))
      if (iv1 > diagram%order) iv1 = diagram%order
      if (diagram%vertexList(iv1)%link(2) == 1) then
         ieh = iv1
         iet = diagram%vertexList(iv1)%plink
      else if (diagram%vertexList(iv1)%link(2) == maxN) then
         iet = iv1
         ieh = diagram%vertexList(iv1)%plink
      else
         return
      end if
      if (ieh < 2 .or. ieh > diagram%order .or. iet < 2 .or. iet > diagram%order) return
      if (diagram%vertexList(ieh)%link(2) /= 1 .or. diagram%vertexList(iet)%link(2) /= maxN) return
      if (diagram%vertexList(iet)%plink /= ieh .or. diagram%vertexList(ieh)%plink /= iet) return
      seen_eh = .false.
      iv = diagram%vertexList(1)%link(3)
      do while (iv /= maxN)
         if (iv == ieh) seen_eh = .true.
         if (iv == iet) exit
         iv = diagram%vertexList(iv)%link(3)
      end do
      if (.not. seen_eh) return
      if (ext_rpos) then
         if (diagram%vertexList(1)%link(3) /= ieh .or. diagram%vertexList(maxN)%link(1) /= iet) return
      end if
      if (ext_rtau) then
         if (diagram%vertexList(ieh)%tau >= 0.5d0 * diagram%vertexList(maxN)%tau) return
         if (diagram%vertexList(iet)%tau <= 0.5d0 * diagram%vertexList(maxN)%tau) return
      end if
      n_extR = n_extR + 1
      order_c = diagram%order
      call eval_state(diagram, .false., lb1, la0, lg1)
      call full_snapshot(diagram)
      call ext_remove_core(diagram, ieh, iet, la1, p_a, p_r)
      acc = exp(max(min(la1 - la0, 700.d0), -700.d0)) * p_a / p_r
      if (urand() < acc) then
         n_extR_acc = n_extR_acc + 1
         accepted = .true.
      else
         call full_restore(diagram)
      end if
   end subroutine ext_remove

   ! Round trip C -> C' (just proposed) -> remove the same pair: must give back log w(C)
   ! and the reverse proposal the add used. The proposed state C' is restored after.
   subroutine ext_roundtrip(diagram, ieh, iet, la0, p_a, p_r)
      type(fynman), intent(inout), target :: diagram
      integer, intent(in) :: ieh, iet
      real(dp), intent(in) :: la0, p_a, p_r
      integer :: iv, ord2, next2
      real(dp) :: la_rt, pa_rt, pr_rt
      ord2 = diagram%order
      next2 = diagram%nph_ext
      do iv = 1, ord2
         call copy_vtex(Nph, dmc_band, nbnd, nsvd, diagram%vertexList(iv), csnap2(iv))
      end do
      call copy_vtex(Nph, dmc_band, nbnd, nsvd, diagram%vertexList(maxN), csnap2(maxN))
      call ext_remove_core(diagram, ieh, iet, la_rt, pa_rt, pr_rt)
      n_rt_ext = n_rt_ext + 1
      rt_ext_la = max(rt_ext_la, abs(la_rt - la0))
      rt_ext_pa = max(rt_ext_pa, abs(pa_rt / p_a - 1.d0))
      rt_ext_pr = max(rt_ext_pr, abs(pr_rt / p_r - 1.d0))
      do iv = 1, ord2
         call copy_vtex(Nph, dmc_band, nbnd, nsvd, csnap2(iv), diagram%vertexList(iv))
      end do
      call copy_vtex(Nph, dmc_band, nbnd, nsvd, csnap2(maxN), diagram%vertexList(maxN))
      diagram%order = ord2
      diagram%nph_ext = next2
   end subroutine ext_roundtrip

   ! Remove the valid external pair (ieh before iet) without an acceptance step.
   ! Returns log w of the new state, the reverse-add proposal p_a and p_r.
   subroutine ext_remove_core(diagram, ieh_in, iet_in, la1, p_a, p_r)
      use multiphonon_update_matrix, only : swap_vertex
      type(fynman), intent(inout), target :: diagram
      integer, intent(in) :: ieh_in, iet_in
      real(dp), intent(out) :: la1, p_a, p_r
      integer :: ieh, iet, t1, t2, iv, ipa, ina, ipb, inb, nu, order_c, qn(3)
      real(dp) :: tmax, pnu(Nph), pnu_fresh(Nph), lam, ft, lb1, lg1, Pq, tauc, tx, pt1, pt2
      logical :: shiftv(maxN), before
      type(vertex), pointer :: veh, vet
      ieh = ieh_in
      iet = iet_in
      order_c = diagram%order
      t1 = order_c - 1
      t2 = order_c
      if (ieh == t2 .and. iet == t1) then
         call swap_vertex(diagram, t1, t2)
      else
         if (ieh /= t1) then
            call swap_vertex(diagram, ieh, t1)
            if (iet == t1) iet = ieh
         end if
         if (iet /= t2) call swap_vertex(diagram, iet, t2)
      end if
      ieh = t1
      iet = t2
      veh => diagram%vertexList(ieh)
      vet => diagram%vertexList(iet)
      if (veh%plink /= iet .or. vet%plink /= ieh) stop 'ext_remove: relabel failed'
      tmax = diagram%vertexList(maxN)%tau
      qn = veh%i_q
      Pq = veh%Pq
      nu = veh%nu
      call mode_probs(veh%gkq, pnu)
      if (mod(n_extR, 10_8) == 0) then
         if (.not. vtmp_init) then
            call CreateVertex(Nph, dmc_band, nbnd, nsvd, vtmp)
            vtmp_init = .true.
         end if
         call copy_vtex(Nph, dmc_band, nbnd, nsvd, veh, vtmp)
         call cal_ek_int(vtmp%i_kin, vtmp%ekin, vtmp%ukin, vtmp%vkin)
         call cal_gkq_vtex_int(vtmp, vtmp%gkq)
         call mode_probs(vtmp%gkq, pnu_fresh)
         gauge_dpnu = max(gauge_dpnu, maxval(abs(pnu_fresh - pnu)))
      end if
      if (ext_mimic) then
         if (veh%link(3) == maxN) then
            tauc = 0.5d0 * tmax
         else
            tauc = min(diagram%vertexList(veh%link(3))%tau, 0.5d0 * tmax)
         end if
         tx = veh%tau
         call nexp(0.d0, tauc, minval(veh%ekin) + veh%wq(nu) - minval(veh%ekout), tx, pt1, .false.)
         if (veh%link(3) == iet) then
            tauc = max(veh%tau, 0.5d0 * tmax)
         else
            tauc = max(diagram%vertexList(vet%link(1))%tau, 0.5d0 * tmax)
         end if
         tx = vet%tau
         call nexp(tauc, tmax, minval(vet%ekin) - vet%wq(nu) - minval(vet%ekout), tx, pt2, .false.)
         ft = pt1 * pt2 * tmax
      else
         lam = max(veh%wq(nu), 0.d0)
         ft = ftau2(0.d0, tmax - vet%tau, tmax - veh%tau, lam)
      end if
      shiftv = .false.
      before = .true.
      iv = diagram%vertexList(1)%link(3)
      do while (iv /= maxN)
         if (iv == ieh) before = .false.
         if (before) shiftv(iv) = .true.
         if (iv == iet) before = .true.
         iv = diagram%vertexList(iv)%link(3)
      end do
      shiftv(ieh) = .false.
      shiftv(iet) = .false.
      ipa = veh%link(1)
      ina = veh%link(3)
      ipb = vet%link(1)
      inb = vet%link(3)
      if (ina == iet) then
         diagram%vertexList(ipa)%link(3) = inb
         diagram%vertexList(inb)%link(1) = ipa
      else
         diagram%vertexList(ipa)%link(3) = ina
         diagram%vertexList(ina)%link(1) = ipa
         diagram%vertexList(ipb)%link(3) = inb
         diagram%vertexList(inb)%link(1) = ipb
      end if
      diagram%order = order_c - 2
      diagram%nph_ext = diagram%nph_ext - 1
      call resync_shift(diagram, qn, shiftv)
      call eval_state(diagram, .false., lb1, la1, lg1)
      p_r = 0.5d0 * 2.d0 / real(order_c - 1, dp)
      p_a = 0.5d0 * (1.d0 / tmax) * Pq * pnu(nu) * ft
   end subroutine ext_remove_core

   subroutine any_move(diagram, accepted)
      type(fynman), intent(inout), target :: diagram
      logical, intent(out) :: accepted
      real(dp) :: t0
      t0 = now()
      if (urand() < 0.5d0) then
         call any_add(diagram, accepted)
      else
         call any_remove(diagram, accepted)
      end if
      t_any = t_any + (now() - t0)
   end subroutine any_move

   subroutine any_add(diagram, accepted)
      type(fynman), intent(inout), target :: diagram
      logical, intent(out) :: accepted
      integer :: ipa, ina, ipb, inb, ia, ib, u, nu, k, qn(3), nspan
      real(dp) :: tmax, tau1, tau2, Pqn, pnu(Nph), lam, L, x, c, ft, la0, la1, lb1, lg1, acc, p_r, p_a
      type(vertex), pointer :: va, vb, vu, vprev
      accepted = .false.
      if (diagram%order + 2 > maxOrder) return
      n_anyA = n_anyA + 1
      tmax = diagram%vertexList(maxN)%tau
      tau1 = urand() * tmax
      ipa = 1
      do while (diagram%vertexList(ipa)%link(3) /= maxN)
         if (diagram%vertexList(diagram%vertexList(ipa)%link(3))%tau >= tau1) exit
         ipa = diagram%vertexList(ipa)%link(3)
      end do
      ina = diagram%vertexList(ipa)%link(3)
      call eval_state(diagram, .false., lb1, la0, lg1)
      call full_snapshot(diagram)
      ia = diagram%order + 1
      ib = diagram%order + 2
      va => diagram%vertexList(ia)
      vb => diagram%vertexList(ib)
      va%tau = tau1
      call copy_in_from_out(va, diagram%vertexList(ipa))
      call sample_q_omp_int(diagram%seed, qn, Pqn, va)
      va%i_q = qn
      va%Pq = Pqn
      va%i_kout = modulo(va%i_kin + qn, nk_svd)
      call cal_ek_int(va%i_kout, va%ekout, va%ukout, va%vkout)
      va%Eoutmin = minval(va%ekout)
      call cal_wq_int(va%i_q, va%wq)
      call cal_gkq_vtex_int(va, va%gkq)
      call cal_gkq_full_vtex_int(va, va%gkq_full)
      call mode_probs(va%gkq, pnu)
      if (sum(pnu) <= 0.d0) then
         call full_restore(diagram)
         return
      end if
      x = urand()
      c = 0.d0
      nu = Nph
      do k = 1, Nph
         c = c + pnu(k)
         if (x < c) then
            nu = k
            exit
         end if
      end do
      va%nu = nu
      lam = max(va%wq(nu), 0.d0)
      L = tmax - tau1
      x = urand()
      if (lam * L < 1.d-8) then
         tau2 = tau1 + x * L
      else
         tau2 = tau1 - log(1.d0 - x * (1.d0 - exp(-lam * L))) / lam
      end if
      if (.not. (tau2 > tau1 .and. tau2 < tmax)) then
         call full_restore(diagram)
         return
      end if
      ft = ftau2(tau1, tau2, tmax, lam)
      ipb = ipa
      do while (diagram%vertexList(ipb)%link(3) /= maxN)
         if (diagram%vertexList(diagram%vertexList(ipb)%link(3))%tau >= tau2) exit
         ipb = diagram%vertexList(ipb)%link(3)
      end do
      inb = diagram%vertexList(ipb)%link(3)
      ! link a after ipa and b after ipb
      if (ipb == ipa) then
         va%link = (/ipa, ib, ib/)
         vb%link = (/ia, ia, ina/)
         diagram%vertexList(ipa)%link(3) = ia
         diagram%vertexList(ina)%link(1) = ib
      else
         va%link = (/ipa, ib, ina/)
         vb%link = (/ipb, ia, inb/)
         diagram%vertexList(ipa)%link(3) = ia
         diagram%vertexList(ina)%link(1) = ia
         diagram%vertexList(ipb)%link(3) = ib
         diagram%vertexList(inb)%link(1) = ib
      end if
      diagram%order = diagram%order + 2
      vb%tau = tau2
      vb%i_q = nk_svd - qn
      vb%Pq = Pqn
      vb%nu = nu
      vb%wq = va%wq
      ! shift the span by +q
      nspan = 1
      vprev => va
      u = va%link(3)
      do while (u /= ib)
         vu => diagram%vertexList(u)
         call copy_in_from_out(vu, vprev)
         vu%i_kout = modulo(vu%i_kin + vu%i_q, nk_svd)
         call cal_ek_int(vu%i_kout, vu%ekout, vu%ukout, vu%vkout)
         vu%Eoutmin = minval(vu%ekout)
         call cal_gkq_vtex_int(vu, vu%gkq)
         call cal_gkq_full_vtex_int(vu, vu%gkq_full)
         nspan = nspan + 1
         vprev => vu
         u = vu%link(3)
      end do
      call copy_in_from_out(vb, vprev)
      call copy_out_from_in(vb, diagram%vertexList(vb%link(3)))
      call cal_gkq_vtex_int(vb, vb%gkq)
      call cal_gkq_full_vtex_int(vb, vb%gkq_full)
      call eval_state(diagram, .false., lb1, la1, lg1)
      p_r = 0.5d0 * 2.d0 / real(diagram%order - 1, dp)
      p_a = 0.5d0 * (1.d0 / tmax) * Pqn * pnu(nu) * ft
      acc = exp(max(min(la1 - la0, 700.d0), -700.d0)) * p_r / p_a
      if (urand() < acc) then
         n_anyA_acc = n_anyA_acc + 1
         sum_spanA = sum_spanA + nspan
         accepted = .true.
      else
         call full_restore(diagram)
      end if
   end subroutine any_add

   subroutine any_remove(diagram, accepted)
      use multiphonon_update_matrix, only : swap_vertex
      type(fynman), intent(inout), target :: diagram
      logical, intent(out) :: accepted
      integer :: iv1, iv2, ia, ib, t1, t2, ipa, ina, ipb, inb, u, nu, nspan, order_c
      real(dp) :: tmax, pnu(Nph), pnu_fresh(Nph), lam, ft, la0, la1, lb1, lg1, acc, p_r, p_a, Pq
      type(vertex), pointer :: va, vb, vu, vprev, vlast
      accepted = .false.
      if (diagram%order < 3) return
      iv1 = 2 + int(urand() * (diagram%order - 1))
      if (iv1 > diagram%order) iv1 = diagram%order
      iv2 = diagram%vertexList(iv1)%link(2)
      if (iv2 == 1 .or. iv2 == maxN) return
      n_anyR = n_anyR + 1
      order_c = diagram%order
      call eval_state(diagram, .false., lb1, la0, lg1)
      call full_snapshot(diagram)
      ! relabel the pair into the last two slots while it is still linked
      if (diagram%vertexList(iv2)%tau < diagram%vertexList(iv1)%tau) then
         ia = iv2; ib = iv1
      else
         ia = iv1; ib = iv2
      end if
      t1 = order_c - 1
      t2 = order_c
      if (ia == t2 .and. ib == t1) then
         call swap_vertex(diagram, t1, t2)
      else
         if (ia /= t1) then
            call swap_vertex(diagram, ia, t1)
            if (ib == t1) ib = ia
         end if
         if (ib /= t2) call swap_vertex(diagram, ib, t2)
      end if
      ia = t1
      ib = t2
      va => diagram%vertexList(ia)
      vb => diagram%vertexList(ib)
      if (va%link(2) /= ib .or. vb%link(2) /= ia) stop 'any_remove: relabel failed'
      tmax = diagram%vertexList(maxN)%tau
      nu = va%nu
      Pq = va%Pq
      call mode_probs(va%gkq, pnu)
      ! diagnostic: the add move would recompute a's matrices with fresh k_out eigenvectors
      if (mod(n_anyR, 10_8) == 0) then
         if (.not. vtmp_init) then
            call CreateVertex(Nph, dmc_band, nbnd, nsvd, vtmp)
            vtmp_init = .true.
         end if
         call copy_vtex(Nph, dmc_band, nbnd, nsvd, va, vtmp)
         call cal_ek_int(vtmp%i_kout, vtmp%ekout, vtmp%ukout, vtmp%vkout)
         call cal_gkq_vtex_int(vtmp, vtmp%gkq)
         call mode_probs(vtmp%gkq, pnu_fresh)
         gauge_dpnu = max(gauge_dpnu, maxval(abs(pnu_fresh - pnu)))
      end if
      lam = max(va%wq(nu), 0.d0)
      ft = ftau2(va%tau, vb%tau, tmax, lam)
      ipa = va%link(1)
      ina = va%link(3)
      ipb = vb%link(1)
      inb = vb%link(3)
      nspan = 1
      if (ina == ib) then
         ! span 1: merge prev(a) -> next(b)
         diagram%vertexList(ipa)%link(3) = inb
         diagram%vertexList(inb)%link(1) = ipa
         if (inb == maxN) then
            if (ipa /= 1) then
               call copy_out_from_in(diagram%vertexList(ipa), diagram%vertexList(maxN))
               call cal_gkq_vtex_int(diagram%vertexList(ipa), diagram%vertexList(ipa)%gkq)
               call cal_gkq_full_vtex_int(diagram%vertexList(ipa), diagram%vertexList(ipa)%gkq_full)
            end if
         else
            call copy_in_from_out(diagram%vertexList(inb), diagram%vertexList(ipa))
            call cal_gkq_vtex_int(diagram%vertexList(inb), diagram%vertexList(inb)%gkq)
            call cal_gkq_full_vtex_int(diagram%vertexList(inb), diagram%vertexList(inb)%gkq_full)
         end if
      else
         ! interior ina .. ipb shifts by -q; the last interior vertex takes next(b)'s in
         vprev => diagram%vertexList(ipa)
         u = ina
         do while (u /= ib)
            vu => diagram%vertexList(u)
            call copy_in_from_out(vu, vprev)
            vu%i_kout = modulo(vu%i_kin + vu%i_q, nk_svd)
            if (vu%link(3) == ib) then
               call copy_out_from_in(vu, diagram%vertexList(inb))
            else
               call cal_ek_int(vu%i_kout, vu%ekout, vu%ukout, vu%vkout)
               vu%Eoutmin = minval(vu%ekout)
            end if
            call cal_gkq_vtex_int(vu, vu%gkq)
            call cal_gkq_full_vtex_int(vu, vu%gkq_full)
            nspan = nspan + 1
            vprev => vu
            u = vu%link(3)
         end do
         vlast => diagram%vertexList(ipb)
         diagram%vertexList(ipa)%link(3) = ina
         diagram%vertexList(ina)%link(1) = ipa
         vlast%link(3) = inb
         diagram%vertexList(inb)%link(1) = ipb
      end if
      diagram%order = order_c - 2
      call eval_state(diagram, .false., lb1, la1, lg1)
      p_r = 0.5d0 * 2.d0 / real(order_c - 1, dp)
      p_a = 0.5d0 * (1.d0 / tmax) * Pq * pnu(nu) * ft
      acc = exp(max(min(la1 - la0, 700.d0), -700.d0)) * p_a / p_r
      if (urand() < acc) then
         n_anyR_acc = n_anyR_acc + 1
         sum_spanR = sum_spanR + nspan
         accepted = .true.
      else
         call full_restore(diagram)
      end if
   end subroutine any_remove

   ! Candidate slow variables for the trace: squared minimum-image grid distance of
   ! the head (periodic) segment momentum from Gamma, and the number of internal lines.
   subroutine slow_vars(diagram, kh2, nint)
      type(fynman), intent(inout), target :: diagram
      integer, intent(out) :: kh2, nint
      integer :: k(3), iv
      k = modulo(diagram%vertexList(1)%i_kout, nk_svd)
      where (k > nk_svd / 2) k = k - nk_svd
      kh2 = sum(k * k)
      nint = 0
      iv = diagram%vertexList(1)%link(3)
      do while (iv /= maxN)
         if (diagram%vertexList(iv)%link(2) /= 1 .and. diagram%vertexList(iv)%link(2) /= maxN) nint = nint + 1
         iv = diagram%vertexList(iv)%link(3)
      end do
      nint = nint / 2
   end subroutine slow_vars

   ! Native-chain measurement with change_q: trace plus self-checks.
   subroutine chq_measure(diagram, stat)
      type(fynman), intent(inout), target :: diagram
      type(diagmc_stat), intent(in) :: stat
      complex(dp) :: dv, dh, dw
      real(dp) :: wf2, lg2, nat_inc, nat_sign, oA, err
      integer :: ns2, kh2, nint
      n_meas = n_meas + 1
      nat_inc = real(stat%Etrue, dp) - e_before
      nat_sign = real(stat%gtrue, dp) - g_before
      call slow_vars(diagram, kh2, nint)
      write(bc_unit, '(2ES24.15,I8,I6,I8,I6)') nat_inc, nat_sign, diagram%order, diagram%nph_ext, kh2, nint
      if (mod(n_meas, 10_8) == 0) then
         call k_check(diagram)
         call g_check(diagram)
         if (diagram%order > 1) then
            call evaluate(diagram, .false., dv, dh, dw, wf2, ns2, lg2)
            oA = real(dh + dw + (wf2 - (diagram%order - 1)) * dv, dp) / abs(real(dv, dp))
            err = abs(oA - nat_inc) / max(1.d0, abs(nat_inc))
            max_o_rel = max(max_o_rel, err)
            if (err > tol) n_bad_o = n_bad_o + 1
            if (abs(sign(1.d0, real(dv, dp)) - nat_sign) > tol) n_bad_sign = n_bad_sign + 1
         end if
      end if
      if (mod(n_meas, 500_8) == 0) call chq_summary()
   end subroutine chq_measure

   subroutine chq_summary()
      integer :: u
      open(newunit=u, file='chq_summary.dat', status='replace', action='write')
      write(u, '(A,F10.4)') 'p_chq ', p_chq
      write(u, '(A,2I14)') 'mv_until_steps_seen ', mv_until, n_mvsteps
      write(u, '(A,3I14)') 'chq_tried_accepted_external ', n_chq, n_chq_acc, n_chq_ext
      write(u, '(A,L2,2I14)') 'ext_on_tried_accepted ', chq_ext_on, n_chq_ext_try, n_chq_ext_acc
      write(u, '(A,I12,4ES12.4)') 'roundtrip_n_dla_dpnu_dg_du ', n_rt, rt_la, rt_pnu, rt_g, rt_u
      write(u, '(A,F10.4,4I14)') 'any_p_add_tried_acc_remove_tried_acc ', p_any, n_anyA, n_anyA_acc, n_anyR, n_anyR_acc
      write(u, '(A,2F12.3,F12.3)') 'any_mean_span_add_rem_acc_time ', sum_spanA / max(1_8, n_anyA_acc), &
           sum_spanR / max(1_8, n_anyR_acc), t_any
      write(u, '(A,2I12)') 'gauge_checks_bad ', n_gchecks, n_gbad
      write(u, '(A,F10.4,4I14,F12.3)') 'ext_p_add_tried_acc_remove_tried_acc_time ', p_ext, n_extA, n_extA_acc, &
           n_extR, n_extR_acc, t_ext
      write(u, '(A,I12,3ES12.4)') 'ext_roundtrip_n_dla_dpa_dpr ', n_rt_ext, rt_ext_la, rt_ext_pa, rt_ext_pr
      write(u, '(A,L2,ES12.4)') 'pnu_fro_max_dpnu_stored_vs_fresh ', pnu_fro, gauge_dpnu
      write(u, '(A,2F12.3)') 'mean_span_tried_accepted ', sum_span_try / max(1_8, n_chq), &
           sum_span_acc / max(1_8, n_chq_acc)
      write(u, '(A,F12.3)') 'time_chq ', t_chq
      write(u, '(A,I12)') 'n_meas ', n_meas
      write(u, '(A,2I12)') 'k_checks_bad ', n_kchecks, n_kbad
      write(u, '(A,I12)') 'n_O_mismatch ', n_bad_o
      write(u, '(A,ES12.4)') 'max_O_rel ', max_o_rel
      write(u, '(A,I12)') 'n_sign_mismatch ', n_bad_sign
      close(u)
      flush(bc_unit)
   end subroutine chq_summary

   ! ---- hooks around update_drive -------------------------------------------------
   subroutine bchain_pre_update(diagram, stat)
      type(fynman), intent(inout), target :: diagram
      type(diagmc_stat), intent(in) :: stat
      logical :: moved
      if (.not. bc_checked) call bc_setup(diagram)
      if (.not. (bc_on .or. mv_on)) return
      if (.not. bc_on) then
         n_mvsteps = n_mvsteps + 1
         if (n_mvsteps > mv_until) return
         if (chq_on) then
            if (urand() < p_chq) call change_q(diagram, moved)
         end if
         if (any_on) then
            if (urand() < p_any) call any_move(diagram, moved)
         end if
         if (ext_on) then
            if (urand() < p_ext) call ext_move(diagram, moved)
         end if
         return
      end if
      if (.not. bc_init) then
         call take_snapshot(diagram)
         call set_current(diagram)
         bc_init = .true.
      end if
      ! B-native add/remove (exact MH for pi_B); native types 1,2 must be off (PA=0).
      if (bc_addrem .and. mod(n_steps, int(blk, 8)) == 0) then
         if (urand() < p_addrem * blk) then
            t_mark = now()
            if (urand() < 0.5d0) then
               call add_B(diagram)
            else
               call remove_B(diagram)
            end if
            t_bmove = t_bmove + (now() - t_mark)
         end if
      end if
      ! change_q is one more pi_A-reversible piece of the composite proposal
      if (chq_on) then
         if (urand() < p_chq) then
            call change_q(diagram, moved)
            if (moved) blk_changed = .true.
         end if
      end if
      if (any_on) then
         if (urand() < p_any) then
            call any_move(diagram, moved)
            if (moved) blk_changed = .true.
         end if
      end if
      if (ext_on) then
         if (urand() < p_ext) then
            call ext_move(diagram, moved)
            if (moved) blk_changed = .true.
         end if
      end if
      acc0 = sum(stat%accept)
      accv0 = stat%accept(1:7)
      t_mark = now()
   end subroutine bchain_pre_update

   subroutine bchain_post_update(diagram, stat)
      type(fynman), intent(inout), target :: diagram
      type(diagmc_stat), intent(in) :: stat
      real(dp) :: r_new, t1, lb1, la1, lg1
      integer :: ty, k
      if (.not. bc_on) return
      t1 = now()
      t_native = t_native + (t1 - t_mark)
      n_steps = n_steps + 1
      if (sum(stat%accept) /= acc0) then
         n_prop = n_prop + 1
         blk_changed = .true.
         ty = 0
         do k = 1, 7
            if (stat%accept(k) /= accv0(k)) ty = k
         end do
         if (ty > 0) prop_t(ty) = prop_t(ty) + 1
      end if
      if (mod(n_steps, int(blk, 8)) /= 0) return
      ! checkpoint: second-stage correction of the composite proposal K_A^blk
      n_blocks = n_blocks + 1
      if (blk_changed) then
         n_blocks_changed = n_blocks_changed + 1
         t1 = now()
         call eval_state(diagram, .true., lb1, la1, lg1)
         r_new = exp(max(min(lb1 - lg1 - la1, 700.d0), -700.d0))
         t_ratio = t_ratio + (now() - t1)
         t1 = now()
         if (urand() < min(1.d0, r_new / r_snap)) then
            n_blocks_acc = n_blocks_acc + 1
            n_acc2 = n_acc2 + 1
            r_snap = r_new
            lwb_cur = lb1; lwa_cur = la1; lgrp_cur = lg1
            call take_snapshot(diagram)
         else
            call restore_snapshot(diagram)
         end if
         t_snap = t_snap + (now() - t1)
         blk_changed = .false.
      end if
      ! Group-invariant Gibbs move; W_B and |G| are group invariants, only w_A changes.
      if (bc_group .and. mod(n_steps, int(refresh_every, 8)) < blk) then
         t1 = now()
         call refresh_modes(diagram)
         call eval_state(diagram, .false., lb1, la1, lg1)
         lwa_cur = la1
         r_snap = exp(max(min(lwb_cur - lgrp_cur - lwa_cur, 700.d0), -700.d0))
         call take_snapshot(diagram)
         t_refresh = t_refresh + (now() - t1)
      end if
      if (mod(n_steps, 1000_8) < blk) call check_cache(diagram)
   end subroutine bchain_post_update

   ! ---- measurement -------------------------------------------------------------
   subroutine bchain_pre_measure(stat)
      type(diagmc_stat), intent(in) :: stat
      if (.not. (bc_on .or. mv_on)) return
      e_before = real(stat%Etrue, dp)
      g_before = real(stat%gtrue, dp)
      t_meas_mark = now()
   end subroutine bchain_pre_measure

   subroutine bchain_post_measure(diagram, stat)
      type(fynman), intent(inout), target :: diagram
      type(diagmc_stat), intent(in) :: stat
      complex(dp) :: fv, fh, fw, dv, dh, dw
      real(dp) :: wf, lg, wf2, lg2, fO, f1, nat_inc, nat_sign, oA, err
      integer :: ns, ns2
      if (.not. (bc_on .or. mv_on)) return
      if (.not. bc_on) then
         call chq_measure(diagram, stat)
         return
      end if
      t_natmeas = t_natmeas + (now() - t_meas_mark)
      t_meas_mark = now()
      n_meas = n_meas + 1
      nat_inc = real(stat%Etrue, dp) - e_before
      nat_sign = real(stat%gtrue, dp) - g_before
      if (diagram%order <= 1) then
         write(bc_unit, '(2ES24.15,I6,ES16.7,2ES24.15,I8)') nat_inc, nat_sign, 0, 0.d0, nat_inc, nat_sign, diagram%order
         return
      end if
      call evaluate(diagram, .false., dv, dh, dw, wf2, ns2, lg2)
      oA = real(dh + dw + (wf2 - (diagram%order - 1)) * dv, dp) / abs(real(dv, dp))
      err = abs(oA - nat_inc) / max(1.d0, abs(nat_inc))
      max_o_rel = max(max_o_rel, err)
      if (err > tol) n_bad_o = n_bad_o + 1
      if (abs(sign(1.d0, real(dv, dp)) - nat_sign) > tol) n_bad_sign = n_bad_sign + 1
      call evaluate(diagram, .true., fv, fh, fw, wf, ns, lg)
      fO = real(fh + fw + (wf - (diagram%order - 1)) * fv, dp) / abs(real(fv, dp))
      f1 = sign(1.d0, real(fv, dp))
      if (ns >= 1 .and. mod(n_meas, 10_8) == 0) call enum_check(diagram)
      sum_seff = sum_seff + ns
      sum_loggroup = sum_loggroup + lg
      write(bc_unit, '(2ES24.15,I6,ES16.7,2ES24.15,I8)') fO, f1, ns, lg, nat_inc, nat_sign, diagram%order
      t_meas = t_meas + (now() - t_meas_mark)
      if (mod(n_meas, 500_8) == 0) then
         call bc_summary()
         if (mv_on) call chq_summary()
      end if
   end subroutine bchain_post_measure

   ! Fold vs explicit enumeration on the first one or two lines of S_eff: the
   ! forced live-mode fold of those lines must equal the sum over all their
   ! live-mode assignments (others fixed), for V and for H + W.
   subroutine enum_check(diagram)
      type(fynman), intent(inout), target :: diagram
      complex(dp) :: fv, fh, fw
      real(dp) :: wf
      integer :: ids(maxN), posof(maxN), la(maxN), lb(maxN), inS(maxN), nv, iv, j, l, l2, nl, ns, partner
      integer :: gl(2), ng, nu1, nu2, keep1, keep2, cnt
      complex(dp) :: sv, sh, sw, dv, dh, dw
      real(dp) :: wf1, lg, err
      integer :: nsx
      logical :: ok
      type(vertex), pointer :: v1a, v1b, v2a, v2b
      nv = 0
      iv = diagram%vertexList(1)%link(3)
      do while (iv /= maxN)
         nv = nv + 1
         ids(nv) = iv
         iv = diagram%vertexList(iv)%link(3)
      end do
      posof = 0
      do j = 1, nv
         posof(ids(j)) = j
      end do
      nl = 0
      do j = 1, nv
         partner = diagram%vertexList(ids(j))%link(2)
         if (partner == 1 .or. partner == maxN) cycle
         if (posof(partner) > j) then
            nl = nl + 1
            la(nl) = j
            lb(nl) = posof(partner)
         end if
      end do
      ns = 0
      call select_S(nl, la, lb, inS, ns)
      ng = 0
      do l2 = 1, ns
         l = inS(l2)
         if (live(diagram%vertexList(ids(la(l)))%wq(diagram%vertexList(ids(la(l)))%nu))) then
            ng = ng + 1
            if (ng <= 2) gl(ng) = l
         end if
      end do
      if (ng < 1) return
      ng = min(ng, 2)
      call evaluate(diagram, .true., fv, fh, fw, wf, nsx, lg, forced_nl=ng, forced=gl(1:ng))
      v1a => diagram%vertexList(ids(la(gl(1))))
      v1b => diagram%vertexList(ids(lb(gl(1))))
      keep1 = v1a%nu
      if (ng == 2) then
         v2a => diagram%vertexList(ids(la(gl(2))))
         v2b => diagram%vertexList(ids(lb(gl(2))))
         keep2 = v2a%nu
      end if
      sv = 0; sh = 0; sw = 0
      cnt = 0
      ! Sum of per-member evaluations (the same lines grouped, current mode only).
      ! The fold normalizes each vertex by its current-mode max|g|, which changes
      ! with the member, so every member is rescaled to the reference normalization.
      do nu1 = 1, Nph
         if (.not. live(v1a%wq(nu1))) cycle
         do nu2 = 1, merge(Nph, 1, ng == 2)
            if (ng == 2) then
               if (.not. live(v2a%wq(nu2))) cycle
            end if
            v1a%nu = nu1; v1b%nu = nu1
            if (ng == 2) then
               v2a%nu = nu2; v2b%nu = nu2
            end if
            call evaluate(diagram, .false., dv, dh, dw, wf1, nsx, lg, forced_nl=ng, forced=gl(1:ng))
            call rescale(dv, dh, dw, wf1)
            sv = sv + dv; sh = sh + dh; sw = sw + dw
            cnt = cnt + 1
         end do
      end do
      v1a%nu = keep1; v1b%nu = keep1
      if (ng == 2) then
         v2a%nu = keep2; v2b%nu = keep2
      end if
      n_enum = n_enum + 1
      err = abs(sv - fv) / abs(fv)
      err = max(err, abs((sh + sw) - (fh + fw)) / max(abs(fh + fw), abs(fv)))
      max_enum_rel = max(max_enum_rel, err)
      if (err > 1.d-7) n_bad_enum = n_bad_enum + 1
   contains
      ! Bring a member evaluation to the normalization of the reference state:
      ! multiply by prod over the grouped vertices of max|g(member nu)| / max|g(keep)|.
      ! The member's own grouped-line phonon energy enters W via the value, and the
      ! fixed term wf1 of the member is subtracted so only grouped-line terms remain.
      subroutine rescale(xv, xh, xw, wf1x)
         complex(dp), intent(inout) :: xv, xh, xw
         real(dp), intent(in) :: wf1x
         real(dp) :: fac
         fac = maxval(abs(v1a%gkq(:, :, v1a%nu))) / maxval(abs(v1a%gkq(:, :, keep1))) &
             * maxval(abs(v1b%gkq(:, :, v1b%nu))) / maxval(abs(v1b%gkq(:, :, keep1)))
         if (ng == 2) then
            fac = fac * maxval(abs(v2a%gkq(:, :, v2a%nu))) / maxval(abs(v2a%gkq(:, :, keep2))) &
                      * maxval(abs(v2b%gkq(:, :, v2b%nu))) / maxval(abs(v2b%gkq(:, :, keep2)))
         end if
         xv = xv * fac
         xh = xh * fac
         xw = xw * fac
         if (wf1x /= wf1x) xw = xw   ! wf1x (fixed term) equals the reference wf by construction
      end subroutine rescale
   end subroutine enum_check

   subroutine bc_summary()
      integer :: u
      open(newunit=u, file='bchain_summary.dat', status='replace', action='write')
      write(u, '(A,L2)') 'grouping ', bc_group
      write(u, '(A,I14)') 'n_steps ', n_steps
      write(u, '(A,I14)') 'n_native_accepted ', n_prop
      write(u, '(A,I14)') 'n_second_stage_accepted ', n_acc2
      write(u, '(A,I6,3I14)') 'block_size_blocks_changed_accepted ', blk, n_blocks, n_blocks_changed, n_blocks_acc
      write(u, '(A,7I12)') 'native_accepted_by_type ', prop_t
      write(u, '(A,7I12)') 'second_stage_accepted_by_type ', acc2_t
      write(u, '(A,4I14)') 'B_add_tried_acc_B_remove_tried_acc ', n_addB, n_addB_acc, n_remB, n_remB_acc
      write(u, '(A,I12)') 'n_meas ', n_meas
      write(u, '(A,I12)') 'n_O_mismatch ', n_bad_o
      write(u, '(A,ES12.4)') 'max_O_rel ', max_o_rel
      write(u, '(A,I12)') 'n_sign_mismatch ', n_bad_sign
      write(u, '(A,I12)') 'n_enum_checks ', n_enum
      write(u, '(A,I12)') 'n_enum_mismatch ', n_bad_enum
      write(u, '(A,ES12.4)') 'max_enum_rel ', max_enum_rel
      write(u, '(A,F12.4)') 'mean_Seff ', sum_seff / max(1_8, n_meas)
      write(u, '(A,F12.4)') 'mean_log_group ', sum_loggroup / max(1_8, n_meas)
      write(u, '(A,I12)') 'n_cache_checks ', n_cache_checks
      write(u, '(A,I12)') 'n_cache_mismatch ', n_cache_bad
      write(u, '(A,ES12.4)') 'max_cache_err ', max_cache_err
      write(u, '(A,7F12.3)') 'time_native_ratio_snap_bmove_refresh_meas_natmeas ', &
           t_native, t_ratio, t_snap, t_bmove, t_refresh, t_meas, t_natmeas
      close(u)
      flush(bc_unit)
   end subroutine bc_summary

end module bchain_mod
