! Future B research: window Rao-Blackwellization of the Luo EZ energy estimator,
! and the exact sign gain of the tiled grouped measure, measured on native chains.
!
! Opt-in with FUTUREB_RB=1; otherwise every entry point returns immediately and
! the native chain, RNG stream and observables are untouched. Measurement-only.
! Any dmc_band, electron or hole, DMC_Method 0 (gkq == gkq_full).
!
! Native weight (update_swap acceptance): D(C) = P(C) * w(C) with
!   w(C) = Re tr[ G_last g_nv ... g_1 G_1 ],   g_v = gkq_v / max|gkq_v|,
!   G_s = diag exp(-(E_b(k_s) - Emin_s) dtau_s)  (native propagator, clamped at e^-100),
!   P(C) = prod_s exp(-Emin_s dtau_s) * prod_lines Dph * prod_v max|gkq_v|  (> 0).
! Native EZ increment (measure_EZ_wfn): N(C)/|w(C)| with
!   N = Re sum_s tr[chain with G_s -> HG_s] + (Wdt - (order-1)) w,
!   HG_s = diag exp(-(E_b-Emin_s) dtau_s) E_b dtau_s.
!
! With prefix products X_s (up to segment s), their H-insertion sums Y_s, and
! suffixes Z_s, V_s (M = Z_s X_s), a window on vertices p..p+3 factorizes as
!   M = R W L,  L = X_p, Lam = Y_p, R = Z_{p+4} G_{p+4}, Rho = V_{p+4} G_{p+4} + Z_{p+4} HG_{p+4},
!   N_H(r) = Re[tr(R W_r Lam) + tr(Rho W_r L) + tr(R W^H_r L)],
! so each of the three local pairings needs only its own 4 vertices and 3 segments.
!
! Per measurement it writes (rb_trace.dat):
!   raw_num raw_den  - native estimator N/|w|, sign(w)
!   rb_num rb_den    - average over all window positions of the fibre conditional
!                      expectation sum_r P_r N_r / sum_r P_r|w_r| (ineligible: raw)
!   n_elig n_win order
!   rho_single       - average over positions of |sum_r P_r w_r| / sum_r P_r|w_r|
!   rho_tiled        - |sum_{C' in G} D| / sum_{C' in G} |D| over the tiled group G(C):
!                      all 3^m re-pairings of the first m = min(#eligible, 8) windows among
!                      positions 1,5,9,..
!   m_tiled m_used   - number of eligible tiled windows, and the number grouped
! Under the native measure E[rho_tiled] = Z_B / Z_A exactly (Lemma 3 partition), so
! <s>_B = <s>_A / E[rho_tiled] is the sign of the grouped measure B.
module rb_window_mod
   use DiagMC
   use pert_param, only : nsvd, rb_method => DMC_Method
   implicit none
   private
   public :: rb_pre, rb_post

   logical, save :: rb_checked = .false., rb_on = .false.
   integer, save :: rb_unit = 0, nb = 0
   real(dp), save :: e_before = 0.d0, g_before = 0.d0
   integer(8), save :: n_meas = 0, n_sign_bad = 0, n_o_bad = 0, n_r0_bad = 0, n_env_bad = 0
   integer(8), save :: n_elig = 0, n_win = 0, n_zero_fibre = 0, n_tiled_skip = 0, n_tiled_bad = 0
   real(dp), save :: max_o_rel = 0.d0, max_r0_rel = 0.d0, max_env_rel = 0.d0, max_tiled_rel = 0.d0
   type(vertex), save :: vs
   integer, parameter :: MAXT = 250, MAXM_ENUM = 8
   complex(dp), allocatable, save :: X(:,:,:), Y(:,:,:), Z(:,:,:), V(:,:,:), gn(:,:,:), TW(:,:,:,:)
   real(dp), allocatable, save :: Gd(:,:), HGd(:,:), TlogP(:,:)
   integer, allocatable, save :: tpos(:), tr0(:)
   real(dp), parameter :: tol = 1.d-8
   ! Phonon-mode grouping diagnostics (FUTUREB_MODES=1, requires FUTUREB_RB=1).
   logical, save :: modes_on = .false., modes_k4 = .false., modes_dump = .false., mode_rb_on = .false.
   integer, save :: mrb_unit = 0
   logical, save :: windows_on = .true., mode_svd_on = .false.
   integer, save :: svd_unit = 0
   integer(8), save :: n_mrb_bad = 0
   real(dp), save :: max_mrb_rel = 0.d0
   integer, save :: dump_unit = 0, dump_every = 10
   integer, save :: mode_unit = 0, k4_every = 20
   ! Momentum grouping diagnostic (FUTUREB_QGROUP=1, every FUTUREB_Q_EVERY-th measurement).
   logical, save :: qgroup_on = .false.
   integer, save :: q_unit = 0, q_every = 50
   integer(8), save :: n_q_bad = 0
   real(dp), save :: max_q_rel = 0.d0
   integer(8), save :: n_mode_bad = 0
   real(dp), save :: max_mode_rel = 0.d0

contains

   subroutine rb_setup()
      character(len=16) :: env
      integer :: st
      rb_checked = .true.
      call get_environment_variable('FUTUREB_RB', env, status=st)
      rb_on = (st == 0 .and. trim(adjustl(env)) == '1')
      if (.not. rb_on) return
      if (rb_method /= 0) stop 'FUTUREB_RB: DMC_Method 0 only'
      nb = dmc_band
      allocate(X(nb,nb,0:maxN), Y(nb,nb,0:maxN), Z(nb,nb,0:maxN), V(nb,nb,0:maxN), gn(nb,nb,maxN))
      allocate(Gd(nb,maxN), HGd(nb,maxN), TW(nb,nb,3,MAXT), TlogP(3,MAXT), tpos(MAXT), tr0(MAXT))
      call CreateVertex(Nph, dmc_band, nbnd, nsvd, vs)
      open(newunit=rb_unit, file='rb_trace.dat', status='replace', action='write')
      write(rb_unit, '(A)') '# raw_num raw_den rb_num rb_den n_elig n_win order rho_single rho_tiled m_tiled m_used'
      call get_environment_variable('FUTUREB_MODES', env, status=st)
      modes_on = (st == 0 .and. trim(adjustl(env)) == '1')
      call get_environment_variable('FUTUREB_RB_WINDOWS', env, status=st)
      if (st == 0 .and. trim(adjustl(env)) == '0') windows_on = .false.
      call get_environment_variable('FUTUREB_MODE_SVD', env, status=st)
      mode_svd_on = (st == 0 .and. trim(adjustl(env)) == '1')
      if (mode_svd_on) then
         open(newunit=svd_unit, file='mode_svd_trace.dat', status='replace', action='write')
         write(svd_unit, '(A)') '# n_lines mean_gamma prod_gamma mean_rho prod_rho mean_rank p1 p2 joint_gamma prod_gamma_pair joint_rho prod_rho_pair'
      end if
      call get_environment_variable('FUTUREB_MODE_RB', env, status=st)
      mode_rb_on = (st == 0 .and. trim(adjustl(env)) == '1')
      if (mode_rb_on) then
         open(newunit=mrb_unit, file='mode_rb_trace.dat', status='replace', action='write')
         write(mrb_unit, '(A)') '# raw_num raw_den mrb_num mrb_den n_lines   (mode Rao-Blackwell, mean over lines)'
      end if
      call get_environment_variable('FUTUREB_MODES_DUMP', env, status=st)
      modes_dump = (st == 0 .and. trim(adjustl(env)) == '1')
      if (modes_dump) then
         open(newunit=dump_unit, file='lines_dump.dat', status='replace', action='write')
         write(dump_unit, '(A)') '# per sampled measurement: "M n_meas n_lines sign" then n_lines rows "la lb rho_line"'
      end if
      call get_environment_variable('FUTUREB_MODES_K4', env, status=st)
      modes_k4 = modes_on .and. (st == 0 .and. trim(adjustl(env)) == '1')
      call get_environment_variable('FUTUREB_MODES_K4_EVERY', env, status=st)
      if (st == 0 .and. len_trim(env) > 0) read(env, *) k4_every
      call get_environment_variable('FUTUREB_QGROUP', env, status=st)
      qgroup_on = (st == 0 .and. trim(adjustl(env)) == '1')
      call get_environment_variable('FUTUREB_Q_EVERY', env, status=st)
      if (st == 0 .and. len_trim(env) > 0) read(env, *) q_every
      if (qgroup_on) then
         open(newunit=q_unit, file='q_trace.dat', status='replace', action='write')
         write(q_unit, '(A)') '# n_meas n_lines line_rank span rho_q rho_mode_same_line'
      end if
      if (modes_on) then
         open(newunit=mode_unit, file='mode_trace.dat', status='replace', action='write')
         write(mode_unit, '(A)') '# n_lines rho_line_mean rho_line_min rho_K2 rho_K3 frac_lines_lt_0.99 rho_l1 rho_l2 rho_l3 rho_K4 prod_all_lines n_laminar prod_laminar pair_rel rho_pi rho_pj rho_pair'
      end if
   end subroutine rb_setup

   subroutine rb_pre(stat)
      type(diagmc_stat), intent(in) :: stat
      if (.not. rb_checked) call rb_setup()
      if (.not. rb_on) return
      e_before = real(stat%Etrue, dp)
      g_before = real(stat%gtrue, dp)
   end subroutine rb_pre

   complex(dp) function ctr(A)
      complex(dp), intent(in) :: A(:,:)
      integer :: i
      ctr = (0.d0, 0.d0)
      do i = 1, size(A, 1)
         ctr = ctr + A(i, i)
      end do
   end function ctr

   subroutine eye(A)
      complex(dp), intent(out) :: A(:,:)
      integer :: i
      A = (0.d0, 0.d0)
      do i = 1, size(A, 1)
         A(i, i) = (1.d0, 0.d0)
      end do
   end subroutine eye

   subroutine rowscale(d, A)
      real(dp), intent(in) :: d(:)
      complex(dp), intent(inout) :: A(:,:)
      integer :: i
      do i = 1, size(A, 1)
         A(i, :) = d(i) * A(i, :)
      end do
   end subroutine rowscale

   subroutine colscale(A, d)
      complex(dp), intent(inout) :: A(:,:)
      real(dp), intent(in) :: d(:)
      integer :: j
      do j = 1, size(A, 2)
         A(:, j) = A(:, j) * d(j)
      end do
   end subroutine colscale

   ! Native propagator / Hpropagator for one segment (band-diagonal).
   subroutine seg_prop(E, dt, G, HG, emin)
      real(dp), intent(in) :: E(:), dt
      real(dp), intent(out) :: G(:), HG(:), emin
      integer :: b
      real(dp) :: el
      emin = minval(E)
      do b = 1, size(E)
         el = E(b) - emin
         G(b) = exp(-min(dt * el, 100.0_dp))
         HG(b) = exp(-dt * el) * (el + emin) * dt
      end do
   end subroutine seg_prop

   subroutine rb_post(diagram, stat)
      type(fynman), intent(inout), target :: diagram
      type(diagmc_stat), intent(in) :: stat
      integer :: ids(maxN), nv, nseg, iv, s, p, j, jj, r, loc0, jb, partner, nwin, nel, r0, ntile, m
      integer :: kin(3), ktmp(3), iq(0:3,3), nuv(0:3), idx(MAXM_ENUM), k, ncomb
      integer, parameter :: pa(3,2) = reshape((/0, 0, 0, 1, 3, 2/), (/3, 2/))  ! slot a (open, close)
      integer, parameter :: pb(3,2) = reshape((/2, 1, 1, 3, 2, 3/), (/3, 2/))  ! slot b (open, close)
      real(dp) :: etrue_inc, g_inc, sgn, w, NH, Nc, Wdt, emin, raw_num, raw_den, rbn, rbd, rhos, rho_t
      real(dp) :: tv(0:3), wqa, wqb, ta, tb, gm, Ssum, lmax, lp, d, sumD, sumA, err
      real(dp) :: wr(3), NHr(3), Nr(3), OW(3), logPr(3), Pr(3)
      real(dp) :: ek(nb), Ej(nb,0:2), Gw(nb), HGw(nb)
      complex(dp) :: A(nb,nb), B(nb,nb), C(nb,nb), L(nb,nb), Lam(nb,nb), Renv(nb,nb), Rho(nb,nb)
      complex(dp) :: T(nb,nb), TH(nb,nb), gr(nb,nb,0:3), grU(nb,nb,0:3), Twin(nb,nb,3)
      complex(dp) :: Cb(nb,nb,MAXM_ENUM), Lf(nb,nb), Rl(nb,nb)
      complex(dp) :: gbuf(nb,nb,Nph), uk(nbnd,nbnd), ukprev(nbnd,nbnd)
      integer :: qa(3), qb(3), nua, nub
      type(vertex), pointer :: vi, vp

      if (.not. rb_on) return
      etrue_inc = real(stat%Etrue, dp) - e_before
      g_inc = real(stat%gtrue, dp) - g_before
      n_meas = n_meas + 1

      nv = 0
      iv = diagram%vertexList(1)%link(3)
      do while (iv /= maxN)
         nv = nv + 1
         ids(nv) = iv
         iv = diagram%vertexList(iv)%link(3)
      end do
      if (nv /= diagram%order - 1) stop 'rb_window: internal vertex count'
      if (nv == 0) then
         write(rb_unit, '(4ES24.15,3I8,2ES24.15,2I6)') etrue_inc, g_inc, etrue_inc, g_inc, 0, 0, &
              diagram%order, 1.d0, 1.d0, 0, 0
         if (mode_rb_on) write(mrb_unit, '(4ES24.15,I6)') etrue_inc, g_inc, etrue_inc, g_inc, 0
         if (modes_on) write(mode_unit, '(I6,10ES20.11,I6,ES20.11,I4,3ES20.11)') 0, 1.d0, 1.d0, 1.d0, 1.d0, 0.d0, &
              1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 0, 1.d0, -1, 1.d0, 1.d0, 1.d0
         return
      end if
      nseg = nv + 1

      ! Segment propagators and normalized vertex matrices, as in measure_EZ_wfn.
      do s = 1, nseg
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
         call seg_prop(vi%ekout(1:nb), vp%tau - vi%tau, Gd(:, s), HGd(:, s), emin)
      end do
      do j = 1, nv
         vi => diagram%vertexList(ids(j))
         gn(:, :, j) = vi%gkq(:, :, vi%nu) / maxval(abs(vi%gkq(:, :, vi%nu)))
      end do

      ! Prefix products X_s and H-insertion sums Y_s.
      call eye(X(:, :, 0))
      Y(:, :, 0) = (0.d0, 0.d0)
      do s = 1, nseg
         A = X(:, :, s - 1)
         B = Y(:, :, s - 1)
         if (s >= 2) then
            A = matmul(gn(:, :, s - 1), A)
            B = matmul(gn(:, :, s - 1), B)
         end if
         X(:, :, s) = A
         call rowscale(Gd(:, s), X(:, :, s))
         C = A
         call rowscale(HGd(:, s), C)
         call rowscale(Gd(:, s), B)
         Y(:, :, s) = B + C
      end do
      ! Suffix products: M = Z_s X_s, V_s = insertions in segments > s.
      call eye(Z(:, :, nseg))
      V(:, :, nseg) = (0.d0, 0.d0)
      do s = nseg, 1, -1
         A = Z(:, :, s)
         call colscale(A, Gd(:, s))
         B = V(:, :, s)
         call colscale(B, Gd(:, s))
         C = Z(:, :, s)
         call colscale(C, HGd(:, s))
         B = B + C
         if (s >= 2) then
            A = matmul(A, gn(:, :, s - 1))
            B = matmul(B, gn(:, :, s - 1))
         end if
         Z(:, :, s - 1) = A
         V(:, :, s - 1) = B
      end do

      w = real(ctr(X(:, :, nseg)), dp)
      NH = real(ctr(Y(:, :, nseg)), dp)
      err = max(abs(real(ctr(Z(:, :, 0)), dp) - w) / abs(w), abs(real(ctr(V(:, :, 0)), dp) - NH) / max(abs(NH), abs(w)))
      max_env_rel = max(max_env_rel, err)
      if (err > tol) n_env_bad = n_env_bad + 1

      Wdt = 0.d0
      do j = 1, nv
         vi => diagram%vertexList(ids(j))
         partner = vi%link(2)
         if (partner /= 1 .and. partner /= maxN) then
            Wdt = Wdt + 0.5d0 * vi%wq(vi%nu) * abs(vi%tau - diagram%vertexList(partner)%tau)
         else
            Wdt = Wdt + vi%wq(vi%nu) * abs(vi%tau - diagram%vertexList(partner)%tau)
         end if
      end do
      Nc = NH + (Wdt - (diagram%order - 1)) * w
      sgn = sign(1.d0, w)
      if (abs(sgn - g_inc) > tol) n_sign_bad = n_sign_bad + 1
      err = abs(Nc / abs(w) - etrue_inc) / max(1.d0, abs(etrue_inc))
      max_o_rel = max(max_o_rel, err)
      if (err > tol) n_o_bad = n_o_bad + 1
      raw_num = Nc / abs(w)
      raw_den = sgn

      nwin = max(nv - 3, 0)
      if (.not. windows_on) nwin = 0
      rbn = 0.d0
      rbd = 0.d0
      rhos = 0.d0
      nel = 0
      ntile = 0
      do p = 1, nwin
         loc0 = -1
         do j = 0, 3
            partner = diagram%vertexList(ids(p + j))%link(2)
            jj = -1
            do r = 0, 3
               if (ids(p + r) == partner) jj = r
            end do
            if (jj < 0) exit
            if (j == 0) loc0 = jj
         end do
         if (jj < 0) then
            rbn = rbn + raw_num
            rbd = rbd + raw_den
            rhos = rhos + 1.d0
            cycle
         end if
         nel = nel + 1
         select case (loc0)
         case (1)
            r0 = 1; jb = 2
         case (3)
            r0 = 2; jb = 1
         case (2)
            r0 = 3; jb = 1
         case default
            stop 'rb_window: bad local pairing'
         end select
         vi => diagram%vertexList(ids(p))
         qa = vi%i_q; nua = vi%nu; wqa = vi%wq(nua)
         vp => diagram%vertexList(ids(p + jb))
         qb = vp%i_q; nub = vp%nu; wqb = vp%wq(nub)
         do j = 0, 3
            tv(j) = diagram%vertexList(ids(p + j))%tau
         end do
         L = X(:, :, p)
         Lam = Y(:, :, p)
         Renv = Z(:, :, p + 4)
         call colscale(Renv, Gd(:, p + 4))
         Rho = V(:, :, p + 4)
         call colscale(Rho, Gd(:, p + 4))
         C = Z(:, :, p + 4)
         call colscale(C, HGd(:, p + 4))
         Rho = Rho + C

         do r = 1, 3
            iq(pa(r, 1), :) = qa;  iq(pa(r, 2), :) = -qa;  nuv(pa(r, 1)) = nua; nuv(pa(r, 2)) = nua
            iq(pb(r, 1), :) = qb;  iq(pb(r, 2), :) = -qb;  nuv(pb(r, 1)) = nub; nuv(pb(r, 2)) = nub
            kin = diagram%vertexList(ids(p))%i_kin
            do j = 0, 3
               vs%i_kin = kin
               vs%i_q = iq(j, :)
               vs%i_kout = kin + iq(j, :)
               vs%nu = nuv(j)
               if (j == 0) then
                  vs%ukin = diagram%vertexList(ids(p))%ukin
               else
                  vs%ukin = ukprev
               end if
               if (j == 3) then
                  vs%ukout = diagram%vertexList(ids(p + 3))%ukout
               else
                  ktmp = kin + iq(j, :)
                  call cal_ek_int(ktmp, ek, uk)
                  vs%ukout = uk
                  ukprev = uk
                  Ej(:, j) = ek
               end if
               kin = kin + iq(j, :)
               call cal_gkq_vtex_int(vs, gbuf)
               grU(:, :, j) = gbuf(:, :, nuv(j))
            end do
            logPr(r) = 0.d0
            do j = 0, 3
               gm = maxval(abs(grU(:, :, j)))
               logPr(r) = logPr(r) + log(gm)
               gr(:, :, j) = grU(:, :, j) / gm
            end do
            T = gr(:, :, 0)
            TH = (0.d0, 0.d0)
            do j = 1, 3
               call seg_prop(Ej(:, j - 1), tv(j) - tv(j - 1), Gw, HGw, emin)
               logPr(r) = logPr(r) - emin * (tv(j) - tv(j - 1))
               A = T
               call rowscale(Gw, T)
               C = A
               call rowscale(HGw, C)
               call rowscale(Gw, TH)
               TH = TH + C
               T = matmul(gr(:, :, j), T)
               TH = matmul(gr(:, :, j), TH)
            end do
            ta = tv(pa(r, 2)) - tv(pa(r, 1))
            tb = tv(pb(r, 2)) - tv(pb(r, 1))
            logPr(r) = logPr(r) + log(Dph(wqa, ta)) + log(Dph(wqb, tb))
            OW(r) = wqa * ta + wqb * tb
            wr(r) = real(ctr(matmul(Renv, matmul(T, L))), dp)
            NHr(r) = real(ctr(matmul(Renv, matmul(T, Lam))) + ctr(matmul(Rho, matmul(T, L))) &
                          + ctr(matmul(Renv, matmul(TH, L))), dp)
            Twin(:, :, r) = T
            if (r == r0) then
               err = abs(wr(r) - w) / abs(w)
               do j = 0, 3
                  vp => diagram%vertexList(ids(p + j))
                  err = max(err, maxval(abs(grU(:, :, j) - vp%gkq(:, :, vp%nu))) / maxval(abs(vp%gkq(:, :, vp%nu))))
               end do
               do j = 0, 2
                  err = max(err, maxval(abs(Ej(:, j) - diagram%vertexList(ids(p + j))%ekout(1:nb))) &
                                 / max(1.d0, maxval(abs(Ej(:, j)))))
               end do
               max_r0_rel = max(max_r0_rel, err)
               if (err > tol) n_r0_bad = n_r0_bad + 1
            end if
         end do

         do r = 1, 3
            Nr(r) = NHr(r) + (Wdt - OW(r0) + OW(r) - (diagram%order - 1)) * wr(r)
         end do
         lmax = maxval(logPr)
         Pr = exp(logPr - lmax)
         Ssum = sum(Pr * abs(wr))
         if (Ssum <= 0.d0) then
            n_zero_fibre = n_zero_fibre + 1
            rbn = rbn + raw_num
            rbd = rbd + raw_den
            rhos = rhos + 1.d0
            cycle
         end if
         rbn = rbn + sum(Pr * Nr) / Ssum
         rbd = rbd + sum(Pr * wr) / Ssum
         rhos = rhos + abs(sum(Pr * wr)) / Ssum
         if (mod(p - 1, 4) == 0 .and. ntile < MAXT) then
            ntile = ntile + 1
            tpos(ntile) = p
            tr0(ntile) = r0
            TW(:, :, :, ntile) = Twin
            TlogP(:, ntile) = logPr - lmax
         end if
      end do
      if (nwin > 0) then
         rbn = rbn / nwin
         rbd = rbd / nwin
         rhos = rhos / nwin
      else
         rbn = raw_num
         rbd = raw_den
         rhos = 1.d0
      end if

      ! Exact tiled-group ratio over all 3^m members (Lemma 4 product structure).
      ! Group = re-pairings of the first min(m, MAXM_ENUM) eligible tiled windows. The
      ! subset depends only on the eligible set E(C), which is invariant under those
      ! re-pairings, so this is still a partition (Lemma 3) and E[rho] = Z_B/Z_A exactly.
      m = min(ntile, MAXM_ENUM)
      if (ntile > MAXM_ENUM) n_tiled_skip = n_tiled_skip + 1
      if (m == 0) then
         rho_t = 1.d0
      else
         Lf = X(:, :, tpos(1))
         Rl = Z(:, :, tpos(m) + 4)
         call colscale(Rl, Gd(:, tpos(m) + 4))
         do j = 1, m - 1
            call eye(C)
            do s = tpos(j) + 4, tpos(j + 1)
               if (s > tpos(j) + 4) C = matmul(gn(:, :, s - 1), C)
               call rowscale(Gd(:, s), C)
            end do
            Cb(:, :, j) = C
         end do
         ncomb = 3**m
         idx(1:m) = 1
         sumD = 0.d0
         sumA = 0.d0
         do k = 1, ncomb
            A = Lf
            lp = 0.d0
            do j = 1, m
               A = matmul(TW(:, :, idx(j), j), A)
               if (j < m) A = matmul(Cb(:, :, j), A)
               lp = lp + TlogP(idx(j), j)
            end do
            A = matmul(Rl, A)
            d = exp(lp) * real(ctr(A), dp)
            sumD = sumD + d
            sumA = sumA + abs(d)
            if (all(idx(1:m) == tr0(1:m))) then
               err = abs(real(ctr(A), dp) - w) / abs(w)
               max_tiled_rel = max(max_tiled_rel, err)
               if (err > tol) n_tiled_bad = n_tiled_bad + 1
            end if
            do j = 1, m
               idx(j) = idx(j) + 1
               if (idx(j) <= 3) exit
               idx(j) = 1
            end do
         end do
         rho_t = abs(sumD) / sumA
      end if

      if (modes_on) call rb_modes(diagram, ids, nv, w)
      if (mode_rb_on) call rb_mode_rb(diagram, ids, nv, w, Nc, Wdt, raw_num, raw_den)
      if (mode_svd_on) call rb_mode_svd(diagram, ids, nv)
      if (qgroup_on .and. mod(n_meas, int(q_every, 8)) == 0) call rb_qgroup(diagram, ids, nv, w)
      n_elig = n_elig + nel
      n_win = n_win + nwin
      write(rb_unit, '(4ES24.15,3I8,2ES24.15,2I6)') raw_num, raw_den, rbn, rbd, nel, nwin, diagram%order, &
           rhos, rho_t, ntile, m
      if (mod(n_meas, 500_8) == 0) call rb_summary()
   end subroutine rb_post


   ! Phonon-mode groups. A group varies the mode index nu of a chosen set of internal
   ! lines over all Nph modes, keeping topology, times and momenta. Lines are chosen
   ! by opening position only (nu-independent), so each rule is a partition and, under
   ! the native measure, E[rho] = Z_B/Z_A with rho = |sum_G D| / sum_G |D|.
   !   rho_mode_line : one line at a time (mean and min over the lines of the sample)
   !   rho_mode_K2/K3: the first 2 / 3 internal lines jointly (Nph^K members)
   !   rho_l1..3     : the first three lines individually (independence test)
   !   rho_K4        : first 4 lines jointly if FUTUREB_MODES_K4=1, on every
   !                   FUTUREB_MODES_K4_EVERY-th measurement (default 20; a fixed,
   !                   state-independent stride); -1 on the others
   !   prod_all_lines: prod over all lines of the single-line rho (independent-line
   !                   estimate of the all-lines mode group; validated against K2..K4)
   !   n_laminar, prod_laminar: greedy non-crossing line set (opening order) and the
   !                   product of its single-line rho; such a set can be mode-summed
   !                   exactly in O(Nph * order) by folding nested lines inside out
   !   pair_rel rho_pi rho_pj rho_pair: independence test on a counter-chosen line
   !                   pair (0 disjoint, 1 nested, 2 crossing): exact joint vs singles
   subroutine rb_modes(diagram, ids, nv, w)
      type(fynman), intent(inout), target :: diagram
      integer, intent(in) :: ids(maxN), nv
      real(dp), intent(in) :: w
      integer :: posof(maxN), la(maxN), lb(maxN), nl, j, s, l, nu, partner, kk, keff, k, ncomb, nsp, i
      integer :: sp(8), lineof(8), idx(4), tmpi, kmax
      real(dp) :: rsum, rmin, dn(Nph), sumd, suma, d, fr, err, dtl(4), rl(maxN), prodall
      real(dp) :: rho2, rho3, rho4, rk, prodlam, rhopair
      integer :: inlam(maxN), nlam, l2, pi_, pj_, prel, sel(2)
      integer, save :: pair_counter = 0
      logical :: ok
      complex(dp) :: Lm(nb,nb), Mm(nb,nb), Rm(nb,nb), Am(nb,nb), Pp(nb,nb,0:8), Rlast(nb,nb)
      type(vertex), pointer :: va, vb, vx

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
      if (nl == 0) then
         write(mode_unit, '(I6,10ES20.11,I6,ES20.11,I4,3ES20.11)') 0, 1.d0, 1.d0, 1.d0, 1.d0, 0.d0, &
              1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 0, 1.d0, -1, 1.d0, 1.d0, 1.d0
         return
      end if

      rsum = 0.d0
      rmin = 1.d0
      fr = 0.d0
      do l = 1, nl
         va => diagram%vertexList(ids(la(l)))
         vb => diagram%vertexList(ids(lb(l)))
         Lm = X(:, :, la(l))
         call eye(Mm)
         do s = la(l) + 1, lb(l)
            if (s > la(l) + 1) Mm = matmul(gn(:, :, s - 1), Mm)
            call rowscale(Gd(:, s), Mm)
         end do
         Rm = Z(:, :, lb(l) + 1)
         call colscale(Rm, Gd(:, lb(l) + 1))
         do nu = 1, Nph
            Am = matmul(Rm, matmul(vb%gkq(:, :, nu), matmul(Mm, matmul(va%gkq(:, :, nu), Lm))))
            dn(nu) = real(ctr(Am), dp) * Dph(va%wq(nu), vb%tau - va%tau)
         end do
         nu = va%nu
         err = abs(dn(nu) / (maxval(abs(va%gkq(:, :, nu))) * maxval(abs(vb%gkq(:, :, nu))) &
                             * Dph(va%wq(nu), vb%tau - va%tau)) - w) / abs(w)
         max_mode_rel = max(max_mode_rel, err)
         if (err > tol) n_mode_bad = n_mode_bad + 1
         if (sum(abs(dn)) > 0.d0) then
            d = abs(sum(dn)) / sum(abs(dn))
         else
            d = 1.d0
         end if
         rl(l) = d
         rsum = rsum + d
         rmin = min(rmin, d)
         if (d < 0.99d0) fr = fr + 1.d0
      end do

      prodall = product(rl(1:nl))
      if (modes_dump .and. mod(n_meas, int(dump_every, 8)) == 0) then
         write(dump_unit, '(A,I10,I6,F6.1)') 'M ', n_meas, nl, sign(1.d0, w)
         do l = 1, nl
            write(dump_unit, '(2I6,ES20.11)') la(l), lb(l), rl(l)
         end do
      end if
      ! Greedy non-crossing (laminar) set by opening order: nu-independent.
      nlam = 0
      prodlam = 1.d0
      do l = 1, nl
         ok = .true.
         do l2 = 1, nlam
            j = inlam(l2)
            if ((la(j) < la(l) .and. la(l) < lb(j) .and. lb(j) < lb(l)) .or. &
                (la(l) < la(j) .and. la(j) < lb(l) .and. lb(l) < lb(j))) then
               ok = .false.
               exit
            end if
         end do
         if (ok) then
            nlam = nlam + 1
            inlam(nlam) = l
            prodlam = prodlam * rl(l)
         end if
      end do
      ! Independence test on a line pair chosen by a state-independent counter rule.
      prel = -1
      rhopair = 1.d0
      pi_ = 1
      pj_ = 1
      if (nl >= 2) then
         pi_ = 1 + mod(pair_counter, nl)
         pj_ = 1 + mod(pi_ + mod(7 * pair_counter + 3, nl - 1), nl)
         if (pj_ == pi_) pj_ = 1 + mod(pi_, nl)
         pair_counter = pair_counter + 1
         if (lb(pi_) < la(pj_) .or. lb(pj_) < la(pi_)) then
            prel = 0                                          ! disjoint
         else if ((la(pi_) < la(pj_) .and. lb(pj_) < lb(pi_)) .or. &
                  (la(pj_) < la(pi_) .and. lb(pi_) < lb(pj_))) then
            prel = 1                                          ! nested
         else
            prel = 2                                          ! crossing
         end if
         sel(1) = pi_
         sel(2) = pj_
         call joint_mode_rho(diagram, ids, sel, 2, la, lb, rhopair)
      end if
      rho2 = 1.d0
      rho3 = 1.d0
      rho4 = -1.d0
      kmax = 3
      if (modes_k4 .and. mod(n_meas, int(k4_every, 8)) == 0) kmax = 4
      do kk = 2, kmax
         keff = min(kk, nl)
         nsp = 2 * keff
         do l = 1, keff
            sp(2 * l - 1) = la(l); lineof(2 * l - 1) = l
            sp(2 * l) = lb(l);     lineof(2 * l) = l
            dtl(l) = diagram%vertexList(ids(lb(l)))%tau - diagram%vertexList(ids(la(l)))%tau
         end do
         do i = 2, nsp
            do j = i, 2, -1
               if (sp(j - 1) > sp(j)) then
                  tmpi = sp(j); sp(j) = sp(j - 1); sp(j - 1) = tmpi
                  tmpi = lineof(j); lineof(j) = lineof(j - 1); lineof(j - 1) = tmpi
               end if
            end do
         end do
         Pp(:, :, 0) = X(:, :, sp(1))
         do j = 1, nsp - 1
            call eye(Mm)
            do s = sp(j) + 1, sp(j + 1)
               if (s > sp(j) + 1) Mm = matmul(gn(:, :, s - 1), Mm)
               call rowscale(Gd(:, s), Mm)
            end do
            Pp(:, :, j) = Mm
         end do
         Rlast = Z(:, :, sp(nsp) + 1)
         call colscale(Rlast, Gd(:, sp(nsp) + 1))
         ncomb = Nph**keff
         idx(1:keff) = 1
         sumd = 0.d0
         suma = 0.d0
         do k = 1, ncomb
            Am = Pp(:, :, 0)
            do j = 1, nsp
               vx => diagram%vertexList(ids(sp(j)))
               Am = matmul(vx%gkq(:, :, idx(lineof(j))), Am)
               if (j < nsp) Am = matmul(Pp(:, :, j), Am)
            end do
            Am = matmul(Rlast, Am)
            d = real(ctr(Am), dp)
            do l = 1, keff
               d = d * Dph(diagram%vertexList(ids(la(l)))%wq(idx(l)), dtl(l))
            end do
            sumd = sumd + d
            suma = suma + abs(d)
            do l = 1, keff
               idx(l) = idx(l) + 1
               if (idx(l) <= Nph) exit
               idx(l) = 1
            end do
         end do
         if (suma > 0.d0) then
            rk = abs(sumd) / suma
         else
            rk = 1.d0
         end if
         if (kk == 2) rho2 = rk
         if (kk == 3) rho3 = rk
         if (kk == 4) rho4 = rk
      end do
      write(mode_unit, '(I6,10ES20.11,I6,ES20.11,I4,3ES20.11)') nl, rsum / nl, rmin, rho2, rho3, fr / nl, rl(1), &
           merge(rl(2), 1.d0, nl >= 2), merge(rl(3), 1.d0, nl >= 3), rho4, prodall, nlam, prodlam, &
           prel, rl(pi_), rl(pj_), rhopair
   end subroutine rb_modes


   ! Momentum group for one internal line: vary its q over the whole nk_svd^3 grid,
   ! keeping topology, times, the mode index and every other line. The line has
   ! opening rank 1 + mod(counter, n_lines), a state-independent choice, so this is
   ! a partition and E[rho_q] = Z_B/Z_A. All segments inside the line are recomputed
   ! (energies, Bloch vectors) and all vertices a..b get new g; L and R are unchanged.
   ! For comparison the same line's mode group rho is written as well.
   subroutine rb_qgroup(diagram, ids, nv, w)
      type(fynman), intent(inout), target :: diagram
      integer, intent(in) :: ids(maxN), nv
      real(dp), intent(in) :: w
      integer, save :: counter = 0
      integer :: posof(maxN), la(maxN), lb(maxN), nl, j, s, partner, rank, a, b, nu, i1, i2, i3, nsp
      integer :: q0(3), q1(3), dq(3), kin(3), ktmp(3), iqv(3), nq, iqc
      real(dp) :: tau_a, tau_b, logp, lmax, err, ek(nb), Gw(nb), HGw(nb), emin, wqn(Nph), dnu(Nph)
      real(dp), allocatable :: dval(:), lval(:)
      complex(dp) :: Lm(nb,nb), Rm(nb,nb), Bq(nb,nb), uk(nbnd,nbnd), ukprev(nbnd,nbnd), gbuf(nb,nb,Nph)
      complex(dp) :: Mm(nb,nb), Am(nb,nb)
      real(dp) :: rho_q, rho_m
      type(vertex), pointer :: va, vb, vj

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
      if (nl == 0) return
      rank = 1 + mod(counter, nl)
      counter = counter + 1
      a = la(rank)
      b = lb(rank)
      va => diagram%vertexList(ids(a))
      vb => diagram%vertexList(ids(b))
      nu = va%nu
      q0 = va%i_q
      tau_a = va%tau
      tau_b = vb%tau
      Lm = X(:, :, a)
      Rm = Z(:, :, b + 1)
      call colscale(Rm, Gd(:, b + 1))

      nq = nk_svd**3
      allocate(dval(nq), lval(nq))
      iqc = 0
      do i1 = 0, nk_svd - 1
      do i2 = 0, nk_svd - 1
      do i3 = 0, nk_svd - 1
         iqc = iqc + 1
         q1 = (/i1, i2, i3/)
         dq = q1 - q0
         logp = 0.d0
         kin = va%i_kin
         do j = a, b
            vj => diagram%vertexList(ids(j))
            if (j == a) then
               iqv = q1
            else if (j == b) then
               iqv = -q1
            else
               iqv = vj%i_q
            end if
            vs%i_kin = kin
            vs%i_q = iqv
            vs%i_kout = kin + iqv
            vs%nu = vj%nu
            if (j == a) then
               vs%ukin = va%ukin
            else
               vs%ukin = ukprev
            end if
            if (j == b) then
               vs%ukout = vb%ukout
            else
               ktmp = kin + iqv
               call cal_ek_int(ktmp, ek, uk)
               vs%ukout = uk
               ukprev = uk
            end if
            call cal_gkq_vtex_int(vs, gbuf)
            if (j == a) then
               Bq = gbuf(:, :, vj%nu)
            else
               Bq = matmul(gbuf(:, :, vj%nu), Bq)
            end if
            if (j < b) then
               call seg_prop(ek, diagram%vertexList(ids(j + 1))%tau - vj%tau, Gw, HGw, emin)
               logp = logp - emin * (diagram%vertexList(ids(j + 1))%tau - vj%tau)
               call rowscale(Gw, Bq)
            end if
            kin = kin + iqv
         end do
         ktmp = q1
         call cal_wq_int(ktmp, wqn)
         logp = logp + log(Dph(wqn(nu), tau_b - tau_a))
         Am = matmul(Rm, matmul(Bq, Lm))
         dval(iqc) = real(ctr(Am), dp)
         lval(iqc) = logp
         if (all(q1 == modulo(q0, nk_svd))) then
            ! self-check against the current weight (normalized units)
            err = 0.d0
            Mm = (1.d0, 0.d0)
            do j = a, b
               vj => diagram%vertexList(ids(j))
               err = err + log(maxval(abs(vj%gkq(:, :, vj%nu))))
            end do
            err = abs(dval(iqc) / exp(err) - w) / abs(w)
            max_q_rel = max(max_q_rel, err)
            if (err > tol) n_q_bad = n_q_bad + 1
         end if
      end do
      end do
      end do
      lmax = maxval(lval)
      dval = dval * exp(lval - lmax)
      rho_q = abs(sum(dval)) / sum(abs(dval))
      deallocate(dval, lval)

      ! Same line, mode group (as in rb_modes).
      Mm = (0.d0, 0.d0)
      call eye(Mm)
      do s = a + 1, b
         if (s > a + 1) Mm = matmul(gn(:, :, s - 1), Mm)
         call rowscale(Gd(:, s), Mm)
      end do
      do j = 1, Nph
         Am = matmul(Rm, matmul(vb%gkq(:, :, j), matmul(Mm, matmul(va%gkq(:, :, j), Lm))))
         dnu(j) = real(ctr(Am), dp) * Dph(va%wq(j), tau_b - tau_a)
      end do
      rho_m = abs(sum(dnu)) / max(sum(abs(dnu)), tiny(1.d0))
      nsp = b - a
      write(q_unit, '(I10,3I6,2ES20.11)') n_meas, nl, rank, nsp, rho_q, rho_m
   end subroutine rb_qgroup


   ! Exact joint mode-group ratio for an arbitrary set of lines (Nph^K members).
   subroutine joint_mode_rho(diagram, ids, sel, keff, la, lb, rho)
      type(fynman), intent(inout), target :: diagram
      integer, intent(in) :: ids(maxN), keff, sel(keff), la(maxN), lb(maxN)
      real(dp), intent(out) :: rho
      integer :: sp(2*keff), lineof(2*keff), idx(keff), nsp, i, j, s, l, k, ncomb, tmpi
      real(dp) :: dtl(keff), sumd, suma, d
      complex(dp) :: Pp(nb,nb,0:2*keff), Mm(nb,nb), Am(nb,nb), Rlast(nb,nb)
      type(vertex), pointer :: vx
      nsp = 2 * keff
      do l = 1, keff
         sp(2 * l - 1) = la(sel(l)); lineof(2 * l - 1) = l
         sp(2 * l) = lb(sel(l));     lineof(2 * l) = l
         dtl(l) = diagram%vertexList(ids(lb(sel(l))))%tau - diagram%vertexList(ids(la(sel(l))))%tau
      end do
      do i = 2, nsp
         do j = i, 2, -1
            if (sp(j - 1) > sp(j)) then
               tmpi = sp(j); sp(j) = sp(j - 1); sp(j - 1) = tmpi
               tmpi = lineof(j); lineof(j) = lineof(j - 1); lineof(j - 1) = tmpi
            end if
         end do
      end do
      Pp(:, :, 0) = X(:, :, sp(1))
      do j = 1, nsp - 1
         call eye(Mm)
         do s = sp(j) + 1, sp(j + 1)
            if (s > sp(j) + 1) Mm = matmul(gn(:, :, s - 1), Mm)
            call rowscale(Gd(:, s), Mm)
         end do
         Pp(:, :, j) = Mm
      end do
      Rlast = Z(:, :, sp(nsp) + 1)
      call colscale(Rlast, Gd(:, sp(nsp) + 1))
      ncomb = Nph**keff
      idx = 1
      sumd = 0.d0
      suma = 0.d0
      do k = 1, ncomb
         Am = Pp(:, :, 0)
         do j = 1, nsp
            vx => diagram%vertexList(ids(sp(j)))
            Am = matmul(vx%gkq(:, :, idx(lineof(j))), Am)
            if (j < nsp) Am = matmul(Pp(:, :, j), Am)
         end do
         Am = matmul(Rlast, Am)
         d = real(ctr(Am), dp)
         do l = 1, keff
            d = d * Dph(diagram%vertexList(ids(la(sel(l))))%wq(idx(l)), dtl(l))
         end do
         sumd = sumd + d
         suma = suma + abs(d)
         do l = 1, keff
            idx(l) = idx(l) + 1
            if (idx(l) <= Nph) exit
            idx(l) = 1
         end do
      end do
      rho = 1.d0
      if (suma > 0.d0) rho = abs(sumd) / suma
   end subroutine joint_mode_rho


   ! Mode Rao-Blackwellization of the native EZ estimator (measurement only).
   ! Line l (opening rank, mode-independent) defines a fibre: all Nph modes of l,
   ! everything else fixed. Under pi_A ~ |Re D|,
   !   E[f_O | fibre_l] = sum_nu N_nu / sum_nu |D_nu|,  E[f_1 | fibre_l] = sum_nu D_nu / sum_nu |D_nu|,
   ! and the mean over lines is unbiased. With M the interior of l (and M_H its
   ! energy insertions), L = X_a, Lam = Y_a, R = Z_{b+1} G_{b+1}, Rho = V_{b+1} G_{b+1} + Z_{b+1} HG_{b+1}:
   !   D_nu  = Re tr(R g_b M g_a L) Dph_nu
   !   NH_nu = Re[tr(R g_b M g_a Lam) + tr(Rho g_b M g_a L) + tr(R g_b M_H g_a L)] Dph_nu
   !   N_nu  = NH_nu + (Wdt - w_cur dtau + w_nu dtau - (order-1)) D_nu
   ! Self-check: the current mode reproduces w and the native numerator Nc.
   subroutine rb_mode_rb(diagram, ids, nv, w, Nc, Wdt, raw_num, raw_den)
      type(fynman), intent(inout), target :: diagram
      integer, intent(in) :: ids(maxN), nv
      real(dp), intent(in) :: w, Nc, Wdt, raw_num, raw_den
      integer :: posof(maxN), j, s, l, nl, partner, a, b, nu
      real(dp) :: dn(Nph), nn(Nph), dtl, dphv, sabs, mnum, mden, fac, err
      complex(dp) :: Lm(nb,nb), Lam(nb,nb), Rm(nb,nb), Rho(nb,nb), Mm(nb,nb), MH(nb,nb), Am(nb,nb), Cm(nb,nb)
      complex(dp) :: T(nb,nb), TH(nb,nb)
      type(vertex), pointer :: va, vb
      posof = 0
      do j = 1, nv
         posof(ids(j)) = j
      end do
      nl = 0
      mnum = 0.d0
      mden = 0.d0
      do j = 1, nv
         partner = diagram%vertexList(ids(j))%link(2)
         if (partner == 1 .or. partner == maxN) cycle
         if (posof(partner) <= j) cycle
         nl = nl + 1
         a = j
         b = posof(partner)
         va => diagram%vertexList(ids(a))
         vb => diagram%vertexList(ids(b))
         dtl = vb%tau - va%tau
         Lm = X(:, :, a)
         Lam = Y(:, :, a)
         Rm = Z(:, :, b + 1)
         call colscale(Rm, Gd(:, b + 1))
         Rho = V(:, :, b + 1)
         call colscale(Rho, Gd(:, b + 1))
         Cm = Z(:, :, b + 1)
         call colscale(Cm, HGd(:, b + 1))
         Rho = Rho + Cm
         call eye(Mm)
         MH = (0.d0, 0.d0)
         do s = a + 1, b
            if (s > a + 1) then
               Mm = matmul(gn(:, :, s - 1), Mm)
               MH = matmul(gn(:, :, s - 1), MH)
            end if
            Am = Mm
            call rowscale(Gd(:, s), Mm)
            Cm = Am
            call rowscale(HGd(:, s), Cm)
            call rowscale(Gd(:, s), MH)
            MH = MH + Cm
         end do
         do nu = 1, Nph
            dphv = Dph(va%wq(nu), dtl)
            T = matmul(vb%gkq(:, :, nu), matmul(Mm, va%gkq(:, :, nu)))
            TH = matmul(vb%gkq(:, :, nu), matmul(MH, va%gkq(:, :, nu)))
            dn(nu) = real(ctr(matmul(Rm, matmul(T, Lm))), dp) * dphv
            nn(nu) = real(ctr(matmul(Rm, matmul(T, Lam))) + ctr(matmul(Rho, matmul(T, Lm))) &
                          + ctr(matmul(Rm, matmul(TH, Lm))), dp) * dphv
            nn(nu) = nn(nu) + (Wdt - va%wq(va%nu) * dtl + va%wq(nu) * dtl - (diagram%order - 1)) * dn(nu)
         end do
         nu = va%nu
         fac = maxval(abs(va%gkq(:, :, nu))) * maxval(abs(vb%gkq(:, :, nu))) * Dph(va%wq(nu), dtl)
         err = max(abs(dn(nu) / fac - w) / abs(w), abs(nn(nu) / fac - Nc) / max(abs(Nc), abs(w)))
         max_mrb_rel = max(max_mrb_rel, err)
         if (err > tol) n_mrb_bad = n_mrb_bad + 1
         sabs = sum(abs(dn))
         if (sabs > 0.d0) then
            mnum = mnum + sum(nn) / sabs
            mden = mden + sum(dn) / sabs
         else
            mnum = mnum + raw_num
            mden = mden + raw_den
         end if
      end do
      if (nl == 0) then
         write(mrb_unit, '(4ES24.15,I6)') raw_num, raw_den, raw_num, raw_den, 0
      else
         write(mrb_unit, '(4ES24.15,I6)') raw_num, raw_den, mnum / nl, mden / nl, nl
      end if
   end subroutine rb_mode_rb


   ! Sign-optimized mode basis diagnostic. For a line (a,b) the mode sum is the
   ! bilinear T_(ik),(lj) = sum_nu Dph_nu g_b(nu)_ik g_a(nu)_lj (rank <= Nph). Its SVD
   ! T = sum_mu s_mu u_mu v_mu^H gives effective modes mu with the same total:
   ! line block = sum_mu s_mu U_mu M W_mu, U_mu = u_mu (i,k), W_mu = conj(v_mu) (l,j).
   ! gamma_l = sum_mu |Re tr(R s U M W L)| / sum_nu |Re tr(R g_b M g_a L) Dph|.
   ! On a single-line fibre E_A[gamma] = Z_A(mu basis)/Z_A exactly; gamma < 1 means
   ! that sampling mu instead of nu has less sign cancellation on that line.
   subroutine rb_mode_svd(diagram, ids, nv)
      type(fynman), intent(inout), target :: diagram
      integer, intent(in) :: ids(maxN), nv
      integer, save :: pair_ctr = 0
      integer :: posof(maxN), la(maxN), lb(maxN), rnk(maxN), j, nl, partner, l, p1, p2, nbb2
      real(dp) :: gam(maxN), rho(maxN), snu, smu, dval, gsum, gprod, rsum, rprod, rk, jn, jm, pn(2), pm(2)
      complex(dp) :: EU(nb, nb, 9, maxN), EW(nb, nb, 9, maxN)
      real(dp) :: sig(9, maxN)
      logical :: dump
      nbb2 = nb * nb
      if (nbb2 > 9) stop 'rb_mode_svd: nb <= 3 only'
      posof = 0
      do j = 1, nv
         posof(ids(j)) = j
      end do
      nl = 0
      do j = 1, nv
         partner = diagram%vertexList(ids(j))%link(2)
         if (partner == 1 .or. partner == maxN) cycle
         if (posof(partner) <= j) cycle
         nl = nl + 1
         la(nl) = j
         lb(nl) = posof(partner)
      end do
      gsum = 0.d0; gprod = 1.d0; rsum = 0.d0; rprod = 1.d0; rk = 0.d0
      do l = 1, nl
         call line_bases(diagram, ids, la(l), lb(l), snu, smu, dval, rnk(l), sig(:, l), &
                         EU(:, :, :, l), EW(:, :, :, l))
         gam(l) = 1.d0
         rho(l) = 1.d0
         if (snu > 0.d0) then
            gam(l) = smu / snu
            rho(l) = abs(dval) / snu
         end if
         gsum = gsum + gam(l); gprod = gprod * gam(l)
         rsum = rsum + rho(l); rprod = rprod * rho(l)
         rk = rk + rnk(l)
      end do
      ! exact two-line basis change on a counter-chosen pair (independence test)
      jn = -1.d0; jm = -1.d0; pn = 1.d0; pm = 1.d0; p1 = 0; p2 = 0
      if (nl >= 2) then
         p1 = 1 + mod(pair_ctr, nl)
         p2 = 1 + mod(p1 + mod(7 * pair_ctr + 3, nl - 1), nl)
         pair_ctr = pair_ctr + 1
         call pair_bases(diagram, ids, la, lb, p1, p2, rnk, sig, EU, EW, jn, jm)
      end if
      if (nl == 0) then
         write(svd_unit, '(I6,4ES20.11,F8.3,2I6,4ES20.11)') 0, 1.d0, 1.d0, 1.d0, 1.d0, 0.d0, 0, 0, -1.d0, -1.d0, -1.d0, -1.d0
      else
         write(svd_unit, '(I6,4ES20.11,F8.3,2I6,4ES20.11)') nl, gsum / nl, gprod, rsum / nl, rprod, rk / nl, &
              p1, p2, jm, merge(gam(max(p1,1)) * gam(max(p2,1)), -1.d0, nl >= 2), jn, &
              merge(rho(max(p1,1)) * rho(max(p2,1)), -1.d0, nl >= 2)
      end if
      dump = modes_dump .and. mod(n_meas, int(dump_every, 8)) == 0
      if (dump) then
         write(dump_unit, '(A,I10,I6,F6.1)') 'G ', n_meas, nl, 0.d0
         do l = 1, nl
            write(dump_unit, '(2I6,2ES20.11)') la(l), lb(l), gam(l), rho(l)
         end do
      end if
   end subroutine rb_mode_svd

   ! One line: nu-basis sum |D|, SVD-basis sum |D|, the signed total, and the SVD
   ! effective matrices (EU at the closer, s*EW at the opener).
   subroutine line_bases(diagram, ids, a, b, snu, smu, dval, rank, sv9, EU, EW)
      type(fynman), intent(inout), target :: diagram
      integer, intent(in) :: ids(maxN), a, b
      real(dp), intent(out) :: snu, smu, dval, sv9(9)
      integer, intent(out) :: rank
      complex(dp), intent(out) :: EU(nb, nb, 9), EW(nb, nb, 9)
      integer :: s, nu, mu, i, k, jj, l2, info, lwork, nbb2
      real(dp) :: dtl, dphv, sv(nb*nb), rwork(8*nb*nb)
      complex(dp) :: Lm(nb,nb), Rm(nb,nb), Mm(nb,nb), Am(nb,nb), Um(nb,nb), Wm(nb,nb)
      complex(dp) :: Tm(nb*nb, nb*nb), Us(nb*nb, nb*nb), VTs(nb*nb, nb*nb), work(8*nb*nb*nb*nb)
      type(vertex), pointer :: va, vb
      nbb2 = nb * nb
      lwork = size(work)
      va => diagram%vertexList(ids(a))
      vb => diagram%vertexList(ids(b))
      dtl = vb%tau - va%tau
      Lm = X(:, :, a)
      Rm = Z(:, :, b + 1)
      call colscale(Rm, Gd(:, b + 1))
      call eye(Mm)
      do s = a + 1, b
         if (s > a + 1) Mm = matmul(gn(:, :, s - 1), Mm)
         call rowscale(Gd(:, s), Mm)
      end do
      Tm = (0.d0, 0.d0)
      snu = 0.d0
      dval = 0.d0
      do nu = 1, Nph
         dphv = Dph(va%wq(nu), dtl)
         Am = matmul(vb%gkq(:, :, nu), matmul(Mm, va%gkq(:, :, nu)))
         snu = snu + abs(real(ctr(matmul(Rm, matmul(Am, Lm))), dp)) * dphv
         dval = dval + real(ctr(matmul(Rm, matmul(Am, Lm))), dp) * dphv
         do k = 1, nb
         do i = 1, nb
            do jj = 1, nb
            do l2 = 1, nb
               Tm(i + (k - 1) * nb, l2 + (jj - 1) * nb) = Tm(i + (k - 1) * nb, l2 + (jj - 1) * nb) &
                    + dphv * vb%gkq(i, k, nu) * va%gkq(l2, jj, nu)
            end do
            end do
         end do
         end do
      end do
      call zgesvd('A', 'A', nbb2, nbb2, Tm, nbb2, sv, Us, nbb2, VTs, nbb2, work, lwork, rwork, info)
      smu = 0.d0
      rank = 0
      sv9 = 0.d0
      do mu = 1, nbb2
         if (sv(mu) <= 1.d-14 * sv(1)) cycle
         rank = rank + 1
         do k = 1, nb
         do i = 1, nb
            Um(i, k) = Us(i + (k - 1) * nb, mu)
         end do
         end do
         do jj = 1, nb
         do l2 = 1, nb
            Wm(l2, jj) = VTs(mu, l2 + (jj - 1) * nb)
         end do
         end do
         EU(:, :, rank) = Um
         EW(:, :, rank) = sv(mu) * Wm
         sv9(rank) = sv(mu)
         Am = matmul(Um, matmul(Mm, EW(:, :, rank)))
         smu = smu + abs(real(ctr(matmul(Rm, matmul(Am, Lm))), dp))
      end do
   end subroutine line_bases

   ! Exact joint sums over two lines: nu x nu (jn = sum|D|) and mu x mu (jm = sum|D|),
   ! both divided by the signed total, i.e. the joint 1/rho-like normalizers; returns
   ! jm/jn relative to the product of single-line gammas via the caller.
   subroutine pair_bases(diagram, ids, la, lb, p1, p2, rnk, sig, EU, EW, jn, jm)
      type(fynman), intent(inout), target :: diagram
      integer, intent(in) :: ids(maxN), la(maxN), lb(maxN), p1, p2, rnk(maxN)
      real(dp), intent(in) :: sig(9, maxN)
      complex(dp), intent(in) :: EU(nb, nb, 9, maxN), EW(nb, nb, 9, maxN)
      real(dp), intent(out) :: jn, jm
      integer :: sp(4), lineof(4), i, j, tmpi, s, k1, k2, nsp
      real(dp) :: dtl1, dtl2, d, sn, sm, tot
      complex(dp) :: Pp(nb, nb, 0:4), Mm(nb, nb), Am(nb, nb), Rl(nb, nb), mat(nb, nb, 4)
      type(vertex), pointer :: vx
      sp = (/ la(p1), lb(p1), la(p2), lb(p2) /)
      lineof = (/ 1, 1, 2, 2 /)
      do i = 2, 4
         do j = i, 2, -1
            if (sp(j - 1) > sp(j)) then
               tmpi = sp(j); sp(j) = sp(j - 1); sp(j - 1) = tmpi
               tmpi = lineof(j); lineof(j) = lineof(j - 1); lineof(j - 1) = tmpi
            end if
         end do
      end do
      nsp = 4
      Pp(:, :, 0) = X(:, :, sp(1))
      do j = 1, nsp - 1
         call eye(Mm)
         do s = sp(j) + 1, sp(j + 1)
            if (s > sp(j) + 1) Mm = matmul(gn(:, :, s - 1), Mm)
            call rowscale(Gd(:, s), Mm)
         end do
         Pp(:, :, j) = Mm
      end do
      Rl = Z(:, :, sp(nsp) + 1)
      call colscale(Rl, Gd(:, sp(nsp) + 1))
      dtl1 = diagram%vertexList(ids(lb(p1)))%tau - diagram%vertexList(ids(la(p1)))%tau
      dtl2 = diagram%vertexList(ids(lb(p2)))%tau - diagram%vertexList(ids(la(p2)))%tau
      sn = 0.d0
      tot = 0.d0
      do k1 = 1, Nph
         do k2 = 1, Nph
            do j = 1, nsp
               vx => diagram%vertexList(ids(sp(j)))
               if (lineof(j) == 1) then
                  mat(:, :, j) = vx%gkq(:, :, k1)
               else
                  mat(:, :, j) = vx%gkq(:, :, k2)
               end if
            end do
            Am = Pp(:, :, 0)
            do j = 1, nsp
               Am = matmul(mat(:, :, j), Am)
               if (j < nsp) Am = matmul(Pp(:, :, j), Am)
            end do
            Am = matmul(Rl, Am)
            d = real(ctr(Am), dp) * Dph(diagram%vertexList(ids(la(p1)))%wq(k1), dtl1) &
                                  * Dph(diagram%vertexList(ids(la(p2)))%wq(k2), dtl2)
            sn = sn + abs(d)
            tot = tot + d
         end do
      end do
      sm = 0.d0
      do k1 = 1, rnk(p1)
         do k2 = 1, rnk(p2)
            do j = 1, nsp
               if (lineof(j) == 1) then
                  if (sp(j) == la(p1)) then
                     mat(:, :, j) = EW(:, :, k1, p1)
                  else
                     mat(:, :, j) = EU(:, :, k1, p1)
                  end if
               else
                  if (sp(j) == la(p2)) then
                     mat(:, :, j) = EW(:, :, k2, p2)
                  else
                     mat(:, :, j) = EU(:, :, k2, p2)
                  end if
               end if
            end do
            Am = Pp(:, :, 0)
            do j = 1, nsp
               Am = matmul(mat(:, :, j), Am)
               if (j < nsp) Am = matmul(Pp(:, :, j), Am)
            end do
            Am = matmul(Rl, Am)
            sm = sm + abs(real(ctr(Am), dp))
         end do
      end do
      jn = abs(tot) / sn          ! joint rho (nu basis, both lines summed)
      jm = sm / sn                 ! joint gamma (mu x mu vs nu x nu)
   end subroutine pair_bases

   subroutine rb_summary()
      integer :: u
      open(newunit=u, file='rb_summary.dat', status='replace', action='write')
      write(u, '(A,I12)') 'n_meas ', n_meas
      write(u, '(A,I12)') 'n_sign_mismatch ', n_sign_bad
      write(u, '(A,I12)') 'n_O_mismatch ', n_o_bad
      write(u, '(A,ES12.4)') 'max_O_rel ', max_o_rel
      write(u, '(A,I12)') 'n_env_mismatch ', n_env_bad
      write(u, '(A,ES12.4)') 'max_env_rel ', max_env_rel
      write(u, '(A,I12)') 'n_r0_mismatch ', n_r0_bad
      write(u, '(A,ES12.4)') 'max_r0_rel ', max_r0_rel
      write(u, '(A,I12)') 'n_tiled_mismatch ', n_tiled_bad
      write(u, '(A,ES12.4)') 'max_tiled_rel ', max_tiled_rel
      write(u, '(A,I12)') 'n_tiled_capped ', n_tiled_skip
      write(u, '(A,I12)') 'n_eligible_windows ', n_elig
      write(u, '(A,I12)') 'n_windows ', n_win
      write(u, '(A,I12)') 'n_zero_fibre ', n_zero_fibre
      write(u, '(A,I12)') 'n_mode_mismatch ', n_mode_bad
      write(u, '(A,I12)') 'n_mode_rb_mismatch ', n_mrb_bad
      write(u, '(A,ES12.4)') 'max_mode_rb_rel ', max_mrb_rel
      write(u, '(A,I12)') 'n_q_mismatch ', n_q_bad
      write(u, '(A,ES12.4)') 'max_q_rel ', max_q_rel
      write(u, '(A,ES12.4)') 'max_mode_rel ', max_mode_rel
      close(u)
      flush(rb_unit)
      if (modes_on) flush(mode_unit)
      if (mode_rb_on) flush(mrb_unit)
      if (mode_svd_on) flush(svd_unit)
      if (qgroup_on) flush(q_unit)
   end subroutine rb_summary

end module rb_window_mod
