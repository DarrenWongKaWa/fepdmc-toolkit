! Excerpt of update_swap from clean P1 binary source (not the full Perturbo file).
! For source-structure audit only. GPL-3 Perturbo-derived.

subroutine update_swap(diagram, stat)
      use DiagMC 
      use linear_da_mod, only: linear_da_hook, linear_da_is_on, &
           linear_da_eval, linear_da_aux_uniform, linear_ld, c2_dump_swap, p1_prop_ell, p1_fixture_on, &
           n_stage1_rej, n_stage2, n_accept

      implicit none 
      type(fynman), target :: diagram
      type(diagmc_stat) :: stat
      integer :: iv1, iv2, ib
      type(vertex), pointer :: vL,vR, vLL, vRR

      real(dp) :: tau1, tau2, tau12,tauL,tauR,tauLL,tauRR, decay, wq1, wq2,Elist(dmc_band), dtau
      real(dp) :: P_accept, P_kchange
   
      real(dp) :: ran
      real(dp) ::  factor 
      type(vertex),pointer :: vLn,vRn
      complex(dp) :: left_env(dmc_band,dmc_band), right_env(dmc_band,dmc_band), A(dmc_band, dmc_band), left_new(dmc_band,dmc_band), mat_old, mat_new 
      complex(dp) :: expH(dmc_band,dmc_band)
      logical :: da_el, da_acc
      real(dp) :: ell_hat, ellR, u1, ell_rev
      integer*8 :: ng_enter, ne_enter
      if(diagram%order<=3)return 

      !randomly choose a pair of vertex 
      call uniform_int_omp(diagram%seed,2,diagram%order,iv1)
      vL => diagram%vertexList(iv1);  iv2 = vL%link(3)
      if(vL%link(2).eq.iv2) return            !the two vertex are not phonon linked 
      if(iv2.eq.maxN)       return            !not the last vertex   
      !count update 
      stat%update(7) = stat%update(7) + 1; call c1_mark_valid(7)

      vR => diagram%vertexList(iv2)
      !the two vertex  that connected directly to the external should not cross 
      if(vL%link(2).eq. 1 .and. vR%link(2) .eq. maxN) return 
      if(vL%link(2).eq. maxN .and. vR%link(2) .eq. 1) return 
      
      vLL => diagram%vertexList(vL%link(1))
      vRR => diagram%vertexList(vR%link(3))
      tauL = vL%tau; tauLL = vLL%tau; 
      tauR = vR%tau; tauRR = vRR%tau

      if(tauL>tauR) stop ' time disordered @ update_swap'

      da_el = .false.
      ell_hat = 0.d0
      ng_enter = c1_ncall(c1_g_prod)+c1_ncall(c1_g_full)+c1_ncall(c1_epc)
      ne_enter = c1_ncall(c1_env_left)+c1_ncall(c1_env_right)+c1_ncall(c1_env_mid)
      if (p1_fixture_on()) call c2_dump_swap(10, iv1, iv2, diagram, 0.d0, 0.d0, 0_8, 0_8, 0.d0, 0.d0, 0.d0, 0.d0)
      if (linear_da_hook()) then
         call c1_mark_hook()
         call linear_da_eval(diagram, iv1, iv2, da_el, ell_hat)
         if (da_el) then
            call c1_mark_eligible()
            call c1_mark_score()
         endif
      endif
      if (linear_da_is_on() .and. da_el) then
         call linear_da_aux_uniform(u1)
         if (log(max(u1, 1.0e-300_dp)) >= min(0.d0, ell_hat)) then
            n_stage1_rej = n_stage1_rej + 1
            if (p1_fixture_on()) call c2_dump_swap(1, iv1, iv2, diagram, ell_hat, 0.d0, &
              c1_ncall(c1_g_prod)+c1_ncall(c1_g_full)+c1_ncall(c1_epc)-ng_enter, &
              c1_ncall(c1_env_left)+c1_ncall(c1_env_right)+c1_ncall(c1_env_mid)-ne_enter, &
              minval(vL%ekout), 0.d0, 0.d0, 0.d0)
            return
         endif
      endif

      vLn => diagram%vertexList(diagram%order+1) 
      vRn => diagram%vertexList(diagram%order+2) 


      vRn%tau = vL%tau
      !q depdent quantity 
      vRn%nu = vR%nu; vRn%i_q = vR%i_q; vRn%Pq = vR%Pq
      vRn%wq = vR%wq; 
      call c0_copy_wq(int(loc(vRn%wq),8),int(loc(vR%wq),8))
      !vRn%uk_f = vR%uk_f; vRn%ukq_f = vR%ukq_f
      !vRn%vq_f = vR%vq_f; vRn%vnq_f = vR%vnq_f
      
      !kin depdent quantity 
      vRn%i_kin = vL%i_kin; vRn%ekin = vL%ekin; vRn%vkin = vL%vkin; vRn%ukin = vL%ukin;  
      !kout depdent quantity 
      vRn%i_kout = vRn%i_kin + vRn%i_q
      call cal_ek_int(vRn%i_kout, vRn%ekout,vRn%ukout, vRn%vkout)
      !new e-ph vertex@vRn 
      call c1_mark_exact()
      call cal_gkq_vtex_int( vRn, vRn%gkq ) 
      call cal_gkq_full_vtex_int( vRn, vRn%gkq_full ) 

      vLn%tau = vR%tau
      !q depdent quantity 
      vLn%nu = vL%nu; vLn%i_q = vL%i_q; vLn%Pq = vL%Pq
      vLn%wq = vL%wq; 
      call c0_copy_wq(int(loc(vLn%wq),8),int(loc(vL%wq),8))

      !kin depdent quantity 
      vLn%i_kin = vRn%i_kout; vLn%ekin = vRn%ekout; vLn%vkin = vRn%vkout; vLn%ukin = vRn%ukout 
      !kout depdent quantity 
      vLn%i_kout = vR%i_kout; vLn%ekout = vR%ekout; vLn%vkout = vR%vkout; vLn%ukout = vR%ukout 
      !new e-ph vertex@vLn
      call cal_gkq_vtex_int( vLn, vLn%gkq ) 
      call cal_gkq_full_vtex_int( vLn, vLn%gkq_full ) 

      !electron propagator ratio 
      tau12 =  vR%tau - vL%tau
      P_kchange = Gel( minval(vRn%ekout),tau12 ) / Gel( minval(vL%ekout),tau12 )

      !ph propagator ratio 
      wq1 = vL%wq(vL%nu); wq2 = vR%wq(vR%nu)
      P_kchange = P_kchange * Dph( wq1, abs( vR%tau - diagram%vertexList(vL%link(2))%tau ) ) /&
                  Dph( wq1, abs( vL%tau - diagram%vertexList(vL%link(2))%tau ) )
      P_kchange = P_kchange * Dph( wq2, abs( vL%tau - diagram%vertexList(vR%link(2))%tau ) ) /&
                  Dph( wq2, abs( vR%tau - diagram%vertexList(vR%link(2))%tau ) )
       
      !matrix term for new e-ph and propagator 
      call left_environment_matrix(diagram, vL%link(1), left_env);
      Elist = vL%ekin - minval(vL%ekin)
      dtau =  vL%tau - vLL%tau
      call propagator(vL%ukin, Elist, dtau, dmc_band, expH)
      left_env = matmul(expH, left_env)
      left_new = left_env

      call right_environment_matrix(diagram,  vR%link(3), right_env)
      Elist = vR%ekout - minval(vR%ekout)
      dtau =  vRR%tau - vR%tau
      call propagator(vR%ukout, Elist, dtau, dmc_band, expH)
      right_env = matmul(right_env, expH) 

      !old mat 
      !e-ph  @vL
      left_env = matmul(vL%gkq(:,:,vL%nu), left_env)
      !propagator : vL -> vR
      Elist = vL%ekout - minval(vL%ekout)
      dtau =  vR%tau - vL%tau
      call propagator(vL%ukout, Elist, dtau, dmc_band, expH)
      left_env = matmul(expH, left_env)
      !e-ph @ vR
      left_env = matmul(vR%gkq(:,:,vR%nu), left_env)

      A = matmul(right_env,left_env)
      mat_old = trace(A, dmc_band)
      if(sample_gt) mat_old = A(ib_sample, ib_sample)

      !new mat 
      !e-ph @vRn
      left_new = matmul(vRn%gkq(:,:,vRn%nu), left_new)
      !propagator : vRn -> vLn 
      Elist = vRn%ekout - minval(vRn%ekout)
      call propagator(vRn%ukout, Elist, dtau, dmc_band, expH)
      left_new = matmul(expH, left_new)
      !e-ph @vLn 
      left_new = matmul(vLn%gkq(:,:,vLn%nu), left_new)

      A = matmul(right_env, left_new)
      mat_new = trace(A, dmc_band)
      if(sample_gt) mat_new = A(ib_sample, ib_sample)

      !for Frohlich case A should be hermition check it 
      !write(*,'(2E20.10)')maxval(abs(A - transpose(conjg(A)))) / maxval(abs(A)), maxval(abs(A))

      factor = real(mat_new) / real(mat_old)
      !write(*,'(3E20.10)')  maxval(abs(left_env)),maxval(abs(right_env)),P_kchange
      P_accept = abs(factor) * P_kchange
      ellR = 0.d0
      if (p1_fixture_on() .and. real(mat_old)==real(mat_old) .and. abs(real(mat_old))>0.d0 .and. &
          real(mat_new)==real(mat_new) .and. abs(real(mat_new))>0.d0) then
         ellR = log(abs(real(mat_new))) - log(abs(real(mat_old)))
         ellR = ellR - (minval(vRn%ekout) - minval(vL%ekout)) * (tauR - tauL)
         ellR = ellR + linear_ld(wq1, tauR - diagram%vertexList(vL%link(2))%tau) - &
              linear_ld(wq1, tauL - diagram%vertexList(vL%link(2))%tau)
         ellR = ellR + linear_ld(wq2, tauL - diagram%vertexList(vR%link(2))%tau) - &
              linear_ld(wq2, tauR - diagram%vertexList(vR%link(2))%tau)
         call c2_dump_swap(3, iv1, iv2, diagram, ell_hat, ellR, &
           c1_ncall(c1_g_prod)+c1_ncall(c1_g_full)+c1_ncall(c1_epc)-ng_enter, &
           c1_ncall(c1_env_left)+c1_ncall(c1_env_right)+c1_ncall(c1_env_mid)-ne_enter, &
           minval(vL%ekout), minval(vRn%ekout), real(mat_old), real(mat_new))
      endif
      call random_number_omp(diagram%seed,ran)

      da_acc = (ran < P_accept)
      if (linear_da_is_on() .and. da_el) then
         if (.not. (real(mat_old) == real(mat_old)) .or. abs(real(mat_old)) == 0.d0) then
            stop 'linear_da: Re(mat_old) is 0 or nonfinite'
         endif
         if (abs(real(mat_new)) == 0.d0) then
            da_acc = .false.
         else
            ellR = log(abs(real(mat_new))) - log(abs(real(mat_old)))
            ellR = ellR - (minval(vRn%ekout) - minval(vL%ekout)) * (tauR - tauL)
            ellR = ellR + linear_ld(wq1, tauR - diagram%vertexList(vL%link(2))%tau) - &
                 linear_ld(wq1, tauL - diagram%vertexList(vL%link(2))%tau)
            ellR = ellR + linear_ld(wq2, tauL - diagram%vertexList(vR%link(2))%tau) - &
                 linear_ld(wq2, tauR - diagram%vertexList(vR%link(2))%tau)
            da_acc = (log(max(ran, 1.0e-300_dp)) < min(0.d0, ellR - ell_hat))
         endif
      endif
      
      if(da_acc) then 
         vLn%link = vL%link;       vRn%link = vR%link
         vLn%plink = vL%plink;     vRn%plink = vR%plink 
         vLn%link(1) = iv2;        vLn%link(3) = vR%link(3) 
         vRn%link(1) = vL%link(1); vRn%link(3) = iv1  
         diagram%vertexList(vL%link(1))%link(3) = iv2
         diagram%vertexList(vR%link(3))%link(1) = iv1
         
         call copy_vtex(Nph, dmc_band, nbnd,nsvd,vLn,diagram%vertexList(iv1))
         call copy_vtex(Nph, dmc_band, nbnd,nsvd,vRn,diagram%vertexList(iv2))
         stat%accept(7) = stat%accept(7) + 1; call c1_mark_commit(7)
         n_accept = n_accept + 1
         if (p1_fixture_on() .and. diagram%vertexList(iv1)%tau <= diagram%vertexList(iv2)%tau) then
            call p1_prop_ell(diagram, iv1, iv2, ell_rev)
            call c2_dump_swap(4, iv1, iv2, diagram, ell_rev, 0.d0, &
              c1_ncall(c1_g_prod)+c1_ncall(c1_g_full)+c1_ncall(c1_epc)-ng_enter, &
              c1_ncall(c1_env_left)+c1_ncall(c1_env_right)+c1_ncall(c1_env_mid)-ne_enter, &
              minval(diagram%vertexList(iv1)%ekout), minval(diagram%vertexList(iv2)%ekout), 0.d0, 0.d0)
         else if (p1_fixture_on()) then
            call p1_prop_ell(diagram, iv2, iv1, ell_rev)
            call c2_dump_swap(4, iv2, iv1, diagram, ell_rev, 0.d0, &
              c1_ncall(c1_g_prod)+c1_ncall(c1_g_full)+c1_ncall(c1_epc)-ng_enter, &
              c1_ncall(c1_env_left)+c1_ncall(c1_env_right)+c1_ncall(c1_env_mid)-ne_enter, &
              minval(diagram%vertexList(iv2)%ekout), minval(diagram%vertexList(iv1)%ekout), 0.d0, 0.d0)
         endif
      else
         n_stage2 = n_stage2 + 1
         if (p1_fixture_on()) call c2_dump_swap(12, iv1, iv2, diagram, ell_hat, ellR, &
           c1_ncall(c1_g_prod)+c1_ncall(c1_g_full)+c1_ncall(c1_epc)-ng_enter, &
           c1_ncall(c1_env_left)+c1_ncall(c1_env_right)+c1_ncall(c1_env_mid)-ne_enter, &
           0.d0, 0.d0, 0.d0, 0.d0)
      end if 

   end subroutine
