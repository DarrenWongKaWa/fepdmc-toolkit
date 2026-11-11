! Fixed Task-5 primitive test driver. Not a generated kernel.
program nk_prim_driver
  use, intrinsic :: iso_fortran_env, only: real64, int64
  use native_kernel_runtime_v1
  implicit none
  integer(int64) :: opcode, n, status, fail_code, kind_sym
  real(real64) :: k_eff, t, dtau, delta, gap, omega, g, L, imag_atol, floor, u, logw
  real(real64) :: sre, sim
  complex(real64) :: z, s, a(2, 2), b(2, 2), c(2, 2)
  character(len=32) :: fail_name

  read(*,*) opcode
  call nk_reset_counters()
  status = NK_OK
  select case (opcode)
  case (1_int64)
    read(*,*) k_eff, t, dtau
    z = nk_electron_scalar(k_eff, t, dtau)
    call write_z(z)
  case (2_int64)
    read(*,*) k_eff, t, dtau, delta, gap
    call nk_electron_twoband(k_eff, t, dtau, delta, gap, c)
    call write_m(c)
  case (3_int64)
    read(*,*) omega, dtau
    z = nk_phonon(omega, dtau)
    call write_z(z)
  case (4_int64)
    read(*,*) g, L, n
    z = nk_prefactor(g, L, n)
    call write_z(z)
  case (5_int64)
    read(*,*) g, L
    call nk_vertex_sigmaz(g, L, c)
    call write_m(c)
  case (6_int64)
    call read_m(a)
    call read_m(b)
    call nk_matmul(a, b, c)
    call write_m(c)
  case (7_int64)
    call read_m(a)
    read(*,*) sre, sim
    s = cmplx(sre, sim, kind=real64)
    call nk_scale_m2(a, s, c)
    call write_m(c)
  case (8_int64)
    call read_m(a)
    call nk_trace_counted(a, z)
    call write_z(z)
  case (9_int64)
    call read_m(a)
    call read_m(b)
    call nk_add_m2(a, b, c)
    call write_m(c)
  case (10_int64)
    read(*,*) sre, sim, imag_atol
    read(*,*) fail_name
    z = cmplx(sre, sim, kind=real64)
    fail_code = NK_INVALID_EXACT_TARGET
    if (trim(fail_name) == "INVALID_CHEAP_WEIGHT") fail_code = NK_INVALID_CHEAP_WEIGHT
    call nk_log_positive_real(z, imag_atol, logw, status, fail_code)
    write(*,'(I0,1X,ES24.16E3)') status, logw
  case (11_int64)
    read(*,*) u, floor
    call nk_log_uniform(u, floor, logw, status)
    write(*,'(I0,1X,ES24.16E3)') status, logw
  case (12_int64)
    read(*,*) u, kind_sym
    call nk_check_log_q(u, kind_sym, logw, status)
    write(*,'(I0,1X,ES24.16E3)') status, logw
  case (13_int64)
    read(*,*) k_eff, t
    write(*,'(ES24.16E3)') nk_xi(k_eff, t)
  case default
    error stop "unknown primitive opcode"
  end select
  write(*,'(I0,1X,I0,1X,I0,1X,I0)') nk_n_electron_twoband, nk_n_vertex_sigmaz, nk_n_matmul, nk_n_trace

contains

  subroutine read_m(m)
    complex(real64), intent(out) :: m(2, 2)
    real(real64) :: r11, i11, r12, i12, r21, i21, r22, i22
    read(*,*) r11, i11, r12, i12, r21, i21, r22, i22
    m(1, 1) = cmplx(r11, i11, kind=real64)
    m(1, 2) = cmplx(r12, i12, kind=real64)
    m(2, 1) = cmplx(r21, i21, kind=real64)
    m(2, 2) = cmplx(r22, i22, kind=real64)
  end subroutine read_m

  subroutine write_z(z)
    complex(real64), intent(in) :: z
    write(*,'(ES24.16E3,1X,ES24.16E3)') real(z, real64), aimag(z)
  end subroutine write_z

  subroutine write_m(m)
    complex(real64), intent(in) :: m(2, 2)
    write(*,'(ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3)') &
      real(m(1, 1), real64), aimag(m(1, 1)), &
      real(m(1, 2), real64), aimag(m(1, 2)), &
      real(m(2, 1), real64), aimag(m(2, 1)), &
      real(m(2, 2), real64), aimag(m(2, 2))
  end subroutine write_m

end program nk_prim_driver
