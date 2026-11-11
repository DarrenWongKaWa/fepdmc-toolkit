! R5-P0 RNG backend shim.
! Official FEP-DMC OpenMP streams call Intel MKL VSL (random_tool.f90).
! This file is used only when MKL headers/libraries are absent.
! It supplies the same Fortran symbols with a xorshift64 uniform generator.
! It does not change Metropolis-Hastings ratios or proposal densities.
! Status: documented compiler/RNG deviation, not a physics rewrite.

module mkl_vsl_type
   implicit none
   type VSL_STREAM_STATE
      integer*8 :: state = 1
      integer :: seed0 = 1
      integer :: brng = 1
   end type VSL_STREAM_STATE
end module mkl_vsl_type

module mkl_vsl
   use mkl_vsl_type
   implicit none
   integer, parameter :: VSL_BRNG_MT19937 = 1
contains

   integer function vslnewstream(stream, brng, seed)
      type(VSL_STREAM_STATE), intent(out) :: stream
      integer, intent(in) :: brng, seed
      integer*8 :: s
      s = int(seed, 8)
      if (s <= 0) s = 1
      stream%state = s
      stream%seed0 = seed
      stream%brng = brng
      vslnewstream = 0
   end function vslnewstream

   integer function vdrnguniform(method, stream, n, r, a, b)
      integer, intent(in) :: method, n
      type(VSL_STREAM_STATE), intent(inout) :: stream
      real*8, intent(out) :: r(*)
      real*8, intent(in) :: a, b
      integer :: i
      real*8 :: u, span
      integer*8 :: x
      integer*8, parameter :: mask52 = 4503599627370495_8
      real*8, parameter :: den52 = 4503599627370496.0d0
      span = b - a
      x = stream%state
      do i = 1, n
         x = ieor(x, ishft(x, 13))
         x = ieor(x, ishft(x, -7))
         x = ieor(x, ishft(x, 17))
         if (x == 0) x = 1
         u = dble(iand(x, mask52)) / den52
         r(i) = a + span * u
      end do
      stream%state = x
      vdrnguniform = 0
      if (method < 0) vdrnguniform = 0
   end function vdrnguniform

end module mkl_vsl
