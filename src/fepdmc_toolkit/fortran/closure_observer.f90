! C1 nested exclusive observer. Does not consume diagram%seed.
! Mode C1_PROFILE=off|coarse|detail is runtime-only; physics path is unchanged.
module closure_observer
  implicit none
  integer, parameter :: c1_prep=1, c1_dispatch=2, c1_meas=3, c1_io_native=4
  integer, parameter :: c1_io_obs=5, c1_other=6, c1_upd1=7, c1_upd7=13
  integer, parameter :: c1_env_left=14, c1_env_right=15, c1_env_mid=16
  integer, parameter :: c1_g_prod=17, c1_g_full=18, c1_epc=19
  integer, parameter :: c1_cal_ek=20, c1_cal_wq=21, c1_prop=22, c1_hprop=23, c1_debug=24
  integer, parameter :: c1_nbucket=24
  integer, parameter :: c1_dp=kind(1.0d0)
  integer, save :: c1_level=0
  logical, save :: c1_inited=.false., c1_event_on=.true.
  integer, save :: c1_update=0, c1_event=0
  integer*8, save :: c1_rate=1, c1_sparse_state=1
  integer, save :: c1_stack_n=0, c1_stack_id(32)
  integer*8, save :: c1_stack_t0(32), c1_ticks(c1_nbucket), c1_incl(c1_nbucket)
  integer*8, save :: c1_count_sel(7), c1_count_valid(7), c1_count_commit(7)
  integer*8, save :: c1_hook=0, c1_eligible=0, c1_score=0, c1_exact=0
  integer*8, save :: c1_n_attempt=0, c1_n_meas=0, c1_n_detail_on=0, c1_n_detail_off=0
  integer*8, save, public :: c1_ncall(c1_nbucket)
  character(len=16), save :: c1_mode='off'
  logical, save :: c1_sparse=.false.
contains
  subroutine c1_init()
    character(len=32) :: env
    integer :: st, nlen
    if (c1_inited) return
    c1_inited=.true.
    c1_level=0
    c1_mode='off'
    call get_environment_variable('C1_PROFILE', env, length=nlen, status=st)
    if (st==0 .and. nlen>0) then
      c1_mode=trim(adjustl(env(1:nlen)))
      if (trim(c1_mode)=='coarse') c1_level=1
      if (trim(c1_mode)=='detail' .or. trim(c1_mode)=='detailed') c1_level=2
    endif
    call get_environment_variable('C1_SPARSE', env, length=nlen, status=st)
    c1_sparse=.false.
    if (st==0 .and. nlen>0) then
      if (env(1:nlen)=='1' .or. env(1:nlen)=='on') c1_sparse=.true.
    endif
    c1_sparse_state=9143311_8
    call get_environment_variable('C1_SPARSE_SEED', env, length=nlen, status=st)
    if (st==0 .and. nlen>0) read(env(1:nlen),*) c1_sparse_state
    if (c1_sparse_state==0_8) c1_sparse_state=1_8
    call system_clock(count_rate=c1_rate)
    if (c1_rate<=0) c1_rate=1
    c1_ticks=0; c1_incl=0; c1_stack_n=0; c1_ncall=0
    c1_count_sel=0; c1_count_valid=0; c1_count_commit=0
    c1_hook=0; c1_eligible=0; c1_score=0; c1_exact=0
    c1_n_attempt=0; c1_n_meas=0; c1_n_detail_on=0; c1_n_detail_off=0
    c1_event_on=.true.
  end subroutine

  logical function c1_active(id)
    integer, intent(in) :: id
    c1_active=.false.
    if (c1_level<=0) return
    if (c1_level==1 .and. id>=c1_env_left) return
    if (c1_level==2 .and. id>=c1_env_left .and. .not.c1_event_on) return
    c1_active=.true.
  end function

  subroutine c1_tic(id)
    integer, intent(in) :: id
    integer :: n
    integer*8 :: t
    if (.not.c1_inited) call c1_init()
    c1_ncall(id)=c1_ncall(id)+1
    if (.not.c1_active(id)) return
    call system_clock(count=t)
    n=c1_stack_n
    if (n>=1) c1_ticks(c1_stack_id(n))=c1_ticks(c1_stack_id(n))+(t-c1_stack_t0(n))
    n=n+1
    if (n>32) error stop 'C1 timer stack overflow'
    c1_stack_n=n
    c1_stack_id(n)=id
    c1_stack_t0(n)=t
    c1_incl(id)=c1_incl(id)-t
  end subroutine

  subroutine c1_toc(id)
    integer, intent(in) :: id
    integer :: n
    integer*8 :: t
    if (.not.c1_active(id)) return
    n=c1_stack_n
    if (n<1) return
    if (c1_stack_id(n)/=id) error stop 'C1 toc bucket mismatch'
    call system_clock(count=t)
    c1_ticks(id)=c1_ticks(id)+(t-c1_stack_t0(n))
    c1_incl(id)=c1_incl(id)+t
    c1_stack_n=n-1
    if (n-1>=1) c1_stack_t0(n-1)=t
  end subroutine

  subroutine c1_begin_attempt(attempt, update)
    integer, intent(in) :: attempt, update
    integer*8 :: x
    real(c1_dp) :: u
    if (.not.c1_inited) call c1_init()
    c1_n_attempt=c1_n_attempt+1
    c1_event=attempt
    c1_update=update
    if (update>=1 .and. update<=7) c1_count_sel(update)=c1_count_sel(update)+1
    c1_event_on=.true.
    if (c1_level==2 .and. c1_sparse) then
      x=c1_sparse_state
      x=ieor(x, ishft(x,13)); x=ieor(x, ishft(x,-7)); x=ieor(x, ishft(x,17))
      if (x==0_8) x=1_8
      c1_sparse_state=x
      u=real(iand(x, (2_8**52-1_8)), c1_dp)/real(2_8**52, c1_dp)
      c1_event_on=(u < 1.0_c1_dp/16.0_c1_dp)
      if (c1_event_on) then
        c1_n_detail_on=c1_n_detail_on+1
      else
        c1_n_detail_off=c1_n_detail_off+1
      endif
    endif
  end subroutine

  subroutine c1_mark_valid(u)
    integer, intent(in) :: u
    if (u>=1 .and. u<=7) c1_count_valid(u)=c1_count_valid(u)+1
  end subroutine

  subroutine c1_mark_commit(u)
    integer, intent(in) :: u
    if (u>=1 .and. u<=7) c1_count_commit(u)=c1_count_commit(u)+1
  end subroutine

  subroutine c1_mark_hook()
    c1_hook=c1_hook+1
  end subroutine
  subroutine c1_mark_eligible()
    c1_eligible=c1_eligible+1
  end subroutine
  subroutine c1_mark_score()
    c1_score=c1_score+1
  end subroutine
  subroutine c1_mark_exact()
    c1_exact=c1_exact+1
  end subroutine
  subroutine c1_mark_meas()
    c1_n_meas=c1_n_meas+1
  end subroutine

  subroutine c1_write()
    integer :: u, ios
    real(c1_dp) :: ex, inc, scale
    character(len=24) :: names(c1_nbucket)
    if (.not.c1_inited .or. c1_level<=0) return
    names = [character(len=24) :: 'prep','dispatch','measurement','io_native','io_observer','other', &
           'update_1','update_2','update_3','update_4','update_5','update_6','update_7', &
           'env_left','env_right','env_mid','g_prod','g_full','epc','cal_ek','cal_wq', &
           'propagator','hpropagator','debug']
    scale=1.0_c1_dp/real(c1_rate,c1_dp)
    open(newunit=u, file='c1_timers.csv', status='replace', action='write', iostat=ios)
    if (ios/=0) return
    write(u,'(a)') 'mode,bucket,exclusive_s,inclusive_s,attempts,measurements,sparse,detail_on,detail_off,note'
    do ios=1,c1_nbucket
      ex=real(c1_ticks(ios),c1_dp)*scale
      inc=real(c1_incl(ios),c1_dp)*scale
      if (c1_level<2 .and. ios>=c1_env_left) then
        write(u,'(a,",",a,",NOT_COMPUTED,NOT_COMPUTED,",i0,",",i0,",",l1,",",i0,",",i0,",uninstrumented_or_mode")') &
          trim(c1_mode), trim(names(ios)), c1_n_attempt, c1_n_meas, c1_sparse, c1_n_detail_on, c1_n_detail_off
      else
        write(u,'(a,",",a,",",es16.8e3,",",es16.8e3,",",i0,",",i0,",",l1,",",i0,",",i0,",exclusive_nested")') &
          trim(c1_mode), trim(names(ios)), ex, inc, c1_n_attempt, c1_n_meas, c1_sparse, c1_n_detail_on, c1_n_detail_off
      endif
    enddo
    close(u)
    open(newunit=u, file='c1_counts.csv', status='replace', action='write')
    write(u,'(a)') 'update,selected,native_valid,commit,hook,eligible,score,exact,attempts,measurements'
    do ios=1,7
      write(u,'(7(i0,","),i0,",",i0,",",i0)') ios, c1_count_sel(ios), c1_count_valid(ios), c1_count_commit(ios), &
        merge(c1_hook,0_8,ios==7), merge(c1_eligible,0_8,ios==7), merge(c1_score,0_8,ios==7), merge(c1_exact,0_8,ios==7), &
        c1_n_attempt, c1_n_meas
    enddo
    close(u)
  end subroutine
end module
