! Solve_master

module solver

  use problem_class, only : problem_type
  use utils, only : pack, unpack
  use logger, only : log_screen

  implicit none
  private

  integer(kind=8), save :: iktotal

  public  :: solve, init_rk45

contains

!=====================================================================
! Master Solver
!
subroutine solve(pb)

  use output, only : write_output
  use my_mpi, only : is_MPI_parallel, finalize_mpi, synchronize_all

  type(problem_type), intent(inout)  :: pb

  ! Synchronise all processes
  if (is_MPI_parallel()) call synchronize_all()

  ! Before the first step, update field and write output (initial state)
  call update_field(pb)
  call write_output(pb)

  iktotal=0
  ! Time loop
  do while (pb%it /= pb%itstop)
    pb%it = pb%it + 1
    ! Do one integration step
    call do_bsstep(pb)
    ! Update field variables
    call update_field(pb)
    ! Check if we need to stop (set pb%itstop = pb%it)
    call check_stop(pb)
    ! Write output (if needed)
    call write_output(pb)
  enddo

  ! Finalise MPI
  if (is_MPI_parallel()) call finalize_mpi()

end subroutine solve



!=====================================================================
! pack, do bs_step and unpack
!
! IMPORTANT NOTE : between pack/unpack pb%v & pb%theta are not up-to-date
! SEISMIC IMPORTANT NOTE: when the CNS model is used, pb%tau is not up-to-date
!
subroutine do_bsstep(pb)

  use derivs_all
  use ode_bs
  use ode_rk45, only: rkf45_d
  use ode_rk45_2, only: rkf45_d2
  use constants, only: FID_SCREEN, SOLVER_TYPE
  use diffusion_solver, only: update_PT_final

  type(problem_type), intent(inout) :: pb

  double precision, dimension(pb%neqs*pb%mesh%nn) :: yt, dydt, yt_scale
  double precision, dimension(pb%neqs*pb%mesh%nn) :: yt_prev
  double precision, dimension(pb%mesh%nn) :: main_var
  integer :: ik, neqs
  character(255) :: msg

  neqs = pb%neqs * pb%mesh%nn

  ! SEISMIC: in the case of the CNS model, solve for tau and not v
  if (pb%i_rns_law == 3) then   ! SEISMIC: CNS model
    main_var = pb%tau
  else  ! SEISMIC: not CNS model (i.e. rate-and-state)
    main_var = pb%v
  endif

  call pack(yt, pb%theta, main_var, pb%sigma, pb%theta2, pb%slip, pb)
  yt_prev = yt

  ! SEISMIC: user-defined switch to use either (1) the Bulirsch-Stoer method, or
  ! the (2) Runge-Kutta-Fehlberg method
  if (SOLVER_TYPE == 0) then
    ! Default value of SOLVER_TYPE has not been altered
    call log_screen("The default solver type (0) has not been altered, and no solver was picked")
    call log_screen("Check the input script and define a solver type > 0")
    stop

  elseif (SOLVER_TYPE == 1) then
    ! Use Bulirsch-Stoer method

    ! this update of derivatives is only needed to set up the scaling (yt_scale)
    call derivs(pb%time,yt,dydt,pb)
    yt_scale=dabs(yt)+dabs(pb%dt_try*dydt)
    ! One step
    call bsstep(yt,dydt,neqs,pb%time,pb%dt_try,pb%acc,yt_scale,pb%dt_did,pb%dt_next,pb,ik)

    ! SEISMIC NOTE: what is happening here?
    if (pb%dt_max >  0.d0) then
      pb%dt_try = min(pb%dt_next,pb%dt_max)
    else
      pb%dt_try = pb%dt_next
    endif

  elseif (SOLVER_TYPE == 2) then
    ! Set-up Runge-Kutta solver

    pb%rk45%iflag = -2 ! Reset to one-step mode each call
    pb%t_prev = pb%time

    100 continue

    ! Call Runge-Kutta solver routine
    call rkf45_d( derivs_rk45, neqs, yt, pb%time, pb%tmax, &
                  pb%acc, 0d0, pb%rk45%iflag, pb%rk45%work, pb%rk45%iwork)

    ! Basic error checking. See description of rkf45_d in ode_rk45.f90 for details
    select case (pb%rk45%iflag)
    case (3)
      call log_screen("RK45 error [3]: relative error tolerance too small")
      call stop_simulation(pb)
    case (4)
      ! write(FID_SCREEN, *) "RK45 warning [4]: integration took more than 3000 derivative evaluations"
      yt = yt_prev
      goto 100
    case (5)
      call log_screen("RK45 error [5]: solution vanished, relative error test is not possible")
      call stop_simulation(pb)
    case (6)
      call log_screen("RK45 error [6]: requested accuracy could not be achieved")
      call stop_simulation(pb)
    case (8)
      call log_screen("RK45 error [8]: invalid input parameters")
      call stop_simulation(pb)
    end select

  elseif (SOLVER_TYPE == 3) then
    ! Set-up Runge-Kutta solver
    pb%t_prev = pb%time
    ! Call Runge-Kutta solver routine
    call rkf45_d2(derivs, yt, pb%time, pb%dt_max, pb%acc, 0d0, pb)
  else
    ! Unknown solver type
    write(msg, *) "Solver type", SOLVER_TYPE, "not recognised"
    call log_screen(msg)
    stop
  endif

  ! Set time step
  pb%dt_did = pb%time - pb%t_prev

  iktotal=ik+iktotal

  call unpack(yt, pb%theta, main_var, pb%sigma, pb%theta2, pb%slip, pb)
  if (pb%features%tp == 1) call update_PT_final(pb%dt_did, pb)

  ! SEISMIC: retrieve the solution for tau in the case of the CNS model, else
  ! retreive the solution for slip velocity
  if (pb%i_rns_law == 3) then
    pb%tau = main_var
  else
    pb%v = main_var
  endif

end subroutine do_bsstep


!=====================================================================
! Update field: slip, tau, potency potency rate, crack,
!
subroutine update_field(pb)

  use friction, only : friction_mu, dtheta_dt
  use friction_cns, only : compute_velocity
  use my_mpi, only: max_allproc, is_MPI_parallel
  use diffusion_solver, only : update_PT_final

  type(problem_type), intent(inout) :: pb

  double precision, dimension(pb%mesh%nn) :: P
  integer :: ivmax

  ! SEISMIC: obtain P at the previous time step
  P = 0d0
  if (pb%features%tp == 1) P = pb%P

  ! SEISMIC: in case of the CNS model, re-compute the slip velocity with
  ! the final value of tau, sigma, and porosity. Otherwise, use the standard
  ! rate-and-state expression to calculate tau as a function of velocity
  if (pb%i_rns_law == 3) then
    pb%v = compute_velocity(pb%tau, pb%sigma-P, pb%theta, pb%theta2, pb)
  else
    pb%tau = (pb%sigma-P) * friction_mu(pb%v,pb%theta,pb) + pb%coh
  endif

  ! Update pb%vmaxglob (required for stopping routine)
  ! Note that ivmax is re-computed globally at output time, so no need
  ! to store this quantity for now.
  ivmax = maxloc(pb%v, 1)
  if (is_MPI_parallel()) then
    call max_allproc(pb%v(ivmax), pb%vmaxglob)
  else
    pb%vmaxglob = pb%v(ivmax)
  endif

end subroutine update_field

!=====================================================================
! check stop:
!
subroutine check_stop(pb)

  use constants, only: FID_SCREEN
  use my_mpi, only: is_MPI_parallel, is_mpi_master, finalize_mpi

  type(problem_type), intent(inout) :: pb

  double precision, save :: vmax_old = 0d0, vmax_older = 0d0
  character(255) :: msg

  if (pb%itstop>0) return

  select case (pb%NSTOP)

   ! STOP if time > tmax
    case (0)
      if (pb%time >= pb%tmax) call stop_simulation(pb)

   ! STOP soon after end of slip localization
    case (1)
      call log_screen("Stop criterion 1 (end of slip localization) is deprecated")
      stop "Terminating..."

   ! STOP 10 ox snapshots after maximum slip rate
    case (2)
      if (pb%it > 2 .and. vmax_old > vmax_older .and. pb%vmaxglob < vmax_old)  &
        pb%itstop = pb%it+10*pb%ox%ntout
      vmax_older = vmax_old
      vmax_old = pb%vmaxglob

   ! STOP at a slip rate threshold (here tmax is threshold velocity)
    case (3)
      if (pb%vmaxglob > pb%tmax) call stop_simulation(pb)

    case default
      write(msg, *) "Stop criterion ", pb%NSTOP, " not implemented"
      call log_screen(msg)
      stop "Terminating..."

  end select

end subroutine check_stop

!=====================================================================
! A "soft" stop of the simulation by letting the solver loop run out
!
subroutine stop_simulation(pb)
  type(problem_type), intent(inout) :: pb

  ! Setting itstop to current iteration will terminate the solver loop
  pb%itstop = pb%it

end subroutine stop_simulation


subroutine init_rk45(pb)

  use problem_class
  use derivs_all
  use ode_rk45, only: rkf45_d
  use constants, only: FID_SCREEN

  type(problem_type), intent(inout) :: pb
  double precision, dimension(pb%neqs*pb%mesh%nn) :: yt
  double precision, dimension(pb%mesh%nn) :: main_var
  integer :: nwork

  call log_screen("Initialising RK45 solver")

  nwork = 3 + 6*pb%neqs*pb%mesh%nn
  pb%rk45%iflag = -1
  allocate(pb%rk45%work(nwork))
  allocate(pb%rk45%iwork(5))

  if (pb%i_rns_law == 3) then   ! SEISMIC: CNS model
    main_var = pb%tau
  else  ! SEISMIC: not CNS model (i.e. rate-and-state)
    main_var = pb%v
  endif

  call pack(yt, pb%theta, main_var, pb%sigma, pb%theta2, pb%slip, pb)

  call rkf45_d( derivs_rk45, pb%neqs*pb%mesh%nn, yt, pb%time, pb%time, &
                pb%acc, 0d0, pb%rk45%iflag, pb%rk45%work, pb%rk45%iwork)

  select case (pb%rk45%iflag)
  case (3)
    call log_screen("RK45 error [3]: relative error tolerance too small")
    stop
  case (4)
    ! call log_screen("RK45 warning [4]: integration took more than 3000 derivative evaluations")
  case (5)
    call log_screen("RK45 error [5]: solution vanished, relative error test is not possible")
    stop
  case (6)
    call log_screen("RK45 error [6]: requested accuracy could not be achieved")
    stop
  case (8)
    call log_screen("RK45 error [8]: invalid input parameters")
    stop
  end select

call log_screen("Finished initialising RK45 solver")

end subroutine init_rk45


end module solver
