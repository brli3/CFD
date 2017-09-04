!*******************************************************************************
program main
!
! This program solves 1-D wave equation:
!
!   d2u/dt2 + c*d2u/dx2 = 0
!
!   u: variable of interest
!   c: wave propagation velocity, assume constant 
!
! for illustration, we can consider a string with fixed ends. u is the 
! vertical displacement. An initial displacement is given.
!
! in domain (x1, x2) over the time interval (t1, t2)
!
! with boundary conditions: u(x1,t) and u(x2,t)
! and initial conditions: u(x,t1) = 0 and du/dt(x,t1) = 0
!
! Using central difference, discrete equation 
! for space index i and time level n reads:
!
!   u(i,n-1) - 2*u(i,n) + u(i,n+1)          u(i-1,n) - 2*u(i,n) + u(i+1,n)
!   ------------------------------- + c^2 * ------------------------------ = 0
!               dt^2                                   dx^2
!
! u(i,n-1) is not available at the first time step and will be obtained from
! the initial condition du/dt. At time step 0:
!
!   du         u(i,n+1) - u(i,n-1)
!  ----(i,n) = ------------------- 
!   dt                 2*dt
!
! Then u(i,n-1) is substituted in the discrete equation at the first time step.
!
!
! Author: Ruipengyu Li
! Modified: 18/04/2017
!
implicit none
integer, parameter :: dp =  selected_real_kind(15)
integer :: i, nx
real(dp), allocatable, dimension(:) :: x, u0, u
real(dp) :: x1, x2, dx, u1, u2, t1, t2
real(dp) :: cfl, dt, c

! number of grid points, check CFL
nx = 21
allocate(x(nx), u0(nx), u(nx))
x = 0.0_dp; u0 = 0.0_dp; u = 0.0_dp
! constant velocity
c = 1.0_dp
! domain coordinate 
x1 = 0.0_dp; x2 = 1.0_dp
! initialize
call initial(x1, x2, x, u0, nx)
dx = x(2) - x(1) ! assume uniform grid
! simulation time (start from zero)
t1 = 0.0_dp; t2 = 0.5_dp
cfl = 1.0_dp
! time step, check CFL
dt = cfl*dx/c
write(*,'(/,a,2(g10.3,1x),/)') 'CFL and Time step:  ', cfl, dt
! B.C.
u1 = 0.0_dp; u2 = 0.0_dp

call fd1d_cds_wave(c, u1, u2, u0, t1, t2, dt, x, u, nx)

! print out x, initial profile and numerical solution at the time t2.
write(*,'(/,a,/)') 'Finite difference analysis of 1-D wave equation.'
write(*,'(t8,3(a,6x),/,3(1x,f9.4))') &
      'x', 'u_init', 'u', & 
      (x(i), u0(i), u(i), i=1,nx)

deallocate(x, u0, u)
stop
end program main
!*******************************************************************************
subroutine initial(xa, xb, x, u0, nx)
! initial condition
implicit none
integer, parameter :: dp = selected_real_kind(15)
integer :: i
integer, intent(in) :: nx
real(dp), intent(in) :: xa, xb
real(dp), intent(out) :: x(nx)
real(dp), intent(out) :: u0(nx)
real(dp) :: dx, umax
! set grid
dx = (xb - xa) / REAL(nx-1, kind=dp)
x(1) = xa
do i=2,nx
  x(i) = x(1) + (i-1)*dx
end do
! initial displacement u(x,0)
umax = 0.05_dp
where(x < 0.7_dp)
  u0 = umax / 0.7_dp * x
else where
  u0 = umax / 0.3_dp * (1.0_dp-x)
end where 
! initial dudt(x1,t1) is 0.
end subroutine initial
!*******************************************************************************
subroutine fd1d_cds_wave(c, u1, u2, u0, t1, t2, dt, x, unp, nx)
!
! Solve 1-D wave equation:
!
!   d2u/dt2 + c*d2u/dx2 = 0
!
!   u: variable of interest
!   c: wave propagation velocity, assume constant 
!
! in domain (x1, x2) over the time interval (t1, t2)
!
! with boundary conditions: u(x1,t) and u(x2,t)
! and initial conditions: u(x,t1) = 0 and du/dt(x,t1) = 0
!
! Using central difference, discrete equation 
! for space index i and time level n reads:
!
!   u(i,n-1) - 2*u(i,n) + u(i,n+1)          u(i-1,n) - 2*u(i,n) + u(i+1,n)
!   ------------------------------- + c^2 * ------------------------------ = 0
!               dt^2                                   dx^2
!
! u(i,n-1) is not available at the first time step and will be obtained from
! the initial condition du/dt. At time step 0:
!
!   du         u(i,n+1) - u(i,n-1)
!  ----(i,n) = ------------------- = 0
!   dt                 2*dt
!
implicit none
integer, parameter :: dp = selected_real_kind(15)
integer :: i, j
integer, intent(in) :: nx
real(dp), intent(in) :: c, u1, u2
real(dp), intent(in) :: t1, t2, dt
real(dp), dimension(nx), intent(in) :: x, u0
real(dp), dimension(nx), intent(out) :: unp
real(dp), allocatable, dimension(:) :: unm, un
real(dp) :: dx, t, cfl 
allocate(unm(nx), un(nx))
unm = u0
un = u0
unp = u0
dx = x(2) - x(1) ! assume uniform grid
cfl = c*dt/dx
t = t1
! u at time level -1.
do i=2,nx-1
  ! from initial condition du/dt(x1,t1)=0
  unm(i) = un(i) + 0.5_dp * cfl**2 * (un(i-1) - 2.0_dp*un(i) + un(i+1))
end do
do while(t <= t2)
  t = t + dt
  do i=2,nx-1
    unp(i) = cfl**2 * (un(i-1) - 2.0_dp*un(i) + un(i+1)) + 2.0_dp*un(i) - unm(i)
  end do
  ! update boundary points
  unp(1) = u1
  unp(nx) = u2
  ! update data
  unm = un
  un = unp
end do
deallocate(unm, un)
end subroutine fd1d_cds_wave
!*******************************************************************************
