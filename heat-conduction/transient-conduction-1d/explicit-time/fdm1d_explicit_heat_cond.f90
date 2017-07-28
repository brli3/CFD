!*******************************************************************************
program main
!
! Solve 1-D heat equation:
!
!   rho*cp*dT/dt = k*d2T/dx2 + qdot(x)
!
!   rho: density
!   cp: specific heat
!   k: thermal conductivity
!      assume constant (one material) for simplicity.
!   qdot: heat generation
!   alpha = k/(rho*cp) : thermal diffusivity

! Discrete equation for space index i and time level n
! using finite difference method in space 
! and explicit method in time integration:
!
!   T(i,n+1) - T(i,n)     k      T(i+1,n)-2T(i,n)+T(i-1,n)    q(i,n)
!   ----------------- = ------ * -------------------------- + ------
!          dt           rho*cp            dx*dx               rho*cp
!
!   CFL = alpha*dt/(dx*dx) is the CFL coefficient.
!   CFL should not greater than 0.5 for stability.
!   This makes the middle coefficient 1-alpha*dt/(dx*dx) positive.
!   In transient heat transfer, Fo = alpha*dt/(dx*dx) is
!   also known as the grid Fourier number.
! 
! Author: Ruipengyu Li
! Modified: 16/04/2017
!
implicit none
integer, parameter :: dp =  selected_real_kind(15)
integer :: i, nx
real(dp) :: k, rho, cp, alpha
real(dp) :: xa, xb, ta, tb, t0, endtime, dt, dx
real(dp) :: cfl
real(dp), allocatable :: x(:), t(:)

! number of grid points
nx = 21
allocate(x(nx), t(nx))

! physical properties
k = 5.0_dp; rho = 1.0_dp; cp = 1.0_dp
alpha = k/rho/cp
! domain size
xa = 0.0_dp; xb = 1.0_dp
dx = (xb - xa) / REAL(nx-1, kind=dp)
! boundary temperatures
ta = 0.0_dp; tb = 0.0_dp
! initial temperature
t0 = 0.0_dp
! simulation time (start from zero)
endtime = 10.0_dp 

! CFL
! Try a value greater than 0.5 if you want to
! see how the solution blows up.
cfl = 0.5
! time step
dt = cfl*dx**2/alpha

write(*,'(/,a,f8.3,/)') 'CFL coef: ', cfl
write(*,'(a,es9.2,/)') 'Time step from CFL: ', dt

call fd1d_explicit_heat_cond(k, rho, cp, xa, xb, ta, tb,  &
                             t0, endtime, dt, x, t, nx)

write(*,'(a,8x,a,/,(3x,f9.2,1x,f10.3))') 'Solution: x', 'T', (x(i), t(i), i=1,nx)

deallocate(x, t)
stop
end program main
!*******************************************************************************
function QDOT(x, time)
! Heat generation term
! Set to zero if there is no generation.
implicit none
integer, parameter :: dp = selected_real_kind(15)
real(dp), intent(in) :: x, time
real(dp) :: QDOT
if(x<0.1) then
  QDOT = 1.0e3_dp
else
  QDOT = 0.0_dp
end if
end function QDOT
!*******************************************************************************
subroutine fd1d_explicit_heat_cond(k, rho, cp, xa, xb, ta, tb, &
                                   t0, endtime, dt, x, t, nx)
! Solve 1-D heat equation:
!
!   rho*cp*dT/dt = k*d2T/dx2 + qdot(x,t)
!
! Forward Euler approx for time
! Central difference scheme for second derivative in space
!
! Discrete equation for space index i and time level n:
!
!   T(i,n+1) - T(i,n)     k      T(i+1,n)-2T(i,n)+T(i-1,n)    q(i,n)
!   ----------------- = ------ * -------------------------- + ------
!          dt           rho*cp            dx*dx               rho*cp
!
!
! The explicit Euler method does not require solving a system
! of linear equations.
! 
implicit none
integer, parameter :: dp = selected_real_kind(15)
integer :: i, nt
integer, intent(in) :: nx
real(dp), intent(in) :: k, rho, cp
real(dp), intent(in) :: t0(nx)
real(dp), intent(in) :: xa, xb, ta, tb, endtime, dt
real(dp), dimension(nx), intent(out) :: x, t
real(dp), allocatable :: tcopy(:)
real(dp), external :: QDOT
real(dp) :: dx, time
allocate(tcopy(nx))
t = t0
tcopy = t0
dx = (xb - xa) / REAL(nx-1, kind=dp)
! set grid
x(1) = xa
do i=2,nx
  x(i) = x(1) + (i-1)*dx
end do
time = 0.0_dp
do while(time<=endtime)
  t(1) = ta
  t(nx) = tb
  do i=2,nx-1
    t(i) = tcopy(i) + (k/rho/cp * (tcopy(i+1) - 2*tcopy(i) + tcopy(i-1)) / &
                       dx**2 + QDOT(x(i),time)/rho/cp) * dt
  end do
tcopy = t
time = time + dt
end do
deallocate(tcopy)
end subroutine fd1d_explicit_heat_cond
!*******************************************************************************
