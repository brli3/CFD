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
! using central differencing in space 
! and implicit method in time integration:
!
!   T(i,n+1) - T(i,n)     k      T(i+1,n+1)-2T(i,n+1)+T(i-1,n+1)    q(i,n+1)
!   ----------------- = ------ * -------------------------------- + --------
!          dt           rho*cp            dx*dx                      rho*cp
!
!   The implicit method is unconditionally stable.
!   You can test it by trying different time step size.
! 
! Author: Ruipengyu Li
! Modified: 16/04/2017
!
implicit none
integer, parameter :: dp =  selected_real_kind(15)
integer :: i, nx
real(dp) :: k, rho, cp
real(dp) :: xa, xb, ta, tb, t0, endtime, dt, dx
real(dp), allocatable :: x(:), t(:)

! number of grid points
nx = 21
allocate(x(nx), t(nx))

! physical properties
k = 5.0_dp; rho = 1.0_dp; cp = 1.0_dp
! domain size
xa = 0.0_dp; xb = 1.0_dp
! boundary temperatures
ta = 0.0_dp; tb = 0.0_dp
! initial temperature
t0 = 0.0_dp
! simulation time (start from zero)
endtime = 10.0_dp 
! time step
dt = 2.0e-4_dp

call fd1d_implicit_heat_cond(k, rho, cp, xa, xb, ta, tb,  &
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
subroutine fd1d_implicit_heat_cond(k, rho, cp, xa, xb, ta, tb, &
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
!   T(i,n+1) - T(i,n)     k      T(i+1,n+1)-2T(i,n+1)+T(i-1,n+1)    q(i,n+1)
!   ----------------- = ------ * -------------------------------- + --------
!          dt           rho*cp            dx*dx                      rho*cp
!
! The implicit Euler method requires solving a system of linear equations.
! The coefficient matrix is tridiagonal.
! 
implicit none
integer, parameter :: dp = selected_real_kind(15)
integer :: i
integer, intent(in) :: nx
real(dp), intent(in) :: k, rho, cp
real(dp), intent(in) :: t0(nx)
real(dp), intent(in) :: xa, xb, ta, tb, endtime, dt
real(dp), dimension(nx), intent(out) :: x, t
real(dp), allocatable :: tcopy(:)
real(dp), allocatable, dimension(:) :: tri_a, tri_b, tri_c, rhs
real(dp), external :: QDOT
real(dp) :: dx, fo, time
allocate(tcopy(nx), tri_a(nx), tri_b(nx), tri_c(nx), rhs(nx))
tri_a = 0.0_dp; tri_b = 0.0_dp; tri_c = 0.0_dp; rhs = 0.0_dp
t = 0.0_dp; tcopy = 0.0_dp
t = t0
tcopy = t0
!write(*,*) t
!write(*,*) tcopy
dx = (xb - xa) / REAL(nx-1, kind=dp)
! set grid
x(1) = xa
do i=2,nx
  x(i) = x(1) + (i-1)*dx
end do
! grid Fourier number
fo = k * dt / rho / cp / dx**2
time = 0.0_dp
do while(time<=endtime)
  ! may be taken out the loop if the matrix is time-independent 
  t(1) = ta
  tri_a(1) = 0.0_dp
  tri_b(1) = 1.0_dp
  tri_c(1) = 0.0_dp
  rhs(1) = t(1)
  do i=2,nx-1
    tri_a(i) = -fo
    tri_b(i) = 1.0_dp + 2.0_dp*fo
    tri_c(i) = -fo
    rhs(i) = tcopy(i) + dt * QDOT(x(i), time + dt) / rho/cp
  end do
  t(nx) = tb
  tri_a(nx) = 0.0_dp
  tri_b(nx) = 1.0_dp
  tri_c(nx) = 0.0_dp
  rhs(nx) = t(nx)
  call thomas(tri_a, tri_b, tri_c, rhs, t, nx)
  tcopy = t
  time = time + dt
end do
deallocate(tcopy, tri_a, tri_b, tri_c, rhs)
end subroutine fd1d_implicit_heat_cond
!*******************************************************************************
subroutine thomas(a, b, c, d, x, n)
! Thomas matrix algorithm
! to solve a tridiagonal matrix system.
! a(i)x(i-1) + b(i)x(i) + c(i)x(i+1) = d(i)
! where a(1) = 0 and c(1) = 0.
! [b1  c1                      ] |x1|       |d1|
! [a2  b2  c2                  ] |x2|       |d2|
! [    a3  b3  c3              ] |x3|       |d3|
! [         .   .   .          ] |..|    =  |..|
! [             .   .   .      ] |..|       |..|
! [            an-1  bn-1  cn-1] |xn-1|     |dn-1|
! [                   an     bn] |xn|       |dn|
! 
! a: subdiagonal
! b: diagonal
! c: superdiagonal 
implicit none
integer, parameter :: dp = selected_real_kind(15)
integer :: i
integer, intent(in) :: n
real(dp), dimension(n), intent(in) :: a, b, c, d
real(dp), dimension(n), intent(out) :: x
real(dp), dimension(n) :: e, f
real(dp) :: temp
! forward sweep
e(1) = c(1) / b(1)
f(1) = d(1) / b(1)
do i=2,n
  temp = b(i) - a(i) * e(i-1)
  e(i) = c(i) / temp 
  f(i) = (d(i) - a(i) * f(i-1)) / temp
end do
! back substitution
x(n) = f(n)
do i=n-1,1,-1
  x(i) = f(i) - e(i) * x(i+1)
end do
end subroutine thomas
!*******************************************************************************
