!*******************************************************************************
program main
! Solve 1-D steady state heat equation:
! d( k(x) dT/dx)dx + qdot(x) = 0
! using finite difference method
! with Dirichlet Boundary condition.
!
! Author: Ruipengyu Li
! Modified: 14/04/2017
implicit none
integer, parameter :: dp =  selected_real_kind(15)
integer :: i
integer :: n
real(dp) :: xa, xb, ta, tb
real(dp), allocatable :: x(:), t(:)
n = 21
allocate(x(n), t(n))
xa = 0.0_dp; xb = 1.0_dp
ta = 0.0_dp; tb = 100.0_dp
call fd1d_steady_heat_cond(xa, xb, ta, tb, x, t, n)
write(*,'(a,/,(2x,f9.2,1x,f10.3))') 'Solution:  x       T', (x(i), t(i), i=1,n)
deallocate(x, t)
stop
end program main
!*******************************************************************************
function K(x)
! Spatial function of thermal conductivity
! Set to one if not needed.
implicit none
integer, parameter :: dp = selected_real_kind(15)
real(dp), intent(in) :: x
real(dp) :: K
if(x < 0.5_dp) then
  k = 0.5_dp
else
  k = 1.0_dp
end if
end function K
!*******************************************************************************
function QDOT(x)
! Function of heat generation term
! Set to zero if there is no generation.
implicit none
integer, parameter :: dp = selected_real_kind(15)
real(dp), intent(in) :: x
real(dp) :: QDOT
QDOT = 0.0_dp
end function QDOT
!*******************************************************************************
subroutine fd1d_steady_heat_cond(xa, xb, ta, tb, x, t, n)
! The PDE is discretized to have n-2 unkowns and equations.
! If we also treat the 2 boundary temperatures as unknowns, 
! then we can write a system of n equations.
implicit none
integer, parameter :: dp = selected_real_kind(15)
integer :: i
integer, intent(in) :: n
real(dp), intent(in) :: xa, xb, ta, tb
real(dp), dimension(n), intent(out) :: x, t
real(dp), allocatable :: tri_a(:), tri_b(:), tri_c(:), rhs(:)
real(dp), external :: K, QDOT
real(dp) :: dx, xm, xp
allocate(tri_a(n), tri_b(n), tri_c(n), rhs(n))
tri_a = 0.0_dp; tri_b = 0.0_dp; tri_c = 0.0_dp
x = 0.0_dp; rhs =  0.0_dp; t = 0.0_dp
dx = (xb - xa) / real(n-1, kind=dp)
x(1) = xa
do i=2,n-1
  x(i) = xa + (i-1)*dx 
end do
x(n) = xb
!do i=1,n
!  x(i) = (REAL(n-i, kind=dp) * xa + REAL(i-1, kind=dp) * xb) &
!         / REAL(n-1, kind=dp)
!end do
tri_a(1) = 0.0_dp
tri_b(1) = 1.0_dp
tri_c(1) = 0.0_dp
rhs(1) = ta
do i=2,n-1
  xm = x(i) - 0.5_dp*dx
  xp = x(i) + 0.5_dp*dx
  tri_a(i) = K(xm) / dx**2
  tri_b(i) = -(K(xm) + K(xp)) / dx**2
  tri_c(i) = K(xp) / dx**2
  rhs(i) = -QDOT(x(i))
end do
tri_a(n) = 0.0_dp
tri_b(n) = 1.0_dp
tri_c(n) = 0.0_dp
rhs(n) = tb
call thomas(tri_a, tri_b, tri_c, rhs, t, n)
deallocate(tri_a, tri_b, tri_c, rhs)
end subroutine fd1d_steady_heat_cond
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
