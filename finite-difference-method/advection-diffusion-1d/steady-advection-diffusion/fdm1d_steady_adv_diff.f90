!*******************************************************************************
program main
! Solve 1-D steady advection-diffusion equation:
!
!   u*dc/dx = D*d2c/dt2
!
!   u: velocity
!   c: scalar variable of interest
!   D: diffusivity (diffusion coefficient)
!  
! using finite difference method
! with boundary condition:
!   c(0) = 0. and c(1) = 1.
! 
! Numerical results are compared with the exact solution.
!
! Author: Ruipengyu Li
! Modified: 14/04/2017
implicit none
integer, parameter :: dp =  selected_real_kind(15)
integer :: i
integer :: n
real(dp) :: xa, xb, ca, cb, u, d, pe
real(dp), allocatable :: x(:), c(:), c_exact(:)
n = 31
allocate(x(n), c(n), c_exact(n))
! domain
xa = 0.0_dp; xb = 1.0_dp
! boundary values
ca = 0.0_dp; cb = 1.0_dp
! velocity and diffusivity
u = 1.0_dp; d = 0.05_dp
! Peclet number
pe = u*(xb-xa) / d

call fd1d_steady_adv_diff(xa, xb, ca, cb, u, d, x, c, n)

! exact solution
c_exact = (EXP(pe*x) - 1.0_dp) / (EXP(pe) - 1.0_dp)

write(*,'(a,5x,a,5x,a,/,(t6,f7.2,1x,f10.4,1x,f10.4))') &
     'Solution:  x', 'c', 'c exact', (x(i), c(i), c_exact(i), i=1,n)

deallocate(x, c, c_exact)
stop
end program main
!*******************************************************************************
subroutine fd1d_steady_adv_diff(xa, xb, ca, cb, u, d, x, c, n)
! Solves 1-D advection-diffusion equation:
!
!   u*dc/dx = D*d2c/dt2
!
! Central difference for first and second derivatives.
!
! The ODE is discretized to have n-2 unkowns and equations.
! If we also treat the 2 boundary temperatures as unknowns, 
! then we can write a system of n equations.
implicit none
integer, parameter :: dp = selected_real_kind(15)
integer :: i
integer, intent(in) :: n
real(dp), intent(in) :: xa, xb, ca, cb, u, d
real(dp), dimension(n), intent(out) :: x, c
real(dp), allocatable :: tri_a(:), tri_b(:), tri_c(:), rhs(:)
real(dp) :: dx
allocate(tri_a(n), tri_b(n), tri_c(n), rhs(n))
tri_a = 0.0_dp; tri_b = 0.0_dp; tri_c = 0.0_dp
x = 0.0_dp; rhs =  0.0_dp; c = 0.0_dp
dx = (xb - xa) / REAL(n-1, kind=dp)
x(1) = xa
do i=2,n-1
  x(i) = xa + (i-1)*dx 
end do
x(n) = xb
tri_a(1) = 0.0_dp
tri_b(1) = 1.0_dp
tri_c(1) = 0.0_dp
rhs(1) = ca
do i=2,n-1
  tri_a(i) = -(0.5_dp*u + d/dx)
  tri_b(i) = 2.0_dp*d / dx
  tri_c(i) = 0.5_dp*u - d/dx
  rhs(i) = 0.0_dp
end do
tri_a(n) = 0.0_dp
tri_b(n) = 1.0_dp
tri_c(n) = 0.0_dp
rhs(n) = cb

call thomas(tri_a, tri_b, tri_c, rhs, c, n)

deallocate(tri_a, tri_b, tri_c, rhs)
end subroutine fd1d_steady_adv_diff
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
