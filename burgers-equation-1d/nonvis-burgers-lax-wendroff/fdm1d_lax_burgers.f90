!*******************************************************************************
program main
!
! This program computes 1-D non-viscous Burgers' equation:
!
!   du/dt + u*du/dx = 0
!
!   u: velocity 
! 
! using Lax-Wendroff scheme.
! The discrete form of equation may be found on references.
! 
! Author: Ruipengyu Li
! Modified: 17/04/2017
!
implicit none
integer, parameter :: dp =  selected_real_kind(15)
integer :: i, nx
real(dp) :: xa, xb, endtime, dt
real(dp), allocatable, dimension(:) :: x, u0, u

! number of grid points, make sure CFL<0.5
nx = 81
allocate(x(nx), u0(nx), u(nx))
x = 0.0_dp; u0 = 0.0_dp; u = 0.0_dp
! domain coordinate 
xa = 0.0_dp; xb = 1.0_dp
! initialize
call initial(xa, xb, x, u0, nx)
! simulation time (start from zero)
endtime = 0.5_dp 
! time step, make sure CFL<0.5
dt = 1.0e-2_dp

! Lax-Wendroff
call fd1d_lax_burg(u0, endtime, dt, x, u, nx)

! print out x, initial profile and results from different schemes.
write(*,'(t8,3(a,5x),/,3(1x,f9.4))') &
      'x', 'u_init', 'u_laxwendroff', & 
      (x(i), u0(i), u(i), i=1,nx)

deallocate(x, u0, u)
stop
end program main
!*******************************************************************************
subroutine initial(xa, xb, x, u_init, nx)
! initial condition
implicit none
integer, parameter :: dp = selected_real_kind(15)
integer :: i
integer, intent(in) :: nx
real(dp), intent(in) :: xa, xb
real(dp), intent(out) :: x(nx)
real(dp), intent(out) :: u_init(nx)
real(dp) :: dx
! set grid
dx = (xb - xa) / REAL(nx-1, kind=dp)
x(1) = xa
do i=2,nx
  x(i) = x(1) + (i-1)*dx
end do
! initialize a square wave
where(x<=0.5_dp)
  u_init = 1.0_dp
else where
  u_init = 0.0_dp
end where 
end subroutine initial 
!*******************************************************************************
subroutine fd1d_lax_burg(u0, endtime, dt, x, u, nx)
!
! Solve 1-D non-viscous Burgers' equation:
!
!   du/dt + u*du/dx = 0
!
! using Lax-Wendroff scheme.
!
! Discrete equation may be found in references
! 
implicit none
integer, parameter :: dp = selected_real_kind(15)
integer :: i
integer, intent(in) :: nx
real(dp), intent(in) :: u0(nx)
real(dp), intent(in) :: x(nx)
real(dp), intent(in) :: endtime, dt
real(dp), dimension(nx), intent(out) :: u
real(dp), allocatable :: uo(:)
real(dp) :: dx, time, m1, m2
allocate(uo(nx))
u = u0
uo = u0
dx = x(2) - x(1) ! assume uniform grid
m1 = 0.5_dp*dt/dx
write(*,'(a,f6.3)') 'dt/dx:', dt/dx
m2 = 0.5_dp * (dt/dx)**2
time = 0.0_dp
do while(time<=endtime)
  do i=2,nx-1
    u(i) = uo(i) - m1*(0.5_dp*uo(i+1)**2-0.5_dp*uo(i-1)**2) + &
           m2*(0.5_dp*(uo(i)+uo(i+1))*(0.5_dp*u(i+1)**2-0.5_dp*u(i)**2) - &
               0.5_dp*(u(i)+u(i-1))*(0.5*u(i)**2-0.5_dp*u(i-1)**2))
  end do
  u(1) = uo(1)
  u(nx) = uo(nx)
  uo = u
  time = time + dt
end do
deallocate(uo)
end subroutine fd1d_lax_burg
!*******************************************************************************
