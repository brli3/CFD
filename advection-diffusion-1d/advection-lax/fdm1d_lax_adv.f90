!*******************************************************************************
program main
!
! Solve 1-D advection equation:
!
!   dc/dt + u*dc/dx = 0
!
!   c: scalar variable of interest
!   u: velocity, assume constant 
! with periodic boundary condition and an initial square wave.
!
! The objective of this program is to show the impoved stablility of 
! the Lax (Friedrichs) method compared with forward in time,
! centred in space (FTCS) scheme.
! Comparison is also made between the Lax and Lax-Wendroff methods.
!
! Discrete equation for space index i and time level n of the Lax method reads:
!
!   c(i,n+1) - 0.5*c(i-1,n) - 0.5*c(i+1,n)        c(i+1,n) - c(i-1,n)
!   -------------------------------------- + u * --------------------- = 0
!                   dt                                   2*dx     
!
! We can show that Lax-Friedrichs scheme stablises the FTCS scheme due to the 
! artificial viscosity, 1/2 term.
! Lax-Wendroff scheme is less dissipative but may suffer from oscillation.
! 
! Author: Ruipengyu Li
! Modified: 17/04/2017
!
implicit none
integer, parameter :: dp =  selected_real_kind(15)
integer :: i, nx
real(dp) :: xa, xb, u, endtime, dt
real(dp), allocatable, dimension(:) :: x, c0, c_lax, c_lax_wen

! number of grid points, make sure CFL<0.5
nx = 81
allocate(x(nx), c0(nx), c_lax(nx), c_lax_wen(nx))
x = 0.0_dp; c0 = 0.0_dp; c_lax = 0.0_dp
! constant velocity
u = 1.0_dp
! domain coordinate 
xa = 0.0_dp; xb = 2.0_dp
! initialize
call initial(xa, xb, x, c0, nx)
! simulation time (start from zero)
endtime = 0.5_dp 
! time step, make sure CFL<0.5
dt = 1.0e-2_dp

! Lax
call fd1d_lax_adv(u, c0, endtime, dt, x, c_lax, nx, 1)
! Lax-Wendroff
call fd1d_lax_adv(u, c0, endtime, dt, x, c_lax_wen, nx, 2)

! print out x, initial profile and results from different schemes.
write(*,'(t8,4(a,5x),/,4(1x,f9.4))') &
      'x', 'c_init', 'c_lax', 'c_lax_wendroff', & 
      (x(i), c0(i), c_lax(i), c_lax_wen(i), i=1,nx)

deallocate(x, c0, c_lax, c_lax_wen)
stop
end program main
!*******************************************************************************
subroutine initial(xa, xb, x, c_init, nx)
! initial condition
implicit none
integer, parameter :: dp = selected_real_kind(15)
integer :: i
integer, intent(in) :: nx
real(dp), intent(in) :: xa, xb
real(dp), intent(out) :: x(nx)
real(dp), intent(out) :: c_init(nx)
real(dp) :: dx
! set grid
dx = (xb - xa) / REAL(nx-1, kind=dp)
x(1) = xa
do i=2,nx
  x(i) = x(1) + (i-1)*dx
end do
! initialize a square wave
where(x>=0.5_dp .and. x<=1.0_dp)
  c_init = 2.0_dp
else where
  c_init = 1.0_dp
end where 
end subroutine initial 
!*******************************************************************************
subroutine fd1d_lax_adv(u, c0, endtime, dt, x, c, nx, scheme)
!
! Solve 1-D advection equation:
!
!   dc/dt + u*dc/dx = 0
!
!   c: scalar variable of interest
!   u: velocity, assume constant 
!
! Discrete equation for space index i and time level n of the Lax method reads:
!
!   c(i,n+1) - 0.5*c(i-1,n) - 0.5*c(i+1,n)        c(i+1,n) - c(i-1,n)
!   -------------------------------------- + u * --------------------- = 0
!                   dt                                   2*dx     
! 
implicit none
integer, parameter :: dp = selected_real_kind(15)
integer :: i
integer, intent(in) :: scheme
integer, intent(in) :: nx
real(dp), intent(in) :: u
real(dp), intent(in) :: c0(nx)
real(dp), intent(in) :: x(nx)
real(dp), intent(in) :: endtime, dt
real(dp), dimension(nx), intent(out) :: c
real(dp), allocatable :: cold(:)
real(dp) :: dx, time, m1, m2
allocate(cold(nx))
c = c0
cold = c0
dx = x(2) - x(1) ! assume uniform grid
m1 = 0.5_dp*u*dt/dx
m2 = 0.5_dp * (u*dt/dx)**2
time = 0.0_dp
do while(time<=endtime)
  select case(scheme)
  case(1) ! Lax
    ! periodic b.c., connect c(1) and c(nx).
    c(1) = (m1+0.5_dp)*cold(nx) + (0.5_dp-m1)*cold(2)
    do i=2,nx-1
      c(i) = (m1+0.5_dp)*cold(i-1) + (0.5_dp-m1)*cold(i+1) 
    end do
    ! periodic b.c., connect c(1) and c(nx).
    c(nx) = (m1+0.5_dp)*cold(nx-1) + (0.5_dp-m1)*cold(1)
  case(2) ! Lax-Wendroff
    c(1) = (m1+m2)*cold(nx) + (1.0_dp-2.0_dp*m2)*c(1) + (m2-m1)*cold(2)
    do i=2,nx-1
      c(i) = (m1+m2)*cold(i-1) + (1.0_dp-2.0_dp*m2)*c(i) + (m2-m1)*cold(i+1)
    end do
    c(nx) = (m1+m2)*cold(nx-1) + (1.0_dp-2.0_dp*m2)*c(nx) + (m2-m1)*cold(1)
  end select
  cold = c
  time = time + dt
end do
deallocate(cold)
end subroutine fd1d_lax_adv
!*******************************************************************************
