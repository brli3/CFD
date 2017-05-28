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
! The objective of this program is to emphasise the importance of 
! spatial bias of the discretization scheme.
!
! For example, discrete equation for space index i and time level n
! using forward in time, centred in space (FTCS) scheme.
!
!   c(i,n+1) - c(i,n)        c(i+1,n)-c(i-1,n)
!   ----------------- + u * ------------------- = 0
!          dt                       2*dx     
!
! We will approximate dc/dx using schemes of 
! backward, forward and centred in space, combined with forward in time.
! i.e. FTBS, FTFS and FTCS
! And the results will be compared against each other.
! 
!
! FTFS fails because the direction of approximation is against that
! of the wave (information) propagation.
! FTCS gives oscillating results because it has no spatial variance.
! 
! Notice that only FTBS has a stable solution for this problem.
! However, the results show that it is dispersive.
! Also play with the number of grid and time step (or CFL number)
! for this scheme.
! 
! Author: Ruipengyu Li
! Modified: 17/04/2017
!
implicit none
integer, parameter :: dp =  selected_real_kind(15)
integer :: i, nx
real(dp) :: xa, xb, u, endtime, dt
real(dp), allocatable, dimension(:) :: x, c0, c_ftbs, c_ftfs, c_ftcs

! number of grid points
nx = 41
allocate(x(nx), c0(nx), c_ftbs(nx), c_ftfs(nx), c_ftcs(nx))
x = 0.0_dp; c0 = 0.0_dp; c_ftbs = 0.0_dp; c_ftfs = 0.0_dp; c_ftcs = 0.0_dp
! constant velocity
u = 1.0_dp
! domain coordinate 
xa = 0.0_dp; xb = 2.0_dp
! initialize
call initial(xa, xb, x, c0, nx)
! simulation time (start from zero)
endtime = 0.5_dp 
! time step
dt = 2.0e-2_dp

! FTBS, stable if CFL<1.0 but dispersive.
call fd1d_adv(u, c0, endtime, dt, x, c_ftbs, nx, 1)
! FTFS, results not bounded.
call fd1d_adv(u, c0, endtime, dt, x, c_ftfs, nx, 2)
! FTCS, oscillating solution.
call fd1d_adv(u, c0, endtime, dt, x, c_ftcs, nx, 3)

! print out x, initial profile and results from different schemes.
write(*,'(t8,5(a,5x),/,5(1x,f9.4))') &
      'x', 'c_init', 'c_ftbs', 'c_ftfs', 'cftcs', & 
      (x(i), c0(i), c_ftbs(i), c_ftfs(i), c_ftcs(i), i=1,nx)

deallocate(x, c0, c_ftbs, c_ftfs, c_ftcs)
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
subroutine fd1d_adv(u, c0, endtime, dt, x, c, nx, scheme)
!
! Solve 1-D advection equation:
!
!   dc/dt + u*dc/dx = 0
!
!   c: scalar variable of interest
!   u: velocity, assume constant 
!
! FTBS, FTFS and FTCS are used.
! For example, discrete equation for space index i and time level n
! using forward in time, centred in space (FTCS) scheme.
!
!   c(i,n+1) - c(i,n)        c(i+1,n)-c(i-1,n)
!   ----------------- + u * ------------------- = 0
!          dt                       2*dx     
!
! The explicit Euler method does not require solving a system
! of linear equations.
! 
implicit none
integer, parameter :: dp = selected_real_kind(15)
integer :: i, j, im1, ip1, ipr
integer, intent(in) :: nx
integer, intent(in) :: scheme
integer, external :: INDEX_WRAP
real(dp), intent(in) :: u
real(dp), intent(in) :: c0(nx)
real(dp), intent(in) :: x(nx)
real(dp), intent(in) :: endtime, dt
real(dp), dimension(nx), intent(out) :: c
real(dp), allocatable :: cold(:)
real(dp) :: dx, time, m 
allocate(cold(nx))
c = c0
cold = c0
dx = x(2) - x(1) ! assume uniform grid
m = u*dt/dx
write(*,'(/,a,f8.3,/)') 'CFL coef u*dt/dx: ', m
time = 0.0_dp
do while(time<=endtime)
  do i=1,nx
    ! periodic b.c., connect c(1) and c(nx).
    im1 = INDEX_WRAP(i-1, 1, nx)
    ip1 = INDEX_WRAP(i+1, 1, nx)
    ipr = INDEX_WRAP(i, 1, nx)
    select case(scheme)
    case(1) ! FTBS
      c(i) = m*cold(im1) + (1-m)*cold(ipr) 
      !c(i) = m*cold(i-1) + (1-m)*cold(i) 
    case(2)
      c(i) = (1+m)*cold(ipr) - m*cold(ip1) 
      !c(i) = (1+m)*cold(i) - m*cold(i+1) 
    case(3)
      c(i) = cold(i) + 0.5_dp*m * (cold(im1) - cold(ip1)) 
      !c(i) = cold(i) + 0.5_dp*m * (cold(i-1) - cold(i+1)) 
    end select
  end do
  cold = c
  time = time + dt
end do
deallocate(cold)
end subroutine fd1d_adv
!*******************************************************************************
function POSMOD(i, j)
! Returns positive modulo of i and j 
implicit none
integer, intent(in) :: i, j
integer :: POSMOD
POSMOD = MOD(i,j)
if(POSMOD<0) POSMOD = POSMOD + ABS(j)
end function POSMOD
!*******************************************************************************
function INDEX_WRAP(i, istart, iend)
! Remap i in the limit of istart and iend
implicit none
integer, intent(in) :: i, istart, iend
integer :: INDEX_WRAP
integer :: imin, imax
integer, external :: POSMOD
imin = MIN(istart, iend)
imax = MAX(istart, iend)
INDEX_WRAP = imin + POSMOD(i-imin, imax-imin+1)
end function INDEX_WRAP
!*******************************************************************************
