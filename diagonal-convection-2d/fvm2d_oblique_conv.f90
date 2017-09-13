!*******************************************************************************
program main 
!--------------------------------------------------------------------------
! This program solves the two-dimensional convection
! equation for a general scalar phi in uniform flow
!
! using finite volume method.
!
! Boundary conditions:
!   prescribed phi at west and south boundaries.
!   outflow conditions or prescribed phi at north and east boundaries.
! 
!                 outlet or phi=1
!          |--------------------------|
!          |                          |
!          |                          |
!          |                          |
!          |                          |
!          |                          |
!    phi=1 |                          | outlet or phi=0
!          |                          |
!          |                          |
!          |                          |
!          |                          |
!          |                          |
!          |                          |
!          ----------------------------
!                     phi=0
!   
! with Cartesian grids and uniform velocity field.
!
! This program demostrates the false(numerical) diffusion caused by UDS.
!
! If CDS is used, the zero values for the main diagonal makes the
! iterative solver fail.
!
! Author: Ruipengyu Li 
! Modified: 05/05/2017
!
! Reference:
!   J. H. Ferziger and M. Peric, Computational Methods for Fluid Dynamics,
!   3rd ed. Springer Berlin Heidelberg, 2001.
!--------------------------------------------------------------------------
implicit none

call fvm2d_conv_oblique()

stop

end program main

!*******************************************************************************
subroutine fvm2d_conv_oblique()
!
! Solves 2D scalar convection equation using finite volume method.
!
implicit none
integer, parameter :: dp = selected_real_kind(15)
character(len=80) :: filename1
integer :: i, j
integer :: ni, nj, nim1, njm1
integer :: nicv, njcv  ! no. of cell centres
integer :: ierr
integer :: east_north_bc
real(dp) :: xmin, xmax, ymin, ymax
real(dp) :: expfx, expfy  ! grid expansion factor
real(dp) :: vel
real(dp) :: dx, dy
real(dp), allocatable, dimension(:) :: x, y, xc, yc
real(dp), allocatable, dimension(:) :: fx, fy  ! interpolation coef
real(dp), allocatable, dimension(:,:) :: u, v  ! velocity at cell face
real(dp), allocatable, dimension(:,:) :: phi
real(dp), allocatable, dimension(:,:) :: ae, aw, an, as, ap, su
real(dp) :: ue, uw, vn, vs
real(dp) :: ge, gw, gn, gs
real(dp) :: ce, cw, cn, cs

filename1 = "velocity_phi.dat"

xmin = 0.0_dp
xmax = 1.0_dp
ymin = 0.0_dp
ymax = 1.0_dp
expfx = 1.0_dp
expfy = 1.0_dp
nicv = 20  ! no. of control volumes
njcv = 20  ! no. of control volumes

vel = 2.0_dp

ni = nicv + 2  ! no. of cell centres
nim1 = ni - 1  ! no. of cell faces
nj = njcv + 2
njm1 = nj - 1

east_north_bc = 1  ! 1: outflow b.c.  2: prescribed phi b.c.

!-----Initialise arrays
allocate(x(1:ni), y(1:nj), xc(1:ni), yc(1:nj), fx(1:ni), fy(1:nj), &
         phi(1:ni,1:nj), &
         u(1:ni,1:nj), v(1:ni,1:nj), &
         ae(1:ni,1:nj), an(1:ni,1:nj), aw(1:ni,1:nj), as(1:ni,1:nj), &
         ap(1:ni,1:nj), su(1:ni,1:nj), &
         stat=ierr)

x(1:ni) = 0.0_dp; y(1:nj) = 0.0_dp
xc(1:ni) = 0.0_dp; yc(1:nj) = 0.0_dp
fx(1:ni) = 0.0_dp; fy(1:nj) = 0.0_dp
phi(1:ni,1:nj) = 0.0_dp
u(1:ni,1:nj) = 0.0_dp; v(1:ni,1:nj) = 0.0_dp
ae(1:ni,1:nj) = 0.0_dp; an(1:ni,1:nj) = 0.0_dp 
aw(1:ni,1:nj) = 0.0_dp; as(1:ni,1:nj) = 0.0_dp
ap(1:ni,1:nj) = 0.0_dp; su(1:ni,1:nj) = 0.0_dp

! x grid size
if (expfx == 1.0_dp) then
  dx = (xmax-xmin) / real(nicv, dp)
else
  dx = (xmax-xmin)*(1.0_dp-expfx)/(1.0_dp-expfx**nicv)
end if
!-----Define x grid (cell faces)
x(1) = xmin
do i=2,nim1
  x(i) = x(i-1) + dx
  dx = dx*expfx
end do
x(ni) = x(nim1)  ! dummy value
!-----Cell centres
xc(1) = x(1)
do i=2,nim1
  xc(i) = 0.5_dp*(x(i)+x(i-1))
end do
xc(ni) = x(nim1)
! y grid size
if (expfy == 1.0_dp) then
  dy = (ymax-ymin) / real(njcv, dp)
else
  dy = (ymax-ymin)*(1.0_dp-expfy)/(1.0_dp-expfy*njcv)
end if
!-----Define y grid
y(1) = ymin
do j=2,njm1
  y(j) = y(j-1) + dy
  dy = dy*expfy
end do
y(nj) = y(njm1)
!-----Cell centres
yc(1) = y(1)
do j=2,njm1
  yc(j) = 0.5_dp*(y(j)+y(j-1))
end do
yc(nj) = y(njm1)
!-----Interpolation factors for CDS, fx=(xe-xP)/(xE-xP)
fx(1) = 0.0_dp
do i=2,nim1
  fx(i) = (x(i)-xc(i))/(xc(i+1)-xc(i))
end do
fy(1) = 0.0_dp
do j=2,njm1
  fy(j) = (y(j)-yc(j))/(yc(j+1)-yc(j))
end do
!-----velocities
u = vel
v = vel
!-----Initialise variable
! west
phi(1,2:njm1) = 1.0_dp
! south
phi(2:nim1,1) = 0.0_dp
! corners
phi(1,1) = 0.5_dp
if(east_north_bc == 2) then  ! prescribed phi
  phi(1:nim1,nj) = 1.0_dp
  phi(ni,1:njm1) = 0.0_dp
  phi(ni,nj) = 0.5_dp
end if
!-----Coefficients
do j=2,njm1
  do i=2,nim1
    ! velocity
    ue = u(i,j)
    uw = u(i-1,j)
    vn = v(i,j)
    vs = v(i,j-1)
    ! mass flux
    ge = ue*(y(j)-y(j-1))
    gw = -uw*(y(j)-y(j-1))  ! negative in terms of outer normal
    gn = vn*(x(i)-x(i-1))
    gs = -vs*(x(i)-x(i-1))
    ce = min(ge, 0.0_dp)
    cw = min(gw, 0.0_dp)
    cn = min(gn, 0.0_dp)
    cs = min(gs, 0.0_dp)
    ! coef matrix
    ae(i,j) = ce
    aw(i,j) = cw
    an(i,j) = cn
    as(i,j) = cs
    ap(i,j) = -(ae(i,j)+aw(i,j)+an(i,j)+as(i,j))
    su(i,j) = 0.0_dp
  end do
end do
!-----West boundary - Dirichlet b.c.
i = 2
do j=2,njm1
  su(i,j) = su(i,j) - aw(i,j)*phi(1,j)
  aw(i,j) = 0.0_dp
end do
!-----East boundary - outflow or prescribed phi
i = nim1
do j=2,njm1
  if (east_north_bc == 1) then  ! outflow
    ap(i,j) = ae(i,j) + ap(i,j)
  else   ! prescribed phi
    su(i,j) = su(i,j) - ae(i,j)*phi(ni,j)
  end if
  ae(i,j) = 0.0_dp
end do
!-----North boundary - outflow or prescribed phi
j = njm1
do i=2,nim1
  if (east_north_bc == 1) then  ! outflow
    ap(i,j) = an(i,j) + ap(i,j)
  else  ! prescribed phi
    su(i,j) = su(i,j) - an(i,j)*phi(i,nj)
  end if
  an(i,j) = 0.0_dp
end do
!-----South boundary - Dirichlet b.c.
j = 2
do i=2,nim1
  su(i,j) = su(i,j) - as(i,j)*phi(i,1)
  as(i,j) = 0.0_dp
end do
call tdma_solver(aw, ae, as, an, ap, su, phi, ni, nj)
! values at outlet planes, zero grad
if (east_north_bc == 1) then
  phi(ni,1:nj) = phi(nim1,1:nj)
  phi(1:ni,nj) = phi(1:ni,njm1)
end if
!-----print out
write(*,'(/,a,/)') '2D Convection of a Step Profile in a Uniform Flow &
                    Oblique to Grid Lines'
if (east_north_bc == 1) then
  write(*,*) 'Outflow conditions for east and north b.c. '
else
  write(*,*) 'Dirichlet conditions for east and north b.c. '
end if
write(*,*) 'UDS used for convection'
write(*,*) 'TDMA used for solver'
if (ni == nj) then
  write(*,'(/,a,/)') 'Profile along the diagonal (1,nj)-(ni,1):'
  do i=1,ni
    write(*,'(1x,2(1x,f7.2))') sqrt((x(i)-x(1))**2 + (y(ni-i)-y(ni))**2), phi(i,ni-i)
  end do
end if
call tecplot_write(x, y, u, v, phi, ni, nj, filename1)
deallocate(x, y, xc, yc, phi, u, v, aw, ae, as, an, ap, su, &
           stat=ierr)
end subroutine fvm2d_conv_oblique
!********************************************************************************
subroutine tecplot_write(x, y, u, v, phi, ni, nj, datafile)
implicit none
integer, parameter :: dp = selected_real_kind(15)  ! Double precision
integer :: i, j, ierr
integer, intent(in) :: ni, nj
character(80), intent(in) :: datafile
real(dp), intent(in) :: x(1:ni), y(1:nj)
real(dp), dimension(1:ni,1:nj) :: u, v, phi

write(*,'(/,a,/)') 'Data file written in Tecplot format'
open(unit=1, file=datafile, status="replace", iostat=ierr)
write(1,*) 'title = "Stagnation point flow 2D - output"'
write(1,*) 'variables = "x", "y", "u", "v", "phi"'
write(1,'(/,a,i4,3x,a,i4,3x,a)') 'zone i=', ni, 'j=', nj, 'f=point'
do j=1,nj
  do i=1,ni
    write(1,'(2x,5(1x,es9.2))') x(i), y(j), u(i,j), v(i,j), phi(i,j)
  end do
end do
close(1)        
end subroutine tecplot_write
!********************************************************************************
subroutine tdma_solver(aw, ae, as, an, ap, su, phi, ni, nj)
! Performs iterative line by line TDMA method for a two dimension problem.
implicit none
integer, parameter :: dp = selected_real_kind(15)  ! Double precision
integer, intent(in) :: ni, nj
real(dp), dimension(1:ni,1:nj), intent(in) :: aw, ae, as, an, ap, su
real(dp), dimension(1:ni,1:nj), intent(inout) :: phi
real(dp) :: res, res1, rsm
real(dp), parameter :: tol = 1.0e-4_dp
integer :: i, j, iter
integer, parameter :: maxit = 100

write(*,*) 'Line by line TDMA solver used.'
do iter=1,maxit
  res = 0.0_dp
  call wesweep(aw, ae, as, an, ap, su, phi, ni, nj)
  do j=2,nj-1
    do i=2,ni-1
      res = res + abs(ae(i,j)*phi(i+1,j) + aw(i,j)*phi(i-1,j) + &
            an(i,j)*phi(i,j+1) + as(i,j)*phi(i,j-1) + ap(i,j)*phi(i,j) - su(i,j))
    end do
  end do
  if (iter == 1) res1 = res
  rsm = res/res1
  write(*,'(a,i4,a,3x,a,es9.2)') 'Iter:', iter, ',', 'RSM = ', rsm
  if (res < 1.0e-10_dp) then
    write(*,*) 'TDMA solver - converged' 
    exit
  elseif (rsm < tol) then 
    write(*,*) 'TDMA solver - converged' 
    exit
  else if (iter == maxit) then
    write(*,*) 'TDMA solver - convergence not reached'
  end if
end do
end subroutine tdma_solver
!********************************************************************************
subroutine wesweep(aw, ae, as, an, ap, su, phi, ni, nj)
! Perform west to east sweeps for 2-D TDMA along vertical lines.
! East and west terms are assumed temperarily known.
implicit none
integer, parameter :: dp = selected_real_kind(15)  ! Double precision
integer :: i, j, ierr
integer, intent(in) :: ni, nj
real(dp), dimension(1:ni,1:nj), intent(in) :: aw, ae, as, an, ap, su
real(dp), dimension(1:ni,1:nj), intent(inout) :: phi
real(dp), allocatable, dimension(:) :: a, b, c, d, x

allocate(a(1:nj), b(1:nj), c(1:nj), d(1:nj), x(1:nj), stat=ierr)
do i=2,ni-1
  do j=2,nj-1
    a(j) = as(i,j)
    b(j) = ap(i,j)
    c(j) = an(i,j)
    d(j) = su(i,j) - aw(i,j)*phi(i-1,j) - ae(i,j)*phi(i+1,j)
  end do
  call tdma(a, b, c, d, x, nj)
  do j=2,nj-1
    phi(i,j) = x(j)
  end do
end do
deallocate(a, b, c, d, x, stat=ierr)
end subroutine wesweep
!********************************************************************************
subroutine tdma(a, b, c, d, x, n)
!  Tri-diagonol system of equations
!|b1 c1                   | | x1 | | d1 |
!|a2 b2 c2                | | x2 | | d2 |
!|   a3 b3 c3             | | x3 | | d3 |
!|      .. .. ..          |*| .. |=| .. |
!|         .. .. ..       | | .. | | .. |
!|          an-1 bn-1 cn-1| |xn-1| |dn-1|
!|                an   bn | | xn | | dn |
!  ith equation in the system:
!  a(i)x(i-1)+b(i)x(i)+c(i)x(i+1)=d(i)  
implicit none 
integer, parameter :: dp = selected_real_kind(15)  ! Double precision
integer :: i, ierr
integer, intent(in) :: n
real(dp), dimension(n), intent(in) :: a, b, c, d
real(dp), dimension(n), intent(out) :: x
real(dp), allocatable, dimension(:) :: p, q
real(dp) :: denom
!-----
allocate(p(1:n), q(1:n), stat=ierr)
p = 0.0_dp
q = 0.0_dp
!-----Forward elimination
! only solve from 2 to ni-1.
! 1 and n are boundary values.
p(2) = c(2) / b(2)
q(2) = d(2) / b(2)
do i=3,n-1
  denom = b(i) - a(i)*p(i-1)
  p(i) = c(i) / denom
  q(i) = (d(i)-a(i)*q(i-1)) / denom
end do
!-----Back substitution
x(n-1) = q(n-1)
do i=n-2,2,-1
  x(i) = q(i) - p(i)*x(i+1)
end do
deallocate(p, q, stat=ierr)
end subroutine tdma
!********************************************************************************
