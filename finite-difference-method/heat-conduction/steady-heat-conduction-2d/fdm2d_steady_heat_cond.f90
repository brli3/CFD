!*******************************************************************************
program main
!
! Solves 2-D steady state heat equation in a rectangular domain
! d(lambda(x,y) dT/dx)/dx + d(lambda(x,y) dT/dy)/dy + qdot = 0
! where lambda(x,y) is the thermal conductivity and
! qdot is the volumetric heat generation term.
! 
! Only Dirichlet boundary condition is implemented.
!
! Author: Ruipengyu Li 
! Modified: 15/04/2017
!
implicit none
integer, parameter :: dp = selected_real_kind(15)
integer :: i, j
integer :: nx, ny
real(dp) :: xmin, xmax, ymin, ymax
real(dp) :: tleft, tright, ttop, tbottom
real(dp), allocatable :: t_temp(:), t(:,:)
! number of grid points
nx = 7
ny = 5 
allocate(t_temp(nx*ny))
allocate(t(nx,ny))
! set up domain
xmin = 0.0_dp
xmax = 0.03_dp
ymin = 0.0_dp
ymax = 0.02_dp
! boundary temperatures
tleft = 300.0_dp
tright = 300.0_dp
tbottom = 300.0_dp
ttop = 300.0_dp

call fdm2d_steady_heat_cond(xmin, xmax, ymin, ymax, tleft, tright, &
                            ttop, tbottom, t_temp, nx, ny)

t = reshape(t_temp, shape(t))
! print out temperature in the spatial coordinate
do j=ny,1,-1
  write(*,'(/,1x,<nx>f7.1,/)') (t(i,j), i=1,nx)
end do
deallocate(t_temp)
deallocate(t)
stop
end program main
!*******************************************************************************
function LAMBDA(x, y)
! Thermal conductivity
implicit none
integer, parameter :: dp = selected_real_kind(15)
real(dp), intent(in) :: x, y
real(dp) ::  LAMBDA
LAMBDA = 20.0_dp
end function LAMBDA
!*******************************************************************************
function QDOT(x, y)
! Volumetric heat generation
implicit none
integer, parameter :: dp = selected_real_kind(15)
real(dp), intent(in) :: x, y
real(dp) :: QDOT
QDOT = 5.0e7_dp
end function QDOT
!*******************************************************************************
subroutine fdm2d_steady_heat_cond(xmin, xmax, ymin, ymax, tleft, tright, &
                                  ttop, tbottom, t, nx, ny)                       
!
! Solve FD linear equations for nx*ny variables including bourndary nodes.
! Boundary nodes provide boundary equations 
! 
! The nodes in the rectangular domain can be numbered by cartesian coordinate
! (i,j) or by a single index k as 
!
!  (ny-1)*nx+1     (ny-1)*nx+2      ...         ny*nx
!      ...             ...          ...          ...
!   (j-1)*nx+1      (j-1)*nx+2   (j-1)*nx+i  (j-1)*nx+nx
!      ...             ...          ...          ...
!     nx+1            nx+2          ...         2*nx
!       1               2           ...          nx
!
implicit none
integer, parameter :: dp = selected_real_kind(15)
integer :: i, j, n
integer, intent(in) :: nx, ny
real(dp), intent(in) :: xmin, xmax, ymin, ymax
real(dp), intent(in) :: tleft, tright, ttop, tbottom
real(dp), intent(out) :: t(nx*ny)
real(dp) :: x(nx), y(ny)
real(dp), allocatable :: a(:,:), t0(:), rhs(:)
real(dp) :: dx, dy, kw, ke, ks, kn, kp
allocate(a(nx*ny,nx*ny), t0(nx*ny), rhs(nx*ny))
a = 0.0_dp; rhs = 0.0_dp; t0 = 0.0_dp
! mesh
call xymesh(xmin, xmax, ymin, ymax, x, y, nx, ny)
! interior nodes
call interior(x, y, a, rhs, nx, ny)
! boundary conditions
call boundary(tleft, tright, ttop, tbottom, x, y, a, rhs, nx, ny)
! linear solver
call gauss_seidel(a, rhs, t, t0, nx*ny)
deallocate(a, rhs, t0)
end subroutine fdm2d_steady_heat_cond
!*******************************************************************************
subroutine xymesh(xmin, xmax, ymin, ymax, x, y, nx, ny)
! Generate x and y vectors
implicit none
integer, parameter :: dp = selected_real_kind(15)
integer :: i, j
integer, intent(in) :: nx, ny
real(dp), intent(in) :: xmin, xmax, ymin, ymax
real(dp), intent(out) :: x(nx), y(ny)
real(dp) :: dx, dy
x = 0.0_dp
y = 0.0_dp
dx = (xmax - xmin) / REAL(nx-1, kind=dp)
dy = (ymax - ymin) / REAL(ny-1, kind=dp)
! set grid
x(1) = xmin
do i=2,nx
  x(i) = x(i) + (i-1)*dx
end do
y(1) = ymin
do j=2,ny
  y(j) = y(j) + (j-1)*dy
end do
end subroutine xymesh
!*******************************************************************************
subroutine interior(x, y, a, rhs, nx, ny)
!
! Set up coefficients of interior nodes from the dicrete FD eqs.
!
! Nodes can be numbered with a single index k = (j-1)*nx + i.
! i and j are column and row indices.
!
!  (ny-1)*nx+1     (ny-1)*nx+2      ...         ny*nx
!      ...             ...          ...          ...
!   (j-1)*nx+1      (j-1)*nx+2   (j-1)*nx+i  (j-1)*nx+nx
!      ...             ...          ...          ...
!     nx+1            nx+2          ...         2*nx
!       1               2           ...          nx
!
implicit none
integer, parameter :: dp = selected_real_kind(15)
integer :: i, j  ! cartesian notation
integer :: k  ! single index notation (row by row from sw to ne)
integer :: kw, ke, ks, kn ! compass notation
integer, intent(in) :: nx, ny
real(dp), intent(in) :: x(nx), y(ny)
real(dp), intent(out) :: rhs(nx*ny)
real(dp), intent(out) :: a(nx*ny, nx*ny)
real(dp), external :: LAMBDA, QDOT
real(dp) :: cdw, cde, cds, cdn ! thermal conductivity
real(dp) :: aw, ae, as, an, ap ! coefficients
real(dp) :: dx, dy
! assuming uniform grid
dx = x(2) - x(1)
dy = y(2) - y(1)
! loop throughout interior domain
do i=2,nx-1
  do j=2,ny-1
    ! compass notation for neighbouring points
    k = (j-1)*nx + i
    kw = k - 1
    ke = k + 1
    ks = k - nx
    kn = k + nx
    ! thermal conductivity at halfway between two points 
    cdw = LAMBDA(0.5_dp*(x(i)+x(i-1)), y(j))
    cde = LAMBDA(0.5_dp*(x(i)+x(i+1)), y(j))
    cds = LAMBDA(x(i), 0.5_dp*(y(j)+y(j-1)))
    cdn = LAMBDA(x(i), 0.5_dp*(y(j)+y(j+1)))
    ! coefficients
    aw = cdw / dx**2
    ae = cde / dx**2
    as = cds / dy**2
    an = cdn / dy**2
    ap = -(aw + ae + as + an)
    ! banded elements in coefficient matrix
    a(k,kw) = aw
    a(k,ke) = ae
    a(k,ks) = as
    a(k,kn) = an
    a(k,k) = ap
    ! right hand side vector
    rhs(k) = -QDOT(x(i), y(j))
  end do
end do
end subroutine interior
!*******************************************************************************
subroutine boundary(tleft, tright, ttop, tbottom, x, y, a, rhs, nx, ny)
!
! Set up boundary conditions.
! Dirichlet B.C. only
!
! Nodes can be numbered with a single index k = (j-1)*nx + i.
! i and j are column and row indices.
!
!  (ny-1)*nx+1     (ny-1)*nx+2      ...         ny*nx
!      ...             ...          ...          ...
!   (j-1)*nx+1      (j-1)*nx+2   (j-1)*nx+i  (j-1)*nx+nx
!      ...             ...          ...          ...
!     nx+1            nx+2          ...         2*nx
!       1               2           ...          nx
!
implicit none
integer, parameter :: dp = selected_real_kind(15)
integer :: i, j  ! cartesian notation
integer :: k  ! single index notation (row by row from sw to ne)
integer :: kw, ke, ks, kn ! compass notation
integer, intent(in) :: nx, ny
real(dp), intent(in) :: tleft, tright, ttop, tbottom
real(dp), intent(in) :: x(nx), y(ny)
real(dp), intent(inout) :: rhs(nx*ny)
real(dp), intent(inout) :: a(nx*ny, nx*ny)
! left
i = 1
do j=2,ny-1
  k = (j-1)*nx + i
  a(k,k) = 1.0_dp
  rhs(k) = tleft
end do
! right
i = nx 
do j=2,ny-1
  k = (j-1)*nx + i
  a(k,k) = 1.0_dp
  rhs(k) = tright
end do
! bottom 
j = 1 
do i=1,nx
  k = (j-1)*nx + i
  a(k,k) = 1.0_dp
  rhs(k) = tbottom
end do
! top
j = ny 
do i=1,nx
  k = (j-1)*nx + i
  a(k,k) = 1.0_dp
  rhs(k) = ttop
end do
end subroutine boundary
!*******************************************************************************
subroutine gauss_seidel(a, b, x, x0, n)
!  Gauss-Seidel iterative method
!  with relaxation 
implicit none 
integer, parameter :: dp = selected_real_kind(15)
integer :: i,j,k
integer, parameter :: imax=500
integer, intent(in) :: n
real(dp), dimension(1:n,1:n), intent(in) :: a
real(dp), dimension(1:n), intent(in) :: x0, b
real(dp), dimension(1:n), intent(out) :: x
real(dp), dimension(1:n) :: x1, x2
real(dp) :: s, dx2, tol=1.0e-7_dp
real(dp), parameter :: rf=0.5
! initial values
x1=x0
x2=x1
do k=1,imax
  do i=1,n
    s=0
    do j=1,n
! update using known values    
      if(j<i) then  
        s=s+a(i,j)*x2(j)
      else if(j>i) then
        s=s+a(i,j)*x1(j)
      end if
    end do
    x2(i)=(b(i)-s)*rf/a(i,i)+(1-rf)*x1(i)
  end do
  dx2=0.0_dp
  do i=1,n
    dx2=dx2+(x1(i)-x2(i))**2
  end do
  dx2=SQRT(dx2)
  if(dx2<tol) exit
  x1=x2
end do
x=x2
end subroutine gauss_seidel 
!*******************************************************************************
