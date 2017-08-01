!--------------------------------------------------------------------------
! This program solves the integral form of a
! unsteady two-dimensional convection-diffusion
! equation for a general scalar phi in a stagnation point flow
!
! using finite volume method.
!
!   
!
! 
!                inlet, phi=0
!          |---------------------------|
!          |                           |
!          |                           |
!          |                           |
!          |                           |
!          |                           |
!phi=phi(y)|wall                       | outlet, d(phi)/dx=0
!          |                           |
!          |                           |
!          |                           |
!          |                           |
!          |                           |
!          |                           |
!          -----------------------------
!            symmetry, d(phi)/dy=0
!   
! with Cartesian grids and known velocity field
!
!
! SIP solver of Stone is used.
!
! Author: Ruipengyu Li 
! Modified: 31/07/2017
!
! Reference:
!   J. H. Ferziger and M. Peric, Computational Methods for Fluid Dynamics,
!   3rd ed. Springer Berlin Heidelberg, 2001.
!--------------------------------------------------------------------------
!*******************************************************************************
module stag_point_mod

implicit none

integer, parameter :: dp = selected_real_kind(15)  ! Double precision

contains

!********************************************************************************
subroutine tecplot_write(x, y, u, v, phi, time, ifile)
! Write out in tecplot format.
implicit none
integer :: i, j, ierr
integer :: ni, nj
integer, intent(in) :: ifile
real(dp), dimension(:), intent(in) :: x, y
real(dp), dimension(:,:), intent(in) :: u, v, phi
real(dp), intent(in) :: time
logical :: first_call = .true.

ni = size(x)
nj = size(y)
write(*,'(/,a,/)') 'Writting in Tecplot format'
if (first_call) then ! print header
  write(ifile,*) 'title = "Stagnation point flow 2D"'
  write(ifile,*) 'variables = "x", "y", "u", "v", "phi"'
  first_call = .false.
end if
write(ifile,'(/,a,es12.5)') 
write(ifile,'(a,2x,a,i3,3x,a,i3,3x,a,3x,a,es12.5)') &
     'zone', 'i=', ni, 'j=', nj, 'f=point', 'solutiontime = ', time
do j = 1, nj
  do i = 1, ni
    write(ifile,'(2x,5(1x,es9.2))') x(i), y(j), u(i,j), v(i,j), phi(i,j)
  end do
end do
end subroutine tecplot_write
!********************************************************************************
subroutine tdma_we(aw, ae, as, an, ap, su, phi)
! Performs iterative line by line TDMA method in a two dimension problem.
implicit none
integer :: ni, nj
integer :: i, j, iter
integer, parameter :: maxit = 1000
real(dp), dimension(:,:), intent(in) :: aw, ae, as, an, ap, su
real(dp), dimension(:,:), intent(out) :: phi
real(dp) :: res, res1, rsm
real(dp), parameter :: tol = 1.0e-4_dp

ni = size(aw(:,1))
nj = size(aw(1,:))

write(*,*) 'Line by line TDMA solver used.'
do iter = 1, maxit
  res = 0.0_dp
  call wesweep(aw, ae, as, an, ap, su, phi)
  do j = 2, nj-1
    do i = 2, ni-1
      res = res + abs(ae(i,j)*phi(i+1,j) + aw(i,j)*phi(i-1,j) + &
                      an(i,j)*phi(i,j+1) + as(i,j)*phi(i,j-1) + &
                      ap(i,j)*phi(i,j) - su(i,j))
    end do
  end do
  if (iter == 1) res1 = res
  rsm = res / res1
  write(*,'(a,i4,a,3x,a,es9.2)') 'Iter:', iter, ',', 'RSM = ', rsm
  if (rsm < tol) then 
    write(*,*) 'TDMA solver - converged' 
    exit
  else if (iter == 1 .and. res1 < 1.0e-10_dp) then
    write(*,*) 'TDMA solver - converged at first iter' 
    exit
  else if (iter == maxit) then
    write(*,*) 'TDMA solver - convergence not reached'
  end if
end do
end subroutine tdma_we
!********************************************************************************
subroutine tdma_sn(aw, ae, as, an, ap, su, phi)
! Performs iterative line by line TDMA method in a two dimension problem.
implicit none
integer :: ni, nj
integer :: i, j, iter
integer, parameter :: maxit = 1000
real(dp), dimension(:,:), intent(in) :: aw, ae, as, an, ap, su
real(dp), dimension(:,:), intent(out) :: phi
real(dp) :: res, res1, rsm
real(dp), parameter :: tol = 1.0e-4_dp

ni = size(aw(:,1))
nj = size(aw(1,:))

write(*,*) 'Line by line TDMA solver used.'
do iter = 1, maxit
  res = 0.0_dp
  call snsweep(aw, ae, as, an, ap, su, phi)
  do j = 2, nj-1
    do i = 2, ni-1
      res = res + abs(ae(i,j)*phi(i+1,j) + aw(i,j)*phi(i-1,j) + &
                      an(i,j)*phi(i,j+1) + as(i,j)*phi(i,j-1) + &
                      ap(i,j)*phi(i,j) - su(i,j))
    end do
  end do
  if (iter == 1) res1 = res
  rsm = res / res1
  write(*,'(a,i4,a,3x,a,es9.2)') 'Iter:', iter, ',', 'RSM = ', rsm
  if (rsm < tol) then 
    write(*,*) 'TDMA solver - converged' 
    exit
  else if (iter == 1 .and. res1 < 1.0e-10_dp) then
    write(*,*) 'TDMA solver - converged at first iter' 
    exit
  else if (iter == maxit) then
    write(*,*) 'TDMA solver - convergence not reached'
  end if
end do
end subroutine tdma_sn
!********************************************************************************
subroutine wesweep(aw, ae, as, an, ap, su, phi)
! Perform west to east sweeps for 2-D TDMA along vertical lines.
! East and west terms are assumed temperarily known.
implicit none
integer :: i, j, ierr
integer :: ni, nj
real(dp), dimension(:,:), intent(in) :: aw, ae, as, an, ap, su
real(dp), dimension(:,:), intent(out) :: phi
real(dp), allocatable, dimension(:) :: a, b, c, d, x

ni = size(aw(:,1))
nj = size(aw(1,:))

allocate(a(1:nj), b(1:nj), c(1:nj), d(1:nj), x(1:nj), stat=ierr)

do i = 2, ni-1
  do j = 2, nj-1
    a(j) = as(i,j)
    b(j) = ap(i,j)
    c(j) = an(i,j)
    d(j) = su(i,j) - aw(i,j)*phi(i-1,j) - ae(i,j)*phi(i+1,j)
  end do
  x(1) = phi(i,1)
  x(ni) = phi(i,ni)
  call tdma(a, b, c, d, x)
  do j = 2, nj-1
    phi(i,j) = x(j)
  end do
end do

deallocate(a, b, c, d, x, stat=ierr)
end subroutine wesweep
!********************************************************************************
subroutine snsweep(aw, ae, as, an, ap, su, phi)
! Perform south to north sweeps for 2-D TDMA along horizontal lines.
! East and west terms are assumed temperarily known.
implicit none
integer :: i, j, ierr
integer :: ni, nj
real(dp), dimension(:,:), intent(in) :: aw, ae, as, an, ap, su
real(dp), dimension(:,:), intent(out) :: phi
real(dp), allocatable, dimension(:) :: a, b, c, d, x

ni = size(aw(:,1))
nj = size(aw(1,:))

allocate(a(1:ni), b(1:ni), c(1:ni), d(1:ni), x(1:ni), stat=ierr)

do j = 2, nj-1
  do i = 2, ni-1
    a(i) = aw(i,j)
    b(i) = ap(i,j)
    c(i) = ae(i,j)
    d(i) = su(i,j) - as(i,j)*phi(i,j-1) - an(i,j)*phi(i,j+1)
  end do
  x(1) = phi(1,j)
  x(ni) = phi(ni,j)
  call tdma(a, b, c, d, x)
  do i = 2, ni-1
    phi(i,j) = x(i)
  end do
end do

deallocate(a, b, c, d, x, stat=ierr)
end subroutine snsweep
!********************************************************************************
subroutine tdma(a, b, c, d, x)
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
integer :: i, ierr
integer :: n
real(dp), dimension(:), intent(in) :: a, b, c, d
real(dp), dimension(:), intent(out) :: x
real(dp), allocatable, dimension(:) :: p, q

n = size(a(:))
allocate(p(1:n), q(1:n), stat=ierr)

!-----Forward elimination
! 1 and n are boundary values.
p(1) = 0.0_dp
q(1) = x(1)
do i = 2, n-1
  p(i) = c(i) / (b(i)-a(i)*p(i-1))
  q(i) = (d(i)-a(i)*q(i-1)) / (b(i)-a(i)*p(i-1))
end do
p(n) = 0.0_dp
q(n) = x(n)
!-----Back substitution
do i = n-1, 2, -1
  x(i) = q(i) - p(i)*x(i+1)
end do
deallocate(p, q, stat=ierr)
end subroutine tdma
!********************************************************************************
subroutine sipsol(aw, ae, as, an, ap, su, phi)
! SIP solver, ILU of Stone (1968)
implicit none
integer :: i, j, ierr
integer :: ni, nj, iter
integer, parameter :: maxit = 1000
real(dp), parameter :: alpha = 0.90_dp
real(dp), parameter :: tol = 1.0e-4_dp
real(dp) :: p1, p2, rsm, resl, res1
real(dp), dimension(:,:), intent(in) :: aw, ae, as, an, ap, su
real(dp), dimension(:,:), intent(out) :: phi
real(dp), allocatable, dimension(:,:) :: lw, ls, lpr, un, ue, res

ni = size(aw(:,1))
nj = size(aw(1,:))
allocate(lw(1:ni,1:nj), ls(1:ni,1:nj), lpr(1:ni,1:nj), &
         un(1:ni,1:nj), ue(1:ni,1:nj), res(1:ni,1:nj), &
         stat=ierr)
lw(:,:) = 0.0_dp; ls(:,:) = 0.0_dp; lpr(:,:) = 0.0_dp
un(:,:) = 0.0_dp; ue(:,:) = 0.0_dp; res(:,:) = 0.0_dp

!-----Calculate coefficients of [L] and [U] matrices
do j = 2, nj-1
  do i = 2, ni-1
    lw(i,j) = aw(i,j) / (1.0_dp + alpha*un(i-1,j))
    ls(i,j) = as(i,j) / (1.0_dp + alpha*ue(i,j-1))
    p1 = alpha * lw(i,j) * un(i-1,j) 
    p2 = alpha * ls(i,j) * ue(i,j-1)
    lpr(i,j) = 1.0_dp / (ap(i,j) + p1 + p2 - lw(i,j)*ue(i-1,j) - &
                         ls(i,j)*un(i,j-1))
    un(i,j) = (an(i,j) - p1) * lpr(i,j)
    ue(i,j) = (ae(i,j) - p2) * lpr(i,j)
  end do
end do
!-----Iterate and calculate residuals
do iter=1,maxit
  resl = 0.0_dp
  do j = 2, nj-1
    do i = 2, ni-1
      res(i,j) = su(i,j) - aw(i,j)*phi(i-1,j) - ae(i,j)*phi(i+1,j) - &
                 an(i,j)*phi(i,j+1) - as(i,j)*phi(i,j-1) - ap(i,j)*phi(i,j) 
      resl = resl + abs(res(i,j))
      res(i,j) = (res(i,j) - ls(i,j)*res(i,j-1) - &
                  lw(i,j)*res(i-1,j)) * lpr(i,j)
    end do
  end do
  if (iter == 1) res1 = resl
  rsm = resl / res1
  ! calculate increment
  do j = nj-1, 2, -1
    do i = ni-1, 2, -1
      res(i,j) = res(i,j) - un(i,j)*res(i,j+1) - ue(i,j)*res(i+1,j)
      phi(i,j) = phi(i,j) + res(i,j)
    end do
  end do
  ! check convergence
  write(*,'(a,i4,a,3x,a,es9.2)') 'Iter:', iter, ',', 'RSM = ', rsm
  if (rsm < tol) then 
    write(*,*) 'SIP solver - converged' 
    exit
  else if (iter == 1 .and. res1 < 1.0e-10_dp) then
    write(*,*) 'SIP solver - converged at first iter' 
    exit
  else if (iter == maxit) then
    write(*,*) 'SIP solver - convergence not reached'
  end if
end do
deallocate(lw, ls, lpr, un, ue, res, stat=ierr)
end subroutine sipsol
!********************************************************************************
end module stag_point_mod
!********************************************************************************
program fvm2d_stag_point
!
! Solves 2D scalar transport equation using finite volume method.
!
use stag_point_mod

implicit none
integer :: i, j, it, nt, pt
integer :: ni, nj, nim1, njm1
integer :: nicv, njcv ! no. of cell centres
integer :: isch, isol, itsch ! schemes
integer :: ierr ! error flag 
integer :: file1
character(len=:), allocatable :: filename1
real(dp) :: xmin, xmax, ymin, ymax
real(dp) :: expfx, expfy  ! grid expansion factor
real(dp) :: den, gam ! density and diff coef
real(dp) :: dx, dy
real(dp) :: ue, uw, vn, vs ! velocities at cell faces
real(dp) :: ge, gw, gn, gs ! mass flow rates
real(dp) :: ce, cw, cn, cs ! convection coeffs
real(dp) :: de, dw, dn, ds ! diffusion coeffs
real(dp) :: time, dt, ct
real(dp) :: fwall ! wall flux
real(dp), allocatable, dimension(:) :: x, y, xc, yc
real(dp), allocatable, dimension(:) :: fx, fy  ! interpolation coef
real(dp), allocatable, dimension(:,:) :: u, v  ! velocity at cell face
real(dp), allocatable, dimension(:,:) :: phi, phio, phioo
real(dp), allocatable, dimension(:,:) :: ae, aw, an, as, ap, su

filename1 = 'velocity_phi.dat'
file1 = 1
open(unit=file1, file=filename1, status='replace', iostat=ierr)

! grid settings
xmin = 0.0_dp
xmax = 1.0_dp
ymin = 0.0_dp
ymax = 1.0_dp
expfx = 1.0_dp
expfy = 1.0_dp
nicv = 20  ! no. of control volumes
njcv = 20  ! no. of control volumes

ni = nicv + 2 ! no. of cell centres inc boundary
nim1 = ni - 1 ! no. of cell faces
nj = njcv + 2
njm1 = nj - 1

den = 1.2_dp ! density
gam = 0.1_dp ! diffusion coef

isch = 2 ! 1:UDS 2:CDS
itsch = 2 ! 1:Explicit Euler 2:Implicit Euler 3:Crank Nicolson
isol = 2 ! 1:TDMA w-e sweep 2:TDMA s-n sweep 3:SIP

dt = 1.0e-0_dp ! time step
nt = 50 ! no. of time step
pt = 10 ! frequency of print

!-----Initialise arrays
allocate(x(1:ni), y(1:nj), xc(1:ni), yc(1:nj), fx(1:ni), fy(1:nj), &
         phi(1:ni, 1:nj), phio(1:ni, 1:nj), phioo(1:ni, 1:nj), &
         u(1:ni,1:nj), v(1:ni,1:nj), & 
         ae(1:ni,1:nj), an(1:ni,1:nj), aw(1:ni,1:nj), as(1:ni,1:nj), &
         ap(1:ni,1:nj), su(1:ni,1:nj), &
         stat=ierr)

if (ierr == 1) then
  write(*,*) 'Allocate returns.'
end if

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
  dx = (xmax-xmin) * (1.0_dp-expfx) / (1.0_dp-expfx**nicv)
end if
!-----Define x grid (cell faces)
x(1) = xmin
do i = 2, nim1
  x(i) = x(i-1) + dx
  dx = dx * expfx
end do
x(ni) = x(nim1)  ! dummy value
!-----Cell centres
xc(1) = x(1)
do i = 2, nim1
  xc(i) = 0.5_dp * (x(i)+x(i-1))
end do
xc(ni) = x(nim1)
! y grid size
if (expfy == 1.0_dp) then
  dy = (ymax-ymin) / real(njcv, dp)
else
  dy = (ymax-ymin) * (1.0_dp-expfy) / (1.0_dp-expfy*njcv)
end if
!-----Define y grid
y(1) = ymin
do j = 2, njm1
  y(j) = y(j-1) + dy
  dy = dy * expfy
end do
y(nj) = y(njm1)
!-----Cell centres
yc(1) = y(1)
do j = 2, njm1
  yc(j) = 0.5_dp * (y(j)+y(j-1))
end do
yc(nj) = y(njm1)
!-----Interpolation factors for CDS, fx=(xe-xP)/(xE-xP)
fx(1) = 0.0_dp
do i = 2, nim1
  fx(i) = (x(i)-xc(i)) / (xc(i+1)-xc(i))
end do
fy(1) = 0.0_dp
do j = 2, njm1
  fy(j) = (y(j)-yc(j)) / (yc(j+1)-yc(j))
end do
! velocities
do j = 1, nj
  do i = 1, ni
    u(i,j) = x(i)
    v(i,j) = -y(j)
  end do
end do
!-----Initialise variable
phio(:,:) = 0.0_dp
phi(2:ni,1:nj) = 0.0_dp
! left wall
phi(1,1:njm1) = 1.0_dp - (yc(1:njm1)-ymin) / (ymax-ymin)
!-----Time loop
time = 0.0_dp
do it = 1, nt
  ap = 0.
  time = time + dt
  ! update solution
  phio(:, :) = phi(:, :)
  ! set boundary vaule
  phi(1,1:njm1) = 1.0_dp - (yc(1:njm1)-ymin) / (ymax-ymin)
  do j = 2, njm1
    do i = 2, nim1
      ! velocity at cell faces
      ue = u(i,j)
      uw = u(i-1,j)
      vn = v(i,j)
      vs = v(i,j-1)
      ! mass flow rate
      ge = den * ue * (y(j)-y(j-1))
      gw = -den * uw * (y(j)-y(j-1))
      gn = den * vn * (x(i)-x(i-1))
      gs = -den * vs * (x(i)-x(i-1))
      ! diffusion term
      de = -gam * (y(j)-y(j-1)) / (xc(i+1)-xc(i))
      dw = -gam * (y(j)-y(j-1)) / (xc(i)-xc(i-1))
      dn = -gam * (x(i)-x(i-1)) / (yc(j+1)-yc(j))
      ds = -gam * (x(i)-x(i-1)) / (yc(j)-yc(j-1))
      ! convection term
      if (isch == 1) then  ! UDS
        ce = min(ge, 0.0_dp)
        cw = min(gw, 0.0_dp)
        cn = min(gn, 0.0_dp)
        cs = min(gs, 0.0_dp)
      else  ! CDS
        ce = ge * fx(i)
        cw = gw * (1.0_dp-fx(i-1))
        cn = gn * fy(j)
        cs = gs * (1.0_dp-fy(j-1))
      end if
      ! coef matrix except ap
      ae(i,j) = ce + de
      aw(i,j) = cw + dw
      an(i,j) = cn + dn
      as(i,j) = cs + ds
      su(i,j) = 0.0_dp ! no source term in equation
    end do
  end do

!-----Boundary conditions:
!   Dirichlet b.c. is included in the linear solver.
!   Neumann b.c. must be specified here.
! East boundary - outflow b.c.
  ae(nim1,2:njm1) = 0.0_dp
! South boundary - symmetry b.c.
  as(2:nim1,2) = 0.0_dp

!-----Time schemes
  ! explicit Euler
  if (itsch == 1) then
    do j = 2, njm1
      do i = 2, nim1
        dx = x(i) - x(i-1)
        ct = den / dt * dx * (y(j) - y(j-1))
        su(i,j) = su(i,j) + ct*phio(i,j) - & 
                 (ae(i,j)*phio(i+1,j) + aw(i,j)*phio(i-1,j) + &
                  an(i,j)*phio(i,j+1) + as(i,j)*phio(i,j-1)) + &
                 (ae(i,j) + aw(i,j) + an(i,j) + as(i,j))*phio(i,j)
        ap(i,j) = ct
        ae(i,j) = 0.0_dp
        aw(i,j) = 0.0_dp
        an(i,j) = 0.0_dp
        as(i,j) = 0.0_dp
      end do
    end do
  ! implicit Euler
  else if (itsch == 2) then
    do j = 2, njm1
      do i = 2, nim1
        dx = x(i) - x(i-1)
        ct = den / dt * dx * (y(j) - y(j-1))
        su(i,j) = su(i,j) + ct*phio(i,j) ! from time discretisation
        ap(i,j) = ct - (ae(i,j) + aw(i,j) + an(i,j) + as(i,j))
      end do
    end do
  ! Crank-Nicolson
  else
    do j = 2, njm1
      do i = 2, nim1
        dx = x(i) - x(i-1)
        ct = den / dt * dx * (y(j) - y(j-1))
        ae(i,j) = 0.5_dp * ae(i,j)
        aw(i,j) = 0.5_dp * aw(i,j)
        an(i,j) = 0.5_dp * an(i,j)
        as(i,j) = 0.5_dp * as(i,j)
        su(i,j) = su(i,j) + ct*phio(i,j) - & 
                 (ae(i,j)*phio(i+1,j) + aw(i,j)*phio(i-1,j) + &
                  an(i,j)*phio(i,j+1) + as(i,j)*phio(i,j-1)) + &
                 (ae(i,j) + aw(i,j) + an(i,j) + as(i,j))*phio(i,j)
        ap(i,j) = ct - (ae(i,j) + aw(i,j) + an(i,j) + as(i,j))
      end do
    end do
  end if
  ! call linear solver
  if (isol == 1) then
    call tdma_we(aw, ae, as, an, ap, su, phi)
  else if (isol == 2) then
    call tdma_we(aw, ae, as, an, ap, su, phi)
  else
    call sipsol(aw, ae, as, an, ap, su, phi)
  end if
  if (mod(it,pt) == 0) then
    ! update outlet and symmetry boundaries
    phi(2:nim1,1) = phi(2:nim1,2)
    phi(ni,1:nj) = phi(nim1,1:nj)
    ! write to file
    call tecplot_write(x, y, u, v, phi, time, file1)
  end if
end do

! West wall heat (scalar) flux
fwall = 0.0_dp
do j = 2, njm1
  fwall = fwall + gam * (y(j)-y(j-1)) * &
          (phi(2,j)-phi(1,j)) / (xc(2)-xc(1))
end do

!-----print out
write(*,'(/,a,/)') '2D Stagnation Point Flow'
if (isch == 1) then
  write(*,*) 'UDS for convection'
else if (isch == 2) then
  write(*,*) 'CDS for convection'
end if
write(*,*) 'CDS for diffusion'
write(*,*) 'SIP solver'
if (itsch == 1) then
  write(*,*) 'Explicit Euler'
else if (itsch == 2) then
  write(*,*) 'Implicit Euler'
else
  write(*,*) 'Crank-Nicolson'
end if
write(*,'(/,a,1x,es11.4)') 'dt = ', dt
write(*,'(a,1x,es11.4)') 'Time = ', time
write(*,'(a,1x,es11.4)') 'Courant = ', dt / dx 
write(*,'(a,1x,es11.4)') 'diff num = ', gam * dt / den / dx**2 
write(*,'(/,a,1x,f10.5)') 'Wall scalar flux =', fwall

deallocate(x, y, xc, yc, phi, phio, phioo, u, v, &
           aw, ae, as, an, ap, su, &
           stat=ierr)

close(file1)         

stop
end program fvm2d_stag_point
!********************************************************************************
