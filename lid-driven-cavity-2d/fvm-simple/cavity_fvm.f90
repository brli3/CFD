!*****************************************************************************
! Lid-Driven Cavity Flow

!   Finite volume method
!   Steady state SIMPLE algorithm
!   Uniform Staggered grid
!   UDS, CDS or Hybrid schemes for advection terms

! Finite volume method is found to converge faster and can handle
! larger Reynolds number than the finite difference method using
! CDS.

! Note that this code only works with uniform mesh.
! All the modules are in one file for the ease of compilation

! Ruipengyu Li
! Modified: 21/08/2017
!*****************************************************************************
module types_mod
implicit none
integer, parameter :: dp = selected_real_kind(15)
end module types_mod

module vars_mod
use types_mod
implicit none
integer :: imon, jmon, ipref, jpref
integer :: nicv, njcv, ni, nj, nim1, njm1
integer :: ierr, itsp, ntsp
integer :: nswpu, nswpv, nswpp
integer :: iadv ! 1:CDS 2:UDS 3:Hybrid
character(len=80) :: errmsg
real(dp) :: xstart, xend, ystart, yend
real(dp) :: urfu, urfv, urfp, error, tol
real(dp) :: resoru, resorv, resorm
real(dp) :: reynolds
real(dp), allocatable, dimension(:) :: x, y, xu, yv, xc, yc
real(dp), allocatable, dimension(:) :: sew, sns, sewu, snsv
real(dp), allocatable, dimension(:,:) :: u, v, p, pp, vis, den
real(dp), allocatable, dimension(:,:) :: uo, vo, po, viso, deno
real(dp), allocatable, dimension(:,:) :: uc, vc
real(dp), allocatable, dimension(:,:) :: ae, aw, an, as, ap, su, sp, du, dv
end module vars_mod

module linsolv_mod
use types_mod
implicit none

contains

subroutine lisolv(istart, jstart, ni, nj, phi)
use vars_mod, only: ap, an, as, ae, aw, su, sp
implicit none
integer, intent (in) :: istart, jstart, ni, nj
real(dp), dimension (ni,nj), intent (out) :: phi
real(dp), dimension (nj) :: a, b, c, d
integer :: i, j, jj, nim1, njm1, jsrm1
real(dp) :: term
nim1 = ni - 1
njm1 = nj - 1
jsrm1 = jstart - 1
a(jsrm1) = 0.0_dp
!-----commence w-e sweep
do i = istart, nim1
  c(jsrm1) = phi(i, jsrm1)
!-----commence s-n traverse
  do j = jstart, njm1
!-----assemble tdma coefficients
    a(j) = an(i, j)
    b(j) = as(i, j)
    c(j) = ae(i, j)*phi(i+1, j) + aw(i, j)*phi(i-1, j) + su(i, j)
    d(j) = ap(i, j)
!-----calculate coefficients of recurrence formula
    term = 1._dp/(d(j)-b(j)*a(j-1))
    a(j) = a(j)*term
    c(j) = (c(j)+b(j)*c(j-1))*term
  end do
!-----obtain new phi"s
  do jj = jstart, njm1
    j = nj + jsrm1 - jj
    phi(i, j) = a(j)*phi(i, j+1) + c(j)
  end do
end do
end subroutine lisolv

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
    c(j) = as(i,j)
    a(j) = ap(i,j)
    b(j) = an(i,j)
    d(j) = su(i,j) + aw(i,j)*phi(i-1,j) + ae(i,j)*phi(i+1,j)
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

subroutine tdma(a, b, c, d, x)
!  Tri-diagonol system of equations
!|a1 -b1                    | | x1 | | d1 |
!|-c2 a2 -b2                | | x2 | | d2 |
!|   -c3 a3 -b3             | | x3 | | d3 |
!|       .. .. ..           |*| .. |=| .. |
!|          .. .. ..        | | .. | | .. |
!|          -cn-1 an-1 -bn-1| |xn-1| |dn-1|
!|                -cn   an  | | xn | | dn |
!  ith equation in the system:
!   a(i)x(i) = b(i)x(i+1) + c(i)x(i-1) + d(i)
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
  p(i) = b(i) / (a(i)-c(i)*p(i-1))
  q(i) = (d(i)+c(i)*q(i-1)) / (a(i)-c(i)*p(i-1))
end do
!-----Back substitution
do i = n-1, 2, -1
  x(i) = p(i)*x(i+1) + q(i)
end do
deallocate(p, q, stat=ierr)
end subroutine tdma

end module linsolv_mod

module calc_mod
use types_mod
use linsolv_mod
implicit none
contains

subroutine calcu()
! u control volume using compass notation
use vars_mod
implicit none
integer :: i, j, n
real(dp) :: arean, areas, areaw, areae, vol
real(dp) :: gee, gww, gne, gnw, gse, gsw, gp, ge, gw, gn, gs
real(dp) :: ce, cw, cn, cs, cp, smp
real(dp) :: vise, visw, visn, viss, de, dw, dn, ds
real(dp) :: dudxe, dudxw, dvdxn, dvdxs
real(dp) :: resor

do j = 2, njm1
  do i = 3, nim1
    arean = sewu(i)
    areas = sewu(i)
    areaw = sns(j)
    areae = sns(j)
    vol = sewu(i)*sns(j)
    ! mass flux at 7 locations surrounding a u-cell
    gee = 0.5_dp*(den(i,j)+den(i+1,j))*u(i+1,j)
    gww = 0.5_dp*(den(i-2,j)+den(i-1,j))*u(i-1,j)
    gne = 0.5_dp*(den(i,j)+den(i,j+1))*v(i,j+1)
    gnw = 0.5_dp*(den(i-1,j)+den(i-1,j+1))*v(i-1,j+1)
    gse = 0.5_dp*(den(i,j-1)+den(i,j))*v(i,j)
    gsw = 0.5_dp*(den(i-1,j-1)+den(i-1,j))*v(i-1,j)
    gp = 0.5_dp*(den(i-1,j)+den(i,j))*u(i,j)
    ! mass flux at u-cell face centres
    ge = 0.5_dp*(gp+gee)
    gw = 0.5_dp*(gww+gp)
    gn = 0.5_dp*(gnw+gne)
    gs = 0.5_dp*(gsw+gse)
    ! convection coef (mass flow rate)
    ce = ge*areae
    cw = gw*areaw
    cn = gn*arean
    cs = gs*areas
    smp = ce - cw + cn - cs
    cp = max(smp, 0.0_dp)
    ! viscosity interpolated at north and south face centres
    vise = vis(i,j)
    visw = vis(i-1,j)
    visn = 0.25_dp*(vis(i-1,j)+vis(i-1,j+1)+vis(i,j+1)+vis(i,j))
    viss = 0.25_dp*(vis(i-1,j-1)+vis(i-1,j)+vis(i,j)+vis(i,j-1))
    ! diffusion coef
    de = vise*areae/(xu(i+1)-xu(i))
    dw = visw*areaw/(xu(i)-xu(i-1))
    dn = visn*arean/(y(j+1)-y(j))
    ds = viss*areas/(y(j)-y(j-1))
    if (iadv == 1) then !CDS
      ae(i,j) = de - 0.5_dp*ce
      aw(i,j) = dw + 0.5_dp*cw
      an(i,j) = dn - 0.5_dp*cn
      as(i,j) = ds + 0.5_dp*cs
    else if (iadv == 2) then ! UDS
      ae(i,j) = de + max(0.0_dp, -ce)
      aw(i,j) = dw + max(0.0_dp, cw)
      an(i,j) = dn + max(0.0_dp, -cn)
      as(i,j) = ds + max(0.0_dp, cs)
    else
      an(i, j) = max(abs(0.5*cn), dn) - 0.5*cn
      as(i, j) = max(abs(0.5*cs), ds) + 0.5*cs
      ae(i, j) = max(abs(0.5*ce), de) - 0.5*ce
      aw(i, j) = max(abs(0.5*cw), dw) + 0.5*cw
    end if
    ! pressure in source term
    su(i,j) = (p(i-1,j)-p(i,j))*0.5_dp*(areae+areaw)
    ! viscous terms in source term
    dudxe = (u(i+1,j)-u(i,j))/sew(i)
    dudxw = (u(i,j)-u(i-1,j))/sew(i-1)
    dvdxn = (v(i,j+1)-v(i-1,j+1))/sewu(i)
    dvdxs = (v(i,j)-v(i-1,j))/sewu(i)
    su(i,j) = su(i,j) + ((vise*dudxe-visw*dudxw)/sewu(i) &
            + (visn*dvdxn-viss*dvdxs)/sns(j))*vol
    su(i,j) = su(i,j) + cp*u(i,j)
    sp(i,j) = -cp
  end do
end do

! boundary

! residual
resoru = 0.0_dp
do j = 2, njm1
  do i = 3, nim1
    ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j) - sp(i,j)
    resor = ap(i,j)*u(i,j) - ae(i,j)*u(i+1,j) - aw(i,j)*u(i-1,j) &
          - an(i,j)*u(i,j+1) -as(i,j)*u(i,j-1) - su(i,j)
    resoru = resoru + abs(resor)
    ! under-relax
    ap(i,j) = ap(i,j)/urfu
    su(i,j) = su(i,j) + (1-urfu)*ap(i,j)*u(i,j)
    ! velocity correction term
    du(i,j) = 0.5*(areae+areaw)/ap(i,j)
  end do
end do

do n = 1, nswpu
  call lisolv(3, 2, ni, nj, u) 
end do
end subroutine calcu

subroutine calcv()
use vars_mod
implicit none
integer :: i, j, n
real(dp) :: arean, areas, areaw, areae, vol
real(dp) :: gne, gnw, gnn, gse, gsw, gss, gp, ge, gw, gn, gs
real(dp) :: ce, cw, cn, cs, cp, smp
real(dp) :: vise, visw, visn, viss, de, dw, dn, ds
real(dp) :: dudye, dudyw, dvdyn, dvdys
real(dp) :: resor

do j = 3, njm1
  do i = 2, nim1
    areae = snsv(j)
    areaw = snsv(j)
    arean = sew(i)
    areas = sew(i)
    vol = sew(i)*snsv(j)
    ! mass flux around v(i,j)
    gne = 0.5_dp*(den(i,j)+den(i+1,j))*u(i+1,j)
    gnw = 0.5_dp*(den(i-1,j)+den(i,j))*u(i,j)
    gnn = 0.5_dp*(den(i,j)+den(i,j+1))*v(i,j+1)
    gse = 0.5_dp*(den(i,j-1)+den(i+1,j-1))*u(i+1,j-1)
    gsw = 0.5_dp*(den(i-1,j-1)+den(i,j-1))*u(i,j-1)
    gss = 0.5_dp*(den(i,j-2)+den(i,j-1))*v(i,j-1)
    gp = 0.5_dp*(den(i,j-1)+den(i,j))*v(i,j)
    ! mass flux at face centre
    ge = 0.5_dp*(gne+gse)
    gw = 0.5_dp*(gnw+gsw)
    gn = 0.5_dp*(gp+gnn)
    gs = 0.5_dp*(gss+gp)
    ! coef
    ce = ge*areae
    cw = gw*areaw
    cn = gn*arean
    cs = gs*areas
    smp = ce - cw + cn - cs
    cp = max(smp, 0.0_dp)
    ! viscosity
    vise = 0.25_dp*(vis(i,j-1)+vis(i,j)+vis(i+1,j)+vis(i+1,j-1))
    visw = 0.25_dp*(vis(i-1,j-1)+vis(i-1,j)+vis(i,j)+vis(i,j-1))
    visn = vis(i,j)
    viss = vis(i,j-1)
    ! diff coef
    de = vise*areae/(x(i+1)-x(i))
    dw = visw*areaw/(x(i)-x(i-1))
    dn = visn*arean/(yv(j+1)-yv(j))
    ds = viss*areas/(yv(j)-yv(j-1))
    if (iadv == 1) then ! CDS
      ae(i,j) = de - 0.5_dp*ce
      aw(i,j) = dw + 0.5_dp*cw
      an(i,j) = dn - 0.5_dp*cn
      as(i,j) = ds + 0.5_dp*cs
    else if (iadv == 2) then ! UDS
      ae(i,j) = de + max(0.0_dp, -ce)
      aw(i,j) = dw + max(0.0_dp, cw)
      an(i,j) = dn + max(0.0_dp, -cn)
      as(i,j) = ds + max(0.0_dp, cs)
    else ! Hybrid
      an(i, j) = max(abs(0.5*cn), dn) - 0.5*cn
      as(i, j) = max(abs(0.5*cs), ds) + 0.5*cs
      ae(i, j) = max(abs(0.5*ce), de) - 0.5*ce
      aw(i, j) = max(abs(0.5*cw), dw) + 0.5*cw
    end if
    su(i,j) = (p(i,j-1)-p(i,j))*0.5_dp*(arean+areas)
    dudye = (u(i+1,j)-u(i+1,j-1))/snsv(j)
    dudyw = (u(i,j)-u(i,j-1))/snsv(j)
    dvdyn = (v(i,j+1)-v(i,j))/sns(j)
    dvdys = (v(i,j)-v(i,j-1))/sns(j-1)
    su(i,j) = su(i,j) + (vise*dudye-visw*dudyw)/sew(i)*vol &
            + (visn*dvdyn-viss*dvdys)/snsv(j)*vol
    sp(i,j) = -cp
    su(i,j) = su(i,j) + cp*v(i,j)
  end do
end do

! boundary

resorv = 0.0_dp
do j = 3, njm1
  do i = 2, nim1
    ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j) - sp(i,j)
    resor = ap(i,j)*v(i,j) - ae(i,j)*v(i+1,j) - aw(i,j)*v(i-1,j) &
          - an(i,j)*v(i,j+1) - as(i,j)*v(i,j-1) - su(i,j)
    resorv = resorv + abs(resor)
    ! under-relax
    ap(i,j) = ap(i,j)/urfv
    su(i,j) = su(i,j) + (1.0_dp-urfv)*ap(i,j)*v(i,j)
    ! p' equation term
    dv(i,j) = 0.5_dp*(arean+areas)/ap(i,j)
  end do
end do

do n = 1, nswpv
  call lisolv(2, 3, ni, nj, v) 
end do

end subroutine calcv

subroutine calcp()
use vars_mod
implicit none
integer :: i, j, n
real(dp) :: arean, areas, areaw, areae, vol
real(dp) :: ge, gw, gn, gs, dene, denw, denn, dens
real(dp) :: ce, cw, cn, cs, cp, smp
real(dp) :: de, dw, dn, ds, ppref

resorm = 0.0_dp
do j = 2, njm1
  do i = 2, nim1
    areae = sns(j)
    areaw = sns(j)
    arean = sew(i)
    areas = sew(i)
    vol = sew(i)*sns(j)
    dene = 0.5_dp*(den(i,j)+den(i+1,j))
    denw = 0.5_dp*(den(i-1,j)+den(i,j))
    denn = 0.5_dp*(den(i,j)+den(i,j+1))
    dens = 0.5_dp*(den(i,j-1)+den(i,j))
    ae(i,j) = dene*areae*du(i+1,j)
    aw(i,j) = denw*areaw*du(i,j)
    an(i,j) = denn*arean*dv(i,j+1)
    as(i,j) = dens*areas*dv(i,j)
    ge = dene*u(i+1,j)
    gw = denw*u(i,j)
    gn = denn*v(i,j+1)
    gs = dens*v(i,j)
    ce = ge*areae
    cw = gw*areaw
    cn = gn*arean
    cs = gs*areas
    smp = ce - cw + cn - cs
    su(i,j) = -smp 
    sp(i,j) = 0.0_dp
    resorm = resorm + abs(smp)
  end do
end do

! boundary

do j = 2, njm1
  do i = 2, nim1
    ! assemble
    ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j) - sp(i,j)
  end do
end do
! solve
pp(:,:) = 0.0_dp
do n = 1, nswpp
  call lisolv(2, 2, ni, nj, pp) 
end do
! update 
ppref = pp(ipref,jpref)
do j = 2, njm1
  do i = 2, nim1
    p(i,j) = p(i,j) + urfp*(pp(i,j)-ppref)
    if (i /= 2) u(i,j) = u(i,j) + du(i,j)*(pp(i-1,j)-pp(i,j))
    if (j /= 2) v(i,j) = v(i,j) + dv(i,j)*(pp(i,j-1)-pp(i,j))
  end do
end do
end subroutine calcp

end module calc_mod

module serv_mod

use types_mod
implicit none

contains

subroutine setgrid(nicv, njcv, xstart, xend, ystart, yend, &
                   x, xu, y, yv, sew, sewu, sns, snsv)
! Cell centred, backward staggered                 
implicit none
integer :: i, j, ni, nj
integer, intent(in) :: nicv, njcv
real(dp), intent(in) :: xstart, xend, ystart, yend
real(dp), dimension(:), intent(out) :: x, xu, y, yv, sew, sewu, sns, snsv
real(dp) :: dx, dy, dxu, dyv

x(:) = 0.0_dp
xu(:) = 0.0_dp
y(:) = 0.0_dp
yv(:) = 0.0_dp
sew(:) = 0.0_dp
sewu(:) = 0.0_dp
sns(:) = 0.0_dp
snsv(:) = 0.0_dp
ni = nicv + 2
nj = njcv + 2
xu(1) = xstart
xu(2) = xstart
dx = (xend-xstart)/nicv
do i = 3, ni
  xu(i) = xu(i-1) + dx
end do
x(1) = xu(2)
do i = 2, ni-1
  x(i) = 0.5_dp*(xu(i)+xu(i+1))
  sew(i) = xu(i+1) - xu(i)
end do
x(ni) = xu(ni)
do i = 3, ni-1
  sewu(i) = x(i) - x(i-1)
end do
yv(1) = ystart
yv(2) = ystart
dy = (yend-ystart)/njcv
do j = 3, nj
  yv(j) = yv(j-1) + dy
end do
y(1) = yv(2)
do j = 2, nj-1
  y(j) = 0.5_dp*(yv(j)+yv(j+1))
  sns(j) = yv(j+1) - yv(j)
end do
y(nj) = yv(nj)
do j = 3, nj-1
  snsv(j) = y(j) - y(j-1)
end do
end subroutine setgrid

subroutine centred_vel(ni, nj, u, v, uc, vc)
implicit none
integer :: i, j, nim1, njm1
integer, intent(in) :: ni, nj
real(dp), dimension(:,:), intent(in) :: u, v
real(dp), dimension(:,:), intent(out) :: uc, vc
nim1 = ni - 1
njm1 = nj - 1
uc(1,:) = u(2,:)
uc(ni,:) = u(ni,:)
vc(:,1) = v(:,2)
vc(:,nj) = v(:,nj)
uc(2:nim1,:) = 0.5_dp*(u(2:nim1,:)+u(3:ni,:))
vc(:,2:njm1) = 0.5_dp*(v(:,2:njm1)+v(:,3:nj))
end subroutine centred_vel

subroutine array_alloc()
use vars_mod
implicit none

allocate(x(1:ni), xu(1:ni), y(1:nj), yv(1:nj), sew(1:ni), &
         sewu(1:ni), sns(1:nj), snsv(1:nj), stat=ierr, errmsg=errmsg)
if (ierr /= 0) write(*,*) 'ALLOACAE ERROR! ', errmsg
x(:) = 0.0_dp
xu(:) = 0.0_dp
y(:) = 0.0_dp
yv(:) = 0.0_dp
sew(:) = 0.0_dp
sewu(:) = 0.0_dp
sns(:) = 0.0_dp
snsv(:) = 0.0_dp

allocate(u(1:ni,1:nj), v(1:ni,1:nj), p(1:ni,1:nj), pp(1:ni, 1:nj), &
         vis(1:ni,1:nj), den(1:ni,1:nj), stat=ierr, errmsg=errmsg)
if (ierr /= 0) write(*,*) 'ALLOCATE ERROR! ', errmsg
u(:,:) = 0.0_dp
v(:,:) = 0.0_dp
p(:,:) = 0.0_dp
pp(:,:) = 0.0_dp
vis(:,:) = 0.0_dp
den(:,:) = 0.0_dp

allocate(uc(1:ni,1:nj), vc(1:ni,1:nj), stat=ierr, errmsg=errmsg)
if (ierr /= 0) write(*,*) 'ALLOCATE ERROR! ', errmsg
uc(:,:) = 0.0_dp
vc(:,:) = 0.0_dp

allocate(ae(1:ni,1:nj), aw(1:ni,1:nj), an(1:ni,1:nj), as(1:ni,1:nj), &
         ap(1:ni,1:nj), su(1:ni,1:nj), sp(1:ni,1:nj), &
         du(1:ni,1:nj), dv(1:ni,1:nj), stat=ierr, errmsg=errmsg)
if (ierr /= 0) write(*,*) 'ALLOCATE ERROR! ', errmsg
ae(:,:) = 0.0_dp
aw(:,:) = 0.0_dp
an(:,:) = 0.0_dp
as(:,:) = 0.0_dp
ap(:,:) = 0.0_dp
su(:,:) = 0.0_dp
sp(:,:) = 0.0_dp
du(:,:) = 0.0_dp
dv(:,:) = 0.0_dp

allocate(uo(1:ni,1:nj), vo(1:ni,1:nj), po(1:ni,1:nj), viso(1:ni,1:nj), &
         deno(1:ni,1:nj), stat=ierr, errmsg=errmsg)
if (ierr /= 0) write(*,*) 'ALLOCATE ERROR! ', errmsg

end subroutine array_alloc

subroutine array_dealloc()
use vars_mod
implicit none

deallocate(x, xu, y, yv, sew, sewu, sns, snsv, stat=ierr)
if (ierr /= 0) write(*,*) 'DEALLOCATE ERROR'
deallocate(u, v, p, vis, den, uc, vc, stat=ierr)
if (ierr /= 0) write(*,*) 'DEALLOCATE ERROR'
deallocate(ae, aw, an, as, ap, su, sp, du, dv, stat=ierr)
if (ierr /= 0) write(*,*) 'DEALLOCATE ERROR'
if (allocated(uo)) deallocate(uo, stat=ierr)
if (ierr /= 0) write(*,*) 'DEALLOCATE ERROR'
if (allocated(vo)) deallocate(vo, stat=ierr)
if (ierr /= 0) write(*,*) 'DEALLOCATE ERROR'
if (allocated(po)) deallocate(po, stat=ierr)
if (ierr /= 0) write(*,*) 'DEALLOCATE ERROR'
if (allocated(viso)) deallocate(viso, stat=ierr)
if (ierr /= 0) write(*,*) 'DEALLOCATE ERROR'
if (allocated(deno)) deallocate(deno, stat=ierr)
if (ierr /= 0) write(*,*) 'DEALLOCATE ERROR'

end subroutine array_dealloc

end module serv_mod

program cavity_fvm
use vars_mod
use serv_mod
use calc_mod
implicit none
integer :: i, j, iprt
integer :: iter, maxit
real(dp) :: utop

!-----Initialize arrays
nicv = 58
njcv = 58
ni = nicv + 2
nj = njcv + 2
nim1 = ni - 1
njm1 = nj - 1

call array_alloc()
!-----Grid
xstart = 0.0_dp
xend = 1.0_dp
ystart = 0.0_dp
yend = 1.0_dp
call setgrid(nicv, njcv, xstart, xend, ystart, yend, &
             x, xu, y, yv, sew, sewu, sns, snsv)
! Fluid properties
utop = 1.0_dp
reynolds = 400.0_dp
den(:,:) = 1000.0_dp
vis(:,:) = den(:,:)*utop*(xend-xstart)/reynolds
u(2:,nj) = utop
! Control parameters
iadv = 3 ! 1: CDS, 2: UDS, 3: Hybrid
maxit = 90000
iprt = 500
urfu = 0.5_dp
urfv = 0.5_dp
urfp = 0.8_dp
nswpu = 3
nswpv = 3
nswpp = 3
ipref = ni/2
jpref = nj/2
imon = ipref
jmon = jpref
do iter = 1, maxit
  ! boundary
  u(2:,nj) = utop
!-----Solve
  call calcu()
  call calcv()
  call calcp()
!-----Convergence
  if (mod(iter, iprt) == 0) then
    write(*,'(a,i5,*(2x,a,es9.2))') &
          'Iter=', iter, 'URes=', resoru, 'VRes=', resorv, 'MRes=', resorm, &
          'Umon=', u(imon,jmon), 'Vmon=', v(imon,jmon), 'Pmon=', p(imon,jmon)
  end if 
  ! this is a crude criterion
  if (max(resoru, resorv, resorm) < 1.0e-10_dp) exit
end do
call centred_vel(ni, nj, u, v, uc, vc)
! Print
write(*,'(/,a)') 'Lid Driven Cavity Flow'
write(*,*) 'Finite Volume Method and SIMPLE Algorithm'
if (iadv == 1) then
  write(*,*) 'CDS for advection'
else if (iadv == 2) then
  write(*,*) 'UDS for advection'
else
  write(*,*) 'Hybrid for advection'
end if
write(*,'(/,2a,i3,3x,a,i3)') 'Grid: ', 'ni = ', nicv+1, 'nj = ', njcv+1
write(*,'(a,1x,es9.2)') 'Re = ', reynolds
write(*,'(a,1x,f5.2)') 'Top u-velocity = ', utop
write(*,'(a,i3,a,i3,a,3(1x,es12.4))') 'x, y and u at (', imon, ',', jmon, &
        ') = ', x(imon), y(jmon), uc(imon,jmon)
! write to tecplot
open(unit=1, file='result.dat', status='replace')
write(1,*) 'title="Lid Driven Cavity Flow"'
write(1,*) 'variables="x", "y", "u", "v", "p"'
write(1,'(a,2x,a,i3,2x,a,i3,2x,a)') 'zone', 'i=', ni, 'j=', nj, 'f=point'
do j = 1, nj
  do i = 1, ni
    write(1,'(*(1x,es14.7))') x(i), y(j), uc(i,j), vc(i,j), p(i,j)
  end do
end do
! write to files
open(unit=2, file='uvel.txt', status='replace')
do j = 1, nj
  write(2,'(*(1x,es14.7))') x(imon), y(j), uc(imon,j)
end do
open(unit=3, file='vvel.txt', status='replace')
do i = 1, ni
  write(3,'(*(1x,es14.7))') x(i), y(jmon), vc(i,jmon)
end do
call array_dealloc()
stop
end program cavity_fvm
