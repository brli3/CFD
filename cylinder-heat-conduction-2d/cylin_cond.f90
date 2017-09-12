!*****************************************************************************
! Cylinder Heat Conduction

!   Finite volume method
!   Steady state, SIMPLE algorithm
!   Uniform, Staggered grid
!   UDS, CDS or Hybrid schemes for advection terms

! Note that this code only deals with cartesian, uniform mesh.

! Ruipengyu Li
! Modified: 04/09/2017
!*****************************************************************************
module vars_mod
implicit none
character(len=80) :: errmsg
integer, parameter :: dp = selected_real_kind(15)
integer, parameter :: max_nq = 10 ! max no. of variables to solve
integer :: nq ! 1=u, 2=v, 3=p, 4=T
integer :: ierr
integer :: ni, nj, nim1, njm1, nim2, njm2
integer :: icrd ! 1:Cartisian 2:Cylindrical 3:Polar
integer :: iadv ! 1:CDS 2:UDS 3:Hybrid 4:Power-law
integer :: imon = 2, jmon = 2, ipref = 1, jpref = 1
integer :: iter = 0, print_iter = 100, max_iter = 10000
integer, dimension(max_nq) :: nsweep = 3 ! number of sweeps for TDMA
logical :: incl_visor = .false. ! viscous source terms in momentum eqs
logical :: incl_falsor = .false. ! false source term in momentum eqs
logical :: llast = .false.
logical :: ltec = .false., ltxt = .false.
logical, dimension(max_nq) :: lsolve = .false.
logical, dimension(max_nq) :: lprint = .false.
real(dp) :: xstart, xend, ystart, yend, rstart
!real(dp) :: diff, flow, acof
real(dp) :: den_ini = 1.0_dp, cp_ini = 1.0_dp ! initial rho and cp
real(dp), dimension(max_nq) :: urf = 0.8_dp ! under-relaxation factor
real(dp), dimension(max_nq) :: resor = 0.0_dp ! residual
real(dp), allocatable, dimension(:) :: x, xu, xdif, sew, sewu, xcvi, xcvip
real(dp), allocatable, dimension(:) :: y, yv, ydif, sns, snsv, ycvr, ycvrs
real(dp), allocatable, dimension(:) :: arx, arxj, arxjp, r, rmn, sx, sxmn
real(dp), allocatable, dimension(:) :: fv, fvp, fx, fxm
real(dp), allocatable, dimension(:) :: fy, fym
real(dp), allocatable, dimension(:,:) :: u, v, p, pp, du, dv, t, uc, vc
real(dp), allocatable, dimension(:,:) :: den, gam, hc
real(dp), allocatable, dimension(:,:) :: ae, aw, an, as, ap, su, sp
end module vars_mod

module serv_mod
use vars_mod, only: dp
implicit none
contains

function DAPEC(flow, diff, iadv)
! calculate D*A(Pe) for aw and as
implicit none
integer, intent(in) :: iadv ! 1:CDS, 2:UDS, 3:Hybrid, 4:Power-law
real(dp), intent(in) :: diff, flow 
real(dp) :: DAPEC 
real(dp) :: temp
if (flow == 0.0_dp) then ! diffusion only
  DAPEC = diff
  return
end if
select case(iadv)
case(1) ! CDS
  DAPEC = diff - 0.5_dp*abs(flow) 
case(2) ! UDS
  DAPEC = diff 
case(3) ! Hybrid
  DAPEC = max(0.0_dp, diff-0.5_dp*abs(flow))
case(4) ! Power-law
  DAPEC = max(0.0_dp, (1.0_dp-0.1_dp*abs(flow/diff))**5)
  DAPEC = DAPEC * diff
case default ! UDS
  DAPEC = diff 
end select
end function DAPEC

function ARIMEAN(phi1, phi2, wf2)
! arithmetic mean
implicit none
real(dp), intent(in) :: phi1, phi2, wf2
real(dp) :: ARIMEAN
ARIMEAN = phi1 * (1 - wf2) + phi2 * wf2
end function ARIMEAN

function HARMEAN(phi1, phi2, wf2)
! harmonic mean
implicit none
real(dp), intent(in) :: phi1, phi2, wf2
real(dp) :: HARMEAN
HARMEAN = phi1 * phi2 / (phi1 * (1.0_dp - wf2) + phi2 * wf2)
end function HARMEAN

end module serv_mod

module case_mod
! case of study
use vars_mod
implicit none
real(dp) :: den0, vis0, gam0, hc0, reyno
real(dp) :: t_hot, t_cold, t_ref, t_mean
real(dp) :: prand, volexp, rayle, qwall, nuss
real(dp) :: gravy
contains

subroutine control()
! mesh input
implicit none 
icrd = 1 ! 1:Cartesian 2:Axisym cylind 3:polar
! domain and grid points
ni = 7
nj = 7
xstart = 0.0_dp
xend = 1.0_dp
ystart = 0.0_dp
yend = 1.0_dp
rstart = 1.0_dp
! control params
lsolve(1:3) = .true.
lsolve(4) = .true.
lprint(1:3) = .true.
lprint(4) = .true.
incl_falsor = .true.
ltec = .true.
ltxt = .false.
iadv = 1 ! 1: CDS, 2: UDS, 3: Hybrid
max_iter = 999999 ! max iteration
print_iter = 500 ! interval for screen print
urf(1) = 0.3_dp 
urf(2) = 0.3_dp
urf(3) = 0.9_dp
urf(4) = 0.6_dp
nsweep(1) = 3
nsweep(2) = 3
nsweep(3) = 3
nsweep(4) = 3
imon = ni / 2
jmon = nj / 2
end subroutine control

subroutine initial()
! initial values. Give unchanged b.c.
implicit none
reyno = 100.0_dp
vis0 = 1.0e-5_dp
den0 = 1.0e0_dp
hc0 = 1.0e3_dp
t_hot = 1.0_dp
t_cold = 0.0_dp
t_ref = 0.0_dp
gravy = -0.981_dp
prand = 0.710_dp
rayle = 1.0e3_dp

gam0 = vis0*hc0/prand
volexp = rayle*vis0**2/den0**2/abs(gravy) &
       /(t_hot-t_cold)/(xend-xstart)**3/prand
den(:,:) = den0
hc(:,:) = hc0
t(:,:) = t_hot
t(ni,:) = t_cold
end subroutine initial

subroutine density()
! specify fluid density as function of temperature
implicit none
end subroutine density

subroutine boundary()
! update boundary
implicit none
t(2:ni,1) = t(2:ni,2) ! adiabatic
t(2:ni,nj) = t(2:ni,njm1) ! adiabatic
end subroutine boundary

subroutine gamsor()
! set and update diffusion coefficients
! store sp and su of variables
! set up additional source terms
implicit none
integer :: i, j
! diffusion coefs
gam(:,:) = vis0
if (nq == 4) gam(:,:) = gam0
if (nq == 4) then 
  gam(:,1) = 0.0_dp
  gam(:,nj) = 0.0_dp
end if
do j = 2, njm1
  do i = 2, nim1
    if (nq == 2) then
      if (j /= 2) then
        ! y buoyancy term
        t_mean = 0.5_dp*(t(i,j)+t(i,j-1))
        su(i,j) = -den0*gravy*volexp*t_mean
      end if
    end if
  end do
end do
end subroutine gamsor

subroutine output()
! output to screen and files.
implicit none
integer :: i, j
if (.not.llast) then
  if (mod(iter, print_iter) == 0) then
    write(*,'(/,a,i6,4(2x,a,es9.2),/,11x,*(2x,a,es12.5))') &
            'Iter=', iter, 'URes=', resor(1), 'VRes=', resor(2), 'MRes=', resor(3), &
            'TRes=', resor(4), 'U=', u(imon,jmon), 'V=', v(imon,jmon), &
            'P=', p(imon,jmon), 'T=', t(imon,jmon)
  end if 
  ! this is a rather crude criterion
  if (maxval(resor) < 1.0e-10_dp .or. iter == max_iter) llast = .true.
else
! write to screen
  qwall = 0.0_dp
  do j = 2, njm1
    qwall = qwall + gam0*(y(j)-y(j-1))*(t(1,j)-t(2,j))/(x(2)-x(1))
  end do
  nuss = qwall / gam0
  write(*,'(/,a,/)') 'Natural Convection in a Square Cavity'
  write(*,*) 'Finite Volume Method, SIMPLE Algorithm'
  if (iadv == 1) then
    write(*,*) 'CDS for advection'
  else if (iadv == 2) then
    write(*,*) 'UDS for advection'
  else
    write(*,*) 'Hybrid for advection'
  end if
  write(*,'(/,2a,i3,3x,a,i3)') 'Grid: ', 'x-cv = ', ni-2, 'y-cv = ', nj-2
  write(*,'(/,a)') 'Values at monitor:'
  write(*,'(*(a,es10.3,3x))') 'x_mon =', x(imon), 'y_mon =', y(jmon)
  write(*,'(*(a,es18.10,3x))') 'u_mon =', uc(imon,jmon), &
          'v_mon =', vc(imon,jmon), 't_mon =', t(imon,jmon)
  write(*,'(/,a,1x,es9.2)') 'Ra = ', rayle
  write(*,'(a,1x,f5.2)') 'Pr = ', vis0*hc0/gam0
  write(*,'(a,1x,f7.4)') 'Nu = ', nuss
  ! write to tecplot
  if (ltec) then
    open(unit=1, file='result.dat', status='replace')
    write(1,*) 'title="Natural Convection in a Square Cavity"'
    write(1,*) 'variables="x", "y", "u", "v", "p", "T"'
    write(1,'(a,2x,a,i3,2x,a,i3,2x,a)') 'zone', 'i=', ni, 'j=', nj, 'f=point'
    do j = 1, nj
      do i = 1, ni
        write(1,'(*(1x,es14.7))') x(i), y(j), uc(i,j), vc(i,j), &
                                  p(i,j), t(i,j)
      end do
    end do
    close(1)
  end if
  ! write to text files
  if (ltxt) then
    open(unit=2, file='uvel.txt', status='replace')
    do j = 1, nj
      write(2,'(*(1x,es14.7))') x(imon), y(j), uc(imon,j)
    end do
    open(unit=3, file='vvel.txt', status='replace')
    do i = 1, ni
      write(3,'(*(1x,es14.7))') x(i), y(jmon), vc(i,jmon)
    end do
    close(2); close(3)
  end if
end if
end subroutine output

end module case_mod
!*****************************************************************************
program main
use vars_mod
use case_mod
implicit none
call control()
call array_alloc()
call grid()
call initial()
do iter = 1, max_iter
  call boundary()
  if (lsolve(1)) then
    !call calcu()
    call calcu_()
    call calcv()
    call calcp()
  end if
  if (lsolve(4)) call calct()
  call output()
  if (llast) exit
end do
call postproc()
call output()
call array_dealloc()
stop
end program main
!*****************************************************************************
subroutine grid()
! Cell centred, backward staggered, uniform 
use vars_mod
implicit none
integer :: i, j
real(dp) :: dx, dy
!                   w                 e
!        W          |        P        |         E
!                   ^                 ^
!      x(i-1)     xu(i)     x(i)    xu(i+1)   x(i+1)

! wf_e = (x_e-x_P)/(x_E-x_P)
! phi_e = phi_E * wf_e + phi_P * (1 - wf_e)
! Uniform grid is used in this program
nim1 = ni - 1
nim2 = ni - 2
njm1 = nj - 1
njm2 = nj - 2
dx = (xend-xstart)/(ni-2)
xu(2) = xstart
do i = 3, ni
  xu(i) = xu(i-1) + dx
end do
sew(2:nim1) = xu(3:ni) - xu(2:nim1)
x(1) = xu(2) ! = 0
x(2:nim1) = 0.5_dp*(xu(2:nim1)+xu(3:ni)) ! cell-centred
x(ni) = xu(ni)
xdif(2:ni) = x(2:ni) - x(1:nim1) ! distance between nodes
sewu(3:nim1) = xdif(3:nim1) ! u-cv
sewu(3) = sewu(3) + xdif(2) ! u-cv at west boundary
sewu(nim1) = sewu(nim1) + xdif(ni) ! u-cv at east boundary
xcvi(3:nim2) = 0.5_dp*sew(3:nim2) ! first half of cv
xcvip(3:nim2) = xcvi(3:nim2) ! second half of cv
xcvip(2) = sew(2) ! west boundary
xcvi(nim1) = sew(nim1) ! east boundary
dy = (yend-ystart)/(nj-2)
yv(2) = ystart
do j = 3, nj
  yv(j) = yv(j-1) + dy
end do
sns(2:njm1) = yv(3:nj) - yv(2:njm1)
y(1) = yv(2) ! = 0
y(2:njm1) = 0.5_dp*(yv(2:njm1)+yv(3:nj)) 
y(nj) = yv(nj)
ydif(2:nj) = y(2:nj) - y(1:njm1)
snsv(3:njm1) = ydif(3:njm1)
snsv(3) = snsv(3) + ydif(2)
snsv(njm1) = snsv(njm1) + ydif(nj)
r(1) = rstart ! starting radius
if (icrd == 1) then ! Cartesian
  rmn(:) = 1.0_dp ! norminal radius
  r(:) = 1.0_dp
else ! Cylind and polar
  do j = 2, nj ! r(1) specified
    r(j) = r(j-1) + ydif(j)
  end do
  rmn(2) = r(1) ! cell face starts from 2
  do j = 3, nj
    rmn(j) = rmn(j-1) + sns(j-1) ! radial position of v(i,j)
  end do
end if
sx(:) = 1.0_dp
sxmn(:) = 1.0_dp
if (icrd == 3) then ! polar
  sx(:) = r(:)
  sxmn(2:) = rmn(2:) ! cell face starts from 2
end if
ycvr(2:njm1) = r(2:njm1)*sns(2:njm1) ! w-e area
arx(2:njm1) = ycvr(2:njm1) ! area normal to x
if (icrd == 3) arx(:) = sns(:)
ycvrs(4:njm2) = 0.5_dp*(r(4:njm2)+r(3:njm2-1))*ydif(4:njm2) ! w-e area of v(i,j)
ycvrs(3) = 0.5_dp*(r(3)+r(1))*snsv(3)
ycvrs(njm1) = 0.5_dp*(r(nj)+r(njm2))*snsv(njm1)
if (icrd == 2) then ! cylind
  arxj(3:njm2) = 0.25_dp*(1.0_dp+rmn(3:njm2)/r(3:njm2))*arx(3:njm2)
  arxjp(3:njm2) = arx(3:njm2) - arxj(3:njm2)
else
  arxj(3:njm2) = 0.5_dp*arx(3:njm2)
  arxjp(3:njm2) = arxj(3:njm2)
end if
arxjp(2) = arx(2) ! south boundary
arxj(njm1) = arx(njm1) ! north boundary
! interpolation factors
fv(3:njm2) = arxjp(3:njm2) / arx(3:njm2)
fvp(3:njm2) = 1.0_dp - fv(3:njm2)
fx(2:ni) = 0.5_dp*sew(1:nim1)/xdif(2:ni)
fxm(2:ni) = 1.0_dp - fx(2:ni)
fy(2:nj) = 0.5_dp*sns(1:njm1)/ydif(2:nj)
fym(2:nj) = 1.0_dp - fy(2:nj)
! initialise
den(:,:) = den_ini
hc(:,:) = cp_ini
! print to check
open(10, file='grid params.txt', status='replace')
write(10,*), 'Coordinates:', icrd
write(10,*), 'x grid'
write(10,'(a,*(f5.2,1x))')'x:     ', x(:)
write(10,'(a,*(f5.2,1x))')'xu:    ', xu(:)
write(10,'(a,*(f5.2,1x))')'xdif:  ', xdif(:)
write(10,'(a,*(f5.2,1x))')'sew:   ', sew(:)
write(10,'(a,*(f5.2,1x))')'sewu:  ', sewu(:)
write(10,'(a,*(f5.2,1x))')'xcvi:  ', xcvi(:)
write(10,'(a,*(f5.2,1x))')'xcvip: ', xcvip(:)

write(10,*), 'y grid'
write(10,'(a,*(f5.2,1x))')'y:     ', y(:)
write(10,'(a,*(f5.2,1x))')'yv:    ', yv(:)
write(10,'(a,*(f5.2,1x))')'ydif:  ', ydif(:)
write(10,'(a,*(f5.2,1x))')'sns:   ', sns(:)
write(10,'(a,*(f5.2,1x))')'snsv:  ', snsv(:)

write(10,*), 'r scale'
write(10,'(a,*(f5.2,1x))')'r:     ', r(:)
write(10,'(a,*(f5.2,1x))')'rmn:   ', rmn(:)
write(10,'(a,*(f5.2,1x))')'sx:    ', sx(:)
write(10,'(a,*(f5.2,1x))')'sxmn:  ', sxmn(:)
write(10,'(a,*(f5.2,1x))')'ycvr:  ', ycvr(:)
write(10,'(a,*(f5.2,1x))')'ycvrs: ', ycvrs(:)
write(10,'(a,*(f5.2,1x))')'arx:   ', arx(:)
write(10,'(a,*(f5.2,1x))')'arxj:  ', arxj(:)
write(10,'(a,*(f5.2,1x))')'arxjp: ', arxjp(:)

write(10,*), 'interp'
write(10,'(a,*(f5.2,1x))')'fx:    ', fx(:)
write(10,'(a,*(f5.2,1x))')'fxm:   ', fxm(:)
write(10,'(a,*(f5.2,1x))')'fv:    ', fv(:)
write(10,'(a,*(f5.2,1x))')'fvp:   ', fvp(:)
write(10,'(a,*(f5.2,1x))')'fy:    ', fy(:)
write(10,'(a,*(f5.2,1x))')'fym:   ', fym(:)
close(10)

end subroutine grid

subroutine calcu()
! u control volume using compass notation
use vars_mod
use case_mod, only: gamsor
implicit none
integer :: i, j, n
real(dp) :: arean, areas, areaw, areae, vol
real(dp) :: gee, gww, gne, gnw, gse, gsw, gp, ge, gw, gn, gs
real(dp) :: ce, cw, cn, cs, cp, smp
real(dp) :: vise, visw, visn, viss, de, dw, dn, ds
real(dp) :: dudxe, dudxw, dvdxn, dvdxs
real(dp) :: reseq
nq = 1
call clearsor()
call gamsor()
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
    vise = gam(i,j)
    visw = gam(i-1,j)
    visn = 0.25_dp*(gam(i-1,j)+gam(i-1,j+1)+gam(i,j+1)+gam(i,j))
    viss = 0.25_dp*(gam(i-1,j-1)+gam(i-1,j)+gam(i,j)+gam(i,j-1))
    ! diffusion coef
    de = vise*areae/(xu(i+1)-xu(i))
    dw = visw*areaw/(xu(i)-xu(i-1))
    dn = visn*arean/(y(j+1)-y(j))
    ds = viss*areas/(y(j)-y(j-1))
    call getanb(i, j, iadv, ce, cw, cn, cs, de, dw, dn, ds, &
                 ae(i,j), aw(i,j), an(i,j), as(i,j))
    su(i,j) = su(i,j)*vol
    sp(i,j) = sp(i,j)*vol
    if (incl_visor) then
      ! viscous terms in source term
      dudxe = (u(i+1,j)-u(i,j))/sew(i)
      dudxw = (u(i,j)-u(i-1,j))/sew(i-1)
      dvdxn = (v(i,j+1)-v(i-1,j+1))/sewu(i)
      dvdxs = (v(i,j)-v(i-1,j))/sewu(i)
      su(i,j) = su(i,j) + ((vise*dudxe-visw*dudxw)/sewu(i) &
              + (visn*dvdxn-viss*dvdxs)/sns(j))*vol
    end if
    if (incl_falsor) then
      ! false source to stabilise
      su(i,j) = su(i,j) + cp*u(i,j)
      sp(i,j) = sp(i,j) - cp
    end if
    ! pressure gradient term
    du(i,j) = vol/xdif(i) ! linear interp for boundary pressure
    su(i,j) = su(i,j) + du(i,j)*(p(i-1,j)-p(i,j))
  end do
end do
resor(1) = 0.0_dp
do j = 2, njm1
  do i = 3, nim1
    ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j) - sp(i,j)
    reseq = ap(i,j)*u(i,j) - ae(i,j)*u(i+1,j) - aw(i,j)*u(i-1,j) &
          - an(i,j)*u(i,j+1) -as(i,j)*u(i,j-1) - su(i,j)
    resor(1) = resor(1) + abs(reseq)
    ap(i,j) = ap(i,j)/urf(1)
    su(i,j) = su(i,j) + (1-urf(1))*ap(i,j)*u(i,j)
    du(i,j) = du(i,j)/ap(i,j)
  end do
end do
do n = 1, nsweep(1)
  call lisolv(3, 2, ni, nj, u) 
end do
end subroutine calcu

subroutine calcu_()
! u control volume using compass notation
use vars_mod
use serv_mod, only: DAPEC, ARIMEAN, HARMEAN
use case_mod, only: gamsor
implicit none
integer :: i, j, n
real(dp) :: den_up, den_ue ! density at u(P) and u(E)
real(dp) :: vol, fl, flm, flp, gm, gmm, flow, diff, acof
real(dp) :: reseq
nq = 1
call clearsor()
call gamsor()
do i = 3, njm1 ! south boundary
  fl = xcvi(i)*v(i,2)*den(i,1)
  flm = xcvip(i-1)*v(i-1,2)*den(i-1,1)
  flow = r(1)*(fl+flm)
  diff = r(1)*(xcvi(i)*gam(i,1)+xcvip(i-1)*gam(i-1,1))/ydif(2)
  acof = DAPEC(flow, diff, iadv) ! get D*A(Pe)
  as(i,2) = acof + max(0.0_dp, flow)
end do
do j = 2, njm1
  ! calculate aw
  flow = arx(j)*u(2,j)*den(1,j) ! west boundary for aw(3,j)
  diff = arx(j)*gam(1,j)/(sew(2)*sx(j))
  acof = DAPEC(flow, diff, iadv)
  aw(3,j) = acof + max(0.0_dp, flow)
  do i = 3, njm1
    if (i == njm1) then ! east boundary for aw(ni,j)
      flow = arx(j)*u(ni,j)*den(ni,j)
      diff = arx(j)*gam(ni,j)/(sew(nim1)*sx(j))
    else
      fl = u(i,j)*(fx(i)*den(i,j)+fxm(i)*den(i-1,j))
      flp = u(i+1,j)*(fx(i+1)*den(i+1,j)+fxm(i+1)*den(i,j))
      flow = arx(j)*0.5_dp*(fl+flp)
      diff = arx(j)*gam(i,j)/(sew(i)*sx(j))
    end if
    acof = DAPEC(flow, diff, iadv)
    aw(i+1,j) = acof + max(0.0_dp, flow)
    ae(i,j) = aw(i+1,j) - flow ! relationship
    if (j == njm1) then
      fl = xcvi(i)*v(i,nj)*den(i,nj)
      flm = xcvip(i-1)*v(i-1,nj)*den(i-1,nj)
      diff = r(nj)*(xcvi(i)*gam(i,nj)+xcvip(i-1)*gam(i-1,nj))/ydif(nj)
    else
      fl = xcvi(i)*v(i,j+1)*(fy(j+1)*den(i,j+1)+fym(j+1)*den(i,j))
      flm = xcvip(i-1)*v(i-1,j+1)*(fy(j+1)*den(i-1,j+1)+fym(j+1)*den(i-1,j))
      gm = gam(i,j)*gam(i,j+1)/(sns(j)*gam(i,j+1)+sns(j+1)*gam(i,j))*xcvi(i)
      gmm = gam(i-1,j)*gam(i-1,j+1)/(sns(j)*gam(i-1,j+1)+sns(j+1)*gam(i-1,j))*xcvip(i-1)
      diff = rmn(j+1)*2.0_dp*(gm+gmm)
    end if
    flow = rmn(j+1)*(fl+flm)
    acof = DAPEC(flow, diff, iadv)
    as(i,j+1) = acof + max(0.0_dp, flow)
    an(i,j) = as(i,j+1) - flow
    vol = ycvr(j)*sewu(i)
    ! source
    su(i,j) = su(i,j)*vol
    sp(i,j) = sp(i,j)*vol
    ! pressure gradient term
    du(i,j) = vol / (xdif(i)*sx(j)) ! linear interp for boundary pressure
    su(i,j) = su(i,j) + du(i,j)*(p(i-1,j)-p(i,j))
  end do
end do
! solve
resor(1) = 0.0_dp
do j = 2, njm1
  do i = 3, nim1
    ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j) - sp(i,j)
    reseq = ap(i,j)*u(i,j) - ae(i,j)*u(i+1,j) - aw(i,j)*u(i-1,j) &
          - an(i,j)*u(i,j+1) -as(i,j)*u(i,j-1) - su(i,j)
    resor(1) = resor(1) + abs(reseq)
    ap(i,j) = ap(i,j) / urf(1)
    su(i,j) = su(i,j) + (1-urf(1))*ap(i,j)*u(i,j)
    du(i,j) = du(i,j) / ap(i,j)
  end do
end do
do n = 1, nsweep(1)
  call lisolv(3, 2, ni, nj, u) 
end do
end subroutine calcu_

subroutine calcv()
use vars_mod
use case_mod, only: gamsor
implicit none
integer :: i, j, n
real(dp) :: arean, areas, areaw, areae, vol
real(dp) :: gne, gnw, gnn, gse, gsw, gss, gp, ge, gw, gn, gs
real(dp) :: ce, cw, cn, cs, cp, smp
real(dp) :: vise, visw, visn, viss, de, dw, dn, ds
real(dp) :: dudye, dudyw, dvdyn, dvdys
real(dp) :: reseq
nq = 2
call clearsor()
call gamsor()
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
    vise = 0.25_dp*(gam(i,j-1)+gam(i,j)+gam(i+1,j)+gam(i+1,j-1))
    visw = 0.25_dp*(gam(i-1,j-1)+gam(i-1,j)+gam(i,j)+gam(i,j-1))
    visn = gam(i,j)
    viss = gam(i,j-1)
    ! diff coef
    de = vise*areae/(x(i+1)-x(i))
    dw = visw*areaw/(x(i)-x(i-1))
    dn = visn*arean/(yv(j+1)-yv(j))
    ds = viss*areas/(yv(j)-yv(j-1))
    call getanb(i, j, iadv, ce, cw, cn, cs, de, dw, dn, ds, &
                 ae(i,j), aw(i,j), an(i,j), as(i,j))
    su(i,j) = su(i,j)*vol
    sp(i,j) = sp(i,j)*vol
    if (incl_visor) then
      ! viscous terms in source 
      dudye = (u(i+1,j)-u(i+1,j-1))/snsv(j)
      dudyw = (u(i,j)-u(i,j-1))/snsv(j)
      dvdyn = (v(i,j+1)-v(i,j))/sns(j)
      dvdys = (v(i,j)-v(i,j-1))/sns(j-1)
      su(i,j) = su(i,j) + (vise*dudye-visw*dudyw)/sew(i)*vol &
              + (visn*dvdyn-viss*dvdys)/snsv(j)*vol
    end if
    if (incl_falsor) then
      su(i,j) = su(i,j) + cp*v(i,j)
      sp(i,j) = sp(i,j) - cp
    end if
    dv(i,j) = vol/ydif(j)
    su(i,j) = su(i,j) + dv(i,j)*(p(i,j-1)-p(i,j))
  end do
end do
resor(2) = 0.0_dp
do j = 3, njm1
  do i = 2, nim1
    ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j) - sp(i,j)
    reseq = ap(i,j)*v(i,j) - ae(i,j)*v(i+1,j) - aw(i,j)*v(i-1,j) &
          - an(i,j)*v(i,j+1) - as(i,j)*v(i,j-1) - su(i,j)
    resor(2) = resor(2) + abs(reseq)
    ap(i,j) = ap(i,j)/urf(2)
    su(i,j) = su(i,j) + (1.0_dp-urf(2))*ap(i,j)*v(i,j)
    dv(i,j) = dv(i,j)/ap(i,j)
  end do
end do
do n = 1, nsweep(2)
  call lisolv(2, 3, ni, nj, v) 
end do
end subroutine calcv

subroutine calcp()
use vars_mod
use case_mod, only: gamsor
implicit none
integer :: i, j, n
real(dp) :: arean, areas, areaw, areae
real(dp) :: ge, gw, gn, gs, dene, denw, denn, dens
real(dp) :: ce, cw, cn, cs, smp
real(dp) :: de, dw, dn, ds, ppref
nq = 3
call clearsor()
call gamsor()
resor(3) = 0.0_dp
do j = 2, njm1
  do i = 2, nim1
    areae = sns(j)
    areaw = sns(j)
    arean = sew(i)
    areas = sew(i)
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
    resor(3) = resor(3) + abs(smp)
  end do
end do
do j = 2, njm1
  do i = 2, nim1
    ! assemble
    ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j) - sp(i,j)
  end do
end do
! solve
pp(:,:) = 0.0_dp
do n = 1, nsweep(3)
  call lisolv(2, 2, ni, nj, pp) 
end do
! update 
ppref = pp(ipref,jpref)
do j = 2, njm1
  do i = 2, nim1
    p(i,j) = p(i,j) + urf(3)*(pp(i,j)-ppref)
    if (i /= 2) u(i,j) = u(i,j) + du(i,j)*(pp(i-1,j)-pp(i,j))
    if (j /= 2) v(i,j) = v(i,j) + dv(i,j)*(pp(i,j-1)-pp(i,j))
  end do
end do
end subroutine calcp

subroutine calct()
use vars_mod
use case_mod
implicit none
integer :: i, j, n
real(dp) :: arean, areas, areaw, areae, vol
real(dp) :: ge, gw, gn, gs, dene, denw, denn, dens
real(dp) :: game, gamw, gamn, gams
real(dp) :: ce, cw, cn, cs, cp, smp
real(dp) :: de, dw, dn, ds, reseq
nq = 4
call clearsor()
call gamsor()
! nominal density: rho*hc
den(:,:) = den(:,:) * hc(:,:)
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
    ge = dene*u(i+1,j)
    gw = denw*u(i,j)
    gn = denn*v(i,j+1)
    gs = dens*v(i,j)
    ce = ge*areae
    cw = gw*areaw
    cn = gn*arean
    cs = gs*areas
    smp = ce - cw + cn - cs
    cp = max(0.0_dp, smp)
    game = 0.5_dp*(gam(i,j)+gam(i+1,j))
    gamw = 0.5_dp*(gam(i-1,j)+gam(i,j))
    gamn = 0.5_dp*(gam(i,j)+gam(i,j+1))
    gams = 0.5_dp*(gam(i,j-1)+gam(i,j))
    de = game*areae/(x(i+1)-x(i))
    dw = gamw*areaw/(x(i)-x(i-1))
    dn = gamn*arean/(y(j+1)-y(j))
    ds = gams*areas/(y(j)-y(j-1))
    call getanb(i, j, iadv, ce, cw, cn, cs, de, dw, dn, ds, &
                 ae(i,j), aw(i,j), an(i,j), as(i,j))
    su(i,j) = su(i,j)*vol + cp*t(i,j)
    sp(i,j) = sp(i,j)*vol - cp
  end do
end do
resor(4) = 0.0_dp
do j = 2, njm1
  do i = 2, nim1
    ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j) - sp(i,j)
    reseq = ap(i,j)*t(i,j) - ae(i,j)*t(i+1,j) - aw(i,j)*t(i-1,j) &
          - an(i,j)*t(i,j+1) - as(i,j)*t(i,j-1) - su(i,j)
    resor(4) = resor(4) + abs(reseq)
    ap(i,j) = ap(i,j)/urf(4)
    su(i,j) = su(i,j) + (1.0_dp-urf(4))*ap(i,j)*t(i,j)
  end do
end do
do n = 1, nsweep(4)
  call lisolv(2, 2, ni, nj, t) 
end do
! original density
den(:,:) = den(:,:) / hc(:,:)
end subroutine calct

subroutine getanb(i, j, iadv, ce, cw, cn, cs, de, dw, dn, ds, &
                  a_e, a_w, a_n, a_s)
! get coefficient                 
use vars_mod, only: dp, ni, nj
implicit none                 
integer, intent(in) :: i, j, iadv
real(dp), intent(in) :: ce, cw, cn, cs, de, dw, dn, ds
real(dp), intent(out) :: a_e, a_w, a_n, a_s
if (iadv == 1) then !CDS
  a_e = de - 0.5_dp*ce
  a_w = dw + 0.5_dp*cw
  a_n = dn - 0.5_dp*cn
  a_s = ds + 0.5_dp*cs
else if (iadv == 2) then ! UDS
  a_e = de + max(0.0_dp, -ce)
  a_w = dw + max(0.0_dp, cw)
  a_n = dn + max(0.0_dp, -cn)
  a_s = ds + max(0.0_dp, cs)
else ! Hybrid
  a_e = max(abs(0.5*ce), de) - 0.5*ce
  a_w = max(abs(0.5*cw), dw) + 0.5*cw
  a_n = max(abs(0.5*cn), dn) - 0.5*cn
  a_s = max(abs(0.5*cs), ds) + 0.5*cs
end if
end subroutine getanb

subroutine postproc()
use vars_mod
use case_mod
implicit none
integer :: i, j
if (lprint(1)) then
  ! centre u vel
  uc(1,:) = u(2,:)
  uc(ni,:) = u(ni,:)
  uc(2:nim1,:) = 0.5_dp*(u(2:nim1,:)+u(3:ni,:))
end if
if (lprint(2)) then
  ! centre v vel
  vc(:,1) = v(:,2)
  vc(:,nj) = v(:,nj)
  vc(:,2:njm1) = 0.5_dp*(v(:,2:njm1)+v(:,3:nj))
end if
if (lprint(3)) then
  ! extrapolate pressure to boundary
  p(1,:) = p(2,:) + (p(2,:)-p(3,:))*xdif(2)/xdif(3)
  p(ni,:) = p(nim1,:) + (p(nim1,:)-p(nim2,:))*xdif(ni)/xdif(nim1)
  p(:,1) = p(:,2) + (p(:,2)-p(:,3))*ydif(2)/ydif(3)
  p(:,nj) = p(:,njm1) + (p(:,njm1)-p(:,njm2))*ydif(nj)/ydif(njm1)
  ! corner pressure
  p(1,1) = p(2,1) + p(1,2) - p(2,2)
  p(1,nj) = p(1,njm1) + p(2,nj) - p(2,njm1)
  p(ni,1) = p(nim1,1) + p(ni,2) - p(nim1,2)
  p(ni,nj) = p(ni,njm1) + p(nim1,nj) - p(nim1,njm1)
  ! with respect to reference point
  p(:,:) = p(:,:) - p(ipref,jpref)
end if
if (lprint(4)) then
  ! corner temperature
  t(1,1) = t(2,1) + t(1,2) - t(2,2)
  t(1,nj) = t(1,njm1) + t(2,nj) - t(2,njm1)
  t(ni,1) = t(nim1,1) + t(ni,2) - t(nim1,2)
  t(ni,nj) = t(ni,njm1) + t(nim1,nj) - t(nim1,njm1)
end if
end subroutine postproc

subroutine lisolv(istart, jstart, ni, nj, phi)
! line by line tdma solver
use vars_mod, only: dp, ap, an, as, ae, aw, su, sp
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
    term = 1.0_dp/(d(j)-b(j)*a(j-1))
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

subroutine clearsor()
! clean source terms because of accumulation
use vars_mod
use case_mod
implicit none
su(:,:) = 0.0_dp
sp(:,:) = 0.0_dp
end subroutine clearsor

subroutine array_alloc()
! allocate arrays
use vars_mod
use case_mod
implicit none
allocate(x(1:ni), xu(1:ni), xdif(1:ni), sew(1:ni), sewu(1:ni), xcvi(1:ni), &
         xcvip(1:ni), fx(1:ni), fxm(1:ni), &
         stat=ierr, errmsg=errmsg)
if (ierr /= 0) write(*,*) 'ALLOACAE ERROR! ', errmsg
x(:) = 0.0_dp; xu(:) = 0.0_dp; xdif(:) = 0.0_dp;  sew(:) = 0.0_dp
sewu(:) = 0.0_dp; xcvi(:) = 0.0_dp; xcvip(:) = 0.0_dp
fx(:) = 0.0_dp; fxm(:) = 0.0_dp

allocate(y(1:nj), yv(1:nj), ydif(1:nj), sns(1:nj), snsv(1:nj), &
         ycvr(1:nj), ycvrs(1:nj), arx(1:nj), arxj(1:nj), arxjp(1:nj), &
         r(1:nj), rmn(1:nj), sx(1:nj), sxmn(1:nj), &
         fy(1:nj), fym(1:nj), fv(1:nj), fvp(1:nj), stat=ierr, errmsg=errmsg)
if (ierr /= 0) write(*,*) 'ALLOACAE ERROR! ', errmsg
y(:) = 0.0_dp; yv(:) = 0.0_dp; ydif(:) = 0.0_dp; sns(:) = 0.0_dp
snsv(:) = 0.0_dp; ycvr(:) = 0.0_dp; ycvrs(:) = 0.0_dp; ycvr(:) = 0.0_dp
ycvrs(:) = 0.0_dp; arx(:) = 0.0_dp; arxj(:) = 0.0_dp; arxjp(:) = 0.0_dp
r(:) = 0.0_dp; rmn(:) = 0.0_dp; sx(:) = 0.0_dp; sxmn(:) = 0.0_dp
fy(:) = 0.0_dp; fym(:) = 0.0_dp; fv(:) = 0.0_dp; fvp(:) = 0.0_dp
allocate(u(1:ni,1:nj), v(1:ni,1:nj), p(1:ni,1:nj), pp(1:ni, 1:nj), &
         den(1:ni,1:nj), du(1:ni,1:nj), dv(1:ni,1:nj), &
         stat=ierr, errmsg=errmsg)
if (ierr /= 0) write(*,*) 'ALLOCATE ERROR! ', errmsg
u(:,:) = 0.0_dp; v(:,:) = 0.0_dp; p(:,:) = 0.0_dp; pp(:,:) = 0.0_dp
den(:,:) = 0.0_dp; du(:,:) = 0.0_dp; dv(:,:) = 0.0_dp

allocate(t(1:ni,1:nj), gam(1:ni,1:nj), hc(1:ni,1:nj), stat=ierr, errmsg=errmsg)
if (ierr /= 0) write(*,*) 'ALLOCATE ERROR! ', errmsg
t(:,:) = 0.0_dp; gam(:,:) = 0.0_dp; hc(:,:) = 0.0_dp

allocate(uc(1:ni,1:nj), vc(1:ni,1:nj), stat=ierr, errmsg=errmsg)
if (ierr /= 0) write(*,*) 'ALLOCATE ERROR! ', errmsg
uc(:,:) = 0.0_dp; vc(:,:) = 0.0_dp

allocate(ae(1:ni,1:nj), aw(1:ni,1:nj), an(1:ni,1:nj), as(1:ni,1:nj), &
         ap(1:ni,1:nj), su(1:ni,1:nj), sp(1:ni,1:nj), &
          stat=ierr, errmsg=errmsg)
if (ierr /= 0) write(*,*) 'ALLOCATE ERROR! ', errmsg
ae(:,:) = 0.0_dp; aw(:,:) = 0.0_dp; an(:,:) = 0.0_dp; as(:,:) = 0.0_dp
ap(:,:) = 0.0_dp; su(:,:) = 0.0_dp; sp(:,:) = 0.0_dp
end subroutine array_alloc

subroutine array_dealloc()
! deallocate arrays
use vars_mod
use case_mod
implicit none
deallocate(x, xu, xdif, sew, sewu, xcvi, xcvip, fv, fvp, fx, fxm, &
           stat=ierr)
if (ierr /= 0) write(*,*) 'DEALLOCATE ERROR'
deallocate(y, yv, ydif, sns, snsv, ycvr, ycvrs, &
           arx, arxj, arxjp, r, rmn, sx, sxmn, &
           fy, fym, stat=ierr)
if (ierr /= 0) write(*,*) 'DEALLOCATE ERROR'
deallocate(u, v, p, pp, den, du, dv, uc, vc, stat=ierr)
if (ierr /= 0) write(*,*) 'DEALLOCATE ERROR'
deallocate(ae, aw, an, as, ap, su, sp, stat=ierr)
if (ierr /= 0) write(*,*) 'DEALLOCATE ERROR'
if (allocated(t)) deallocate(t, stat=ierr)
if (ierr /= 0) write(*,*) 'DEALLOCATE ERROR'
if (allocated(gam)) deallocate(gam, stat=ierr)
if (ierr /= 0) write(*,*) 'DEALLOCATE ERROR'
if (allocated(hc)) deallocate(hc, stat=ierr)
if (ierr /= 0) write(*,*) 'DEALLOCATE ERROR'
end subroutine array_dealloc
