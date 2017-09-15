!*****************************************************************************
! Conduction in Hollow Cylinder with Heat Generation

!                                Adiabatic
!                ________________________________________
!               |                |           |           |
!               |                |S=200-0.1T |           |
!               |                |    k=4    |           |
!               |                |___________|           |
!               |                                        | <-
!  T = T0(y+1)  |                                        | <- q'' = 50
!               |          T0 = 20                       | <-
!               |          k = 2                         |
!               |                                        |
!               |                                        |
!               |                                        |
!               |________________________________________| 
!                              -> -> -> -> ->      
!                              h = 50, Tf = 5
!                             

!  Length = 2.0
!  Thickness = 1.0
!  Inner radius = 1.0
!  Thermal conductivity = 2.0
!  initial temperature = 20
!  A heating region with a different conductivity (S=80-2T, k=1)

!  B.C.: 
!       left - prescribed temperature
!       right - constant heat flux
!       top - adiabatic
!       bottom - convective heat flux

!   Finite volume method, steady state simulation
!   Solves temperature only. 
!   Problem is specified in module user.

! Ruipengyu Li
! Modified: 14/09/2017
!*****************************************************************************
module vars
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
integer :: iter = 0, screen_iter = 100, max_iter = 10000
integer, dimension(max_nq) :: nsweep = 3 ! number of sweeps for TDMA
logical :: llast = .false.
logical :: ltec = .false., ltxt = .false.
logical, dimension(max_nq) :: lsolve = .false.
logical, dimension(max_nq) :: lprint = .false.
real(dp) :: xstart, xend, ystart, yend, rstart
real(dp) :: den_ini = 1.0_dp, cp_ini = 1.0_dp ! initial rho and cp
real(dp), dimension(max_nq) :: urf = 0.8_dp ! under-relaxation factor
real(dp), dimension(max_nq) :: resor = 0.0_dp ! residual
real(dp), allocatable, dimension(:) :: x, xu, xdif, sew, sewu, xcvi, xcvip
real(dp), allocatable, dimension(:) :: y, yv, ydif, sns, snsv, ycvr, ycvrs
real(dp), allocatable, dimension(:) :: arx, arxj, arxjp, r, rmn, sx, sxmn
real(dp), allocatable, dimension(:) :: fv, fvp, fx, fxm
real(dp), allocatable, dimension(:) :: fy, fym
real(dp), allocatable, dimension(:,:) :: u, v, p, pp, du, dv, t, uc, vc
real(dp), allocatable, dimension(:,:) :: den, gam, cp
real(dp), allocatable, dimension(:,:) :: ae, aw, an, as, ap, su, sp
end module vars

module user
! case of study
use vars
implicit none
integer :: ibnd = 2 ! boundary treatment 1:boundary value update 
!                                        2:additional source term
integer :: is1, is2, js1, js2 ! region of heating
real(dp) :: t0, gam0, gam1, htc, tf, qin, qwall
contains

subroutine control()
! mesh input
implicit none 
icrd = 2 ! 1:Cartesian 2:Axisym cylind 3:polar
! solution domain
ni = 52
nj = 52
xstart = 0.0_dp
xend = 2.0_dp
ystart = 0.0_dp
yend = 1.0_dp
rstart = 1.0_dp
! control params
lsolve(1:3) = .false.
lsolve(4) = .true.
lprint(1:3) = .false.
lprint(4) = .true.
ltec = .true.
ltxt = .false.
iadv = 1 ! 1:CDS 2:UDS 3:Hybrid 4:Power-law
max_iter = 19999 ! max iteration
screen_iter = 500 ! interval for screen print
urf(1) = 0.3_dp 
urf(2) = 0.3_dp
urf(3) = 0.9_dp
urf(4) = 1.0_dp
nsweep(1) = 3
nsweep(2) = 3
nsweep(3) = 3
nsweep(4) = 3
imon = ni / 2
jmon = nj / 2
end subroutine control

subroutine initial()
! initial values and unchanged b.c.
implicit none
integer :: i, j
is1 = ni/2 - 10
is2 = ni/2 + 10
js2 = nj
js1 = js2 - 20
gam0 = 3.0_dp ! thermal conductivity of the domain
gam1 = 1.0_dp ! thermal conductivity of the heating material
qin = 50.0_dp ! heat input at south boundary 
htc = 50.0_dp ! convective heat transfer coef of cooling fluid
t0 = 20.0_dp ! initial temperature
tf = 10.0_dp ! temperature of cooling fluid
t(2:,:) = t0
t(1,:) = t0 * (1.0_dp + y(:)) ! temperature at left boundary 
end subroutine initial 

subroutine density() 
! specify fluid density as function of temperature 
implicit none
end subroutine density

subroutine boundary() 
! update boundary
implicit none 
if (ibnd == 2 .and. .not.llast) return ! boundary only update at last
t(2:nim1,nj) = t(2:nim1,njm1) ! north
t(2:nim1,1) = (gam0*t(2:nim1,2) + htc*ydif(2)*tf) / (gam0 + htc*ydif(2)) ! south
t(ni,2:njm1) = t(nim1,2:njm1) + qin * xdif(ni) / gam0 ! east
end subroutine boundary 

subroutine gamsor() 
! set and update diffusion coefficients
! store sp and su of variables 
! set up additional source terms
! diffusion coefs 
implicit none
integer :: i, j
gam(:,:) = gam0
gam(is1:is2,js1:js2) = gam1
su(is1:is2,js1:js2) = 200.0_dp ! linear source term
sp(is1:is2,js1:js2) = -0.1_dp ! heat generation
if (ibnd == 1) return ! no further treatment, boundary value is updated
gam(2:nim1,nj) = 0.0_dp ! north, adiabatic
do i = 2, nim1 ! south
  su(i,2) = su(i,2) + r(1)/arx(2) * tf / (ydif(2)/gam(i,2) + 1.0_dp/htc)
  sp(i,2) = sp(i,2) - r(1)/arx(2) / (ydif(2)/gam(i,2) + 1.0_dp/htc)
  gam(i,1) = 0.0_dp
end do
su(nim1,2:njm1) = su(nim1,2:njm1) + qin / sew(nim1) ! east +q*A/V
gam(ni,2:njm1) = 0.0_dp ! east
end subroutine gamsor

subroutine output()
! output to screen and files.
implicit none
integer :: i, j
if (.not.llast) then
  if (mod(iter, screen_iter) == 0) then
    write(*,'(/,a,i6,4(1x,a,es10.3),/,11x,*(1x,a,es13.6))') &
            'Iter=', iter, 'URes=', resor(1), 'VRes=', resor(2), 'MRes=', resor(3), &
            'TRes=', resor(4), 'U=', u(imon,jmon), 'V=', v(imon,jmon), &
            'P=', p(imon,jmon), 'T=', t(imon,jmon)
  end if 
  ! this is a rather crude convergence criterion
  if (maxval(resor) < 1.0e-10_dp .or. iter == max_iter-1) llast = .true.
else
! write to screen
  qwall = 0.0_dp
  do j = 2, njm1
    qwall = qwall + gam(1,j)*(y(j)-y(j-1))*(t(1,j)-t(2,j))/(x(2)-x(1))
  end do
  write(*,'(/,a,/)') 'Conduction in Hollow Cylinder'
  write(*,*) 'Finite Volume Method'
  if (any(lsolve(:3))) then ! if flow is solved
    if (iadv == 1) then
      write(*,*) 'CDS for advection'
    else if (iadv == 3) then
      write(*,*) 'Hybrid for advection'
    else if (iadv == 4) then
      write(*,*) 'Power-law for advection'
    else 
      write(*,*) 'UDS for advection'
    end if
  end if
  write(*,'(/,2a,i3,3x,a,i3)') 'Grid: ', 'x-cv = ', ni-2, 'y-cv = ', nj-2
  write(*,'(/,a)') 'Values at monitor:'
  write(*,'(*(a,es10.3,3x))') 'x_mon =', x(imon), 'y_mon =', y(jmon)
  write(*,'(*(a,es18.10,3x))') 'u_mon =', uc(imon,jmon), &
          'v_mon =', vc(imon,jmon), 't_mon =', t(imon,jmon)
  write(*,'(a,1x,es12.5)') 'Qwall_west = ', qwall
  ! write to tecplot
  if (ltec) then
    open(unit=1, file='result.dat', status='replace')
    write(1,*) 'title="Natural Convection in a Square Cavity"'
    write(1,*) 'variables="x", "y", "Cond", "T"'
    write(1,'(a,2x,a,i3,2x,a,i3,2x,a)') 'zone', 'i=', ni, 'j=', nj, 'f=point'
    do j = 1, nj
      do i = 1, ni
        write(1,'(*(1x,es14.7))') x(i), y(j), &
                                  gam(i,j), t(i,j)
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

end module user

module serv
use vars, only: dp
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

end module serv
!*****************************************************************************
program main
! main program
use vars
use user
implicit none
call control()
call array_init()
call grid()
call initial()
do iter = 1, max_iter
  call boundary()
  if (lsolve(1)) then
    call calcu()
    call calcv()
    call calcp()
  end if
  if (lsolve(4)) call calct 
  if (llast) exit
  call output()
end do
call postproc()
call output()
call array_close()
stop
end program main
!*****************************************************************************
subroutine grid()
! Cell centred, backward staggered, uniform 
use vars
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
! initialise density and specific heat
den(:,:) = den_ini
cp(:,:) = cp_ini
end subroutine grid

subroutine calcu()
! u control volume using compass notation
use vars
use serv, only: DAPEC, ARIMEAN, HARMEAN
use user, only: gamsor
implicit none
integer :: i, j, n
real(dp) :: vis_ne, vis_nw ! viscosity at north corners of u-cv
real(dp) :: vol, fl, flm, flp, flow, diff, acof
real(dp) :: reseq
nq = 1
call clrsor()
call gamsor()
do i = 3, njm1 ! south boundary for as(i,2)
  fl = xcvi(i) * v(i,2) * den(i,1)
  flm = xcvip(i-1) * v(i-1,2) * den(i-1,1)
  flow = r(1) * (fl+flm)
  diff = r(1) * (xcvi(i)*gam(i,1)+xcvip(i-1)*gam(i-1,1)) / ydif(2)
  acof = DAPEC(flow, diff, iadv) ! get D*A(Pe)
  as(i,2) = acof + max(0.0_dp, flow)
end do
resor(1) = 0.0_dp
do j = 2, njm1
  ! calculate aw
  flow = arx(j) * u(2,j) * den(1,j) ! west boundary for aw(3,j)
  diff = arx(j) * gam(1,j) / (sew(2)*sx(j))
  acof = DAPEC(flow, diff, iadv)
  aw(3,j) = acof + max(0.0_dp, flow)
  do i = 3, njm1
    if (i == njm1) then ! east boundary for aw(ni,j)
      flow = arx(j) * u(ni,j) *den(ni,j)
      diff = arx(j) * gam(ni,j) / (sew(nim1)*sx(j))
    else
      fl = u(i,j) * ARIMEAN(den(i-1,j), den(i,j), fx(i))
      flp = u(i+1,j) * ARIMEAN(den(i,j), den(i+1,j), fx(i+1))
      flow = arx(j) * 0.5_dp * (fl+flp)
      diff = arx(j) * gam(i,j) / (sew(i)*sx(j))
    end if
    acof = DAPEC(flow, diff, iadv)
    aw(i+1,j) = acof + max(0.0_dp, flow)
    ae(i,j) = aw(i+1,j) - flow ! relationship
    if (j == njm1) then
      fl = xcvi(i) * v(i,nj) * den(i,nj)
      flm = xcvip(i-1) * v(i-1,nj) * den(i-1,nj)
      diff = r(nj) * (xcvi(i)*gam(i,nj)+xcvip(i-1)*gam(i-1,nj)) / ydif(nj)
    else
      fl = xcvi(i) * v(i,j+1) * ARIMEAN(den(i,j), den(i,j+1), fy(j+1))
      flm = xcvip(i-1) * v(i-1,j+1) * ARIMEAN(den(i-1,j), den(i-1,j+1), fy(j+1))
      vis_ne = HARMEAN(gam(i,j), gam(i,j+1), fy(j+1))
      vis_nw = HARMEAN(gam(i-1,j), gam(i-1,j+1), fy(j+1))
      diff = rmn(j+1) * (vis_ne*xcvi(i)+vis_nw*xcvip(i-1)) / ydif(j+1)
    end if
    flow = rmn(j+1) * (fl + flm)
    acof = DAPEC(flow, diff, iadv)
    as(i,j+1) = acof + max(0.0_dp, flow)
    an(i,j) = as(i,j+1) - flow
    vol = ycvr(j) * sewu(i)
    ! source
    su(i,j) = su(i,j) * vol
    sp(i,j) = sp(i,j) * vol
    ! pressure gradient term
    du(i,j) = vol / (xdif(i)*sx(j)) ! linear interp for boundary pressure
    su(i,j) = su(i,j) + du(i,j)*(p(i-1,j)-p(i,j))
    ! assemble and under-relax
    ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j) - sp(i,j)
    reseq = ap(i,j)*u(i,j) - ae(i,j)*u(i+1,j) - aw(i,j)*u(i-1,j) &
          - an(i,j)*u(i,j+1) -as(i,j)*u(i,j-1) - su(i,j)
    resor(1) = resor(1) + abs(reseq)
    ap(i,j) = ap(i,j) / urf(1)
    su(i,j) = su(i,j) + (1-urf(1))*ap(i,j)*u(i,j)
    du(i,j) = du(i,j) / ap(i,j) ! for pressure correction
  end do
end do
do n = 1, nsweep(1)
  call tdma(3, 2, ni, nj, u) 
end do
end subroutine calcu

subroutine calcv()
use vars
use serv, only: DAPEC, ARIMEAN, HARMEAN
use user, only: gamsor
implicit none
integer :: i, j, n
real(dp) :: vis_ne, vis_se, gm, gmm
real(dp) :: area, vol, fl, flm, flp, flow, diff, acof
real(dp) :: reseq
nq = 2
call clrsor()
call gamsor()
do i = 2, nim1
  ! calculate as(i,3) from south boundary
  area = r(1) * sew(i)
  flow = area * v(i,2) * den(i,1)
  diff = area * gam(i,1) * sns(2)
  acof = DAPEC(flow, diff, iadv)
  as(i,3) = acof + max(0.0_dp, flow)
end do
resor(2) = 0.0_dp
do j = 3, njm1
  fl = arxj(j) * u(2,j) * den(1,j)
  flm = arxjp(j-1) * u(2,j-1) * den(1,j-1)
  flow = fl + flm
  diff = (arxj(j)*gam(1,j) + arxjp(j-1)*gam(1,j-1)) / (xdif(2)*sxmn(j))
  acof = DAPEC(flow, diff, iadv)
  aw(2,j) = acof + max(0.0_dp, flow)
  do i = 2, nim1
    if (i == nim1) then ! east bound for aw(i+1)
      fl = arxj(j) * u(ni,j) * den(ni,j)
      flm = arxjp(j-1) * u(ni,j-1) * den(ni,j-1)
      diff = (arxj(j)*gam(ni,j) + arxjp(j-1)*gam(ni,j-1)) / (xdif(ni)*sxmn(j))
    else
      fl = arxj(j) * u(i+1,j) * ARIMEAN(den(i,j), den(i+1,j), fx(i+1))
      flm = arxjp(j-1) * u(i+1,j-1) * ARIMEAN(den(i,j-1), den(i+1,j-1), fx(i+1))
      vis_ne = HARMEAN(gam(i,j), gam(i+1,j), fx(i+1))
      vis_se = HARMEAN(gam(i,j-1), gam(i+1,j-1), fx(i+1))
      diff = (vis_ne*arxj(j) + vis_se*arxjp(j-1)) / (xdif(i+1)*sxmn(j))
    end if
    flow = fl + flm
    acof = DAPEC(flow, diff, iadv)
    aw(i+1,j) = acof + max(0.0_dp, flow)
    ae(i,j) = aw(i+1,j) - flow
    if (j == njm1) then
      area = r(nj) * sew(i) ! north bound
      flow = area * v(i,nj) * den(i,nj)
      diff = area * gam(i,nj) / sns(njm1)
    else
      area = r(j) * sew(i) ! north face of v-cv
      fl = v(i,j) * rmn(j)*sew(i) * ARIMEAN(den(i,j-1), den(i,j), fy(j))
      flp = v(i,j+1) * rmn(j+1)*sew(i) * ARIMEAN(den(i,j), den(i,j+1), fy(j+1))
      flow = ARIMEAN(fl, flp, fvp(j))
      diff = area * gam(i,j) / sns(j)
    end if
    acof = DAPEC(flow, diff, iadv)
    as(i,j+1) = acof + max(0.0_dp, flow)
    an(i,j) = as(i,j+1) - flow
    vol = ycvrs(j) * sew(i) ! v-cv volume
    ! source
    su(i,j) = su(i,j) * vol
    sp(i,j) = sp(i,j) * vol
    dv(i,j) = vol / ydif(j)
    su(i,j) = su(i,j) + dv(i,j)*(p(i,j-1)-p(i,j))
    ! assemble and under-relax
    ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j) - sp(i,j)
    reseq = ap(i,j)*v(i,j) - ae(i,j)*v(i+1,j) - aw(i,j)*v(i-1,j) &
          - an(i,j)*v(i,j+1) - as(i,j)*v(i,j-1) - su(i,j)
    resor(2) = resor(2) + abs(reseq)
    ap(i,j) = ap(i,j) / urf(2)
    su(i,j) = su(i,j) + (1.0_dp-urf(2))*ap(i,j)*v(i,j)
    dv(i,j) = dv(i,j) / ap(i,j)
  end do
end do
do n = 1, nsweep(2)
  call tdma(2, 3, ni, nj, v) 
end do
end subroutine calcv

subroutine calcp()
! solve pressure correction equation in main cv
! SIMPLE algorithm
use vars
use serv, only: DAPEC, ARIMEAN, HARMEAN
use user, only: gamsor
implicit none
integer :: i, j, n
real(dp) :: areaew, arean, areas, vol ! face areas and cell volume
real(dp) :: dene, denw, denn, dens ! face densities
real(dp) :: ce, cw, cn, cs ! convection coefs
real(dp) :: ppref ! reference pressure correction
nq = 3
call clrsor()
call gamsor()
do j = 2, njm1
  do i = 2, nim1
    vol = ycvr(j) * sew(i)
    sp(i,j) = sp(i,j) * vol
    su(i,j) = su(i,j) * vol ! extra sources
  end do
end do
do i = 2, nim1 ! as(i,2)
  areas = rmn(2) * sew(i) ! area of south bound
  cs = den(i,1) * v(i,2) * areas
  su(i,2) = su(i,2) + cs
  as(i,2) = 0.0_dp ! adiabatic-like
end do
resor(3) = 0.0_dp
do j = 2, njm1 ! aw(2,j)
  areaew = arx(j)
  cw = den(1,j) * u(2,j) * areaew
  su(2,j) = su(2,j) + cw
  aw(2,j) = 0.0_dp
  do i = 2, nim1
    if (i == nim1) then ! ae(nim1,j)
      ce = den(ni,j) * u(ni,j) * areaew
      su(i,j) = su(i,j) - ce
      ae(i,j) = 0.0_dp
    else
      dene = ARIMEAN(den(i,j), den(i+1,j), fx(i+1))
      ce = den(i+1,j) * u(i+1,j) * areaew
      su(i,j) = su(i,j) - ce
      ae(i,j) = dene * areaew * du(i+1,j)
      su(i+1,j) = su(i+1,j) + ce
      aw(i+1,j) = ae(i,j) ! relationship
    end if
    arean = rmn(j+1) * sew(i)
    if (j == njm1) then ! an(i,njm1)
      cn = den(i,nj) * v(i,nj) * arean
      su(i,j) = su(i,j) - cn
      an(i,j) = 0.0_dp
    else
      denn = ARIMEAN(den(i,j), den(i,j+1), fy(j+1))
      cn = denn * v(i,j+1) * arean
      su(i,j) = su(i,j) - cn
      an(i,j) = denn * arean * dv(i,j+1)
      su(i,j+1) = su(i,j+1) + cn
      as(i,j+1) = an(i,j) ! relationship
    end if
    ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j) - sp(i,j)
    resor(3) = resor(3) + abs(su(i,j))
    pp(i,j) = 0.0_dp ! initial pressure correction field
  end do
end do
do n = 1, nsweep(3)
  call tdma(2, 2, ni, nj, pp) 
end do
! update u, v, p
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
use vars
use serv, only: DAPEC, ARIMEAN, HARMEAN
use user, only: gamsor
implicit none
integer :: i, j, n
real(dp) :: areaew, arean, areas, vol
real(dp) :: dene, denw, denn, dens
real(dp) :: game, gamw, gamn, gams
real(dp) :: flow, diff, acof
real(dp) :: reseq
nq = 4
call clrsor()
call gamsor()
if (lsolve(4)) den(:,:) = den(:,:) * cp(:,:) ! nominal density: rho*hc
do i = 2, nim1 ! as(i,2)
  areas = r(1) * sew(i) 
  flow = den(i,1) * v(i,2) * areas
  diff = gam(i,1) * areas / ydif(2)
  acof = DAPEC(flow, diff, iadv)
  as(i,2) = acof + max(0.0_dp, flow)
end do
resor(4) = 0.0_dp
do j = 2, njm1
  areaew = arx(j)
  denw = den(1,j)
  flow = denw * u(2,j) * areaew
  diff = gam(1,j) * areaew / (xdif(2)*sx(j))
  acof = DAPEC(flow, diff, iadv)
  aw(2,j) = acof + max(0.0_dp, flow)
  do i = 2, nim1
    if (i == nim1) then
      dene = den(nj,j)
      flow = dene * u(ni,j) * areaew
      diff = gam(ni,j) * areaew / (xdif(ni)*sx(j))
    else
      dene = ARIMEAN(den(i,j), den(i+1,j), fx(i+1))
      flow = dene * u(i+1,j) * areaew
      game = HARMEAN(gam(i,j), gam(i+1,j), fx(i+1))
      diff = game * areaew / (xdif(i+1)*sx(j))
    end if
    acof = DAPEC(flow, diff, iadv)
    aw(i+1,j) = acof + max(0.0_dp, flow)
    ae(i,j) = aw(i+1,j) - flow
    arean = rmn(j+1) * sew(i)
    if (j == njm1) then
      denn = den(i,nj)
      flow = denn * v(i,nj) * arean
      diff = gam(i,nj) * arean / ydif(nj)
    else
      denn = ARIMEAN(den(i,j), den(i,j+1), fy(j+1))
      flow = denn * v(i,j+1) * arean
      gamn = HARMEAN(gam(i,j), gam(i,j+1), fy(j+1))
      diff = gamn * arean / ydif(j+1)
    end if
    acof = DAPEC(flow, diff, iadv)
    as(i,j+1) = acof + max(0.0_dp, flow)
    an(i,j) = as(i,j+1) - flow
    vol = ycvr(j) * sew(i)
    ! source
    sp(i,j) = sp(i,j) * vol
    su(i,j) = su(i,j) * vol
    ! assemble and under-relax
    ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j) - sp(i,j)
    reseq = ap(i,j)*t(i,j) - ae(i,j)*t(i+1,j) - aw(i,j)*t(i-1,j) &
          - an(i,j)*t(i,j+1) - as(i,j)*t(i,j-1) - su(i,j)
    resor(4) = resor(4) + abs(reseq)
    ap(i,j) = ap(i,j)/urf(4)
    su(i,j) = su(i,j) + (1.0_dp-urf(4))*ap(i,j)*t(i,j)
  end do
end do
do n = 1, nsweep(4)
  call tdma(2, 2, ni, nj, t) 
end do
if (lsolve(4)) den(:,:) = den(:,:) / cp(:,:) ! recover density
end subroutine calct

subroutine postproc()
use vars
implicit none
integer :: i, j
if (lprint(1)) then ! centre u vel
  uc(1,:) = u(2,:)
  uc(ni,:) = u(ni,:)
  uc(2:nim1,:) = 0.5_dp*(u(2:nim1,:)+u(3:ni,:))
end if
if (lprint(2)) then ! centre v vel
  vc(:,1) = v(:,2)
  vc(:,nj) = v(:,nj)
  vc(:,2:njm1) = 0.5_dp*(v(:,2:njm1)+v(:,3:nj))
end if
if (lprint(3)) then ! pressure
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
if (lprint(4)) then ! corner temperature
  t(1,1) = t(2,1) + t(1,2) - t(2,2)
  t(1,nj) = t(1,njm1) + t(2,nj) - t(2,njm1)
  t(ni,1) = t(nim1,1) + t(ni,2) - t(nim1,2)
  t(ni,nj) = t(ni,njm1) + t(nim1,nj) - t(nim1,njm1)
end if
end subroutine postproc

subroutine tdma(istart, jstart, ni, nj, phi)
! line by line tdma solver
use vars, only: dp, ap, an, as, ae, aw, su, sp
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
end subroutine tdma

subroutine clrsor()
! clean source terms for accumulation
use vars, only: dp, sp, su
implicit none
sp(:,:) = 0.0_dp
su(:,:) = 0.0_dp
end subroutine clrsor

subroutine array_init()
! allocate arrays
use vars
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

allocate(t(1:ni,1:nj), gam(1:ni,1:nj), cp(1:ni,1:nj), stat=ierr, errmsg=errmsg)
if (ierr /= 0) write(*,*) 'ALLOCATE ERROR! ', errmsg
t(:,:) = 0.0_dp; gam(:,:) = 0.0_dp; cp(:,:) = 0.0_dp

allocate(uc(1:ni,1:nj), vc(1:ni,1:nj), stat=ierr, errmsg=errmsg)
if (ierr /= 0) write(*,*) 'ALLOCATE ERROR! ', errmsg
uc(:,:) = 0.0_dp; vc(:,:) = 0.0_dp

allocate(ae(1:ni,1:nj), aw(1:ni,1:nj), an(1:ni,1:nj), as(1:ni,1:nj), &
         ap(1:ni,1:nj), su(1:ni,1:nj), sp(1:ni,1:nj), &
          stat=ierr, errmsg=errmsg)
if (ierr /= 0) write(*,*) 'ALLOCATE ERROR! ', errmsg
ae(:,:) = 0.0_dp; aw(:,:) = 0.0_dp; an(:,:) = 0.0_dp; as(:,:) = 0.0_dp
ap(:,:) = 0.0_dp; su(:,:) = 0.0_dp; sp(:,:) = 0.0_dp
end subroutine array_init

subroutine array_close()
! deallocate arrays
use vars
implicit none
deallocate(x, xu, xdif, sew, sewu, xcvi, xcvip, fv, fvp, fx, fxm, &
           stat=ierr)
if (ierr /= 0) write(*,*) 'DEALLOCATE ERROR'
deallocate(y, yv, ydif, sns, snsv, ycvr, ycvrs, &
           arx, arxj, arxjp, r, rmn, sx, sxmn, &
           fy, fym, stat=ierr)
if (ierr /= 0) write(*,*) 'DEALLOCATE ERROR'
deallocate(u, v, p, pp, den, gam, du, dv, uc, vc, stat=ierr)
if (ierr /= 0) write(*,*) 'DEALLOCATE ERROR'
deallocate(ae, aw, an, as, ap, su, sp, stat=ierr)
if (ierr /= 0) write(*,*) 'DEALLOCATE ERROR'
if (allocated(t)) deallocate(t, stat=ierr)
if (ierr /= 0) write(*,*) 'DEALLOCATE ERROR'
if (allocated(cp)) deallocate(cp, stat=ierr)
if (ierr /= 0) write(*,*) 'DEALLOCATE ERROR'
end subroutine array_close
