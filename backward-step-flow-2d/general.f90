!**********************************************************************
!  2-d turbulent recirculating flows 
!  General subroutines for teacht
!  Ruipengyu Li April 2016
!**********************************************************************
module general
!common
  integer :: it, jt, ni, nj, nim1, njm1, nim2, njm2
  integer :: niter, iread, ifile
  real :: great
!headings
  character (len=6) :: hedu, hedv, hedp, hedt, hedk
  character (len=6) :: hedd, hedm, heda, hedb
!u-velocity
  integer :: nswpu
  real :: resoru, urfu
  real, allocatable :: dxepu(:), dxpwu(:), sewu(:)
!v-velocity
  integer :: nswpv
  real :: resorv, urfv
  real, allocatable :: dynpv(:), dypsv(:), snsv(:), rcv(:)
!pressure correction
  integer :: ipref, jpref, nswpp
  real :: resorm, urfp
  real, allocatable :: du(:, :), dv(:, :)
!ten
  integer :: nswpk
  real :: resork, urfk
!tdis
  integer :: nswpd
  real :: resore, urfe
!variables
  real, allocatable :: u(:, :), v(:, :), p(:, :), pp(:, :)
  real, allocatable :: te(:, :), ed(:, :)
  real, allocatable :: uc(:, :), vc(:, :)
!geometry
  integer :: indcos
  real, allocatable :: sns(:), sew(:), dxep(:), dxpw(:), dynp(:), dyps(:)
  real, allocatable :: x(:), y(:), xu(:), yv(:), r(:), rv(:)
!fluid properties
  real :: urfvis, viscos, densit, prandt
  real, allocatable :: den(:, :), vis(:, :)
!case 
  real :: uin, tein, edin, flowin, alamda
  real :: rsmall, rlarge, al1, al2
  real :: jstep, istep, jstp1, jstm1, istp1, istm1
  real :: jexit, iexit, jexp1, jexm1, iexp1, iexm1
!turbulence 
  real :: cd, cmu, c1, c2, cappa, elog, pred, prte
  real, allocatable :: gen(:, :), sukd(:, :), spkd(:, :)
!wall function 
  real, allocatable :: yplusn(:), xplusw(:), taun(:), tauw(:)
  real, allocatable :: ypluss(:), xpluse(:), taus(:), taue(:)
!coefficients of discretised equations
  real, allocatable :: ap(:, :), an(:, :), as(:, :), ae(:, :), aw(:, :)
  real, allocatable :: su(:, :), sp(:, :)
!heat transfer
  integer :: nswpt
  real :: resort, urft, prandl, pfun
  real :: tin, twalln, twalls, twalle, twallw
  real, allocatable :: t(:, :), gamh(:, :)
!criteria
  logical incalu, incalv, incalp, incalk, incald
  logical incalm, incala, incalb, incalt, inpro
end module general
!**********************************************************************

!**********************************************************************     
subroutine init()
!----------------------------------------------------------------------
!calculate coordinates of velocity locations and cell dimensions
!----------------------------------------------------------------------
  use general
  implicit none
  integer :: i, j
  do j = 1, nj
    r(j) = y(j)
    if (indcos==1) r(j) = 1.0
  end do
  dxpw(1) = 0.0
  dxep(ni) = 0.0
  do i = 1, nim1
    dxep(i) = x(i+1) - x(i)
    dxpw(i+1) = dxep(i)
  end do
  dyps(1) = 0.0
  dynp(nj) = 0.0
  do j = 1, njm1
    dynp(j) = y(j+1) - y(j)
    dyps(j+1) = dynp(j)
  end do
  sew(1) = 0.0
  sew(ni) = 0.0
  do i = 2, nim1
    sew(i) = 0.5*(dxep(i)+dxpw(i))
  end do
  sns(1) = 0.0
  sns(nj) = 0.0
  do j = 2, njm1
    sns(j) = 0.5*(dynp(j)+dyps(j))
  end do
  xu(1) = 0.0
  do i = 2, ni
    xu(i) = 0.5*(x(i)+x(i-1))
  end do
  dxpwu(1) = 0.0
  dxpwu(2) = 0.0
  dxepu(1) = 0.0
  dxepu(ni) = 0.0
  do i = 2, nim1
    dxepu(i) = xu(i+1) - xu(i)
    dxpwu(i+1) = dxepu(i)
  end do
  sewu(1) = 0.0
  sewu(2) = 0.0
  do i = 3, nim1
    sewu(i) = 0.5*(dxepu(i)+dxpwu(i))
  end do
  yv(1) = 0.0
  rv(1) = 0.0
  do j = 2, nj
    rv(j) = 0.5*(r(j)+r(j-1))
    rcv(j) = 0.5*(rv(j)+rv(j-1))
    yv(j) = 0.5*(y(j)+y(j-1))
  end do
  dypsv(1) = 0.0
  dypsv(2) = 0.0
  dynpv(nj) = 0.0
  do j = 2, njm1
    dynpv(j) = yv(j+1) - yv(j)
    dypsv(j+1) = dynpv(j)
  end do
  snsv(1) = 0.0
  snsv(2) = 0.0
  snsv(nj) = 0.0
  do j = 3, njm1
    snsv(j) = 0.5*(dynpv(j)+dypsv(j))
  end do

!----------------------------------------------------------------------
! set variables to zero except den() and vis()
!----------------------------------------------------------------------
  do i = 1, ni
    do j = 1, nj
      u(i, j) = 0.0
      v(i, j) = 0.0
      p(i, j) = 0.0
      pp(i, j) = 0.0
      t(i, j) = tin
      te(i, j) = 0.0
      ed(i, j) = 0.0
      den(i, j) = densit
      vis(i, j) = viscos
      gamh(i, j) = viscos/prandl
      du(i, j) = 0.0
      dv(i, j) = 0.0
      su(i, j) = 0.0
      sp(i, j) = 0.0
    end do
  end do
end subroutine init
!**********************************************************************

!**********************************************************************
subroutine props()
!----------------------------------------------------------------------
! calculate fluid properties
! calculate effective viscocity
! note the provision for under-relaxation (urfvis)
! exchange coefficients are calculated in calcte and calced
!----------------------------------------------------------------------
  use general
  implicit none
  integer :: i, j
  real :: visold
  do i = 2, nim1
    vis(i, 1) = vis(i, 2)
    do j = 2, njm1
      visold = vis(i, j)
      if (ed(i,j)==0.) then
        vis(i, j) = viscos
      else
        vis(i, j) = den(i, j)*te(i, j)**2*cmu/ed(i, j) + viscos
      end if
!-----under-relax viscosity
      vis(i, j) = urfvis*vis(i, j) + (1.-urfvis)*visold
      gamh(i, j) = viscos/prandl + (vis(i,j)-viscos)/prandt
    end do
  end do
end subroutine props
!**********************************************************************

!**********************************************************************
subroutine print(istart, jstart, ni, nj, it, jt, x, y, phi, head, ifile)
!-----printing field values of variables

  dimension phi(it, jt), x(1), y(1), store(50), head(6)
  dimension f(7), f4(11)
  data f/4h(1h , 4h,a6,, 4hi3, , 4h11i , 4h10, , 4h7x, , 4ha6) /
  data f4/4h 1i , 4h 2i , 4h 3i , 4h 4i , 4h 5i , 4h 6i , 4h 7i , 4h 8i , 4h 9i , 4h10i , 4h11i /
  data hi, hy/4hi = , 4hy = /
  ista = istart - 11
  ista = ista + 11
  iend = ista + 10
  iend = min0(ni, iend)
  iel = max0(iend-ista, 1)
  f(4) = f4(iel)

  open (unit=4, file='output/uvel.dat')
  open (unit=5, file='output/vvel.dat')
  open (unit=6, file='output/press.dat')
  open (unit=7, file='output/temp.dat')
  open (unit=8, file='output/tke.dat')
  open (unit=9, file='output/ed.dat')
  open (unit=10, file='output/visc.dat')

  if (ista/=istart) write (ifile, 115)

  do jj = jstart, nj
    j = jstart + nj - jj
    do i = istart, ni
      store(i) = phi(i, j)
    end do
    write (ifile, 113) y(j), (store(i), i=istart, ni)
  end do
  write (ifile, 114)(x(i), i=istart, ni)
  return
  110 format ('0', 20('*-'), 6a6, 20('-*'))
  112 format ('  j')
  113 format (f7.5, 1p20e10.2)
  114 format ('0.00000 ', f8.5, 20f10.5)
  115 format (///)
end subroutine print
!**********************************************************************


!**********************************************************************
subroutine tecprint(istart, jstart, ni, nj, it, jt, x, y, phi, head, ifile)
  implicit none
  integer :: i, j, istart, jstart, ifile, ni, nj, it, jt
  integer :: nim1, njm1
  character (len=6) :: head
  real :: x(ni), y(nj)
  real :: phi(ni, nj)
  real :: store(ni*nj)
  open (unit=4, file='output/uvel.dat')
  open (unit=5, file='output/vvel.dat')
  open (unit=6, file='output/press.dat')
  open (unit=7, file='output/temp.dat')
  open (unit=8, file='output/tke.dat')
  open (unit=9, file='output/ed.dat')
  open (unit=10, file='output/visc.dat')
  nim1 = ni - 1
  njm1 = nj - 1
  write (ifile, 120) nim1, njm1
  do j = jstart, nj
    do i = istart, ni
      store(i) = phi(i, j)
      write (ifile, 121) x(i), y(j), store(i)
    end do
  end do
  120 format (/'ZONE', 3x, 'I=', i3, 3x, 'j=', i3, 3x, 'F=POINT')
  121 format (1x, es10.3, 3x, es10.3, 3x, es10.3)
end subroutine tecprint
!**********************************************************************

!**********************************************************************
subroutine velprint(istart, jstart, ni, nj, it, jt, x, y, uvel, vvel, ifile)
  implicit none
  integer :: i, j, istart, jstart, ifile, ni, nj, it, jt
  integer :: nim1, njm1
  character (len=6) :: head
  real :: x(ni), y(nj)
  real :: uvel(ni, nj), vvel(ni, nj)
  real :: store2(ni*nj), store3(ni*nj)
  open (unit=14, file='output/vel.dat')
  nim1 = ni - 1 
  njm1 = nj - 1
  write (ifile, 120) nim1, njm1
  do j = jstart, nj
    do i = istart, ni
      store2(i) = uvel(i, j)
      store3(i) = vvel(i, j)
      write (ifile, 121) x(i), y(j), store2(i), store3(i)
    end do
  end do
  120 format (/'ZONE', 3x, 'I=', i3, 3x, 'j=', i3, 3x, 'F=POINT')
  121 format (1x, es10.3, 3x, es10.3, 3x, es10.3, 3x, es10.3)
end subroutine velprint
  
!**********************************************************************

!**********************************************************************
subroutine lisolv(istart, jstart, ni, nj, it, jt, phi)
  use general, only: ap, an, as, ae, aw, su, sp
  implicit none
  integer, intent (in) :: istart, jstart, ni, nj, it, jt
  real, dimension (it, jt), intent (out) :: phi
  real, dimension (2*it) :: a, b, c, d
  integer :: i, j, jj, nim1, njm1, jsrm1
  real :: term
  nim1 = ni - 1
  njm1 = nj - 1
  jsrm1 = jstart - 1
  a(jsrm1) = 0.0
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
      term = 1./(d(j)-b(j)*a(j-1))
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
!**********************************************************************     

!**********************************************************************
subroutine centred_vels()
!----------------------------------------------------------------------
!  obtain velocity values at cell centre
!----------------------------------------------------------------------
  use general
  implicit none
  integer :: i, j
!-----interpolate u-velocity	
  do i = 3, nim1
    do j = 2, nj
      uc(i, j) = u(i, j) + 0.5*(u(i+1,j)-u(i,j))
    end do
  end do
!-----interpolate v-velocity
  do i = 2, ni
    do j = 3, njm1
      vc(i, j) = v(i, j) + 0.5*(v(i,j)-v(i,j+1))
    end do
  end do
end subroutine centred_vels
!**********************************************************************

!**********************************************************************
subroutine closefile()
  use general, only: ifile
  implicit none
  do ifile = 4, 10
    close (unit=ifile)
  end do
  close (unit=14)
end subroutine closefile
!**********************************************************************

!**********************************************************************
subroutine writeout()
  use general
  implicit none
  call centred_vels
  if (incalu) call tecprint(2, 2, ni, nj, it, jt, xu, y, u, hedu, 4)
  if (incalv) call tecprint(2, 2, ni, nj, it, jt, x, yv, v, hedv, 5)
  if (incalp) call tecprint(2, 2, ni, nj, it, jt, x, y, p, hedp, 6)
  if (incalt) call tecprint(2, 2, ni, nj, it, jt, x, y, t, hedt, 7)
  if (incalk) call tecprint(2, 2, ni, nj, it, jt, x, y, te, hedk, 8)
  if (incald) call tecprint(2, 2, ni, nj, it, jt, x, y, ed, hedd, 9)
  if (incalu .and. incalv) call velprint(2, 2, ni, nj, it, jt, x, y, &
                                         uc, vc, 14)
  call closefile
end subroutine writeout
!**********************************************************************

