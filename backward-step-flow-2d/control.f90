!**********************************************************************
!  2-d planar or cylindrical turbulent flows
!  Teach-t in Fortran 95
!  Modified by Ruipengyu Li
!  April 2016 
!**********************************************************************      
subroutine contro()

  use general
  use promod
  use calc

!chapter  0  0  0  0  0  0  0  0  preliminaries  0  0  0  0  0  0  0  0

  great = 1.e30
  istep = 2
  niter = 1
!----------------------------------------------------------------------
!  modify the numbers for different grid
!  i-index of maximum dimension of dependent variables 
!  match the number when plot in column   
  it = 100
  jt = 100
  ni = 100
  nj = 100
!----------------------------------------------------------------------
  nswpu = 3
  nswpv = 3
  nswpp = 5
  nswpk = 3
  nswpd = 3
  nswpt = 3
!  i/o from files 
  open (unit=1, file='input/head.in')
  open (unit=2, file='output/prop.out')
  open (unit=3, file='output/stress.out')
  open (unit=11, file='output/kplus.dat')
  open (unit=12, file='output/lplus.dat')
  open (unit=13, file='output/residual.out')

  do iread = 1, 3
    read (1, *)
  end do
  read (1, 010) hedu, hedv, hedp, hedt, hedk, hedd, hedm, heda, hedb
!
!chapter  1  1  1  1  1  parameters and control indices  1  1  1  1  1  1
!
!-----allocate arrays
  allocate (u(it,jt), v(it,jt), p(it,jt), pp(it,jt), du(it,jt), &
            dv(it,jt))
  allocate (te(it,jt), ed(it,jt), gen(it,jt), sukd(it,jt), spkd(it,jt))
  allocate (uc(it,jt), vc(it,jt))
  allocate (den(it,jt), vis(it,jt), t(it,jt), gamh(it,jt))
  allocate (ap(it,jt), an(it,jt), as(it,jt), ae(it,jt), aw(it,jt), &
            su(it,jt), sp(it,jt))
  allocate (dxepu(2*it), dxpwu(2*jt), sewu(2*it))
  allocate (dynpv(2*jt), dypsv(2*jt), snsv(2*jt), rcv(2*jt))
  allocate (sns(2*jt), sew(2*it), dxep(2*it), dxpw(2*it), dynp(2*jt), &
            dyps(2*jt))
  allocate (x(2*it), y(2*jt), xu(2*it), yv(2*jt), r(2*jt), rv(2*jt))
  allocate (yplusn(2*jt), xplusw(2*it), taun(2*jt), tauw(2*it))
  allocate (ypluss(2*jt), xpluse(2*it), taus(2*jt), taue(2*it))
!-----grid
  nim1 = ni - 1
  njm1 = nj - 1
  indcos = 2
  jstep = 8
  istp1 = istep + 1
  istm1 = istep - 1
  jstp1 = jstep + 1
  jstm1 = jstep - 1
  rlarge = 0.0124
  altot = 16.0*rlarge
  epsx = 1.05 !1.2
  sumx = 0.5*epsx**(ni-4) + (epsx**(ni-3)-1.)/(epsx-1.) + 0.5
  dx = altot/sumx
  x(1) = -0.5*dx
  x(2) = -x(1)
  do i = 3, nim1
    x(i) = x(i-1) + dx
    dx = epsx*dx
  end do
  x(ni) = x(nim1) - x(ni-2) + x(nim1)
  al1 = 0.5*(x(istep)+x(istm1))
  al2 = altot - al1
  dy = rlarge/float(nj-2)
  y(1) = -0.5*dy
  do j = 2, nj
    y(j) = y(j-1) + dy
  end do
  rsmall = 0.5*(y(jstep)+y(jstp1))
!-----dependent variable selection
  incalu = .true.
  incalv = .true.
  incalp = .true.
  incalk = .true.
  incald = .true.
  inpro = .true.
  incalt = .true.
!-----fluid properties
  densit = 1000.0
  prandl = 7.0
!-----turbulence constants
  cmu = 0.09
  cd = 1.00
  c1 = 1.44
  c2 = 1.92
  cappa = .4187
  elog = 9.793
  pred = cappa*cappa/(c2-c1)/(cmu**.5)
  prte = 1.0
  prandt = 0.9
  pfun = prandl/prandt
  pfun = 9.24*(pfun**0.75-1.0)*(1.0+0.28*exp(-0.007*pfun))
!-----boundary values
  uin = 0.5
  ularge = uin*(rsmall/rlarge)**2
  tin = 0.0
  twalln = 100.0
  turbin = .03
  tein = turbin*uin**2
  alamda = 0.005
  edin = tein**1.5/(alamda*rlarge)
  viscos = 0.8e-3
!-----pressure calculation
  ipref = 2
  jpref = 2
!-----program control and monitor
  maxit = 500
  imon = 10
  jmon = 10
  urfu = 0.5
  urfv = 0.5
  urfp = 1.0
  urfe = 0.7
  urfk = 0.7
  urft = 1.0
  urfvis = 0.7
!  number of iterations to be completed for each intermediate step  
  indpri = 2
  sormax = 0.005
!
!chapter  2  2  2  2  2  2  initial operations  2  2  2  2  2  2  2  2  2
!
!-----calculate geometrical quantities and set variables to zero
  call init
!-----initialise variable fields
  flowin = 0.0
  arden = 0.0
  do j = 2, jstep
    u(2, j) = uin
    t(1, j) = tin
    te(1, j) = tein
    ed(1, j) = edin
    arden = 0.5*(den(1,j)+den(2,j))*r(j)*sns(j)
    flowin = flowin + arden*u(2, j)
  end do
  jfin = jstp1
  do i = 2, ni
    if (i>=istp1) jfin = nj
    factor = (yv(jstp1)*rv(jstp1))/(yv(jfin)*rv(jfin))
    jend = jfin - 1
    do j = 2, jend
      te(i, j) = tein
      ed(i, j) = edin
      u(i, j) = uin*factor
    end do
  end do
  do i = 2, nim1
    t(i, nj) = twalln
!  yplus is calculated in routine modte   
!  initial yplus
    yplusn(i) = 11.0
  end do
  do j = jstep, nj
    t(1, j) = tin
    xplusw(j) = 11.0
    if (j==jstep) xplusw(j) = 0.0
  end do
  call props
!-----initial output
  write (2, 210)
  write (2, 211)
  write (2, 220) uin
  write (2, 221) tin
  rsdrl = rsmall/rlarge
  re = uin*densit*2.0*rsmall/viscos
  write (2, 230) re
  write (2, 223) prandl
  write (2, 240) rsdrl
  write (2, 260) densit
  write (2, 250) viscos
  write (2, 222) twalln
  call writeout
!.......calculate residual sources normalization factors................
  flowin = 0.0
  xmonin = 0.0
  do j = 2, njm1
    arden = 0.5*(den(1,j)+den(2,j))*r(j)*sns(j)
    flowin = flowin + arden*u(2, j)
    xmonin = xmonin + arden*u(2, j)*u(2, j)
  end do
  snormt = flowin*abs(tin-twalln)
  resort = 0.0
!
!chapter  3  3  3  3  3  3  3  iteration loop  3  3  3  3  3  3  3  3  3
!
  write (13, 310) imon, jmon
  do
    if (incalu) call calcu
    if (incalv) call calcv
    if (incalp) call calcp
    if (incalk) call calcte
    if (incald) call calced
    if (incalt) call calct
!-----update fluid properities
    if (inpro) call props
!-----intermediate output
    resorm = resorm/flowin
    resoru = resoru/xmonin
    resorv = resorv/xmonin
    resort = resort/snormt
    dummy = t(imon, jmon)
    write (13, 311) niter, resoru, resorv, resorm, resort, resork, &
                   resore, u(imon, jmon), v(imon, jmon), p(imon, jmon), &
                   dummy, te(imon, njm1), ed(imon, njm1)
    if (mod(niter,indpri)/=0) then
      sorce = max(resorm, resoru, resorv, resort)
    else
      call writeout
      write (13, 310) imon, jmon
    end if
!-----termination tests  
    if (niter==999 .and. sorce>1.0e4*sormax) exit
    if (niter==maxit) exit
    if (sorce<=sormax) exit
    niter = niter + 1
  end do
!
!chapter  4  4  4  4  4  4  final operations and output  4  4  4  4  4  4
!
  call writeout
!-----calculation of non dimensional turbulence energy and length scale
  do i = 2, nim1
    do j = 2, njm1
      su(i, j) = te(i, j)*den(i, j)/abs(taun(i))
      sp(i, j) = te(i, j)**1.5/ed(i, j)/rlarge
    end do
  end do
  call tecprint(2, 2, ni, nj, it, jt, x, y, su, heda, 11)
  call tecprint(2, 2, ni, nj, it, jt, x, y, sp, hedb, 12)
!-----calculation of shear-stress coeff. and nusselt number 
! along large duct wall
  j = njm1
  xlref = rlarge - rsmall
  tref = tin
  delt = abs(tin-twalln)
  cpp = 4.190e3
  tk = 6.0e-1
  dyn = yv(nj) - y(njm1)
  write (3, 402)
  do i = istep, nim1
    ssc = taun(i)/(densit*ularge*ularge)
    xuh = xu(i)/(rlarge-rsmall)
    xh = x(i)/xlref
    uplusn = alog(elog*yplusn(i))/cappa
    hfluxn = den(i, j)*cpp*sqrt(te(i,j))*(tref-t(i,njm1))/prandt/(uplusn+pfun)
    if (yplusn(i)<=11.63) hfluxn = tk*abs(twalln-t(i,njm1))/dyn
    hfluxn = abs(hfluxn)
    anun = hfluxn*xlref/(tk*delt)
    write (3, 403) i, xuh, ssc, anun, xh
  end do
  stop

!-----format statements
  010 format (6a6)
  210 format(1h0,47x,'kase t1 turbulent flow through a sudden &
             enlargement'////)
  211 format(1h0,50x,'lesson 8 - - solution of energy equation'////)
  220 format(//1h0,15x,'inlet fluid velocity ',t60,1h=,3x,1pe11.3)
  221 format(1h0,15x,'inlet fluid temperature ',t60,1h=,3x,1pe11.3)
  222 format(///1h0,15x,'thermal boundary conditions are - - -'//&
      		1h0,25x,'prescribed pipe wall temperature = ',1pe11.3//&
      		1h0,25x,'adiabatic step wall '//) 
  223 format(1h0,15x,'prandtl number',t60,1h=,3x,1pe11.3)
  230 format(1h0,15x,'reynolds number ',t60,1h=,3x,1pe11.3) 
  240 format(1h0,15x,'diameter ratio ',t60,1h=,3x,1pe11.3)
  250 format(1h0,15x,'laminar viscosity ',t60,1h=,3x,1pe11.3)
  260 format(1h0,15x,'fluid density ',t60,1h=,3x,1pe11.3)
310 format(1h0,'iter   ','i---------------absolute residual source sum&
      s---------------i   i-------field values at monitoring location','&
      (',i2,',',i2,')','--------i' / 2x,'no.',3x,'umom',6x,'vmom',6x,'ma&
      ss',6x,'ener',6x,'tkin',6x,'disp',10x,'u',9x,'v',9x,'p',9x,'t',9x,&
      'k',9x,'d'/)
311   format(1h ,i3,4x,1p6e10.3,3x,1p6e10.3)
  402 format(///5x,1hi,7x,5hxu(i),6x,10hs.s.coeff.,'nusselt no. ',&
      5x,'x(i)')
  403 format(/5x,i5,3(1pe11.3),2x,1pe11.3)
end subroutine contro
