!**********************************************************************	
!  Problem modification 
!  Boundary implementation
!  May add additional walls, extend flow domain,etc.
!**********************************************************************	
module promod
  use general
  real :: cdterm, sqrtk, tmult, arden, ardent, flow, uinc, area
  real :: gt, gencou, genres, diterm
  real :: dudy, dvdx, denu, denv, dyn
  real :: xp, yp
  real :: xplusa, yplusa
  real :: uplus, vplus
contains
!chapter  0  0  0  0  0  0  0  preliminaries  0  0  0  0  0  0  0  0  0

!chapter  1  1  1  1  1  1  1  1  properties  1  1  1  1  1  1  1  1  1
!**********************************************************************
  subroutine modpro()
    implicit none
!-----no modifications for this problem 
  end subroutine modpro
!**********************************************************************
!chapter  2  2  2  2  2  2  2  2  u momentum  2  2  2  2  2  2  2  2  2
!**********************************************************************
  subroutine modu()
    implicit none
    integer :: i, j
    real :: area, vol, term
!-----top wall
!  inserts shear force along cylindrical wall 
    cdterm = cmu**0.25
!  distance from north wall (cell face) 
    yp = yv(nj) - y(njm1)
    j = njm1
    do i = 3, nim1
      sqrtk = sqrt(0.5*(te(i,j)+te(i-1,j)))
      denu = 0.5*(den(i,j)+den(i-1,j))
!  n.b. yplusn is set to 11 initially
!  yplus is calculated in the k-equation part
!  considering turbulence generation == dissipation 
      yplusa = 0.5*(yplusn(i)+yplusn(i-1))
!  select the appropriate stress equation according to y+  
      if (yplusa<=11.63) then
        tmult = viscos/yp
!  coefficient of wall shear stress expression
!  n.b. elog == 9.793
      else
        tmult = denu*cdterm*sqrtk*cappa/log(elog*yplusa)
      end if
!  north wall shear stress  
      taun(i) = -tmult*u(i, j)
!  north wall shear force   
!  insert via source term sp(i,j)    
!  sp term is coefficient of u
!  only sp term since u-velocity at wall is zero
      sp(i, j) = sp(i, j) - tmult*sewu(i)*rv(nj)
!  the usual shear expression is suppressed by setting an(i,j)=0    
      an(i, j) = 0.0
    end do
!  values for the edges
    taun(2) = taun(3)
    taun(ni) = taun(nim1)
!-----side wall
!  no treatment for velocity normal to the wall
!  set normal fluxes to zero
    do j = jstp1, njm1
      aw(3, j) = 0.0
    end do
!-----symmetry  axis
!  set normal fluxes to zero
    do i = 1, ni
!  symmetric boundary lies on scalar cell face 
!  between scalar points (i,1) and (i,2)
      u(i, 1) = u(i, 2)
      as(i, 2) = 0.0
    end do
!-----outlet
    ardent = 0.0
    flow = 0.0
    do j = 2, njm1
!  area of east/west wall times density of fluid
      arden = 0.5*(den(nim1,j)+den(nim1-1,j))*r(j)*sns(j)
!  sum at cross-section
      ardent = ardent + arden
!  mass flow rate   
      flow = flow + arden*u(nim1, j)
    end do
!  uniform increment of u-velocity at outlet of flow domain  
    uinc = (flowin-flow)/ardent
!  obtain u at outlet by overall mass balance
    do j = 2, njm1
      u(ni, j) = u(nim1, j) + uinc
    end do
  end subroutine modu
!**********************************************************************
!chapter  3  3  3  3  3  3  3  3  v momentum  3  3  3  3  3  3  3  3  3
!**********************************************************************
  subroutine modv()
    implicit none
    integer :: i, j
    real :: area, vol, term
!-----side wall
    cdterm = cmu**0.25
    xp = x(istep) - xu(istep)
    i = istep
    do j = jstp1, njm1
      sqrtk = sqrt(0.5*(te(i,j)+te(i,j-1)))
      denv = 0.5*(den(i,j)+den(i,j-1))
      xplusa = 0.5*(xplusw(j)+xplusw(j-1))
      if (xplusa<=11.63) then
        tmult = viscos/xp
      else
        tmult = denv*cdterm*sqrtk*cappa/log(elog*xplusa)
      end if
      tauw(j) = -tmult*v(i, j)
      sp(i, j) = sp(i, j) - tmult*snsv(j)*rv(j)
      aw(i, j) = 0.0
    end do
    tauw(jstep) = tauw(jstp1)
    tauw(nj) = tauw(njm1)
!-----top wall
    do i = istep, nim1
      an(i, njm1) = 0.0
    end do
!-----symmetry  axis
    do i = 2, nim1
      as(i, 3) = 0.0
    end do
!-----outlet
!  no condition needed if the pipe is sufficiently long
!  from corners or steps
  end subroutine modv
!**********************************************************************
!chapter  4  4  4  4  4  4  pressure correction  4  4  4  4  4  4  4  4
!**********************************************************************
  subroutine modp()
    implicit none
!  no modification needed
!  pressure is guessed and evaluated
!  only pressure difference matters  
  end subroutine modp
!**********************************************************************
!chapter  5  5  5  5  5  5  5  thermal energy  5  5  5  5  5  5  5  5  5
!**********************************************************************
  subroutine modt()
    implicit none
    integer :: i, j
    real :: area, vol, term
!-----top wall (constant temperature)
    cdterm = cmu**0.25
    j = njm1
!  distance between first grid and wall
    dyn = yv(nj) - y(njm1)
    do i = 2, nim1
      an(i, j) = 0.0
      area = rv(j+1)*sew(i)
      if (yplusn(i)<=11.63) then
        gt = viscos/(prandl*dyn)
      else
        uplus = log(elog*yplusn(i))/cappa
        gt = den(i, j)*cdterm*sqrt(te(i,j))/(prandt*(uplus+pfun))
      end if
      term = gt*area
!  heat rate is incorporated through two source terms
!  why added to two terms?
!  separate the equation in terms of t and twall(!=0) 
      su(i, j) = su(i, j) + term*twalln
      sp(i, j) = sp(i, j) - term
    end do
!-----side wall (adiabatic)
    i = 2
    do j = jstp1, njm1
      aw(i, j) = 0.0
    end do
!-----symmetry axis 
    j = 2
    do i = 2, nim1
      as(i, j) = 0.0
    end do
  end subroutine modt
!**********************************************************************
!chapter 6 6 6 6 6 6 6 turbulent kinetic energy
!**********************************************************************
  subroutine modte()
    implicit none
    integer :: i, j
    real :: area, vol, term
!-----top wall
    cdterm = cmu**0.25
    yp = yv(nj) - y(njm1)
    j = njm1
    do i = 2, nim1
      denu = den(i, j)
      sqrtk = sqrt(te(i,j))
      vol = r(j)*sns(j)*sew(i)
!  modify generation term  
!  part of generation term modified in terms of wall shear stress
      gencou = 0.5*(abs(taun(i+1)*u(i+1,j))+abs(taun(i)*u(i,j)))/yp
!  y+ formula consider k production == dissipation
      yplusn(i) = denu*sqrtk*cdterm*yp/viscos
!  u-velocity evaluated at cell centre
!  average of nearby values  
      dudy = ((u(i,j)+u(i+1,j)+u(i,j+1)+u(i+1,j+1))/4.&
           -(u(i,j)+u(i+1,j)+u(i,j-1)+u(i+1,j-1))/4.)/sns(j)
!  replace term in original generation     
      genres = gen(i, j) - vis(i, j)*dudy**2
      gen(i, j) = genres + gencou
!  second term in k-source    
      if (yplusn(i)<=11.63) then
!  u+ == y+    
        diterm = den(i, j)*(cmu**.75)*sqrtk*yplusn(i)/yp
      else
!  u+ == log(ey+)/cappa   
        diterm = den(i, j)*(cmu**.75)*sqrtk*log(elog*yplusn(i))/(cappa*yp)
      end if
!  incorporated through source su and sp    
!  volumetric generation
      su(i, j) = gen(i, j)*vol + sukd(i, j)
      sp(i, j) = -diterm*vol + spkd(i, j)
!  no flux contributes to wall
!  turbulence energy goes to zero at a wall     
      an(i, j) = 0.0
    end do
!-----side wall
    xp = x(istep) - xu(istep)
    i = istep
    do j = jstp1, njm1
      denv = den(i, j)
      sqrtk = sqrt(te(i,j))
      vol = r(j)*sns(j)*sew(i)
      xplusw(j) = denv*sqrtk*cdterm*xp/viscos
      gencou = 0.5*(abs(tauw(j+1)*v(i,j+1))+abs(tauw(j)*v(i,j)))/xp
      dvdx = ((v(i,j)+v(i,j+1)+v(i+1,j)+v(i+1,j+1))/4. &
           -(v(i,j)+v(i,j+1)+v(i-1,j)+v(i-1,j+1))/4.)/sew(i)
      genres = gen(i, j) - vis(i, j)*dvdx**2
      gen(i, j) = genres + gencou
      if (xplusw(j)<=11.63) then
        diterm = den(i, j)*(cmu**.75)*sqrtk*xplusw(j)/xp
      else
        diterm = den(i, j)*(cmu**.75)*sqrtk*log(elog*xplusw(j))/(cappa*xp)
      end if
      su(i, j) = sukd(i, j) + gen(i, j)*vol
      sp(i, j) = spkd(i, j) - diterm*vol
      aw(i, j) = 0.0
    end do
!-----symmetry  axis
    j = 2
    do i = 2, nim1
      te(i, 1) = te(i, 2)
      dudy = ((u(i,j)+u(i+1,j)+u(i,j+1)+u(i+1,j+1))/4. &
           -(u(i,j)+u(i+1,j)+u(i,j-1)+u(i+1,j-1))/4.)/sns(j)
      vol = r(j)*sns(j)*sew(i)
      gen(i, j) = gen(i, j) - vis(i, j)*dudy**2
      su(i, j) = sukd(i, j) + gen(i, j)*vol
      as(i, 2) = 0.0
    end do
  end subroutine modte
!**********************************************************************
!chapter  7 7 7 7 7 7 7 7 dissipation
!**********************************************************************
  subroutine moded()
!  energy dissipation is fixed by equilibrium relations
    implicit none
    integer :: i, j
    real :: term
!-----top wall
    yp = yv(nj) - y(njm1)
    j = njm1
!  unlike k, epsilon reaches its highest value (much higher than in a
!  free stream) at a wall  
!  it's difficult to modify as(i,j) in such cases
!  equilibrium relation
    term = (cmu**.75)/(cappa*yp)
    do i = 2, nim1
!  overwrite through source coefficients 
!  effective replacement
!  sp == -10^30
!  su == epsilon*10^30  
      su(i, j) = great*term*te(i, j)**1.5
      sp(i, j) = -great
    end do
!-----side wall
    xp = x(istep) - xu(istep)
    i = istep
    term = (cmu**.75)/(cappa*xp)
    njm2 = nj - 2
    do j = jstp1, njm2
      su(i, j) = great*term*te(i, j)**1.5
      sp(i, j) = -great
    end do
!-----symmetry  axis
    do i = 2, nim1
      as(i, 2) = 0.0
    end do
  end subroutine moded
!**********************************************************************
end module promod
!**********************************************************************

