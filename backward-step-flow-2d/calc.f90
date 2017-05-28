!**********************************************************************	
!  Calculate variable profiles 
!  Ruipengyu Li April-2016				
!**********************************************************************	
module calc
  use general
  use promod
  real :: arean, areas, areaew
  real :: gn, gnw, gs, gsw, ge, gp, gw, gse
  real :: cn, cs, ce, cw
  real :: dn, ds, de, dw
  real :: visn, viss, vise, visw, vdr
  real :: smp, cp, cpo, denn, dens, dene, denw
  real :: gamp, gamm, gamn, gams, game, gamw
  real :: dudyp, dudym, dvdym, dvdyp, dudxp, dudxm, &
          dvdxm, dvdxp, dudx, dvdy
  real :: sorvol, resor, rgamp, rgamm, ppref
contains
!**********************************************************************
  subroutine calcu()
    implicit none
    integer :: i, j, n
    real :: vol

!chapter  0  0  0  0  0  0  0  0  preliminaries  0  0  0  0  0  0  0  0

!chapter  1  1  1  1  1  1  assembly of coefficients  1  1  1  1  1  1  1
    do i = 3, nim1
      do j = 2, njm1
!-----compute areas and volume
!  u-cell
        arean = rv(j+1)*sewu(i)
        areas = rv(j)*sewu(i)
        areaew = r(j)*sns(j)
        vol = r(j)*sewu(i)*sns(j)
!-----calculate convection coefficients 
!  interpolation of density  
!  mass flux of scalar cell faces
        gn = 0.5*(den(i,j+1)+den(i,j))*v(i, j+1)
        gnw = 0.5*(den(i-1,j)+den(i-1,j+1))*v(i-1, j+1)
        gs = 0.5*(den(i,j-1)+den(i,j))*v(i, j)
        gsw = 0.5*(den(i-1,j)+den(i-1,j-1))*v(i-1, j)
        ge = 0.5*(den(i+1,j)+den(i,j))*u(i+1, j)
!  mass flux at location of velocity
        gp = 0.5*(den(i,j)+den(i-1,j))*u(i, j)
        gw = 0.5*(den(i-1,j)+den(i-2,j))*u(i-1, j)
!  convection coefficients at cell face centre of u-cells
        cn = 0.5*(gn+gnw)*arean
        cs = 0.5*(gs+gsw)*areas
        ce = 0.5*(ge+gp)*areaew
        cw = 0.5*(gp+gw)*areaew
!-----calculate diffusion coefficients
!  evaluate at u-cell 
        visn = 0.25*(vis(i,j)+vis(i,j+1)+vis(i-1,j)+vis(i-1,j+1))
        viss = 0.25*(vis(i,j)+vis(i,j-1)+vis(i-1,j)+vis(i-1,j-1))
        dn = visn*arean/dynp(j)
        ds = viss*areas/dyps(j)
        de = vis(i, j)*areaew/dxepu(i)
        dw = vis(i-1, j)*areaew/dxpwu(i)
!-----calculate coefficients of source terms
!  net outflow from control volume
!  local mass sink
        smp = cn - cs + ce - cw
!  maximum of zero and net outflow
        cp = max(0.0, smp)
        cpo = cp
!-----assemble main coefficients
!  hybrid differencing scheme
        an(i, j) = max(abs(0.5*cn), dn) - 0.5*cn
        as(i, j) = max(abs(0.5*cs), ds) + 0.5*cs
        ae(i, j) = max(abs(0.5*ce), de) - 0.5*ce
        aw(i, j) = max(abs(0.5*cw), dw) + 0.5*cw
        du(i, j) = areaew
!  coefficient c of linearized source treatment
!  pressure as source term
        su(i, j) = cpo*u(i, j) + du(i, j)*(p(i-1,j)-p(i,j))
!  coefficient b of linearized source treatment
        sp(i, j) = -cp
!  viscous diffusion in source term
!  u-vel gradient at node (i,j)
        dudxp = (u(i+1,j)-u(i,j))/sew(i)
!  u-vel gradient at node (i-1,j)
        dudxm = (u(i,j)-u(i-1,j))/sew(i-1)
        su(i, j) = (vis(i,j)*dudxp-vis(i-1,j)*dudxm)/sewu(i)*vol + su(i, j)
!  viscosity at midpoint of upstream wall of cell
        gamp = 0.25*(vis(i,j)+vis(i-1,j)+vis(i,j+1)+vis(i-1,j+1))
!  v-vel gradient at midpoint of north wall of u-cell
        dvdxp = rv(j+1)*(v(i,j+1)-v(i-1,j+1))/dxep(i)
!  viscosity at midpoint of downstream wall of cell
        gamm = 0.25*(vis(i,j)+vis(i-1,j)+vis(i,j-1)+vis(i-1,j-1))
!  v-vel gradient at midpoint of south wall of u-cell
        dvdxm = rv(j)*(v(i,j)-v(i-1,j))/dxep(i)
!  finial su value including viscous terms
        su(i, j) = su(i, j) + (gamp*dvdxp-gamm*dvdxm)/sns(j)/r(j)*vol
      end do
    end do
!chapter  2  2  2  2  2  2  2  problem modifications  2  2  2  2  2  2  2
!  boundary condition
    call modu
!
!chapter  3  final coeff. assembly and residual source calculation  3  3
!
    resoru = 0.0
    do i = 3, nim1
      do j = 2, njm1
        ap(i, j) = an(i, j) + as(i, j) + ae(i, j) + aw(i, j) - sp(i, j)
!  coefficient of velocity correction term for u-vel
        du(i, j) = du(i, j)/ap(i, j)
!  residual source for individual control volume
        resor = an(i, j)*u(i, j+1) + as(i, j)*u(i, j-1) + ae(i, j)*u(i+1, j) &
              + aw(i, j)*u(i-1, j) - ap(i, j)*u(i, j) + su(i, j)
!  volume
        vol = r(j)*sew(i)*sns(j)
        sorvol = great*vol
        if (-sp(i,j)>0.5*sorvol) resor = resor/sorvol
!  sum of residual sources within calculation domain for u-equation
        resoru = resoru + abs(resor)
!-----under-relaxation
!  lbl formula
        ap(i, j) = ap(i, j)/urfu
        su(i, j) = su(i, j) + (1.-urfu)*ap(i, j)*u(i, j)
        du(i, j) = du(i, j)*urfu
      end do
    end do
!
!chapter  4  4  4  solution of difference equation  4  4  4  4  4  4  4
!
    do n = 1, nswpu
!  under-relaxation is performed 
!  implicitly through coefficient modification
      call lisolv(3, 2, ni, nj, it, jt, u)
    end do
  end subroutine calcu
!**********************************************************************

!**********************************************************************      
  subroutine calcv()
    implicit none
    integer :: i, j, n
    real :: vol

!chapter  0  0  0  0  0  0  0  0  preliminaries  0  0  0  0  0  0  0  0

!chapter  1  1  1  1  1  1  assembly of coefficients  1  1  1  1  1  1  1
    do i = 2, nim1
      do j = 3, njm1
!-----compute areas and volume
        arean = rcv(j+1)*sew(i)
        areas = rcv(j)*sew(i)
        areaew = rv(j)*snsv(j)
        vol = rv(j)*sew(i)*snsv(j)
!-----calculate convection coefficients 
        gn = 0.5*(den(i,j+1)+den(i,j))*v(i, j+1)
        gp = 0.5*(den(i,j)+den(i,j-1))*v(i, j)
        gs = 0.5*(den(i,j-1)+den(i,j-2))*v(i, j-1)
        ge = 0.5*(den(i+1,j)+den(i,j))*u(i+1, j)
        gse = 0.5*(den(i,j-1)+den(i+1,j-1))*u(i+1, j-1)
        gw = 0.5*(den(i,j)+den(i-1,j))*u(i, j)
        gsw = 0.5*(den(i,j-1)+den(i-1,j-1))*u(i, j-1)
        cn = 0.5*(gn+gp)*arean
        cs = 0.5*(gp+gs)*areas
        ce = 0.5*(ge+gse)*areaew
        cw = 0.5*(gw+gsw)*areaew
!-----calculate diffusion coefficients
        vise = 0.25*(vis(i,j)+vis(i+1,j)+vis(i,j-1)+vis(i+1,j-1))
        visw = 0.25*(vis(i,j)+vis(i-1,j)+vis(i,j-1)+vis(i-1,j-1))
        dn = vis(i, j)*arean/dynpv(j)
        ds = vis(i, j-1)*areas/dypsv(j)
        de = vise*areaew/dxep(i)
        dw = visw*areaew/dxpw(i)
!-----calculate coefficients of source terms
        smp = cn - cs + ce - cw
        cp = max(0.0, smp)
        cpo = cp
!-----assemble main coefficients
        an(i, j) = max(abs(0.5*cn), dn) - 0.5*cn
        as(i, j) = max(abs(0.5*cs), ds) + 0.5*cs
        ae(i, j) = max(abs(0.5*ce), de) - 0.5*ce
        aw(i, j) = max(abs(0.5*cw), dw) + 0.5*cw
        dv(i, j) = 0.5*(arean+areas)
        su(i, j) = cpo*v(i, j) + dv(i, j)*(p(i,j-1)-p(i,j))
        sp(i, j) = -cp
        if (indcos==2) sp(i, j) = sp(i, j) - vis(i, j)*vol/rv(j)**2
        if (indcos==2) sp(i, j) = sp(i, j) - vis(i, j)*vol/rv(j)**2
        dudyp = (u(i+1,j)-u(i+1,j-1))/dyps(j)
        gamp = 0.25*(vis(i,j)+vis(i+1,j)+vis(i,j-1)+vis(i+1,j-1))
        gamm = 0.25*(vis(i,j)+vis(i-1,j)+vis(i,j-1)+vis(i-1,j-1))
        dudym = (u(i,j)-u(i,j-1))/dyps(j)
        su(i, j) = su(i, j) + (gamp*dudyp-gamm*dudym)/sew(i)*vol
        dvdyp = (v(i,j+1)-v(i,j))/sns(j)
        rgamp = vis(i, j)*r(j)
        dvdym = (v(i,j)-v(i,j-1))/sns(j-1)
        rgamm = vis(i, j-1)*r(j-1)
        su(i, j) = su(i, j) + (rgamp*dvdyp-rgamm*dvdym)/(r(j)*sns(j))*vol
      end do
    end do
!
!chapter  2  2  2  2  2  2  2  problem modifications  2  2  2  2  2  2  2
!
    call modv
!
!chapter  3  final coeff. assembly and residual source calculation  3  3
!
    resorv = 0.0
    do i = 2, nim1
      do j = 3, njm1
        ap(i, j) = an(i, j) + as(i, j) + ae(i, j) + aw(i, j) - sp(i, j)
        dv(i, j) = dv(i, j)/ap(i, j)
        resor = an(i, j)*v(i, j+1) + as(i, j)*v(i, j-1) + ae(i, j)*v(i+1, j) &
              + aw(i, j)*v(i-1, j) - ap(i, j)*v(i, j) + su(i, j)
        vol = r(j)*sew(i)*sns(j)
        sorvol = great*vol
        if (-sp(i,j)>0.5*sorvol) resor = resor/sorvol
        resorv = resorv + abs(resor)
!-----under-relaxation
        ap(i, j) = ap(i, j)/urfv
        su(i, j) = su(i, j) + (1.-urfv)*ap(i, j)*v(i, j)
        dv(i, j) = dv(i, j)*urfv
      end do
    end do
!
!chapter  4  4  4  solution of difference equation  4  4  4  4  4  4  4
!
    do n = 1, nswpv
      call lisolv(2, 3, ni, nj, it, jt, v)
    end do
  end subroutine calcv
!**********************************************************************

!**********************************************************************      
  subroutine calcp()
    implicit none
    integer :: i, j, n
    real :: vol

!chapter  0  0  0  0  0  0  0  0  preliminaries  0  0  0  0  0  0  0  0
!
    resorm = 0.0
!chapter  1  1  1  1  1  1  assembly of coefficients  1  1  1  1  1  1  1

    do i = 2, nim1
      do j = 2, njm1
!  scalar cell
!-----compute areas and volume
        arean = rv(j+1)*sew(i)
        areas = rv(j)*sew(i)
        areaew = r(j)*sns(j)
!-----calculate coefficients
        denn = 0.5*(den(i,j)+den(i,j+1))
        dens = 0.5*(den(i,j)+den(i,j-1))
        dene = 0.5*(den(i,j)+den(i+1,j))
        denw = 0.5*(den(i,j)+den(i-1,j))
        an(i, j) = denn*arean*dv(i, j+1)
        as(i, j) = dens*areas*dv(i, j)
        ae(i, j) = dene*areaew*du(i+1, j)
        aw(i, j) = denw*areaew*du(i, j)
!-----calculate source terms
        cn = denn*v(i, j+1)*arean
        cs = dens*v(i, j)*areas
        ce = dene*u(i+1, j)*areaew
        cw = denw*u(i, j)*areaew
        smp = cn - cs + ce - cw
        sp(i, j) = 0.0
        su(i, j) = -smp
!-----compute sum of absolute mass sources
        resorm = resorm + abs(smp)
      end do
    end do
!
!chapter  2  2  2  2  2  2  2  problem modifications  2  2  2  2  2  2  2
!
    call modp
!
!chapter  3  3  3  3  3  final coefficient assembly  3  3  3  3  3  3  3
!
    do i = 2, nim1
      do j = 2, njm1
        ap(i, j) = an(i, j) + as(i, j) + ae(i, j) + aw(i, j) - sp(i, j)
      end do
    end do
!
!chapter  4  4  4  4  4  solution of difference equations  4  4  4  4  4
!
!chapter  5  5  5  5  correct velocities and pressure  5  5  5  5  5  5
    do n = 1, nswpp
      call lisolv(2, 2, ni, nj, it, jt, pp)
    end do
!
!-----velocities
    do i = 2, nim1
      do j = 2, njm1
!  velocity correction formula
        if (i/=2) u(i, j) = u(i, j) + du(i, j)*(pp(i-1,j)-pp(i,j))
        if (j/=2) v(i, j) = v(i, j) + dv(i, j)*(pp(i,j-1)-pp(i,j))
      end do
    end do
!-----pressures (with provision for under-relaxation)
    ppref = pp(ipref, jpref)
    do i = 2, nim1
      do j = 2, njm1
        p(i, j) = p(i, j) + urfp*(pp(i,j)-ppref)
        pp(i, j) = 0.0
      end do
    end do
  end subroutine calcp
!**********************************************************************

!**********************************************************************
  subroutine calct()
    implicit none
    integer :: i, j, n
    real :: vol
!
!chapter  0  0  0  0  0  0  0  preliminaries  0  0  0  0  0  0  0
!
!
!chapter  1  1  1  1  1  1  assembly of coefficients  1  1  1  1  1  1 
!
    do i = 2, nim1
      do j = 2, njm1
!-----compute areas and volume
!  scalar cell
        arean = rv(j+1)*sew(i)
        areas = rv(j)*sew(i)
        areaew = r(j)*sns(j)
        vol = r(j)*sns(j)*sew(i)
!-----calculate convection coefficients 
        gn = 0.5*(den(i,j)+den(i,j+1))*v(i, j+1)
        gs = 0.5*(den(i,j)+den(i,j-1))*v(i, j)
        ge = 0.5*(den(i,j)+den(i+1,j))*u(i+1, j)
        gw = 0.5*(den(i,j)+den(i-1,j))*u(i, j)
        cn = gn*arean
        cs = gs*areas
        ce = ge*areaew
        cw = gw*areaew
!-----calculate diffusion coefficients
        gamn = 0.5*(gamh(i,j)+gamh(i,j+1))
        gams = 0.5*(gamh(i,j)+gamh(i,j-1))
        game = 0.5*(gamh(i,j)+gamh(i+1,j))
        gamw = 0.5*(gamh(i,j)+gamh(i-1,j))
        dn = gamn*arean/dynp(j)
        ds = gams*areas/dyps(j)
        de = game*areaew/dxep(i)
        dw = gamw*areaew/dxpw(i)
!-----source terms
        smp = cn - cs + ce - cw
        cp = max(0.0, smp)
        cpo = cp
!-----assemble main coefficients
        an(i, j) = max(abs(0.5*cn), dn) - 0.5*cn
        as(i, j) = max(abs(0.5*cs), ds) + 0.5*cs
        ae(i, j) = max(abs(0.5*ce), de) - 0.5*ce
        aw(i, j) = max(abs(0.5*cw), dw) + 0.5*cw
        su(i, j) = cpo*t(i, j)
        sp(i, j) = -cp
      end do
    end do
!
!chapter  2  2  2  2  2  2  problem modifications  2  2  2  2  2  2
!
    call modt
!chapter 3 final coefficient assembly and residual source calculation 3
!
    resort = 0.0
    do i = 2, nim1
      do j = 2, njm1
        ap(i, j) = an(i, j) + as(i, j) + ae(i, j) + aw(i, j) - sp(i, j)
        resor = an(i, j)*t(i, j+1) + as(i, j)*t(i, j-1) + ae(i, j)*t(i+1, j) &
              + aw(i, j)*t(i-1, j) - ap(i, j)*t(i, j) + su(i, j)
        vol = r(j)*sew(i)*sns(j)
        sorvol = great*vol
        if (-sp(i,j)>0.5*sorvol) resor = resor/sorvol
        resort = resort + abs(resor)
!-----under-relaxation
        ap(i, j) = ap(i, j)/urft
        su(i, j) = su(i, j) + (1.0-urft)*ap(i, j)*t(i, j)
      end do
    end do
!
!chapter  4  4  4  4  4  solution of difference equations  4  4  4  4  4
!
    do n = 1, nswpt
      call lisolv(2, 2, ni, nj, it, jt, t)
    end do
  end subroutine calct

!**********************************************************************      
  subroutine calcte()
    implicit none
    integer :: i, j, n
    real :: vol
!
!chapter  0  0  0  0  0  0  0  preliminaries  0  0  0  0  0  0  0
!

!
!chapter  1  1  1  1  1  1  assembly of coefficients  1  1  1  1  1  1 
!
!  sigma_k
    prte = 1.0
    do i = 2, nim1
      do j = 2, njm1
!-----compute areas and volume
        arean = rv(j+1)*sew(i)
        areas = rv(j)*sew(i)
        areaew = r(j)*sns(j)
        vol = r(j)*sns(j)*sew(i)
!-----calculate convection coefficients 
        gn = 0.5*(den(i,j)+den(i,j+1))*v(i, j+1)
        gs = 0.5*(den(i,j)+den(i,j-1))*v(i, j)
        ge = 0.5*(den(i,j)+den(i+1,j))*u(i+1, j)
        gw = 0.5*(den(i,j)+den(i-1,j))*u(i, j)
        cn = gn*arean
        cs = gs*areas
        ce = ge*areaew
        cw = gw*areaew
!-----calculate diffusion coefficients
        gamn = 0.5*(vis(i,j)+vis(i,j+1))/prte
        gams = 0.5*(vis(i,j)+vis(i,j-1))/prte
        game = 0.5*(vis(i,j)+vis(i+1,j))/prte
        gamw = 0.5*(vis(i,j)+vis(i-1,j))/prte
        dn = gamn*arean/dynp(j)
        ds = gams*areas/dyps(j)
        de = game*areaew/dxep(i)
        dw = gamw*areaew/dxpw(i)
!-----source terms
        smp = cn - cs + ce - cw
        cp = max(0.0, smp)
        cpo = cp
        dudx = (u(i+1,j)-u(i,j))/sew(i)
        dvdy = (v(i,j+1)-v(i,j))/sns(j)
        dudy = ((u(i,j)+u(i+1,j)+u(i,j+1)+u(i+1,j+1))/4. &
             -(u(i,j)+u(i+1,j)+u(i,j-1)+u(i+1,j-1))/4.)/sns(j)
        dvdx = ((v(i,j)+v(i,j+1)+v(i+1,j)+v(i+1,j+1))/4. &
             -(v(i,j)+v(i,j+1)+v(i-1,j)+v(i-1,j+1))/4.)/sew(i)
!  generation of te
        gen(i, j) = (2.*(dudx**2+dvdy**2)+(dudy+dvdx)**2)*vis(i, j)
        vdr = v(i, j)/rv(j)
        if (j==2) vdr = 0.0
        if (indcos==2) gen(i, j) = gen(i, j) + &
                                   vis(i, j)*0.5*(vdr+v(i,j+1)/rv(j+1))**2
!-----assemble main coefficients
        an(i, j) = max(abs(0.5*cn), dn) - 0.5*cn
        as(i, j) = max(abs(0.5*cs), ds) + 0.5*cs
        ae(i, j) = max(abs(0.5*ce), de) - 0.5*ce
        aw(i, j) = max(abs(0.5*cw), dw) + 0.5*cw
        su(i, j) = cpo*te(i, j)
        sukd(i, j) = su(i, j)
!  c in source term
        su(i, j) = su(i, j) + gen(i, j)*vol
        sp(i, j) = -cp
        spkd(i, j) = sp(i, j)
!  b in source term
        sp(i, j) = sp(i, j) - cd*cmu*den(i, j)**2*te(i, j)*vol/vis(i, j)
      end do
    end do
!
!chapter  2  2  2  2  2  2  problem modifications  2  2  2  2  2  2
!
    call modte
!
!chapter 3 final coefficient assembly and residual source calculation 3
!
    resork = 0.0
    do i = 2, nim1
      do j = 2, njm1
        ap(i, j) = an(i, j) + as(i, j) + ae(i, j) + aw(i, j) - sp(i, j)
        resor = an(i, j)*te(i, j+1) + as(i, j)*te(i, j-1) + ae(i, j)*te(i+1, j) &
              + aw(i, j)*te(i-1, j) - ap(i, j)*te(i, j) + su(i, j)
        vol = r(j)*sew(i)*sns(j)
        sorvol = great*vol
        if (-sp(i,j)>0.5*sorvol) resor = resor/sorvol
        resork = resork + abs(resor)
!-----under-relaxation
        ap(i, j) = ap(i, j)/urfk
        su(i, j) = su(i, j) + (1.-urfk)*ap(i, j)*te(i, j)
      end do
    end do
!
!chapter  4  4  4  4  4  solution of difference equations  4  4  4  4  4
!
    do n = 1, nswpk
      call lisolv(2, 2, ni, nj, it, jt, te)
    end do
  end subroutine calcte
!**********************************************************************


!**********************************************************************
  subroutine calced()
    implicit none
    integer :: i, j, n
    real :: vol
!
!chapter  0  0  0  0  0  0  0  preliminaries  0  0  0  0  0  0  0
!

!chapter  1  1  1  1  1  1  assembly of coefficients  1  1  1  1  1  1 
!
    do i = 2, nim1
      do j = 2, njm1
!-----compute areas and volume
        arean = rv(j+1)*sew(i)
        areas = rv(j)*sew(i)
        areaew = r(j)*sns(j)
        vol = r(j)*sns(j)*sew(i)
!-----calculate convection coefficients 
        gn = 0.5*(den(i,j)+den(i,j+1))*v(i, j+1)
        gs = 0.5*(den(i,j)+den(i,j-1))*v(i, j)
        ge = 0.5*(den(i,j)+den(i+1,j))*u(i+1, j)
        gw = 0.5*(den(i,j)+den(i-1,j))*u(i, j)
        cn = gn*arean
        cs = gs*areas
        ce = ge*areaew
        cw = gw*areaew
!-----calculate diffusion coefficients
        gamn = 0.5*(vis(i,j)+vis(i,j+1))/pred
        gams = 0.5*(vis(i,j)+vis(i,j-1))/pred
        game = 0.5*(vis(i,j)+vis(i+1,j))/pred
        gamw = 0.5*(vis(i,j)+vis(i-1,j))/pred
        dn = gamn*arean/dynp(j)
        ds = gams*areas/dyps(j)
        de = game*areaew/dxep(i)
        dw = gamw*areaew/dxpw(i)
!-----source terms
        smp = cn - cs + ce - cw
        cp = max(0.0, smp)
        cpo = cp
!-----assemble main coefficients
        an(i, j) = max(abs(0.5*cn), dn) - 0.5*cn
        as(i, j) = max(abs(0.5*cs), ds) + 0.5*cs
        ae(i, j) = max(abs(0.5*ce), de) - 0.5*ce
        aw(i, j) = max(abs(0.5*cw), dw) + 0.5*cw
        su(i, j) = cpo*ed(i, j)
        sukd(i, j) = su(i, j)
!  c in source term (expressed in te)
!  25.04.2016 Ruipengyu Li
!  effective viscosity is used instead of the turbulent viscosity
!  the original source term in the epsilon equation is 
!  su(i, j) = su(i, j) + c1*ed(i, j)*gen(i, j)*vol/te(i, j)
!  
        su(i,j)=su(i,j)+c1*cmu*gen(i,j)*vol*den(i,j)*te(i,j)/vis(i,j)
        sp(i, j) = -cp
        spkd(i, j) = sp(i, j)
!  b in source term
        sp(i, j) = sp(i, j) - c2*den(i, j)*ed(i, j)*vol/te(i, j)
      end do
    end do
!
!chapter  2  2  2  2  2  2  problem modifications  2  2  2  2  2  2
!
    call moded
!
!chapter 3 final coefficient assembly and residual source calculation 3
    resore = 0.0
!
    do i = 2, nim1
      do j = 2, njm1
        ap(i, j) = an(i, j) + as(i, j) + ae(i, j) + aw(i, j) - sp(i, j)
        resor = an(i, j)*ed(i, j+1) + as(i, j)*ed(i, j-1) + ae(i, j)*ed(i+1, j)&
              + aw(i, j)*ed(i-1, j) - ap(i, j)*ed(i, j) + su(i, j)
        vol = r(j)*sns(j)*sew(i)
        sorvol = great*vol
        if (-sp(i,j)>0.5*sorvol) resor = resor/sorvol
        resore = resore + abs(resor)
!-----under-relaxation
        ap(i, j) = ap(i, j)/urfe
        su(i, j) = su(i, j) + (1.-urfe)*ap(i, j)*ed(i, j)
      end do
    end do

!chapter  4  4  4  4  4  solution of difference equations  4  4  4  4  4

    do n = 1, nswpd
      call lisolv(2, 2, ni, nj, it, jt, ed)
    end do
  end subroutine calced
!**********************************************************************     
end module calc

