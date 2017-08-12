!*****************************************************************************
! Lid-driven cavity flow
!
!   Finite difference method
!   Cell-centred grid
!   Implicit and steady state SIMPLE are implemented
!
!   Note that the steady state SIMPLE only converges under a very
!   coarse mesh (21*21) and low Re (10) in this particular case.
!   
!   Implicit SIMPLE can deal with higher Re (hundreds) but time step and 
!   underrelaxation factor have to be sufficiently small. And the 
!   computational time may be large.
!
!   If both transient and steady state simulations can converge, the
!   Implicit SIMPLE will take more iterations than the steady state one.
!
! Ruipengyu Li
! Modified: 11/08/2017
!*****************************************************************************
module cavity_mod

implicit none

integer, parameter :: dp = selected_real_kind(15)

contains
!*****************************************************************************
subroutine sipsol(i1, i2, j1, j2, aw, ae, as, an, ap, su, phi)
! SIP solver, ILU of Stone (1968)
implicit none
integer :: i, j, ierr, iter
integer, intent(in) :: i1, i2, j1, j2
integer, parameter :: maxit = 3000
real(dp), parameter :: alpha = 0.80_dp
real(dp), parameter :: tol = 1.0e-6_dp
real(dp) :: p1, p2, rsm, resl, res1
real(dp), dimension(:,:), intent(in) :: aw, ae, as, an, ap, su
real(dp), dimension(:,:), intent(out) :: phi
real(dp), allocatable, dimension(:,:) :: lw, ls, lpr, un, ue, res

allocate(lw(i1:i2,j1:j2), ls(i1:i2,j1:j2), lpr(i1:i2,j1:j2), &
         un(i1:i2,j1:j2), ue(i1:i2,j1:j2), res(i1:i2,j1:j2), &
         stat=ierr)
lw(:,:) = 0.0_dp; ls(:,:) = 0.0_dp; lpr(:,:) = 0.0_dp
un(:,:) = 0.0_dp; ue(:,:) = 0.0_dp; res(:,:) = 0.0_dp

!-----Calculate coefficients of [L] and [U] matrices
do j = j1+1, j2-1
  do i = i1+1, i2-1
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
do iter = 1, maxit
  resl = 0.0_dp
  do j = j1+1, j2-1
    do i = i1+1, i2-1
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
  do j = j2-1, j1+1, -1
    do i = i2-1, i1+1, -1
      res(i,j) = res(i,j) - un(i,j)*res(i,j+1) - ue(i,j)*res(i+1,j)
      phi(i,j) = phi(i,j) + res(i,j)
    end do
  end do
  ! check convergence
  if (rsm < tol) then 
    exit
  else if (iter == 1 .and. res1 < 1.0e-10_dp) then
    exit
  else if (iter == maxit) then
    !write(*,*) 'SIP solver - convergence not reached'
  end if
end do
deallocate(lw, ls, lpr, un, ue, res, stat=ierr)
end subroutine sipsol
!*****************************************************************************
end module cavity_mod
!*****************************************************************************
program main
use cavity_mod
implicit none

integer :: i, j, ni, nj
integer :: itsp, ntsp, iter
integer :: ierr
integer :: imon, jmon ! monitor
character(len=80) :: msg
logical :: steadysim = .false.
real(dp) :: dx, dy, dt 
real(dp) :: utop, re ! Reynolds no.
real(dp) :: urfu, urfv, urfp
real(dp) :: resnorm, error, tol ! for convergence
real(dp) :: vel1, vel2 ! for advection term uv
real(dp) :: rex, rey, pterm ! some terms in coefficients
real(dp), allocatable, dimension(:) :: x, y
real(dp), allocatable, dimension(:,:) :: uo, vo, po ! uncorrected
real(dp), allocatable, dimension(:,:) :: pp ! correction
real(dp), allocatable, dimension(:,:) :: u, v, p ! corrected
real(dp), allocatable, dimension(:,:) :: uc, vc, pc ! grid values
real(dp), allocatable, dimension(:,:) :: ae, aw, an, as, ap, su ! grid values

ni = 11
nj = 11
imon = ni/2
jmon = nj/2
dt = 1.0e-4_dp
ntsp = 999999
utop = 1.0_dp
re = 10.0_dp
error = 0.0_dp
tol = 1.0e-15_dp
urfu = 1.0_dp
urfv = 1.0_dp
urfp = 0.5_dp ! under-relaxation factor for pressure
steadysim = .false.

! Under current grid layout, u and v have one node 
! less than p in y and x direction, respectively
allocate(u(1:ni+1,1:nj+1), v(1:ni+1,1:nj+1), p(1:ni+1,1:nj+1), &
         uo(1:ni+1,1:nj+1), vo(1:ni+1,1:nj+1), po(1:ni+1,1:nj+1), &
         ae(1:ni+1,1:nj+1), aw(1:ni+1,1:nj+1), an(1:ni+1,1:nj+1), &
         as(1:ni+1,1:nj+1), ap(1:ni+1,1:nj+1), su(1:ni+1,1:nj+1), &
         pp(1:ni+1,1:nj+1), x(1:ni), y(1:nj), &
         uc(1:ni,1:nj), vc(1:ni,1:nj), pc(1:ni,1:nj), &
         stat=ierr, errmsg=msg)
u(:,:) = 0.0_dp      
v(:,:) = 0.0_dp    
p(:,:) = 0.0_dp 
uo(:,:) = 0.0_dp  
vo(:,:) = 0.0_dp      
po(:,:) = 0.0_dp     
pp(:,:) = 0.0_dp 
uc(:,:) = 0.0_dp     
vc(:,:) = 0.0_dp 
pc(:,:) = 0.0_dp      
ae(:,:) = 0.0_dp
aw(:,:) = 0.0_dp
an(:,:) = 0.0_dp
as(:,:) = 0.0_dp
ap(:,:) = 0.0_dp
su(:,:) = 0.0_dp
x(:) = 0.0_dp      
y(:) = 0.0_dp      

dx = 1.0_dp / (ni-1)
dy = 1.0_dp / (nj-1)

! Grid
do i = 1, ni
  x(i) = (i-1) * dx
end do
do j = 1, nj
  y(j) = (j-1) * dy
end do

! Initial condition
uo(1:ni,1:nj-1) = 0.0_dp
uo(1:ni,nj:nj+1) = 1.0_dp
vo(1:ni+1,1:nj) = 0.0_dp
po(1:ni+1,1:nj+1) = 0.0_dp

rex = 1.0_dp/re/dx**2
rey = 1.0_dp/re/dy**2
pterm = 1.0_dp/dt + 2.0_dp*(rex+rey)
do itsp = 1, ntsp
!-----u momentum
  do j = 2, nj
    do i = 2, ni-1
      vel1 = 0.5_dp*(vo(i,j)+vo(i+1,j))
      vel2 = 0.5_dp*(vo(i,j-1)+vo(i+1,j-1))
      ae(i,j) = -0.5_dp*uo(i+1,j)/dx + rex
      aw(i,j) = 0.5_dp*uo(i-1,j)/dx + rex
      an(i,j) = -0.5_dp*vel1/dy + rey
      as(i,j) = 0.5_dp*vel2/dy + rey
      if (steadysim) then
        su(i,j) = (po(i,j)-po(i+1,j))/dx
        ap(i,j) = 2.0_dp*(rex+rey)
      else
        su(i,j) = uo(i,j)/dt + (po(i,j)-po(i+1,j))/dx
        ap(i,j) = pterm
      end if
      ! under-relax
      su(i,j) = su(i,j) + (1-urfu)*ap(i,j)/urfu*uo(i,j)
      ap(i,j) = ap(i,j)/urfu
    end do
  end do
  call sipsol(1, ni, 1, nj+1, -aw, -ae, -as, -an, ap, su, u)
!-----v momentum
  do j = 2, nj-1
    do i = 2, ni
      vel1 = 0.5_dp*(uo(i,j)+uo(i,j+1))
      vel2 = 0.5_dp*(uo(i-1,j)+uo(i-1,j+1))
      ae(i,j) = -0.5_dp*vel1/dx + rex
      aw(i,j) = 0.5_dp*vel2/dx + rex
      an(i,j) = -0.5_dp*vo(i,j+1)/dy + rey
      as(i,j) = -0.5_dp*vo(i,j-1)/dy + rey
      if (steadysim) then
        su(i,j) = (po(i,j)-po(i,j+1))/dy
        ap(i,j) = 2.0_dp*(rex+rey)
      else
        su(i,j) = vo(i,j)/dt + (po(i,j)-po(i,j+1))/dy
        ap(i,j) = pterm
      end if
      ! under-relax
      su(i,j) = su(i,j) + (1-urfv)*ap(i,j)/urfv*vo(i,j)
      ap(i,j) = ap(i,j)/urfv
    end do
  end do
  call sipsol(1, ni+1, 1, nj, -aw, -ae, -as, -an, ap, su, v)
!-----p' from continuity
  pp(:,:) = 0.0_dp
  do j = 2, nj
    do i = 2, ni
      if (steadysim) then
        ae(i,j) = 1.0_dp/(2.0_dp*(rex+rey))/dx**2 
        aw(i,j) = 1.0_dp/(2.0_dp*(rex+rey))/dx**2 
        an(i,j) = 1.0_dp/(2.0_dp*(rex+rey))/dy**2 
        as(i,j) = 1.0_dp/(2.0_dp*(rex+rey))/dy**2 
      else
        ae(i,j) = 1.0_dp/pterm/dx**2 
        aw(i,j) = 1.0_dp/pterm/dx**2 
        an(i,j) = 1.0_dp/pterm/dy**2 
        as(i,j) = 1.0_dp/pterm/dy**2 
      end if
      ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j)
      su(i,j) = (u(i-1,j)-u(i,j))/dx + (v(i,j-1)-v(i,j))/dy
    end do
  end do
  call sipsol(1, ni+1, 1, nj+1, -aw, -ae, -as, -an, ap, su, pp)
!-----correct u, v and p 
  do j = 2, nj
    do i = 2, ni
      ! under-relax
      p(i,j) =  p(i,j) + urfp*pp(i,j)
    end do
  end do
  do j = 2, nj
    do i = 2, ni-1
      u(i,j) = u(i,j) - (pp(i+1,j)-pp(i,j))/pterm/dx
    end do
  end do
  do j = 2, nj-1
    do i = 2, ni
      v(i,j) = v(i,j) - (pp(i,j+1)-pp(i,j))/pterm/dy
    end do
  end do
!-----boundary condition
  u(1:ni,1) = -u(1:ni,2)
  u(1:ni,nj+1) = 2.0_dp*utop - u(1:ni,nj)
  u(1,2:nj) = 0.0_dp
  u(ni,2:nj) = 0.0_dp

  v(1,1:nj) = -v(2,1:nj)
  v(ni+1,1:nj) = -v(ni,1:nj)
  v(2:ni,1) = 0.0_dp
  v(2:ni,nj) = 0.0_dp

  p(1,1:ni+1) = p(2,1:ni+1)
  p(ni+1,1:nj+1) = p(ni,1:nj+1)
  p(2:ni,nj+1) = p(2:ni,nj)
  p(2:ni,1) = p(2:ni,2)
!-----convergence
  error = 0.0_dp
  do j = 2, nj
    do i = 2, ni
      ! calculate norm of continuity residual and flux over the central plane.
      resnorm = ((u(i,j)-u(i-1,j))/dx + (v(i,j)-v(i,j-1))/dy)**2
      error = error + resnorm
    end do
  end do
  error = sqrt(error)! RMS as a criterion for convergence
  if (mod(itsp, 1000) == 0) &
  write(*,'(a,i6,2(2x,a,es9.2))') &
          'Iter=', itsp, 'MassErr=', error, 'u(mon)=', u(imon,jmon)
  if (error < tol) exit
  uo(:,:) = u(:,:)
  vo(:,:) = v(:,:)
  po(:,:) = p(:,:)
end do
if (error > 1.0e-2_dp) then
  write(*,*) 'Convergence failed!'
  stop
end if
! values at grid points
do j = 1, nj
  do i = 1, ni
    uc(i,j) = 0.5_dp*(u(i,j)+u(i,j+1))
    vc(i,j) = 0.5_dp*(v(i,j)+v(i+1,j))
    pc(i,j) = 0.25_dp*(p(i,j)+p(i+1,j)+p(i,j+1)+p(i+1,j+1)) 
  end do
end do
! screen write out
write(*,'(/,a,/)') 'Lid Driven Cavity Flow'
write(*,*) '  Finite difference method'
write(*,*) '  CDS for advection and diffusion terms'
if (steadysim) then
  write(*,*) '  Steady state SIMPLE'
else
  write(*,*) '  Implicit Euler in time with SIMPLE'
end if
write(*,'(/,2a,i3,3x,a,i3)') 'Grid: ', 'ni = ', ni, 'nj = ', nj
write(*,'(a,1x,es9.2)') 'Re = ', re
write(*,'(a,1x,f5.2)') 'Top u-velocity = ', utop
write(*,'(a,i3,a,i3,a,3(1x,es12.4))') 'x, y and u at(', imon, ',', jmon, &
        ') = ', x(imon), y(jmon), uc(imon,jmon)

! write to tecplot
open(unit=1, file='result.dat', status='replace')
write(1,*) 'title="Lid Driven Cavity Flow"'
write(1,*) 'variables="x", "y", "u", "v", "p"'
write(1,'(a,2x,a,i3,2x,a,i3,2x,a)') 'zone', 'i=', ni, 'j=', nj, 'f=point'
do j = 1, nj
  do i = 1, ni
    write(1,'(*(1x,es14.7))') x(i), y(j), uc(i,j), vc(i,j), pc(i,j)
  end do
end do
open(unit=2, file='uvel.txt', status='replace')
do j = 1, nj
  write(2,'(*(1x,es14.7))') x(imon), y(j), uc(imon,j)
end do
open(unit=3, file='vvel.txt', status='replace')
do i = 1, ni
  write(3,'(*(1x,es14.7))') x(i), y(jmon), vc(i,jmon)
end do

deallocate(u, v, p, pp, uc, vc, pc, stat=ierr)

close(1) 
close(2) 
close(3) 

stop
end program main
!*****************************************************************************
