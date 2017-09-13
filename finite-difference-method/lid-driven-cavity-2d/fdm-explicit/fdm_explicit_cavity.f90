! Lid-driven cavity flow
!
! A short program to demenstrate SIMPLE algorithm
!
!   Finite difference method
!   CDS for advection and diffusion terms
!   Explicit Euler time marching (iteration)
!   SIMPLE algorithm
!
!   Note that this method requires a significant amount of computational
!   time for relatively large Reynolds number and fine mesh.
!
! Ruipengyu Li
! Modified: 11/08/2017
!*****************************************************************************
module cavity_mod

implicit none

integer, parameter :: dp = selected_real_kind(15)

contains

subroutine sipsol(i1, i2, j1, j2, aw, ae, as, an, ap, su, phi)
! SIP solver, ILU of Stone (1968)
implicit none
integer :: i, j, ierr, iter
integer, intent(in) :: i1, i2, j1, j2
integer, parameter :: maxit = 600
real(dp), parameter :: alpha = 0.88_dp
real(dp), parameter :: tol = 1.0e-5_dp
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

end module cavity_mod
!*****************************************************************************
program main
use cavity_mod 
implicit none

integer :: i, j, ni, nj, imon, jmon
integer :: itsp, ntsp
integer :: ierr
character(len=80) :: msg
real(dp) :: dx, dy, dt 
real(dp) :: utop, re ! Reynolds no.
real(dp) :: urfp, resnorm, error, tol
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
ntsp = 500000
utop = 1.0_dp
re = 10.0_dp
dt = 1.0e-4_dp
urfp = 0.8_dp ! under-relaxation factor for pressure
tol = 1.0e-15_dp

allocate(u(1:ni,1:nj+1), v(1:ni+1,1:nj), p(1:ni+1,1:nj+1), &
         uo(1:ni,1:nj+1), vo(1:ni+1,1:nj), po(1:ni+1,1:nj+1), &
         pp(1:ni+1,1:nj+1), &
         uc(1:ni,1:nj), vc(1:ni,1:nj), pc(1:ni,1:nj), &
         x(1:ni), y(1:nj), &
         stat=ierr, errmsg=msg)
         
allocate(ae(1:ni+1,1:nj+1), aw(1:ni+1,1:nj+1), an(1:ni+1,1:nj+1), &
         as(1:ni+1,1:nj+1), ap(1:ni+1,1:nj+1), su(1:ni+1,1:nj+1), &
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
x(:) = 0.0_dp      
y(:) = 0.0_dp      

ae(:,:) = 0.0_dp
aw(:,:) = 0.0_dp
an(:,:) = 0.0_dp
as(:,:) = 0.0_dp
ap(:,:) = 0.0_dp
su(:,:) = 0.0_dp

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
po(1:ni+1,1:nj+1) = 1.0_dp

do itsp = 1, ntsp
!-----u momentum
  do j = 2, nj
    do i = 2, ni-1
      u(i,j) = -0.5_dp*(uo(i+1,j)*uo(i+1,j) - uo(i-1,j)*uo(i-1,j))/dx &
              - 0.25_dp*((uo(i,j)+uo(i,j+1))*(vo(i,j)+vo(i+1,j)) &
              -          (uo(i,j)+uo(i,j-1))*(vo(i,j-1)+vo(i+1,j-1)))/dy &
              - (po(i+1,j)-po(i,j))/dx &
              + ((uo(i-1,j)-2.0_dp*uo(i,j)+uo(i+1,j))/dx**2 &
              +  (uo(i,j-1)-2.0_dp*uo(i,j)+uo(i,j+1))/dy**2)/re
      u(i,j) = uo(i,j) + dt*u(i,j)
    end do
  end do
!-----v momentum
  do j = 2, nj-1
    do i = 2, ni
      v(i,j) = -0.25_dp*((uo(i,j)+uo(i,j+1))*(vo(i,j)+vo(i+1,j)) &
              -          (uo(i-1,j)+uo(i-1,j+1))*(vo(i,j)+vo(i-1,j)))/dx &
              -0.5_dp*(vo(i,j+1)*vo(i,j+1) - vo(i,j-1)*vo(i,j-1))/dy &
              - (po(i,j+1)-po(i,j))/dy &
              + ((vo(i-1,j)-2.0_dp*vo(i,j)+vo(i+1,j))/dx**2 &
              +  (vo(i,j-1)-2.0_dp*vo(i,j)+vo(i,j+1))/dy**2)/re
      v(i,j) = vo(i,j) + dt*v(i,j)
    end do
  end do
!-----p' from continuity
  pp(:,:) = 0.0_dp
  do j = 2, nj
    do i = 2, ni
      ! coefficients for p' equation
      ae(i,j) = dt/dx**2
      aw(i,j) = dt/dx**2
      an(i,j) = dt/dy**2
      as(i,j) = dt/dy**2
      su(i,j) = -(u(i,j)-u(i-1,j))/dx - (v(i,j)-v(i,j-1))/dy
      ap(i,j) = 2*(dt/dx**2 + dt/dy**2)
    end do
  end do
  call sipsol(1, ni+1, 1, nj+1, -aw, -ae, -as, -an, ap, su, pp)
!-----correct u, v and p 
  do j = 2, nj
    do i = 2, ni
      p(i,j) =  p(i,j) + urfp*pp(i,j)
    end do
  end do
  do j = 2, nj
    do i = 2, ni-1
      u(i,j) = u(i,j) - dt/dx*(pp(i+1,j)-pp(i,j))
    end do
  end do
  do j = 2, nj-1
    do i = 2, ni
      v(i,j) = v(i,j) - dt/dy*(pp(i,j+1)-pp(i,j))
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
  error = sqrt(error)
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
write(*,*) '  Explicit in time, SIMPLE algorithm'
write(*,'(/,2a,i3,3x,a,i3)') 'Grid: ', 'ni = ', ni, 'nj = ', nj
write(*,'(a,1x,es9.2)') 'Re = ', re
write(*,'(a,1x,es9.2)') 'dt = ', dt
write(*,'(a,1x,es9.2)') 'under-relax pressure = ', urfp
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
