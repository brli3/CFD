! Lid-driven cavity flow
!
! A short program to demenstrate SIMPLE algorithm
!
!   Finite difference method
!   Cell-centred grid
!   Explicit Euler time marching (iteration)
!   SIMPLE algorithm

! Ruipengyu Li
! Modified: 07/08/2017
!*****************************************************************************
program main

implicit none

integer, parameter :: dp = selected_real_kind(15)
integer :: i, j, ni, nj
integer :: itsp, ntsp
integer :: ierr
character(len=80) :: msg
real(dp) :: dx, dy, dt 
real(dp) :: utop, re ! Reynolds no.
real(dp) :: udfp, resm, error, tol
real(dp), allocatable, dimension(:) :: x, y
real(dp), allocatable, dimension(:,:) :: uo, vo, po ! uncorrected
real(dp), allocatable, dimension(:,:) :: pp ! correction
real(dp), allocatable, dimension(:,:) :: u, v, p ! corrected
real(dp), allocatable, dimension(:,:) :: uc, vc, pc ! grid values

ni = 101
nj = 101
ntsp = 100000
utop = 1.0_dp
re = 500.0_dp
error = 0.0_dp
tol = 1.0e-5_dp

allocate(u(1:ni,1:nj+1), v(1:ni+1,1:nj), p(1:ni+1,1:nj+1), &
         uo(1:ni,1:nj+1), vo(1:ni+1,1:nj), po(1:ni+1,1:nj+1), &
         pp(1:ni+1,1:nj+1), &
         uc(1:ni,1:nj), vc(1:ni,1:nj), pc(1:ni,1:nj), &
         x(1:ni), y(1:nj), &
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

dx = 1.0_dp / (ni-1)
dy = 1.0_dp / (nj-1)
dt = 1.0e-4_dp
udfp = 0.8_dp ! under-relaxation factor for pressure

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
      pp(i,j) = (dt/dx**2*(pp(i+1,j)+pp(i-1,j)) &
               + dt/dy**2*(pp(i,j+1)+pp(i,j-1)) &
               - (u(i,j)-u(i-1,j))/dx-(v(i,j)-v(i,j-1))/dy) &
               / (2*(dt/dx**2 + dt/dy**2))
    end do
  end do
!-----correct u, v and p 
  do j = 2, nj
    do i = 2, ni
      p(i,j) =  p(i,j) + udfp*pp(i,j)
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
      resm = (u(i,j)-u(i-1,j))/dx + (v(i,j)-v(i,j-1))/dy
      error = error + resm**2
    end do
  end do
  error = sqrt(error) ! RMS as a criterion for convergence
  if (mod(itsp, 500) == 0) &
  write(*,'(a,i6,2x,a,es9.2)') 'Iter=', itsp, 'Error=', error
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
write(*,*) '  Cell-centred staggered grid'
write(*,*) '  SIMPLE algorithm'
write(*,'(/,2a,i3,3x,a,i3)') 'Grid: ', 'ni = ', ni, 'nj = ', nj
write(*,'(a,1x,es9.2)') 'Re = ', re
write(*,'(a,1x,f5.2)') 'Top u-velocity = ', utop

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
close(1) 

deallocate(u, v, p, pp, uc, vc, pc, stat=ierr)

stop
end program main
!*****************************************************************************
