program tec
! This program is to demenstrate the data format for Tecplot.
! Ruipengyu Li
! Modified 24/12/2017

implicit none
integer :: i, j, ni, nj, ierr
integer, parameter :: dp = selected_real_kind(15)
real(dp) :: dx, dy
real(dp), allocatable, dimension(:) x, y, p, t

ni = 5
nj = 5
allocate(x(ni), y(nj), p(ni,nj), t(ni,nj), status=ierr)
if (ierr /= 0) print*, 'ERROR! Allocate arrays'
x(:) = 0.0_dp
y(:) = 0.0_dp
p(:,:) = 0.0_dp
t(:,:) = 0.0_dp
dx = 1.0_dp / real(ni-1)
dy = 1.0_dp / real(nj-1)
do i = 2, ni
  x(i) = x(i-1) + dx
end do
do j = 2, nj
  y(j) = y(j-1) + dy
end do
do j = 1, nj
  do i = 1, ni
    p(i,j) = x(i) * y(j)
    t(i,j) = x(i) + y(j)
  end do
end do
open(unit=1, file='output.dat', status='replace', iostat=ierr)
if (ierr /= 0) print*, 'ERROR! Open file'
100 format(





close(1)
stop
end program tec
