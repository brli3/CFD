program tec
! This program is to demenstrate the data format
! of ASCII files for Tecplot.
! Ruipengyu Li
! Modified 25/12/2017

implicit none
integer :: i, j, ni, nj, ierr
integer, parameter :: dp = selected_real_kind(15)
real(dp) :: dx, dy, dt, time
real(dp), allocatable :: x(:), y(:), p(:,:), t(:,:)
logical :: zone1 = .true.

ni = 3
nj = 3
allocate(x(ni), y(nj), p(ni,nj), t(ni,nj), stat=ierr)
if (ierr /= 0) print*, 'ERROR! Allocate array'
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

! Create a single line of text
open(unit=1, file='xy_line.dat', status='replace', iostat=ierr)
if (ierr /= 0) print*, 'ERROR! Open file 1'
write(1,100) min(ni,nj)
100 format('TITLE="Simple Line Text"', /, 'VARIABLES="X" "Y"', /, &
           'ZONE', 1x, /, 'I=', I3, 1x, 'DATAPACKING=POINT')
do i = 1, min(ni,nj)
  write(1,'(*(es11.2))') x(i), y(i)
end do

! Create a grid file
open(unit=2, file='grid.dat', status='replace', iostat=ierr)
if (ierr /= 0) print*, 'ERROR! Open file 2'
write(2,101) ni, nj, 1
101 format('TITLE="Grid File 2D"', /, 'FILETYPE=GRID', 1x, 'VARIABLES="X" "Y"', /, &
           'ZONE', /, 'I=', I3, 1x, 'J=', I3, 1x, 'K=', I3, /, &
           'ZONETYPE=Ordered', 1x, 'DATAPACKING=BLOCK')
do j = 1, nj         
  do i = 1, ni
    write(2,'(es11.2)') x(i)
  end do
end do
do j = 1, nj
  do i = 1, ni
    write(2,'(es11.2)') y(j)
  end do
end do

! Create a solution file
! N.B. load grid and solution files on Tecplot at the same time
open(unit=3, file='solu.dat', status='replace', iostat=ierr)
if (ierr /= 0) print*, 'ERROR! Open file 3'
write(3,102) ni, nj, 1
102 format('TITLE="Solution File 2D"', /, 'FILETYPE=SOLUTION', 1x, 'VARIABLES="Pressure"', /, &
           'ZONE', /, 'I=', I3, 1x, 'J=', I3, 1x, 'K=', I3, /, &
           'ZONETYPE=Ordered', 1x, 'DATAPACKING=BLOCK')
do j = 1, nj
  do i = 1, ni
    p(i,j) = x(i) * y(j)
    write(3,'(es11.2)') p(i,j)
  end do
end do

! Create a full transient data file
open(unit=4, file='trans.dat', status='replace', iostat=ierr)
if (ierr /= 0) print*, 'ERROR! Open file 4'
write(4,103)
time = 0.0_dp
dt = 0.2_dp
do while (time <= 1.0_dp)
  write(4,104) ni, nj, 1, time
  do j = 1, nj
    do i = 1, ni
      write(4,'(es11.2)') x(i)
    end do
  end do
  do j = 1, nj
    do i = 1, ni
      write(4,'(es11.2)') y(j)
    end do
  end do
  do j = 1, nj
    do i = 1, ni
      t(i,j) = (x(i) + y(j)) * time
      write(4,'(es11.2)') t(i,j)
    end do
  end do
  time = time + dt
end do
103 format('TITLE="Full Transient Data File"', /, 'FILETYPE=FULL', /, &
           'VARIABLES="X", "Y", "Temperatue"', /)
104 format(/, 'ZONE', 1x 'I=', I3, 1x, 'J=', I3, 1x, 'K=', I3, /, &
           'SOLUTIONTIME=', es11.2, /, 'ZONETYPE=Ordered', 1x, 'DATAPACKING=BLOCK')

! Create transient solution file to be used with a transient grid file
! N.B. SOLUTIONTIME cannot be specified for current Tecplot version
open(unit=5, file='grid_trans.dat', status='replace', iostat=ierr)
if (ierr /= 0) print*, 'ERROR! Open file 5'
open(unit=6, file='solu_trans.dat', status='replace', iostat=ierr)
if (ierr /= 0) print*, 'ERROR! Open file 6'
write(5,105)
write(6,106)
time = 0.0_dp
dt = 0.2_dp
do while (time <= 1.0_dp)
  write(5,107) ni, nj, 1
  write(6,107) ni, nj, 1
  do j = 1, nj
    do i = 1, ni
      write(5,'(es11.2)') x(i)
    end do
  end do
  do j = 1, nj
    do i = 1, ni
      write(5,'(es11.2)') y(j)
    end do
  end do
  do j = 1, nj
    do i = 1, ni
      t(i,j) = (x(i) + y(j)) * time
      write(6,'(es11.2)') t(i,j)
    end do
  end do
  time = time + dt
end do
105 format('TITLE="Transient Grid File"', /, 'FILETYPE=GRID', /, &
           'VARIABLES="X","Y"', /)
106 format('TITLE="Transient Solution File"', /, 'FILETYPE=SOLUTION', /, &
           'VARIABLES="Temperatue"', /)
107 format(/, 'ZONE', 1x 'I=', I3, 1x, 'J=', I3, 1x, 'K=', I3, /, &
           /, 'ZONETYPE=Ordered', 1x, 'DATAPACKING=BLOCK')


! Create full transient data file with shared xy locations
open(unit=7, file='trans_sharexy.dat', status='replace', iostat=ierr)
if (ierr /= 0) print*, 'ERROR! Open file 7'
write(7,103)
time = 0.0_dp
dt = 0.2_dp
do while (time <= 1.0_dp)
  if (zone1) then
    write(7,104) ni, nj, 1, time
    do j = 1, nj
      do i = 1, ni
        write(7,'(es11.2)') x(i)
      end do
    end do
    do j = 1, nj
      do i = 1, ni
        write(7,'(es11.2)') y(j)
      end do
    end do
    zone1 = .false.
  else
    write(7,108) ni, nj, 1, time
  end if
  do j = 1, nj
    do i = 1, ni
      t(i,j) = (x(i) + y(j)) * time
      write(7,'(es11.2)') t(i,j)
    end do
  end do
  time = time + dt
end do
108 format(/, 'ZONE', 1x 'I=', I3, 1x, 'J=', I3, 1x, 'K=', I3, /, &
           'SOLUTIONTIME=', es11.2, /, 'VARSHARELIST=([1-2]=1)', /, &
           'ZONETYPE=Ordered', 1x, 'DATAPACKING=BLOCK')

close(1)
close(2)
close(3)
close(4)
close(5)
close(6)
close(7)
stop
end program tec
