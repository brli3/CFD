program main
! Demo for velocity moving time average
! Useful for unsteady simulation
! Ruipengyu Li 15/12/2017
implicit none
integer, parameter :: dp = selected_real_kind(15)
integer :: i, ni
real :: dx, t, dt, t_end, t0, dt_av
real, allocatable, dimension(:) :: x, u, uav, ufluc 

open(unit=1, file='u_point.txt', status='replace')
open(unit=2, file='uav_point.txt', status='replace')

ni = 11
allocate(x(ni), u(ni), uav(ni), ufluc(ni))
u(:) = 0.0_dp
uav(:) = 0.0_dp
ufluc(:) = 0.0_dp
x(1) = 0.0_dp
dx = 1.0_dp / (ni-1)
do i = 2, ni
  x(i) = x(i-1) + dx
end do
write(*,*) 'Initial coordinate and velocity:'
do i = 1, ni
  u(i) = (x(i)-0.5_dp)**2 
  write(*,'(5f7.3)') x(i), u(i)
end do

t = 0.0_dp
t0 = t
t_end = 10.0_dp
dt = 0.1_dp
dt_av = 0.2_dp

write(*,'(/,t2,2(a,f6.2,5x),/)') 'dt =', dt, 'dt_av =', dt_av

do while (t < t_end)
  t = t + dt 
  write(*,'(2(a,f6.2,5x))') 't =', t, 't0 =', t0
  call random_number(u(:))
  u(:) = exp(u(:)) * abs(t - t_end/2.) 
  write(1,'(2(1x,es14.7))') t, u(ni/2)
  uav(:) = uav(:) + u(:) * dt
  if (t - t0 >= dt_av) then
    uav(:) = uav(:) / (t - t0)
    ufluc(:) = u(:) - uav(:)
    write(*,'(4(2x,a))') '  x', '    u', '   uav', ' ufluc'
    do i = 1, ni
      write(*,'(5f7.3)') x(i), u(i), uav(i), ufluc(i)
    end do
    write(2,'(2(1x,es14.7))') t, uav(ni/2)
    uav(:) = 0.0_dp
    t0 = t
    print *
  end if
  call random_number(dt)
  dt = 0.1_dp * dt
  if (t + dt > t_end) dt = t_end - t
end do

deallocate(x, u, uav, ufluc)
close(1)
close(2)
stop
end program main
