!--------------------------------------------------------------------------
! This program solves one-dimensional convection-diffusion
! equation:
!
!   d (rho*u*phi) / dx = d (Gamma*d(phi)/dx) / dx
!
! by finite difference method
! with Dirichlet boundary conditions on both ends.
!
! The discrete equation can be written as
!   
!  ae(i)*phi(i+1) + aw(i)*phi(i-1) + ap(i)*phi(i) = su(i) 
!
! Discretization schemes:
!   diffusion term - Cetral Differencing Scheme (CDS)
!   convection term - Upwind Differencing Scheme (UDS) or CDS
!
! We can demenstrate the oscillations using CDS and false diffusion using UDS.
! And how they are reduced with grid refinement.
!
! Oscillations related to CDS can be reduced by ensuring that 
! the local Pe is not larger than 2 in regions of high gradient change of phi.
! This can be demenstrated using a non-uniform grid.
!
! Finally, we can prove the order of accuracy of numerical schemes by
! calculating the average error with reduction of grid spacing
!
! Author: Ruipengyu Li 
! Modified: 22/04/2017
!
! Reference:
!   J. H. Ferziger and M. Peric, Computational Methods for Fluid Dynamics,
!   3rd ed. Springer Berlin Heidelberg, 2001.
!--------------------------------------------------------------------------
module conv_diff_mod

implicit none

integer, parameter :: dp = selected_real_kind(15) ! Double precision

contains

subroutine tdma(a, b, c, d, x)
!  Tri-diagonol system of equations
!|a1 -b1                    | | x1 | | d1 |
!|-c2 a2 -b2                | | x2 | | d2 |
!|   -c3 a3 -b3             | | x3 | | d3 |
!|       .. .. ..           |*| .. |=| .. |
!|          .. .. ..        | | .. | | .. |
!|          -cn-1 an-1 -bn-1| |xn-1| |dn-1|
!|                -cn   an  | | xn | | dn |
!  ith equation in the system:
!   a(i)x(i) = b(i)x(i+1) + c(i)x(i-1) + d(i)
implicit none 
integer :: i, ierr
integer :: n
real(dp), dimension(:), intent(in) :: a, b, c, d
real(dp), dimension(:), intent(out) :: x
real(dp), allocatable, dimension(:) :: p, q
n = size(a(:))
allocate(p(1:n), q(1:n), stat=ierr)
!-----Forward elimination
! 1 and n are boundary values.
p(1) = 0.0_dp
q(1) = x(1)
do i = 2, n-1
  p(i) = b(i) / (a(i)-c(i)*p(i-1))
  q(i) = (d(i)+c(i)*q(i-1)) / (a(i)-c(i)*p(i-1))
end do
!-----Back substitution
do i = n-1, 2, -1
  x(i) = p(i)*x(i+1) + q(i)
end do
deallocate(p, q, stat=ierr)
end subroutine tdma

function PHI_EXACT(x, pe, l, phi_i, phi_o)
! Exact solution of 1-D advection-diffusion equation.
implicit none
real(dp), intent(in) :: x, pe, l, phi_i, phi_o
real(dp) :: PHI_EXACT

PHI_EXACT = phi_i + (phi_o-phi_i)*(EXP(x*pe/l)-1.0_dp)/(EXP(pe)-1.0_dp)

end function PHI_EXACT

end module conv_diff_mod

program fd1dscd
! Convection term: use upwind or central difference scheme.
! Diffusion term: use central difference scheme.
use conv_diff_mod 
implicit none
integer :: i, ni, nim1
integer :: ierr, csm
real(dp) :: den, u, gam, pe, phi_in, phi_out, err
real(dp) :: xmin, xmax, dx, expf
real(dp) :: cw, ce, dw, de
real(dp), allocatable, dimension(:) :: x, phi, phi_ex
real(dp), allocatable, dimension(:) :: aw, ae, ap, su
!-----Initialize variables
ni = 11 
nim1 = ni - 1
den = 1.0_dp  ! density
u = 1.0_dp  ! velocity, assume constant
gam = 0.02_dp  ! diffusion coefficient, assume constant
phi_in = 0.0_dp  ! left boundary value
phi_out = 1.0_dp ! right boundary value
xmin = 0.0_dp
xmax = 1.0_dp  ! length of domain
expf = 1.0_dp  ! grid expansion factor
csm = 1  ! convection scheme: 1.UDS 2.CDS
!-----Initialize arrays
allocate(x(ni), phi(ni), phi_ex(ni), stat=ierr)
allocate(aw(ni), ae(ni), ap(ni), su(ni), stat=ierr)
aw = 0.0_dp
ae = 0.0_dp
ap = 0.0_dp
su = 0.0_dp
x = 0.0_dp
phi = 0.0_dp
phi_ex = 0.0_dp
!-----Define grid
if (expf==1.0_dp) then
  ! uniform grid
  dx = (xmax-xmin) / (ni-1)
else
  ! non-uniform grid
  ! a little math here: need to know sum of a geometric series
  ! to get the first dx.
  dx = (xmax-xmin) * (1.0_dp-expf) / (1.0_dp-expf**(ni-1))
end if
x(1) = xmin
do i = 2, ni
  x(i) = x(i-1) + dx
  dx = dx * expf
end do
!-----Boundary values
phi(1) = phi_in
phi(ni) = phi_out
!-----Compute coefficients of algebraic equations at each point
do i = 2, nim1
  if (csm==1) then
!-----Upwind advection
    ce = min(den*u, 0.0_dp) / (x(i+1)-x(i))
    cw = -max(den*u, 0.0_dp) / (x(i)-x(i-1))
  else
!-----Central difference advection
    ce = den*u / (x(i+1)-x(i-1))
    cw = -den*u / (x(i+1)-x(i-1)) 
  end if
!-----Central difference diffusion
  de = -2.0_dp*gam / ((x(i+1)-x(i-1)) * (x(i+1)-x(i)))
  dw = -2.0_dp*gam / ((x(i+1)-x(i-1)) * (x(i)-x(i-1)))
!-----Assemble coefficient
  ae(i) = (ce + de)
  aw(i) = (cw + dw)
  ap(i) = aw(i) + ae(i)
end do 
!-----Boundary nodes treatment 
! included in the solver
!su(2) = su(2) + aw(2)*phi(1)
!aw(2) = 0.0_dp
!su(nim1) = su(nim1) + ae(nim1)*phi(ni)
!ae(nim1) = 0.0_dp
!-----Solve
call tdma(ap, ae, aw, su, phi)
!-----Exact solution
pe = den*u*xmax/gam
err = 0.0_dp
do i = 1, ni
  phi_ex(i) = PHI_EXACT(x(i), pe, xmax-xmin, phi_in, phi_out)
  err = err + ABS(phi(i) - phi_ex(i))
end do
err = err / ni
!-----Write out
write(*,*) '   One-dimension convection-diffusion problem'
write(*,*)
write(*,*) '   Peclet number: ', pe
write(*,*) '   Error avg: ', err
! coefficients
write(*,'(/,t8,4(a2,13x))') 'aw', 'ae', 'ap', 'su'
do i = 1, ni
  write(*,'(4(f15.3))') aw(i), ae(i), ap(i), su(i)
end do
! phi
write(*,'(/,t8,a2,10x,a3,6x,a9,5x,a5)') 'x', 'phi', &
                                        'phi_exact', 'Error'
do i = 1, ni
  write(*,'(4(f12.5))') x(i), phi(i), phi_ex(i), &
                        phi_ex(i)-phi(i)
end do
deallocate(aw, ae, ap, su)
deallocate(x, phi, phi_ex)
end program fd1dscd
!********************************************************************************
