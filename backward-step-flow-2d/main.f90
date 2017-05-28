!**********************************************************************
!  teach-t (modified in F90)
!  Ruipengyu Li  
!  April-2016
!**********************************************************************
program teach 
!----------------------------------------------------------------------
! solves conservation equations for heat, mass, momentum, etc.
! uses primitive variables with the velocities and pressure derived
! from SIMPLE algorithm and all equations are solved by line by line 
! method of TDMA.
! solve relevant conservation equations by means of a hybrid scheme.
! covers steady, incompressible, turbulent, plane or axisymetric flows.
!----------------------------------------------------------------------
   implicit none 
      call contro()
end program teach 
!**********************************************************************
