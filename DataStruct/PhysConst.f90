!----------------------------------------------------------------
! (c) Copyright, 2001 by the Regents of the University of California.
! PhysConstclass: physical constants class in CONSTANTS module of 
!                 DATA STRUCTURE layer.
! Version: 1.0
! Author: Ji Qiang, LANL, 3/7/01
! Description: This class defines the physical constant parameters used
!              in the simulation.
! Comments:
!----------------------------------------------------------------
      module PhysConstclass
!        use mpistub
        implicit none
      
        !physical parameters and constants ---------------------
        double precision :: Pi
        double precision :: Clight !speed of light in vaccum
        double precision :: Scxl !length scale
        double precision :: Rad2deg !conversion factor from radian to degree
        double precision :: Epsilon0 !permittivity of vacuum
        double precision :: Scfreq !time scale
        !<<<<<<<<<<< Kilean <<<<<<<<<<<<<<<<<
        double precision,parameter :: proton_mass = 938272081d0, &
        &electron_mass = 510998.94d0,&
        &eV_2_kg = 1.782662d-36
        
      contains
        subroutine construct_PhysConst(freq)
        double precision, intent(in) :: freq

        Clight = 299792458.0d0 
        Pi = 2*dasin(1.0d0)
        Scxl = Clight/(2*pi*freq)
        Rad2deg = 180.0d0/Pi
        Epsilon0 = 8.854187817d-12
        Scfreq = freq

        end subroutine construct_PhysConst
 
      end module PhysConstclass
