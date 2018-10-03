!----------------------------------------------------------------
! (c) Copyright, 2001 by the Regents of the University of California.
! NumConstclass: numerical constants class in CONSTANTS module of
!                 DATA STRUCTURE layer.
! Version: 1.0
! Author: Ji Qiang, LANL, 3/7/01
! Description: This class defines the maximum size for numerical parameters
!              in the simulation.
! Comments: 
!----------------------------------------------------------------
      module NumConstclass
        implicit none
      
        !numerical parameters --------------------------
        integer, parameter :: Pdim = 6 !phase dimension
        !maximum # of particles on a single processor.
        integer, parameter :: Nplcmax = 5000000
        !maximum # of grids in x,y,z dimension on a single processor.
        integer, parameter :: Nxlcmax = 128
        integer, parameter :: Nylcmax = 128
        integer, parameter :: Nzlcmax = 128
        !maximum # of beam line elements
        integer, parameter :: Nblemtmax = 5000
        !maximum # of drift space
        integer, parameter :: Ndriftmax = 2400
        !maximum # of quadrupoles
        integer, parameter :: Nquadmax = 1000
        integer, parameter :: Ncfmax = 500
        !maximum # of dipoles
        integer, parameter :: Ndipolemax = 100
        !maximum # of rf gaps
        integer, parameter :: Ncclmax = 1000
        integer, parameter :: Nccdtlmax = 1000
        integer, parameter :: Ndtlmax = 1000
        integer, parameter :: Nscmax = 1000
        !maximum # of beam position monitors
        integer, parameter :: Nbpmmax = 1000
        !maximum # of magnetic solenoid
        integer, parameter :: Nslmax = 1000
        !maximum # of magnetic solenoid with RF field
        integer, parameter :: Nslrfmax = 1000
        !maximum # of nonlinear insert elements
        integer, parameter :: Nnllmax = 1000

      end module NumConstclass
