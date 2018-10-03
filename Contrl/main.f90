! Copyright: Lawrence Berkeley Laboratory.
      program main
      use AccSimulatorclass
      implicit none
      include 'mpif.h'
      double precision :: time

      call construct_AccSimulator(time)
      call run_AccSimulator()
      call destruct_AccSimulator(time)

      end program main
