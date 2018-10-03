module parallel
!      include 'mpif.h'
! # of processors, processor id
      integer ,save :: nvp=1,idproc=0
! current communicator, MPI datatypes
      integer ,save :: lworld,mreal,mcplx,m2real,mntgr
! MPI operators
      integer ,save :: mpisum, mpimax, mpimaxloc
! the MPI implementation may require more than one error code
      integer ,save ,dimension(MPI_STATUS_SIZE) :: mpistat

   contains 

      subroutine init_parallel
      logical flag

!      call MPI_INITIALIZED(flag,mpistat)
!      if(.not.flag)then
!        call MPI_INIT(mpistat)
!      endif
      lworld=MPI_COMM_WORLD
! how many processors?
      call MPI_COMM_SIZE(lworld,nvp,mpistat)
! get processor id
      call MPI_COMM_RANK(lworld,idproc,mpistat)
      print*,"idproc in init: ",idproc
! in the future, these can be set differently for other precision:
      mreal=MPI_DOUBLE_PRECISION
      mcplx=MPI_DOUBLE_COMPLEX
      m2real=MPI_2DOUBLE_PRECISION
      mntgr=MPI_INTEGER
!     mntgr=MPI_2INTEGER
      mpisum=MPI_SUM
      mpimax=MPI_MAX
      mpimaxloc=MPI_MAXLOC
      end subroutine init_parallel
!
      subroutine end_parallel
      include 'mpif.h'
      logical flag
      
      call MPI_INITIALIZED(flag,mpistat)
      if(flag)then
        call MPI_BARRIER(lworld,mpistat)
        call MPI_FINALIZE(mpistat)
      endif
      end subroutine end_parallel

end module parallel
