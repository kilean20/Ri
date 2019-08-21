!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! BeamBunchclass: Charged beam bunch class in Beam module of APPLICATION 
!                 layer.
! Version: beta
! Author: Kilean Hwang
! Description: Read  beam pipe info from external file
! Comments: 
!----------------------------------------------------------------
module PipeInfoClass
  type, private :: PipeInfoType
    double precision, allocatable, dimension(:) :: s
    double precision, allocatable, dimension(:) :: x
    double precision, allocatable, dimension(:) :: y
    integer :: n
  end type PipeInfoType
  type (PipeInfoType), public :: PipeInfo
contains


subroutine readPipeInfo()
  implicit none
  include 'mpif.h'
  integer :: iUnit,eastat,i,myrank,ierr
  logical :: file_open
  
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank,ierr);
  if(myrank == 0) then
    iUnit = 2132161
    file_open = .true. 
    DO while ( file_open ) 
      iUnit = iUnit + 1 
      inquire(iUnit, opened = file_open ) 
    end DO
  
    PipeInfo%n = 0
    open(iUnit,file='pipeinfo.in',status='old',action='read')
    loop1 : DO
      READ(iUnit,*,iostat=eastat)
      IF (eastat < 0) THEN
        EXIT loop1
      ELSE IF (eastat > 0) THEN
        STOP 'IO-error'
      ENDIF
      PipeInfo%n = PipeInfo%n+1
    END DO loop1
    close(iUnit)
  endif
  call MPI_Bcast(PipeInfo%n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  allocate(PipeInfo%s(PipeInfo%n),PipeInfo%x(PipeInfo%n),PipeInfo%y(PipeInfo%n))
  
  if(myrank == 0) then
    open(iUnit,file='pipeinfo.in',status='old',action='read')
    loop2 : DO i=1,PipeInfo%n
      READ(iUnit,*,iostat=eastat) PipeInfo%s(i), PipeInfo%x(i), PipeInfo%y(i)
    END DO loop2
    close(iUnit)

  endif 
  call MPI_Bcast(PipeInfo%s,PipeInfo%n,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(PipeInfo%x,PipeInfo%n,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(PipeInfo%y,PipeInfo%n,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
  
end subroutine readPipeInfo


subroutine getPipeInfo(z,piperad,piperad2)
  implicit none
  double precision, intent(in) :: z
  double precision, intent(out) :: piperad, piperad2
  integer :: i
  
  loop1 : DO i=PipeInfo%n,1,-1
    if ( z > PipeInfo%s(i) ) then
      piperad = PipeInfo%x(i)
      piperad2 = PipeInfo%y(i)
      exit
    endif
  END DO loop1
  
end subroutine getPipeInfo

end module