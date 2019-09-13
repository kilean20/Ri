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
  integer, private, parameter :: rectangular_=1, elliptic_=2
  type, private :: PipeInfoType
    double precision, allocatable, dimension(:) :: s
    double precision, allocatable, dimension(:) :: x
    double precision, allocatable, dimension(:) :: y
    integer, allocatable, dimension(:) :: shape
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
  allocate(PipeInfo%s(PipeInfo%n),PipeInfo%x(PipeInfo%n),PipeInfo%y(PipeInfo%n),PipeInfo%shape(PipeInfo%n))
  
  if(myrank == 0) then
    open(iUnit,file='pipeinfo.in',status='old',action='read')
    loop2 : DO i=1,PipeInfo%n
      READ(iUnit,*,iostat=eastat) PipeInfo%s(i), PipeInfo%x(i), PipeInfo%y(i), PipeInfo%shape(i)
    END DO loop2
    close(iUnit)

  endif 
  call MPI_Bcast(PipeInfo%s,PipeInfo%n,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(PipeInfo%x,PipeInfo%n,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(PipeInfo%y,PipeInfo%n,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
  
end subroutine readPipeInfo


subroutine getPipeInfo(z, &
                       pipe_shape0,pipe_x0,pipe_y0,pipe_z0, &
                       pipe_shape1,pipe_x1,pipe_y1,pipe_z1)
  implicit none
  double precision, intent(in) :: z
  integer, intent(out) :: pipe_shape0,pipe_shape1
  double precision, intent(out) :: pipe_x0,pipe_y0,pipe_z0,pipe_x1,pipe_y1,pipe_z1
  integer :: i
  
  loop1 : DO i=1,PipeInfo%n
    if ( z > PipeInfo%s(i) ) then
      pipe_z0 = PipeInfo%s(i)
      pipe_x0 = PipeInfo%x(i)
      pipe_y0 = PipeInfo%y(i)
      pipe_z1 = PipeInfo%s(i+1)
      pipe_x1 = PipeInfo%x(i+1)
      pipe_y1 = PipeInfo%y(i+1)
      pipe_shape0 = PipeInfo%shape(i)
      pipe_shape1 = PipeInfo%shape(i+1)
      exit
    endif
  END DO loop1
  
end subroutine getPipeInfo


subroutine getLossInfo(x,y,z,&
                       pipe_shape0,pipe_x0,pipe_y0,pipe_z0, &
                       pipe_shape1,pipe_x1,pipe_y1,pipe_z1, &
                       flagLost)
  implicit none
  double precision, intent(in) :: x,y,z,pipe_x0,pipe_y0,pipe_z0,pipe_x1,pipe_y1,pipe_z1
  integer, intent(in) :: pipe_shape0,pipe_shape1
  logical, intent(out) :: flagLost
  double precision :: r,rx,ry,pipe_r,pipe_r0,pipe_r1
  r = sqrt(x**2 + y**2)
  rx = x/pipe_x0
  ry = y/pipe_y0
  if(pipe_shape0==rectangular_) then
    if(abs(rx) >= abs(ry)) then
      pipe_r0 = r*(pipe_x0/abs(x))
    else
      pipe_r0 = r*(pipe_y0/abs(y))
    endif
  else
    pipe_r0 = r/sqrt(rx**2 + ry**2)
  endif
  rx = x/pipe_x1
  ry = y/pipe_y1
  if(pipe_shape0==rectangular_) then
    if(abs(rx) >= abs(ry)) then
      pipe_r1 = r*(pipe_x0/abs(x))
    else
      pipe_r1 = r*(pipe_y0/abs(y))
    endif
  else
    pipe_r1 = r/sqrt(rx**2 + ry**2)
  endif
  pipe_r = (pipe_r1*(z-pipe_z0) + pipe_r0*(pipe_z1-z))/(pipe_z1-pipe_z0)
  if(r > pipe_r) then
    flagLost = .True.
  else
    flagLost = .False.
  endif
end subroutine

end module