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
  type :: PipeInfoType
    double precision, allocatable, dimension(:) :: s
    double precision, allocatable, dimension(:) :: x
    double precision, allocatable, dimension(:) :: y
    integer, allocatable, dimension(:) :: shape
    integer :: n
  end type PipeInfoType
contains


subroutine readPipeInfo(PipeInfo)
  implicit none
  include 'mpif.h'
  type(PipeInfoType), intent(out) :: PipeInfo
  integer :: iUnit,eastat,i,myrank,ierr
  double precision :: tmp
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
  !print*, 'myrank,PipeInfo%n=',myrank,PipeInfo%n
  allocate(PipeInfo%s(PipeInfo%n),PipeInfo%x(PipeInfo%n),PipeInfo%y(PipeInfo%n),PipeInfo%shape(PipeInfo%n))
  
  if(myrank == 0) then
    open(iUnit,file='pipeinfo.in',status='old',action='read')
    loop2 : DO i=1,PipeInfo%n
      READ(iUnit,*,iostat=eastat) PipeInfo%s(i), PipeInfo%x(i), PipeInfo%y(i), tmp
       PipeInfo%shape(i) = tmp
    END DO loop2
    close(iUnit)

  endif 
  call MPI_Bcast(PipeInfo%s,PipeInfo%n,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(PipeInfo%x,PipeInfo%n,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(PipeInfo%y,PipeInfo%n,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(PipeInfo%shape,PipeInfo%n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  
  
  if(myrank == 0) then
    print*, "---- pipe_info: s, x, y, shape ----"
    loop3 : DO i=1,PipeInfo%n
      print*, PipeInfo%s(i), PipeInfo%x(i), PipeInfo%y(i), PipeInfo%shape(i)
    END DO loop3
  endif
end subroutine readPipeInfo


subroutine getPipeInfo(PipeInfo,z,&
                       pipe_shape0,pipe_x0,pipe_y0,pipe_z0, &
                       pipe_shape1,pipe_x1,pipe_y1,pipe_z1)
  implicit none
  type(PipeInfoType), intent(in) :: PipeInfo
  double precision, intent(in) :: z
  integer, intent(out) :: pipe_shape0,pipe_shape1
  double precision, intent(out) :: pipe_x0,pipe_y0,pipe_z0,pipe_x1,pipe_y1,pipe_z1
  integer :: i
  
  
  loop1 : DO i=2,PipeInfo%n
    if ( z <= PipeInfo%s(i) ) then
      pipe_z0 = PipeInfo%s(i-1)
      pipe_x0 = PipeInfo%x(i-1)
      pipe_y0 = PipeInfo%y(i-1)
      pipe_z1 = PipeInfo%s(i)
      pipe_x1 = PipeInfo%x(i)
      pipe_y1 = PipeInfo%y(i)
      pipe_shape0 = PipeInfo%shape(i-1)
      pipe_shape1 = PipeInfo%shape(i)
      exit
    endif
  END DO loop1  
end subroutine getPipeInfo


subroutine getLossInfo(PipeInfo,x,y,z,&
                       pipe_shape0,pipe_x0,pipe_y0,pipe_z0, &
                       pipe_shape1,pipe_x1,pipe_y1,pipe_z1, &
                       flagLost)
  implicit none
  type(PipeInfoType), intent(in) :: PipeInfo
  double precision, intent(in) :: x,y,z,pipe_x0,pipe_y0,pipe_z0,pipe_x1,pipe_y1,pipe_z1
  integer, intent(in) :: pipe_shape0,pipe_shape1
  logical, intent(out) :: flagLost
  double precision :: r,rx,ry,pipe_r,pipe_r0,pipe_r1
  
  flagLost = .False.
  
  if(pipe_shape0==pipe_shape1) then
    if(pipe_shape0==rectangular_) then
      if(pipe_x0 == pipe_x1 .and. pipe_y0 == pipe_y1) then
        if (abs(x) >= pipe_x0 .or. abs(y) >= pipe_y0) then
          flagLost = .True.
        endif
      else 
        pipe_r0 = (pipe_x1*(z-pipe_z0) + pipe_x0*(pipe_z1-z))/(pipe_z1-pipe_z0)
        pipe_r1 = (pipe_y1*(z-pipe_z0) + pipe_y0*(pipe_z1-z))/(pipe_z1-pipe_z0)
        if (abs(x) >= pipe_r0 .or. abs(y) >= pipe_r1) then
          flagLost = .True.
        endif
      endif
    else if(pipe_shape0==elliptic_) then 
      if(pipe_x0 == pipe_x1 .and. pipe_y0 == pipe_y1) then
        if( ((x/pipe_x0)**2 + (y/pipe_y0)**2) >= 1.0) then
          flagLost = .True.
        endif
      else
        pipe_r0 = (pipe_x1*(z-pipe_z0) + pipe_x0*(pipe_z1-z))/(pipe_z1-pipe_z0)
        pipe_r1 = (pipe_y1*(z-pipe_z0) + pipe_y0*(pipe_z1-z))/(pipe_z1-pipe_z0)
        if( ((x/pipe_r0)**2 + (y/pipe_r1)**2) >= 1.0) then
          flagLost = .True.
        endif
      endif
    endif
  else
    !print*, "pipe_shape0!=pipe_shape1,pipe_z0,1 and shape0,1=",pipe_z0,pipe_z1,pipe_shape0,pipe_shape1
    flagLost = .False.
  endif
end subroutine getLossInfo

end module