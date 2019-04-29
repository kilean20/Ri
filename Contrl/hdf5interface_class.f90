! Interface bewteen IMPACTZ and hdf5io module for OpenPMD
! modification: 4/26/2019 Kilean
!  @ for impact particle data

module hdf5_interface_class
  use parallel_class
  use hdf5io_class

  implicit none

  type(parallel),target :: p
  class(parallel), pointer :: pp => null()
  type(hdf5file) :: file_hdf5
  real, dimension(:,:), allocatable :: arr
  real, dimension(:,:,:), allocatable :: arr3d
  real, dimension(:,:), allocatable :: part
  integer :: HDF5noff, HDF5nyp, HDF5Out,HDF5FieldOutFlag=0,HDF5iteration

  contains

  subroutine init_hdf5_interface(ny)
    integer, intent(in) :: ny
    ! Initialize MPI environment for HDF5
    call p%new()

    pp => p

    if (mod(ny,pp%getnvp()) /= 0) then
       print *, 'The number of processor cannot divde nz=', ny
       call p%del()
       call exit
    else
       HDF5nyp = ny/pp%getnvp()
    end if

    HDF5noff = pp%getidproc()*HDF5nyp
    !print*,HDF5nyp,HDF5noff
  end subroutine init_hdf5_interface

  subroutine hdf5_particle_output(BptPointer,Nplocal,nfile,iteration,mass,pName,normalF)
    double precision, dimension(:,:), intent(in) :: BptPointer
    double precision, intent(in) :: mass !kg
    character(len=*), intent(in) :: pName
    double precision, optional, intent(in) :: normalF(6)
    integer, intent(in) :: Nplocal,nfile,iteration
    integer :: ierr

    ! Initialize the file data for particle postion x
    call file_hdf5%new(iter=iteration,&
                &particleName=trim(pName),&
                &unitDimension=(/1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/),&
                &chiter=str(nfile),&
                &records='position',&
                &component='x')

    ! Write position x
    if(present(normalF)) then
      call pwpart(pp,file_hdf5,BptPointer(1,:)*normalF(1),Nplocal,ierr)
    else
      call pwpart(pp,file_hdf5,BptPointer(1,:),Nplocal,ierr)
    endif
    
     ! Initialize the file data for particle momentum x
    call file_hdf5%new(&
                &unitDimension=(/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/),&
                &records='moments',&
                &component='x')
    !print*,BptPointer(2,1:6)

    ! Write momentum x
    if(present(normalF)) then
      call pwpart(pp,file_hdf5,BptPointer(2,:)*normalF(2),Nplocal,ierr)
    else
      call pwpart(pp,file_hdf5,BptPointer(2,:),Nplocal,ierr)
    endif
    
    ! Initialize the file data for particle postion y
    call file_hdf5%new(&
                &unitDimension=(/1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/),&
                &records='position',&
                &component='y')

    ! Write position y
    if(present(normalF)) then
      call pwpart(pp,file_hdf5,BptPointer(3,:)*normalF(3),Nplocal,ierr)
    else
      call pwpart(pp,file_hdf5,BptPointer(3,:),Nplocal,ierr)
    endif
    
     ! Initialize the file data for particle momentum y
    call file_hdf5%new(&
                &unitDimension=(/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/),&
                &records='moments',&
                &component='y')

    ! Write momentum y
    if(present(normalF)) then
      call pwpart(pp,file_hdf5,BptPointer(4,:)*normalF(4),Nplocal,ierr)
    else
      call pwpart(pp,file_hdf5,BptPointer(4,:),Nplocal,ierr)
    endif
    
    ! Initialize the file data for particle postion z
    call file_hdf5%new(&
                &unitDimension=(/1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/),&
                &records='position',&
                &component='z')

    ! Write position t (or position z)
    if(present(normalF)) then
      call pwpart(pp,file_hdf5,BptPointer(5,:)*normalF(5),Nplocal,ierr)
    else
      call pwpart(pp,file_hdf5,BptPointer(5,:),Nplocal,ierr)
    endif

     ! Initialize the file data for particle momentum z
    call file_hdf5%new(&
                &unitDimension=(/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/),&
                &records='moments',&
                &component='z')

    ! Write momentum pt (or momentum pz)
    if(present(normalF)) then
      call pwpart(pp,file_hdf5,BptPointer(6,:)*normalF(6),Nplocal,ierr)
    else
      call pwpart(pp,file_hdf5,BptPointer(6,:),Nplocal,ierr)
    endif

    ! Initialize the file data for particle mass
    call file_hdf5%new(&
                &unitDimension=(/0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0/),&
                &records='mass',&
                &component='')
    call pwpart(pp,file_hdf5,mass,ierr)
  end subroutine hdf5_particle_output

  subroutine hdf5_field3d_output(rho,nx,ny,nz,iteration)
    double precision, dimension(:,:,:), intent(in) :: rho
    integer, intent(in) :: nx,ny,nz,iteration
    integer :: ierr

    ! Initialize the file data for 3D mesh rho3d
    call file_hdf5%new(iter=iteration,&
                &axisLabels=(/'x','y','z'/),&
                &gridSpacing=(/1.0,1.0,1.0/),&
                &gridGlobalOffset=(/0.0d0,0.0d0,0.0d0/),&
                &position=(/0.0,0.0,0.0/),&
                &unitDimension=(/-3.d0,0.d0,1.d0,1.d0,0.d0,0.d0,0.d0/),&
                &records='field1')

    ! Write rho3d
    !print*,nx,ny,nz,HDF5nyp,HDF5noff
    !if (mod(ny,pp%getnvp()) /= 0) then
    !   print *, 'The number of processor cannot divde nz=', ny
    !   call p%del()
    !   call exit
    !else
    !   HDF5nyp = ny/pp%getnvp()
    !end if

    HDF5noff = pp%getidproc()*ny

    call pwfield(pp,file_hdf5,rho(:,:,:),(/nx,ny*pp%getnvp(),nz/),(/nx,ny,nz/),&
    		&(/0,HDF5noff,0/),ierr)

  end subroutine hdf5_field3d_output

end module hdf5_interface_class

