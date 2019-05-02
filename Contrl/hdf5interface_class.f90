! Interface bewteen IMPACTZ and hdf5io module for OpenPMD
! modification: 4/26/2019 Kilean
!               - particle data 
!                 - openPMD format
!                 - simple array format for parallel read-in from Impact
!               - TBT

module hdf5_interface_class
  use hdf5io_class

  implicit none

  type(hdf5file), private :: file_hdf5
  integer :: HDF5noff, HDF5nyp, HDF5Out,HDF5FieldOutFlag=0

  contains

  subroutine openPMD_particle_output(BptPointer,Nplocal,nfile,iteration,&
                                 &mass,pName,samplePeriod,normalF)
    double precision, dimension(:,:), intent(in) :: BptPointer
    double precision, intent(in) :: mass !kg
    character(len=*), intent(in) :: pName
    double precision, optional, intent(in) :: normalF(6)
    integer, intent(in) :: Nplocal,nfile,iteration,samplePeriod
    integer :: ierr

!    ! Initialize the file data for particle postion x
    call file_hdf5%new(iter=iteration,&
                      &particleName=trim(pName),&
                      &filename = 'openPMD.'//trim(str(nfile))//'.h5',&
                      &unitDimension=(/1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/),&
                      &records='position',&
                      &component='x')
    if(samplePeriod>1) then
      ! Write position x
      if(present(normalF)) then
        call pwpart(file_hdf5,BptPointer(1,1:Nplocal:samplePeriod)*normalF(1),Nplocal/samplePeriod,ierr)
      else
        call pwpart(file_hdf5,BptPointer(1,1:Nplocal:samplePeriod),Nplocal/samplePeriod,ierr)
      endif
      
       ! Initialize the file data for particle momentum x
      call file_hdf5%new(&
                  &unitDimension=(/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/),&
                  &records='moments',&
                  &component='x')
      !print*,BptPointer(2,1:6)

      ! Write momentum x
      if(present(normalF)) then
        call pwpart(file_hdf5,BptPointer(2,1:Nplocal:samplePeriod)*normalF(2),Nplocal/samplePeriod,ierr)
      else
        call pwpart(file_hdf5,BptPointer(2,1:Nplocal:samplePeriod),Nplocal/samplePeriod,ierr)
      endif
      
      ! Initialize the file data for particle postion y
      call file_hdf5%new(&
                  &unitDimension=(/1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/),&
                  &records='position',&
                  &component='y')

      ! Write position y
      if(present(normalF)) then
        call pwpart(file_hdf5,BptPointer(3,1:Nplocal:samplePeriod)*normalF(3),Nplocal/samplePeriod,ierr)
      else
        call pwpart(file_hdf5,BptPointer(3,1:Nplocal:samplePeriod),Nplocal/samplePeriod,ierr)
      endif
      
       ! Initialize the file data for particle momentum y
      call file_hdf5%new(&
                  &unitDimension=(/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/),&
                  &records='moments',&
                  &component='y')

      ! Write momentum y
      if(present(normalF)) then
        call pwpart(file_hdf5,BptPointer(4,1:Nplocal:samplePeriod)*normalF(4),Nplocal/samplePeriod,ierr)
      else
        call pwpart(file_hdf5,BptPointer(4,1:Nplocal:samplePeriod),Nplocal/samplePeriod,ierr)
      endif
      
      ! Initialize the file data for particle postion z
      call file_hdf5%new(&
                  &unitDimension=(/1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/),&
                  &records='position',&
                  &component='z')

      ! Write position t (or position z)
      if(present(normalF)) then
        call pwpart(file_hdf5,BptPointer(5,1:Nplocal:samplePeriod)*normalF(5),Nplocal/samplePeriod,ierr)
      else
        call pwpart(file_hdf5,BptPointer(5,1:Nplocal:samplePeriod),Nplocal/samplePeriod,ierr)
      endif

       ! Initialize the file data for particle momentum z
      call file_hdf5%new(&
                  &unitDimension=(/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/),&
                  &records='moments',&
                  &component='z')

      ! Write momentum pt (or momentum pz)
      if(present(normalF)) then
        call pwpart(file_hdf5,BptPointer(6,1:Nplocal:samplePeriod)*normalF(6),Nplocal/samplePeriod,ierr)
      else
        call pwpart(file_hdf5,BptPointer(6,1:Nplocal:samplePeriod),Nplocal/samplePeriod,ierr)
      endif
    else
      ! Write position x
      if(present(normalF)) then
        call pwpart(file_hdf5,BptPointer(1,1:Nplocal)*normalF(1),Nplocal,ierr)
      else
        call pwpart(file_hdf5,BptPointer(1,1:Nplocal),Nplocal,ierr)
      endif
      
       ! Initialize the file data for particle momentum x
      call file_hdf5%new(&
                  &unitDimension=(/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/),&
                  &records='moments',&
                  &component='x')
      !print*,BptPointer(2,1:6)

      ! Write momentum x
      if(present(normalF)) then
        call pwpart(file_hdf5,BptPointer(2,1:Nplocal)*normalF(2),Nplocal,ierr)
      else
        call pwpart(file_hdf5,BptPointer(2,1:Nplocal),Nplocal,ierr)
      endif
      
      ! Initialize the file data for particle postion y
      call file_hdf5%new(&
                  &unitDimension=(/1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/),&
                  &records='position',&
                  &component='y')

      ! Write position y
      if(present(normalF)) then
        call pwpart(file_hdf5,BptPointer(3,1:Nplocal)*normalF(3),Nplocal,ierr)
      else
        call pwpart(file_hdf5,BptPointer(3,1:Nplocal),Nplocal,ierr)
      endif
      
       ! Initialize the file data for particle momentum y
      call file_hdf5%new(&
                  &unitDimension=(/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/),&
                  &records='moments',&
                  &component='y')

      ! Write momentum y
      if(present(normalF)) then
        call pwpart(file_hdf5,BptPointer(4,1:Nplocal)*normalF(4),Nplocal,ierr)
      else
        call pwpart(file_hdf5,BptPointer(4,1:Nplocal),Nplocal,ierr)
      endif
      
      ! Initialize the file data for particle postion z
      call file_hdf5%new(&
                  &unitDimension=(/1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/),&
                  &records='position',&
                  &component='z')

      ! Write position t (or position z)
      if(present(normalF)) then
        call pwpart(file_hdf5,BptPointer(5,1:Nplocal)*normalF(5),Nplocal,ierr)
      else
        call pwpart(file_hdf5,BptPointer(5,1:Nplocal),Nplocal,ierr)
      endif

       ! Initialize the file data for particle momentum z
      call file_hdf5%new(&
                  &unitDimension=(/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/),&
                  &records='moments',&
                  &component='z')

      ! Write momentum pt (or momentum pz)
      if(present(normalF)) then
        call pwpart(file_hdf5,BptPointer(6,1:Nplocal)*normalF(6),Nplocal,ierr)
      else
        call pwpart(file_hdf5,BptPointer(6,1:Nplocal),Nplocal,ierr)
      endif
    endif

    ! Initialize the file data for particle mass
    call file_hdf5%new(&
                &unitDimension=(/0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0/),&
                &records='mass',&
                &component='')
    call pwpart(file_hdf5,mass,ierr)
  end subroutine openPMD_particle_output
  
  
  subroutine hdf5_particle_output(nfile,pData,npt)
!=======================================================
!  assumes 'mat' is column split among mpi tasks
!=======================================================
  implicit none
  integer, intent(in) :: nfile,npt
  double precision, intent(in) :: pData(9,npt)

    call hdf5_write_matrix_double('partcl.'//trim(str(nfile))//'.h5',&
                             &pData,9,npt)

  end subroutine hdf5_particle_output



end module hdf5_interface_class

