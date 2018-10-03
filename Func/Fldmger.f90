!----------------------------------------------------------------
! (c) Copyright, 2001 by the Regents of the University of California.
! Fldmgerclass: Exchange guard grid data between neighboring processors 
! class in Communication module of FUNCTION layer.
! Version: 1.0
! Author: Ji Qiang, LANL, 7/13/01
! Description: This class defines the functions to sum up the particle
!              contribution from neighboring processor domain, exchange
!              the potential, exchange the field for interpolation between
!              neighboring processors.
! Comments:
!----------------------------------------------------------------
        module Fldmgerclass
          use Timerclass
          use Pgrid2dclass
        contains
!----------------------------------------------------------------
! neighboring grid communication for the 3D open boundary conditions
        ! sum up the contributions from neighboring guard cells.
        subroutine guardsum1_Fldmger(x,innx,inny,innz,grid) 
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx, inny, innz
        double precision, dimension(innx,inny,innz), intent(inout) :: x
        type (Pgrid2d), intent(in) :: grid
        double precision, dimension(innx,innz) :: sendbuf1, recvbuf1
        double precision, dimension(innx,inny) :: sendbuf2, recvbuf2
        double precision, dimension(innz)  :: tmp1,sendtmp
        integer :: left,right,bottom,top,msid,&
                   ierr,i,j,k,nx,ny
        integer :: myid,myidx,myidy,comm2d,commcol,commrow,totnp,npy,&
                   npx
        integer, dimension(2) :: tmpcoord
        integer status(MPI_STATUS_SIZE)
        double precision :: t0
        integer :: nsend,jadd,kadd
    
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)
        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        if(myidx.ne.(npx-1)) then
          right = myidx + 1
        else
          right = MPI_PROC_NULL
        endif
        if(myidx.ne.0) then
          left = myidx - 1
        else
          left = MPI_PROC_NULL
        endif 
        if(npx.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif

        if(myidy.ne.npy-1) then
          top = myidy + 1
        else
          top = MPI_PROC_NULL
        endif
        if(myidy.ne.0) then
          bottom = myidy -1
        else
          bottom = MPI_PROC_NULL
        endif
        if(npy.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif

        if(npy.gt.1) then

        nsend = innz*innx
          call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,1,commcol,&
                         msid,ierr)
          do k = 1, innz
            do i = 1, innx
              sendbuf1(i,k) = x(i,inny,k)
            enddo
          enddo
          call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,1,commcol,&
                        ierr)
          call MPI_WAIT(msid,status,ierr)
        if(myidy.ne.0) then
          do k = 1, innz
            do i = 1, innx
              x(i,2,k) = x(i,2,k) + recvbuf1(i,k)
            enddo
          enddo
        endif
          
          call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,2,commcol,&
                         msid,ierr)
          do k = 1, innz
            do i = 1, innx
              sendbuf1(i,k) = x(i,1,k)
            enddo
          enddo
          call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,2,commcol,&
                        ierr)
          call MPI_WAIT(msid,status,ierr)
        if(myidy.ne.(npy-1)) then
          do k = 1, innz
            do i = 1, innx
              x(i,inny-1,k) = x(i,inny-1,k) + recvbuf1(i,k)
            enddo
          enddo
        endif

        endif

        if(npx.gt.1) then
          nsend = innx*(inny-2*jadd)
          call MPI_IRECV(recvbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,left,3,commrow,&
                         msid,ierr)
          do j = 1+jadd, inny-jadd
            do i = 1, innx
              sendbuf2(i,j-jadd) = x(i,j,innz)
            enddo
          enddo
          call MPI_SEND(sendbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,right,3,commrow,&
                        ierr)
          call MPI_WAIT(msid,status,ierr)
          if(myidx.ne.0) then
          do j = 1+jadd, inny-jadd
            do i = 1, innx
              x(i,j,2) = x(i,j,2) + recvbuf2(i,j-jadd)
            enddo
          enddo
          endif

          call MPI_IRECV(recvbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,right,4,commrow,&
                         msid,ierr)
          do j = 1+jadd, inny-jadd
            do i = 1, innx
              sendbuf2(i,j-jadd) = x(i,j,1)
            enddo
          enddo
          call MPI_SEND(sendbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,left,4,commrow,&
                        ierr)
          call MPI_WAIT(msid,status,ierr)
          if(myidx.ne.(npx-1)) then
          do j = 1+jadd, inny-jadd
            do i = 1, innx
              x(i,j,innz-1) = x(i,j,innz-1) + recvbuf2(i,j-jadd)
            enddo
          enddo
          endif
        endif

        call MPI_BARRIER(comm2d,ierr)
        t_guardsum = t_guardsum + elapsedtime_Timer(t0)

        end subroutine guardsum1_Fldmger

        ! exchange grid information between neighboring guard cell
        ! to calculate E from phi. 
        subroutine guardexch1_Fldmger(x,innx,inny,innz,grid)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx, inny, innz
        double precision, dimension(innx,inny,innz), intent(inout) :: x
        type (Pgrid2d), intent(in) :: grid
        double precision, dimension(innx,innz) :: sendbuf1, recvbuf1
        double precision, dimension(innx,inny) :: sendbuf2, recvbuf2
        integer :: left,right,bottom,top
        integer :: leftx,rightx,bottomy,topy
        integer :: myid,myidx,myidy,comm2d,commcol,commrow,totnp,npx,&
                   npy,msid
        integer, dimension(2) :: tmpcoord
        integer :: i,j,k,ierr
        integer status(MPI_STATUS_SIZE)
        double precision :: t0
        integer :: nsend,jadd,kadd
        
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        !define left,right,top,and bottom for cyclic shift purpose.
        if(myidx.eq.0) then
          left = npx - 1
          right = myidx + 1
        else if(myidx.eq.(npx-1)) then
          left = myidx - 1
          right = 0
        else
          left = myidx - 1
          right = myidx + 1
        endif
        if(npx.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif

        if(myidy.eq.0) then
          bottom = npy - 1
          top = myidy + 1
        else if(myidy.eq.(npy-1)) then
          bottom = myidy - 1
          top = 0
        else
          bottom = myidy - 1
          top = myidy + 1
        endif
        if(npy.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif

        if(npy.gt.1) then

        do k = 1+kadd, innz-kadd
          do i = 1, innx
            sendbuf1(i,k-kadd) = x(i,inny-1,k)
          enddo
        enddo
        nsend = (innz-2*kadd)*innx
        call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,1,&
                      commcol,msid,ierr)
        call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,1,&
                      commcol,ierr)
        call MPI_WAIT(msid,status,ierr)
        do k = 1+kadd, innz-kadd
          do i = 1, innx
            x(i,1,k) = recvbuf1(i,k-kadd)
          enddo
        enddo
          
        do k = 1+kadd, innz-kadd
          do i = 1, innx
            sendbuf1(i,k-kadd) = x(i,2,k)
          enddo
        enddo
        call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,2,&
                      commcol,msid,ierr)
        call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,2,&
                      commcol,ierr)
        call MPI_WAIT(msid,status,ierr)
        do k = 1+kadd, innz-kadd
          do i = 1, innx
            x(i,inny,k) = recvbuf1(i,k-kadd)
          enddo
        enddo

        endif

        if(npx.gt.1) then

        nsend = innx*(inny-2*jadd)
        do k = 1+jadd, inny-jadd
          do j = 1, innx
            sendbuf2(j,k-jadd) = x(j,k,innz-1)
          enddo
        enddo
        call MPI_IRECV(recvbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,left,3,&
                      commrow,msid,ierr)
        call MPI_SEND(sendbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,right,3,&
                      commrow,ierr)
        call MPI_WAIT(msid,status,ierr)
        do k = 1+jadd, inny-jadd
          do j = 1, innx
            x(j,k,1) = recvbuf2(j,k-jadd)
          enddo
        enddo
          
        do k = 1+jadd, inny-jadd
          do j = 1, innx
            sendbuf2(j,k-jadd) = x(j,k,2)
          enddo
        enddo
        call MPI_IRECV(recvbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,right,4,&
                      commrow,msid,ierr)
        call MPI_SEND(sendbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,left,4,&
                      commrow,ierr)
        call MPI_WAIT(msid,status,ierr)
        do k = 1+jadd, inny-jadd
          do j = 1, innx
            x(j,k,innz) = recvbuf2(j,k-jadd)
          enddo
        enddo
 
        endif

        call MPI_BARRIER(comm2d,ierr)
        t_guardexch = t_guardexch + elapsedtime_Timer(t0)

        end subroutine guardexch1_Fldmger

!----------------------------------------------------------------
! neighboring grid communication for the 2D open, longitudinal periodic
! boundary conditions
        ! sum up the contributions from neighboring guard cells.
        subroutine guardsum2_Fldmger(x,innx,inny,innz,grid) 
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx, inny, innz
        double precision, dimension(innx,inny,innz), intent(inout) :: x
        type (Pgrid2d), intent(in) :: grid
        double precision, dimension(innx,innz) :: sendbuf1, recvbuf1
        double precision, dimension(innx,inny) :: sendbuf2, recvbuf2
        double precision, dimension(innz)  :: tmp1,sendtmp
        integer :: left,right,bottom,top,msid,&
                   ierr,i,j,k,nx,ny
        integer :: myid,myidx,myidy,comm2d,commcol,commrow,totnp,npy,&
                   npx
        integer, dimension(2) :: tmpcoord
        integer status(MPI_STATUS_SIZE)
        double precision :: t0,tmp
        integer :: nsend,jadd,kadd
    
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)
        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        if(myidx.ne.(npx-1)) then
          right = myidx + 1
        else
          right = 0
        endif
        if(myidx.ne.0) then
          left = myidx - 1
        else
          left = npx - 1
        endif 
        if(npx.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif

        if(myidy.ne.npy-1) then
          top = myidy + 1
        else
          top = MPI_PROC_NULL
        endif
        if(myidy.ne.0) then
          bottom = myidy -1
        else
          bottom = MPI_PROC_NULL
        endif
        if(npy.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif

        if(npy.gt.1) then

        nsend = innz*innx
          call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,1,commcol,&
                         msid,ierr)
          do k = 1, innz
            do i = 1, innx
              sendbuf1(i,k) = x(i,inny,k)
            enddo
          enddo
          call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,1,commcol,&
                        ierr)
          call MPI_WAIT(msid,status,ierr)
        if(myidy.ne.0) then
          do k = 1, innz
            do i = 1, innx
              x(i,2,k) = x(i,2,k) + recvbuf1(i,k)
            enddo
          enddo
        endif
          
          call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,2,commcol,&
                         msid,ierr)
          do k = 1, innz
            do i = 1, innx
              sendbuf1(i,k) = x(i,1,k)
            enddo
          enddo
          call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,2,commcol,&
                        ierr)
          call MPI_WAIT(msid,status,ierr)
        if(myidy.ne.(npy-1)) then
          do k = 1, innz
            do i = 1, innx
              x(i,inny-1,k) = x(i,inny-1,k) + recvbuf1(i,k)
            enddo
          enddo
        endif

        endif

        if(npx.gt.1) then
          if(myidx.eq.(npx-1)) then
            do j = 1+jadd, inny-jadd
              do i = 1, innx
                x(i,j,innz) = x(i,j,innz-1)
              enddo
            enddo
          endif
          if(myidx.eq.0) then
            do j = 1+jadd, inny-jadd
              do i = 1, innx
                x(i,j,1) = x(i,j,2)
              enddo
            enddo
          endif
          nsend = innx*(inny-2*jadd)
          call MPI_IRECV(recvbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,left,3,commrow,&
                         msid,ierr)
          do j = 1+jadd, inny-jadd
            do i = 1, innx
              sendbuf2(i,j-jadd) = x(i,j,innz)
            enddo
          enddo
          call MPI_SEND(sendbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,right,3,commrow,&
                        ierr)
          call MPI_WAIT(msid,status,ierr)
          do j = 1+jadd, inny-jadd
            do i = 1, innx
              x(i,j,2) = x(i,j,2) + recvbuf2(i,j-jadd)
            enddo
          enddo

          call MPI_IRECV(recvbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,right,4,commrow,&
                         msid,ierr)
          do j = 1+jadd, inny-jadd
            do i = 1, innx
              sendbuf2(i,j-jadd) = x(i,j,1)
            enddo
          enddo
          call MPI_SEND(sendbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,left,4,commrow,&
                        ierr)
          call MPI_WAIT(msid,status,ierr)
          do j = 1+jadd, inny-jadd
            do i = 1, innx
              x(i,j,innz-1) = x(i,j,innz-1) + recvbuf2(i,j-jadd)
            enddo
          enddo
        else
          do j = 1+jadd, inny-jadd
            do i = 1, innx
              tmp = x(i,j,innz)
              x(i,j,innz) = x(i,j,innz) + x(i,j,1)
              x(i,j,1) = x(i,j,1) + tmp
            enddo
          enddo
        endif

        call MPI_BARRIER(comm2d,ierr)
        t_guardsum = t_guardsum + elapsedtime_Timer(t0)

        end subroutine guardsum2_Fldmger

        ! exchange grid information between neighboring guard cell
        ! to calculate E from phi. 
        subroutine guardexch2_Fldmger(x,innx,inny,innz,grid)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx, inny, innz
        double precision, dimension(innx,inny,innz), intent(inout) :: x
        type (Pgrid2d), intent(in) :: grid
        double precision, dimension(innx,innz) :: sendbuf1, recvbuf1
        double precision, dimension(innx,inny) :: sendbuf2, recvbuf2
        integer :: left,right,bottom,top
        integer :: leftx,rightx,bottomy,topy
        integer :: myid,myidx,myidy,comm2d,commcol,commrow,totnp,npx,&
                   npy,msid
        integer, dimension(2) :: tmpcoord
        integer :: i,j,k,ierr
        integer status(MPI_STATUS_SIZE)
        double precision :: t0
        integer :: nsend,jadd,kadd
        
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        !define left,right,top,and bottom for cyclic shift purpose.
        if(myidx.eq.0) then
          left = npx - 1
          right = myidx + 1
        else if(myidx.eq.(npx-1)) then
          left = myidx - 1
          right = 0
        else
          left = myidx - 1
          right = myidx + 1
        endif
        if(npx.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif

        if(myidy.eq.0) then
          bottom = npy - 1
          top = myidy + 1
        else if(myidy.eq.(npy-1)) then
          bottom = myidy - 1
          top = 0
        else
          bottom = myidy - 1
          top = myidy + 1
        endif
        if(npy.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif

        if(npy.gt.1) then

        do k = 1+kadd, innz-kadd
          do i = 1, innx
            sendbuf1(i,k-kadd) = x(i,inny-1,k)
          enddo
        enddo
        nsend = (innz-2*kadd)*innx
        call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,1,&
                      commcol,msid,ierr)
        call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,1,&
                      commcol,ierr)
        call MPI_WAIT(msid,status,ierr)
        do k = 1+kadd, innz-kadd
          do i = 1, innx
            x(i,1,k) = recvbuf1(i,k-kadd)
          enddo
        enddo
          
        do k = 1+kadd, innz-kadd
          do i = 1, innx
            sendbuf1(i,k-kadd) = x(i,2,k)
          enddo
        enddo
        call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,2,&
                      commcol,msid,ierr)
        call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,2,&
                      commcol,ierr)
        call MPI_WAIT(msid,status,ierr)
        do k = 1+kadd, innz-kadd
          do i = 1, innx
            x(i,inny,k) = recvbuf1(i,k-kadd)
          enddo
        enddo

        endif

        if(npx.gt.1) then

        nsend = innx*(inny-2*jadd)
        !apply longitudinal periodic BC for potential.
        if(myidx.ne.(npx-1)) then
          do k = 1+jadd, inny-jadd
            do j = 1, innx
              sendbuf2(j,k-jadd) = x(j,k,innz-1)
            enddo
          enddo
        else
          do k = 1+jadd, inny-jadd
            do j = 1, innx
              sendbuf2(j,k-jadd) = x(j,k,innz-2)
            enddo
          enddo
        endif
        call MPI_IRECV(recvbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,left,3,&
                      commrow,msid,ierr)
        call MPI_SEND(sendbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,right,3,&
                      commrow,ierr)
        call MPI_WAIT(msid,status,ierr)
        do k = 1+jadd, inny-jadd
          do j = 1, innx
              x(j,k,1) = recvbuf2(j,k-jadd)
          enddo
        enddo
          
        !apply longitudinal periodic BC for potential.
        if(myidx.ne.0) then
          do k = 1+jadd, inny-jadd
            do j = 1, innx
              sendbuf2(j,k-jadd) = x(j,k,2)
            enddo
          enddo
        else
          do k = 1+jadd, inny-jadd
            do j = 1, innx
              sendbuf2(j,k-jadd) = x(j,k,3)
            enddo
          enddo
        endif
        call MPI_IRECV(recvbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,right,4,&
                      commrow,msid,ierr)
        call MPI_SEND(sendbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,left,4,&
                      commrow,ierr)
        call MPI_WAIT(msid,status,ierr)
        do k = 1+jadd, inny-jadd
          do j = 1, innx
            x(j,k,innz) = recvbuf2(j,k-jadd)
          enddo
        enddo

        !apply longitudinal periodic BC for potential.
        nsend = inny*innx
        if(myidx.eq.0) then
          do k = 1, inny
            do j = 1, innx
              sendbuf2(j,k) = x(j,k,2)
            enddo
          enddo
          call MPI_SEND(sendbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,left,4,&
                        commrow,ierr)
        else if(myidx.eq.(npx-1)) then
          call MPI_IRECV(recvbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,right,4,&
                         commrow,msid,ierr)
          call MPI_WAIT(msid,status,ierr)
          do k = 1, inny
            do j = 1, innx
              x(j,k,innz-1) = recvbuf2(j,k)
            enddo
          enddo
        endif

        else

          do k = 1, inny
            do j = 1, innx
              x(j,k,innz) = x(j,k,1)
            enddo
          enddo
 
        endif

        call MPI_BARRIER(comm2d,ierr)
        t_guardexch = t_guardexch + elapsedtime_Timer(t0)

        end subroutine guardexch2_Fldmger
!----------------------------------------------------------------
! neighboring communcition for the round pipe in cylindrincal coord. with
! longitudinal open BC.
        subroutine guardsum3_Fldmger(x,innx,inny,innz,grid) 
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx, inny, innz
        double precision, dimension(innx,inny,innz), intent(inout) :: x
        type (Pgrid2d), intent(in) :: grid
        double precision, dimension(innx,innz) :: sendbuf1, recvbuf1
        !double precision, dimension(innx,inny-2) :: sendbuf2, recvbuf2
        double precision, allocatable, dimension(:,:) :: sendbuf2, recvbuf2
        double precision, dimension(innz)  :: tmp1,sendtmp
        integer :: left,right,bottom,top,msid,ierr,i,j,k,nx,ny
        integer :: myid,myidx,myidy,comm2d,commcol,commrow,totnp,npy,&
                   npx,jadd,kadd
        integer, dimension(2) :: tmpcoord
        integer status(MPI_STATUS_SIZE)
        double precision :: t0,tmp
        integer :: nsend
    
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)
        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        if(npy.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif
        if(npx.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif

        if(myidx.ne.(npx-1)) then
          right = myidx + 1
        else
          right = MPI_PROC_NULL
        endif
        if(myidx.ne.0) then
          left = myidx - 1
        else
          left = MPI_PROC_NULL
        endif 

        if(myidy.ne.npy-1) then
          top = myidy + 1
        else
          top = 0
        endif
        if(myidy.ne.0) then
          bottom = myidy - 1
        else
          bottom = npy - 1
        endif

        if(npy.gt.1) then
          if(myidy.eq.(npy-1)) then
            do k = 1, innz
              do i = 1, innx
                x(i,inny,k) = x(i,inny-1,k)
              enddo
            enddo
          endif
          if(myidy.eq.0) then
            do k = 1, innz
              do i = 1, innx
                x(i,1,k) = x(i,2,k)
              enddo
            enddo
          endif
          nsend = innz*innx
          call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,1,commcol,&
                         msid,ierr)
          do k = 1, innz
            do i = 1, innx
              sendbuf1(i,k) = x(i,inny,k)
            enddo
          enddo
          call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,1,commcol,&
                        ierr)
          call MPI_WAIT(msid,status,ierr)
          do k = 1, innz
            do i = 1, innx
              x(i,2,k) = x(i,2,k) + recvbuf1(i,k)
            enddo
          enddo
          
          call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,2,commcol,&
                         msid,ierr)
          do k = 1, innz
            do i = 1, innx
              sendbuf1(i,k) = x(i,1,k)
            enddo
          enddo
          call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,2,commcol,&
                        ierr)
          call MPI_WAIT(msid,status,ierr)
          do k = 1, innz
            do i = 1, innx
              x(i,inny-1,k) = x(i,inny-1,k) + recvbuf1(i,k)
            enddo
          enddo
        else
          do k = 1, innz
            do i = 1, innx
              tmp = x(i,inny,k)
              x(i,inny,k) = x(i,inny,k) + x(i,1,k)
              x(i,1,k) = x(i,1,k) + tmp
            enddo
          enddo
        endif

        if(npx.gt.1) then
          allocate(sendbuf2(innx,inny-2*jadd))
          allocate(recvbuf2(innx,inny-2*jadd))

          nsend = innx*(inny-2*jadd)
          call MPI_IRECV(recvbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,left,3,commrow,&
                         msid,ierr)
          do j = 1+jadd, inny-jadd
            do i = 1, innx
              sendbuf2(i,j-jadd) = x(i,j,innz)
            enddo
          enddo
          call MPI_SEND(sendbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,right,3,commrow,&
                        ierr)
          call MPI_WAIT(msid,status,ierr)
          if(myidx.ne.0) then
          do j = 1+jadd, inny-jadd
            do i = 1, innx
              x(i,j,2) = x(i,j,2) + recvbuf2(i,j-jadd)
            enddo
          enddo
          endif

          call MPI_IRECV(recvbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,right,4,commrow,&
                         msid,ierr)
          do j = 1+jadd, inny-jadd
            do i = 1, innx
              sendbuf2(i,j-jadd) = x(i,j,1)
            enddo
          enddo
          call MPI_SEND(sendbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,left,4,commrow,&
                        ierr)
          call MPI_WAIT(msid,status,ierr)
          if(myidx.ne.(npx-1)) then
          do j = 1+jadd, inny-jadd
            do i = 1, innx
              x(i,j,innz-1) = x(i,j,innz-1) + recvbuf2(i,j-jadd)
            enddo
          enddo
          endif
          
          deallocate(sendbuf2)
          deallocate(recvbuf2)
        endif

        call MPI_BARRIER(comm2d,ierr)
        t_guardsum = t_guardsum + elapsedtime_Timer(t0)

        end subroutine guardsum3_Fldmger

        subroutine guardexch3_Fldmger(x,innx,inny,innz,grid)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx, inny, innz
        double precision, dimension(innx,inny,innz), intent(inout) :: x
        type (Pgrid2d), intent(in) :: grid
        double precision, dimension(innx,innz) :: sendbuf1, recvbuf1
!        double precision, dimension(innx,inny-2) :: sendbuf2, recvbuf2
!        double precision, allocatable, dimension(:,:) :: sendbuf1, recvbuf1
        double precision, allocatable, dimension(:,:) :: sendbuf2, recvbuf2
        integer :: left,right,bottom,top
        integer :: leftx,rightx,bottomy,topy
        integer :: myid,myidx,myidy,comm2d,commcol,commrow,totnp,npx,&
                   npy,msid
        integer, dimension(2) :: tmpcoord
        integer :: i,j,k,ierr
        integer status(MPI_STATUS_SIZE)
        double precision :: t0
        integer :: nsend,jadd,kadd
        
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        !define left,right,top,and bottom for cyclic shift purpose.
        if(myidx.eq.0) then
          left = npx - 1
          right = myidx + 1
        else if(myidx.eq.(npx-1)) then
          left = myidx - 1
          right = 0
        else
          left = myidx - 1
          right = myidx + 1
        endif

        if(myidy.eq.0) then
          bottom = npy - 1
          top = myidy + 1
        else if(myidy.eq.(npy-1)) then
          bottom = myidy - 1
          top = 0
        else
          bottom = myidy - 1
          top = myidy + 1
        endif
        if(npy.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif
        if(npx.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif

!        if(npy.gt.1) then
!        endif

        if(npx.gt.1) then

        allocate(sendbuf2(innx,inny-2*jadd))
        allocate(recvbuf2(innx,inny-2*jadd))
        nsend = innx*(inny-2*jadd)
          do k = 1+jadd, inny-jadd
            do j = 1, innx
              sendbuf2(j,k-jadd) = x(j,k,innz-1)
            enddo
          enddo
        call MPI_IRECV(recvbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,left,3,&
                      commrow,msid,ierr)
        call MPI_SEND(sendbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,right,3,&
                      commrow,ierr)
        call MPI_WAIT(msid,status,ierr)
        do k = 1+jadd, inny-jadd
          do j = 1, innx
              x(j,k,1) = recvbuf2(j,k-jadd)
          enddo
        enddo
          
        do k = 1+jadd, inny-jadd
           do j = 1, innx
             sendbuf2(j,k-jadd) = x(j,k,2)
           enddo
        enddo
        call MPI_IRECV(recvbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,right,4,&
                      commrow,msid,ierr)
        call MPI_SEND(sendbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,left,4,&
                      commrow,ierr)
        call MPI_WAIT(msid,status,ierr)
        do k = 1+jadd, inny-jadd
          do j = 1, innx
            x(j,k,innz) = recvbuf2(j,k-jadd)
          enddo
        enddo

        deallocate(sendbuf2)
        deallocate(recvbuf2)
 
        endif

        if(npy.gt.1) then

!        allocate(sendbuf1(innx,innz-2*kadd))
!        allocate(recvbuf1(innx,innz-2*kadd))
        ! periodic BC
        if(myidy.ne.(npy-1)) then
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              sendbuf1(i,k-kadd) = x(i,inny-1,k)
            enddo
          enddo
        else
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              sendbuf1(i,k-kadd) = x(i,inny-2,k)
            enddo
          enddo
        endif
        nsend = (innz-2*kadd)*innx
        call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,1,&
                      commcol,msid,ierr)
        call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,1,&
                      commcol,ierr)
        call MPI_WAIT(msid,status,ierr)
        do k = 1+kadd, innz-kadd
          do i = 1, innx
            x(i,1,k) = recvbuf1(i,k-kadd)
          enddo
        enddo
          
        !periodic BC
        if(myidy.ne.0) then
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              sendbuf1(i,k-kadd) = x(i,2,k)
            enddo
          enddo
        else
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              sendbuf1(i,k-kadd) = x(i,3,k)
            enddo
          enddo
        endif
        call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,2,&
                      commcol,msid,ierr)
        call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,2,&
                      commcol,ierr)
        call MPI_WAIT(msid,status,ierr)
        do k = 1+kadd, innz-kadd
          do i = 1, innx
            x(i,inny,k) = recvbuf1(i,k-kadd)
          enddo
        enddo

          nsend = innx*innz
          if(myidy.eq.0) then
            do k = 1, innz
              do i = 1, innx
                sendbuf1(i,k) = x(i,2,k)
              enddo
            enddo
            call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,2,&
                          commcol,ierr)
          else if(myidy.eq.(npy-1)) then
            call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,2,&
                           commcol,msid,ierr)
            call MPI_WAIT(msid,status,ierr)
            do k = 1, innz
              do i = 1, innx
                x(i,inny-1,k) = recvbuf1(i,k)
              enddo
            enddo
          endif

!        deallocate(sendbuf1)
!        deallocate(recvbuf1)
        else
          do k = 1, innz
            do i = 1, innx
              x(i,inny,k) = x(i,1,k)
            enddo
          enddo
        endif

        call MPI_BARRIER(comm2d,ierr)
        t_guardexch = t_guardexch + elapsedtime_Timer(t0)

        end subroutine guardexch3_Fldmger
!----------------------------------------------------------------
! neighboring communcition for the round pipe in cylindrincal coord. with
! longitudinal periodic BC.
        subroutine guardsum4_Fldmger(x,innx,inny,innz,grid) 
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx, inny, innz
        double precision, dimension(innx,inny,innz), intent(inout) :: x
        type (Pgrid2d), intent(in) :: grid
        double precision, dimension(innx,innz) :: sendbuf1, recvbuf1
        !double precision, dimension(innx,inny-2) :: sendbuf2, recvbuf2
        double precision, allocatable, dimension(:,:) :: sendbuf2, recvbuf2
        double precision, dimension(innz)  :: tmp1,sendtmp
        integer :: left,right,bottom,top,msid,ierr,i,j,k,nx,ny
        integer :: myid,myidx,myidy,comm2d,commcol,commrow,totnp,npy,&
                   npx,jadd,kadd
        integer, dimension(2) :: tmpcoord
        integer status(MPI_STATUS_SIZE)
        double precision :: t0,tmp
        integer :: nsend
    
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)
        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        if(npy.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif
        if(npx.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif

        if(myidx.ne.(npx-1)) then
          right = myidx + 1
        else
          right = 0
        endif
        if(myidx.ne.0) then
          left = myidx - 1
        else
          left = npx - 1
        endif 

        if(myidy.ne.npy-1) then
          top = myidy + 1
        else
          top = 0
        endif
        if(myidy.ne.0) then
          bottom = myidy - 1
        else
          bottom = npy - 1
        endif

        if(npy.gt.1) then
          if(myidy.eq.(npy-1)) then
            do k = 1, innz
              do i = 1, innx
                x(i,inny,k) = x(i,inny-1,k)
              enddo
            enddo
          endif
          if(myidy.eq.0) then
            do k = 1, innz
              do i = 1, innx
                x(i,1,k) = x(i,2,k)
              enddo
            enddo
          endif
          nsend = innz*innx
          call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,1,commcol,&
                         msid,ierr)
          do k = 1, innz
            do i = 1, innx
              sendbuf1(i,k) = x(i,inny,k)
            enddo
          enddo
          call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,1,commcol,&
                        ierr)
          call MPI_WAIT(msid,status,ierr)
          do k = 1, innz
            do i = 1, innx
              x(i,2,k) = x(i,2,k) + recvbuf1(i,k)
            enddo
          enddo
          
          call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,2,commcol,&
                         msid,ierr)
          do k = 1, innz
            do i = 1, innx
              sendbuf1(i,k) = x(i,1,k)
            enddo
          enddo
          call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,2,commcol,&
                        ierr)
          call MPI_WAIT(msid,status,ierr)
          do k = 1, innz
            do i = 1, innx
              x(i,inny-1,k) = x(i,inny-1,k) + recvbuf1(i,k)
            enddo
          enddo
        else
          do k = 1, innz
            do i = 1, innx
              tmp = x(i,inny,k)
              x(i,inny,k) = x(i,inny,k) + x(i,1,k)
              x(i,1,k) = x(i,1,k) + tmp
            enddo
          enddo
        endif

        if(npx.gt.1) then
          allocate(sendbuf2(innx,inny-2*jadd))
          allocate(recvbuf2(innx,inny-2*jadd))

          if(myidx.eq.(npx-1)) then
            do j = 1+jadd, inny-jadd
              do i = 1, innx
                x(i,j,innz) = x(i,j,innz-1)
              enddo
            enddo
          endif
          if(myidx.eq.0) then
            do j = 1+jadd, inny-jadd
              do i = 1, innx
                x(i,j,1) = x(i,j,2)
              enddo
            enddo
          endif
          nsend = innx*(inny-2*jadd)
          call MPI_IRECV(recvbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,left,3,commrow,&
                         msid,ierr)
          do j = 1+jadd, inny-jadd
            do i = 1, innx
              sendbuf2(i,j-jadd) = x(i,j,innz)
            enddo
          enddo
          call MPI_SEND(sendbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,right,3,commrow,&
                        ierr)
          call MPI_WAIT(msid,status,ierr)
          do j = 1+jadd, inny-jadd
            do i = 1, innx
              x(i,j,2) = x(i,j,2) + recvbuf2(i,j-jadd)
            enddo
          enddo

          call MPI_IRECV(recvbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,right,4,commrow,&
                         msid,ierr)
          do j = 1+jadd, inny-jadd
            do i = 1, innx
              sendbuf2(i,j-jadd) = x(i,j,1)
            enddo
          enddo
          call MPI_SEND(sendbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,left,4,commrow,&
                        ierr)
          call MPI_WAIT(msid,status,ierr)
          do j = 1+jadd, inny-jadd
            do i = 1, innx
              x(i,j,innz-1) = x(i,j,innz-1) + recvbuf2(i,j-jadd)
            enddo
          enddo
          
          deallocate(sendbuf2)
          deallocate(recvbuf2)
        else
          do j = 1, inny
            do i = 1, innx
              tmp = x(i,j,innz)
              x(i,j,innz) = x(i,j,innz) + x(i,j,1)
              x(i,j,1) = x(i,j,1) + tmp
            enddo
          enddo
        endif

        call MPI_BARRIER(comm2d,ierr)
        t_guardsum = t_guardsum + elapsedtime_Timer(t0)

        end subroutine guardsum4_Fldmger

        subroutine guardexch4_Fldmger(x,innx,inny,innz,grid)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx, inny, innz
        double precision, dimension(innx,inny,innz), intent(inout) :: x
        type (Pgrid2d), intent(in) :: grid
        double precision, dimension(innx,innz) :: sendbuf1, recvbuf1
        double precision, dimension(innx,inny) :: sendbuf2, recvbuf2
!        double precision, allocatable, dimension(:,:) :: sendbuf1, recvbuf1
!        double precision, allocatable, dimension(:,:) :: sendbuf2, recvbuf2
        integer :: left,right,bottom,top
        integer :: leftx,rightx,bottomy,topy
        integer :: myid,myidx,myidy,comm2d,commcol,commrow,totnp,npx,&
                   npy,msid
        integer, dimension(2) :: tmpcoord
        integer :: i,j,k,ierr
        integer status(MPI_STATUS_SIZE)
        double precision :: t0
        integer :: nsend,jadd,kadd
        
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        !define left,right,top,and bottom for cyclic shift purpose.
        if(myidx.eq.0) then
          left = npx - 1
          right = myidx + 1
        else if(myidx.eq.(npx-1)) then
          left = myidx - 1
          right = 0
        else
          left = myidx - 1
          right = myidx + 1
        endif

        if(myidy.eq.0) then
          bottom = npy - 1
          top = myidy + 1
        else if(myidy.eq.(npy-1)) then
          bottom = myidy - 1
          top = 0
        else
          bottom = myidy - 1
          top = myidy + 1
        endif
        if(npy.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif
        if(npx.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif

        if(npy.gt.1) then

!        allocate(sendbuf1(innx,innz-2*kadd))
!        allocate(recvbuf1(innx,innz-2*kadd))
        ! periodic BC
        if(myidy.ne.(npy-1)) then
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              sendbuf1(i,k-kadd) = x(i,inny-1,k)
            enddo
          enddo
        else
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              sendbuf1(i,k-kadd) = x(i,inny-2,k)
            enddo
          enddo
        endif
        nsend = (innz-2*kadd)*innx
        call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,1,&
                      commcol,msid,ierr)
        call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,1,&
                      commcol,ierr)
        call MPI_WAIT(msid,status,ierr)
        do k = 1+kadd, innz-kadd
          do i = 1, innx
            x(i,1,k) = recvbuf1(i,k-kadd)
          enddo
        enddo
          
        !periodic BC
        if(myidy.ne.0) then
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              sendbuf1(i,k-kadd) = x(i,2,k)
            enddo
          enddo
        else
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              sendbuf1(i,k-kadd) = x(i,3,k)
            enddo
          enddo
        endif
        call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,2,&
                      commcol,msid,ierr)
        call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,2,&
                      commcol,ierr)
        call MPI_WAIT(msid,status,ierr)
        do k = 1+kadd, innz-kadd
          do i = 1, innx
            x(i,inny,k) = recvbuf1(i,k-kadd)
          enddo
        enddo

        if(myidy.eq.0) then
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              sendbuf1(i,k-kadd) = x(i,2,k)
            enddo
          enddo
          call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,2,&
                        commcol,ierr)
        else if(myidy.eq.(npy-1)) then
          call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,2,&
                         commcol,msid,ierr)
          call MPI_WAIT(msid,status,ierr)
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              x(i,inny-1,k) = recvbuf1(i,k-kadd)
            enddo
          enddo
        endif
   
!        deallocate(sendbuf1)
!        deallocate(recvbuf1)

        else
          do k = 1+kadd,innz-kadd
            do i = 1, innx 
              x(i,inny,k) = x(i,1,k)
            enddo
          enddo
        endif

        if(npx.gt.1) then

!        allocate(sendbuf2(innx,inny-2*jadd))
!        allocate(recvbuf2(innx,inny-2*jadd))
        nsend = innx*(inny-2*jadd)
        !apply longitudinal periodic BC for potential.
        if(myidx.ne.(npx-1)) then
          do k = 1+jadd, inny-jadd
            do j = 1, innx
              sendbuf2(j,k-jadd) = x(j,k,innz-1)
            enddo
          enddo
        else
          do k = 1+jadd, inny-jadd
            do j = 1, innx
              sendbuf2(j,k-jadd) = x(j,k,innz-2)
            enddo
          enddo
        endif
        call MPI_IRECV(recvbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,left,3,&
                      commrow,msid,ierr)
        call MPI_SEND(sendbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,right,3,&
                      commrow,ierr)
        call MPI_WAIT(msid,status,ierr)
        do k = 1+jadd, inny-jadd
          do j = 1, innx
              x(j,k,1) = recvbuf2(j,k-jadd)
          enddo
        enddo
          
        !apply longitudinal periodic BC for potential.
        if(myidx.ne.0) then
          do k = 1+jadd, inny-jadd
            do j = 1, innx
              sendbuf2(j,k-jadd) = x(j,k,2)
            enddo
          enddo
        else
          do k = 1+jadd, inny-jadd
            do j = 1, innx
              sendbuf2(j,k-jadd) = x(j,k,3)
            enddo
          enddo
        endif
        call MPI_IRECV(recvbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,right,4,&
                      commrow,msid,ierr)
        call MPI_SEND(sendbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,left,4,&
                      commrow,ierr)
        call MPI_WAIT(msid,status,ierr)
        do k = 1+jadd, inny-jadd
          do j = 1, innx
            x(j,k,innz) = recvbuf2(j,k-jadd)
          enddo
        enddo

        !apply longitudinal periodic BC for potential.
        nsend = inny*innx
        if(myidx.eq.0) then
          do k = 1, inny
            do j = 1, innx
              sendbuf2(j,k) = x(j,k,2)
            enddo
          enddo
          call MPI_SEND(sendbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,left,4,&
                        commrow,ierr)
        else if(myidx.eq.(npx-1)) then
          call MPI_IRECV(recvbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,right,4,&
                         commrow,msid,ierr)
          call MPI_WAIT(msid,status,ierr)
          do k = 1, inny
            do j = 1, innx
              x(j,k,innz-1) = recvbuf2(j,k)
            enddo
          enddo
        endif

!        deallocate(sendbuf2)
!        deallocate(recvbuf2)

        else
          do k = 1, inny
            do j = 1, innx
              x(j,k,innz) = x(j,k,1)
            enddo
          enddo
        endif

        call MPI_BARRIER(comm2d,ierr)
        t_guardexch = t_guardexch + elapsedtime_Timer(t0)

        end subroutine guardexch4_Fldmger
!----------------------------------------------------------------
        ! exchange the information for interpolation 
        ! from neighboring guard cells.
        subroutine boundint4_Fldmger(x1,x2,x3,innx,inny,innz,grid) 
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx, inny, innz
        double precision, dimension(innx,inny,innz), intent(inout) :: x1,x2,x3
        type (Pgrid2d), intent(in) :: grid
!        double precision, dimension(innx,innz-2,3) :: sendbuf1, recvbuf1
        double precision, allocatable, dimension(:,:,:) :: sendbuf1, recvbuf1
        double precision, dimension(innx,inny,3) :: sendbuf2, recvbuf2
        integer :: left,right,bottom,top,msid,&
                   ierr,i,j,k,nx,ny
        integer :: myid,myidx,myidy,comm2d,commcol,commrow,totnp,&
                   npx,npy
        integer status(MPI_STATUS_SIZE)
        double precision :: t0
        integer :: nsend,jadd,kadd
    
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)
        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)
        if(npy.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif
        if(npx.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif

        top = myidy + 1
        bottom = myidy - 1
        left = myidx - 1
        right = myidx + 1

        ny = inny - 2*jadd

! This could be modified to improve speed
        if(npy.gt.1) then

        nsend = innx*(innz-2*kadd)*3
        allocate(recvbuf1(innx,innz-2*kadd,3))
        allocate(sendbuf1(innx,innz-2*kadd,3))

        if(myidy.ne.0) then
          call MPI_IRECV(recvbuf1(1,1,1),nsend,MPI_DOUBLE_PRECISION,bottom,1,commcol,&
                        msid,ierr)
        endif
        if(myidy.ne.(npy-1)) then
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              sendbuf1(i,k-kadd,1) = x1(i,inny-1,k)
            enddo
          enddo
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              sendbuf1(i,k-kadd,2) = x2(i,inny-1,k)
            enddo
          enddo
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              sendbuf1(i,k-kadd,3) = x3(i,inny-1,k)
            enddo
          enddo
          call MPI_SEND(sendbuf1(1,1,1),nsend,MPI_DOUBLE_PRECISION,top,1,&
                        commcol,ierr)
        endif
        if(myidy.ne.0) then
          call MPI_WAIT(msid,status,ierr)
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              x1(i,1,k) = recvbuf1(i,k-kadd,1)
            enddo
          enddo
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              x2(i,1,k) = recvbuf1(i,k-kadd,2)
            enddo
          enddo
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              x3(i,1,k) = recvbuf1(i,k-kadd,3)
            enddo
          enddo
        endif

        if(myidy.ne.(npy-1)) then
          call MPI_IRECV(recvbuf1(1,1,1),nsend,MPI_DOUBLE_PRECISION,top,2,commcol,&
                        msid,ierr)
        endif
        if(myidy.ne.0) then
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              sendbuf1(i,k-kadd,1) = x1(i,2,k)
            enddo
          enddo
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              sendbuf1(i,k-kadd,2) = x2(i,2,k)
            enddo
          enddo
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              sendbuf1(i,k-kadd,3) = x3(i,2,k)
            enddo
          enddo
          call MPI_SEND(sendbuf1(1,1,1),nsend,MPI_DOUBLE_PRECISION,bottom,2,commcol,&
                        ierr)
        endif
        if(myidy.ne.(npy-1)) then
          call MPI_WAIT(msid,status,ierr)
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              x1(i,inny,k) = recvbuf1(i,k-kadd,1)
            enddo
          enddo
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              x2(i,inny,k) = recvbuf1(i,k-kadd,2)
            enddo
          enddo
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              x3(i,inny,k) = recvbuf1(i,k-kadd,3)
            enddo
          enddo
        endif

        deallocate(recvbuf1)
        deallocate(sendbuf1)

        endif

        if(npx.gt.1) then

        nsend = innx*inny*3
        if(myidx.ne.0) then
          call MPI_IRECV(recvbuf2(1,1,1),nsend,MPI_DOUBLE_PRECISION,left,3,commrow,&
                         msid,ierr)
        endif
        if(myidx.ne.(npx-1)) then
          do j = 1, inny
            do k = 1, innx
              sendbuf2(k,j,1) = x1(k,j,innz-1)
            enddo
          enddo
          do j = 1, inny
            do k = 1, innx
              sendbuf2(k,j,2) = x2(k,j,innz-1)
            enddo
          enddo
          do j = 1, inny
            do k = 1, innx
              sendbuf2(k,j,3) = x3(k,j,innz-1)
            enddo
          enddo
          call MPI_SEND(sendbuf2(1,1,1),nsend,MPI_DOUBLE_PRECISION,right,3,commrow,&
                        ierr)
        endif
        if(myidx.ne.0) then
          call MPI_WAIT(msid,status,ierr)
          do j = 1, inny
            do k = 1, innx
              x1(k,j,1) = recvbuf2(k,j,1)
            enddo
          enddo
          do j = 1, inny
            do k = 1, innx
              x2(k,j,1) = recvbuf2(k,j,2)
            enddo
          enddo
          do j = 1, inny
            do k = 1, innx
              x3(k,j,1) = recvbuf2(k,j,3)
            enddo
          enddo
        endif

        if(myidx.ne.(npx-1)) then
          call MPI_IRECV(recvbuf2(1,1,1),nsend,MPI_DOUBLE_PRECISION,right,4,commrow,&
                         msid,ierr)
        endif
        if(myidx.ne.0) then
          do j = 1, inny
            do k = 1, innx
              sendbuf2(k,j,1) = x1(k,j,2)
            enddo
          enddo
          do j = 1, inny
            do k = 1, innx
              sendbuf2(k,j,2) = x2(k,j,2)
            enddo
          enddo
          do j = 1, inny
            do k = 1, innx
              sendbuf2(k,j,3) = x3(k,j,2)
            enddo
          enddo
          call MPI_SEND(sendbuf2(1,1,1),nsend,MPI_DOUBLE_PRECISION,left,4,commrow, &
                        ierr)
        endif
        if(myidx.ne.(npx-1)) then
          call MPI_WAIT(msid,status,ierr)
          do j = 1, inny
            do k = 1, innx
              x1(k,j,innz) = recvbuf2(k,j,1)
            enddo
          enddo
          do j = 1, inny
            do k = 1, innx
              x2(k,j,innz) = recvbuf2(k,j,2)
            enddo
          enddo
          do j = 1, inny
            do k = 1, innx
              x3(k,j,innz) = recvbuf2(k,j,3)
            enddo
          enddo
        endif

        endif

        call MPI_BARRIER(comm2d,ierr)
        t_boundint = t_boundint + elapsedtime_Timer(t0)

        end subroutine boundint4_Fldmger


        !exchange the information for interpolation using CIC.
        subroutine boundint_Fldmger(x,innx,inny,innz,grid) 
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx, inny, innz
        double precision, dimension(innx,inny,innz), intent(inout) :: x
        type (Pgrid2d), intent(in) :: grid
        double precision, dimension(innx,innz-2) :: sendbuf1, recvbuf1
        double precision, dimension(innx,inny) :: sendbuf2, recvbuf2
        integer :: left,right,bottom,top,msid,&
                   ierr,i,j,k,nx,ny
        integer :: myid,myidx,myidy,comm2d,commcol,commrow,totnp,&
                   npx,npy
        integer status(MPI_STATUS_SIZE)
        double precision :: t0
        integer :: nsend
    
        call starttime_Timer(t0)
    
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)
        call getsize_Pgrid2d(grid,totnp,npy,npx)

        top = myidy + 1
        bottom = myidy - 1
        left = myidx - 1
        right = myidx + 1

        ny = inny - 2

! This could be modified to improve speed
        nsend = innx*(innz-2)
        if(myidy.ne.0) then
          call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,1,commcol,&
                        msid,ierr)
        endif
        if(myidy.ne.(npy-1)) then
          do k = 2, innz-1
            do i = 1, innx
              sendbuf1(i,k-1) = x(i,inny-1,k)
            enddo
          enddo

          call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,1,&
                        commcol,ierr)
        endif
        if(myidy.ne.0) then
          call MPI_WAIT(msid,status,ierr)
          do k = 2, innz-1
            do i = 1, innx
              x(i,1,k) = recvbuf1(i,k-1)
            enddo
          enddo
        endif

        if(myidy.ne.(npy-1)) then
          call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,2,commcol,&
                        msid,ierr)
        endif
        if(myidy.ne.0) then
          do k = 2, innz-1
            do i = 1, innx
              sendbuf1(i,k-1) = x(i,2,k)
            enddo
          enddo
          call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,2,commcol,&
                        ierr)
        endif
        if(myidy.ne.(npy-1)) then
          call MPI_WAIT(msid,status,ierr)
          do k = 2, innz-1
            do i = 1, innx
              x(i,inny,k) = recvbuf1(i,k-1)
            enddo
          enddo
        endif

        nsend = innx*inny
        if(myidx.ne.0) then
          call MPI_IRECV(recvbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,left,3,commrow,&
                         msid,ierr)
        endif
        if(myidx.ne.(npx-1)) then
          do j = 1, inny
            do k = 1, innx
              sendbuf2(k,j) = x(k,j,innz-1)
            enddo
          enddo
          call MPI_SEND(sendbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,right,3,commrow,&
                        ierr)
        endif
        if(myidx.ne.0) then
          call MPI_WAIT(msid,status,ierr)
          do j = 1, inny
            do k = 1, innx
              x(k,j,1) = recvbuf2(k,j)
            enddo
          enddo
        endif

        if(myidx.ne.(npx-1)) then
          call MPI_IRECV(recvbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,right,4,commrow,&
                         msid,ierr)
        endif
        if(myidx.ne.0) then
          do j = 1, inny
            do k = 1, innx
              sendbuf2(k,j) = x(k,j,2)
            enddo
          enddo
          call MPI_SEND(sendbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,left,4,commrow, &
                        ierr)
        endif
        if(myidx.ne.(npx-1)) then
          call MPI_WAIT(msid,status,ierr)
          do j = 1, inny
            do k = 1, innx
              x(k,j,innz) = recvbuf2(k,j)
            enddo
          enddo
        endif

        t_boundint = t_boundint + elapsedtime_Timer(t0)

        end subroutine boundint_Fldmger

        !move the local distributed array xin(nxlc,nylc,nzlc) to
        !a new local distributed array xout(nxnewlc,nynewlc,nzlc).
        subroutine fieldmv1bak_Fldmger(xin,nxlc,nylc,nzlc,nygb,&
          myid,comm,np,lgrid,rgrid,dispold,nxnewlc,nynewlc,xout)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nxlc,nylc,nzlc,nygb
        integer, intent(in) :: myid,comm,np,lgrid,rgrid,nxnewlc,nynewlc
        integer, dimension(0:np-1), intent(in) :: dispold
        double precision, dimension(nxlc,nylc,nzlc), intent(in) :: xin
        double precision, dimension(nxnewlc,nynewlc,nzlc), intent(out) :: xout
        double precision, allocatable, dimension(:,:,:) :: left,right
        double precision, allocatable, dimension(:,:,:) :: templeft,tempright
        double precision, dimension(nxnewlc,nynewlc,nzlc) :: contleft,contright
        integer, dimension(2) :: ileft, iright, jleft, jright
        integer, dimension(nylc) :: msk
        integer, dimension(nygb) :: mmsk
        integer status(MPI_STATUS_SIZE), ierr, msid
        integer :: i,gbi,il,ir,in,nleft0,nright0,myleft,myright
        integer :: iout, totout, ii, iir, iil, j, k

        ileft(1) = 0
        iright(1) = 0
        ileft(2) = myid
        iright(2) = myid
        msk = 2
        do i = 1, nylc
          gbi = dispold(myid) + i
          if(gbi.lt.lgrid) then
            ileft(1) = ileft(1) + 1
            msk(i) = -1
          else if((gbi.le.nygb).and.(gbi.gt.rgrid)) then
            iright(1) = iright(1) + 1
            msk(i) = 1
          else if(gbi.le.nygb) then
            msk(i) = 0
          else
          endif
        enddo

        if(ileft(1).gt.0) then
        allocate(left(nxnewlc,ileft(1),nzlc))
        endif
        if(iright(1).gt.0) then
        allocate(right(nxnewlc,iright(1),nzlc))
        endif
        do k = 1, nzlc
          il = 0
          ir = 0
        do i = 1, nylc
          if(msk(i).eq.(-1)) then
            il = il + 1
            do j = 1, nxnewlc
              if(j.le.nxlc) then
                left(j,il,k) = xin(j,i,k)
              else
                left(j,il,k) = 0.0
              endif
            enddo
          else if(msk(i).eq.1) then
            ir = ir + 1
            do j = 1, nxnewlc
              if(j.le.nxlc) then
                right(j,ir,k) = xin(j,i,k)
              else
                right(j,ir,k) = 0.0
              endif
            enddo
          else
          endif
        enddo
        enddo
 
        nleft0 = 0
        nright0 = 0
        myleft = myid - 1
        myright = myid + 1

        do

        jleft = 0
        jright = 0
        if(myid.ne.0) then
          call MPI_SEND(ileft,2,MPI_INTEGER,myleft,0,comm,ierr)
        endif
        if(myid.ne.(np-1)) then
          call MPI_RECV(jleft,2,MPI_INTEGER,myright,0,comm,status,&
                        ierr)
        endif

        if(myid.ne.(np-1)) then
          call MPI_SEND(iright,2,MPI_INTEGER,myright,0,comm,ierr)
        endif
        if(myid.ne.0) then
          call MPI_RECV(jright,2,MPI_INTEGER,myleft,0,comm,status,&
                        ierr)
        endif

        !send outgoing particles to left neibhoring processor.
        if(jleft(1).gt.0) then
          allocate(templeft(nxnewlc,jleft(1),nzlc))
        endif
        if(myid.ne.0) then
          ileft(1) = nxnewlc*ileft(1)*nzlc
          if(ileft(1).gt.0) then 
          call MPI_SEND(left(1,1,1),ileft(1),MPI_DOUBLE_PRECISION,myleft,0,comm,&
                        ierr)
          endif
          ileft(1) = ileft(1)/(nxnewlc*nzlc)
        endif
        if(myid.ne.(np-1)) then
          jleft(1) = nxnewlc*jleft(1)*nzlc
          if(jleft(1).gt.0) then
          call MPI_RECV(templeft(1,1,1),jleft(1),MPI_DOUBLE_PRECISION,myright, &
                        0,comm,status,ierr)
          endif
          jleft(1) = jleft(1)/(nxnewlc*nzlc) 
        endif
        
        if(ileft(1).gt.0) then
          deallocate(left)
        endif
        ileft(1) = 0
        ileft(2) = jleft(2)
        mmsk = 2
        do i = 1, jleft(1)
          gbi = dispold(jleft(2)) + i
          if(gbi.lt.lgrid) then
            mmsk(i) = -1
            ileft(1) = ileft(1) + 1
          endif
        enddo 
        if(ileft(1).gt.0) then
          allocate(left(nxnewlc,ileft(1),nzlc))
        endif
        do k = 1, nzlc
          iil = 0
          ii = 0
        do i = 1,jleft(1)
          if(mmsk(i).eq.(-1)) then
            iil = iil + 1
            left(:,iil,k) = templeft(:,i,k)
          else
            ii = ii + 1
            contleft(:,ii+nleft0,k) = templeft(:,i,k)
          endif
        enddo
        enddo
        nleft0 = nleft0 + (jleft(1)-ileft(1))
        if(jleft(1).gt.0) then
          deallocate(templeft)
        endif

        !send outgoing particles to right neibhoring processor.
        if(jright(1).gt.0) then
          allocate(tempright(nxnewlc,jright(1),nzlc))
        endif
        if(myid.ne.(np-1)) then
          iright(1) = nxnewlc*iright(1)*nzlc
          if(iright(1).gt.0) then
          call MPI_SEND(right(1,1,1),iright(1),MPI_DOUBLE_PRECISION,myright,&
                        0,comm,ierr)
          endif
          iright(1) = iright(1)/(nxnewlc*nzlc)
        endif
        if(myid.ne.0) then
          jright(1) = nxnewlc*jright(1)*nzlc
          if(jright(1).gt.0) then
          call MPI_RECV(tempright(1,1,1),jright(1),MPI_DOUBLE_PRECISION,myleft,&
                        0,comm,status,ierr)
          endif
          jright(1) = jright(1)/(nxnewlc*nzlc)
        endif
     
        if(jright(1).gt.0) then
          allocate(templeft(nxnewlc,nright0,nzlc))
          do k = 1, nzlc
            do j = 1, nright0
              do i = 1, nxnewlc
                templeft(i,j,k) = contright(i,j,k)
              enddo
            enddo
          enddo
        endif
        if(iright(1).gt.0) then
          deallocate(right)
        endif
        iright(1) = 0
        iright(2) = jright(2)
        mmsk = 2
        do i = 1, jright(1)
          gbi = dispold(jright(2)+1) - jright(1) + i 
          if(gbi.gt.rgrid) then
            mmsk(i) = 1
            iright(1) = iright(1) + 1
          endif
        enddo
        if(iright(1).gt.0) then
          allocate(right(nxnewlc,iright(1),nzlc))
        endif
        do k = 1, nzlc
          iir = 0
          ii = 0
        do i = 1, jright(1)
          if(mmsk(i).eq.1) then
            iir = iir + 1
            right(:,iir,k) =  tempright(:,i,k)
          else
            ii = ii + 1
            contright(:,ii,k) = tempright(:,i,k)
!            contright(:,ii+nright0,k) = tempright(:,i,k)
          endif
        enddo
        enddo
        if(jright(1).gt.0) then
          do k = 1, nzlc
            do j = 1, nright0
              do i = 1, nxnewlc
                contright(i,j+ii,k) = templeft(i,j,k)
              enddo
            enddo
          enddo
          deallocate(templeft)
        endif
        nright0 = nright0 + (jright(1)-iright(1))
        if(jright(1).gt.0) then
          deallocate(tempright)
        endif

        iout = iright(1) + ileft(1)
        call MPI_ALLREDUCE(iout,totout,1,MPI_INTEGER,MPI_SUM, &
                           comm,ierr)
        if(totout.eq.0) then
          exit
        endif

        enddo
            
        if(ileft(1).gt.0) then
          deallocate(left)
        endif
        if(iright(1).gt.0) then
          deallocate(right)     
        endif

        do k = 1, nzlc
        do i = 1, nright0
          do j = 1, nxnewlc
            xout(j,i,k) = contright(j,i,k)        
          enddo
        enddo
        enddo
        do k = 1, nzlc
          in = 0
        do i = 1, nylc
          if(msk(i).eq.0) then
            in = in + 1
            do j = 1, nxnewlc
              if(j.le.nxlc) then
                xout(j,nright0+in,k) = xin(j,i,k)
              else
                xout(j,nright0+in,k) = 0.0
              endif
            enddo
          endif
        enddo
        enddo
        in = in + nright0
        do k = 1, nzlc
        do i = 1, nleft0
          do j = 1, nxnewlc
            xout(j,in+i,k) = contleft(j,i,k)
          enddo
        enddo
        enddo
        in = in + nleft0
        do k = 1, nzlc
        do i = in+1, nynewlc
          do j = 1, nxnewlc
            xout(j,i,k) = 0.0
          enddo
        enddo
        enddo

        end subroutine fieldmv1bak_Fldmger

        !move the local distributed array xin(nxlc,nylc,nzlc) to
        !a new local distributed array xout(nxlc,nylc,nznewlc).
        subroutine fieldmv2bak_Fldmger(xin,nxlc,nylc,nzlc,nzgb,&
          myid,comm,np,lgrid,rgrid,dispold,nznewlc,xout)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nxlc,nylc,nzlc,nzgb
        integer, intent(in) :: myid,comm,np,lgrid,rgrid,nznewlc
        integer, dimension(0:np-1), intent(in) :: dispold
        double precision, dimension(nxlc,nylc,nzlc), intent(in) :: xin
        double precision, dimension(nxlc,nylc,nznewlc), intent(out) :: xout
        double precision, allocatable, dimension(:,:,:) :: left,right
        double precision, allocatable, dimension(:,:,:) :: templeft,tempright
        double precision, dimension(nxlc,nylc,nznewlc) :: contleft,contright
        integer, dimension(2) :: ileft, iright, jleft, jright
        integer, dimension(nzlc) :: msk
        integer, dimension(nzgb) ::  mmsk
        integer status(MPI_STATUS_SIZE), ierr, msid
        integer :: i,gbi,il,ir,in,nleft0,nright0,myleft,myright
        integer :: iout, totout, ii, iir, iil, j, k

        ileft(1) = 0
        iright(1) = 0
        ileft(2) = myid
        iright(2) = myid
        msk = 2
        do i = 1, nzlc
          gbi = dispold(myid) + i
          if(gbi.lt.lgrid) then
            ileft(1) = ileft(1) + 1
            msk(i) = -1
          else if((gbi.le.nzgb).and.(gbi.gt.rgrid)) then
            iright(1) = iright(1) + 1
            msk(i) = 1
          else if(gbi.le.nzgb) then
            msk(i) = 0
          else
          endif
        enddo

        if(ileft(1).gt.0) then
          allocate(left(nxlc,nylc,ileft(1)))
        endif
        if(iright(1).gt.0) then
          allocate(right(nxlc,nylc,iright(1)))
        endif
        do k = 1, nylc
          il = 0
          ir = 0
        do i = 1, nzlc
          if(msk(i).eq.(-1)) then
            il = il + 1
            do j = 1, nxlc
              left(j,k,il) = xin(j,k,i)
            enddo
          else if(msk(i).eq.1) then
            ir = ir + 1
            do j = 1, nxlc
              right(j,k,ir) = xin(j,k,i)
            enddo
          else
          endif
        enddo
        enddo
 
        nleft0 = 0
        nright0 = 0
        myleft = myid - 1
        myright = myid + 1

        do

!        ij = ij + 1

        jleft = 0
        jright = 0
        if(myid.ne.0) then
          call MPI_SEND(ileft,2,MPI_INTEGER,myleft,0,comm,ierr)
        endif
        if(myid.ne.(np-1)) then
          call MPI_RECV(jleft,2,MPI_INTEGER,myright,0,comm,status,&
                        ierr)
        endif

        if(myid.ne.(np-1)) then
          call MPI_SEND(iright,2,MPI_INTEGER,myright,0,comm,ierr)
        endif
        if(myid.ne.0) then
          call MPI_RECV(jright,2,MPI_INTEGER,myleft,0,comm,status,&
                        ierr)
        endif

        !send outgoing particles to left neibhoring processor.
        if(jleft(1).gt.0) then
          allocate(templeft(nxlc,nylc,jleft(1)))
        endif
        if(myid.ne.0) then
          ileft(1) = nxlc*ileft(1)*nylc
          if(ileft(1).gt.0) then
          call MPI_SEND(left(1,1,1),ileft(1),MPI_DOUBLE_PRECISION,myleft,0,comm,&
                        ierr)
          endif 
          ileft(1) = ileft(1)/(nxlc*nylc)
        endif
        if(myid.ne.(np-1)) then
          jleft(1) = nxlc*jleft(1)*nylc
          if(jleft(1).gt.0) then
          call MPI_RECV(templeft(1,1,1),jleft(1),MPI_DOUBLE_PRECISION,myright, &
                        0,comm,status,ierr)
          endif
          jleft(1) = jleft(1)/(nxlc*nylc) 
        endif
        
        if(ileft(1).gt.0) then
          deallocate(left)
        endif
        ileft(1) = 0
        ileft(2) = jleft(2)
        mmsk = 2
        do i = 1, jleft(1)
          gbi = dispold(jleft(2)) + i
          if(gbi.lt.lgrid) then
            mmsk(i) = -1
            ileft(1) = ileft(1) + 1
          endif
        enddo 
        if(ileft(1).gt.0) then
          allocate(left(nxlc,nylc,ileft(1)))
        endif
        do k = 1, nylc
          iil = 0
          ii = 0
        do i = 1,jleft(1)
          if(mmsk(i).eq.(-1)) then
            iil = iil + 1
            do j = 1, nxlc
              left(j,k,iil) = templeft(j,k,i)
            enddo
          else
            ii = ii + 1
            do j = 1, nxlc
              contleft(j,k,ii+nleft0) = templeft(j,k,i)
            enddo
          endif
        enddo
        enddo
        nleft0 = nleft0 + (jleft(1)-ileft(1))
        if(jleft(1).gt.0) then
          deallocate(templeft)
        endif

        !send outgoing particles to right neibhoring processor.
        if(jright(1).gt.0) then
          allocate(tempright(nxlc,nylc,jright(1)))
        endif
        if(myid.ne.(np-1)) then
          iright(1) = nxlc*iright(1)*nylc
          if(iright(1).gt.0) then
          call MPI_SEND(right(1,1,1),iright(1),MPI_DOUBLE_PRECISION,myright,&
                        0,comm,ierr)
          endif
          iright(1) = iright(1)/(nxlc*nylc)
        endif
        if(myid.ne.0) then
          jright(1) = nxlc*jright(1)*nylc
          if(jright(1).gt.0) then
          call MPI_RECV(tempright(1,1,1),jright(1),MPI_DOUBLE_PRECISION,myleft,&
                        0,comm,status,ierr)
          endif
          jright(1) = jright(1)/(nxlc*nylc)
        endif
     
        if(jright(1).gt.0) then
          allocate(templeft(nxlc,nylc,nright0))
          do k = 1, nright0
            do j = 1, nylc
              do i = 1, nxlc
                templeft(i,j,k) = contright(i,j,k)
              enddo
            enddo
          enddo
        endif
        if(iright(1).gt.0) then
          deallocate(right)
        endif
        iright(1) = 0
        iright(2) = jright(2)
        mmsk = 2
        do i = 1, jright(1)
          gbi = dispold(jright(2)+1)-jright(1)+i
          if(gbi.gt.rgrid) then
            mmsk(i) = 1
            iright(1) = iright(1) + 1
          endif
        enddo
        if(iright(1).gt.0) then
          allocate(right(nxlc,nylc,iright(1)))
        endif
        do k = 1, nylc
          iir = 0
          ii = 0
        do i = 1, jright(1)
          if(mmsk(i).eq.1) then
            iir = iir + 1
            right(:,k,iir) =  tempright(:,k,i)
          else
            ii = ii + 1
!            contright(:,k,ii+nright0) = tempright(:,k,i)
            contright(:,k,ii) = tempright(:,k,i)
          endif
        enddo
        enddo
        if(jright(1).gt.0) then
          do k = 1, nright0
            do j = 1, nylc
              do i = 1, nxlc
                contright(i,j,k+ii) = templeft(i,j,k) 
              enddo
            enddo
          enddo
          deallocate(templeft)
        endif
        nright0 = nright0 + (jright(1)-iright(1))
        if(jright(1).gt.0) then
          deallocate(tempright)
        endif

        iout = iright(1) + ileft(1)
        call MPI_ALLREDUCE(iout,totout,1,MPI_INTEGER,MPI_SUM, &
                           comm,ierr)
        if(totout.eq.0) then
          exit
        endif

        enddo
            
        if(ileft(1).gt.0) then
          deallocate(left)
        endif
        if(iright(1).gt.0) then
          deallocate(right)     
        endif

        do k = 1, nylc
        do i = 1, nright0
          do j = 1, nxlc
            xout(j,k,i) = contright(j,k,i)        
          enddo
        enddo
        enddo
        do k = 1, nylc
          in = 0
        do i = 1, nzlc
          if(msk(i).eq.0) then
            in = in + 1
            do j = 1, nxlc
              xout(j,k,nright0+in) = xin(j,k,i)
            enddo
          endif
        enddo
        enddo
        in = in + nright0
        do k = 1, nylc
        do i = 1, nleft0
          do j = 1, nxlc
            xout(j,k,in+i) = contleft(j,k,i)
          enddo
        enddo
        enddo
        in = in + nleft0
        do k = 1, nylc
        do i = in+1, nznewlc
          do j = 1, nxlc
            xout(j,k,i) = 0.0
          enddo
        enddo
        enddo

        end subroutine fieldmv2bak_Fldmger

        !move the local distributed array xin(nxlc,nylc,nzlc) to
        !a new local distributed array xout(nxnewlc,nynewlc,nzlc).
        subroutine fieldmv1_Fldmger(xin,nxlc,nylc,nzlc,nygb,&
          myid,comm,np,lgrid,rgrid,dispold,nxnewlc,nynewlc,xout)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nxlc,nylc,nzlc,nygb
        integer, intent(in) :: myid,comm,np,lgrid,rgrid,nxnewlc,nynewlc
        integer, dimension(0:np-1), intent(in) :: dispold
        double precision, dimension(nxlc,nylc,nzlc), intent(in) :: xin
        double precision, dimension(nxnewlc,nynewlc,nzlc), intent(out) :: xout
        double precision, allocatable, dimension(:,:,:) :: left,right
        double precision, allocatable, dimension(:,:,:) :: templeft,tempright
        double precision, dimension(nxnewlc,nynewlc,nzlc) :: contleft,contright
        integer, dimension(2) :: ileft, iright, jleft, jright
        integer, dimension(nylc) :: msk
        integer, dimension(nygb) :: mmsk
        integer status(MPI_STATUS_SIZE), ierr, msid
        integer :: i,gbi,il,ir,in,nleft0,nright0,myleft,myright
        integer :: iout, totout, ii, iir, iil, j, k
        integer statarry(MPI_STATUS_SIZE,4), req(4)

        ileft(1) = 0
        iright(1) = 0
        ileft(2) = myid
        iright(2) = myid
        msk = 2
        do i = 1, nylc
          gbi = dispold(myid) + i
          if(gbi.lt.lgrid) then
            ileft(1) = ileft(1) + 1
            msk(i) = -1
          else if((gbi.le.nygb).and.(gbi.gt.rgrid)) then
            iright(1) = iright(1) + 1
            msk(i) = 1
          else if(gbi.le.nygb) then
            msk(i) = 0
          else
          endif
        enddo

        allocate(left(nxnewlc,ileft(1),nzlc))
        allocate(right(nxnewlc,iright(1),nzlc))
        do k = 1, nzlc
          il = 0
          ir = 0
        do i = 1, nylc
          if(msk(i).eq.(-1)) then
            il = il + 1
            do j = 1, nxnewlc
              if(j.le.nxlc) then
                left(j,il,k) = xin(j,i,k)
              else
                left(j,il,k) = 0.0
              endif
            enddo
          else if(msk(i).eq.1) then
            ir = ir + 1
            do j = 1, nxnewlc
              if(j.le.nxlc) then
                right(j,ir,k) = xin(j,i,k)
              else
                right(j,ir,k) = 0.0
              endif
            enddo
          else
          endif
        enddo
        enddo
 
        call MPI_BARRIER(comm,ierr)

        nleft0 = 0
        nright0 = 0
        if(myid.ne.0) then
          myleft = myid - 1
        else
          myleft = MPI_PROC_NULL
        endif
        if(myid.ne.(np-1)) then
          myright = myid + 1
        else
          myright = MPI_PROC_NULL
        endif

        do

        jleft = 0
        jright = 0
        
        call MPI_IRECV(jleft,2,MPI_INTEGER,myright,0,comm,req(1),&
                        ierr)
        call MPI_IRECV(jright,2,MPI_INTEGER,myleft,0,comm,req(2),&
                        ierr)
        call MPI_ISEND(ileft,2,MPI_INTEGER,myleft,0,comm,req(3),ierr)
        call MPI_ISEND(iright,2,MPI_INTEGER,myright,0,comm,req(4),ierr)
        call MPI_WAITALL(4,req,statarry,ierr)

        !send outgoing particles to left neibhoring processor.
        allocate(templeft(nxnewlc,jleft(1),nzlc))
        ileft(1) = nxnewlc*ileft(1)*nzlc
        jleft(1) = nxnewlc*jleft(1)*nzlc
        call MPI_IRECV(templeft(1,1,1),jleft(1),MPI_DOUBLE_PRECISION,myright, &
                      0,comm,msid,ierr)
        call MPI_SEND(left(1,1,1),ileft(1),MPI_DOUBLE_PRECISION,myleft,0,comm,&
                      ierr)
        call MPI_WAIT(msid,status,ierr)
        ileft(1) = ileft(1)/(nxnewlc*nzlc)
        jleft(1) = jleft(1)/(nxnewlc*nzlc) 
        
        deallocate(left)
        ileft(1) = 0
        ileft(2) = jleft(2)
        mmsk = 2
        do i = 1, jleft(1)
          gbi = dispold(jleft(2)) + i
          if(gbi.lt.lgrid) then
            mmsk(i) = -1
            ileft(1) = ileft(1) + 1
          endif
        enddo 
        allocate(left(nxnewlc,ileft(1),nzlc))
        do k = 1, nzlc
          iil = 0
          ii = 0
        do i = 1,jleft(1)
          if(mmsk(i).eq.(-1)) then
            iil = iil + 1
            left(:,iil,k) = templeft(:,i,k)
          else
            ii = ii + 1
            contleft(:,ii+nleft0,k) = templeft(:,i,k)
          endif
        enddo
        enddo
        nleft0 = nleft0 + (jleft(1)-ileft(1))
        deallocate(templeft)

        !send outgoing particles to right neibhoring processor.
        allocate(tempright(nxnewlc,jright(1),nzlc))
        iright(1) = nxnewlc*iright(1)*nzlc
        jright(1) = nxnewlc*jright(1)*nzlc
        call MPI_IRECV(tempright(1,1,1),jright(1),MPI_DOUBLE_PRECISION,myleft,&
                      0,comm,msid,ierr)
        call MPI_SEND(right(1,1,1),iright(1),MPI_DOUBLE_PRECISION,myright,&
                      0,comm,ierr)
        call MPI_WAIT(msid,status,ierr)
        iright(1) = iright(1)/(nxnewlc*nzlc)
        jright(1) = jright(1)/(nxnewlc*nzlc)
     
        allocate(templeft(nxnewlc,nright0,nzlc))
          do k = 1, nzlc
            do j = 1, nright0
              do i = 1, nxnewlc
                templeft(i,j,k) = contright(i,j,k)
              enddo
            enddo
          enddo
        deallocate(right)
        iright(1) = 0
        iright(2) = jright(2)
        mmsk = 2
        do i = 1, jright(1)
          gbi = dispold(jright(2)+1) - jright(1) + i 
          if(gbi.gt.rgrid) then
            mmsk(i) = 1
            iright(1) = iright(1) + 1
          endif
        enddo
        allocate(right(nxnewlc,iright(1),nzlc))
        do k = 1, nzlc
          iir = 0
          ii = 0
        do i = 1, jright(1)
          if(mmsk(i).eq.1) then
            iir = iir + 1
            right(:,iir,k) =  tempright(:,i,k)
          else
            ii = ii + 1
            contright(:,ii,k) = tempright(:,i,k)
!            contright(:,ii+nright0,k) = tempright(:,i,k)
          endif
        enddo
        enddo
          do k = 1, nzlc
            do j = 1, nright0
              do i = 1, nxnewlc
                contright(i,j+ii,k) = templeft(i,j,k)
              enddo
            enddo
          enddo
        deallocate(templeft)
        nright0 = nright0 + (jright(1)-iright(1))
        deallocate(tempright)

        iout = iright(1) + ileft(1)
        call MPI_ALLREDUCE(iout,totout,1,MPI_INTEGER,MPI_SUM, &
                           comm,ierr)
        if(totout.eq.0) then
          exit
        endif

        enddo
            
        deallocate(left)
        deallocate(right)     

        do k = 1, nzlc
        do i = 1, nright0
          do j = 1, nxnewlc
            xout(j,i,k) = contright(j,i,k)        
          enddo
        enddo
        enddo
        do k = 1, nzlc
          in = 0
        do i = 1, nylc
          if(msk(i).eq.0) then
            in = in + 1
            do j = 1, nxnewlc
              if(j.le.nxlc) then
                xout(j,nright0+in,k) = xin(j,i,k)
              else
                xout(j,nright0+in,k) = 0.0
              endif
            enddo
          endif
        enddo
        enddo
        in = in + nright0
        do k = 1, nzlc
        do i = 1, nleft0
          do j = 1, nxnewlc
            xout(j,in+i,k) = contleft(j,i,k)
          enddo
        enddo
        enddo
        in = in + nleft0
        do k = 1, nzlc
        do i = in+1, nynewlc
          do j = 1, nxnewlc
            xout(j,i,k) = 0.0
          enddo
        enddo
        enddo

        end subroutine fieldmv1_Fldmger

        !move the local distributed array xin(nxlc,nylc,nzlc) to
        !a new local distributed array xout(nxlc,nylc,nznewlc).
        subroutine fieldmv2_Fldmger(xin,nxlc,nylc,nzlc,nzgb,&
          myid,comm,np,lgrid,rgrid,dispold,nznewlc,xout)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nxlc,nylc,nzlc,nzgb
        integer, intent(in) :: myid,comm,np,lgrid,rgrid,nznewlc
        integer, dimension(0:np-1), intent(in) :: dispold
        double precision, dimension(nxlc,nylc,nzlc), intent(in) :: xin
        double precision, dimension(nxlc,nylc,nznewlc), intent(out) :: xout
        double precision, allocatable, dimension(:,:,:) :: left,right
        double precision, allocatable, dimension(:,:,:) :: templeft,tempright
        double precision, dimension(nxlc,nylc,nznewlc) :: contleft,contright
        integer, dimension(2) :: ileft, iright, jleft, jright
        integer, dimension(nzlc) :: msk
        integer, dimension(nzgb) ::  mmsk
        integer status(MPI_STATUS_SIZE), ierr, msid
        integer :: i,gbi,il,ir,in,nleft0,nright0,myleft,myright
        integer :: iout, totout, ii, iir, iil, j, k
        integer statarry(MPI_STATUS_SIZE,4), req(4)

        ileft(1) = 0
        iright(1) = 0
        ileft(2) = myid
        iright(2) = myid
        msk = 2
        do i = 1, nzlc
          gbi = dispold(myid) + i
          if(gbi.lt.lgrid) then
            ileft(1) = ileft(1) + 1
            msk(i) = -1
          else if((gbi.le.nzgb).and.(gbi.gt.rgrid)) then
            iright(1) = iright(1) + 1
            msk(i) = 1
          else if(gbi.le.nzgb) then
            msk(i) = 0
          else
          endif
        enddo

        allocate(left(nxlc,nylc,ileft(1)))
        allocate(right(nxlc,nylc,iright(1)))
        do k = 1, nylc
          il = 0
          ir = 0
        do i = 1, nzlc
          if(msk(i).eq.(-1)) then
            il = il + 1
            do j = 1, nxlc
              left(j,k,il) = xin(j,k,i)
            enddo
          else if(msk(i).eq.1) then
            ir = ir + 1
            do j = 1, nxlc
              right(j,k,ir) = xin(j,k,i)
            enddo
          else
          endif
        enddo
        enddo
 
        call MPI_BARRIER(comm,ierr)

        nleft0 = 0
        nright0 = 0
!        myleft = myid - 1
!        myright = myid + 1
        if(myid.ne.0) then
          myleft = myid - 1
        else
          myleft = MPI_PROC_NULL
        endif
        if(myid.ne.(np-1)) then
          myright = myid + 1
        else
          myright = MPI_PROC_NULL
        endif

        do

        jleft = 0
        jright = 0
        call MPI_IRECV(jleft,2,MPI_INTEGER,myright,0,comm,req(1),&
                      ierr)
        call MPI_IRECV(jright,2,MPI_INTEGER,myleft,0,comm,req(2),&
                      ierr)
        call MPI_ISEND(ileft,2,MPI_INTEGER,myleft,0,comm,req(3),ierr)
        call MPI_ISEND(iright,2,MPI_INTEGER,myright,0,comm,req(4),ierr)
        call MPI_WAITALL(4,req,statarry,ierr)

        !send outgoing particles to left neibhoring processor.
        allocate(templeft(nxlc,nylc,jleft(1)))
        ileft(1) = nxlc*ileft(1)*nylc
        jleft(1) = nxlc*jleft(1)*nylc
        call MPI_IRECV(templeft(1,1,1),jleft(1),MPI_DOUBLE_PRECISION,myright, &
                      0,comm,msid,ierr)
        call MPI_SEND(left(1,1,1),ileft(1),MPI_DOUBLE_PRECISION,myleft,0,comm,&
                      ierr)
        call MPI_WAIT(msid,status,ierr)
        ileft(1) = ileft(1)/(nxlc*nylc)
        jleft(1) = jleft(1)/(nxlc*nylc) 
        
        deallocate(left)
        ileft(1) = 0
        ileft(2) = jleft(2)
        mmsk = 2
        do i = 1, jleft(1)
          gbi = dispold(jleft(2)) + i
          if(gbi.lt.lgrid) then
            mmsk(i) = -1
            ileft(1) = ileft(1) + 1
          endif
        enddo 
        allocate(left(nxlc,nylc,ileft(1)))
        do k = 1, nylc
          iil = 0
          ii = 0
        do i = 1,jleft(1)
          if(mmsk(i).eq.(-1)) then
            iil = iil + 1
            do j = 1, nxlc
              left(j,k,iil) = templeft(j,k,i)
            enddo
          else
            ii = ii + 1
            do j = 1, nxlc
              contleft(j,k,ii+nleft0) = templeft(j,k,i)
            enddo
          endif
        enddo
        enddo
        nleft0 = nleft0 + (jleft(1)-ileft(1))
        deallocate(templeft)

        !send outgoing particles to right neibhoring processor.
        allocate(tempright(nxlc,nylc,jright(1)))
        iright(1) = nxlc*iright(1)*nylc
        jright(1) = nxlc*jright(1)*nylc
        call MPI_IRECV(tempright(1,1,1),jright(1),MPI_DOUBLE_PRECISION,myleft,&
                      0,comm,msid,ierr)
        call MPI_SEND(right(1,1,1),iright(1),MPI_DOUBLE_PRECISION,myright,&
                      0,comm,ierr)
        call MPI_WAIT(msid,status,ierr)
        iright(1) = iright(1)/(nxlc*nylc)
        jright(1) = jright(1)/(nxlc*nylc)
     
        allocate(templeft(nxlc,nylc,nright0))
          do k = 1, nright0
            do j = 1, nylc
              do i = 1, nxlc
                templeft(i,j,k) = contright(i,j,k)
              enddo
            enddo
          enddo
        deallocate(right)
        iright(1) = 0
        iright(2) = jright(2)
        mmsk = 2
        do i = 1, jright(1)
          gbi = dispold(jright(2)+1)-jright(1)+i
          if(gbi.gt.rgrid) then
            mmsk(i) = 1
            iright(1) = iright(1) + 1
          endif
        enddo
          allocate(right(nxlc,nylc,iright(1)))
        do k = 1, nylc
          iir = 0
          ii = 0
        do i = 1, jright(1)
          if(mmsk(i).eq.1) then
            iir = iir + 1
            right(:,k,iir) =  tempright(:,k,i)
          else
            ii = ii + 1
!            contright(:,k,ii+nright0) = tempright(:,k,i)
            contright(:,k,ii) = tempright(:,k,i)
          endif
        enddo
        enddo
        do k = 1, nright0
          do j = 1, nylc
            do i = 1, nxlc
              contright(i,j,k+ii) = templeft(i,j,k) 
            enddo
          enddo
        enddo
        deallocate(templeft)
        nright0 = nright0 + (jright(1)-iright(1))
        deallocate(tempright)

        iout = iright(1) + ileft(1)
        call MPI_ALLREDUCE(iout,totout,1,MPI_INTEGER,MPI_SUM, &
                           comm,ierr)
        if(totout.eq.0) then
          exit
        endif

        enddo
            
        deallocate(left)
        deallocate(right)     

        do k = 1, nylc
        do i = 1, nright0
          do j = 1, nxlc
            xout(j,k,i) = contright(j,k,i)        
          enddo
        enddo
        enddo
        do k = 1, nylc
          in = 0
        do i = 1, nzlc
          if(msk(i).eq.0) then
            in = in + 1
            do j = 1, nxlc
              xout(j,k,nright0+in) = xin(j,k,i)
            enddo
          endif
        enddo
        enddo
        in = in + nright0
        do k = 1, nylc
        do i = 1, nleft0
          do j = 1, nxlc
            xout(j,k,in+i) = contleft(j,k,i)
          enddo
        enddo
        enddo
        in = in + nleft0
        do k = 1, nylc
        do i = in+1, nznewlc
          do j = 1, nxlc
            xout(j,k,i) = 0.0
          enddo
        enddo
        enddo

        end subroutine fieldmv2_Fldmger

      end module Fldmgerclass
