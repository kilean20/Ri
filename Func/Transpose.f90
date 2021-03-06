!----------------------------------------------------------------
! (c) Copyright, 2001 by the Regents of the University of California.
! Transposeclass: 2D and 3D array parallel transpose class in Linear  
!                 Algebra module of FUNCTION layer.
! Version: 1.0
! Author: Ji Qiang, LANL, 3/7/01
! Description: This class defines a tranpose class which contains 2D and 3D 
!              tranpose functions. (Both arrays are distributed on 2D 
!              Cartisian processor array). 
! Comments:
!----------------------------------------------------------------
      module Transposeclass
      use Timerclass

      contains
      ! Subroutine Trans2d: get the transpose of 2D double precision 
      ! array with uneven distribution.
      subroutine trans2dold_TRANSP(ny,nx,nysizex,nsizex,xin,tempmtr,np,stab,&
                         rtab,comm)
      implicit none
      include 'mpif.h'
      integer, intent(in) :: ny,nx,nysizex,nsizex,np,comm
      double complex, dimension(ny,nsizex), intent(in) :: xin
      integer, dimension(0:np-1), intent(in) :: stab,rtab
      double complex, dimension(nx,nysizex), intent(out) :: tempmtr
      double complex, dimension(ny*nsizex) :: sendbuf
      double complex, dimension(nx*nysizex) :: recvbuf
      integer, dimension(0:np-1) :: sendcount,recvcount,senddisp, &
                                recvdisp,rtabdisp 
      integer :: tempi,ierr,i,j,i0,j0,k,iblock
      double precision :: t0,t1

      call starttime_Timer(t0)

      if(np.ne.1) then

      ! put 2d array (ny*nsizex) into 1d buffer.
      do i = 1, ny
        do j = 1, nsizex
          sendbuf(j+(i-1)*nsizex) = xin(i,j)
        enddo
      enddo

!      do i = 0, np-1
!        sendcount(i) = stab(i)*nsizex
!      enddo
      sendcount = stab*nsizex
      senddisp(0) = 0
      do i = 1, np-1
        senddisp(i) = senddisp(i-1) + sendcount(i-1)
      enddo

!      do i = 0, np-1
!        recvcount(i) = rtab(i)*nysizex
!      enddo
      recvcount = rtab*nysizex
      recvdisp(0) = 0
      do i = 1, np-1
        recvdisp(i) = recvdisp(i-1) + recvcount(i-1)
      enddo

      rtabdisp = recvdisp/nysizex

!      call starttime_Timer(t1)
      ! do all-to-all scatter.
      call MPI_ALLTOALLV(sendbuf,sendcount,senddisp,MPI_DOUBLE_COMPLEX,&
           recvbuf,recvcount,recvdisp,MPI_DOUBLE_COMPLEX,comm,ierr)
!      t_init = t_init + elapsedtime_Timer(t1)

!      print*,"after all to all"
      ! reshape the 1d receiving buffer into 2d array (nx*nysizex).
!      do i = 1, nx*nysizex
!        iblock = 0
!        do j = 0, np-1
!          if(i.gt.recvdisp(j)) then
!            iblock = iblock + 1
!          endif
!        enddo
       
!        k = i - recvdisp(iblock-1)

!        i0 = (k-1)/(recvcount(iblock-1)/nysizex) + 1
!        j0 = k - (i0-1)*(recvcount(iblock-1)/nysizex)
!        tempi = i0
!        i0 = j0
!        j0 = tempi
!        i0 = (recvdisp(iblock-1)/nysizex)+i0
!        tempmtr(i0,j0) = recvbuf(i)
!      enddo
      do i = 0, np-1
        i0 = rtabdisp(i)
        j0 = recvdisp(i)
        do j = 1, nysizex
          do k = 1, rtab(i)
            tempmtr(i0+k,j) = recvbuf(j0+(j-1)*rtab(i)+k)
          enddo
        enddo
      enddo

      else

      do j = 1, nsizex
        do i = 1, ny
          tempmtr(j,i) = xin(i,j)
        enddo
      enddo

      endif
    
      t_transp = t_transp + elapsedtime_Timer(t0)

      end subroutine trans2dold_TRANSP

      ! Subroutine Trans2d: get the transpose of 2D double precision array with
      ! uneven distribution.
      subroutine trans2d_TRANSP(ny,nx,nysizex,nsizex,xin,tempmtr,np,stab,&
                         rtab,comm)
      use Timerclass
      implicit none
      include 'mpif.h'
      integer, intent(in) :: ny,nx,nysizex,nsizex,np,comm
      double complex, dimension(ny,nsizex), intent(in) :: xin
      integer, dimension(0:np-1), intent(in) :: stab,rtab
      double complex, dimension(nx,nysizex), intent(out) :: tempmtr
      double complex, dimension(ny*nsizex) :: sendbuf
      double complex, dimension(nx*nysizex) :: recvbuf
      integer, dimension(0:np-1) :: sendcount,recvcount,senddisp, &
                                recvdisp,rtabdisp,stabdisp 
      integer :: tempi,ierr,i,j,i0,j0,k,iblock
      double precision :: t0,t1

      call starttime_Timer(t0)

      if(np.ne.1) then
       
!      do i = 0, np-1
!        sendcount(i) = stab(i)*nsizex
!      enddo
      sendcount = stab*nsizex
      senddisp(0) = 0
      do i = 1, np-1
        senddisp(i) = senddisp(i-1) + sendcount(i-1)
      enddo

      stabdisp = senddisp/nsizex

      ! put 2d array (ny*nsizex) into 1d buffer.
!      do i = 1, ny
!        do j = 1, nsizex
!          sendbuf(j+(i-1)*nsizex) = xin(i,j)
!        enddo
!      enddo
       do i = 0, np-1 
         i0 = senddisp(i)
         j0 = stabdisp(i)
         do j = 1, nsizex
           do k = 1, stab(i)
             sendbuf(k+(j-1)*stab(i)+i0) = xin(k+j0,j)
           enddo
         enddo
       enddo

!      do i = 0, np-1
!        recvcount(i) = rtab(i)*nysizex
!      enddo
      recvcount = rtab*nysizex
      recvdisp(0) = 0
      do i = 1, np-1
        recvdisp(i) = recvdisp(i-1) + recvcount(i-1)
      enddo

      rtabdisp = recvdisp/nysizex

!      t_enlarge = t_enlarge + elapsedtime_Timer(t0)

      call MPI_BARRIER(comm,ierr)
!      call starttime_Timer(t1)
      ! do all-to-all scatter.
      call MPI_ALLTOALLV(sendbuf,sendcount,senddisp,MPI_DOUBLE_COMPLEX,&
           recvbuf,recvcount,recvdisp,MPI_DOUBLE_COMPLEX,comm,ierr)
      call MPI_BARRIER(comm,ierr)
!      t_init = t_init + elapsedtime_Timer(t1)

!      call starttime_Timer(t1)
!      print*,"after all to all"
      ! reshape the 1d receiving buffer into 2d array (nx*nysizex).
!      do i = 1, nx*nysizex
!        iblock = 0
!        do j = 0, np-1
!          if(i.gt.recvdisp(j)) then
!            iblock = iblock + 1
!          endif
!        enddo
       
!        k = i - recvdisp(iblock-1)

!        i0 = (k-1)/(recvcount(iblock-1)/nysizex) + 1
!        j0 = k - (i0-1)*(recvcount(iblock-1)/nysizex)
!        tempi = i0
!        i0 = j0
!        j0 = tempi
!        i0 = (recvdisp(iblock-1)/nysizex)+i0
!        tempmtr(i0,j0) = recvbuf(i)
!      enddo

!      do i = 0, np-1
!        i0 = rtabdisp(i)
!        j0 = recvdisp(i)
!        do j = 1, nysizex
!          do k = 1, rtab(i)
!            tempmtr(i0+k,j) = recvbuf(j0+(j-1)*rtab(i)+k)
!          enddo
!        enddo
!      enddo

!      do i = 0, np-1
!        i0 = rtabdisp(i)
!        j0 = recvdisp(i)
!        do k = 1, rtab(i)
!          do j = 1, nysizex
!            tempmtr(i0+k,j) = recvbuf(j0+(k-1)*nysizex+j)
!          enddo
!        enddo
!      enddo

      do i = 0, np-1
        i0 = rtabdisp(i)
        j0 = recvdisp(i)
        do j = 1, nysizex
          do k = 1, rtab(i)
            tempmtr(i0+k,j) = recvbuf(j0+(k-1)*nysizex+j)
          enddo
        enddo
      enddo

      else

      do j = 1, nsizex
        do i = 1, ny
          tempmtr(j,i) = xin(i,j)
        enddo
      enddo

      endif

!      t_shrink = t_shrink + elapsedtime_Timer(t1)
    
      t_transp = t_transp + elapsedtime_Timer(t0)

      end subroutine trans2d_TRANSP

      ! Subroutine Trans2d: get the transpose of 2D double precision array with
      ! uneven distribution.
      subroutine trans2d2_TRANSP(ny,nx,nysizex,nsizex,xin,xout,np,stab,&
                         rtab,comm,myid)
      use Timerclass
      implicit none
      include 'mpif.h'
      integer, intent(in) :: ny,nx,nysizex,nsizex,np,comm,myid
      double complex, dimension(ny,nsizex), intent(in) :: xin
      integer, dimension(0:np-1), intent(in) :: stab,rtab
      double complex, dimension(nx,nysizex), intent(out) :: xout
      double complex, allocatable, dimension(:,:) :: s
      double complex, allocatable, dimension(:,:) :: t
      integer, dimension(0:np-1) :: senddisp, &
                                recvdisp 
      integer :: ierr,i,j,k,kstrt,ks,ir,msid,sendnum,recnum,joff,koff
      integer status(MPI_STATUS_SIZE)
      double precision :: t0

      call starttime_Timer(t0)

      if(np.ne.1) then

        senddisp(0) = 0
        do i = 1, np-1
          senddisp(i) = senddisp(i-1) + stab(i-1)
        enddo
       
        recvdisp(0) = 0
        do i = 1, np-1
          recvdisp(i) = recvdisp(i-1) + rtab(i-1)
        enddo

        kstrt = myid + 1
        ks = myid -1
        do i = 1, np
          ir = i - (1+ks)
          if(ir.lt.1) ir = ir + np
          ! post receive
          allocate(t(nysizex,rtab(ir-1)))
          recnum = rtab(ir-1)*nysizex
          if(ir.ne.kstrt) call MPI_IRECV(t,recnum,MPI_DOUBLE_COMPLEX,&
                          ir-1,i,comm,msid,ierr)
          ! send data
          allocate(s(stab(ir-1),nsizex))
          joff = senddisp(ir-1)
          do k = 1, nsizex
            do j = 1, stab(ir-1)
              s(j,k) = xin(j+joff,k)
            enddo
          enddo
          ! copy data to oneself directly.
          if(ir.eq.kstrt) then
            do k = 1, nsizex
              do j = 1, stab(ir-1)
                t(j,k) = s(j,k)
              enddo
            enddo
          else
            sendnum = stab(ir-1)*nsizex
            call MPI_SEND(s,sendnum,MPI_DOUBLE_COMPLEX,ir-1,i,comm,ierr)
          endif
          !receive data
          if(ir.ne.kstrt) call MPI_WAIT(msid,status,ierr)
          koff = recvdisp(ir-1)
          do k = 1, rtab(ir-1)
            do j = 1, nysizex
              xout(k+koff,j) = t(j,k)
            enddo
          enddo
          deallocate(s)
          deallocate(t)
        enddo  

      else

      do j = 1, nsizex
        do i = 1, ny
          xout(j,i) = xin(i,j)
        enddo
      enddo

      endif

      t_transp = t_transp + elapsedtime_Timer(t0)

      end subroutine trans2d2_TRANSP

      ! Subroutine Trans3d2: get the transpose of 3D double precision array with
      ! uneven distribution.
      subroutine trans3d2_TRANSP(ny,nx,nysizex,nsizex,xin,xout,np,stab,&
                         rtab,comm,myid,nz)
      use Timerclass
      implicit none
      include 'mpif.h'
      integer, intent(in) :: ny,nx,nysizex,nsizex,np,comm,myid,nz
      double complex, dimension(ny,nsizex,nz), intent(in) :: xin
      integer, dimension(0:np-1), intent(in) :: stab,rtab
      double complex, dimension(nx,nysizex,nz), intent(out) :: xout
      double complex, allocatable, dimension(:,:,:) :: s
      double complex, allocatable, dimension(:,:,:) :: t
      integer, dimension(0:np-1) :: senddisp, &
                                recvdisp 
      integer :: ierr,i,j,k,kstrt,ks,ir,msid,sendnum,recnum,joff,koff
      integer status(MPI_STATUS_SIZE)
      integer :: k0
      double precision :: t0

!      call MPI_BARRIER(comm,ierr)
      call starttime_Timer(t0)

      if(np.ne.1) then

        senddisp(0) = 0
        do i = 1, np-1
          senddisp(i) = senddisp(i-1) + stab(i-1)
        enddo
       
        recvdisp(0) = 0
        do i = 1, np-1
          recvdisp(i) = recvdisp(i-1) + rtab(i-1)
        enddo

        kstrt = myid + 1
        ks = myid -1
        do i = 1, np
          ir = i - (1+ks)
          if(ir.lt.1) ir = ir + np
          ! post receive
          allocate(t(nysizex,rtab(ir-1),nz))
          recnum = rtab(ir-1)*nysizex*nz
          if(ir.ne.kstrt) call MPI_IRECV(t,recnum,MPI_DOUBLE_COMPLEX,&
                          ir-1,i,comm,msid,ierr)
          ! send data
          allocate(s(stab(ir-1),nsizex,nz))
          joff = senddisp(ir-1)
          do k0 = 1, nz
          do k = 1, nsizex
            do j = 1, stab(ir-1)
              s(j,k,k0) = xin(j+joff,k,k0)
            enddo
          enddo
          enddo
          ! copy data to oneself directly.
          if(ir.eq.kstrt) then
            do k0 = 1, nz
            do k = 1, nsizex
              do j = 1, stab(ir-1)
                t(j,k,k0) = s(j,k,k0)
              enddo
            enddo
            enddo
          else
            sendnum = stab(ir-1)*nsizex*nz
            call MPI_SEND(s,sendnum,MPI_DOUBLE_COMPLEX,ir-1,i,comm,ierr)
          endif
          !receive data
          if(ir.ne.kstrt) call MPI_WAIT(msid,status,ierr)
          koff = recvdisp(ir-1)
          do k0 = 1, nz
          do k = 1, rtab(ir-1)
            do j = 1, nysizex
              xout(k+koff,j,k0) = t(j,k,k0)
            enddo
          enddo
          enddo
          deallocate(s)
          deallocate(t)
        enddo  

      else

      do k0 = 1, nz
      do j = 1, nsizex
        do i = 1, ny
          xout(j,i,k0) = xin(i,j,k0)
        enddo
      enddo
      enddo

      endif

!      call MPI_BARRIER(comm,ierr)
      t_transp = t_transp + elapsedtime_Timer(t0)

      end subroutine trans3d2_TRANSP

      ! Subroutine Trans3d3: get the transpose of 3D double precision array with
      ! uneven distribution.
      subroutine trans3d3_TRANSP(nx,nsizey,nsizez,nsizexz,xin,xout,np,stab,&
                         rtab,comm,myid,nz)
      use Timerclass
      implicit none
      include 'mpif.h'
      integer, intent(in) :: nx,nsizey,nsizez,nsizexz,np,comm,myid,nz
      double complex, dimension(nx,nsizey,nsizez), intent(in) :: xin
      integer, dimension(0:np-1), intent(in) :: stab,rtab
      double complex, dimension(nz,nsizey,nsizexz), intent(out) :: xout
      double complex, allocatable, dimension(:,:,:) :: s
      double complex, allocatable, dimension(:,:,:) :: t
      integer, dimension(0:np-1) :: senddisp, &
                                recvdisp 
      integer :: ierr,i,j,k,kstrt,ks,ir,msid,sendnum,recnum,joff,koff
      integer status(MPI_STATUS_SIZE)
      integer :: k0
      double precision :: t0

!      call MPI_BARRIER(comm,ierr)
      call starttime_Timer(t0)

      if(np.ne.1) then

        senddisp(0) = 0
        do i = 1, np-1
          senddisp(i) = senddisp(i-1) + stab(i-1)
        enddo
       
        recvdisp(0) = 0
        do i = 1, np-1
          recvdisp(i) = recvdisp(i-1) + rtab(i-1)
        enddo

        kstrt = myid + 1
        ks = myid -1
        do i = 1, np
          ir = i - (1+ks)
          if(ir.lt.1) ir = ir + np
          ! post receive
          allocate(t(nsizexz,nsizey,rtab(ir-1)))
          recnum = rtab(ir-1)*nsizey*nsizexz
          if(ir.ne.kstrt) call MPI_IRECV(t,recnum,MPI_DOUBLE_COMPLEX,&
                          ir-1,i,comm,msid,ierr)
          ! send data
          allocate(s(stab(ir-1),nsizey,nsizez))
          joff = senddisp(ir-1)
          do k0 = 1, nsizez
          do k = 1, nsizey
            do j = 1, stab(ir-1)
              s(j,k,k0) = xin(j+joff,k,k0)
            enddo
          enddo
          enddo
          ! copy data to oneself directly.
          if(ir.eq.kstrt) then
            do k0 = 1, nsizez
            do k = 1, nsizey
              do j = 1, stab(ir-1)
                t(j,k,k0) = s(j,k,k0)
              enddo
            enddo
            enddo
          else
            sendnum = stab(ir-1)*nsizey*nsizez
            call MPI_SEND(s,sendnum,MPI_DOUBLE_COMPLEX,ir-1,i,comm,ierr)
          endif
          !receive data
          if(ir.ne.kstrt) call MPI_WAIT(msid,status,ierr)
          koff = recvdisp(ir-1)
          do j = 1, rtab(ir-1) 
            do k = 1, nsizey
              do k0 = 1, nsizexz
                xout(j+koff,k,k0) = t(k0,k,j)
              enddo
            enddo
          enddo
          deallocate(s)
          deallocate(t)
        enddo  

      else

      do k0 = 1, nsizez
      do j = 1, nsizey
        do i = 1, nx
          xout(k0,j,i) = xin(i,j,k0)
        enddo
      enddo
      enddo

      endif

!      call MPI_BARRIER(comm,ierr)
      t_transp = t_transp + elapsedtime_Timer(t0)

      end subroutine trans3d3_TRANSP

      ! Subroutine Trans3d: get the transpose of 3D double precision array with
      ! uneven distribution.
      subroutine trans3d_TRANSP(ny,nx,nysizex,nsizex,xin,tempmtr,np,stab,&
                         rtab,comm,nzlcal)
      use Timerclass
      implicit none
      include 'mpif.h'
      integer, intent(in) :: ny,nx,nysizex,nsizex,np,comm,nzlcal
      double complex, dimension(ny,nsizex,nzlcal), intent(in) :: xin
      integer, dimension(0:np-1), intent(in) :: stab,rtab
      double complex, dimension(nx,nysizex,nzlcal), intent(out) :: tempmtr
      double complex, dimension(ny*nsizex*nzlcal) :: sendbuf
      double complex, dimension(nx*nysizex*nzlcal) :: recvbuf
      integer, dimension(0:np-1) :: sendcount,recvcount,senddisp, &
                                recvdisp,rtabdisp,stabdisp 
      integer :: tempi,ierr,i,j,i0,j0,k,iblock
      double precision :: t0,t1
      integer :: ntmp,m

!      call MPI_BARRIER(comm,ierr)
      call starttime_Timer(t0)

      ntmp = nsizex*nzlcal

!      print*,"trans3np: ",np,ntmp,nsizex,nzlcal
      if(np.ne.1) then
       
!      do i = 0, np-1
!        sendcount(i) = stab(i)*ntmp
!      enddo
      sendcount = stab*ntmp
      senddisp(0) = 0
      do i = 1, np-1
        senddisp(i) = senddisp(i-1) + sendcount(i-1)
      enddo

      stabdisp = senddisp/ntmp

      ! put 3d array (ny*nsizex) into 1d buffer.
!      do i = 1, ny
!        do j = 1, nsizex
!          sendbuf(j+(i-1)*nsizex) = xin(i,j)
!        enddo
!      enddo
       do i = 0, np-1 
         i0 = senddisp(i)
         j0 = stabdisp(i)
         do m = 1, nzlcal
         do j = 1, nsizex
           do k = 1, stab(i)
             sendbuf(k+(j-1)*stab(i)+(m-1)*nsizex*stab(i)+i0) = &
             xin(k+j0,j,m)
           enddo
         enddo
         enddo
       enddo

       ntmp = nysizex*nzlcal
!      do i = 0, np-1
!        recvcount(i) = rtab(i)*ntmp
!      enddo
      recvcount = rtab*ntmp
      recvdisp(0) = 0
      do i = 1, np-1
        recvdisp(i) = recvdisp(i-1) + recvcount(i-1)
      enddo

      rtabdisp = recvdisp/ntmp

!      t_enlarge = t_enlarge + elapsedtime_Timer(t0)

!      call MPI_BARRIER(comm,ierr)
!      call starttime_Timer(t1)
      ! do all-to-all scatter.
      call MPI_ALLTOALLV(sendbuf,sendcount,senddisp,MPI_DOUBLE_COMPLEX,&
           recvbuf,recvcount,recvdisp,MPI_DOUBLE_COMPLEX,comm,ierr)
!      call MPI_BARRIER(comm,ierr)
!      t_init = t_init + elapsedtime_Timer(t1)

!      call starttime_Timer(t1)
!      print*,"after all to all"
      ! reshape the 1d receiving buffer into 2d array (nx*nysizex).
!      do i = 1, nx*nysizex
!        iblock = 0
!        do j = 0, np-1
!          if(i.gt.recvdisp(j)) then
!            iblock = iblock + 1
!          endif
!        enddo
       
!        k = i - recvdisp(iblock-1)

!        i0 = (k-1)/(recvcount(iblock-1)/nysizex) + 1
!        j0 = k - (i0-1)*(recvcount(iblock-1)/nysizex)
!        tempi = i0
!        i0 = j0
!        j0 = tempi
!        i0 = (recvdisp(iblock-1)/nysizex)+i0
!        tempmtr(i0,j0) = recvbuf(i)
!      enddo

!      do i = 0, np-1
!        i0 = rtabdisp(i)
!        j0 = recvdisp(i)
!        do j = 1, nysizex
!          do k = 1, rtab(i)
!            tempmtr(i0+k,j) = recvbuf(j0+(j-1)*rtab(i)+k)
!          enddo
!        enddo
!      enddo

!      do i = 0, np-1
!        i0 = rtabdisp(i)
!        j0 = recvdisp(i)
!        do k = 1, rtab(i)
!          do j = 1, nysizex
!            tempmtr(i0+k,j) = recvbuf(j0+(k-1)*nysizex+j)
!          enddo
!        enddo
!      enddo

      do i = 0, np-1
        i0 = rtabdisp(i)
        j0 = recvdisp(i)
        do m = 1, nzlcal
        do j = 1, nysizex
          do k = 1, rtab(i)
            tempmtr(i0+k,j,m) = & 
            recvbuf(j0+(k-1)*nysizex+j+(m-1)*nysizex*rtab(i))
          enddo
        enddo
        enddo
      enddo

      else

      do m = 1, nzlcal
      do j = 1, nsizex
        do i = 1, ny
          tempmtr(j,i,m) = xin(i,j,m)
        enddo
      enddo
      enddo

      endif

!      t_shrink = t_shrink + elapsedtime_Timer(t1)
    
!      call MPI_BARRIER(comm,ierr)
      t_transp = t_transp + elapsedtime_Timer(t0)

      end subroutine trans3d_TRANSP

      ! Subroutine Trans3d3: get the transpose of 3D double precision array with
      ! uneven distribution.
      subroutine trans3d3r_TRANSP(nx,nsizey,nsizez,nsizexz,xin,xout,np,stab,&
                         rtab,comm,myid,nz)
      use Timerclass
      implicit none
      include 'mpif.h'
      integer, intent(in) :: nx,nsizey,nsizez,nsizexz,np,comm,myid,nz
      double precision, dimension(nx,nsizey,nsizez), intent(in) :: xin
      integer, dimension(0:np-1), intent(in) :: stab,rtab
      double precision, dimension(nz,nsizey,nsizexz), intent(out) :: xout
      double precision, allocatable, dimension(:,:,:) :: s
      double precision, allocatable, dimension(:,:,:) :: t
      integer, dimension(0:np-1) :: senddisp, &
                                recvdisp 
      integer :: ierr,i,j,k,kstrt,ks,ir,msid,sendnum,recnum,joff,koff
      integer status(MPI_STATUS_SIZE)
      integer :: k0
      double precision :: t0

!      call MPI_BARRIER(comm,ierr)
      call starttime_Timer(t0)

      if(np.ne.1) then

        senddisp(0) = 0
        do i = 1, np-1
          senddisp(i) = senddisp(i-1) + stab(i-1)
        enddo
       
        recvdisp(0) = 0
        do i = 1, np-1
          recvdisp(i) = recvdisp(i-1) + rtab(i-1)
        enddo

        kstrt = myid + 1
        ks = myid -1
        do i = 1, np
          ir = i - (1+ks)
          if(ir.lt.1) ir = ir + np
          ! post receive
          allocate(t(nsizexz,nsizey,rtab(ir-1)))
          recnum = rtab(ir-1)*nsizey*nsizexz
          if(ir.ne.kstrt) call MPI_IRECV(t,recnum,MPI_DOUBLE_PRECISION,&
                          ir-1,i,comm,msid,ierr)
          ! send data
          allocate(s(stab(ir-1),nsizey,nsizez))
          joff = senddisp(ir-1)
          do k0 = 1, nsizez
          do k = 1, nsizey
            do j = 1, stab(ir-1)
              s(j,k,k0) = xin(j+joff,k,k0)
            enddo
          enddo
          enddo
          ! copy data to oneself directly.
          if(ir.eq.kstrt) then
            do k0 = 1, nsizez
            do k = 1, nsizey
              do j = 1, stab(ir-1)
                t(j,k,k0) = s(j,k,k0)
              enddo
            enddo
            enddo
          else
            sendnum = stab(ir-1)*nsizey*nsizez
            call MPI_SEND(s,sendnum,MPI_DOUBLE_PRECISION,ir-1,i,comm,ierr)
          endif
          !receive data
          if(ir.ne.kstrt) call MPI_WAIT(msid,status,ierr)
          koff = recvdisp(ir-1)
          do j = 1, rtab(ir-1) 
            do k = 1, nsizey
              do k0 = 1, nsizexz
                xout(j+koff,k,k0) = t(k0,k,j)
              enddo
            enddo
          enddo
          deallocate(s)
          deallocate(t)
        enddo  

      else

      do k0 = 1, nsizez
      do j = 1, nsizey
        do i = 1, nx
          xout(k0,j,i) = xin(i,j,k0)
        enddo
      enddo
      enddo

      endif

!      call MPI_BARRIER(comm,ierr)
      t_transp = t_transp + elapsedtime_Timer(t0)

      end subroutine trans3d3r_TRANSP

      ! Subroutine Trans3d: get the transpose of 3D double precision array with
      ! uneven distribution.
      subroutine trans3dr_TRANSP(ny,nx,nysizex,nsizex,xin,tempmtr,np,stab,&
                         rtab,comm,nzlcal)
      use Timerclass
      implicit none
      include 'mpif.h'
      integer, intent(in) :: ny,nx,nysizex,nsizex,np,comm,nzlcal
      double precision, dimension(ny,nsizex,nzlcal), intent(in) :: xin
      integer, dimension(0:np-1), intent(in) :: stab,rtab
      double precision, dimension(nx,nysizex,nzlcal), intent(out) :: tempmtr
      double precision, dimension(ny*nsizex*nzlcal) :: sendbuf
      double precision, dimension(nx*nysizex*nzlcal) :: recvbuf
      integer, dimension(0:np-1) :: sendcount,recvcount,senddisp, &
                                recvdisp,rtabdisp,stabdisp 
      integer :: tempi,ierr,i,j,i0,j0,k,iblock
      double precision :: t0,t1
      integer :: ntmp,m

!      call MPI_BARRIER(comm,ierr)
      call starttime_Timer(t0)

      ntmp = nsizex*nzlcal

      if(np.ne.1) then
       
!      do i = 0, np-1
!        sendcount(i) = stab(i)*ntmp
!      enddo
      sendcount = stab*ntmp
      senddisp(0) = 0
      do i = 1, np-1
        senddisp(i) = senddisp(i-1) + sendcount(i-1)
      enddo

      stabdisp = senddisp/ntmp

      ! put 3d array (ny*nsizex) into 1d buffer.
!      do i = 1, ny
!        do j = 1, nsizex
!          sendbuf(j+(i-1)*nsizex) = xin(i,j)
!        enddo
!      enddo
       do i = 0, np-1 
         i0 = senddisp(i)
         j0 = stabdisp(i)
         do m = 1, nzlcal
         do j = 1, nsizex
           do k = 1, stab(i)
             sendbuf(k+(j-1)*stab(i)+(m-1)*nsizex*stab(i)+i0) = &
             xin(k+j0,j,m)
           enddo
         enddo
         enddo
       enddo

       ntmp = nysizex*nzlcal
!      do i = 0, np-1
!        recvcount(i) = rtab(i)*ntmp
!      enddo
      recvcount = rtab*ntmp
      recvdisp(0) = 0
      do i = 1, np-1
        recvdisp(i) = recvdisp(i-1) + recvcount(i-1)
      enddo

      rtabdisp = recvdisp/ntmp

!      t_enlarge = t_enlarge + elapsedtime_Timer(t0)

!      call MPI_BARRIER(comm,ierr)
!      call starttime_Timer(t1)
      ! do all-to-all scatter.
      call MPI_ALLTOALLV(sendbuf,sendcount,senddisp,MPI_DOUBLE_PRECISION,&
           recvbuf,recvcount,recvdisp,MPI_DOUBLE_PRECISION,comm,ierr)
!      call MPI_BARRIER(comm,ierr)
!      t_init = t_init + elapsedtime_Timer(t1)

!      call starttime_Timer(t1)
!      print*,"after all to all"
      ! reshape the 1d receiving buffer into 2d array (nx*nysizex).
!      do i = 1, nx*nysizex
!        iblock = 0
!        do j = 0, np-1
!          if(i.gt.recvdisp(j)) then
!            iblock = iblock + 1
!          endif
!        enddo
       
!        k = i - recvdisp(iblock-1)

!        i0 = (k-1)/(recvcount(iblock-1)/nysizex) + 1
!        j0 = k - (i0-1)*(recvcount(iblock-1)/nysizex)
!        tempi = i0
!        i0 = j0
!        j0 = tempi
!        i0 = (recvdisp(iblock-1)/nysizex)+i0
!        tempmtr(i0,j0) = recvbuf(i)
!      enddo

!      do i = 0, np-1
!        i0 = rtabdisp(i)
!        j0 = recvdisp(i)
!        do j = 1, nysizex
!          do k = 1, rtab(i)
!            tempmtr(i0+k,j) = recvbuf(j0+(j-1)*rtab(i)+k)
!          enddo
!        enddo
!      enddo

!      do i = 0, np-1
!        i0 = rtabdisp(i)
!        j0 = recvdisp(i)
!        do k = 1, rtab(i)
!          do j = 1, nysizex
!            tempmtr(i0+k,j) = recvbuf(j0+(k-1)*nysizex+j)
!          enddo
!        enddo
!      enddo

      do i = 0, np-1
        i0 = rtabdisp(i)
        j0 = recvdisp(i)
        do m = 1, nzlcal
        do j = 1, nysizex
          do k = 1, rtab(i)
            tempmtr(i0+k,j,m) = & 
            recvbuf(j0+(k-1)*nysizex+j+(m-1)*nysizex*rtab(i))
          enddo
        enddo
        enddo
      enddo

      else

      do m = 1, nzlcal
      do j = 1, nsizex
        do i = 1, ny
          tempmtr(j,i,m) = xin(i,j,m)
        enddo
      enddo
      enddo

      endif

!      t_shrink = t_shrink + elapsedtime_Timer(t1)
    
!      call MPI_BARRIER(comm,ierr)
      t_transp = t_transp + elapsedtime_Timer(t0)

      end subroutine trans3dr_TRANSP

      end module Transposeclass
