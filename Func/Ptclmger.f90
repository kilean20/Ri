!----------------------------------------------------------------
! (c) Copyright, 2003 by the Regents of the University of California.
! Ptclmgerclass: Particles moving manager class in Communication module 
!                of FUNCTION layer.
! Version: 2.0
! Author: Ji Qiang, LBNL, 7/24/03
! Description: This class defines functions to transport particles to 
!              their local compuatation processor domain through an
!              iterative neighboring processor communcation process.
! Comments: We have added 3 more attributes to the particle array:
!           x,px,y,py,t,pt,charge/mass,charge weight,id
!----------------------------------------------------------------
        module Ptclmgerclass
          use Timerclass
          use Pgrid2dclass
        contains
        ! move particles from one processor to 4 neighboring processors.
        subroutine ptsmv_Ptclmger(Ptsl,Nptlocal,grid,pdim,npmax,lcrange)
        implicit none
        include 'mpif.h'
        type (Pgrid2d), intent(in) :: grid
        integer, intent(inout) :: Nptlocal
        integer, intent(in) :: pdim,npmax
        double precision, pointer, dimension(:,:) :: Ptsl
        double precision, dimension(:),intent(in) :: lcrange
!        double precision, dimension(pdim,npmax) :: Ptsl
        double precision, allocatable, dimension(:,:) :: left,right,up,down
        double precision, allocatable, dimension(:,:) :: temp1,recv
        integer :: myid,myidx,myidy,totnp,npy,npx, &
                   comm2d,commcol,commrow
        integer :: ileft,iright,iup,idown,iupright,iupleft,&
                   idownleft,idownright
        integer :: jleft,jright,jup,jdown,jupright,jupleft,&
                   jdownleft,jdownright
        integer :: myleft,myright,myup,mydown,myupright,myupleft,&
                   mydwnleft,mydwnright
        integer :: nsmall,i,j,numpts,ic
        integer, dimension(2) :: tmpcoord
        logical, allocatable, dimension(:) :: msk,mmsk
        integer :: numbuf,nmv,nmv0,nout,ii,totnmv,ij
        integer :: msid,ierr
        integer status(MPI_STATUS_SIZE) 
        integer statarry(MPI_STATUS_SIZE,4), req(4)
        double precision :: t0,t1

!        call starttime_Timer(t1)

        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        if(Nptlocal.gt.480) then
          nsmall = Nptlocal/16
        else
          nsmall = Nptlocal
        endif
        allocate(left(9,nsmall))
        allocate(right(9,nsmall))
        allocate(up(9,nsmall))
        allocate(down(9,nsmall))

!        left = 0.0
!        right = 0.0
!        up = 0.0
!        down = 0.0
        ileft = 0
        iright = 0
        iup = 0
        idown = 0

        allocate(msk(Nptlocal))
!        msk = .true.

        do i = 1, Nptlocal
          msk(i) = .true.
          if(Ptsl(5,i).le.lcrange(5)) then
            if(myidx.ne.0) then
              ileft = ileft + 1
              left(:,ileft) = Ptsl(:,i)
              msk(i) = .false.
              if(mod(ileft,nsmall).eq.0) then
                allocate(temp1(9,ileft))
                do j = 1,ileft
                  temp1(:,j) = left(:,j)
                enddo
                deallocate(left)
                allocate(left(9,ileft+nsmall))
                do j = 1,ileft
                  left(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            else
              if(Ptsl(3,i).gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = Ptsl(:,i)
                msk(i) = .false.
                if(mod(iup,nsmall).eq.0) then
                  allocate(temp1(9,iup))
                  do j = 1, iup
                    temp1(:,j) = up(:,j)
                  enddo
                  deallocate(up)
                  allocate(up(9,iup+nsmall))
                  do j = 1, iup
                    up(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else if(Ptsl(3,i).le.lcrange(3)) then
                if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = Ptsl(:,i)
                msk(i) = .false.
                if(mod(idown,nsmall).eq.0) then
                  allocate(temp1(9,idown))
                  do j = 1, idown
                    temp1(:,j) = down(:,j)
                  enddo
                  deallocate(down)
                  allocate(down(9,idown+nsmall))
                  do j = 1, idown
                    down(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else
              endif
            endif
          else if(Ptsl(5,i).gt.lcrange(6)) then
            if(myidx.ne.(npx-1)) then
              iright = iright + 1
              right(:,iright) = Ptsl(:,i)
              msk(i) = .false.
              if(mod(iright,nsmall).eq.0) then
                allocate(temp1(9,iright))
                do j =1, iright
                  temp1(:,j) = right(:,j)
                enddo
                deallocate(right)
                allocate(right(9,iright+nsmall))
                do j =1, iright
                  right(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            else
              if(Ptsl(3,i).gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = Ptsl(:,i)
                msk(i) = .false.
                if(mod(iup,nsmall).eq.0) then
                  allocate(temp1(9,iup))
                  do j = 1, iup
                    temp1(:,j) = up(:,j)
                  enddo
                  deallocate(up)
                  allocate(up(9,iup+nsmall))
                  do j = 1, iup
                    up(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else if(Ptsl(3,i).le.lcrange(3)) then
                if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = Ptsl(:,i)
                msk(i) = .false.
                if(mod(idown,nsmall).eq.0) then
                  allocate(temp1(9,idown))
                  do j = 1, idown
                    temp1(:,j) = down(:,j)
                  enddo
                  deallocate(down)
                  allocate(down(9,idown+nsmall))
                  do j = 1, idown
                    down(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else
              endif
            endif
          else if(Ptsl(3,i).gt.lcrange(4)) then
            if(myidy.ne.(npy-1)) then
              iup = iup + 1
              up(:,iup) = Ptsl(:,i)
              msk(i) = .false.
              if(mod(iup,nsmall).eq.0) then
                allocate(temp1(9,iup))
                do j = 1, iup
                  temp1(:,j) = up(:,j)
                enddo
                deallocate(up)
                allocate(up(9,iup+nsmall))
                do j = 1, iup
                  up(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            endif
          else if(Ptsl(3,i).le.lcrange(3)) then
            if(myidy.ne.0) then
              idown = idown + 1
              down(:,idown) = Ptsl(:,i)
              msk(i) = .false.
              if(mod(idown,nsmall).eq.0) then
                allocate(temp1(9,idown))
                do j = 1, idown
                  temp1(:,j) = down(:,j)
                enddo
                deallocate(down)
                allocate(down(9,idown+nsmall))
                do j = 1, idown
                  down(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            endif
          else
          endif

        enddo

!        t_ptsmv1 = t_ptsmv1 + elapsedtime_Timer(t1)
        nmv0 = 0
        nout = ileft+iright+iup+idown
!        print*,"in comm: ",Nptlocal,nout
        allocate(recv(9,nmv0+1))
        allocate(temp1(9,nmv0+1))
        ij = 0
        call MPI_BARRIER(comm2d,ierr)
        if(myidx.ne.(npx-1)) then
          myright = myidx + 1
        else
          myright = MPI_PROC_NULL
        endif
        if(myidx.ne.0) then
          myleft = myidx - 1
        else
          myleft = MPI_PROC_NULL
        endif 

        if(myidy.ne.npy-1) then
          myup = myidy + 1
        else
          myup = MPI_PROC_NULL
        endif
        if(myidy.ne.0) then
          mydown = myidy -1
        else
          mydown = MPI_PROC_NULL
        endif

        do

        ij = ij + 1
        jleft = 0
        jright = 0

        call MPI_IRECV(jleft,1,MPI_INTEGER,myright,0,commrow,req(1),&
                      ierr)
        call MPI_IRECV(jright,1,MPI_INTEGER,myleft,0,commrow,req(2),&
                        ierr)
        call MPI_ISEND(ileft,1,MPI_INTEGER,myleft,0,commrow,req(3),&
                       ierr)
        call MPI_ISEND(iright,1,MPI_INTEGER,myright,0,commrow,req(4),&
                       ierr)
        call MPI_WAITALL(4,req,statarry,ierr) 

!        call MPI_BARRIER(commrow,ierr)
!        if(myid.eq.0) then
!          print*,"pass 1:"
!        endif

        jup = 0
        jdown = 0

        call MPI_IRECV(jdown,1,MPI_INTEGER,myup,0,commcol,req(1),&
                      ierr)
        call MPI_IRECV(jup,1,MPI_INTEGER,mydown,0,commcol,req(2),&
                      ierr)
        call MPI_ISEND(idown,1,MPI_INTEGER,mydown,0,commcol,req(3),&
                       ierr)
        call MPI_ISEND(iup,1,MPI_INTEGER,myup,0,commcol,req(4),&
                       ierr)
        call MPI_WAITALL(4,req,statarry,ierr) 

        numbuf = jleft+jright+jup+jdown 
        
!        call MPI_BARRIER(commcol,ierr)
!        if(myid.eq.0) then
!          print*,"pass 2:"
!        endif

        deallocate(recv)
        allocate(recv(9,numbuf+nmv0+1))
        do i = 1, nmv0
          recv(:,i) = temp1(:,i)
        enddo
        deallocate(temp1)

        !send outgoing particles to left neibhoring processor.
        !allocate(temp1(9,jleft))
        ! temp1 = 0.0
        jleft = 9*jleft
        ileft = 9*ileft
        !call MPI_IRECV(temp1(1,1),jleft,MPI_DOUBLE_PRECISION,myright,&
        !               0,commrow,msid,ierr)
        call MPI_IRECV(recv(1,nmv0+1),jleft,MPI_DOUBLE_PRECISION,myright,&
                       0,commrow,msid,ierr)
        call MPI_SEND(left(1,1),ileft,MPI_DOUBLE_PRECISION,myleft,&
                      0,commrow,ierr)
        call MPI_WAIT(msid,status,ierr) 
        ileft = ileft/9
        jleft = jleft/9
        !do i = 1, jleft
        !  recv(:,i+nmv0) = temp1(:,i)
        !enddo
        nmv0 = nmv0+jleft
        !deallocate(temp1)
        
        !send outgoing particles to right neibhoring processor.
        !allocate(temp1(9,jright))
        !temp1 = 0.0
        jright = 9*jright
        iright = 9*iright
        !call MPI_IRECV(temp1(1,1),jright,MPI_DOUBLE_PRECISION,myleft,&
        !                0,commrow,msid,ierr)
        call MPI_IRECV(recv(1,nmv0+1),jright,MPI_DOUBLE_PRECISION,myleft,&
                        0,commrow,msid,ierr)
        call MPI_SEND(right(1,1),iright,MPI_DOUBLE_PRECISION,myright,&
                      0,commrow,ierr)
        call MPI_WAIT(msid,status,ierr) 
        iright = iright/9
        jright = jright/9
        !do i = 1, jright
        !  recv(:,i+nmv0) = temp1(:,i)
        !enddo
        nmv0 = nmv0 + jright
        !deallocate(temp1)

!        call MPI_BARRIER(commrow,ierr)
!        if(myid.eq.0) then
!          print*,"pass 3:"
!        endif

        !send outgoing particles to down neibhoring processor.
        !allocate(temp1(9,jdown))
        !temp1 = 0.0
        jdown = 9*jdown
        idown = 9*idown
        !call MPI_IRECV(temp1(1,1),jdown,MPI_DOUBLE_PRECISION,myup,&
        !                0,commcol,msid,ierr)
        call MPI_IRECV(recv(1,nmv0+1),jdown,MPI_DOUBLE_PRECISION,myup,&
                        0,commcol,msid,ierr)
        call MPI_SEND(down(1,1),idown,MPI_DOUBLE_PRECISION,mydown,&
                        0,commcol,ierr)
        call MPI_WAIT(msid,status,ierr)
        idown = idown/9
        jdown = jdown/9
        !do i = 1, jdown
        !  recv(:,i+nmv0) = temp1(:,i)
        !enddo
        nmv0 = nmv0 + jdown
        !deallocate(temp1)

        !send outgoing particles to up neibhoring processor.
        !allocate(temp1(9,jup))
        !temp1 = 0.0
        jup = 9*jup
        iup = 9*iup
        !call MPI_IRECV(temp1(1,1),jup,MPI_DOUBLE_PRECISION,mydown,&
        !              0,commcol,msid,ierr)
        call MPI_IRECV(recv(1,nmv0+1),jup,MPI_DOUBLE_PRECISION,mydown,&
                      0,commcol,msid,ierr)
        call MPI_SEND(up(1,1),iup,MPI_DOUBLE_PRECISION,myup,&
                      0,commcol,ierr)
        call MPI_WAIT(msid,status,ierr)
        iup = iup/9
        jup = jup/9
        !do i = 1, jup
        !  recv(:,i+nmv0) = temp1(:,i)
        !enddo
        nmv0 = nmv0 + jup
        !deallocate(temp1)
 
!        call MPI_BARRIER(commcol,ierr)
!        if(myid.eq.0) then
!          print*,"pass 4:"
!        endif

        deallocate(up)
        deallocate(down)
        deallocate(left)
        deallocate(right)

        if(numbuf.gt.480) then
          nsmall = numbuf/16
        else
          nsmall = numbuf
        endif
        allocate(left(9,nsmall))
        allocate(right(9,nsmall))
        allocate(up(9,nsmall))
        allocate(down(9,nsmall))
        allocate(mmsk(numbuf))
        ileft = 0
        iright = 0
        iup = 0
        idown = 0
        nmv0 = nmv0 - numbuf
        do i = 1, numbuf
          mmsk(i) = .true.
          ii = i+nmv0
          if(recv(5,ii).le.lcrange(5)) then
            if(myidx.ne.0) then
              ileft = ileft + 1
              left(:,ileft) = recv(:,ii)
              mmsk(i) = .false.
              if(mod(ileft,nsmall).eq.0) then
                allocate(temp1(9,ileft))
                do j = 1,ileft
                  temp1(:,j) = left(:,j)
                enddo
                deallocate(left)
                allocate(left(9,ileft+nsmall))
                do j = 1,ileft
                  left(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            else
              if(recv(3,ii).gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = recv(:,ii)
                mmsk(i) = .false.
                if(mod(iup,nsmall).eq.0) then
                  allocate(temp1(9,iup))
                  do j = 1, iup
                    temp1(:,j) = up(:,j)
                  enddo
                  deallocate(up)
                  allocate(up(9,iup+nsmall))
                  do j = 1, iup
                    up(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else if(recv(3,ii).le.lcrange(3)) then
                if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = recv(:,ii)
                mmsk(i) = .false.
                if(mod(idown,nsmall).eq.0) then
                  allocate(temp1(9,idown))
                  do j = 1, idown
                    temp1(:,j) = down(:,j)
                  enddo
                  deallocate(down)
                  allocate(down(9,idown+nsmall))
                  do j = 1, idown
                    down(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else
              endif
            endif
          else if(recv(5,ii).gt.lcrange(6)) then
            if(myidx.ne.(npx-1)) then
              iright = iright + 1
              right(:,iright) = recv(:,ii)
              mmsk(i) = .false.
              if(mod(iright,nsmall).eq.0) then
                allocate(temp1(9,iright))
                do j =1, iright
                  temp1(:,j) = right(:,j)
                enddo
                deallocate(right)
                allocate(right(9,iright+nsmall))
                do j =1, iright
                  right(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            else
              if(recv(3,ii).gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = recv(:,ii)
                mmsk(i) = .false.
                if(mod(iup,nsmall).eq.0) then
                  allocate(temp1(9,iup))
                  do j = 1, iup
                    temp1(:,j) = up(:,j)
                  enddo
                  deallocate(up)
                  allocate(up(9,iup+nsmall))
                  do j = 1, iup
                    up(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else if(recv(3,ii).le.lcrange(3)) then
                if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = recv(:,ii)
                mmsk(i) = .false.
                if(mod(idown,nsmall).eq.0) then
                  allocate(temp1(9,idown))
                  do j = 1, idown
                    temp1(:,j) = down(:,j)
                  enddo
                  deallocate(down)
                  allocate(down(9,idown+nsmall))
                  do j = 1, idown
                    down(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else
              endif
            endif
          else if(recv(3,ii).gt.lcrange(4)) then
            if(myidy.ne.(npy-1)) then
              iup = iup + 1
              up(:,iup) = recv(:,ii)
              mmsk(i) = .false.
              if(mod(iup,nsmall).eq.0) then
                allocate(temp1(9,iup))
                do j = 1, iup
                  temp1(:,j) = up(:,j)
                enddo
                deallocate(up)
                allocate(up(9,iup+nsmall))
                do j = 1, iup
                  up(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            endif
          else if(recv(3,ii).le.lcrange(3)) then
            if(myidy.ne.0) then
              idown = idown + 1
              down(:,idown) = recv(:,ii)
              mmsk(i) = .false.
              if(mod(idown,nsmall).eq.0) then
                allocate(temp1(9,idown))
                do j = 1, idown
                  temp1(:,j) = down(:,j)
                enddo
                deallocate(down)
                allocate(down(9,idown+nsmall))
                do j = 1, idown
                  down(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            endif
          else
          endif
        enddo
        nmv = ileft+iright+idown+iup
        call MPI_ALLREDUCE(nmv,totnmv,1,MPI_INTEGER,MPI_SUM, &
                           comm2d,ierr)
        if(totnmv.eq.0) then
          deallocate(mmsk)
          deallocate(up)
          deallocate(down)
          deallocate(left)
          deallocate(right)
          nmv0 = nmv0 + numbuf - nmv
          exit
        endif

        ic = 0
        allocate(temp1(9,nmv0+numbuf-nmv))
        do i = 1, nmv0
          temp1(:,i) = recv(:,i)
        enddo
        do i = 1, numbuf
          ii = i + nmv0
          if(mmsk(i)) then
            ic = ic + 1
            temp1(:,ic+nmv0) = recv(:,ii)
          endif
        enddo
        deallocate(mmsk)
        nmv0 = nmv0 + numbuf - nmv

!        print*," loop ", ij,numbuf,nmv,nmv0
!        call MPI_BARRIER(comm2d,ierr)

        enddo

        !copy the remaining local particles into a temporary array.
        numpts = Nptlocal-nout
        allocate(temp1(9,numpts))
        ic = 0
        do i = 1, Nptlocal
          if(msk(i)) then
            ic = ic + 1
            do j = 1, 9
              temp1(j,ic) = Ptsl(j,i)
            enddo
          endif
        enddo

        deallocate(msk)

!        call MPI_BARRIER(comm2d,ierr)
        !recopy the remaining local particles back to Ptsl which has 
        !a new size now.
        Nptlocal = numpts+nmv0 
        deallocate(Ptsl)
        allocate(Ptsl(9,Nptlocal))
        do i = 1, numpts
          do j = 1, 9
            Ptsl(j,i) = temp1(j,i)
          enddo
        enddo
        deallocate(temp1)
        do i = 1, nmv0
          ii = i + numpts
          do j = 1, 9
            Ptsl(j,ii) = recv(j,i)
          enddo
        enddo

        deallocate(recv)

        call MPI_BARRIER(comm2d,ierr)
        t_ptsmv = t_ptsmv + elapsedtime_Timer(t0)

        end subroutine ptsmv_ptclmger

        subroutine ptsmv1_ptclmger(Ptsl,Nptlocal,grid,pdim,npmax,lcrange)
        implicit none
        include 'mpif.h'
        type (Pgrid2d), intent(in) :: grid
        integer, intent(inout) :: Nptlocal
        integer, intent(in) :: pdim,npmax
        double precision, pointer, dimension(:,:) :: Ptsl
!        double precision, dimension(pdim,npmax) :: Ptsl
        double precision, dimension(:),intent(in) :: lcrange
        double precision, allocatable, dimension(:,:) :: left,right,up,down
        double precision, allocatable, dimension(:,:) :: temp1,recv
        integer :: myid,myidx,myidy,totnp,npy,npx, &
                   comm2d,commcol,commrow
        integer :: ileft,iright,iup,idown,iupright,iupleft,&
                   idownleft,idownright
        integer :: jleft,jright,jup,jdown,jupright,jupleft,&
                   jdownleft,jdownright
        integer :: myleft,myright,myup,mydown,myupright,myupleft,&
                   mydwnleft,mydwnright
        integer status(MPI_STATUS_SIZE)
        integer :: nsmall,i,j,numpts,ic,ierr
        integer, dimension(2) :: tmpcoord
        logical, allocatable, dimension(:) :: msk,mmsk
        double precision :: t0,t1
        integer :: numbuf,nmv,nmv0,nout,ii,totnmv,ij

        call starttime_Timer(t0)
!        call starttime_Timer(t1)

        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        if(Nptlocal.gt.480) then
          nsmall = Nptlocal/16
        else
          nsmall = Nptlocal
        endif
        allocate(left(9,nsmall))
        allocate(right(9,nsmall))
        allocate(up(9,nsmall))
        allocate(down(9,nsmall))

        left = 0.0
        right = 0.0
        up = 0.0
        down = 0.0
        ileft = 0
        iright = 0
        iup = 0
        idown = 0

        allocate(msk(Nptlocal))
        msk = .true.

        do i = 1, Nptlocal
          if(Ptsl(5,i).le.lcrange(5)) then
            if(myidx.ne.0) then
              ileft = ileft + 1
              left(:,ileft) = Ptsl(:,i)
              msk(i) = .false.
              if(mod(ileft,nsmall).eq.0) then
                allocate(temp1(9,ileft))
                do j = 1,ileft
                  temp1(:,j) = left(:,j)
                enddo
                deallocate(left)
                allocate(left(9,ileft+nsmall))
                do j = 1,ileft
                  left(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            else
              if(Ptsl(3,i).gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = Ptsl(:,i)
                msk(i) = .false.
                if(mod(iup,nsmall).eq.0) then
                  allocate(temp1(9,iup))
                  do j = 1, iup
                    temp1(:,j) = up(:,j)
                  enddo
                  deallocate(up)
                  allocate(up(9,iup+nsmall))
                  do j = 1, iup
                    up(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else if(Ptsl(3,i).le.lcrange(3)) then
                if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = Ptsl(:,i)
                msk(i) = .false.
                if(mod(idown,nsmall).eq.0) then
                  allocate(temp1(9,idown))
                  do j = 1, idown
                    temp1(:,j) = down(:,j)
                  enddo
                  deallocate(down)
                  allocate(down(9,idown+nsmall))
                  do j = 1, idown
                    down(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else
              endif
            endif
          else if(Ptsl(5,i).gt.lcrange(6)) then
            if(myidx.ne.(npx-1)) then
              iright = iright + 1
              right(:,iright) = Ptsl(:,i)
              msk(i) = .false.
              if(mod(iright,nsmall).eq.0) then
                allocate(temp1(9,iright))
                do j =1, iright
                  temp1(:,j) = right(:,j)
                enddo
                deallocate(right)
                allocate(right(9,iright+nsmall))
                do j =1, iright
                  right(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            else
              if(Ptsl(3,i).gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = Ptsl(:,i)
                msk(i) = .false.
                if(mod(iup,nsmall).eq.0) then
                  allocate(temp1(9,iup))
                  do j = 1, iup
                    temp1(:,j) = up(:,j)
                  enddo
                  deallocate(up)
                  allocate(up(9,iup+nsmall))
                  do j = 1, iup
                    up(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else if(Ptsl(3,i).le.lcrange(3)) then
                if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = Ptsl(:,i)
                msk(i) = .false.
                if(mod(idown,nsmall).eq.0) then
                  allocate(temp1(9,idown))
                  do j = 1, idown
                    temp1(:,j) = down(:,j)
                  enddo
                  deallocate(down)
                  allocate(down(9,idown+nsmall))
                  do j = 1, idown
                    down(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else
              endif
            endif
          else if(Ptsl(3,i).gt.lcrange(4)) then
            if(myidy.ne.(npy-1)) then
              iup = iup + 1
              up(:,iup) = Ptsl(:,i)
              msk(i) = .false.
              if(mod(iup,nsmall).eq.0) then
                allocate(temp1(9,iup))
                do j = 1, iup
                  temp1(:,j) = up(:,j)
                enddo
                deallocate(up)
                allocate(up(9,iup+nsmall))
                do j = 1, iup
                  up(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            endif
          else if(Ptsl(3,i).le.lcrange(3)) then
            if(myidy.ne.0) then
              idown = idown + 1
              down(:,idown) = Ptsl(:,i)
              msk(i) = .false.
              if(mod(idown,nsmall).eq.0) then
                allocate(temp1(9,idown))
                do j = 1, idown
                  temp1(:,j) = down(:,j)
                enddo
                deallocate(down)
                allocate(down(9,idown+nsmall))
                do j = 1, idown
                  down(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            endif
          else
          endif

        enddo

!        t_ptsmv1 = t_ptsmv1 + elapsedtime_Timer(t1)
        nmv0 = 0
        nout = ileft+iright+iup+idown
!        print*,"in comm: ",Nptlocal,nout
        allocate(recv(9,nmv0+1))
        allocate(temp1(9,nmv0+1))
        ij = 0

        do

        ij = ij + 1
        myleft = myidx - 1
        myright = myidx + 1
        jleft = 0
        jright = 0

        if(myidx.ne.0) then
          call MPI_SEND(ileft,1,MPI_INTEGER,myleft,0,commrow,ierr)
        endif
        if(myidx.ne.(npx-1)) then
          call MPI_RECV(jleft,1,MPI_INTEGER,myright,0,commrow,status,&
                        ierr)
        endif

        if(myidx.ne.(npx-1)) then
          call MPI_SEND(iright,1,MPI_INTEGER,myright,0,commrow,ierr)
        endif
        if(myidx.ne.0) then
          call MPI_RECV(jright,1,MPI_INTEGER,myleft,0,commrow,status,&
                        ierr)
        endif

        myup = myidy + 1
        mydown = myidy - 1
        jup = 0
        jdown = 0

        if(myidy.ne.0) then
          call MPI_SEND(idown,1,MPI_INTEGER,mydown,0,commcol,ierr)
        endif
        if(myidy.ne.(npy-1)) then
          call MPI_RECV(jdown,1,MPI_INTEGER,myup,0,commcol,status,&
                        ierr)
        endif
        
        if(myidy.ne.(npy-1)) then
          call MPI_SEND(iup,1,MPI_INTEGER,myup,0,commcol,ierr)
        endif
        if(myidy.ne.0) then
          call MPI_RECV(jup,1,MPI_INTEGER,mydown,0,commcol,status,&
                        ierr)
        endif

        numbuf = jleft+jright+jup+jdown 
        
        deallocate(recv)
        allocate(recv(9,numbuf+nmv0+1))
        do i = 1, nmv0
          recv(:,i) = temp1(:,i)
        enddo
        deallocate(temp1)

        !send outgoing particles to left neibhoring processor.
        allocate(temp1(9,jleft))
        temp1 = 0.0
        if(myidx.ne.0) then
          ileft = 9*ileft
          call MPI_SEND(left(1,1),ileft,MPI_DOUBLE_PRECISION,myleft,&
                        0,commrow,ierr)
          ileft = ileft/9
        endif
        if(myidx.ne.(npx-1)) then
          jleft = 9*jleft
          call MPI_RECV(temp1(1,1),jleft,MPI_DOUBLE_PRECISION,myright,&
                        0,commrow,&
                        status,ierr)
          jleft = jleft/9
        endif
        do i = 1, jleft
          recv(:,i+nmv0) = temp1(:,i)
        enddo
        nmv0 = nmv0+jleft
        deallocate(temp1)
        
!        call MPI_BARRIER(comm2d,ierr)
        !send outgoing particles to right neibhoring processor.
        allocate(temp1(9,jright))
        temp1 = 0.0
        if(myidx.ne.(npx-1)) then
          iright = 9*iright
          call MPI_SEND(right(1,1),iright,MPI_DOUBLE_PRECISION,myright,&
                        0,commrow,&
                        ierr)
          iright = iright/9
        endif
        if(myidx.ne.0) then
          jright = 9*jright
          call MPI_RECV(temp1(1,1),jright,MPI_DOUBLE_PRECISION,myleft,&
                        0,commrow,&
                        status,ierr)
          jright = jright/9
        endif
        do i = 1, jright
          recv(:,i+nmv0) = temp1(:,i)
        enddo
        nmv0 = nmv0 + jright
        deallocate(temp1)

!        call MPI_BARRIER(comm2d,ierr)
        !send outgoing particles to down neibhoring processor.
        allocate(temp1(9,jdown))
        temp1 = 0.0
        if(myidy.ne.0) then
          idown = 9*idown
          call MPI_SEND(down(1,1),idown,MPI_DOUBLE_PRECISION,mydown,&
                        0,commcol,&
                        ierr)
          idown = idown/9
        endif
        if(myidy.ne.(npy-1)) then
          jdown = 9*jdown
          call MPI_RECV(temp1(1,1),jdown,MPI_DOUBLE_PRECISION,myup,&
                        0,commcol,&
                        status,ierr)
          jdown = jdown/9
        endif
        do i = 1, jdown
          recv(:,i+nmv0) = temp1(:,i)
        enddo
        nmv0 = nmv0 + jdown
        deallocate(temp1)

!        call MPI_BARRIER(comm2d,ierr)
        !send outgoing particles to up neibhoring processor.
        allocate(temp1(9,jup))
        temp1 = 0.0
        if(myidy.ne.(npy-1)) then
          iup = 9*iup
          call MPI_SEND(up(1,1),iup,MPI_DOUBLE_PRECISION,myup,&
                        0,commcol,&
                        ierr)
          iup = iup/9
        endif
        if(myidy.ne.0) then
          jup = 9*jup
          call MPI_RECV(temp1(1,1),jup,MPI_DOUBLE_PRECISION,mydown,&
                        0,commcol,&
                        status,ierr)
          jup = jup/9
        endif
        do i = 1, jup
          recv(:,i+nmv0) = temp1(:,i)
        enddo
        nmv0 = nmv0 + jup
        deallocate(temp1)
 
        deallocate(up)
        deallocate(down)
        deallocate(left)
        deallocate(right)

        if(numbuf.gt.480) then
          nsmall = numbuf/16
        else
          nsmall = numbuf
        endif
        allocate(left(9,nsmall))
        allocate(right(9,nsmall))
        allocate(up(9,nsmall))
        allocate(down(9,nsmall))
        allocate(mmsk(numbuf))
        mmsk = .true.
        ileft = 0
        iright = 0
        iup = 0
        idown = 0
        nmv0 = nmv0 - numbuf
        do i = 1, numbuf
          ii = i+nmv0
          if(recv(5,ii).le.lcrange(5)) then
            if(myidx.ne.0) then
              ileft = ileft + 1
              left(:,ileft) = recv(:,ii)
              mmsk(i) = .false.
              if(mod(ileft,nsmall).eq.0) then
                allocate(temp1(9,ileft))
                do j = 1,ileft
                  temp1(:,j) = left(:,j)
                enddo
                deallocate(left)
                allocate(left(9,ileft+nsmall))
                do j = 1,ileft
                  left(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            else
              if(recv(3,ii).gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = recv(:,ii)
                mmsk(i) = .false.
                if(mod(iup,nsmall).eq.0) then
                  allocate(temp1(9,iup))
                  do j = 1, iup
                    temp1(:,j) = up(:,j)
                  enddo
                  deallocate(up)
                  allocate(up(9,iup+nsmall))
                  do j = 1, iup
                    up(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else if(recv(3,ii).le.lcrange(3)) then
                if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = recv(:,ii)
                mmsk(i) = .false.
                if(mod(idown,nsmall).eq.0) then
                  allocate(temp1(9,idown))
                  do j = 1, idown
                    temp1(:,j) = down(:,j)
                  enddo
                  deallocate(down)
                  allocate(down(9,idown+nsmall))
                  do j = 1, idown
                    down(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else
              endif
            endif
          else if(recv(5,ii).gt.lcrange(6)) then
            if(myidx.ne.(npx-1)) then
              iright = iright + 1
              right(:,iright) = recv(:,ii)
              mmsk(i) = .false.
              if(mod(iright,nsmall).eq.0) then
                allocate(temp1(9,iright))
                do j =1, iright
                  temp1(:,j) = right(:,j)
                enddo
                deallocate(right)
                allocate(right(9,iright+nsmall))
                do j =1, iright
                  right(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            else
              if(recv(3,ii).gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = recv(:,ii)
                mmsk(i) = .false.
                if(mod(iup,nsmall).eq.0) then
                  allocate(temp1(9,iup))
                  do j = 1, iup
                    temp1(:,j) = up(:,j)
                  enddo
                  deallocate(up)
                  allocate(up(9,iup+nsmall))
                  do j = 1, iup
                    up(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else if(recv(3,ii).le.lcrange(3)) then
                if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = recv(:,ii)
                mmsk(i) = .false.
                if(mod(idown,nsmall).eq.0) then
                  allocate(temp1(9,idown))
                  do j = 1, idown
                    temp1(:,j) = down(:,j)
                  enddo
                  deallocate(down)
                  allocate(down(9,idown+nsmall))
                  do j = 1, idown
                    down(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else
              endif
            endif
          else if(recv(3,ii).gt.lcrange(4)) then
            if(myidy.ne.(npy-1)) then
              iup = iup + 1
              up(:,iup) = recv(:,ii)
              mmsk(i) = .false.
              if(mod(iup,nsmall).eq.0) then
                allocate(temp1(9,iup))
                do j = 1, iup
                  temp1(:,j) = up(:,j)
                enddo
                deallocate(up)
                allocate(up(9,iup+nsmall))
                do j = 1, iup
                  up(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            endif
          else if(recv(3,ii).le.lcrange(3)) then
            if(myidy.ne.0) then
              idown = idown + 1
              down(:,idown) = recv(:,ii)
              mmsk(i) = .false.
              if(mod(idown,nsmall).eq.0) then
                allocate(temp1(9,idown))
                do j = 1, idown
                  temp1(:,j) = down(:,j)
                enddo
                deallocate(down)
                allocate(down(9,idown+nsmall))
                do j = 1, idown
                  down(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            endif
          else
          endif
        enddo
        nmv = ileft+iright+idown+iup
        call MPI_ALLREDUCE(nmv,totnmv,1,MPI_INTEGER,MPI_SUM, &
                           comm2d,ierr)
        if(totnmv.eq.0) then
          deallocate(mmsk)
          deallocate(up)
          deallocate(down)
          deallocate(left)
          deallocate(right)
          nmv0 = nmv0 + numbuf - nmv
          exit
        endif

        ic = 0
        allocate(temp1(9,nmv0+numbuf-nmv))
        do i = 1, nmv0
          temp1(:,i) = recv(:,i)
        enddo
        do i = 1, numbuf
          ii = i + nmv0
          if(mmsk(i)) then
            ic = ic + 1
            temp1(:,ic+nmv0) = recv(:,ii)
          endif
        enddo
        deallocate(mmsk)
        nmv0 = nmv0 + numbuf - nmv

!        print*," loop ", ij,numbuf,nmv,nmv0

        enddo

        !copy the remaining local particles into a temporary array.
        numpts = Nptlocal-nout
        allocate(temp1(9,numpts))
        ic = 0
        do i = 1, Nptlocal
          if(msk(i)) then
            ic = ic + 1
            temp1(:,ic) = Ptsl(:,i)
          endif
        enddo

        deallocate(msk)

!        call MPI_BARRIER(comm2d,ierr)
        !recopy the remaining local particles back to Ptsl which has 
        !a new size now.
        Nptlocal = numpts+nmv0 
        deallocate(Ptsl)
        allocate(Ptsl(9,Nptlocal))
        do i = 1, numpts
          Ptsl(:,i) = temp1(:,i)
        enddo
        deallocate(temp1)
        do i = 1, nmv0
          ii = i + numpts
          Ptsl(:,ii) = recv(:,i)
        enddo

        deallocate(recv)

        t_ptsmv = t_ptsmv + elapsedtime_Timer(t0)

        end subroutine ptsmv1_ptclmger

        ! move particles from one processor to 4 neighboring processors.
        subroutine ptsmv2_ptclmger(Ptsl,Nptlocal,grid,pdim,npmax,lcrange)
        implicit none
        include 'mpif.h'
        type (Pgrid2d), intent(in) :: grid
        integer, intent(inout) :: Nptlocal
        integer, intent(in) :: pdim,npmax
        double precision, pointer, dimension(:,:) :: Ptsl
!        double precision, dimension(pdim,npmax) :: Ptsl
        double precision, dimension(:),intent(in) :: lcrange
        integer, parameter :: nptmv = 100000
        double precision, dimension(9,3*nptmv) :: left,right,up,down
        double precision, allocatable, dimension(:,:) :: temp1,recv
        integer :: myid,myidx,myidy,totnp,npy,npx, &
                   comm2d,commcol,commrow
        integer :: ileft,iright,iup,idown,iupright,iupleft,&
                   idownleft,idownright
        integer :: jleft,jright,jup,jdown,jupright,jupleft,&
                   jdownleft,jdownright
        integer :: myleft,myright,myup,mydown,myupright,myupleft,&
                   mydwnleft,mydwnright
        integer :: nsmall,i,j,numpts,ic
        integer, dimension(2) :: tmpcoord
        logical, dimension(Nptlocal) :: msk
        logical, allocatable, dimension(:) :: mmsk
        integer :: numbuf,nmv,nmv0,nout,ii,totnmv,ij
        integer :: msid,ierr
        integer status(MPI_STATUS_SIZE) 
        integer statarry(MPI_STATUS_SIZE,4), req(4)
        integer :: flag,Nptlocal0,nst,nout0,iileft,iiright,iiup,iidown
        integer :: totflag
        double precision :: t0

        call starttime_Timer(t0)

        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)
        if(myidx.ne.(npx-1)) then
          myright = myidx + 1
        else
          myright = MPI_PROC_NULL
        endif
        if(myidx.ne.0) then
          myleft = myidx - 1
        else
          myleft = MPI_PROC_NULL
        endif 

        if(myidy.ne.npy-1) then
          myup = myidy + 1
        else
          myup = MPI_PROC_NULL
        endif
        if(myidy.ne.0) then
          mydown = myidy -1
        else
          mydown = MPI_PROC_NULL
        endif

!        call MPI_BARRIER(comm2d,ierr)

        flag = 0
        Nptlocal0 = Nptlocal
        nout0 = 0

        do 

        ileft = 0
        iright = 0
        iup = 0
        idown = 0
        iileft = 0
        iiright = 0
        iiup = 0
        iidown = 0
        do i = 1, Nptlocal0 - nout0
          msk(i) = .true.
          if(Ptsl(5,i).le.lcrange(5)) then
            if(myidx.ne.0) then
              iileft = iileft + 1
              if(iileft.le.nptmv) then
                ileft = ileft + 1
                left(:,ileft) = Ptsl(:,i)
                msk(i) = .false.
              endif
            else
              if(Ptsl(3,i).gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                  iiup = iiup + 1
                  if(iiup.le.nptmv) then
                    iup = iup + 1
                    up(:,iup) = Ptsl(:,i)
                    msk(i) = .false.
                  endif
                endif
              else if(Ptsl(3,i).le.lcrange(3)) then
                if(myidy.ne.0) then
                  iidown = iidown + 1
                  if(iidown.le.nptmv) then
                    idown = idown + 1
                    down(:,idown) = Ptsl(:,i)
                    msk(i) = .false.
                  endif
                endif
              else
              endif
            endif
          else if(Ptsl(5,i).gt.lcrange(6)) then
            if(myidx.ne.(npx-1)) then
              iiright = iiright + 1
              if(iiright.le.nptmv) then
                iright = iright + 1
                right(:,iright) = Ptsl(:,i)
                msk(i) = .false.
              endif
            else
              if(Ptsl(3,i).gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                  iiup = iiup + 1
                  if(iiup.le.nptmv) then
                    iup = iup + 1
                    up(:,iup) = Ptsl(:,i)
                    msk(i) = .false.
                  endif
                endif
              else if(Ptsl(3,i).le.lcrange(3)) then
                if(myidy.ne.0) then
                  iidown = iidown + 1
                  if(iidown.le.nptmv) then
                    idown = idown + 1
                    down(:,idown) = Ptsl(:,i)
                    msk(i) = .false.
                  endif
                endif
              else
              endif
            endif
          else if(Ptsl(3,i).gt.lcrange(4)) then
            if(myidy.ne.(npy-1)) then
              iiup = iiup + 1
              if(iiup.le.nptmv) then
                iup = iup + 1
                up(:,iup) = Ptsl(:,i)
                msk(i) = .false.
              endif
            endif
          else if(Ptsl(3,i).le.lcrange(3)) then
            if(myidy.ne.0) then
              iidown = iidown + 1
              if(iidown.le.nptmv) then
                idown = idown + 1
                down(:,idown) = Ptsl(:,i)
                msk(i) = .false.
              endif
            endif
          else
          endif

        enddo

        if((iileft.gt.nptmv).or.(iiright.gt.nptmv).or.(iiup.gt.nptmv) &
            .or.(iidown.gt.nptmv)) then
           flag = 1
        else
           flag = 0
        endif

        nmv0 = 0
        nout = ileft+iright+iup+idown
        allocate(recv(9,nmv0+1))
        allocate(temp1(9,nmv0+1))
        ij = 0
        call MPI_BARRIER(comm2d,ierr)

        do

        ij = ij + 1
        jleft = 0
        jright = 0

        call MPI_IRECV(jleft,1,MPI_INTEGER,myright,0,commrow,req(1),&
                      ierr)
        call MPI_IRECV(jright,1,MPI_INTEGER,myleft,0,commrow,req(2),&
                        ierr)
        call MPI_ISEND(ileft,1,MPI_INTEGER,myleft,0,commrow,req(3),&
                       ierr)
        call MPI_ISEND(iright,1,MPI_INTEGER,myright,0,commrow,req(4),&
                       ierr)
        call MPI_WAITALL(4,req,statarry,ierr) 

!        call MPI_BARRIER(commrow,ierr)
!        if(myid.eq.0) then
!          print*,"pass 1:"
!        endif

        jup = 0
        jdown = 0

        call MPI_IRECV(jdown,1,MPI_INTEGER,myup,0,commcol,req(1),&
                      ierr)
        call MPI_IRECV(jup,1,MPI_INTEGER,mydown,0,commcol,req(2),&
                      ierr)
        call MPI_ISEND(idown,1,MPI_INTEGER,mydown,0,commcol,req(3),&
                       ierr)
        call MPI_ISEND(iup,1,MPI_INTEGER,myup,0,commcol,req(4),&
                       ierr)
        call MPI_WAITALL(4,req,statarry,ierr) 

        numbuf = jleft+jright+jup+jdown 
        
!        call MPI_BARRIER(commcol,ierr)
!        if(myid.eq.0) then
!          print*,"pass 2:"
!        endif

        deallocate(recv)
        allocate(recv(9,numbuf+nmv0+1))
        do i = 1, nmv0
          recv(:,i) = temp1(:,i)
        enddo
        deallocate(temp1)

        nst = nmv0 + 1
        !send outgoing particles to left neibhoring processor.
        jleft = 9*jleft
        ileft = 9*ileft
        call MPI_IRECV(recv(1,nst),jleft,MPI_DOUBLE_PRECISION,myright,&
                       0,commrow,msid,ierr)
        call MPI_SEND(left(1,1),ileft,MPI_DOUBLE_PRECISION,myleft,&
                      0,commrow,ierr)
        call MPI_WAIT(msid,status,ierr) 
        ileft = ileft/9
        jleft = jleft/9
        nmv0 = nmv0+jleft
        
        nst = nmv0 + 1
        !send outgoing particles to right neibhoring processor.
        jright = 9*jright
        iright = 9*iright
        call MPI_IRECV(recv(1,nst),jright,MPI_DOUBLE_PRECISION,myleft,&
                        0,commrow,msid,ierr)
        call MPI_SEND(right(1,1),iright,MPI_DOUBLE_PRECISION,myright,&
                      0,commrow,ierr)
        call MPI_WAIT(msid,status,ierr) 
        iright = iright/9
        jright = jright/9
        nmv0 = nmv0 + jright

!        call MPI_BARRIER(commrow,ierr)
!        if(myid.eq.0) then
!          print*,"pass 3:"
!        endif

        nst = nmv0 + 1
        !send outgoing particles to down neibhoring processor.
        jdown = 9*jdown
        idown = 9*idown
        call MPI_IRECV(recv(1,nst),jdown,MPI_DOUBLE_PRECISION,myup,&
                        0,commcol,msid,ierr)
        call MPI_SEND(down(1,1),idown,MPI_DOUBLE_PRECISION,mydown,&
                        0,commcol,ierr)
        call MPI_WAIT(msid,status,ierr)
        idown = idown/9
        jdown = jdown/9
        nmv0 = nmv0 + jdown

        nst = nmv0 + 1
        !send outgoing particles to up neibhoring processor.
        jup = 9*jup
        iup = 9*iup
        call MPI_IRECV(recv(1,nst),jup,MPI_DOUBLE_PRECISION,mydown,&
                      0,commcol,msid,ierr)
        call MPI_SEND(up(1,1),iup,MPI_DOUBLE_PRECISION,myup,&
                      0,commcol,ierr)
        call MPI_WAIT(msid,status,ierr)
        iup = iup/9
        jup = jup/9
        nmv0 = nmv0 + jup
 
!        call MPI_BARRIER(commcol,ierr)
!        if(myid.eq.0) then
!          print*,"pass 4:"
!        endif

        allocate(mmsk(numbuf))
        ileft = 0
        iright = 0
        iup = 0
        idown = 0
        nmv0 = nmv0 - numbuf
        do i = 1, numbuf
          mmsk(i) = .true.
          ii = i+nmv0
          if(recv(5,ii).le.lcrange(5)) then
            if(myidx.ne.0) then
              ileft = ileft + 1
              left(:,ileft) = recv(:,ii)
              mmsk(i) = .false.
            else
              if(recv(3,ii).gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = recv(:,ii)
                mmsk(i) = .false.
                endif
              else if(recv(3,ii).le.lcrange(3)) then
                if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = recv(:,ii)
                mmsk(i) = .false.
                endif
              else
              endif
            endif
          else if(recv(5,ii).gt.lcrange(6)) then
            if(myidx.ne.(npx-1)) then
              iright = iright + 1
              right(:,iright) = recv(:,ii)
              mmsk(i) = .false.
            else
              if(recv(3,ii).gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = recv(:,ii)
                mmsk(i) = .false.
                endif
              else if(recv(3,ii).le.lcrange(3)) then
                if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = recv(:,ii)
                mmsk(i) = .false.
                endif
              else
              endif
            endif
          else if(recv(3,ii).gt.lcrange(4)) then
            if(myidy.ne.(npy-1)) then
              iup = iup + 1
              up(:,iup) = recv(:,ii)
              mmsk(i) = .false.
            endif
          else if(recv(3,ii).le.lcrange(3)) then
            if(myidy.ne.0) then
              idown = idown + 1
              down(:,idown) = recv(:,ii)
              mmsk(i) = .false.
            endif
          else
          endif
        enddo
        nmv = ileft+iright+idown+iup
        call MPI_ALLREDUCE(nmv,totnmv,1,MPI_INTEGER,MPI_SUM, &
                           comm2d,ierr)
        if(totnmv.eq.0) then
          nmv0 = nmv0 + numbuf - nmv
          deallocate(mmsk)
          exit
        endif

        ic = 0
        allocate(temp1(9,nmv0+numbuf-nmv))
        do i = 1, nmv0
          temp1(:,i) = recv(:,i)
        enddo
        do i = 1, numbuf
          ii = i + nmv0
          if(mmsk(i)) then
            ic = ic + 1
            temp1(:,ic+nmv0) = recv(:,ii)
          endif
        enddo
        nmv0 = nmv0 + numbuf - nmv
        deallocate(mmsk)

!        call MPI_BARRIER(comm2d,ierr)

        enddo

        !copy the remaining local particles into a temporary array.
        numpts = Nptlocal-nout
        allocate(temp1(9,numpts))
        ic = 0
        do i = 1, Nptlocal0-nout0
          if(msk(i)) then
            ic = ic + 1
            do j = 1, 9
              temp1(j,ic) = Ptsl(j,i)
            enddo
          endif
        enddo
        do i = Nptlocal0-nout0+1, Nptlocal
          ii = i-nout
          do j = 1, 9
            temp1(j,ii) = Ptsl(j,i)
          enddo
        enddo
 
!        call MPI_BARRIER(comm2d,ierr)
        !recopy the remaining local particles back to Ptsl which has 
        !a new size now.
        Nptlocal = numpts+nmv0 
        deallocate(Ptsl)
        allocate(Ptsl(9,Nptlocal))
        do i = 1, numpts
          do j = 1, 9
            Ptsl(j,i) = temp1(j,i)
          enddo
        enddo
        deallocate(temp1)
        do i = 1, nmv0
          ii = i + numpts
          do j = 1, 9
            Ptsl(j,ii) = recv(j,i)
          enddo
        enddo

        deallocate(recv)

        nout0 = nout0 + nout

          call MPI_ALLREDUCE(flag,totflag,1,MPI_INTEGER,MPI_SUM, &
                           comm2d,ierr)
          if(totflag.eq.0) exit

        enddo

        t_ptsmv = t_ptsmv + elapsedtime_Timer(t0)

        end subroutine ptsmv2_ptclmger

        subroutine ptsmv2perd_ptclmger(Ptsl,Nptlocal,grid,pdim,npmax,lcrange,&
                                       perd)
        implicit none
        include 'mpif.h'
        type (Pgrid2d), intent(in) :: grid
        integer, intent(inout) :: Nptlocal
        integer, intent(in) :: pdim,npmax
        double precision, pointer, dimension(:,:) :: Ptsl
!        double precision, dimension(pdim,npmax) :: Ptsl
        double precision, dimension(:),intent(in) :: lcrange
        double precision, intent(in) :: perd
        integer, parameter :: nptmv = 10000
        double precision, dimension(9,3*nptmv) :: left,right,up,down
        double precision, allocatable, dimension(:,:) :: temp1,recv
        integer :: myid,myidx,myidy,totnp,npy,npx, &
                   comm2d,commcol,commrow
        integer :: ileft,iright,iup,idown,iupright,iupleft,&
                   idownleft,idownright
        integer :: jleft,jright,jup,jdown,jupright,jupleft,&
                   jdownleft,jdownright
        integer :: myleft,myright,myup,mydown,myupright,myupleft,&
                   mydwnleft,mydwnright
        integer :: nsmall,i,j,numpts,ic
        integer, dimension(2) :: tmpcoord
        logical, dimension(Nptlocal) :: msk
        logical, allocatable, dimension(:) :: mmsk
        integer :: numbuf,nmv,nmv0,nout,ii,totnmv,ij
        integer :: msid,ierr
        integer status(MPI_STATUS_SIZE) 
        integer statarry(MPI_STATUS_SIZE,4), req(4)
        integer :: flag,Nptlocal0,nst,nout0,iileft,iiright,iiup,iidown
        integer :: totflag
        double precision :: t0
        integer :: tmpkz
        double precision :: tmpz,halfperd

        call starttime_Timer(t0)

        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)
        if(myidx.ne.(npx-1)) then
          myright = myidx + 1
        else
          myright = MPI_PROC_NULL
        endif
        if(myidx.ne.0) then
          myleft = myidx - 1
        else
          myleft = MPI_PROC_NULL
        endif 

        if(myidy.ne.npy-1) then
          myup = myidy + 1
        else
          myup = MPI_PROC_NULL
        endif
        if(myidy.ne.0) then
          mydown = myidy -1
        else
          mydown = MPI_PROC_NULL
        endif

!        call MPI_BARRIER(comm2d,ierr)

        halfperd = perd/2
        flag = 0
        Nptlocal0 = Nptlocal
        nout0 = 0

        do 

        ileft = 0
        iright = 0
        iup = 0
        idown = 0
        iileft = 0
        iiright = 0
        iiup = 0
        iidown = 0
        do i = 1, Nptlocal0 - nout0
          msk(i) = .true.
          tmpkz = Ptsl(5,i)/halfperd
          tmpz = mod(Ptsl(5,i),halfperd) - mod(tmpkz,2)*halfperd
!          if(Ptsl(5,i).le.lcrange(5)) then
          if(tmpz.le.lcrange(5)) then
            if(myidx.ne.0) then
              iileft = iileft + 1
              if(iileft.le.nptmv) then
                ileft = ileft + 1
                left(:,ileft) = Ptsl(:,i)
                msk(i) = .false.
              endif
            else
              if(Ptsl(3,i).gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                  iiup = iiup + 1
                  if(iiup.le.nptmv) then
                    iup = iup + 1
                    up(:,iup) = Ptsl(:,i)
                    msk(i) = .false.
                  endif
                endif
              else if(Ptsl(3,i).le.lcrange(3)) then
                if(myidy.ne.0) then
                  iidown = iidown + 1
                  if(iidown.le.nptmv) then
                    idown = idown + 1
                    down(:,idown) = Ptsl(:,i)
                    msk(i) = .false.
                  endif
                endif
              else
              endif
            endif
          !else if(Ptsl(5,i).gt.lcrange(6)) then
          else if(tmpz.gt.lcrange(6)) then
            if(myidx.ne.(npx-1)) then
              iiright = iiright + 1
              if(iiright.le.nptmv) then
                iright = iright + 1
                right(:,iright) = Ptsl(:,i)
                msk(i) = .false.
              endif
            else
              if(Ptsl(3,i).gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                  iiup = iiup + 1
                  if(iiup.le.nptmv) then
                    iup = iup + 1
                    up(:,iup) = Ptsl(:,i)
                    msk(i) = .false.
                  endif
                endif
              else if(Ptsl(3,i).le.lcrange(3)) then
                if(myidy.ne.0) then
                  iidown = iidown + 1
                  if(iidown.le.nptmv) then
                    idown = idown + 1
                    down(:,idown) = Ptsl(:,i)
                    msk(i) = .false.
                  endif
                endif
              else
              endif
            endif
          else if(Ptsl(3,i).gt.lcrange(4)) then
            if(myidy.ne.(npy-1)) then
              iiup = iiup + 1
              if(iiup.le.nptmv) then
                iup = iup + 1
                up(:,iup) = Ptsl(:,i)
                msk(i) = .false.
              endif
            endif
          else if(Ptsl(3,i).le.lcrange(3)) then
            if(myidy.ne.0) then
              iidown = iidown + 1
              if(iidown.le.nptmv) then
                idown = idown + 1
                down(:,idown) = Ptsl(:,i)
                msk(i) = .false.
              endif
            endif
          else
          endif

        enddo

        if((iileft.gt.nptmv).or.(iiright.gt.nptmv).or.(iiup.gt.nptmv) &
            .or.(iidown.gt.nptmv)) then
           flag = 1
        else
           flag = 0
        endif

        nmv0 = 0
        nout = ileft+iright+iup+idown
        allocate(recv(9,nmv0+1))
        allocate(temp1(9,nmv0+1))
        ij = 0
        call MPI_BARRIER(comm2d,ierr)

        do

        ij = ij + 1
        jleft = 0
        jright = 0

        call MPI_IRECV(jleft,1,MPI_INTEGER,myright,0,commrow,req(1),&
                      ierr)
        call MPI_IRECV(jright,1,MPI_INTEGER,myleft,0,commrow,req(2),&
                        ierr)
        call MPI_ISEND(ileft,1,MPI_INTEGER,myleft,0,commrow,req(3),&
                       ierr)
        call MPI_ISEND(iright,1,MPI_INTEGER,myright,0,commrow,req(4),&
                       ierr)
        call MPI_WAITALL(4,req,statarry,ierr) 

!        call MPI_BARRIER(commrow,ierr)
!        if(myid.eq.0) then
!          print*,"pass 1:"
!        endif

        jup = 0
        jdown = 0

        call MPI_IRECV(jdown,1,MPI_INTEGER,myup,0,commcol,req(1),&
                      ierr)
        call MPI_IRECV(jup,1,MPI_INTEGER,mydown,0,commcol,req(2),&
                      ierr)
        call MPI_ISEND(idown,1,MPI_INTEGER,mydown,0,commcol,req(3),&
                       ierr)
        call MPI_ISEND(iup,1,MPI_INTEGER,myup,0,commcol,req(4),&
                       ierr)
        call MPI_WAITALL(4,req,statarry,ierr) 

        numbuf = jleft+jright+jup+jdown 
        
!        call MPI_BARRIER(commcol,ierr)
!        if(myid.eq.0) then
!          print*,"pass 2:"
!        endif

        deallocate(recv)
        allocate(recv(9,numbuf+nmv0+1))
        do i = 1, nmv0
          recv(:,i) = temp1(:,i)
        enddo
        deallocate(temp1)

        nst = nmv0 + 1
        !send outgoing particles to left neibhoring processor.
        jleft = 9*jleft
        ileft = 9*ileft
        call MPI_IRECV(recv(1,nst),jleft,MPI_DOUBLE_PRECISION,myright,&
                       0,commrow,msid,ierr)
        call MPI_SEND(left(1,1),ileft,MPI_DOUBLE_PRECISION,myleft,&
                      0,commrow,ierr)
        call MPI_WAIT(msid,status,ierr) 
        ileft = ileft/9
        jleft = jleft/9
        nmv0 = nmv0+jleft
        
        nst = nmv0 + 1
        !send outgoing particles to right neibhoring processor.
        jright = 9*jright
        iright = 9*iright
        call MPI_IRECV(recv(1,nst),jright,MPI_DOUBLE_PRECISION,myleft,&
                        0,commrow,msid,ierr)
        call MPI_SEND(right(1,1),iright,MPI_DOUBLE_PRECISION,myright,&
                      0,commrow,ierr)
        call MPI_WAIT(msid,status,ierr) 
        iright = iright/9
        jright = jright/9
        nmv0 = nmv0 + jright

!        call MPI_BARRIER(commrow,ierr)
!        if(myid.eq.0) then
!          print*,"pass 3:"
!        endif

        nst = nmv0 + 1
        !send outgoing particles to down neibhoring processor.
        jdown = 9*jdown
        idown = 9*idown
        call MPI_IRECV(recv(1,nst),jdown,MPI_DOUBLE_PRECISION,myup,&
                        0,commcol,msid,ierr)
        call MPI_SEND(down(1,1),idown,MPI_DOUBLE_PRECISION,mydown,&
                        0,commcol,ierr)
        call MPI_WAIT(msid,status,ierr)
        idown = idown/9
        jdown = jdown/9
        nmv0 = nmv0 + jdown

        nst = nmv0 + 1
        !send outgoing particles to up neibhoring processor.
        jup = 9*jup
        iup = 9*iup
        call MPI_IRECV(recv(1,nst),jup,MPI_DOUBLE_PRECISION,mydown,&
                      0,commcol,msid,ierr)
        call MPI_SEND(up(1,1),iup,MPI_DOUBLE_PRECISION,myup,&
                      0,commcol,ierr)
        call MPI_WAIT(msid,status,ierr)
        iup = iup/9
        jup = jup/9
        nmv0 = nmv0 + jup
 
!        call MPI_BARRIER(commcol,ierr)
!        if(myid.eq.0) then
!          print*,"pass 4:"
!        endif

        allocate(mmsk(numbuf))
        ileft = 0
        iright = 0
        iup = 0
        idown = 0
        nmv0 = nmv0 - numbuf
        do i = 1, numbuf
          mmsk(i) = .true.
          ii = i+nmv0
          tmpkz = recv(5,ii)/halfperd
          tmpz = mod(recv(5,ii),halfperd) - mod(tmpkz,2)*halfperd
          !if(recv(5,ii).le.lcrange(5)) then
          if(tmpz.le.lcrange(5)) then
            if(myidx.ne.0) then
              ileft = ileft + 1
              left(:,ileft) = recv(:,ii)
              mmsk(i) = .false.
            else
              if(recv(3,ii).gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = recv(:,ii)
                mmsk(i) = .false.
                endif
              else if(recv(3,ii).le.lcrange(3)) then
                if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = recv(:,ii)
                mmsk(i) = .false.
                endif
              else
              endif
            endif
          !else if(recv(5,ii).gt.lcrange(6)) then
          else if(tmpz.gt.lcrange(6)) then
            if(myidx.ne.(npx-1)) then
              iright = iright + 1
              right(:,iright) = recv(:,ii)
              mmsk(i) = .false.
            else
              if(recv(3,ii).gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = recv(:,ii)
                mmsk(i) = .false.
                endif
              else if(recv(3,ii).le.lcrange(3)) then
                if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = recv(:,ii)
                mmsk(i) = .false.
                endif
              else
              endif
            endif
          else if(recv(3,ii).gt.lcrange(4)) then
            if(myidy.ne.(npy-1)) then
              iup = iup + 1
              up(:,iup) = recv(:,ii)
              mmsk(i) = .false.
            endif
          else if(recv(3,ii).le.lcrange(3)) then
            if(myidy.ne.0) then
              idown = idown + 1
              down(:,idown) = recv(:,ii)
              mmsk(i) = .false.
            endif
          else
          endif
        enddo
        nmv = ileft+iright+idown+iup
        call MPI_ALLREDUCE(nmv,totnmv,1,MPI_INTEGER,MPI_SUM, &
                           comm2d,ierr)
        if(totnmv.eq.0) then
          nmv0 = nmv0 + numbuf - nmv
          deallocate(mmsk)
          exit
        endif

        ic = 0
        allocate(temp1(9,nmv0+numbuf-nmv))
        do i = 1, nmv0
          temp1(:,i) = recv(:,i)
        enddo
        do i = 1, numbuf
          ii = i + nmv0
          if(mmsk(i)) then
            ic = ic + 1
            temp1(:,ic+nmv0) = recv(:,ii)
          endif
        enddo
        nmv0 = nmv0 + numbuf - nmv
        deallocate(mmsk)

!        call MPI_BARRIER(comm2d,ierr)

        enddo

        !copy the remaining local particles into a temporary array.
        numpts = Nptlocal-nout
        allocate(temp1(9,numpts))
        ic = 0
        do i = 1, Nptlocal0-nout0
          if(msk(i)) then
            ic = ic + 1
            do j = 1, 9
              temp1(j,ic) = Ptsl(j,i)
            enddo
          endif
        enddo
        do i = Nptlocal0-nout0+1, Nptlocal
          ii = i-nout
          do j = 1, 9
            temp1(j,ii) = Ptsl(j,i)
          enddo
        enddo
 
!        call MPI_BARRIER(comm2d,ierr)
        !recopy the remaining local particles back to Ptsl which has 
        !a new size now.
        Nptlocal = numpts+nmv0 
        deallocate(Ptsl)
        allocate(Ptsl(9,Nptlocal))
        do i = 1, numpts
          do j = 1, 9
            Ptsl(j,i) = temp1(j,i)
          enddo
        enddo
        deallocate(temp1)
        do i = 1, nmv0
          ii = i + numpts
          do j = 1, 9
            Ptsl(j,ii) = recv(j,i)
          enddo
        enddo

        deallocate(recv)

        nout0 = nout0 + nout

          call MPI_ALLREDUCE(flag,totflag,1,MPI_INTEGER,MPI_SUM, &
                           comm2d,ierr)
          if(totflag.eq.0) exit

        enddo

        t_ptsmv = t_ptsmv + elapsedtime_Timer(t0)

        end subroutine ptsmv2perd_ptclmger


        ! move particles from one processor to 4 neighboring processors.
        subroutine ptsmv3_ptclmger(Ptsl,Nptlocal,grid,pdim,npmax,lcrange)
        implicit none
        include 'mpif.h'
        type (Pgrid2d), intent(in) :: grid
        integer, intent(inout) :: Nptlocal
        integer, intent(in) :: pdim,npmax
        double precision, pointer, dimension(:,:) :: Ptsl
!        double precision, dimension(pdim,npmax) :: Ptsl
        double precision, dimension(:),intent(in) :: lcrange
        integer, parameter :: nptmv = 50000
        double precision, dimension(9,nptmv) :: sender
        double precision, allocatable, dimension(:,:) :: temp1,recv
        integer :: myid,myidx,myidy,totnp,npy,npx, &
                   comm2d,commcol,commrow
        integer :: myleft,myright,myup,mydown
        integer :: ileft,iright,iup,idown,jleft,jright,jup,jdown
        integer :: flag,nst,nout0,iileft,iiright,iiup,iidown
        integer :: i,j,k,numpts,ic,msid,ierr
        integer :: numbuf,nmv,nmv0,nout,ii,totnmv,ij,iout,jin
        integer, allocatable, dimension(:) :: ms1,ms2
        integer status(MPI_STATUS_SIZE) 
        integer statarry(MPI_STATUS_SIZE,2), req(2)
        logical, allocatable, dimension(:) :: mmsk,msk
        double precision :: t0

        call starttime_Timer(t0)

        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)
        if(myidx.ne.(npx-1)) then
          myright = myidx + 1
        else
          myright = MPI_PROC_NULL
        endif
        if(myidx.ne.0) then
          myleft = myidx - 1
        else
          myleft = MPI_PROC_NULL
        endif 

        if(myidy.ne.npy-1) then
          myup = myidy + 1
        else
          myup = MPI_PROC_NULL
        endif
        if(myidy.ne.0) then
          mydown = myidy -1
        else
          mydown = MPI_PROC_NULL
        endif

!        call MPI_BARRIER(comm2d,ierr)

        iileft = 0
        iiright = 0
        allocate(msk(Nptlocal))
        allocate(ms1(Nptlocal))
        allocate(ms2(Nptlocal))
        do i = 1, Nptlocal
          msk(i) = .true.
          if(Ptsl(5,i).le.lcrange(5)) then
            if(myidx.ne.0) then
              iileft = iileft + 1
              ms1(iileft) = i
              msk(i) = .false.
            endif
          else if(Ptsl(5,i).gt.lcrange(6)) then
            if(myidx.ne.(npx-1)) then
              iiright = iiright + 1
              ms2(iiright) = i
              msk(i) = .false.
            endif
          else
          endif
        enddo

!        call MPI_BARRIER(comm2d,ierr)
        nout0 = 0
        nmv0 = 0
        allocate(recv(9,nmv0+1))
        do

        if(iileft.gt.nptmv) then
           flag = 1
           iout = nptmv
        else
           flag = 0
           iout = iileft
        endif
        nout = iout
        do i = 1, iout
          k = ms1(i+nout0)
          do j = 1, 9
            sender(j,i) = Ptsl(j,k)
          enddo
        enddo

        allocate(temp1(9,nmv0+1))
        do i = 1, nmv0
          temp1(:,i) = recv(:,i)
        enddo

        do
          jin = 0
          call MPI_IRECV(jin,1,MPI_INTEGER,myright,0,commrow,req(1),&
                        ierr)
          call MPI_ISEND(iout,1,MPI_INTEGER,myleft,0,commrow,req(2),&
                         ierr)
          call MPI_WAITALL(2,req,statarry,ierr) 

          numbuf = jin
          deallocate(recv)
          allocate(recv(9,numbuf+nmv0+1))
          do i = 1, nmv0
            recv(:,i) = temp1(:,i)
          enddo
          deallocate(temp1)

          nst = nmv0 + 1
          !send outgoing particles to left neibhoring processor.
          jin = 9*jin
          iout = 9*iout
          call MPI_IRECV(recv(1,nst),jin,MPI_DOUBLE_PRECISION,myright,&
                         0,commrow,msid,ierr)
          call MPI_SEND(sender(1,1),iout,MPI_DOUBLE_PRECISION,myleft,&
                        0,commrow,ierr)
          call MPI_WAIT(msid,status,ierr) 
          iout = iout/9
          jin = jin/9
          nmv0 = nmv0+jin
        
          allocate(mmsk(numbuf))
          iout = 0
          nmv0 = nmv0 - numbuf
          do i = 1, numbuf
            mmsk(i) = .true.
            ii = i+nmv0
            if(recv(5,ii).le.lcrange(5)) then
              if(myidx.ne.0) then
                iout = iout + 1
                sender(:,iout) = recv(:,ii)
                mmsk(i) = .false.
              endif
            endif
          enddo
          nmv = iout
          call MPI_ALLREDUCE(nmv,totnmv,1,MPI_INTEGER,MPI_SUM, &
                             comm2d,ierr)
          if(totnmv.eq.0) then
            nmv0 = nmv0 + numbuf - nmv
            deallocate(mmsk)
            exit
          endif

          ic = 0
          allocate(temp1(9,nmv0+numbuf-nmv))
          do i = 1, nmv0
            temp1(:,i) = recv(:,i)
          enddo
          do i = 1, numbuf
            ii = i + nmv0
            if(mmsk(i)) then
              ic = ic + 1
              temp1(:,ic+nmv0) = recv(:,ii)
            endif
          enddo
          nmv0 = nmv0 + numbuf - nmv
          deallocate(mmsk)
        enddo

        if(flag.eq.1) then
          nout0 = nout0 + nout
          iileft = iileft - nptmv
        else
          exit
        endif

        enddo

        ! send the particles moving to the right neighoring PE.
        nout0 = 0
        do

        if(iiright.gt.nptmv) then
           flag = 1
           iout = nptmv
        else
           flag = 0
           iout = iiright
        endif
        nout = iout
        do i = 1, iout
          k = ms2(i+nout0)
          do j = 1, 9
            sender(j,i) = Ptsl(j,k)
          enddo
        enddo

        allocate(temp1(9,nmv0+1))
        do i = 1, nmv0
          temp1(:,i) = recv(:,i)
        enddo

        do
          jin = 0
          call MPI_IRECV(jin,1,MPI_INTEGER,myleft,0,commrow,req(1),&
                        ierr)
          call MPI_ISEND(iout,1,MPI_INTEGER,myright,0,commrow,req(2),&
                         ierr)
          call MPI_WAITALL(2,req,statarry,ierr) 

          numbuf = jin
          deallocate(recv)
          allocate(recv(9,numbuf+nmv0+1))
          do i = 1, nmv0
            recv(:,i) = temp1(:,i)
          enddo
          deallocate(temp1)

          nst = nmv0 + 1
          !send outgoing particles to left neibhoring processor.
          jin = 9*jin
          iout = 9*iout
          call MPI_IRECV(recv(1,nst),jin,MPI_DOUBLE_PRECISION,myleft,&
                         0,commrow,msid,ierr)
          call MPI_SEND(sender(1,1),iout,MPI_DOUBLE_PRECISION,myright,&
                        0,commrow,ierr)
          call MPI_WAIT(msid,status,ierr) 
          iout = iout/9
          jin = jin/9
          nmv0 = nmv0+jin
        
          allocate(mmsk(numbuf))
          iout = 0
          nmv0 = nmv0 - numbuf
          do i = 1, numbuf
            mmsk(i) = .true.
            ii = i+nmv0
            if(recv(5,ii).gt.lcrange(6)) then
              if(myidx.ne.0) then
                iout = iout + 1
                sender(:,iout) = recv(:,ii)
                mmsk(i) = .false.
              endif
            endif
          enddo
          nmv = iout
          call MPI_ALLREDUCE(nmv,totnmv,1,MPI_INTEGER,MPI_SUM, &
                             comm2d,ierr)
          if(totnmv.eq.0) then
            nmv0 = nmv0 + numbuf - nmv
            deallocate(mmsk)
            exit
          endif

          ic = 0
          allocate(temp1(9,nmv0+numbuf-nmv))
          do i = 1, nmv0
            temp1(:,i) = recv(:,i)
          enddo
          do i = 1, numbuf
            ii = i + nmv0
            if(mmsk(i)) then
              ic = ic + 1
              temp1(:,ic+nmv0) = recv(:,ii)
            endif
          enddo
          nmv0 = nmv0 + numbuf - nmv
          deallocate(mmsk)
        enddo

        if(flag.eq.1) then
          nout0 = nout0 + nout
          iiright = iiright - nptmv
        else
          exit
        endif

        enddo

        numpts = Nptlocal - iileft - iiright
        allocate(temp1(9,numpts))
        ic = 0
        do i = 1, Nptlocal
          if(msk(i)) then
            ic = ic + 1
            do j = 1, 9
              temp1(j,ic) = Ptsl(j,i)
            enddo
          endif
        enddo
        deallocate(msk)
        deallocate(Ptsl)
        Nptlocal = numpts + nmv0
        allocate(Ptsl(9,Nptlocal))
        do i = 1, numpts
          do j = 1, 9
            Ptsl(j,i) = temp1(j,i)
          enddo
        enddo
        deallocate(temp1)
        do i = 1, nmv0
          ii = i + numpts
          do j = 1, 9
            Ptsl(j,ii) = recv(j,i)
          enddo
        enddo
        deallocate(recv)
        deallocate(ms1)
        deallocate(ms2)
        allocate(msk(Nptlocal))
        allocate(ms1(Nptlocal))
        allocate(ms2(Nptlocal))

        iiup = 0
        iidown = 0
        do i = 1, Nptlocal
          msk(i) = .true.
          if(Ptsl(3,i).gt.lcrange(4)) then
            if(myidy.ne.(npy-1)) then
              iiup = iiup + 1
              ms1(iiup) = i
              msk(i) = .false.
            endif
          else if(Ptsl(3,i).le.lcrange(3)) then
            if(myidy.ne.0) then
              iidown = iidown + 1
              ms2(iidown) = i
              msk(i) = .false.
            endif
          else
          endif
        enddo

!        call MPI_BARRIER(comm2d,ierr)
        nout0 = 0
        nmv0 = 0
        allocate(recv(9,nmv0+1))
        do

        if(iiup.gt.nptmv) then
           flag = 1
           iout = nptmv
        else
           flag = 0
           iout = iiup
        endif
        nout = iout
        do i = 1, iout
          k = ms1(i+nout0)
          do j = 1, 9
            sender(j,i) = Ptsl(j,k)
          enddo
        enddo

        allocate(temp1(9,nmv0+1))
        do i = 1, nmv0
          temp1(:,i) = recv(:,i)
        enddo

        do
          jin = 0
          call MPI_IRECV(jin,1,MPI_INTEGER,mydown,0,commcol,req(1),&
                        ierr)
          call MPI_ISEND(iout,1,MPI_INTEGER,myup,0,commcol,req(2),&
                         ierr)
          call MPI_WAITALL(2,req,statarry,ierr) 

          numbuf = jin
          deallocate(recv)
          allocate(recv(9,numbuf+nmv0+1))
          do i = 1, nmv0
            recv(:,i) = temp1(:,i)
          enddo
          deallocate(temp1)

          nst = nmv0 + 1
          !send outgoing particles to left neibhoring processor.
          jin = 9*jin
          iout = 9*iout
          call MPI_IRECV(recv(1,nst),jin,MPI_DOUBLE_PRECISION,mydown,&
                         0,commcol,msid,ierr)
          call MPI_SEND(sender(1,1),iout,MPI_DOUBLE_PRECISION,myup,&
                        0,commcol,ierr)
          call MPI_WAIT(msid,status,ierr) 
          iout = iout/9
          jin = jin/9
          nmv0 = nmv0+jin
        
          allocate(mmsk(numbuf))
          iout = 0
          nmv0 = nmv0 - numbuf
          do i = 1, numbuf
            mmsk(i) = .true.
            ii = i+nmv0
            if(recv(3,ii).gt.lcrange(4)) then
              if(myidy.ne.(npy-1)) then
                iout = iout + 1
                sender(:,iout) = recv(:,ii)
                mmsk(i) = .false.
              endif
            endif
          enddo
          nmv = iout
          call MPI_ALLREDUCE(nmv,totnmv,1,MPI_INTEGER,MPI_SUM, &
                             comm2d,ierr)
          if(totnmv.eq.0) then
            nmv0 = nmv0 + numbuf - nmv
            deallocate(mmsk)
            exit
          endif

          ic = 0
          allocate(temp1(9,nmv0+numbuf-nmv))
          do i = 1, nmv0
            temp1(:,i) = recv(:,i)
          enddo
          do i = 1, numbuf
            ii = i + nmv0
            if(mmsk(i)) then
              ic = ic + 1
              temp1(:,ic+nmv0) = recv(:,ii)
            endif
          enddo
          nmv0 = nmv0 + numbuf - nmv
          deallocate(mmsk)
        enddo

        if(flag.eq.1) then
          nout0 = nout0 + nout
          iiup = iiup - nptmv
        else
          exit
        endif

        enddo

        ! send the particles moving to the down neighoring PE.
        nout0 = 0
        do

        if(iidown.gt.nptmv) then
           flag = 1
           iout = nptmv
        else
           flag = 0
           iout = iidown
        endif
        nout = iout
        do i = 1, iout
          k = ms2(i+nout0)
          do j = 1, 9
            sender(j,i) = Ptsl(j,k)
          enddo
        enddo

        allocate(temp1(9,nmv0+1))
        do i = 1, nmv0
          temp1(:,i) = recv(:,i)
        enddo

        do
          jin = 0
          call MPI_IRECV(jin,1,MPI_INTEGER,myup,0,commcol,req(1),&
                        ierr)
          call MPI_ISEND(iout,1,MPI_INTEGER,mydown,0,commcol,req(2),&
                         ierr)
          call MPI_WAITALL(2,req,statarry,ierr) 

          numbuf = jin
          deallocate(recv)
          allocate(recv(9,numbuf+nmv0+1))
          do i = 1, nmv0
            recv(:,i) = temp1(:,i)
          enddo
          deallocate(temp1)

          nst = nmv0 + 1
          !send outgoing particles to left neibhoring processor.
          jin = 9*jin
          iout = 9*iout
          call MPI_IRECV(recv(1,nst),jin,MPI_DOUBLE_PRECISION,myup,&
                         0,commcol,msid,ierr)
          call MPI_SEND(sender(1,1),iout,MPI_DOUBLE_PRECISION,mydown,&
                        0,commcol,ierr)
          call MPI_WAIT(msid,status,ierr) 
          iout = iout/9
          jin = jin/9
          nmv0 = nmv0+jin
        
          allocate(mmsk(numbuf))
          iout = 0
          nmv0 = nmv0 - numbuf
          do i = 1, numbuf
            mmsk(i) = .true.
            ii = i+nmv0
            if(recv(3,ii).le.lcrange(3)) then
              if(myidy.ne.0) then
                iout = iout + 1
                sender(:,iout) = recv(:,ii)
                mmsk(i) = .false.
              endif
            endif
          enddo
          nmv = iout
          call MPI_ALLREDUCE(nmv,totnmv,1,MPI_INTEGER,MPI_SUM, &
                             comm2d,ierr)
          if(totnmv.eq.0) then
            nmv0 = nmv0 + numbuf - nmv
            deallocate(mmsk)
            exit
          endif

          ic = 0
          allocate(temp1(9,nmv0+numbuf-nmv))
          do i = 1, nmv0
            temp1(:,i) = recv(:,i)
          enddo
          do i = 1, numbuf
            ii = i + nmv0
            if(mmsk(i)) then
              ic = ic + 1
              temp1(:,ic+nmv0) = recv(:,ii)
            endif
          enddo
          nmv0 = nmv0 + numbuf - nmv
          deallocate(mmsk)
        enddo

        if(flag.eq.1) then
          nout0 = nout0 + nout
          iidown = iidown - nptmv
        else
          exit
        endif

        enddo

        numpts = Nptlocal - iiup - iidown
        allocate(temp1(9,numpts))
        ic = 0
        do i = 1, Nptlocal
          if(msk(i)) then
            ic = ic + 1
            do j = 1, 9
              temp1(j,ic) = Ptsl(j,i)
            enddo
          endif
        enddo
        deallocate(msk)
        deallocate(Ptsl)
        Nptlocal = numpts + nmv0
        allocate(Ptsl(9,Nptlocal))
        do i = 1, numpts
          do j = 1, 9
            Ptsl(j,i) = temp1(j,i)
          enddo
        enddo
        deallocate(temp1)
        do i = 1, nmv0
          ii = i + numpts
          do j = 1, 9
            Ptsl(j,ii) = recv(j,i)
          enddo
        enddo
        deallocate(recv)
        deallocate(ms1)
        deallocate(ms2)

        t_ptsmv = t_ptsmv + elapsedtime_Timer(t0)

        end subroutine ptsmv3_ptclmger

        ! move particles from one processor to 4 neighboring processors
        ! in Cylindrical Coordinate.
        subroutine ptsmv4r_Ptclmger(Ptsl,Nptlocal,grid,pdim,npmax,lcrange)
        implicit none
        include 'mpif.h'
        type (Pgrid2d), intent(in) :: grid
        integer, intent(inout) :: Nptlocal
        integer, intent(in) :: pdim,npmax
        double precision, pointer, dimension(:,:) :: Ptsl
        double precision, dimension(:),intent(in) :: lcrange
!        double precision, dimension(pdim,npmax) :: Ptsl
        double precision, allocatable, dimension(:,:) :: left,right,up,down
        double precision, allocatable, dimension(:,:) :: temp1,recv
        integer :: myid,myidx,myidy,totnp,npy,npx, &
                   comm2d,commcol,commrow
        integer :: ileft,iright,iup,idown,iupright,iupleft,&
                   idownleft,idownright
        integer :: jleft,jright,jup,jdown,jupright,jupleft,&
                   jdownleft,jdownright
        integer :: myleft,myright,myup,mydown,myupright,myupleft,&
                   mydwnleft,mydwnright
        integer :: nsmall,i,j,numpts,ic
        integer, dimension(2) :: tmpcoord
        logical, allocatable, dimension(:) :: msk,mmsk
        integer :: numbuf,nmv,nmv0,nout,ii,totnmv,ij
        integer :: msid,ierr
        integer status(MPI_STATUS_SIZE) 
        integer statarry(MPI_STATUS_SIZE,4), req(4)
        double precision :: t0,t1,ri,thi,pi

!        call starttime_Timer(t1)

        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        if(Nptlocal.gt.480) then
          nsmall = Nptlocal/16
        else
          nsmall = Nptlocal
        endif
        allocate(left(9,nsmall))
        allocate(right(9,nsmall))
        allocate(up(9,nsmall))
        allocate(down(9,nsmall))

!        left = 0.0
!        right = 0.0
!        up = 0.0
!        down = 0.0
        ileft = 0
        iright = 0
        iup = 0
        idown = 0

        allocate(msk(Nptlocal))
!        msk = .true.

        pi = 2*asin(1.0)

        do i = 1, Nptlocal
          ri = sqrt(Ptsl(1,i)*Ptsl(1,i)+Ptsl(3,i)*Ptsl(3,i))
          if(Ptsl(1,i).gt.0.0) then
            if(Ptsl(3,i).gt.0.0) then
              thi = asin(Ptsl(3,i)/ri)
            else
              thi = 2*pi+asin(Ptsl(3,i)/ri)
            endif
          else
            thi = pi - asin(Ptsl(3,i)/ri)
          endif
          msk(i) = .true.
          if(Ptsl(5,i).le.lcrange(5)) then
            if(myidx.ne.0) then
              ileft = ileft + 1
              left(:,ileft) = Ptsl(:,i)
              msk(i) = .false.
              if(mod(ileft,nsmall).eq.0) then
                allocate(temp1(9,ileft))
                do j = 1,ileft
                  temp1(:,j) = left(:,j)
                enddo
                deallocate(left)
                allocate(left(9,ileft+nsmall))
                do j = 1,ileft
                  left(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            else
              if(thi.gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = Ptsl(:,i)
                msk(i) = .false.
                if(mod(iup,nsmall).eq.0) then
                  allocate(temp1(9,iup))
                  do j = 1, iup
                    temp1(:,j) = up(:,j)
                  enddo
                  deallocate(up)
                  allocate(up(9,iup+nsmall))
                  do j = 1, iup
                    up(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else if(thi.le.lcrange(3)) then
                if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = Ptsl(:,i)
                msk(i) = .false.
                if(mod(idown,nsmall).eq.0) then
                  allocate(temp1(9,idown))
                  do j = 1, idown
                    temp1(:,j) = down(:,j)
                  enddo
                  deallocate(down)
                  allocate(down(9,idown+nsmall))
                  do j = 1, idown
                    down(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else
              endif
            endif
          else if(Ptsl(5,i).gt.lcrange(6)) then
            if(myidx.ne.(npx-1)) then
              iright = iright + 1
              right(:,iright) = Ptsl(:,i)
              msk(i) = .false.
              if(mod(iright,nsmall).eq.0) then
                allocate(temp1(9,iright))
                do j =1, iright
                  temp1(:,j) = right(:,j)
                enddo
                deallocate(right)
                allocate(right(9,iright+nsmall))
                do j =1, iright
                  right(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            else
              if(thi.gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = Ptsl(:,i)
                msk(i) = .false.
                if(mod(iup,nsmall).eq.0) then
                  allocate(temp1(9,iup))
                  do j = 1, iup
                    temp1(:,j) = up(:,j)
                  enddo
                  deallocate(up)
                  allocate(up(9,iup+nsmall))
                  do j = 1, iup
                    up(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else if(thi.le.lcrange(3)) then
                if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = Ptsl(:,i)
                msk(i) = .false.
                if(mod(idown,nsmall).eq.0) then
                  allocate(temp1(9,idown))
                  do j = 1, idown
                    temp1(:,j) = down(:,j)
                  enddo
                  deallocate(down)
                  allocate(down(9,idown+nsmall))
                  do j = 1, idown
                    down(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else
              endif
            endif
          else if(thi.gt.lcrange(4)) then
            if(myidy.ne.(npy-1)) then
              iup = iup + 1
              up(:,iup) = Ptsl(:,i)
              msk(i) = .false.
              if(mod(iup,nsmall).eq.0) then
                allocate(temp1(9,iup))
                do j = 1, iup
                  temp1(:,j) = up(:,j)
                enddo
                deallocate(up)
                allocate(up(9,iup+nsmall))
                do j = 1, iup
                  up(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            endif
          else if(thi.le.lcrange(3)) then
            if(myidy.ne.0) then
              idown = idown + 1
              down(:,idown) = Ptsl(:,i)
              msk(i) = .false.
              if(mod(idown,nsmall).eq.0) then
                allocate(temp1(9,idown))
                do j = 1, idown
                  temp1(:,j) = down(:,j)
                enddo
                deallocate(down)
                allocate(down(9,idown+nsmall))
                do j = 1, idown
                  down(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            endif
          else
          endif

        enddo

!        t_ptsmv1 = t_ptsmv1 + elapsedtime_Timer(t1)
        nmv0 = 0
        nout = ileft+iright+iup+idown
!        print*,"in comm: ",Nptlocal,nout
        allocate(recv(9,nmv0+1))
        allocate(temp1(9,nmv0+1))
        ij = 0
        call MPI_BARRIER(comm2d,ierr)
        if(myidx.ne.(npx-1)) then
          myright = myidx + 1
        else
          myright = MPI_PROC_NULL
        endif
        if(myidx.ne.0) then
          myleft = myidx - 1
        else
          myleft = MPI_PROC_NULL
        endif 

        if(myidy.ne.npy-1) then
          myup = myidy + 1
        else
          myup = MPI_PROC_NULL
        endif
        if(myidy.ne.0) then
          mydown = myidy -1
        else
          mydown = MPI_PROC_NULL
        endif

        do

        ij = ij + 1
        jleft = 0
        jright = 0

        call MPI_IRECV(jleft,1,MPI_INTEGER,myright,0,commrow,req(1),&
                      ierr)
        call MPI_IRECV(jright,1,MPI_INTEGER,myleft,0,commrow,req(2),&
                        ierr)
        call MPI_ISEND(ileft,1,MPI_INTEGER,myleft,0,commrow,req(3),&
                       ierr)
        call MPI_ISEND(iright,1,MPI_INTEGER,myright,0,commrow,req(4),&
                       ierr)
        call MPI_WAITALL(4,req,statarry,ierr) 

!        call MPI_BARRIER(commrow,ierr)
!        if(myid.eq.0) then
!          print*,"pass 1:"
!        endif

        jup = 0
        jdown = 0

        call MPI_IRECV(jdown,1,MPI_INTEGER,myup,0,commcol,req(1),&
                      ierr)
        call MPI_IRECV(jup,1,MPI_INTEGER,mydown,0,commcol,req(2),&
                      ierr)
        call MPI_ISEND(idown,1,MPI_INTEGER,mydown,0,commcol,req(3),&
                       ierr)
        call MPI_ISEND(iup,1,MPI_INTEGER,myup,0,commcol,req(4),&
                       ierr)
        call MPI_WAITALL(4,req,statarry,ierr) 

        numbuf = jleft+jright+jup+jdown 
        
!        call MPI_BARRIER(commcol,ierr)
!        if(myid.eq.0) then
!          print*,"pass 2:"
!        endif

        deallocate(recv)
        allocate(recv(9,numbuf+nmv0+1))
        do i = 1, nmv0
          recv(:,i) = temp1(:,i)
        enddo
        deallocate(temp1)

        !send outgoing particles to left neibhoring processor.
        !allocate(temp1(9,jleft))
        ! temp1 = 0.0
        jleft = 9*jleft
        ileft = 9*ileft
        !call MPI_IRECV(temp1(1,1),jleft,MPI_DOUBLE_PRECISION,myright,&
        !               0,commrow,msid,ierr)
        call MPI_IRECV(recv(1,nmv0+1),jleft,MPI_DOUBLE_PRECISION,myright,&
                       0,commrow,msid,ierr)
        call MPI_SEND(left(1,1),ileft,MPI_DOUBLE_PRECISION,myleft,&
                      0,commrow,ierr)
        call MPI_WAIT(msid,status,ierr) 
        ileft = ileft/9
        jleft = jleft/9
        !do i = 1, jleft
        !  recv(:,i+nmv0) = temp1(:,i)
        !enddo
        nmv0 = nmv0+jleft
        !deallocate(temp1)
        
        !send outgoing particles to right neibhoring processor.
        !allocate(temp1(9,jright))
        !temp1 = 0.0
        jright = 9*jright
        iright = 9*iright
        !call MPI_IRECV(temp1(1,1),jright,MPI_DOUBLE_PRECISION,myleft,&
        !                0,commrow,msid,ierr)
        call MPI_IRECV(recv(1,nmv0+1),jright,MPI_DOUBLE_PRECISION,myleft,&
                        0,commrow,msid,ierr)
        call MPI_SEND(right(1,1),iright,MPI_DOUBLE_PRECISION,myright,&
                      0,commrow,ierr)
        call MPI_WAIT(msid,status,ierr) 
        iright = iright/9
        jright = jright/9
        !do i = 1, jright
        !  recv(:,i+nmv0) = temp1(:,i)
        !enddo
        nmv0 = nmv0 + jright
        !deallocate(temp1)

!        call MPI_BARRIER(commrow,ierr)
!        if(myid.eq.0) then
!          print*,"pass 3:"
!        endif

        !send outgoing particles to down neibhoring processor.
        !allocate(temp1(9,jdown))
        !temp1 = 0.0
        jdown = 9*jdown
        idown = 9*idown
        !call MPI_IRECV(temp1(1,1),jdown,MPI_DOUBLE_PRECISION,myup,&
        !                0,commcol,msid,ierr)
        call MPI_IRECV(recv(1,nmv0+1),jdown,MPI_DOUBLE_PRECISION,myup,&
                        0,commcol,msid,ierr)
        call MPI_SEND(down(1,1),idown,MPI_DOUBLE_PRECISION,mydown,&
                        0,commcol,ierr)
        call MPI_WAIT(msid,status,ierr)
        idown = idown/9
        jdown = jdown/9
        !do i = 1, jdown
        !  recv(:,i+nmv0) = temp1(:,i)
        !enddo
        nmv0 = nmv0 + jdown
        !deallocate(temp1)

        !send outgoing particles to up neibhoring processor.
        !allocate(temp1(9,jup))
        !temp1 = 0.0
        jup = 9*jup
        iup = 9*iup
        !call MPI_IRECV(temp1(1,1),jup,MPI_DOUBLE_PRECISION,mydown,&
        !              0,commcol,msid,ierr)
        call MPI_IRECV(recv(1,nmv0+1),jup,MPI_DOUBLE_PRECISION,mydown,&
                      0,commcol,msid,ierr)
        call MPI_SEND(up(1,1),iup,MPI_DOUBLE_PRECISION,myup,&
                      0,commcol,ierr)
        call MPI_WAIT(msid,status,ierr)
        iup = iup/9
        jup = jup/9
        !do i = 1, jup
        !  recv(:,i+nmv0) = temp1(:,i)
        !enddo
        nmv0 = nmv0 + jup
        !deallocate(temp1)
 
!        call MPI_BARRIER(commcol,ierr)
!        if(myid.eq.0) then
!          print*,"pass 4:"
!        endif

        deallocate(up)
        deallocate(down)
        deallocate(left)
        deallocate(right)

        if(numbuf.gt.480) then
          nsmall = numbuf/16
        else
          nsmall = numbuf
        endif
        allocate(left(9,nsmall))
        allocate(right(9,nsmall))
        allocate(up(9,nsmall))
        allocate(down(9,nsmall))
        allocate(mmsk(numbuf))
        ileft = 0
        iright = 0
        iup = 0
        idown = 0
        nmv0 = nmv0 - numbuf
        do i = 1, numbuf
          mmsk(i) = .true.
          ii = i+nmv0
          ri = sqrt(recv(1,ii)*recv(1,ii)+recv(3,ii)*recv(3,ii))
          if(recv(1,ii).gt.0.0) then
            if(recv(3,ii).gt.0.0) then
              thi = asin(recv(3,ii)/ri)
            else
              thi = 2*pi+asin(recv(3,ii)/ri)
            endif
          else
            thi = pi - asin(recv(3,ii)/ri)
          endif
          if(recv(5,ii).le.lcrange(5)) then
            if(myidx.ne.0) then
              ileft = ileft + 1
              left(:,ileft) = recv(:,ii)
              mmsk(i) = .false.
              if(mod(ileft,nsmall).eq.0) then
                allocate(temp1(9,ileft))
                do j = 1,ileft
                  temp1(:,j) = left(:,j)
                enddo
                deallocate(left)
                allocate(left(9,ileft+nsmall))
                do j = 1,ileft
                  left(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            else
              if(thi.gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = recv(:,ii)
                mmsk(i) = .false.
                if(mod(iup,nsmall).eq.0) then
                  allocate(temp1(9,iup))
                  do j = 1, iup
                    temp1(:,j) = up(:,j)
                  enddo
                  deallocate(up)
                  allocate(up(9,iup+nsmall))
                  do j = 1, iup
                    up(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else if(thi.le.lcrange(3)) then
                if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = recv(:,ii)
                mmsk(i) = .false.
                if(mod(idown,nsmall).eq.0) then
                  allocate(temp1(9,idown))
                  do j = 1, idown
                    temp1(:,j) = down(:,j)
                  enddo
                  deallocate(down)
                  allocate(down(9,idown+nsmall))
                  do j = 1, idown
                    down(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else
              endif
            endif
          else if(recv(5,ii).gt.lcrange(6)) then
            if(myidx.ne.(npx-1)) then
              iright = iright + 1
              right(:,iright) = recv(:,ii)
              mmsk(i) = .false.
              if(mod(iright,nsmall).eq.0) then
                allocate(temp1(9,iright))
                do j =1, iright
                  temp1(:,j) = right(:,j)
                enddo
                deallocate(right)
                allocate(right(9,iright+nsmall))
                do j =1, iright
                  right(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            else
              if(thi.gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = recv(:,ii)
                mmsk(i) = .false.
                if(mod(iup,nsmall).eq.0) then
                  allocate(temp1(9,iup))
                  do j = 1, iup
                    temp1(:,j) = up(:,j)
                  enddo
                  deallocate(up)
                  allocate(up(9,iup+nsmall))
                  do j = 1, iup
                    up(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else if(thi.le.lcrange(3)) then
                if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = recv(:,ii)
                mmsk(i) = .false.
                if(mod(idown,nsmall).eq.0) then
                  allocate(temp1(9,idown))
                  do j = 1, idown
                    temp1(:,j) = down(:,j)
                  enddo
                  deallocate(down)
                  allocate(down(9,idown+nsmall))
                  do j = 1, idown
                    down(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else
              endif
            endif
          else if(thi.gt.lcrange(4)) then
            if(myidy.ne.(npy-1)) then
              iup = iup + 1
              up(:,iup) = recv(:,ii)
              mmsk(i) = .false.
              if(mod(iup,nsmall).eq.0) then
                allocate(temp1(9,iup))
                do j = 1, iup
                  temp1(:,j) = up(:,j)
                enddo
                deallocate(up)
                allocate(up(9,iup+nsmall))
                do j = 1, iup
                  up(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            endif
          else if(thi.le.lcrange(3)) then
            if(myidy.ne.0) then
              idown = idown + 1
              down(:,idown) = recv(:,ii)
              mmsk(i) = .false.
              if(mod(idown,nsmall).eq.0) then
                allocate(temp1(9,idown))
                do j = 1, idown
                  temp1(:,j) = down(:,j)
                enddo
                deallocate(down)
                allocate(down(9,idown+nsmall))
                do j = 1, idown
                  down(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            endif
          else
          endif
        enddo
        nmv = ileft+iright+idown+iup
        call MPI_ALLREDUCE(nmv,totnmv,1,MPI_INTEGER,MPI_SUM, &
                           comm2d,ierr)
        if(totnmv.eq.0) then
          deallocate(mmsk)
          deallocate(up)
          deallocate(down)
          deallocate(left)
          deallocate(right)
          nmv0 = nmv0 + numbuf - nmv
          exit
        endif

        ic = 0
        allocate(temp1(9,nmv0+numbuf-nmv))
        do i = 1, nmv0
          temp1(:,i) = recv(:,i)
        enddo
        do i = 1, numbuf
          ii = i + nmv0
          if(mmsk(i)) then
            ic = ic + 1
            temp1(:,ic+nmv0) = recv(:,ii)
          endif
        enddo
        deallocate(mmsk)
        nmv0 = nmv0 + numbuf - nmv

!        print*," loop ", ij,numbuf,nmv,nmv0
!        call MPI_BARRIER(comm2d,ierr)

        enddo

        !copy the remaining local particles into a temporary array.
        numpts = Nptlocal-nout
        allocate(temp1(9,numpts))
        ic = 0
        do i = 1, Nptlocal
          if(msk(i)) then
            ic = ic + 1
            do j = 1, 9
              temp1(j,ic) = Ptsl(j,i)
            enddo
          endif
        enddo

        deallocate(msk)

!        call MPI_BARRIER(comm2d,ierr)
        !recopy the remaining local particles back to Ptsl which has 
        !a new size now.
        Nptlocal = numpts+nmv0 
        deallocate(Ptsl)
        allocate(Ptsl(9,Nptlocal))
        do i = 1, numpts
          do j = 1, 9
            Ptsl(j,i) = temp1(j,i)
          enddo
        enddo
        deallocate(temp1)
        do i = 1, nmv0
          ii = i + numpts
          do j = 1, 9
            Ptsl(j,ii) = recv(j,i)
          enddo
        enddo

        deallocate(recv)

        call MPI_BARRIER(comm2d,ierr)
        t_ptsmv = t_ptsmv + elapsedtime_Timer(t0)

        end subroutine ptsmv4r_Ptclmger

        ! move particles from one processor to 4 neighboring processors.
        subroutine ptsmv5r_Ptclmger(Ptsl,Nptlocal,grid,pdim,npmax,lcrange)
        implicit none
        include 'mpif.h'
        type (Pgrid2d), intent(in) :: grid
        integer, intent(inout) :: Nptlocal
        integer, intent(in) :: pdim,npmax
        double precision, pointer, dimension(:,:) :: Ptsl
!        double precision, dimension(pdim,npmax) :: Ptsl
        double precision, dimension(:),intent(in) :: lcrange
        integer, parameter :: nptmv = 100000
        double precision, dimension(9,3*nptmv) :: left,right,up,down
        double precision, allocatable, dimension(:,:) :: temp1,recv
        integer :: myid,myidx,myidy,totnp,npy,npx, &
                   comm2d,commcol,commrow
        integer :: ileft,iright,iup,idown,iupright,iupleft,&
                   idownleft,idownright
        integer :: jleft,jright,jup,jdown,jupright,jupleft,&
                   jdownleft,jdownright
        integer :: myleft,myright,myup,mydown,myupright,myupleft,&
                   mydwnleft,mydwnright
        integer :: nsmall,i,j,numpts,ic
        integer, dimension(2) :: tmpcoord
        logical, dimension(Nptlocal) :: msk
        logical, allocatable, dimension(:) :: mmsk
        integer :: numbuf,nmv,nmv0,nout,ii,totnmv,ij
        integer :: msid,ierr
        integer status(MPI_STATUS_SIZE) 
        integer statarry(MPI_STATUS_SIZE,4), req(4)
        integer :: flag,Nptlocal0,nst,nout0,iileft,iiright,iiup,iidown
        integer :: totflag
        double precision :: t0,ri,thi,pi

        call starttime_Timer(t0)

        pi = 2.0*asin(1.0)

        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)
        if(myidx.ne.(npx-1)) then
          myright = myidx + 1
        else
          myright = MPI_PROC_NULL
        endif
        if(myidx.ne.0) then
          myleft = myidx - 1
        else
          myleft = MPI_PROC_NULL
        endif 

        if(myidy.ne.npy-1) then
          myup = myidy + 1
        else
          myup = MPI_PROC_NULL
        endif
        if(myidy.ne.0) then
          mydown = myidy -1
        else
          mydown = MPI_PROC_NULL
        endif

!        call MPI_BARRIER(comm2d,ierr)

        flag = 0
        Nptlocal0 = Nptlocal
        nout0 = 0

        do 

        ileft = 0
        iright = 0
        iup = 0
        idown = 0
        iileft = 0
        iiright = 0
        iiup = 0
        iidown = 0
        do i = 1, Nptlocal0 - nout0
          ri = sqrt(Ptsl(1,i)*Ptsl(1,i)+Ptsl(3,i)*Ptsl(3,i))
          if(Ptsl(1,i).gt.0.0) then
            if(Ptsl(3,i).gt.0.0) then
              thi = asin(Ptsl(3,i)/ri)
            else
              thi = 2*pi+asin(Ptsl(3,i)/ri)
            endif
          else
            thi = pi - asin(Ptsl(3,i)/ri)
          endif
          msk(i) = .true.
          if(Ptsl(5,i).le.lcrange(5)) then
            if(myidx.ne.0) then
              iileft = iileft + 1
              if(iileft.le.nptmv) then
                ileft = ileft + 1
                left(:,ileft) = Ptsl(:,i)
                msk(i) = .false.
              endif
            else
              if(thi.gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                  iiup = iiup + 1
                  if(iiup.le.nptmv) then
                    iup = iup + 1
                    up(:,iup) = Ptsl(:,i)
                    msk(i) = .false.
                  endif
                endif
              else if(thi.le.lcrange(3)) then
                if(myidy.ne.0) then
                  iidown = iidown + 1
                  if(iidown.le.nptmv) then
                    idown = idown + 1
                    down(:,idown) = Ptsl(:,i)
                    msk(i) = .false.
                  endif
                endif
              else
              endif
            endif
          else if(Ptsl(5,i).gt.lcrange(6)) then
            if(myidx.ne.(npx-1)) then
              iiright = iiright + 1
              if(iiright.le.nptmv) then
                iright = iright + 1
                right(:,iright) = Ptsl(:,i)
                msk(i) = .false.
              endif
            else
              if(thi.gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                  iiup = iiup + 1
                  if(iiup.le.nptmv) then
                    iup = iup + 1
                    up(:,iup) = Ptsl(:,i)
                    msk(i) = .false.
                  endif
                endif
              else if(thi.le.lcrange(3)) then
                if(myidy.ne.0) then
                  iidown = iidown + 1
                  if(iidown.le.nptmv) then
                    idown = idown + 1
                    down(:,idown) = Ptsl(:,i)
                    msk(i) = .false.
                  endif
                endif
              else
              endif
            endif
          else if(thi.gt.lcrange(4)) then
            if(myidy.ne.(npy-1)) then
              iiup = iiup + 1
              if(iiup.le.nptmv) then
                iup = iup + 1
                up(:,iup) = Ptsl(:,i)
                msk(i) = .false.
              endif
            endif
          else if(thi.le.lcrange(3)) then
            if(myidy.ne.0) then
              iidown = iidown + 1
              if(iidown.le.nptmv) then
                idown = idown + 1
                down(:,idown) = Ptsl(:,i)
                msk(i) = .false.
              endif
            endif
          else
          endif

        enddo

        if((iileft.gt.nptmv).or.(iiright.gt.nptmv).or.(iiup.gt.nptmv) &
            .or.(iidown.gt.nptmv)) then
           flag = 1
        else
           flag = 0
        endif

        nmv0 = 0
        nout = ileft+iright+iup+idown
        allocate(recv(9,nmv0+1))
        allocate(temp1(9,nmv0+1))
        ij = 0
        call MPI_BARRIER(comm2d,ierr)

        do

        ij = ij + 1
        jleft = 0
        jright = 0

        call MPI_IRECV(jleft,1,MPI_INTEGER,myright,0,commrow,req(1),&
                      ierr)
        call MPI_IRECV(jright,1,MPI_INTEGER,myleft,0,commrow,req(2),&
                        ierr)
        call MPI_ISEND(ileft,1,MPI_INTEGER,myleft,0,commrow,req(3),&
                       ierr)
        call MPI_ISEND(iright,1,MPI_INTEGER,myright,0,commrow,req(4),&
                       ierr)
        call MPI_WAITALL(4,req,statarry,ierr) 

!        call MPI_BARRIER(commrow,ierr)
!        if(myid.eq.0) then
!          print*,"pass 1:"
!        endif

        jup = 0
        jdown = 0

        call MPI_IRECV(jdown,1,MPI_INTEGER,myup,0,commcol,req(1),&
                      ierr)
        call MPI_IRECV(jup,1,MPI_INTEGER,mydown,0,commcol,req(2),&
                      ierr)
        call MPI_ISEND(idown,1,MPI_INTEGER,mydown,0,commcol,req(3),&
                       ierr)
        call MPI_ISEND(iup,1,MPI_INTEGER,myup,0,commcol,req(4),&
                       ierr)
        call MPI_WAITALL(4,req,statarry,ierr) 

        numbuf = jleft+jright+jup+jdown 
        
!        call MPI_BARRIER(commcol,ierr)
!        if(myid.eq.0) then
!          print*,"pass 2:"
!        endif

        deallocate(recv)
        allocate(recv(9,numbuf+nmv0+1))
        do i = 1, nmv0
          recv(:,i) = temp1(:,i)
        enddo
        deallocate(temp1)

        nst = nmv0 + 1
        !send outgoing particles to left neibhoring processor.
        jleft = 9*jleft
        ileft = 9*ileft
        call MPI_IRECV(recv(1,nst),jleft,MPI_DOUBLE_PRECISION,myright,&
                       0,commrow,msid,ierr)
        call MPI_SEND(left(1,1),ileft,MPI_DOUBLE_PRECISION,myleft,&
                      0,commrow,ierr)
        call MPI_WAIT(msid,status,ierr) 
        ileft = ileft/9
        jleft = jleft/9
        nmv0 = nmv0+jleft
        
        nst = nmv0 + 1
        !send outgoing particles to right neibhoring processor.
        jright = 9*jright
        iright = 9*iright
        call MPI_IRECV(recv(1,nst),jright,MPI_DOUBLE_PRECISION,myleft,&
                        0,commrow,msid,ierr)
        call MPI_SEND(right(1,1),iright,MPI_DOUBLE_PRECISION,myright,&
                      0,commrow,ierr)
        call MPI_WAIT(msid,status,ierr) 
        iright = iright/9
        jright = jright/9
        nmv0 = nmv0 + jright

!        call MPI_BARRIER(commrow,ierr)
!        if(myid.eq.0) then
!          print*,"pass 3:"
!        endif

        nst = nmv0 + 1
        !send outgoing particles to down neibhoring processor.
        jdown = 9*jdown
        idown = 9*idown
        call MPI_IRECV(recv(1,nst),jdown,MPI_DOUBLE_PRECISION,myup,&
                        0,commcol,msid,ierr)
        call MPI_SEND(down(1,1),idown,MPI_DOUBLE_PRECISION,mydown,&
                        0,commcol,ierr)
        call MPI_WAIT(msid,status,ierr)
        idown = idown/9
        jdown = jdown/9
        nmv0 = nmv0 + jdown

        nst = nmv0 + 1
        !send outgoing particles to up neibhoring processor.
        jup = 9*jup
        iup = 9*iup
        call MPI_IRECV(recv(1,nst),jup,MPI_DOUBLE_PRECISION,mydown,&
                      0,commcol,msid,ierr)
        call MPI_SEND(up(1,1),iup,MPI_DOUBLE_PRECISION,myup,&
                      0,commcol,ierr)
        call MPI_WAIT(msid,status,ierr)
        iup = iup/9
        jup = jup/9
        nmv0 = nmv0 + jup
 
!        call MPI_BARRIER(commcol,ierr)
!        if(myid.eq.0) then
!          print*,"pass 4:"
!        endif

        allocate(mmsk(numbuf))
        ileft = 0
        iright = 0
        iup = 0
        idown = 0
        nmv0 = nmv0 - numbuf
        do i = 1, numbuf
          ii = i+nmv0
          ri = sqrt(recv(1,ii)*recv(1,ii)+recv(3,ii)*recv(3,ii))
          if(recv(1,ii).gt.0.0) then
            if(recv(3,ii).gt.0.0) then
              thi = asin(recv(3,ii)/ri)
            else
              thi = 2*pi+asin(recv(3,ii)/ri)
            endif
          else
            thi = pi - asin(recv(3,ii)/ri)
          endif
          mmsk(i) = .true.
          if(recv(5,ii).le.lcrange(5)) then
            if(myidx.ne.0) then
              ileft = ileft + 1
              left(:,ileft) = recv(:,ii)
              mmsk(i) = .false.
            else
              if(thi.gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = recv(:,ii)
                mmsk(i) = .false.
                endif
              else if(thi.le.lcrange(3)) then
                if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = recv(:,ii)
                mmsk(i) = .false.
                endif
              else
              endif
            endif
          else if(recv(5,ii).gt.lcrange(6)) then
            if(myidx.ne.(npx-1)) then
              iright = iright + 1
              right(:,iright) = recv(:,ii)
              mmsk(i) = .false.
            else
              if(thi.gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = recv(:,ii)
                mmsk(i) = .false.
                endif
              else if(thi.le.lcrange(3)) then
                if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = recv(:,ii)
                mmsk(i) = .false.
                endif
              else
              endif
            endif
          else if(thi.gt.lcrange(4)) then
            if(myidy.ne.(npy-1)) then
              iup = iup + 1
              up(:,iup) = recv(:,ii)
              mmsk(i) = .false.
            endif
          else if(thi.le.lcrange(3)) then
            if(myidy.ne.0) then
              idown = idown + 1
              down(:,idown) = recv(:,ii)
              mmsk(i) = .false.
            endif
          else
          endif
        enddo
        nmv = ileft+iright+idown+iup
        call MPI_ALLREDUCE(nmv,totnmv,1,MPI_INTEGER,MPI_SUM, &
                           comm2d,ierr)
        if(totnmv.eq.0) then
          nmv0 = nmv0 + numbuf - nmv
          deallocate(mmsk)
          exit
        endif

        ic = 0
        allocate(temp1(9,nmv0+numbuf-nmv))
        do i = 1, nmv0
          temp1(:,i) = recv(:,i)
        enddo
        do i = 1, numbuf
          ii = i + nmv0
          if(mmsk(i)) then
            ic = ic + 1
            temp1(:,ic+nmv0) = recv(:,ii)
          endif
        enddo
        nmv0 = nmv0 + numbuf - nmv
        deallocate(mmsk)

!        call MPI_BARRIER(comm2d,ierr)

        enddo

        !copy the remaining local particles into a temporary array.
        numpts = Nptlocal-nout
        allocate(temp1(9,numpts))
        ic = 0
        do i = 1, Nptlocal0-nout0
          if(msk(i)) then
            ic = ic + 1
            do j = 1, 9
              temp1(j,ic) = Ptsl(j,i)
            enddo
          endif
        enddo
        do i = Nptlocal0-nout0+1, Nptlocal
          ii = i-nout
          do j = 1, 9
            temp1(j,ii) = Ptsl(j,i)
          enddo
        enddo
 
!        call MPI_BARRIER(comm2d,ierr)
        !recopy the remaining local particles back to Ptsl which has 
        !a new size now.
        Nptlocal = numpts+nmv0 
        deallocate(Ptsl)
        allocate(Ptsl(9,Nptlocal))
        do i = 1, numpts
          do j = 1, 9
            Ptsl(j,i) = temp1(j,i)
          enddo
        enddo
        deallocate(temp1)
        do i = 1, nmv0
          ii = i + numpts
          do j = 1, 9
            Ptsl(j,ii) = recv(j,i)
          enddo
        enddo

        deallocate(recv)

        nout0 = nout0 + nout

          call MPI_ALLREDUCE(flag,totflag,1,MPI_INTEGER,MPI_SUM, &
                           comm2d,ierr)
          if(totflag.eq.0) exit

        enddo

        t_ptsmv = t_ptsmv + elapsedtime_Timer(t0)

        end subroutine ptsmv5r_Ptclmger


      end module Ptclmgerclass
