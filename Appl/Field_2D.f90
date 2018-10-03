!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! FieldQuantclass: 3D field quantity class in Field module of APPLICATION
!                 layer.
! Version: beta
! Author: Ji Qiang 
! Description: This class defines a 3-D field quantity in the accelerator.
!              The field quantity can be updated at each step. 
! Comments:
!----------------------------------------------------------------
      module FieldQuantclass
        use Timerclass
        use CompDomclass
        use Pgrid2dclass
        use FFTclass
        use Besselclass
        use Transposeclass
        use PhysConstclass
        use Dataclass
        type FieldQuant
!          private
          !# of mesh points in x and y directions.
          integer :: Nx,Ny,Nz,Nxlocal,Nylocal,Nzlocal
          !# field quantity array.
          double precision, pointer, dimension(:,:,:) :: FieldQ
        end type FieldQuant
      contains
        !Initialize field class.
        subroutine construct_FieldQuant(this,innx,inny,innz,geom,grid) 
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(out) :: this
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer, intent(in) :: innx, inny, innz 
        integer :: myid, myidx, myidy, nptot,nproccol,nprocrow
        integer, allocatable, dimension(:,:,:) :: LocalTable

        call getsize_Pgrid2d(grid,nptot,nproccol,nprocrow) 
        call getpost_Pgrid2d(grid,myid,myidy,myidx)

        allocate(LocalTable(2,0:nprocrow-1,0:nproccol-1))
        call getlctabnm_CompDom(geom,LocalTable)

        this%Nx = innx
        this%Ny = inny
        this%Nz = innz
        this%Nxlocal = innx
        if(nproccol.gt.1) then
          this%Nylocal = LocalTable(2,myidx,myidy)+2
        else
          this%Nylocal = LocalTable(2,myidx,myidy)
        endif
        if(nprocrow.gt.1) then
          this%Nzlocal = LocalTable(1,myidx,myidy)+2
        else
          this%Nzlocal = LocalTable(1,myidx,myidy)
        endif
  
        allocate(this%FieldQ(this%Nxlocal,this%Nylocal,this%Nzlocal))
        this%FieldQ = 0.0

        deallocate(LocalTable)

        end subroutine construct_FieldQuant
   
        ! set field quantity.
        subroutine set_FieldQuant(this,innx,inny,innz,geom,grid,&
                                  nprx,npry) 
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(inout) :: this
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer, intent(in) :: innx, inny, innz, nprx, npry 
        integer :: myid, myidx, myidy
        integer, dimension(2,0:nprx-1,0:npry-1)::LocalTable

        call getpost_Pgrid2d(grid,myid,myidy,myidx)

        call getlctabnm_CompDom(geom,LocalTable)

        this%Nx = innx
        this%Ny = inny
        this%Nz = innz
        this%Nxlocal = innx
        if(npry.gt.1) then
          this%Nylocal = LocalTable(2,myidx,myidy)+2
        else
          this%Nylocal = LocalTable(2,myidx,myidy)
        endif
        if(nprx.gt.1) then
          this%Nzlocal = LocalTable(1,myidx,myidy)+2
        else
          this%Nzlocal = LocalTable(1,myidx,myidy)
        endif
        this%Nxlocal = innx
  
        deallocate(this%FieldQ) 
        allocate(this%FieldQ(this%Nxlocal,this%Nylocal,this%Nzlocal)) 
        this%FieldQ = 0.0

        end subroutine set_FieldQuant

!----------------------------------------------------------------------
! update potential (solving Possion's equation) with 3D isolated 
! boundary conditions.
        subroutine update3O_FieldQuant(this,source,fldgeom,grid,nxlc,&
          nylc,nzlc,nprocrow,nproccol,nylcr,nzlcr)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(in) :: fldgeom
        double precision, dimension(nxlc,nylc,nzlc), intent(in) :: source
        type (Pgrid2d), intent(in) :: grid
        integer, intent(in) :: nxlc,nylc,nzlc,&
                               nprocrow,nproccol,nylcr,nzlcr
        type (FieldQuant), intent(inout) :: this
        double precision, dimension(3) :: msize
        double precision :: hx, hy, hz, temp
        integer, dimension(0:nprocrow-1) :: ypzstable,pztable
        integer, dimension(0:nproccol-1) :: xpystable,pytable
        double precision , dimension(nxlc,nylcr,nzlcr) :: rho
        integer :: myid,myidz,myidy,&
                   comm2d,commcol,commrow
        integer :: nxpylc2,nypzlc2
        integer :: nsxy1,nsxy2,nsyz1,nsyz2
        integer :: i,j,k,inxglb,inyglb,inzglb,innx,inny,innz
        integer, dimension(2,0:nprocrow-1,0:nproccol-1)::LocalTable
        integer, dimension(3) :: glmshnm
        integer :: jadd,kadd,ierr
        double precision :: t0

        call getpost_Pgrid2d(grid,myid,myidy,myidz)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        call getlctabnm_CompDom(fldgeom,LocalTable)
        
        do i = 0, nprocrow-1
          pztable(i) = LocalTable(1,i,0)
        enddo
        do i = 0, nproccol-1
          pytable(i) = LocalTable(2,0,i)
        enddo

        call getmnum_CompDom(fldgeom,glmshnm)
        inxglb = glmshnm(1) 
        inyglb = glmshnm(2)
        inzglb = glmshnm(3)

        innz = nzlcr
        inny = nylcr
        innx = nxlc

        if(nprocrow.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif
        if(nproccol.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif
        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              rho(i,j,k) = source(i,j+jadd,k+kadd)
            enddo
          enddo
        enddo
        
        ! +1 is from the real to complex fft.
        nsxy1 = (inxglb+1)/nproccol
        nsxy2 = (inxglb+1) - nproccol*nsxy1
        do i = 0, nproccol-1
          if(i.le.(nsxy2-1)) then
            xpystable(i) = nsxy1+1
          else
            xpystable(i) = nsxy1
          endif
        enddo

        nsyz1 = 2*inyglb/nprocrow
        nsyz2 = 2*inyglb - nprocrow*nsyz1
        do i = 0, nprocrow-1
          if(i.le.(nsyz2-1)) then
            ypzstable(i) = nsyz1+1
          else
            ypzstable(i) = nsyz1
          endif
        enddo

        nxpylc2 = xpystable(myidy)
        nypzlc2 = ypzstable(myidz)

        call getmsize_CompDom(fldgeom,msize)
        hx = msize(1)
        hy = msize(2)
        hz = msize(3)

        ! Open boundary conditions!
        call openBC3D(innx,inny,innz,rho,hx,hy,hz, &
        nxpylc2,nypzlc2,myidz,myidy,nprocrow,nproccol,commrow,commcol,&
        comm2d,pztable,pytable,ypzstable,xpystable,&
        inxglb,inyglb,inzglb)

        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              this%FieldQ(i,j+jadd,k+kadd) = rho(i,j,k)
            enddo
          enddo
        enddo

        call MPI_BARRIER(comm2d,ierr)
        t_field = t_field + elapsedtime_Timer(t0)

        end subroutine update3O_FieldQuant

        ! Solving Poisson's equation with open BCs.
        subroutine openBC3D(innx,inny,innz,rho,hx,hy,hz,&
           nxpylc2,nypzlc2,myidz,myidy,npz,npy,commrow,commcol,comm2d, &
           pztable,pytable,ypzstable,xpystable,inxglb,&
           inyglb,inzglb)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx,inny,innz,inxglb,inyglb,inzglb
        integer, intent(in) :: nxpylc2,nypzlc2
        integer, intent(in) :: myidz,myidy,npz,npy,commrow,commcol,&
                               comm2d
        integer, dimension(0:npz-1), intent(in) :: pztable,ypzstable
        integer, dimension(0:npy-1), intent(in) :: pytable,xpystable
        double precision, dimension(innx,inny,innz), intent(inout) :: rho
        double precision, intent(in) :: hx, hy, hz
        double precision :: scalex,scaley,scalez
        integer :: i,j,k,n1,n2,n3,nylc22,nzlc22
        double complex, allocatable, dimension(:,:,:) :: rho2out
        double complex, allocatable, dimension(:,:,:) :: grn
        integer :: ginny,ginnz,ierr

        n1 = 2*inxglb
        n2 = 2*inyglb
        n3 = 2*inzglb

        nylc22 = nxpylc2
        nzlc22 = nypzlc2

        scalex = 1.0
        scaley = 1.0
        scalez = 1.0

        allocate(rho2out(n3,nylc22,nzlc22))

        call fft3d1_FFT(n1,n2,n3,innz,inny,nylc22,nzlc22,&
                1,scalex,scaley,scalez,rho,ypzstable,pztable,&
        xpystable,pytable,npz,commrow,npy,commcol,comm2d,myidz,myidy, &
        rho2out)

        !c compute FFT of the Green function on the grid:
        ! here the +1 is from the unsymmetry of green function
        ! on double-sized grid.
        if(myidz.eq.(npz-1)) then
           ginnz = innz + 1
        else
           ginnz = innz
        endif
        if(myidy.eq.(npy-1)) then
           ginny = inny + 1
        else
           ginny = inny
        endif
        allocate(grn(n3,nylc22,nzlc22))
        !call greenf1(inxglb,inyglb,inzglb,ginnz,ginny,nylc22,nzlc22, &
        !       hx,hy,hz,myidz,npz,commrow,myidy,npy,commcol,comm2d,&
        !          ypzstable,pztable,xpystable,pytable,grn)
        call greenf1tIntnew(inxglb,inyglb,inzglb,ginnz,ginny,nylc22,nzlc22, &
               hx,hy,hz,myidz,npz,commrow,myidy,npy,commcol,comm2d,&
                  ypzstable,pztable,xpystable,pytable,grn)

        ! multiply transformed charge density and transformed Green 
        ! function:
        do k = 1, nzlc22
          do j = 1, nylc22
            do i = 1, n3
              rho2out(i,j,k) = rho2out(i,j,k)*grn(i,j,k)
            enddo
          enddo
        enddo

        deallocate(grn)

        ! inverse FFT:
        scalex = 1.0/float(n1)
        scaley = 1.0/float(n2)
        scalez = 1.0/float(n3)
        call invfft3d1_FFT(n3,n2,n1,nylc22,nzlc22,inny,innz,&
               -1,scalex,scaley,scalez,rho2out,pztable,ypzstable,&
        pytable,xpystable,npz,commrow,npy,commcol,comm2d,myidz,myidy,&
        rho)

        deallocate(rho2out)

        end subroutine openBC3D

        ! green function for extended array.
        subroutine greenf1(nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz,&
                  hx,hy,hz,myidx,npx,commrow,myidy,npy,commcol,comm2d,&
                   xstable,xrtable,ystable,yrtable,grnout)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz
        integer, intent(in) :: myidx,myidy,npx,npy,commrow,commcol,&
                               comm2d
        integer, dimension(0:npx-1),intent(in) :: xstable,xrtable
        integer, dimension(0:npy-1),intent(in) :: ystable,yrtable
        double precision, intent(in) :: hx, hy, hz
        double complex, intent (out), &
                       dimension (2*nz,nsizexy,nsizeyz) :: grnout
        integer, dimension(0:npx-1) :: gxrtable
        integer, dimension(0:npy-1) :: gyrtable
        integer :: i,j,k,kk,ii,jj,iii,jjj,kkk,n1,n2,n3,ns1,ns2
        integer :: nblockx,nblocky,ksign,nxx,ierr
        double precision, dimension (nx+1,nsizey,nsizez) :: grn
        double precision :: scalex,scaley,scalez
        double precision :: t0
        double precision, dimension(2*nx,nsizey) :: tmp1
        double complex, dimension(nx+1,nsizey) :: tmp10
        double complex, dimension(2*ny,nsizexy) :: tmp2
        double complex, dimension(2*nz,nsizexy) :: tmp3
        double complex, allocatable, dimension(:,:,:) :: x1
        double complex, allocatable, dimension(:,:,:) :: x0

        call starttime_Timer(t0)

        gxrtable = xrtable
        gxrtable(npx-1) = xrtable(npx-1) + 1
        nblockx = 0
        do i = 0, myidx-1
          nblockx = nblockx + gxrtable(i)
        enddo

        gyrtable = yrtable
        gyrtable(npy-1) = yrtable(npy-1) + 1
        nblocky = 0
        do i = 0, myidy-1
          nblocky = nblocky + gyrtable(i)
        enddo

        do k = 1, nsizez
          do j = 1, nsizey
            do i = 1, nx+1
              jj = j + nblocky
              kk = k + nblockx
                kkk = kk - 1
                jjj = jj - 1
                iii = i - 1
              if((iii*iii+jjj*jjj+kkk*kkk) .ne. 0) then
               grn(i,j,k)=1.0/sqrt((hx*iii)**2+(hy*jjj)**2+(hz*kkk)**2)
              endif
            enddo
          enddo
        enddo
        if((myidx.eq.0).and.(myidy.eq.0)) then
          if(nsizez.gt.1) then
            grn(1,1,1) = grn(1,1,2)
          else
            grn(1,1,1) = 1.0
          endif
        endif

        scalex = 1.0
        scaley = 1.0
        scalez = 1.0
        n1 = 2*nx
        n2 = 2*ny
        n3 = 2*nz
        ksign = 1

        nxx = n1/2 + 1
        !FFTs along y and z dimensions.
        allocate(x0(n1/2+1,nsizey,nsizez))
        do k = 1, nsizez
          do j = 1, nsizey 
            do i = 1, n1/2+1
              tmp1(i,j) = grn(i,j,k)
            enddo
            do i = n1/2+2, n1
              tmp1(i,j) = grn(n1-i+2,j,k)
            enddo
          enddo

          ! FFTs along z dimensions:
          call fftrclocal_FFT(ksign,scalex,tmp1,n1,nsizey,tmp10)

          do j = 1, nsizey 
            do i = 1, n1/2+1
              x0(i,j,k) = tmp10(i,j)
            enddo
          enddo
        enddo

        allocate(x1(n2/2+1,nsizexy,nsizez))

        ! FFTs along y dimensions:
!        call MPI_BARRIER(commcol,ierr)
        call trans3d_TRANSP(nxx,n2/2+1,nsizexy,nsizey,x0,x1,npy,&
                     ystable,gyrtable,commcol,nsizez)
        deallocate(x0)
        allocate(x0(n2,nsizexy,nsizez))

        do k = 1, nsizez
          do j = 1, nsizexy 
            do i = 1, n2/2+1
              tmp2(i,j) = x1(i,j,k)
            enddo
            do i = n2/2+2,n2
              tmp2(i,j) = x1(n2-i+2,j,k)
            enddo
          enddo

          call fftlocal_FFT(ksign,scaley,tmp2,n2,nsizexy)

          do j = 1, nsizexy 
            do i = 1, n2
              x0(i,j,k) = tmp2(i,j) 
            enddo
          enddo
        enddo

        deallocate(x1)
        allocate(x1(n3/2+1,nsizexy,nsizeyz))
!        call MPI_BARRIER(commcol,ierr)
        call trans3d3_TRANSP(n2,nsizexy,nsizez,nsizeyz,x0,x1,npx,&
                      xstable,gxrtable,commrow,myidx,n3/2+1)
        deallocate(x0)

        do k = 1, nsizeyz
          do j = 1, nsizexy
            do i = 1, n3/2+1
              tmp3(i,j) = x1(i,j,k) 
            enddo
            do i = n3/2+2,n3
              tmp3(i,j) = x1(n3-i+2,j,k)
            enddo
          enddo

          call fftlocal_FFT(ksign,scalez,tmp3,n3,nsizexy)

          do j = 1, nsizexy
            do i = 1, n3
              grnout(i,j,k) = tmp3(i,j)
            enddo
          enddo
        enddo

        deallocate(x1)

        t_greenf = t_greenf + elapsedtime_Timer(t0)

        end subroutine greenf1

        ! green function for extended array.
        subroutine greenf1tIntnew(nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz,&
                  hx,hy,hz,myidx,npx,commrow,myidy,npy,commcol,comm2d,&
                   xstable,xrtable,ystable,yrtable,grnout)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz
        integer, intent(in) :: myidx,myidy,npx,npy,commrow,commcol,&
                               comm2d
        integer, dimension(0:npx-1),intent(in) :: xstable,xrtable
        integer, dimension(0:npy-1),intent(in) :: ystable,yrtable
        double precision, intent(in) :: hx, hy, hz
        double complex, intent (out), &
                       dimension (2*nz,nsizexy,nsizeyz) :: grnout
        integer, dimension(0:npx-1) :: gxrtable
        integer, dimension(0:npy-1) :: gyrtable
        integer :: i,j,k,kk,ii,jj,iii,jjj,kkk,n1,n2,n3,ns1,ns2
        integer :: nblockx,nblocky,ksign,nxx,ierr
        double precision, dimension (nx+1,nsizey,nsizez) :: grn
        double precision, dimension (nx+2,nsizey+1,nsizez+1) :: grntmp
        double precision :: scalex,scaley,scalez
        double precision :: t0
        double precision, dimension(2*nx,nsizey) :: tmp1
        double complex, dimension(nx+1,nsizey) :: tmp10
        double complex, dimension(2*ny,nsizexy) :: tmp2
        double complex, dimension(2*nz,nsizexy) :: tmp3
        double complex, allocatable, dimension(:,:,:) :: x1
        double complex, allocatable, dimension(:,:,:) :: x0
        double precision :: rr,aa,bb,cc,dd,ee,ff,ss
        double complex :: gg,gg2
        double complex :: ggrr
        double precision, dimension(2) :: xx,yy,zz
        double precision, dimension(3) :: vv
        integer :: n,i0,j0,k0
        double precision :: recfourpi

        recfourpi = 1.0/(8.0*asin(1.0))

        call starttime_Timer(t0)

!        if(myidx.eq.0 .and. myidy.eq.0) then
!          print*,"into integrated Green function....."
!        endif
        gxrtable = xrtable
        gxrtable(npx-1) = xrtable(npx-1) + 1
        nblockx = 0
        do i = 0, myidx-1
          nblockx = nblockx + gxrtable(i)
        enddo

        gyrtable = yrtable
        gyrtable(npy-1) = yrtable(npy-1) + 1
        nblocky = 0
        do i = 0, myidy-1
          nblocky = nblocky + gyrtable(i)
        enddo

        do k0 = 1, nsizez+1
          do j0 = 1, nsizey+1
            do i0 = 1, nx+2
              jj = j0 + nblocky
              kk = k0 + nblockx
              kkk = kk - 1
              jjj = jj - 1
              iii = i0 - 1

              vv(1) = iii*hx-hx/2
              vv(2) = jjj*hy-hy/2
              vv(3) = kkk*hz-hz/2
            
                rr = sqrt(vv(1)**2+vv(2)**2+vv(3)**2)
                aa = vv(1)**2*vv(3)+(vv(2)**2+vv(3)**2)*vv(3) + &
                            vv(1)*vv(3)*rr
                bb = vv(1)**2*vv(2) + vv(1)*vv(2)*rr
                cc = vv(1)*(vv(1)**2+vv(2)**2)+vv(3)*vv(1)*(vv(3)+rr)
                dd = vv(3)*vv(2)*(vv(3) + rr)
                ee = vv(1)**2*vv(2)+vv(2)*(vv(2)**2+vv(3)*(vv(3)+rr))
                ff = vv(1)*vv(3)*(vv(3)+rr)
                ss = 4*vv(2)*vv(3)*log(vv(1)+rr) + &
                        4*vv(1)*vv(3)*log(vv(2)+rr) + &
                        4*vv(1)*vv(2)*log(vv(3)+rr)

                gg2 = cmplx(0.0,vv(3)**2)*log(cmplx(aa**2-bb**2,2*aa*bb)/ &
                        (aa**2+bb**2) ) + cmplx(0.0,vv(1)**2)*&
                        log(cmplx(cc**2-dd**2,2*cc*dd )/(cc**2+dd**2))+&
                        cmplx(0.0,vv(2)**2)*log(cmplx(ee**2-ff**2,2*ee*ff)/ &
                        (ee**2+ff**2) ) + ss
               
              !grntmp(i0,j0,k0) = recfourpi*real(gg2)/(4*hx*hy*hz) !wrong in Z code
              grntmp(i0,j0,k0) = real(gg2)/(4*hx*hy*hz)

            enddo
          enddo
        enddo

        do k0 = 1, nsizez
          do j0 = 1, nsizey
            do i0 = 1, nx+1
              grn(i0,j0,k0) = grntmp(i0+1,j0+1,k0+1)-grntmp(i0,j0+1,k0+1)-grntmp(i0+1,j0,k0+1)+&
                              grntmp(i0,j0,k0+1)-grntmp(i0+1,j0+1,k0)+grntmp(i0,j0+1,k0)+&
                              grntmp(i0+1,j0,k0)-grntmp(i0,j0,k0)
            enddo
          enddo
        enddo

!              xx(1) = iii*hx-hx/2
!              xx(2) = iii*hx+hx/2
!              yy(1) = jjj*hy-hy/2
!              yy(2) = jjj*hy+hy/2
!              zz(1) = kkk*hz-hz/2
!              zz(2) = kkk*hz+hz/2
!       
!              !find the integrated Green function.
!              n = 0
!              do k = 1, 2
!                do j = 1, 2
!                  do i = 1, 2
!                    n = n+1
!                    rr(n) = sqrt(xx(i)**2+yy(j)**2+zz(k)**2)
!                    aa(n) = xx(i)**2*zz(k)+(yy(j)**2+zz(k)**2)*zz(k) + &
!                            xx(i)*zz(k)*rr(n)
!                    bb(n) = xx(i)**2*yy(j) + xx(i)*yy(j)*rr(n)
!                    cc(n) = xx(i)*(xx(i)**2+yy(j)**2)+zz(k)*xx(i)*(zz(k)+rr(n))
!                    dd(n) = zz(k)*yy(j)*(zz(k) + rr(n))
!                    ee(n) = xx(i)**2*yy(j)+yy(j)*(yy(j)**2+zz(k)*(zz(k)+rr(n)))
!                    ff(n) = xx(i)*zz(k)*(zz(k)+rr(n))
!                    ss(n) = 4*yy(j)*zz(k)*log(xx(i)+rr(n)) + &
!                            4*xx(i)*zz(k)*log(yy(j)+rr(n)) + &
!                            4*xx(i)*yy(j)*log(zz(k)+rr(n))
!                    gg(n) = cmplx(0.0,zz(k)**2)*log(cmplx(aa(n)**2-bb(n)**2,2*aa(n)*bb(n))/ &
!                            (aa(n)**2+bb(n)**2) ) + cmplx(0.0,xx(i)**2)*&
!                            log(cmplx(cc(n)**2-dd(n)**2,2*cc(n)*dd(n) )/(cc(n)**2+dd(n)**2))+&
!                            cmplx(0.0,yy(j)**2)*log(cmplx(ee(n)**2-ff(n)**2,2*ee(n)*ff(n))/ &
!                            (ee(n)**2+ff(n)**2) )
!                    gg2(n) = ss(n) +  gg(n)
!                  enddo
!                enddo
!              enddo
!              ggrr = (-gg2(1)+gg2(2)+gg2(3)-gg2(4)+gg2(5)-gg2(6)-gg2(7)+gg2(8))/4
!              grn(i0,j0,k0) = real(ggrr)/(hx*hy*hz)
!            enddo
!          enddo
!        enddo


        if((myidx.eq.0).and.(myidy.eq.0)) then
          if(nsizez.gt.1) then
            grn(1,1,1) = grn(1,1,2)
          else
            grn(1,1,1) = 1.0
          endif
        endif

        scalex = 1.0
        scaley = 1.0
        scalez = 1.0
        n1 = 2*nx
        n2 = 2*ny
        n3 = 2*nz
        ksign = 1

        nxx = n1/2 + 1
        !FFTs along y and z dimensions.
        allocate(x0(n1/2+1,nsizey,nsizez))
        do k = 1, nsizez
          do j = 1, nsizey 
            do i = 1, n1/2+1
              tmp1(i,j) = grn(i,j,k)
            enddo
            do i = n1/2+2, n1
              tmp1(i,j) = grn(n1-i+2,j,k)
            enddo
          enddo

          ! FFTs along z dimensions:
          call fftrclocal_FFT(ksign,scalex,tmp1,n1,nsizey,tmp10)

          do j = 1, nsizey 
            do i = 1, n1/2+1
              x0(i,j,k) = tmp10(i,j)
            enddo
          enddo
        enddo

        allocate(x1(n2/2+1,nsizexy,nsizez))

        ! FFTs along y dimensions:
!        call MPI_BARRIER(commcol,ierr)
        call trans3d_TRANSP(nxx,n2/2+1,nsizexy,nsizey,x0,x1,npy,&
                     ystable,gyrtable,commcol,nsizez)
        deallocate(x0)
        allocate(x0(n2,nsizexy,nsizez))

        do k = 1, nsizez
          do j = 1, nsizexy 
            do i = 1, n2/2+1
              tmp2(i,j) = x1(i,j,k)
            enddo
            do i = n2/2+2,n2
              tmp2(i,j) = x1(n2-i+2,j,k)
            enddo
          enddo

          call fftlocal_FFT(ksign,scaley,tmp2,n2,nsizexy)

          do j = 1, nsizexy 
            do i = 1, n2
              x0(i,j,k) = tmp2(i,j) 
            enddo
          enddo
        enddo

        deallocate(x1)
        allocate(x1(n3/2+1,nsizexy,nsizeyz))
!        call MPI_BARRIER(commcol,ierr)
        call trans3d3_TRANSP(n2,nsizexy,nsizez,nsizeyz,x0,x1,npx,&
                      xstable,gxrtable,commrow,myidx,n3/2+1)
        deallocate(x0)

        do k = 1, nsizeyz
          do j = 1, nsizexy
            do i = 1, n3/2+1
              tmp3(i,j) = x1(i,j,k) 
            enddo
            do i = n3/2+2,n3
              tmp3(i,j) = x1(n3-i+2,j,k)
            enddo
          enddo

          call fftlocal_FFT(ksign,scalez,tmp3,n3,nsizexy)

          do j = 1, nsizexy
            do i = 1, n3
              grnout(i,j,k) = tmp3(i,j)
            enddo
          enddo
        enddo

        deallocate(x1)

        t_greenf = t_greenf + elapsedtime_Timer(t0)

        end subroutine greenf1tIntnew


!--------------------------------------------------------------------
! Field potential from 2D transverse open, 1D longidutinal periodic 
        !update potential (solving Possion's equation) with isolated 
        ! boundary conditions.
        subroutine update2O1P_FieldQuant(this,source,fldgeom,grid,nxlc,&
          nylc,nzlc,nprocrow,nproccol,nylcr,nzlcr)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(in) :: fldgeom
        double precision, dimension(nxlc,nylc,nzlc), intent(in) :: source
        type (Pgrid2d), intent(in) :: grid
        integer, intent(in) :: nxlc,nylc,nzlc,&
                               nprocrow,nproccol,nylcr,nzlcr
        type (FieldQuant), intent(inout) :: this
        double precision, dimension(3) :: msize
        double precision :: hx, hy, hz, temp
        integer, dimension(0:nprocrow-1) :: ypzstable,pztable
        integer, dimension(0:nproccol-1) :: xpystable,pytable
        double precision , dimension(nxlc,nylcr,nzlcr) :: rho
        integer :: myid,myidz,myidy,&
                   comm2d,commcol,commrow
        integer :: nxpylc2,nypzlc2
        integer :: nsxy1,nsxy2,nsyz1,nsyz2
        integer :: i,j,k,inxglb,inyglb,inzglb,innx,inny,innz
        integer, dimension(2,0:nprocrow-1,0:nproccol-1)::LocalTable
        integer, dimension(3) :: glmshnm
        integer :: jadd,kadd,ierr
        double precision :: t0

        call getpost_Pgrid2d(grid,myid,myidy,myidz)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        call getlctabnm_CompDom(fldgeom,LocalTable)
        
        do i = 0, nprocrow-1
          pztable(i) = LocalTable(1,i,0)
        enddo
        pztable(nprocrow-1) = pztable(nprocrow-1) - 1 !periodic BC.
        do i = 0, nproccol-1
          pytable(i) = LocalTable(2,0,i)
        enddo

        call getmnum_CompDom(fldgeom,glmshnm)
        inxglb = glmshnm(1) 
        inyglb = glmshnm(2)
        inzglb = glmshnm(3)-1  !periodic BC.

        innz = nzlcr
        inny = nylcr
        innx = nxlc

        if(nprocrow.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif
        if(nproccol.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif
        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              rho(i,j,k) = source(i,j+jadd,k+kadd)
            enddo
          enddo
        enddo
        
        ! +1 is from the double precision to complex fft.
        nsxy1 = (inxglb+1)/nproccol
        nsxy2 = (inxglb+1) - nproccol*nsxy1
        do i = 0, nproccol-1
          if(i.le.(nsxy2-1)) then
            xpystable(i) = nsxy1+1
          else
            xpystable(i) = nsxy1
          endif
        enddo

        nsyz1 = 2*inyglb/nprocrow
        nsyz2 = 2*inyglb - nprocrow*nsyz1
        do i = 0, nprocrow-1
          if(i.le.(nsyz2-1)) then
            ypzstable(i) = nsyz1+1
          else
            ypzstable(i) = nsyz1
          endif
        enddo

        nxpylc2 = xpystable(myidy)
        nypzlc2 = ypzstable(myidz)

        call getmsize_CompDom(fldgeom,msize)
        hx = msize(1)
        hy = msize(2)
        hz = msize(3)

        ! Open boundary conditions!
        call open2perd1(innx,inny,innz,rho,hx,hy,hz, &
        nxpylc2,nypzlc2,myidz,myidy,nprocrow,nproccol,commrow,commcol,&
        comm2d,pztable,pytable,ypzstable,xpystable,&
        inxglb,inyglb,inzglb)

        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              this%FieldQ(i,j+jadd,k+kadd) = rho(i,j,k)
            enddo
          enddo
        enddo
        !wrong peiodic condition, works only if nprocrow = 1
        !if(myidz.eq.(nprocrow-1)) then
        !  do j = 1, inny
        !    do i = 1, innx
        !      this%FieldQ(i,j+jadd,innz+1+kadd) = rho(i,j,1+kadd)
        !    enddo
        !  enddo
        !endif

        !open(8,file="phir",status="unknown")
        !do i = 1, innx
        !    write(8,1000)(i-1)*hx,this%FieldQ(i,1,1),&
        !    this%FieldQ(i,1,innz/2),this%FieldQ(i,1,innz)
        !enddo
        !close(8)
        !open(9,file="phithe",status="unknown")
        !do j = 1, inny+1
        !  write(9,1000)(j-1)*hy,this%FieldQ(2,j,2),&
        !  this%FieldQ(3,j,2),this%FieldQ(4,j,2)
        !enddo
        !close(9)
        !open(10,file="phiz",status="unknown")
        !do k = 1, innz+1
        !  write(10,1000)(k-1)*hz,this%FieldQ(1,1,k),&
        !  this%FieldQ(innx/2,1,k),this%FieldQ(innx,1,k)
        !enddo
        !close(10)

1000    format(4(1x,e13.6))

        call MPI_BARRIER(comm2d,ierr)
        t_field = t_field + elapsedtime_Timer(t0)

        end subroutine update2O1P_FieldQuant

        ! Solving Poisson's equation with open BCs.
        subroutine open2perd1(innx,inny,innz,rho,hx,hy,hz,&
           nxpylc2,nypzlc2,myidz,myidy,npz,npy,commrow,commcol,comm2d, &
           pztable,pytable,ypzstable,xpystable,inxglb,&
           inyglb,inzglb)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx,inny,innz,inxglb,inyglb,inzglb
        integer, intent(in) :: nxpylc2,nypzlc2
        integer, intent(in) :: myidz,myidy,npz,npy,commrow,commcol,&
                               comm2d
        integer, dimension(0:npz-1), intent(in) :: pztable,ypzstable
        integer, dimension(0:npy-1), intent(in) :: pytable,xpystable
        double precision, dimension(innx,inny,innz), intent(inout) :: rho
        double precision, intent(in) :: hx, hy, hz
        !double precision, dimension(innx,inny,innz) :: rhosum
        double precision :: scalex,scaley,scalez,sumtest,sumtest2
        integer :: i,j,k,n1,n2,n3,nylc22,nzlc22
        double complex, allocatable, dimension(:,:,:) :: rho2out
        double complex, allocatable, dimension(:,:,:) :: grn
        integer :: ginny,ginnz,ierr

!        rhosum = rho

        n1 = 2*inxglb
        n2 = 2*inyglb
        n3 = inzglb    ! periodic boundary condition.

        nylc22 = nxpylc2
        nzlc22 = nypzlc2

        scalex = 1.0
        scaley = 1.0
        scalez = 1.0

        allocate(rho2out(n3,nylc22,nzlc22))

        call fft3d2_FFT(n1,n2,n3,innz,inny,nylc22,nzlc22,&
                1,scalex,scaley,scalez,rho,ypzstable,pztable,&
        xpystable,pytable,npz,commrow,npy,commcol,comm2d,myidz,myidy, &
        rho2out)

        !c compute FFT of the Green function on the grid:
        ! here the +1 is from the unsymmetry of green function
        ! on double-sized grid.
        ginnz = innz  !periodic boundary condition.
        if(myidy.eq.(npy-1)) then
           ginny = inny + 1
        else
           ginny = inny
        endif
        allocate(grn(n3,nylc22,nzlc22))
        call greenf2(inxglb,inyglb,inzglb,ginnz,ginny,nylc22,nzlc22, &
               hx,hy,hz,myidz,npz,commrow,myidy,npy,commcol,comm2d,&
                  ypzstable,pztable,xpystable,pytable,grn)

        ! multiply transformed charge density and transformed Green 
        ! function:
        do k = 1, nzlc22
          do j = 1, nylc22
            do i = 1, n3
              rho2out(i,j,k) = rho2out(i,j,k)*grn(i,j,k)
            enddo
          enddo
        enddo

        deallocate(grn)

        ! inverse FFT:
        scalex = 1.0/float(n1)
        scaley = 1.0/float(n2)
        scalez = 1.0/float(n3)
        call invfft3d2_FFT(n3,n2,n1,nylc22,nzlc22,inny,innz,&
               -1,scalex,scaley,scalez,rho2out,pztable,ypzstable,&
        pytable,xpystable,npz,commrow,npy,commcol,comm2d,myidz,myidy,&
        rho)

!        sumtest = 0.0
!        sumtest2 = 0.0
!        do k = 1, innz
!          do j = 1, inny
!            do i = 1, innx
!              sumtest2 = sumtest2 + rhosum(i,j,k)
!              sumtest = sumtest + abs(rho(i,j,k)-rhosum(i,j,k))
!            enddo
!          enddo
!        enddo
!        print*,"sumtest-fft: ",sumtest,sumtest2

        deallocate(rho2out)

        end subroutine open2perd1

        ! green function for extended array.
        subroutine greenf2(nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz,&
                  hx,hy,hz,myidx,npx,commrow,myidy,npy,commcol,comm2d,&
                   xstable,xrtable,ystable,yrtable,grnout)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz
        integer, intent(in) :: myidx,myidy,npx,npy,commrow,commcol,&
                               comm2d
        integer, dimension(0:npx-1),intent(in) :: xstable,xrtable
        integer, dimension(0:npy-1),intent(in) :: ystable,yrtable
        double precision, intent(in) :: hx, hy, hz
        double complex, intent (out), &
                       dimension (nz,nsizexy,nsizeyz) :: grnout
        integer, dimension(0:npx-1) :: gxrtable
        integer, dimension(0:npy-1) :: gyrtable
        integer :: i,j,k,kk,ii,jj,iii,jjj,kkk,n1,n2,n3,ns1,ns2
        integer :: nblockx,nblocky,ksign,nxx,ierr,kz
        double precision, dimension (nx+1,nsizey,nsizez) :: grn
        double precision :: scalex,scaley,scalez
        double precision :: t0
        double precision, dimension(2*nx,nsizey) :: tmp1
        double complex, dimension(nx+1,nsizey) :: tmp10
        double complex, dimension(2*ny,nsizexy) :: tmp2
        double complex, dimension(nz,nsizexy) :: tmp3
        double complex, allocatable, dimension(:,:,:) :: x1
        double complex, allocatable, dimension(:,:,:) :: x0
        double precision :: btlmda,sumtmp,tmptmp

        call starttime_Timer(t0)

        btlmda = nz*hz
!        print*,"btlmda: ",btlmda
        gxrtable = xrtable  !periodic BC
        nblockx = 0
        do i = 0, myidx-1
          nblockx = nblockx + gxrtable(i)
        enddo

        gyrtable = yrtable
        gyrtable(npy-1) = yrtable(npy-1) + 1
        nblocky = 0
        do i = 0, myidy-1
          nblocky = nblocky + gyrtable(i)
        enddo
        sumtmp = 0.0
        do k = 1, nsizez
          do j = 1, nsizey
            do i = 1, nx+1
              jj = j + nblocky
              kk = k + nblockx
              if(kk.le.nz/2) then
                kkk = kk - 1
              else
                kkk = kk - nz - 1
              endif
              jjj = jj - 1
              iii = i - 1
              grn(i,j,k) = 0.0
              do kz = -8, 8
                tmptmp=(hx*iii)**2+(hy*jjj)**2+(hz*kkk+kz*btlmda)**2
                if(tmptmp .gt. 1.0e-16) then
                 sumtmp=1.0/sqrt(tmptmp)
                endif
                grn(i,j,k) = grn(i,j,k) + sumtmp
              enddo
            enddo
          enddo
        enddo
        if((myidx.eq.0).and.(myidy.eq.0)) then
          if(nsizez.gt.1) then
            grn(1,1,1) = grn(1,1,2)
          else
            grn(1,1,1) = 1.0
          endif
        endif

        scalex = 1.0
        scaley = 1.0
        scalez = 1.0
        n1 = 2*nx
        n2 = 2*ny
        n3 = nz   !periodic boundary condition.
        ksign = 1

        nxx = n1/2 + 1
        !FFTs along x and y dimensions.
        allocate(x0(n1/2+1,nsizey,nsizez))
        do k = 1, nsizez
          do j = 1, nsizey 
            do i = 1, n1/2+1
              tmp1(i,j) = grn(i,j,k)
            enddo
            do i = n1/2+2, n1
              tmp1(i,j) = grn(n1-i+2,j,k)
            enddo
          enddo

          ! FFTs along x dimensions:
          call fftrclocal_FFT(ksign,scalex,tmp1,n1,nsizey,tmp10)

          do j = 1, nsizey 
            do i = 1, n1/2+1
              x0(i,j,k) = tmp10(i,j)
            enddo
          enddo
        enddo

        allocate(x1(n2/2+1,nsizexy,nsizez))

        ! FFTs along y dimensions:
!        call MPI_BARRIER(commcol,ierr)
        call trans3d_TRANSP(nxx,n2/2+1,nsizexy,nsizey,x0,x1,npy,&
                     ystable,gyrtable,commcol,nsizez)
        deallocate(x0)
        allocate(x0(n2,nsizexy,nsizez))

        do k = 1, nsizez
          do j = 1, nsizexy 
            do i = 1, n2/2+1
              tmp2(i,j) = x1(i,j,k)
            enddo
            do i = n2/2+2,n2
              tmp2(i,j) = x1(n2-i+2,j,k)
            enddo
          enddo

          call fftlocal_FFT(ksign,scaley,tmp2,n2,nsizexy)

          do j = 1, nsizexy 
            do i = 1, n2
              x0(i,j,k) = tmp2(i,j) 
            enddo
          enddo
        enddo

        deallocate(x1)
        allocate(x1(n3,nsizexy,nsizeyz))
!        call MPI_BARRIER(commcol,ierr)
        call trans3d3_TRANSP(n2,nsizexy,nsizez,nsizeyz,x0,x1,npx,&
                      xstable,gxrtable,commrow,myidx,n3)
        deallocate(x0)

        do k = 1, nsizeyz
          do j = 1, nsizexy
            do i = 1, n3
              tmp3(i,j) = x1(i,j,k) 
            enddo
          enddo

          call fftlocal_FFT(ksign,scalez,tmp3,n3,nsizexy)

          do j = 1, nsizexy
            do i = 1, n3
              grnout(i,j,k) = tmp3(i,j)
            enddo
          enddo
        enddo

        deallocate(x1)

        t_greenf = t_greenf + elapsedtime_Timer(t0)

        end subroutine greenf2
!-------------------------------------------------------------------------

        subroutine update3_FieldQuant(this,source,fldgeom,grid,nxlc,&
          nylc,nzlc,nprocrow,nproccol,nylcr,nzlcr,&
          besscoef,bessnorm,gml,modth,nmod)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(in) :: fldgeom
        type (FieldQuant), intent(inout) :: this
        type (Pgrid2d), intent(in) :: grid
        double precision, dimension(nxlc,nylc,nzlc), intent(in) :: source
        double precision, dimension(nxlc,nxlc,nmod), intent(in) :: besscoef
        double precision, dimension(nxlc,nmod), intent(in) :: bessnorm,gml
        integer, dimension(nylcr), intent(in) :: modth
        integer, intent(in) :: nxlc,nylc,nzlc,nmod,&
                               nprocrow,nproccol,nylcr,nzlcr
        double precision, dimension(3) :: msize
        double precision :: hx, hy, hz, temp
        integer, dimension(0:nprocrow-1) :: pztable,xpzstable,pzdisp
        integer, dimension(0:nproccol-1) :: xpystable,pytable
        double precision , dimension(nxlc,nylcr,nzlcr) :: rho
        integer :: myid,myidz,myidy,comm2d,commcol,commrow
        integer :: nxpylc,nxpzlc,nsxy1,nsxy2,nsxz1,nsxz2
        integer :: i,j,k,inxglb,inyglb,inzglb,innx,inny,innz
        integer, dimension(2,0:nprocrow-1,0:nproccol-1)::LocalTable
        integer, dimension(3) :: glmshnm
        integer :: jadd,kadd,ierr,ksign,kk
        double precision :: t0,pi,fourpi,testsum,totsum,testrho,totrho,rad,scaley
        double precision :: sec0,sec1,sec,tval,sumtest
        double precision :: sumrho0,sumrho0tot,sumrho1,sumrho1tot,sumrho2,&
        sumrho2tot,sumrho3,sumrho3tot,sumrho4,sumrho4tot,sumrho5,sumrho5tot
        double precision :: sumbess,sumbesstot,sumbn,sumbntot,sumod,sumodtot

        call getpost_Pgrid2d(grid,myid,myidy,myidz)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        call getlctabnm_CompDom(fldgeom,LocalTable)
        
        do i = 0, nprocrow-1
          pztable(i) = LocalTable(1,i,0)
        enddo
        pztable(nprocrow-1) = pztable(nprocrow-1)
        pzdisp(0) = 0
        do i = 1, nprocrow-1
          pzdisp(i) = pzdisp(i-1)+pztable(i-1)
        enddo
        do i = 0, nproccol-1
          pytable(i) = LocalTable(2,0,i)
        enddo
        pytable(nproccol-1) = pytable(nproccol-1) - 1 !periodic BC.

        call getmnum_CompDom(fldgeom,glmshnm)
        inxglb = glmshnm(1) 
        inyglb = glmshnm(2)-1  !periodic BC.
        inzglb = glmshnm(3)

        innz = nzlcr
        inny = nylcr
        innx = nxlc

        !used for transpose between serial x and parallel y.
        nsxy1 = inxglb/nproccol
        nsxy2 = inxglb - nproccol*nsxy1
        do i = 0, nproccol-1
          if(i.le.(nsxy2-1)) then
            xpystable(i) = nsxy1+1
          else
            xpystable(i) = nsxy1
          endif
        enddo
        nxpylc = xpystable(myidy)

        !used for transpose between serial x and parallel z.
        nsxz1 = inxglb/nprocrow  !periodic BC
        nsxz2 = inxglb - nprocrow*nsxz1
        do i = 0, nprocrow-1
          if(i.le.(nsxz2-1)) then
            xpzstable(i) = nsxz1+1
          else
            xpzstable(i) = nsxz1
          endif
        enddo
        nxpzlc = xpzstable(myidz)

        call getmsize_CompDom(fldgeom,msize)
        hx = msize(1)
        hy = msize(2)
        hz = msize(3)

        if(nprocrow.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif
        if(nproccol.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif

        do k = 1, innz
          do j = 1, inny
            rho(1,j,k) = source(1,j+jadd,k+kadd)/(0.25*hx*hx*hy*hz)
            do i = 2, innx-1
              rho(i,j,k) = source(i,j+jadd,k+kadd)/((i-1)*hx*hx*hy*hz)*&
                           1.0
            enddo
            rho(innx,j,k) = source(innx,j+jadd,k+kadd)/&
            (0.25*(2*innx-2.5)*hx*hx*hy*hz)
          enddo
        enddo

        rad = (innx-1)*hx
        pi = 2*asin(1.0)
        !do k = 1, innz
        !  kk = k + pzdisp(myidz)
        !  do j = 1, inny
        !    do i = 1, innx
        !      if((i-1)*hx.le.0.5*rad) then
        !         !rho(i,j,k) = (i-1)*hx*(cos((j-1)*hy))**2*24.0/&
        !         !             (pi*rad**3)
        !         !rho(i,j,k) = 10.0
        !         rho(i,j,k) = exp(-10*((i-1)*hx+(kk-1)*hz)**2)
        !      else
        !         rho(i,j,k) = 0.0
        !      endif
        !    enddo
        !  enddo
        !enddo

!        sumrho0 = 0.0
!        do k = 1,innz
!          do j = 1, inny
!            do i = 1, innx
!              sumrho0 = sumrho0 + rho(i,j,k)
!            enddo
!          enddo
!        enddo
!        call MPI_ALLREDUCE(sumrho0,sumrho0tot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        print*,"sumrho0: ",sumrho0tot,innx,inny,innz,myidz,myidy

        !sec0 = MPI_WTIME()
        ksign = 1
        scaley = 1.0 
        call fft3d3_FFT(inxglb,inyglb,innz,inny,nxpylc,&
        ksign,scaley,rho,xpystable,pytable,nproccol,commcol,comm2d)

!        sumrho1 = 0.0
!        do k = 1,innz
!          do j = 1, inny
!            do i = 1, innx
!              sumrho1 = sumrho1 + rho(i,j,k)
!            enddo
!          enddo
!        enddo
!        call MPI_ALLREDUCE(sumrho1,sumrho1tot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        print*,"sumrho1: ",sumrho1tot

        call Bessr_Bessel(inxglb,inny,innz,rho,besscoef,&
                         bessnorm,hx,modth,myidy,nmod)

!        sumrho2 = 0.0
!        do k = 1,innz
!          do j = 1, inny
!            do i = 1, innx
!              sumrho2 = sumrho2 + rho(i,j,k)
!            enddo
!          enddo
!        enddo
!        call MPI_ALLREDUCE(sumrho2,sumrho2tot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        print*,"sumrho2: ",sumrho2tot
!        sumbess = 0.0
!        sumbn = 0.0
!        if(myidy.eq.0) then
!
!        do k = 1,nmod
!          do j = 1, inxglb
!            sumbn = sumbn + bessnorm(j,k)
!            do i = 1, inxglb
!              sumbess = sumbess + besscoef(i,j,k)
!            enddo
!          enddo
!        enddo
!
!        else
!
!        do k = 1,nmod
!          do j = 1, inxglb
!            sumbn = sumbn + bessnorm(j,k)
!            do i = 1, inxglb
!              sumbess = sumbess + besscoef(i,j,k)
!            enddo
!          enddo
!        enddo
!
!        endif
!
!        sumod = 0.0
!        do i = 1, inny
!          sumod = sumod + modth(i)
!        enddo
!        call MPI_ALLREDUCE(sumbess,sumbesstot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        call MPI_ALLREDUCE(sumbn,sumbntot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        call MPI_ALLREDUCE(sumod,sumodtot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        print*,"sumbess: ",sumbesstot,sumbntot,sumodtot,nmod

        call Gaussz3(inxglb,inzglb,innz,nxpzlc,inny,rho,&
        xpzstable,pztable,nprocrow,commrow,myidz,myidy,gml,modth,hz,&
        nmod)

!        sumrho3 = 0.0
!        do k = 1,innz
!          do j = 1, inny
!            do i = 1, innx
!              sumrho3 = sumrho3 + rho(i,j,k)
!            enddo
!          enddo
!        enddo
!        call MPI_ALLREDUCE(sumrho3,sumrho3tot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        print*,"sumrho3: ",sumrho3tot

        call BessrI_Bessel(inxglb,inny,innz,rho,besscoef,modth,&
                               myidy,nmod)

!        sumrho4 = 0.0
!        do k = 1,innz
!          do j = 1, inny
!            do i = 1, innx
!              sumrho4 = sumrho4 + rho(i,j,k)
!            enddo
!          enddo
!        enddo
!        call MPI_ALLREDUCE(sumrho4,sumrho4tot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        print*,"sumrho4: ",sumrho4tot

        ksign = -1
        scaley = 1.0/float(inyglb) 
        call fft3d3_FFT(inxglb,inyglb,innz,inny,nxpylc,&
        ksign,scaley,rho,xpystable,pytable,nproccol,commcol,comm2d)
        !sec1 = MPI_WTIME()
        !sec = sec1 - sec0
        !call MPI_ALLREDUCE(sec,tval,1,MPI_REAL,MPI_MAX, &
        !                     MPI_COMM_WORLD,ierr)
        !print*,"Field time: ",tval

!        sumrho5 = 0.0
!        do k = 1,innz
!          do j = 1, inny
!            do i = 1, innx
!              sumrho5 = sumrho5 + rho(i,j,k)
!            enddo
!          enddo
!        enddo
!        call MPI_ALLREDUCE(sumrho5,sumrho5tot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        print*,"sumrho5: ",sumrho5tot

        !testrho = 0.0
        !do k = 1, innz
        !  do j = 1, inny
        !    do i = 1, innx
        !      testrho = testrho + rho(i,j,k)
        !    enddo
        !  enddo
        !enddo
        !call MPI_ALLREDUCE(testrho,totrho,1,MPI_REAL,MPI_SUM, &
        !                     MPI_COMM_WORLD,ierr)
        !print*,"totrho: ",totrho
        
        fourpi = 4*pi
        ! 1/(4*pi) is due to the multiplication of 4pi in Rob's code.
        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              this%FieldQ(i,j+jadd,k+kadd) = rho(i,j,k)*fourpi
            enddo
          enddo
        enddo
        !wrong peiodic condition, works only if nprocrow = 1
!        if(myidy.eq.(nproccol-1)) then
!          do k = 1, innz
!            do i = 1, innx
!              this%FieldQ(i,inny+1+jadd,k+kadd)=rho(i,1,k)*fourpi
!            enddo
!          enddo
!        endif
!        sumtest = 0.0
!        do k = 1, innz
!          do j = 1, inny
!            do i = 1, innx
!              sumtest = sumtest + rho(i,j,k)
!            enddo
!          enddo
!        enddo
!        print*,"sumtest, potential: ",sumtest

        call MPI_BARRIER(comm2d,ierr)
        t_field = t_field + elapsedtime_Timer(t0)

        end subroutine update3_FieldQuant

        subroutine Gaussz3(nx,nz,nsizez,nsizexz,nsizey,x,&
        xstable,xrtable,nprocrow,commrow,myidx,myidy,gm,modth,hz,nmod)
        implicit none
        include 'mpif.h'
        integer,intent(in) :: nx,nz,nsizez,nsizexz,nsizey,nmod
        integer,intent(in) :: nprocrow,commrow
        integer,intent(in) :: myidx,myidy
        integer,dimension(0:nprocrow-1),intent(in) :: xstable,xrtable
        double precision, dimension(nx,nsizey,nsizez), intent(inout) :: x
        double precision, dimension(nx,nmod), intent(in) :: gm
        double precision, intent(in) :: hz
        integer, dimension(nsizey), intent(in) :: modth
        integer,dimension(0:nprocrow-1) :: rdisp
        integer :: i,j,k,jm
        double precision :: t0,tmp1,tmp2
        double precision, dimension(nz) :: aa,cc
        double precision, dimension(nz,nsizey,nsizexz) :: bb
        double precision, allocatable, dimension(:,:,:) :: x0

        call starttime_Timer(t0)

        rdisp(0) = 0
        do i = 1, nprocrow-1
          rdisp(i) = rdisp(i-1) + xstable(i-1)
        enddo

        allocate(x0(nz,nsizey,nsizexz))
        !transpose z direction into serial direction:
        call trans3d3r_TRANSP(nx,nsizey,nsizez,nsizexz,x,x0,nprocrow,&
                      xstable,xrtable,commrow,myidx,nz)

        do i = 1, nz -1
          aa(i) = 1.0
        enddo
        aa(nz) = 0.0
        cc(1) = 0.0
        do i = 2, nz 
          cc(i) = 1.0
        enddo
        do j = 1, nsizey
          if(myidy.ne.0) then
            jm = modth(j)-modth(1) + 1
          else
!            jm = modth(j)-modth(1) + 1
!            if(j.eq.2) jm = j
            if(j.le.2) then
              jm = j
            else
              jm = modth(j)-modth(1) + 2
            endif
          endif
          do k = 1, nsizexz
            tmp1 = hz*hz*gm(k+rdisp(myidx),jm)*gm(k+rdisp(myidx),jm)
            tmp2 = exp(-hz*gm(k+rdisp(myidx),jm))
            bb(1,j,k) = -2 - tmp1 + tmp2
            do i = 2, nz-1
              bb(i,j,k) = -2 - tmp1
            enddo
            bb(nz,j,k) = -2 - tmp1 + tmp2
          enddo
        enddo
       
        !eliminate the lower off-diagnal terms.
        do k = 1, nsizexz
          do j = 1, nsizey
            do i = 2, nz
              bb(i,j,k) = bb(i,j,k) - aa(i-1)*cc(i)/bb(i-1,j,k)
            enddo
          enddo
        enddo
        !the negative sign due to the "-" in the Poisson's equation. 
        do k = 1, nsizexz
          do j = 1, nsizey
            do i = 1, nz
              x0(i,j,k) = -x0(i,j,k)*hz*hz
            enddo
          enddo
        enddo

        !forward substitute for Gauss Elemination.
        do k = 1, nsizexz
          do j = 1, nsizey
            do i = 2, nz
              x0(i,j,k) = x0(i,j,k) - x0(i-1,j,k)*cc(i)/bb(i-1,j,k)
            enddo
          enddo
        enddo

        !backward substitute
        do k = 1, nsizexz
          do j = 1, nsizey
            x0(nz,j,k) = x0(nz,j,k)/bb(nz,j,k)
            do i = nz-1, 1, -1
              x0(i,j,k) = (x0(i,j,k) - x0(i+1,j,k)*aa(i))/bb(i,j,k)
            enddo
          enddo
        enddo

        call trans3d3r_TRANSP(nz,nsizey,nsizexz,nsizez,x0,x,nprocrow,&
                      xrtable,xstable,commrow,myidx,nx)

        deallocate(x0)

        t_fft2dhpf = t_fft2dhpf + elapsedtime_Timer(t0)

        return
        end subroutine Gaussz3

!-------------------------------------------------------------------------
! solving the 3D Poisson equation in the cylindrical coordinate system
! subject to the Dirchet boundary conditions at both ends, z=0, z=Zn.
! here, we have used a Fourier expansion along theta and a Bessel expansion
! also radial direction.
! we need to incorporate the boundary conditions into the source term
! at the both ends.
        subroutine update32_FieldQuant(this,source,fldgeom,grid,nxlc,&
          nylc,nzlc,nprocrow,nproccol,nylcr,nzlcr,&
          besscoef,bessnorm,gml,modth,nmod,phi1,phin)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(in) :: fldgeom
        type (FieldQuant), intent(inout) :: this
        type (Pgrid2d), intent(in) :: grid
        double precision, dimension(nxlc,nylc,nzlc), intent(in) :: source
        double precision, dimension(nxlc,nxlc,nmod), intent(in) :: besscoef
        double precision, dimension(nxlc,nmod), intent(in) :: bessnorm,gml
        double precision, dimension(nxlc,nylcr), intent(in) :: phi1,phin
        integer, dimension(nylcr), intent(in) :: modth
        integer, intent(in) :: nxlc,nylc,nzlc,nmod,&
                               nprocrow,nproccol,nylcr,nzlcr
        double precision, dimension(3) :: msize
        double precision :: hx, hy, hz, temp
        integer, dimension(0:nprocrow-1) :: pztable,xpzstable,pzdisp
        integer, dimension(0:nproccol-1) :: xpystable,pytable
        double precision , dimension(nxlc,nylcr,nzlcr) :: rho
        integer :: myid,myidz,myidy,comm2d,commcol,commrow
        integer :: nxpylc,nxpzlc,nsxy1,nsxy2,nsxz1,nsxz2
        integer :: i,j,k,inxglb,inyglb,inzglb,innx,inny,innz
        integer, dimension(2,0:nprocrow-1,0:nproccol-1)::LocalTable
        integer, dimension(3) :: glmshnm
        integer :: jadd,kadd,ierr,ksign,kk
        double precision :: t0,pi,fourpi,testsum,totsum,testrho,totrho,rad,scaley
        double precision :: sec0,sec1,sec,tval,sumtest
        double precision :: sumrho0,sumrho0tot,sumrho1,sumrho1tot,sumrho2,&
        sumrho2tot,sumrho3,sumrho3tot,sumrho4,sumrho4tot,sumrho5,sumrho5tot
        double precision :: sumbess,sumbesstot,sumbn,sumbntot,sumod,sumodtot

        call getpost_Pgrid2d(grid,myid,myidy,myidz)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        call getlctabnm_CompDom(fldgeom,LocalTable)
        
        do i = 0, nprocrow-1
          pztable(i) = LocalTable(1,i,0)
        enddo
        pztable(nprocrow-1) = pztable(nprocrow-1)
        pzdisp(0) = 0
        do i = 1, nprocrow-1
          pzdisp(i) = pzdisp(i-1)+pztable(i-1)
        enddo
        do i = 0, nproccol-1
          pytable(i) = LocalTable(2,0,i)
        enddo
        pytable(nproccol-1) = pytable(nproccol-1) - 1 !periodic BC.

        call getmnum_CompDom(fldgeom,glmshnm)
        inxglb = glmshnm(1) 
        inyglb = glmshnm(2)-1  !periodic BC.
        inzglb = glmshnm(3)

        innz = nzlcr
        inny = nylcr
        innx = nxlc

        !used for transpose between serial x and parallel y.
        nsxy1 = inxglb/nproccol
        nsxy2 = inxglb - nproccol*nsxy1
        do i = 0, nproccol-1
          if(i.le.(nsxy2-1)) then
            xpystable(i) = nsxy1+1
          else
            xpystable(i) = nsxy1
          endif
        enddo
        nxpylc = xpystable(myidy)

        !used for transpose between serial x and parallel z.
        nsxz1 = inxglb/nprocrow  !periodic BC
        nsxz2 = inxglb - nprocrow*nsxz1
        do i = 0, nprocrow-1
          if(i.le.(nsxz2-1)) then
            xpzstable(i) = nsxz1+1
          else
            xpzstable(i) = nsxz1
          endif
        enddo
        nxpzlc = xpzstable(myidz)

        call getmsize_CompDom(fldgeom,msize)
        hx = msize(1)
        hy = msize(2)
        hz = msize(3)

        if(nprocrow.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif
        if(nproccol.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif

        do k = 1, innz
          do j = 1, inny
            rho(1,j,k) = source(1,j+jadd,k+kadd)/(0.125*hx*hx*hy*hz)
            do i = 2, innx-1
              rho(i,j,k) = source(i,j+jadd,k+kadd)/((i-1)*hx*hx*hy*hz)*&
                           1.0
            enddo
            rho(innx,j,k) = source(innx,j+jadd,k+kadd)/&
            (0.25*(2*innx-2.5)*hx*hx*hy*hz)
          enddo
        enddo

        !set Dirchet boundary conditions.
        if(myidz.eq.0) then
          do j = 1, inny
            do i = 1, innx
              rho(i,j,1) = phi1(i,j)
            enddo
          enddo
        endif
        if(myidz.eq.nprocrow) then
          do j = 1, inny
            do i = 1, innx
              rho(i,j,innz) = phin(i,j)
            enddo
          enddo
        endif

        rad = (innx-1)*hx
        pi = 2*asin(1.0)
        !do k = 1, innz
        !  kk = k + pzdisp(myidz)
        !  do j = 1, inny
        !    do i = 1, innx
        !      if((i-1)*hx.le.0.5*rad) then
        !         !rho(i,j,k) = (i-1)*hx*(cos((j-1)*hy))**2*24.0/&
        !         !             (pi*rad**3)
        !         !rho(i,j,k) = 10.0
        !         rho(i,j,k) = exp(-10*((i-1)*hx+(kk-1)*hz)**2)
        !      else
        !         rho(i,j,k) = 0.0
        !      endif
        !    enddo
        !  enddo
        !enddo

!        sumrho0 = 0.0
!        do k = 1,innz
!          do j = 1, inny
!            do i = 1, innx
!              sumrho0 = sumrho0 + rho(i,j,k)
!            enddo
!          enddo
!        enddo
!        call MPI_ALLREDUCE(sumrho0,sumrho0tot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        print*,"sumrho0: ",sumrho0tot,innx,inny,innz,myidz,myidy

        !sec0 = MPI_WTIME()
        ksign = 1
        scaley = 1.0 
        call fft3d3_FFT(inxglb,inyglb,innz,inny,nxpylc,&
        ksign,scaley,rho,xpystable,pytable,nproccol,commcol,comm2d)

!        sumrho1 = 0.0
!        do k = 1,innz
!          do j = 1, inny
!            do i = 1, innx
!              sumrho1 = sumrho1 + rho(i,j,k)
!            enddo
!          enddo
!        enddo
!        call MPI_ALLREDUCE(sumrho1,sumrho1tot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        print*,"sumrho1: ",sumrho1tot

        call Bessr_Bessel(inxglb,inny,innz,rho,besscoef,&
                         bessnorm,hx,modth,myidy,nmod)

!        sumrho2 = 0.0
!        do k = 1,innz
!          do j = 1, inny
!            do i = 1, innx
!              sumrho2 = sumrho2 + rho(i,j,k)
!            enddo
!          enddo
!        enddo
!        call MPI_ALLREDUCE(sumrho2,sumrho2tot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        print*,"sumrho2: ",sumrho2tot
!        sumbess = 0.0
!        sumbn = 0.0
!        if(myidy.eq.0) then
!
!        do k = 1,nmod
!          do j = 1, inxglb
!            sumbn = sumbn + bessnorm(j,k)
!            do i = 1, inxglb
!              sumbess = sumbess + besscoef(i,j,k)
!            enddo
!          enddo
!        enddo
!
!        else
!
!        do k = 1,nmod
!          do j = 1, inxglb
!            sumbn = sumbn + bessnorm(j,k)
!            do i = 1, inxglb
!              sumbess = sumbess + besscoef(i,j,k)
!            enddo
!          enddo
!        enddo
!
!        endif
!
!        sumod = 0.0
!        do i = 1, inny
!          sumod = sumod + modth(i)
!        enddo
!        call MPI_ALLREDUCE(sumbess,sumbesstot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        call MPI_ALLREDUCE(sumbn,sumbntot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        call MPI_ALLREDUCE(sumod,sumodtot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        print*,"sumbess: ",sumbesstot,sumbntot,sumodtot,nmod

        call Gaussz3(inxglb,inzglb,innz,nxpzlc,inny,rho,&
        xpzstable,pztable,nprocrow,commrow,myidz,myidy,gml,modth,hz,&
        nmod)

!        sumrho3 = 0.0
!        do k = 1,innz
!          do j = 1, inny
!            do i = 1, innx
!              sumrho3 = sumrho3 + rho(i,j,k)
!            enddo
!          enddo
!        enddo
!        call MPI_ALLREDUCE(sumrho3,sumrho3tot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        print*,"sumrho3: ",sumrho3tot

        call BessrI_Bessel(inxglb,inny,innz,rho,besscoef,modth,&
                               myidy,nmod)

!        sumrho4 = 0.0
!        do k = 1,innz
!          do j = 1, inny
!            do i = 1, innx
!              sumrho4 = sumrho4 + rho(i,j,k)
!            enddo
!          enddo
!        enddo
!        call MPI_ALLREDUCE(sumrho4,sumrho4tot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        print*,"sumrho4: ",sumrho4tot

        ksign = -1
        scaley = 1.0/float(inyglb) 
        call fft3d3_FFT(inxglb,inyglb,innz,inny,nxpylc,&
        ksign,scaley,rho,xpystable,pytable,nproccol,commcol,comm2d)
        !sec1 = MPI_WTIME()
        !sec = sec1 - sec0
        !call MPI_ALLREDUCE(sec,tval,1,MPI_REAL,MPI_MAX, &
        !                     MPI_COMM_WORLD,ierr)
        !print*,"Field time: ",tval

!        sumrho5 = 0.0
!        do k = 1,innz
!          do j = 1, inny
!            do i = 1, innx
!              sumrho5 = sumrho5 + rho(i,j,k)
!            enddo
!          enddo
!        enddo
!        call MPI_ALLREDUCE(sumrho5,sumrho5tot,1,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,MPI_COMM_WORLD,ierr)
!        print*,"sumrho5: ",sumrho5tot

        !testrho = 0.0
        !do k = 1, innz
        !  do j = 1, inny
        !    do i = 1, innx
        !      testrho = testrho + rho(i,j,k)
        !    enddo
        !  enddo
        !enddo
        !call MPI_ALLREDUCE(testrho,totrho,1,MPI_REAL,MPI_SUM, &
        !                     MPI_COMM_WORLD,ierr)
        !print*,"totrho: ",totrho
        
        fourpi = 4*pi
        ! 1/(4*pi) is due to the multiplication of 4pi in Rob's code.
        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              this%FieldQ(i,j+jadd,k+kadd) = rho(i,j,k)*fourpi
            enddo
          enddo
        enddo
        !wrong peiodic condition, works only if nprocrow = 1
!        if(myidy.eq.(nproccol-1)) then
!          do k = 1, innz
!            do i = 1, innx
!              this%FieldQ(i,inny+1+jadd,k+kadd)=rho(i,1,k)*fourpi
!            enddo
!          enddo
!        endif
!        sumtest = 0.0
!        do k = 1, innz
!          do j = 1, inny
!            do i = 1, innx
!              sumtest = sumtest + rho(i,j,k)
!            enddo
!          enddo
!        enddo
!        print*,"sumtest, potential: ",sumtest

        call MPI_BARRIER(comm2d,ierr)
        t_field = t_field + elapsedtime_Timer(t0)

        end subroutine update32_FieldQuant

        !Solving the tradiagonal matrix using Gaussian elimination.
        !here, the boundary conditions are applied at both ends,
        !i.e. i = 1, i = Nz. Hence, the x in the global coordinate
        !z = 0, z= Z_n has to contain the boundary condition.
        subroutine Gaussz32(nx,nz,nsizez,nsizexz,nsizey,x,&
        xstable,xrtable,nprocrow,commrow,myidx,myidy,gm,modth,hz,nmod)
        implicit none
        include 'mpif.h'
        integer,intent(in) :: nx,nz,nsizez,nsizexz,nsizey,nmod
        integer,intent(in) :: nprocrow,commrow
        integer,intent(in) :: myidx,myidy
        integer,dimension(0:nprocrow-1),intent(in) :: xstable,xrtable
        double precision, dimension(nx,nsizey,nsizez), intent(inout) :: x
        double precision, dimension(nx,nmod), intent(in) :: gm
        double precision, intent(in) :: hz
        integer, dimension(nsizey), intent(in) :: modth
        integer,dimension(0:nprocrow-1) :: rdisp
        integer :: i,j,k,jm
        double precision :: t0,tmp1,tmp2
        double precision, dimension(nz) :: aa,cc
        double precision, dimension(nz,nsizey,nsizexz) :: bb
        double precision, allocatable, dimension(:,:,:) :: x0

        call starttime_Timer(t0)

        rdisp(0) = 0
        do i = 1, nprocrow-1
          rdisp(i) = rdisp(i-1) + xstable(i-1)
        enddo

        allocate(x0(nz,nsizey,nsizexz))
        !transpose z direction into serial direction:
        call trans3d3r_TRANSP(nx,nsizey,nsizez,nsizexz,x,x0,nprocrow,&
                      xstable,xrtable,commrow,myidx,nz)

        do i = 1, nz -1
          aa(i) = 1.0
        enddo
        aa(nz) = 0.0
        cc(1) = 0.0
        do i = 2, nz 
          cc(i) = 1.0
        enddo
        do j = 1, nsizey
          if(myidy.ne.0) then
            jm = modth(j)-modth(1) + 1
          else
!            jm = modth(j)-modth(1) + 1
!            if(j.eq.2) jm = j
            if(j.le.2) then
              jm = j
            else
              jm = modth(j)-modth(1) + 2
            endif
          endif
          do k = 1, nsizexz
            tmp1 = hz*hz*gm(k+rdisp(myidx),jm)*gm(k+rdisp(myidx),jm)
            tmp2 = 0.0
            bb(1,j,k) = -2 - tmp1 
            do i = 2, nz-1
              bb(i,j,k) = -2 - tmp1
            enddo
            bb(nz,j,k) = -2 - tmp1
          enddo
        enddo
       
        !eliminate the lower off-diagnal terms.
        do k = 1, nsizexz
          do j = 1, nsizey
            do i = 3, nz-1
              bb(i,j,k) = bb(i,j,k) - aa(i-1)*cc(i)/bb(i-1,j,k)
            enddo
          enddo
        enddo
        !the negative sign due to the "-" in the Poisson's equation. 
        !put in the Dirchet boundary condition at both ends
        do k = 1, nsizexz
          do j = 1, nsizey
            x0(2,j,k) = -x0(2,j,k)*hz*hz -x0(1,j,k)
            do i = 3, nz-2
              x0(i,j,k) = -x0(i,j,k)*hz*hz
            enddo
            x0(nz-1,j,k) = -x0(nz-1,j,k)*hz*hz -x0(nz,j,k)
          enddo
        enddo

        !forward substitute for Gauss Elemination.
        do k = 1, nsizexz
          do j = 1, nsizey
            do i = 3, nz-1
              x0(i,j,k) = x0(i,j,k) - x0(i-1,j,k)*cc(i)/bb(i-1,j,k)
            enddo
          enddo
        enddo

        !backward substitute
        do k = 1, nsizexz
          do j = 1, nsizey
            x0(nz-1,j,k) = x0(nz-1,j,k)/bb(nz-1,j,k)
            do i = nz-2, 2, -1
              x0(i,j,k) = (x0(i,j,k) - x0(i+1,j,k)*aa(i))/bb(i,j,k)
            enddo
          enddo
        enddo

        call trans3d3r_TRANSP(nz,nsizey,nsizexz,nsizez,x0,x,nprocrow,&
                      xrtable,xstable,commrow,myidx,nx)

        deallocate(x0)

        t_fft2dhpf = t_fft2dhpf + elapsedtime_Timer(t0)

        return
        end subroutine Gaussz32

!-------------------------------------------------------------------------

!-------------------------------------------------------------------------

        subroutine update4_FieldQuant(this,source,fldgeom,grid,nxlc,&
          nylc,nzlc,nprocrow,nproccol,nylcr,nzlcr,&
          gblam)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(in) :: fldgeom
        double precision, dimension(nxlc,nylc,nzlc), intent(in) :: source
        type (Pgrid2d), intent(in) :: grid
        integer, intent(in) :: nxlc,nylc,nzlc,&
                               nprocrow,nproccol,nylcr,nzlcr
        type (FieldQuant), intent(inout) :: this
        double precision, intent(in) :: gblam
        double precision, dimension(3) :: msize
        double precision :: hx, hy, hz, temp
        integer, dimension(0:nprocrow-1) :: ypzstable,pztable,xpzstable
        integer, dimension(0:nproccol-1) :: xpystable,pytable,zpystable
        double precision , dimension(nxlc,nylcr,nzlcr) :: rho
        integer :: myid,myidz,myidy,&
                   comm2d,commcol,commrow
        integer :: nxpylc,nypzlc,nxpzlc,nzpylc
        integer :: nsxy1,nsxy2,nsyz1,nsyz2
        integer :: i,j,k,inxglb,inyglb,inzglb,innx,inny,innz
        integer, dimension(2,0:nprocrow-1,0:nproccol-1)::LocalTable
        integer, dimension(3) :: glmshnm
        integer :: jadd,kadd,ierr
        double precision :: t0,pi,fourpi,testsum,totsum,testrho,totrho

        call getpost_Pgrid2d(grid,myid,myidy,myidz)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        call getlctabnm_CompDom(fldgeom,LocalTable)
        
        do i = 0, nprocrow-1
          pztable(i) = LocalTable(1,i,0)
        enddo
        pztable(nprocrow-1) = pztable(nprocrow-1) - 1 !periodic BC.
        do i = 0, nproccol-1
          pytable(i) = LocalTable(2,0,i)
        enddo
        pytable(nproccol-1) = pytable(nproccol-1) - 1 !periodic BC.

        call getmnum_CompDom(fldgeom,glmshnm)
        inxglb = glmshnm(1) 
        inyglb = glmshnm(2)-1  !periodic BC.
        inzglb = glmshnm(3)-1  !periodic BC.

        innz = nzlcr
        inny = nylcr
        innx = nxlc

        nsxy1 = inxglb/nproccol
        nsxy2 = inxglb - nproccol*nsxy1
        do i = 0, nproccol-1
          if(i.le.(nsxy2-1)) then
            xpystable(i) = nsxy1+1
          else
            xpystable(i) = nsxy1
          endif
        enddo
        nxpylc = xpystable(myidy)

        ! /2 + 1 is from the double precision to complex fft.
        nsyz1 = (inyglb/2+1)/nprocrow  !periodic BC
        nsyz2 = (inyglb/2+1) - nprocrow*nsyz1
        do i = 0, nprocrow-1
          if(i.le.(nsyz2-1)) then
            ypzstable(i) = nsyz1+1
          else
            ypzstable(i) = nsyz1
          endif
        enddo
        nypzlc = ypzstable(myidz)

        nsxy1 = inxglb/nprocrow
        nsxy2 = inxglb - nprocrow*nsxy1
        do i = 0, nprocrow-1
          if(i.le.(nsxy2-1)) then
            xpzstable(i) = nsxy1+1
          else
            xpzstable(i) = nsxy1
          endif
        enddo
        nxpzlc = xpzstable(myidz)

        nsyz1 = inzglb/nproccol  !periodic BC
        nsyz2 = inzglb - nproccol*nsyz1
        do i = 0, nproccol-1
          if(i.le.(nsyz2-1)) then
            zpystable(i) = nsyz1+1
          else
            zpystable(i) = nsyz1
          endif
        enddo
        nzpylc = zpystable(myidy)

        call getmsize_CompDom(fldgeom,msize)
        !hx = msize(1)
        hx = msize(1)
        hy = msize(2)
        hz = msize(3)
        !if(myid.eq.0) then
        !  print*,"hx: ",hx,hy,hz,(inxglb-1)*hx
        !endif

        if(nprocrow.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif
        if(nproccol.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif

        !from the "source(i,j,k): # of particles on the grid" to get the
        !density on the grid.
        do k = 1, innz
          do j = 1, inny
            rho(1,j,k) = source(1,j+jadd,k+kadd)/(0.25*hx*hx*hy*hz)
            do i = 2, innx-1
              rho(i,j,k) = source(i,j+jadd,k+kadd)/((i-1)*hx*hx*hy*hz)*&
                           1.0
            enddo
            rho(innx,j,k) = source(innx,j+jadd,k+kadd)/&
            (0.25*(2*innx-2.5)*hx*hx*hy*hz)
          enddo
        enddo

        ! Finite boundary conditions!
        call finiteBC(innx,inny,innz,rho,hx,hy,hz,nxpylc,nypzlc,&
        nxpzlc,nzpylc,myidz,myidy,nprocrow,nproccol,commrow,commcol,&
        comm2d,pztable,pytable,ypzstable,xpystable,xpzstable,zpystable,&
        inxglb,inyglb,inzglb,gblam)

        pi = 2*asin(1.0d0)
        fourpi = 4*pi
        ! 1/(4*pi) is due to the multiplication of 4pi in Rob's code.
        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              this%FieldQ(i,j+jadd,k+kadd) = rho(i,j,k)*fourpi
            enddo
          enddo
        enddo
        !wrong peiodic condition, works only if nprocrow = 1
!        if(myidz.eq.(nprocrow-1)) then
!          do j = 1, inny
!            do i = 1, innx
!              this%FieldQ(i,j+jadd,innz+1+kadd)=rho(i,j,1)*fourpi
!            enddo
!          enddo
!        endif
!        if(myidy.eq.(nproccol-1)) then
!         do k = 1, innz
!            do i = 1, innx
!              this%FieldQ(i,inny+1+jadd,k+kadd)=rho(i,1,k)*fourpi
!            enddo
!          enddo
!        endif
!        if((nprocrow.eq.1).and.(nproccol.eq.1)) then
!          do i = 1, innx
!            this%FieldQ(i,inny+1,innz+1) = rho(i,1,1)*fourpi
!          enddo
!        endif

        !open(8,file="phir",status="unknown")
        !do i = 1, innx
        !    write(8,1000)(i-1)*hx,this%FieldQ(i,1,1),&
        !    this%FieldQ(i,1,innz/2),this%FieldQ(i,1,innz)
        !enddo
        !close(8)
        !open(9,file="phithe",status="unknown")
        !do j = 1, inny+1
        !  write(9,1000)(j-1)*hy,this%FieldQ(2,j,2),&
        !  this%FieldQ(3,j,2),this%FieldQ(4,j,2)
        !enddo
        !close(9)
        !open(10,file="phiz",status="unknown")
        !do k = 1, innz+1
        !  write(10,1000)(k-1)*hz,this%FieldQ(1,1,k),&
        !  this%FieldQ(innx/2,1,k),this%FieldQ(innx,1,k)
        !enddo
        !close(10)

1000    format(4(1x,e13.6))

        call MPI_BARRIER(comm2d,ierr)
        t_field = t_field + elapsedtime_Timer(t0)

        end subroutine update4_FieldQuant

        ! Solving Poisson's equation with open BCs.
        subroutine finiteBC(innx,inny,innz,rho,hx,hy,hz,nxpylc,nypzlc,&
        nxpzlc,nzpylc,myidz,myidy,npz,npy,commrow,commcol,comm2d,&
        pztable,pytable,ypzstable,xpystable,xpzstable,zpystable,inxglb,&
        inyglb,inzglb,gblam)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx,inny,innz,inxglb,inyglb,inzglb
        integer, intent(in) :: nxpylc,nypzlc,nxpzlc,nzpylc
        integer, intent(in) :: myidz,myidy,npz,npy,commrow,commcol,&
                               comm2d
        integer, dimension(0:npz-1), intent(in) :: pztable,ypzstable,&
                                                   xpzstable
        integer, dimension(0:npy-1), intent(in) :: pytable,xpystable,&
                                                   zpystable
        double precision, dimension(innx,inny,innz), intent(inout) :: rho
        double precision, intent(in) :: hx, hy, hz, gblam
        double precision :: scalex,scaley,scalez,pi
        integer :: i,j,k,kk,nzz
        double complex, dimension(innx,inny,innz) :: rhofft
        integer, dimension(0:npz-1) :: zdisp
        integer, dimension(0:npy-1) :: ydisp
        double precision, dimension(innx-1) :: aa,cc
        double precision, dimension(innx-1,inny,innz) :: bb
!        double precision, dimension(innx,inny,innz) :: testrho
!        double precision :: diffrho

        zdisp(0) = 0
        do i = 1, npz-1
          zdisp(i) = zdisp(i-1) + pztable(i-1)
        enddo
        ydisp(0) = 0
        do i = 1, npy-1
          ydisp(i) = ydisp(i-1) + pytable(i-1)
        enddo

        scalex = 1.0
        scaley = 1.0
        scalez = 1.0

        !print*,"before fft3d4 ",myidy,inxglb,inyglb,inzglb,innz,inny
        !testrho = rho
        call fft3d4_FFT(inxglb,inyglb,inzglb,innz,inny,nxpylc,&
        nypzlc,1,scalex,scaley,scalez,rho,ypzstable,pztable,xpystable,&
        pytable,npz,commrow,npy,commcol,comm2d,myidz,myidy,rhofft)
        !print*,"after fft3d4 ",myidy,nxpylc,nypzlc

        aa(1) = 2.0/(0.5*hx*hx)
        cc(1) = 0.0
        do i = 2, innx-1
          aa(i) = 1.0/(hx*hx) + 0.5/((i-1)*hx*hx)
          cc(i) = 1.0/(hx*hx) - 0.5/((i-1)*hx*hx)
        enddo

        nzz = inzglb/2 + 1
        pi = 2*asin(1.0)
        do k = 1, innz
          kk = k + zdisp(myidz)
          do j = 1, inny
            !bb(1,j,k) = -2.0/(0.5*hx*hx) - &
            !            (2*pi*(k+zdisp(myidz)-1)/gblam)**2
            if(kk.le.nzz) then
              bb(1,j,k) = -2.0/(0.5*hx*hx) - &
                        (2*pi*(kk-1)/gblam)**2
            else
              bb(1,j,k) = -2.0/(0.5*hx*hx) - &
                        (2*pi*(kk-1-inzglb)/gblam)**2
            endif
            do i = 2, innx-1
              !bb(i,j,k) = -2.0/(hx*hx)-(j+ydisp(myidy)-1)**2/&
              !((i-1)*(i-1)*hx*hx)-(2*pi*(k+zdisp(myidz)-1)/gblam)**2
              if(kk.le.nzz) then
                bb(i,j,k) = -2.0/(hx*hx)-(j+ydisp(myidy)-1)**2/&
                ((i-1)*(i-1)*hx*hx)-(2*pi*(kk-1)&
                /gblam)**2
              else
                bb(i,j,k) = -2.0/(hx*hx)-(j+ydisp(myidy)-1)**2/&
                ((i-1)*(i-1)*hx*hx)-(2*pi*(kk-1-inzglb)&
                /gblam)**2
              endif
            enddo
          enddo
        enddo

        do k = 1, innz
          do j = 1, inny
            do i = 2, innx-1
              bb(i,j,k) = bb(i,j,k) - aa(i-1)*cc(i)/bb(i-1,j,k)
            enddo
          enddo
        enddo
        !the negative sign due to the "-" in the Poisson's equation.
        do k = 1, innz
          do j = 1, inny
            do i = 1, innx-2
              rhofft(i,j,k) = - rhofft(i,j,k)
            enddo
            rhofft(innx,j,k) = 0.0
            rhofft(innx-1,j,k) = -rhofft(innx-1,j,k)-aa(innx-1)*&
                                 rhofft(innx,j,k)
          enddo
        enddo

        !forward substitute for Gauss Elemination.
        do k = 1, innz
          do j = 1, inny
            rhofft(1,j,k) = rhofft(1,j,k)
            do i = 2, innx-1
              rhofft(i,j,k) = rhofft(i,j,k) &
                               -rhofft(i-1,j,k)*cc(i)/bb(i-1,j,k)
            enddo
          enddo
        enddo

        !backward substitute.
        do k = 1, innz
          do j = 1, inny
            rhofft(innx-1,j,k) = rhofft(innx-1,j,k)/bb(innx-1,j,k)
            do i = innx-2, 1, -1
              rhofft(i,j,k) = (rhofft(i,j,k)-aa(i)*rhofft(i+1,j,k))/&
                               bb(i,j,k)
            enddo
          enddo
        enddo

        !print*,"before invfft3d4 ",myidy,inxglb,inyglb,inzglb,innz,inny
        !print*,"before invfft3d4 ",myidy,sum(real(rhofft)),sum(imag(rhofft)),nxpzlc,nzpylc
        ! inverse FFT:
        scalex = 1.0/float(inxglb)
        scaley = 1.0/float(inyglb)
        scalez = 1.0/float(inzglb)
        call invfft3d4_FFT(inxglb,inyglb,inzglb,innz,inny,nxpzlc,&
        nzpylc,-1,scalex,scaley,scalez,rhofft,xpzstable,pztable,&
        zpystable,pytable,npz,commrow,npy,commcol,comm2d,myidz,myidy,&
        rho)
        !print*,"after invfft3d4 ",myidy,nxpylc,nypzlc
        
        !diffrho = 0.0
        !do k = 1, innz
        !  do j = 1, inny
        !    do i = 1, innx
        !      diffrho = diffrho + abs(testrho(i,j,k)-rho(i,j,k))
        !    enddo
        !  enddo
        !enddo
        !print*,"diffrho: ",diffrho

        end subroutine finiteBC

!-------------------------------------------------------------------------
        subroutine update5_FieldQuant(this,source,fldgeom,grid,nxlc,&
          nylc,nzlc,nprocrow,nproccol,nylcr,nzlcr)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(in) :: fldgeom
        type (FieldQuant), intent(inout) :: this
        type (Pgrid2d), intent(in) :: grid
        double precision, dimension(nxlc,nylc,nzlc), intent(in) :: source
        integer, intent(in) :: nxlc,nylc,nzlc,&
                               nprocrow,nproccol,nylcr,nzlcr
        double precision, dimension(3) :: msize
        double precision :: hx, hy, hz, temp
        integer, dimension(0:nprocrow-1) :: pztable,xpzstable,zdisp
        integer, dimension(0:nproccol-1) :: xpystable,pytable,ydisp
        double precision , dimension(nxlc-1,nylcr,nzlcr) :: rho
        real*8, dimension(nxlc-1) :: tmprho
        integer :: myid,myidz,myidy,comm2d,commcol,commrow
        integer :: nxpylc,nxpzlc,nsxy1,nsxy2,nsxz1,nsxz2
        integer :: i,j,k,inxglb,inyglb,inzglb,innx,inny,innz
        integer, dimension(2,0:nprocrow-1,0:nproccol-1)::LocalTable
        integer, dimension(3) :: glmshnm
        integer :: jadd,kadd,ierr,ksign,kk,jj
        double precision :: t0,pi,fourpi,testsum,totsum,testrho,totrho,&
                xrad,yrad,scaley
        !double precision , dimension(nxlc-1,nylcr,nzlcr) :: testt
        !double precision, dimension(nxlc) :: betan
        double precision :: length

        call getpost_Pgrid2d(grid,myid,myidy,myidz)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        call getlctabnm_CompDom(fldgeom,LocalTable)
        
        do i = 0, nprocrow-1
          pztable(i) = LocalTable(1,i,0)
        enddo
        zdisp(0) = 0
        do i = 0, nprocrow-1
          zdisp(i) = zdisp(i-1)+pztable(i-1)
        enddo
        do i = 0, nproccol-1
          pytable(i) = LocalTable(2,0,i)
        enddo
        pytable(nproccol-1) = pytable(nproccol-1) - 1 !periodic BC.
        ydisp(0) = 0
        do i = 1, nproccol-1
          ydisp(i) = ydisp(i-1) + pytable(i-1)
        enddo

        call getmnum_CompDom(fldgeom,glmshnm)
        inxglb = glmshnm(1)-1  !periodic BC
        inyglb = glmshnm(2)-1  !periodic BC.
        inzglb = glmshnm(3)

        innz = nzlcr
        inny = nylcr
        innx = nxlc-1  !periodic BC along x.

        !used for transpose between serial x and parallel y.
        nsxy1 = inxglb/nproccol
        nsxy2 = inxglb - nproccol*nsxy1
        do i = 0, nproccol-1
          if(i.le.(nsxy2-1)) then
            xpystable(i) = nsxy1+1
          else
            xpystable(i) = nsxy1
          endif
        enddo
        nxpylc = xpystable(myidy)

        !used for transpose between serial x and parallel z.
        nsxz1 = inxglb/nprocrow  !periodic BC
        nsxz2 = inxglb - nprocrow*nsxz1
        do i = 0, nprocrow-1
          if(i.le.(nsxz2-1)) then
            xpzstable(i) = nsxz1+1
          else
            xpzstable(i) = nsxz1
          endif
        enddo
        nxpzlc = xpzstable(myidz)

        call getmsize_CompDom(fldgeom,msize)
        hx = msize(1)
        hy = msize(2)
        hz = msize(3)
        xrad = inxglb*hx
        yrad = inyglb*hy

        if(nprocrow.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif
        if(nproccol.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif

        do k = 1, innz
          do j = 1, inny
            rho(1,j,k) = source(1,j+jadd,k+kadd)/(hx*hy*hz)
            do i = 2, innx-1
              rho(i,j,k) = source(i,j+jadd,k+kadd)/(hx*hy*hz)
            enddo
            rho(innx,j,k) = source(innx,j+jadd,k+kadd)/&
            (hx*hy*hz)
          enddo
        enddo
        pi = 2*asin(1.0)
        !do k = 1, innz
        !  kk = k + zdisp(myidz)
        !  do j = 1, inny
        !    jj = j + ydisp(myidy)
        !    do i = 1, innx
        !      if((i.gt.innx/4).and.(i.lt.3*innx/4).and. &
        !         (jj.gt.inyglb/4).and.(jj.lt.3*inyglb/4) ) then
        !         !rho(i,j,k) = (i-1)*hx*(cos((j-1)*hy))**2*24.0/&
        !         !             (pi*rad**3)
        !         !rho(i,j,k) = 10.0
        !         rho(i,j,k) = exp(-10*((i-1)*hx+(kk-1)*hz)**2)
        !!      else
        !         rho(i,j,k) = 0.0
        !      endif
        !    enddo
        !  enddo
        !enddo
!
        !testt = rho

        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              tmprho(i) = rho(i,j,k)
            enddo
            call sinft(tmprho,innx)
            do i = 1, innx
              rho(i,j,k) = tmprho(i)
            enddo
          enddo
        enddo
          
        scaley = 1.0 
        call fft3d5_FFT(inxglb,inyglb,innz,inny,nxpylc,&
        scaley,rho,xpystable,pytable,nproccol,commcol,comm2d)

        call Gaussz5(inxglb,inzglb,innz,nxpzlc,inny,rho,&
        xpzstable,pztable,nprocrow,nproccol,commrow,myidz,myidy,ydisp,&
        hz,xrad,yrad)

        scaley = 2.0/float(inyglb) 
        call fft3d5_FFT(inxglb,inyglb,innz,inny,nxpylc,&
        scaley,rho,xpystable,pytable,nproccol,commcol,comm2d)

        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              tmprho(i) = rho(i,j,k)
            enddo
            call sinft(tmprho,innx)
            do i = 1, innx
              rho(i,j,k) = tmprho(i)*2/float(innx)
            enddo
          enddo
        enddo

        !testrho = 0.0
        !do k = 1, innz
        !  do j = 1, inny
        !    do i = 1, innx
        !      !testrho = testrho + abs(rho(i,j,k)-testt(i,j,k))
        !      testrho = testrho + rho(i,j,k)
        !    enddo
        !  enddo
        !enddo
        !call MPI_ALLREDUCE(testrho,totrho,1,MPI_REAL,MPI_SUM, &
        !                     MPI_COMM_WORLD,ierr)
        !!print*,"totrho: ",totrho
        !print*,"testrho: ",testrho
       
        fourpi = 4*pi
        ! 1/(4*pi) is due to the multiplication of 4pi in Rob's code.
        do k = 1, innz
          do j = 1, inny
            this%FieldQ(1,j+jadd,k+kadd) = 0.0
            do i = 2, innx
              this%FieldQ(i,j+jadd,k+kadd) = rho(i,j,k)*fourpi
            enddo
            this%FieldQ(innx+1,j+jadd,k+kadd) = 0.0
          enddo
        enddo
        !wrong peiodic condition, works only if nprocrow = 1
        if(myidy.eq.0) then
          do k = 1, innz
            do i = 1, innx+1
              this%FieldQ(i,1+jadd,k+kadd)= 0.0
            enddo
          enddo
        endif
        if(myidy.eq.(nproccol-1)) then
          do k = 1, innz
            do i = 1, innx+1
              this%FieldQ(i,inny+1+jadd,k+kadd)= 0.0
            enddo
          enddo
        endif

        call MPI_BARRIER(comm2d,ierr)
        t_field = t_field + elapsedtime_Timer(t0)

        end subroutine update5_FieldQuant

        subroutine Gaussz5(nx,nz,nsizez,nsizexz,nsizey,x,&
        xstable,xrtable,nprocrow,nproccol,commrow,myidx,myidy,ydisp,&
        hz,xa,ya)
        implicit none
        include 'mpif.h'
        integer,intent(in) :: nx,nz,nsizez,nsizexz,nsizey
        integer,intent(in) :: nprocrow,commrow,nproccol
        integer,intent(in) :: myidx,myidy
        integer,dimension(0:nprocrow-1),intent(in) :: xstable,xrtable
        double precision, dimension(nx,nsizey,nsizez), intent(inout) :: x
        double precision, intent(in) :: hz,xa,ya
        integer, dimension(0:nproccol-1), intent(in) :: ydisp
        integer,dimension(0:nprocrow-1) :: xdisp
        integer :: i,j,k,jm,kk,jj
        double precision :: t0,tmp1,tmp2,pi
        double precision, dimension(nz) :: aa,cc
        double precision, dimension(nz,nsizey,nsizexz) :: bb
        double precision, allocatable, dimension(:,:,:) :: x0

        call starttime_Timer(t0)

        pi = 2.0*asin(1.0)
        xdisp(0) = 0
        do i = 1, nprocrow-1
          xdisp(i) = xdisp(i-1) + xstable(i-1)
        enddo

        allocate(x0(nz,nsizey,nsizexz))
        !transpose z direction into serial direction:
        call trans3d3r_TRANSP(nx,nsizey,nsizez,nsizexz,x,x0,nprocrow,&
                      xstable,xrtable,commrow,myidx,nz)

        do i = 1, nz -1
          aa(i) = 1.0
        enddo
        aa(nz) = 0.0
        cc(1) = 0.0
        do i = 2, nz 
          cc(i) = 1.0
        enddo
        do k = 1, nsizexz
          kk = xdisp(myidx) + k
          do j = 1, nsizey
            jj = ydisp(myidy) + j
            tmp1 = hz*hz*(((jj-1)*pi/ya)**2+((kk-1)*pi/xa)**2)
            tmp2 = exp(-sqrt(tmp1))
            bb(1,j,k) = -2 - tmp1 + tmp2
            do i = 2, nz-1
              bb(i,j,k) = -2 - tmp1
            enddo
            bb(nz,j,k) = -2 - tmp1 + tmp2
          enddo
        enddo
       
        !eliminate the lower off-diagnal terms.
        do k = 1, nsizexz
          do j = 1, nsizey
            do i = 2, nz
              bb(i,j,k) = bb(i,j,k) - aa(i-1)*cc(i)/bb(i-1,j,k)
            enddo
          enddo
        enddo
        !the negative sign due to the "-" in the Poisson's equation. 
        do k = 1, nsizexz
          do j = 1, nsizey
            do i = 1, nz
              x0(i,j,k) = -x0(i,j,k)*hz*hz
            enddo
          enddo
        enddo

        if(myidx.ne.0 .or. myidy.ne.0) then
        !forward substitute for Gauss Elemination.
          do k = 1, nsizexz
          do j = 1, nsizey
            do i = 2, nz
              x0(i,j,k) = x0(i,j,k) - x0(i-1,j,k)*cc(i)/bb(i-1,j,k)
            enddo
          enddo
        enddo

        !backward substitute
        do k = 1, nsizexz
          do j = 1, nsizey
            x0(nz,j,k) = x0(nz,j,k)/bb(nz,j,k)
            do i = nz-1, 1, -1
              x0(i,j,k) = (x0(i,j,k) - x0(i+1,j,k)*aa(i))/bb(i,j,k)
            enddo
          enddo
        enddo

        else

        !forward substitute for Gauss Elemination.
        do k = 2, nsizexz
          do j = 2, nsizey
            do i = 2, nz
              x0(i,j,k) = x0(i,j,k) - x0(i-1,j,k)*cc(i)/bb(i-1,j,k)
            enddo
          enddo
        enddo

        !backward substitute
        do k = 2, nsizexz
          do j = 2, nsizey
            x0(nz,j,k) = x0(nz,j,k)/bb(nz,j,k)
            do i = nz-1, 1, -1
              x0(i,j,k) = (x0(i,j,k) - x0(i+1,j,k)*aa(i))/bb(i,j,k)
            enddo
          enddo
        enddo

        endif

        call trans3d3r_TRANSP(nz,nsizey,nsizexz,nsizez,x0,x,nprocrow,&
                      xrtable,xstable,commrow,myidx,nx)

        deallocate(x0)

        t_fft2dhpf = t_fft2dhpf + elapsedtime_Timer(t0)

        return

        end subroutine Gaussz5
!-------------------------------------------------------------------------

        subroutine update6_FieldQuant(this,source,fldgeom,grid,nxlc,&
          nylc,nzlc,nprocrow,nproccol,nylcr,nzlcr)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(in) :: fldgeom
        type (FieldQuant), intent(inout) :: this
        type (Pgrid2d), intent(in) :: grid
        double precision, dimension(nxlc,nylc,nzlc), intent(in) :: source
        integer, intent(in) :: nxlc,nylc,nzlc,&
                               nprocrow,nproccol,nylcr,nzlcr
        double precision, dimension(3) :: msize
        double precision :: hx, hy, hz, temp
        integer, dimension(0:nprocrow-1) :: pztable,xpzstable,zdisp
        integer, dimension(0:nproccol-1) :: xpystable,pytable,ydisp
        double precision , dimension(nxlc-1,nylcr,nzlcr) :: rho
        real*8, dimension(nxlc-1) :: tmprho
        integer :: myid,myidz,myidy,comm2d,commcol,commrow
        integer :: nxpylc,nxpzlc,nsxy1,nsxy2,nsxz1,nsxz2
        integer :: i,j,k,inxglb,inyglb,inzglb,innx,inny,innz
        integer, dimension(2,0:nprocrow-1,0:nproccol-1)::LocalTable
        integer, dimension(3) :: glmshnm
        integer :: jadd,kadd,ierr,ksign,jj,kk
        double precision :: t0,pi,fourpi,testsum,totsum,testrho,totrho,&
                xrad,yrad,scaley,zrad
        !double precision , dimension(nxlc-1,nylcr,nzlcr) :: testt
        !double precision, dimension(nxlc) :: betan
        double precision :: length

        call getpost_Pgrid2d(grid,myid,myidy,myidz)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        call getlctabnm_CompDom(fldgeom,LocalTable)
        
        do i = 0, nprocrow-1
          pztable(i) = LocalTable(1,i,0)
        enddo
        pztable(nprocrow-1) = pztable(nprocrow-1) - 1 !periodic BC.
        zdisp(0) = 0
        do i = 0, nprocrow-1
          zdisp(i) = zdisp(i-1)+pztable(i-1)
        enddo
        do i = 0, nproccol-1
          pytable(i) = LocalTable(2,0,i)
        enddo
        pytable(nproccol-1) = pytable(nproccol-1) - 1 !periodic BC.
        ydisp(0) = 0
        do i = 1, nproccol-1
          ydisp(i) = ydisp(i-1) + pytable(i-1)
        enddo

        call getmnum_CompDom(fldgeom,glmshnm)
        inxglb = glmshnm(1)-1  !periodic BC
        inyglb = glmshnm(2)-1  !periodic BC.
        inzglb = glmshnm(3)-1  !periodic BC.

        innz = nzlcr
        inny = nylcr
        innx = nxlc-1  !periodic BC along x.

        !used for transpose between serial x and parallel y.
        nsxy1 = inxglb/nproccol
        nsxy2 = inxglb - nproccol*nsxy1
        do i = 0, nproccol-1
          if(i.le.(nsxy2-1)) then
            xpystable(i) = nsxy1+1
          else
            xpystable(i) = nsxy1
          endif
        enddo
        nxpylc = xpystable(myidy)

        !used for transpose between serial x and parallel z.
        nsxz1 = inxglb/nprocrow  !periodic BC
        nsxz2 = inxglb - nprocrow*nsxz1
        do i = 0, nprocrow-1
          if(i.le.(nsxz2-1)) then
            xpzstable(i) = nsxz1+1
          else
            xpzstable(i) = nsxz1
          endif
        enddo
        nxpzlc = xpzstable(myidz)

        call getmsize_CompDom(fldgeom,msize)
        hx = msize(1)
        hy = msize(2)
        hz = msize(3)
        xrad = inxglb*hx
        yrad = inyglb*hy
        zrad = inzglb*hz

        if(nprocrow.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif
        if(nproccol.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif

        do k = 1, innz
          do j = 1, inny
            rho(1,j,k) = source(1,j+jadd,k+kadd)/(hx*hy*hz)
            do i = 2, innx-1
              rho(i,j,k) = source(i,j+jadd,k+kadd)/(hx*hy*hz)
            enddo
            rho(innx,j,k) = source(innx,j+jadd,k+kadd)/&
            (hx*hy*hz)
          enddo
        enddo
        pi = 2*asin(1.0d0)
        !do k = 1, innz
        !  kk = k + zdisp(myidz)
        !  do j = 1, inny
        !    jj = j + ydisp(myidy)
        !    do i = 1, innx
        !      if((i.gt.innx/4).and.(i.lt.3*innx/4).and. &
        !         (jj.gt.inyglb/4).and.(jj.lt.3*inyglb/4) ) then
        !         !rho(i,j,k) = (i-1)*hx*(cos((j-1)*hy))**2*24.0/&
        !         !             (pi*rad**3)
        !         rho(i,j,k) = 10.0
        !         !rho(i,j,k) = exp(-10*((i-1)*hx+(k-1)*hz)**2)
        !      else
        !         rho(i,j,k) = 0.0
        !      endif
        !    enddo
        !  enddo
        !enddo
        !testrho = 0.0
        !do k = 1, innz
        !  do j = 1, inny
        !    do i = 1, innx
        !      testrho = testrho + rho(i,j,k)*hx*hy*hz
        !    enddo
        !  enddo
        !enddo
        !print*,"testrho-first: ",testrho

        !testt = rho

        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              tmprho(i) = rho(i,j,k)
            enddo
            call sinft(tmprho,innx)
            do i = 1, innx
              rho(i,j,k) = tmprho(i)
            enddo
          enddo
        enddo
          
        scaley = 1.0 
        call fft3d6_FFT(inxglb,inyglb,innz,inny,nxpylc,&
        scaley,rho,xpystable,pytable,nproccol,commcol,comm2d)

        call Gaussz6(inxglb,inzglb,innz,nxpzlc,inny,rho,&
        xpzstable,pztable,nprocrow,nproccol,commrow,myidz,myidy,ydisp,&
        xrad,yrad,zrad)

        scaley = 2.0/float(inyglb) 
        call fft3d6_FFT(inxglb,inyglb,innz,inny,nxpylc,&
        scaley,rho,xpystable,pytable,nproccol,commcol,comm2d)

        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              tmprho(i) = rho(i,j,k)
            enddo
            call sinft(tmprho,innx)
            do i = 1, innx
              rho(i,j,k) = tmprho(i)*2/float(innx)
            enddo
          enddo
        enddo

        !testrho = 0.0
        !do k = 1, innz
        !  do j = 1, inny
        !    do i = 1, innx
        !      !testrho = testrho + abs(rho(i,j,k)-testt(i,j,k))
        !      testrho = testrho + rho(i,j,k)
        !    enddo
        !  enddo
        !enddo
        !print*,"testrho: ",testrho
        !call MPI_ALLREDUCE(testrho,totrho,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
        !                     MPI_COMM_WORLD,ierr)
        !print*,"totrho: ",totrho,testrho
       
        fourpi = 4*pi
        ! 1/(4*pi) is due to the multiplication of 4pi in Rob's code.
        do k = 1, innz
          do j = 1, inny
            this%FieldQ(1,j+jadd,k+kadd) = 0.0
            do i = 2, innx
              this%FieldQ(i,j+jadd,k+kadd) = rho(i,j,k)*fourpi
            enddo
            this%FieldQ(innx+1,j+jadd,k+kadd) = 0.0
          enddo
        enddo
        !wrong peiodic condition, works only if nprocrow = 1
        if(myidy.eq.0) then
          do k = 1, innz
            do i = 1, innx+1
              this%FieldQ(i,1+jadd,k+kadd)= 0.0
            enddo
          enddo
        endif
        if(myidy.eq.(nproccol-1)) then
          do k = 1, innz
            do i = 1, innx+1
              this%FieldQ(i,inny+1+jadd,k+kadd)= 0.0
            enddo
          enddo
        endif
        if(myidz.eq.(nprocrow-1)) then
          if(myidy.ne.(nproccol-1)) then
            do j = 1, inny
              do i = 1, innx+1
                this%FieldQ(i,j+jadd,innz+1+kadd)=rho(i,j,1)*fourpi
              enddo
            enddo
          else
            do j = 1, inny+1
              do i = 1, innx+1
                this%FieldQ(i,j+jadd,innz+1+kadd)=rho(i,j,1)*fourpi
              enddo
            enddo
          endif
        endif
        
        call MPI_BARRIER(comm2d,ierr)
        t_field = t_field + elapsedtime_Timer(t0)

        end subroutine update6_FieldQuant

        subroutine Gaussz6(nx,nz,nsizez,nsizexz,nsizey,x,&
        xstable,xrtable,nprocrow,nproccol,commrow,myidx,myidy,ydisp,&
        xa,ya,za)
        implicit none
        include 'mpif.h'
        integer,intent(in) :: nx,nz,nsizez,nsizexz,nsizey
        integer,intent(in) :: nprocrow,commrow,nproccol
        integer,intent(in) :: myidx,myidy
        integer,dimension(0:nprocrow-1),intent(in) :: xstable,xrtable
        double precision, dimension(nx,nsizey,nsizez), intent(inout) :: x
        double precision, intent(in) :: xa,ya,za
        integer, dimension(0:nproccol-1), intent(in) :: ydisp
        integer,dimension(0:nprocrow-1) :: xdisp
        integer :: i,j,k,jm,kk,jj,ii,ksign
        double precision :: t0,tmp2,pi,scalez
        double precision, dimension(nz,nsizey) :: tmp1,tmp10
        double precision, allocatable, dimension(:,:,:) :: x0
        integer :: jadd,kadd

        call starttime_Timer(t0)

        pi = 2.0*asin(1.0)
        xdisp(0) = 0
        do i = 1, nprocrow-1
          xdisp(i) = xdisp(i-1) + xstable(i-1)
        enddo

        allocate(x0(nz,nsizey,nsizexz))
        !transpose z direction into serial direction:
        call trans3d3r_TRANSP(nx,nsizey,nsizez,nsizexz,x,x0,nprocrow,&
                      xstable,xrtable,commrow,myidx,nz)

        if(myidx.ne.0) then
          kadd = 0
        else
          kadd = 1
        endif
        if(myidy.ne.0) then
          jadd = 0
        else
          jadd = 1
        endif

        ksign = 1
        scalez = 1.0
        pi = 2*asin(1.0)
        do k = 1+kadd, nsizexz
          kk = xdisp(myidx) + k - 1
          do j = 1, nsizey
            do i = 1, nz
              tmp1(i,j) = x0(i,j,k)
            enddo
          enddo

          call fftrclocal_FFT(ksign,scalez,tmp1,nz,nsizey,tmp10)

          do j = 1+jadd, nsizey
            jj = ydisp(myidy) + j - 1

            x0(1,j,k) = tmp10(1,j)/((kk*pi/xa)**2+(jj*pi/ya)**2)
            x0(2,j,k) = tmp10(2,j)/((kk*pi/xa)**2+(jj*pi/ya)**2+&
                          (2*(nz/2)*pi/za)**2)
            !x0(1,j,k) = tmp10(1,j)
            !x0(2,j,k) = tmp10(2,j)
            do i = 3, nz
              ii = (i+1)/2 - 1
              !no negative sign is due to -rho in Poisson's equation.
              x0(i,j,k) = tmp10(i,j)/((kk*pi/xa)**2+(jj*pi/ya)**2+&
                          (2*ii*pi/za)**2)
              !x0(i,j,k) = tmp10(i,j)
            enddo
          enddo
        enddo

        ksign = -1
        scalez = 1.0/float(nz)
        !scalez = 2.0/float(nz)
        do k = 1, nsizexz
          do j = 1, nsizey
            do i = 1, nz
              tmp1(i,j) = x0(i,j,k)
            enddo
          enddo

          call fftcrlocal_FFT(ksign,scalez,tmp1,nz,nsizey,tmp10)

          do j = 1, nsizey
            do i = 1, nz
              x0(i,j,k) = tmp10(i,j) !using machine lib mfftcrlocal
              !x0(i,j,k) = tmp10(i,j)*scalez*2 !using NR lib
            enddo
          enddo
        enddo

        call trans3d3r_TRANSP(nz,nsizey,nsizexz,nsizez,x0,x,nprocrow,&
                      xrtable,xstable,commrow,myidx,nx)

        deallocate(x0)

        t_fft2dhpf = t_fft2dhpf + elapsedtime_Timer(t0)

        return
        end subroutine Gaussz6

!----------------------------------------------------------------------
! update potential (solving 2D Possion's equation) with 2D open 
! boundary conditions.
        subroutine update7_FieldQuant(this,source,fldgeom,grid,nxlc,&
          nylc,nzlc,nprocrow,nproccol,nylcr,nzlcr)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(in) :: fldgeom
        double precision, dimension(nxlc,nylc,nzlc), intent(in) :: source
        type (Pgrid2d), intent(in) :: grid
        integer, intent(in) :: nxlc,nylc,nzlc,&
                               nprocrow,nproccol,nylcr,nzlcr
        type (FieldQuant), intent(inout) :: this
        double precision, dimension(3) :: msize
        double precision :: hx, hy, hz, temp, rhoavgz
        integer, dimension(0:nprocrow-1) :: ypzstable,pztable
        integer, dimension(0:nproccol-1) :: xpystable,pytable
        double precision , dimension(nxlc,nylcr,nzlcr) :: rho
        integer :: myid,myidz,myidy,&
                   comm2d,commcol,commrow
        integer :: nxpylc2,nypzlc2
        integer :: nsxy1,nsxy2,nsyz1,nsyz2
        integer :: i,j,k,inxglb,inyglb,inzglb,innx,inny,innz
        integer, dimension(2,0:nprocrow-1,0:nproccol-1)::LocalTable
        integer, dimension(3) :: glmshnm
        integer :: jadd,kadd,ierr
        double precision :: t0

        call getpost_Pgrid2d(grid,myid,myidy,myidz)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        call getlctabnm_CompDom(fldgeom,LocalTable)
        
        do i = 0, nprocrow-1
          pztable(i) = LocalTable(1,i,0)
        enddo
        do i = 0, nproccol-1
          pytable(i) = LocalTable(2,0,i)
        enddo

        call getmnum_CompDom(fldgeom,glmshnm)
        inxglb = glmshnm(1) 
        inyglb = glmshnm(2)
        inzglb = glmshnm(3)

        innz = nzlcr
        inny = nylcr
        innx = nxlc

        if(nprocrow.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif
        if(nproccol.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif
        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              rho(i,j,k) = source(i,j+jadd,k+kadd)
            enddo
          enddo
        enddo
        
        ! +1 is from the real to complex fft.
        nsxy1 = (inxglb+1)/nproccol
        nsxy2 = (inxglb+1) - nproccol*nsxy1
        do i = 0, nproccol-1
          if(i.le.(nsxy2-1)) then
            xpystable(i) = nsxy1+1
          else
            xpystable(i) = nsxy1
          endif
        enddo

        nsyz1 = 2*inyglb/nprocrow
        nsyz2 = 2*inyglb - nprocrow*nsyz1
        do i = 0, nprocrow-1
          if(i.le.(nsyz2-1)) then
            ypzstable(i) = nsyz1+1
          else
            ypzstable(i) = nsyz1
          endif
        enddo

        nxpylc2 = xpystable(myidy)
        nypzlc2 = ypzstable(myidz)

        call getmsize_CompDom(fldgeom,msize)
        hx = msize(1)
        hy = msize(2)
        hz = msize(3)

!        print*, 'Just before grid flatten:'
!        do k = 1, innz
!           write(*,*) 8,8,k,rho(8,8,k)
!        enddo

        ! Flatten the 3D grid to produce 2D grid values:
        do i = 1, innx
          do j = 1, inny
            rhoavgz = sum(rho(i,j,:))/dble(innz)
            do k = 1, innz
              rho(i,j,k) = rhoavgz  !Avg over z to produce 2D density
            enddo
          enddo
        enddo        

!        print*, 'Just after grid flatten:'
!        print*, 'innx,inny,innz = ',innx,inny,innz
!        do k = 1, innz
!           write(*,*) 8,8,k,rho(8,8,k)
!        enddo

        ! Open 2D boundary conditions!
!        call openBC3D(innx,inny,innz,rho,hx,hy,hz, &
        call openBC2D(innx,inny,innz,rho,hx,hy,hz, &
        nxpylc2,nypzlc2,myidz,myidy,nprocrow,nproccol,commrow,commcol,&
        comm2d,pztable,pytable,ypzstable,xpystable,&
        inxglb,inyglb,inzglb)

!        print*, 'Just after openBC2D:'

        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              this%FieldQ(i,j+jadd,k+kadd) = rho(i,j,k)
            enddo
          enddo
        enddo

        call MPI_BARRIER(comm2d,ierr)
        t_field = t_field + elapsedtime_Timer(t0)

        end subroutine update7_FieldQuant


        ! Solving Poisson's equation with open BCs.
        subroutine openBC2D(innx,inny,innz,rho,hx,hy,hz,&
           nxpylc2,nypzlc2,myidz,myidy,npz,npy,commrow,commcol,comm2d, &
           pztable,pytable,ypzstable,xpystable,inxglb,&
           inyglb,inzglb)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx,inny,innz,inxglb,inyglb,inzglb
        integer, intent(in) :: nxpylc2,nypzlc2
        integer, intent(in) :: myidz,myidy,npz,npy,commrow,commcol,&
                               comm2d
        integer, dimension(0:npz-1), intent(in) :: pztable,ypzstable
        integer, dimension(0:npy-1), intent(in) :: pytable,xpystable
        double precision, dimension(innx,inny,innz), intent(inout) :: rho
        double precision, intent(in) :: hx, hy, hz
        double precision :: scalex,scaley,scalez
        integer :: i,j,k,n1,n2,n3,nylc22,nzlc22
        double complex, allocatable, dimension(:,:,:) :: rho2out
        double complex, allocatable, dimension(:,:,:) :: grn
        integer :: ginny,ginnz,ierr

!        print*, 'Inside openBC2D.'

        n1 = 2*inxglb
        n2 = 2*inyglb
        n3 = 2*inzglb

        nylc22 = nxpylc2
        nzlc22 = nypzlc2

        scalex = 1.0
        scaley = 1.0
        scalez = 1.0

        allocate(rho2out(n3,nylc22,nzlc22))

        call fft3d1_FFT(n1,n2,n3,innz,inny,nylc22,nzlc22,&
                1,scalex,scaley,scalez,rho,ypzstable,pztable,&
        xpystable,pytable,npz,commrow,npy,commcol,comm2d,myidz,myidy, &
        rho2out)

        !c compute FFT of the Green function on the grid:
        ! here the +1 is from the unsymmetry of green function
        ! on double-sized grid.
        if(myidz.eq.(npz-1)) then
           ginnz = innz + 1
        else
           ginnz = innz
        endif
        if(myidy.eq.(npy-1)) then
           ginny = inny + 1
        else
           ginny = inny
        endif
        allocate(grn(n3,nylc22,nzlc22))
!        print*, 'Just before greenf2D.'
        call greenf2D(inxglb,inyglb,inzglb,ginnz,ginny,nylc22,nzlc22, &
               hx,hy,hz,myidz,npz,commrow,myidy,npy,commcol,comm2d,&
                  ypzstable,pztable,xpystable,pytable,grn)
!        print*, 'Just after greenf2D.'

        ! multiply transformed charge density and transformed Green 
        ! function:
!        print*, 'After forward FFTs:'
!        print*, 'nzlc22,nylc22,n3=',nzlc22,nylc22,n3
!        do i = 1, n3
!           write(*,10) i,8,8,real(rho2out(i,8,8)),aimag(rho2out(i,8,8)),&
!           real(grn(i,8,8)),aimag(grn(i,8,8))
!        enddo

        do k = 1, nzlc22
          do j = 1, nylc22
            do i = 1, n3
!              write(*,10) i,8,8,real(rho2out(i,8,8)),aimag(rho2out(i,8,8)),&
!                      real(grn(i,8,8)),aimag(grn(i,8,8))
              rho2out(i,j,k) = rho2out(i,j,k)*grn(i,j,k)
            enddo
          enddo
        enddo
10      format(3(1x,I2),4(1x,g32.16))

        deallocate(grn)

        ! inverse FFT:
        scalex = 1.0/float(n1)
        scaley = 1.0/float(n2)
        scalez = 1.0/float(n3)
!        print*, 'Just before inverse FFT'
        call invfft3d1_FFT(n3,n2,n1,nylc22,nzlc22,inny,innz,&
               -1,scalex,scaley,scalez,rho2out,pztable,ypzstable,&
        pytable,xpystable,npz,commrow,npy,commcol,comm2d,myidz,myidy,&
        rho)
!        print*, 'Just after inverse FFT'

        deallocate(rho2out)

        end subroutine openBC2D

        ! green function for extended array.
        subroutine greenf2D(nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz,&
                  hx,hy,hz,myidx,npx,commrow,myidy,npy,commcol,comm2d,&
                   xstable,xrtable,ystable,yrtable,grnout)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz
        integer, intent(in) :: myidx,myidy,npx,npy,commrow,commcol,&
                               comm2d
        integer, dimension(0:npx-1),intent(in) :: xstable,xrtable
        integer, dimension(0:npy-1),intent(in) :: ystable,yrtable
        double precision, intent(in) :: hx, hy, hz
        double complex, intent (out), &
                       dimension (2*nz,nsizexy,nsizeyz) :: grnout
        integer, dimension(0:npx-1) :: gxrtable
        integer, dimension(0:npy-1) :: gyrtable
        integer :: i,j,k,kk,ii,jj,iii,jjj,kkk,n1,n2,n3,ns1,ns2
        integer :: nblockx,nblocky,ksign,nxx,ierr
        double precision, dimension (nx+1,nsizey,nsizez) :: grn
        double precision :: scalex,scaley,scalez
        double precision :: t0,fudge
        double precision, dimension(2*nx,nsizey) :: tmp1
        double complex, dimension(nx+1,nsizey) :: tmp10
        double complex, dimension(2*ny,nsizexy) :: tmp2
        double complex, dimension(2*nz,nsizexy) :: tmp3
        double complex, allocatable, dimension(:,:,:) :: x1
        double complex, allocatable, dimension(:,:,:) :: x0

        call starttime_Timer(t0)

        gxrtable = xrtable
        gxrtable(npx-1) = xrtable(npx-1) + 1
        nblockx = 0
        do i = 0, myidx-1
          nblockx = nblockx + gxrtable(i)
        enddo

        gyrtable = yrtable
        gyrtable(npy-1) = yrtable(npy-1) + 1
        nblocky = 0
        do i = 0, myidy-1
          nblocky = nblocky + gyrtable(i)
        enddo

!  Best so far:  fudge = 1.31d0
!        fudge = 1.38d0 (big)
        fudge = 1.37d0

        do k = 1, nsizez
          do j = 1, nsizey
            do i = 1, nx+1
              jj = j + nblocky
              kk = k + nblockx
                kkk = kk - 1
                jjj = jj - 1
                iii = i - 1
              if((iii*iii+jjj*jjj).ne.0) then
               grn(i,j,k)=-fudge*log((hx*iii)**2+(hy*jjj)**2)   !2D case
              endif
            enddo
          enddo
        enddo
        if((myidx.eq.0).and.(myidy.eq.0)) then
          if(nsizez.gt.1) then
            grn(1,1,1) = grn(1,1,2)
          else
            grn(1,1,1) = 1.0
          endif
        endif

        scalex = 1.0
        scaley = 1.0
        scalez = 1.0
        n1 = 2*nx
        n2 = 2*ny
        n3 = 2*nz
        ksign = 1

        nxx = n1/2 + 1
        !FFTs along y and z dimensions.
        allocate(x0(n1/2+1,nsizey,nsizez))
        do k = 1, nsizez
          do j = 1, nsizey 
            do i = 1, n1/2+1
              tmp1(i,j) = grn(i,j,k)
            enddo
            do i = n1/2+2, n1
              tmp1(i,j) = grn(n1-i+2,j,k)
            enddo
          enddo

          ! FFTs along z dimensions:
          call fftrclocal_FFT(ksign,scalex,tmp1,n1,nsizey,tmp10)

          do j = 1, nsizey 
            do i = 1, n1/2+1
              x0(i,j,k) = tmp10(i,j)
            enddo
          enddo
        enddo

        allocate(x1(n2/2+1,nsizexy,nsizez))

        ! FFTs along y dimensions:
!        call MPI_BARRIER(commcol,ierr)
        call trans3d_TRANSP(nxx,n2/2+1,nsizexy,nsizey,x0,x1,npy,&
                     ystable,gyrtable,commcol,nsizez)
        deallocate(x0)
        allocate(x0(n2,nsizexy,nsizez))

        do k = 1, nsizez
          do j = 1, nsizexy 
            do i = 1, n2/2+1
              tmp2(i,j) = x1(i,j,k)
            enddo
            do i = n2/2+2,n2
              tmp2(i,j) = x1(n2-i+2,j,k)
            enddo
          enddo

          call fftlocal_FFT(ksign,scaley,tmp2,n2,nsizexy)

          do j = 1, nsizexy 
            do i = 1, n2
              x0(i,j,k) = tmp2(i,j) 
            enddo
          enddo
        enddo

        deallocate(x1)
        allocate(x1(n3/2+1,nsizexy,nsizeyz))
!        call MPI_BARRIER(commcol,ierr)
        call trans3d3_TRANSP(n2,nsizexy,nsizez,nsizeyz,x0,x1,npx,&
                      xstable,gxrtable,commrow,myidx,n3/2+1)
        deallocate(x0)

        do k = 1, nsizeyz
          do j = 1, nsizexy
            do i = 1, n3/2+1
              tmp3(i,j) = x1(i,j,k) 
            enddo
            do i = n3/2+2,n3
              tmp3(i,j) = x1(n3-i+2,j,k)
            enddo
          enddo

          call fftlocal_FFT(ksign,scalez,tmp3,n3,nsizexy)

          do j = 1, nsizexy
            do i = 1, n3
              grnout(i,j,k) = tmp3(i,j)
            enddo
          enddo
        enddo

        deallocate(x1)

        t_greenf = t_greenf + elapsedtime_Timer(t0)

        end subroutine greenf2D

        ! green function for extended array.
        subroutine greenfInt2D(nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz,&
                  hx,hy,hz,myidx,npx,commrow,myidy,npy,commcol,comm2d,&
                   xstable,xrtable,ystable,yrtable,grnout)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz
        integer, intent(in) :: myidx,myidy,npx,npy,commrow,commcol,&
                               comm2d
        integer, dimension(0:npx-1),intent(in) :: xstable,xrtable
        integer, dimension(0:npy-1),intent(in) :: ystable,yrtable
        double precision, intent(in) :: hx, hy, hz
        double complex, intent (out), &
                       dimension (2*nz,nsizexy,nsizeyz) :: grnout
        integer, dimension(0:npx-1) :: gxrtable
        integer, dimension(0:npy-1) :: gyrtable
        integer :: i,j,k,kk,ii,jj,iii,jjj,kkk,n1,n2,n3,ns1,ns2
        integer :: nblockx,nblocky,ksign,nxx,ierr
        double precision, dimension (nx+1,nsizey,nsizez) :: grn
        double precision, dimension (nx+2,nsizey+1,nsizez+1) :: grntmp
        double precision :: scalex,scaley,scalez
        double precision :: t0
        double precision, dimension(2*nx,nsizey) :: tmp1
        double complex, dimension(nx+1,nsizey) :: tmp10
        double complex, dimension(2*ny,nsizexy) :: tmp2
        double complex, dimension(2*nz,nsizexy) :: tmp3
        double complex, allocatable, dimension(:,:,:) :: x1
        double complex, allocatable, dimension(:,:,:) :: x0
        double precision :: rr,aa,bb,cc,dd,ee,ff,ss
        double complex :: gg,gg2
        double complex :: ggrr
        double precision, dimension(2) :: xx,yy,zz
        double precision, dimension(3) :: vv
        integer :: n,i0,j0,k0
        double precision :: recfourpi

        recfourpi = 1.0/(8.0*asin(1.0))

        call starttime_Timer(t0)

!        if(myidx.eq.0 .and. myidy.eq.0) then
!          print*,"into integrated Green function....."
!        endif
        gxrtable = xrtable
        gxrtable(npx-1) = xrtable(npx-1) + 1
        nblockx = 0
        do i = 0, myidx-1
          nblockx = nblockx + gxrtable(i)
        enddo

        gyrtable = yrtable
        gyrtable(npy-1) = yrtable(npy-1) + 1
        nblocky = 0
        do i = 0, myidy-1
          nblocky = nblocky + gyrtable(i)
        enddo

        do k0 = 1, nsizez+1
          do j0 = 1, nsizey+1
            do i0 = 1, nx+2
              jj = j0 + nblocky
              kk = k0 + nblockx
              kkk = kk - 1
              jjj = jj - 1
              iii = i0 - 1

              vv(1) = iii*hx-hx/2
              vv(2) = jjj*hy-hy/2
              vv(3) = kkk*hz-hz/2
            
                rr = sqrt(vv(1)**2+vv(2)**2+vv(3)**2)
                aa = vv(1)**2*vv(3)+(vv(2)**2+vv(3)**2)*vv(3) + &
                            vv(1)*vv(3)*rr
                bb = vv(1)**2*vv(2) + vv(1)*vv(2)*rr
                cc = vv(1)*(vv(1)**2+vv(2)**2)+vv(3)*vv(1)*(vv(3)+rr)
                dd = vv(3)*vv(2)*(vv(3) + rr)
                ee = vv(1)**2*vv(2)+vv(2)*(vv(2)**2+vv(3)*(vv(3)+rr))
                ff = vv(1)*vv(3)*(vv(3)+rr)
                ss = 4*vv(2)*vv(3)*log(vv(1)+rr) + &
                        4*vv(1)*vv(3)*log(vv(2)+rr) + &
                        4*vv(1)*vv(2)*log(vv(3)+rr)

                gg2 = cmplx(0.0,vv(3)**2)*log(cmplx(aa**2-bb**2,2*aa*bb)/ &
                        (aa**2+bb**2) ) + cmplx(0.0,vv(1)**2)*&
                        log(cmplx(cc**2-dd**2,2*cc*dd )/(cc**2+dd**2))+&
                        cmplx(0.0,vv(2)**2)*log(cmplx(ee**2-ff**2,2*ee*ff)/ &
                        (ee**2+ff**2) ) + ss
               
              !grntmp(i0,j0,k0) = recfourpi*real(gg2)/(4*hx*hy*hz) !wrong in Z code
              grntmp(i0,j0,k0) = real(gg2)/(4*hx*hy*hz)

            enddo
          enddo
        enddo

        do k0 = 1, nsizez
          do j0 = 1, nsizey
            do i0 = 1, nx+1
              grn(i0,j0,k0) = grntmp(i0+1,j0+1,k0+1)-grntmp(i0,j0+1,k0+1)-grntmp(i0+1,j0,k0+1)+&
                              grntmp(i0,j0,k0+1)-grntmp(i0+1,j0+1,k0)+grntmp(i0,j0+1,k0)+&
                              grntmp(i0+1,j0,k0)-grntmp(i0,j0,k0)
            enddo
          enddo
        enddo

!              xx(1) = iii*hx-hx/2
!              xx(2) = iii*hx+hx/2
!              yy(1) = jjj*hy-hy/2
!              yy(2) = jjj*hy+hy/2
!              zz(1) = kkk*hz-hz/2
!              zz(2) = kkk*hz+hz/2
!       
!              !find the integrated Green function.
!              n = 0
!              do k = 1, 2
!                do j = 1, 2
!                  do i = 1, 2
!                    n = n+1
!                    rr(n) = sqrt(xx(i)**2+yy(j)**2+zz(k)**2)
!                    aa(n) = xx(i)**2*zz(k)+(yy(j)**2+zz(k)**2)*zz(k) + &
!                            xx(i)*zz(k)*rr(n)
!                    bb(n) = xx(i)**2*yy(j) + xx(i)*yy(j)*rr(n)
!                    cc(n) = xx(i)*(xx(i)**2+yy(j)**2)+zz(k)*xx(i)*(zz(k)+rr(n))
!                    dd(n) = zz(k)*yy(j)*(zz(k) + rr(n))
!                    ee(n) = xx(i)**2*yy(j)+yy(j)*(yy(j)**2+zz(k)*(zz(k)+rr(n)))
!                    ff(n) = xx(i)*zz(k)*(zz(k)+rr(n))
!                    ss(n) = 4*yy(j)*zz(k)*log(xx(i)+rr(n)) + &
!                            4*xx(i)*zz(k)*log(yy(j)+rr(n)) + &
!                            4*xx(i)*yy(j)*log(zz(k)+rr(n))
!                    gg(n) = cmplx(0.0,zz(k)**2)*log(cmplx(aa(n)**2-bb(n)**2,2*aa(n)*bb(n))/ &
!                            (aa(n)**2+bb(n)**2) ) + cmplx(0.0,xx(i)**2)*&
!                            log(cmplx(cc(n)**2-dd(n)**2,2*cc(n)*dd(n) )/(cc(n)**2+dd(n)**2))+&
!                            cmplx(0.0,yy(j)**2)*log(cmplx(ee(n)**2-ff(n)**2,2*ee(n)*ff(n))/ &
!                            (ee(n)**2+ff(n)**2) )
!                    gg2(n) = ss(n) +  gg(n)
!                  enddo
!                enddo
!              enddo
!              ggrr = (-gg2(1)+gg2(2)+gg2(3)-gg2(4)+gg2(5)-gg2(6)-gg2(7)+gg2(8))/4
!              grn(i0,j0,k0) = real(ggrr)/(hx*hy*hz)
!            enddo
!          enddo
!        enddo


        if((myidx.eq.0).and.(myidy.eq.0)) then
          if(nsizez.gt.1) then
            grn(1,1,1) = grn(1,1,2)
          else
            grn(1,1,1) = 1.0
          endif
        endif

        scalex = 1.0
        scaley = 1.0
        scalez = 1.0
        n1 = 2*nx
        n2 = 2*ny
        n3 = 2*nz
        ksign = 1

        nxx = n1/2 + 1
        !FFTs along y and z dimensions.
        allocate(x0(n1/2+1,nsizey,nsizez))
        do k = 1, nsizez
          do j = 1, nsizey 
            do i = 1, n1/2+1
              tmp1(i,j) = grn(i,j,k)
            enddo
            do i = n1/2+2, n1
              tmp1(i,j) = grn(n1-i+2,j,k)
            enddo
          enddo

          ! FFTs along z dimensions:
          call fftrclocal_FFT(ksign,scalex,tmp1,n1,nsizey,tmp10)

          do j = 1, nsizey 
            do i = 1, n1/2+1
              x0(i,j,k) = tmp10(i,j)
            enddo
          enddo
        enddo

        allocate(x1(n2/2+1,nsizexy,nsizez))

        ! FFTs along y dimensions:
!        call MPI_BARRIER(commcol,ierr)
        call trans3d_TRANSP(nxx,n2/2+1,nsizexy,nsizey,x0,x1,npy,&
                     ystable,gyrtable,commcol,nsizez)
        deallocate(x0)
        allocate(x0(n2,nsizexy,nsizez))

        do k = 1, nsizez
          do j = 1, nsizexy 
            do i = 1, n2/2+1
              tmp2(i,j) = x1(i,j,k)
            enddo
            do i = n2/2+2,n2
              tmp2(i,j) = x1(n2-i+2,j,k)
            enddo
          enddo

          call fftlocal_FFT(ksign,scaley,tmp2,n2,nsizexy)

          do j = 1, nsizexy 
            do i = 1, n2
              x0(i,j,k) = tmp2(i,j) 
            enddo
          enddo
        enddo

        deallocate(x1)
        allocate(x1(n3/2+1,nsizexy,nsizeyz))
!        call MPI_BARRIER(commcol,ierr)
        call trans3d3_TRANSP(n2,nsizexy,nsizez,nsizeyz,x0,x1,npx,&
                      xstable,gxrtable,commrow,myidx,n3/2+1)
        deallocate(x0)

        do k = 1, nsizeyz
          do j = 1, nsizexy
            do i = 1, n3/2+1
              tmp3(i,j) = x1(i,j,k) 
            enddo
            do i = n3/2+2,n3
              tmp3(i,j) = x1(n3-i+2,j,k)
            enddo
          enddo

          call fftlocal_FFT(ksign,scalez,tmp3,n3,nsizexy)

          do j = 1, nsizexy
            do i = 1, n3
              grnout(i,j,k) = tmp3(i,j)
            enddo
          enddo
        enddo

        deallocate(x1)

        t_greenf = t_greenf + elapsedtime_Timer(t0)

        end subroutine greenfInt2D


!-------------------------------------------------------------------------
        subroutine setval_FieldQuant(this,i,j,k,value)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(out) :: this
        integer, intent(in) :: i, j, k
        double precision, intent(in) :: value

        this%FieldQ(i,j,k) = value

        end subroutine setval_FieldQuant

        double precision function get_FieldQuant(this,i,j,k)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(in) :: this
        integer, intent(in) :: i, j, k

        get_FieldQuant = this%FieldQ(i,j,k)

        end function get_FieldQuant

        subroutine getglb_FieldQuant(this,temp)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(in) :: this
        type (FieldQuant), intent(out) :: temp
        integer :: i, j, k, lcnz,lcny,lcnx
        double precision :: value

        lcnx = this%Nxlocal
        lcny = this%Nylocal
        lcnz = this%Nzlocal
    
        do k = 1, lcnz
          do j = 1, lcny
            do i = 1, lcnx
              value = get_FieldQuant(this,i,j,k)
              call setval_FieldQuant(temp,i,j,k,value)
            enddo
          enddo 
        enddo

        end subroutine getglb_FieldQuant

        subroutine getlcgrid_FieldQuant(this,nxlc,nylc,nzlc)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(in) :: this
        integer, intent(out) :: nxlc,nylc,nzlc

        nxlc = this%Nxlocal
        nylc = this%Nylocal
        nzlc = this%Nzlocal

        end subroutine

        subroutine destruct_FieldQuant(this)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(out) :: this

        deallocate(this%FieldQ) 

        end subroutine destruct_FieldQuant

!-------------------------------------------------------------------------
       !longitudinal and transverse wakefield using analytical formulae
       subroutine wakefieldana_FieldQuant(Nz,xwakez,ywakez,recvdensz,exwake,eywake,ezwake,&
                                       hz,aa,gg,leng,flagbtw)
       implicit none
       include 'mpif.h'
       integer, intent(in) :: Nz,flagbtw
       double precision, intent(in) :: hz, aa, gg, leng
       double precision, dimension(Nz,2), intent(in) :: recvdensz
       double precision, dimension(Nz), intent(in) :: xwakez,ywakez
       double precision, dimension(Nz), intent(out) :: exwake,ezwake,eywake
       double precision, dimension(2*Nz,1) :: densz2n,densz2nout,&
                         greenwake,greenwakeout
       integer :: kz,twonz,one,ksign,kkzz,i
       double precision :: scale,zz00,pilc,Z0,alpha1,alpha,zz
       double precision :: coef1,coef2,densconst,offset,zzmax
       real*8 :: tmptmp,pbtw1,pbtw2,pbtw3,pbtw4
       real*8 :: leng1,leng2

       pilc = 2*asin(1.0)
  
       do kz = 1, Nz
          densz2n(kz,1) = recvdensz(kz,1) 
       enddo
       do kz = Nz + 1, 2*Nz
          densz2n(kz,1) = 0.0d0
       enddo

       twonz = 2*Nz
       one = 1
       ksign = 1
       scale = 1.0d0
       call fftrclocal2_FFT(ksign,scale,densz2n,twonz,one,densz2nout)

       !longitudinal wakefield function
       !the following longitudinal wakefield function is from
       !"Short-range dipole wakefields in accelerating structures for the NLC",
       !K. L. F. Bane, SLAC-PUB-9663, 2003.
       pilc = 2*asin(1.0)
       alpha1 = 0.4648
       alpha = 1.-alpha1*sqrt(gg/leng)-(1.-2*alpha1)*(gg/leng)
       !slac wake
       !zz00 = gg*(aa/(alpha*leng))**2/8 
       !fermi wake
       zz00 = 0.41*(abs(aa/leng))**0.8*(gg/leng)**1.6*aa
       Z0 = 120*pilc
       !greenwake(1,1) = 0.0
       zz = 0.0d0
       !parameter for Elettra BTW
       pbtw1 = 1226.0d0
       pbtw2 = 3.0d-4
       pbtw3 = 0.494
       pbtw4 = 494.0

       !here, leng1 and leng2 are the length factors.
       !Elegant uses RF module length instead of cavity length.  
       leng1 = 1.0500439034821d0
       leng2 = 1.33408842738323d0
       if(flagbtw.eq.1) then !for Backward TWS
         !The 0.5 is from S+
         greenwake(1,1) = 0.5d12*(pbtw1*exp(-sqrt(zz/pbtw2))+&
            pbtw3*2*(sqrt(zz+0.5*hz))/(hz)+pbtw4*sqrt(zz))
         do kz = 2, Nz+1
           zz = (kz-1)*hz
           greenwake(kz,1) = 1.0d12*(pbtw1*exp(-sqrt(zz/pbtw2))+&
             pbtw3*2*(sqrt(zz+0.5*hz)-sqrt(zz-0.5*hz))/(hz)+&
             pbtw4*sqrt(zz))
         enddo
       else if(flagbtw.eq.2) then !for Tesla standing wave structure
         !The 0.5 is from S+
         greenwake(1,1) = leng1*0.5d12*38.1*(1.165*exp(-sqrt(zz/3.65d-3))-&
            0.165)
         do kz = 2, Nz+1
           zz = (kz-1)*hz
           greenwake(kz,1) = leng1*1.0d12*38.1*(1.165*exp(-sqrt(zz/3.65d-3))-&
             0.165)
         enddo
       else if(flagbtw.eq.3) then !for Tesla 3rd harm. standing wave structure
         !The 0.5 is from S+
         !greenwake(1,1) = leng2*0.5d12*130.0d0*(1.075*exp(-sqrt(zz/2.25d-3))-&
         !   0.075)
         !3/18/08 changed following Sasha's new formulae
         greenwake(1,1) = leng2*0.5d12*130.0d0*(1.02*exp(-sqrt(zz/2.0d-3))-&
            0.02)
         do kz = 2, Nz+1
           zz = (kz-1)*hz
!           greenwake(kz,1) = leng2*1.0d12*130.0d0*(1.075*exp(-sqrt(zz/2.25d-3))-&
!             0.075)
           greenwake(kz,1) = leng2*1.0d12*130.0d0*(1.02*exp(-sqrt(zz/2.0d-3))-&
             0.02)
         enddo
       else !for TWS
         !The 0.5 is from S+
         greenwake(1,1) = 0.5d0*Z0*Clight*exp(-sqrt(zz/zz00))/(pilc*aa*aa)
         !do kz = 1, Nz+1
         do kz = 2, Nz+1
           zz = (kz-1)*hz
           greenwake(kz,1) = Z0*Clight*exp(-sqrt(zz/zz00))/(pilc*aa*aa)
         enddo
       endif
       do kz = Nz+2, twonz
         greenwake(kz,1) = greenwake(twonz-kz+2,1)
       enddo
       do kz = 2, Nz
         greenwake(kz,1) = 0.0
       enddo

       call fftrclocal2_FFT(ksign,scale,greenwake,twonz,one,&
             greenwakeout)

       !do kz = 1, twonz
       !  greenwake(kz,1) = densz2nout(kz,1)*greenwakeout(kz,1)
       !enddo
       do kz = 1, 2
          greenwake(kz,1) = densz2nout(kz,1)*greenwakeout(kz,1)
       enddo
       do kz = 2, twonz/2
          greenwake(2*kz-1,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz-1,1)-&
                            densz2nout(2*kz,1)*greenwakeout(2*kz,1)
          greenwake(2*kz,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz,1)+&
                            densz2nout(2*kz,1)*greenwakeout(2*kz-1,1)
       enddo

       scale = 1.0d0/twonz
       ksign = -1
       call fftcrlocal2_FFT(ksign,scale,greenwake,twonz,one,&
            greenwakeout)

       !the "-" sign is from definition of longitudinal wake function, see Wangler's book
       do kz = 1, Nz
         ezwake(kz) = -greenwakeout(kz,1)*hz
       enddo
       
       if((flagbtw.eq.2) .or. (flagbtw.eq.3)) then !no transverse wake for standing wave cavity so far
         exwake = 0.0d0
         eywake = 0.0d0
         goto 100
       endif

       !The following is for comparison purpose
       !calculate the short range longitudinal wakefield by direct summation
!       do kz = 1, Nz
!         densz2n(kz,1) = 0.5*Z0*Clight/(pilc*aa*aa)*recvdensz(kz,1)*hz
!         do kkzz = kz+1, Nz
!           zz = (kkzz-kz)*hz
!       !    !densz2n(kz,1) = densz2n(kz,1) + 4*Z0*Clight*zz00*(1.0-(1.+sqrt(zz/zz00))*&
!       !    !                 exp(-sqrt(zz/zz00)))/(pilc*aa*aa*aa*aa)*recvdensz(kkzz,1)*&
!       !    !                 recvdensz(kkzz,2)*hz
!       !    !densz2n(kz,1) = densz2n(kz,1) + coef1*(zz/zz00)*recvdensz(kkzz,1)* &
!       !    !                                       recvdensz(kkzz,2)*hz
!       !    !densz2n(kz,1) = densz2n(kz,1) + coef1*(1.0-(1.+sqrt(zz/zz00))*&
!       !    !                 exp(-sqrt(zz/zz00)))*recvdensz(kkzz,1)*&
!       !    !                 recvdensz(kkzz,2)*hz
!           densz2n(kz,1) = densz2n(kz,1) + Z0*Clight*exp(-sqrt(zz/zz00))/(pilc*aa*aa)*& 
!                            recvdensz(kkzz,1)*hz
!       
!         enddo
!       enddo
    
       !write(17,*)erwake
       !write(17,*)densz2n(1:Nz,1)
       !call flush_(17)
!       do i = 1, Nz
!         zz = (i-1)*hz
!         tmptmp = Z0*Clight*exp(-sqrt(zz/zz00))/(pilc*aa*aa)
!         write(19,*)zz,-densz2n(i,1),ezwake(i),recvdensz(i,1),tmptmp
!       enddo
       !exwake(:) = densz2n(1:Nz,1)

!------------------------------------------------------------------------------
       !calculate the transverse wakefield effects
       do kz = 1, Nz
          densz2n(kz,1) = recvdensz(kz,1)*xwakez(kz)
       enddo
       do kz = Nz + 1, 2*Nz
          densz2n(kz,1) = 0.0d0
       enddo
 
       twonz = 2*Nz
       one = 1
       ksign = 1
       !scale = 1./(twonz)
       scale = 1.0d0
       call fftrclocal2_FFT(ksign,scale,densz2n,twonz,one,densz2nout)
 
       !tranverse wakefield
       !the following transverse wakefield function is from
       !"Short-range dipole wakefields in accelerating structures for the NLC",
       !K. L. F. Bane, SLAC-PUB-9663, 2003. here, 1.1 is an average as suggested       !on page 11 for cell 45.
       !for LCLS slac
       !zz00 = 1.1*0.169*aa*(aa/leng)**1.17*(gg/aa)**0.38
       !for Fermi Elettra
       zz00 = 1.0*0.169*aa*(aa/leng)**1.17*(gg/aa)**0.38
       coef1 = 4*Z0*Clight*zz00/(pilc*aa*aa*aa*aa)
       !densconst = 0.5e-6
       !offset = 0.001
       !zzmax = 0.002
       !coef2 = coef1/zz00*densconst*offset
       !print*,"aa wake: ",aa,leng,gg,zz00,Z0,coef1,coef2,offset,densconst
       ! parameter for Elettra BTW linac  
       pbtw1 = 2.8d4
       pbtw2 = 1.2d-4
       pbtw3 = 1.2d-4
       pbtw4 = 1.4d4
       if(flagbtw.eq.1) then !for Backward TWS
         greenwake(1,1) = 0.0
         do kz = 1, Nz+1
           zz = (kz-1)*hz
           greenwake(kz,1) = 1.0d12*(pbtw1*(1.0-(1.0+sqrt(zz/pbtw2))*exp(-sqrt(zz/pbtw3))) + &
                             pbtw4*sqrt(zz))
         enddo
       else ! for TWS
         greenwake(1,1) = 0.0
         do kz = 1, Nz+1
           zz = (kz-1)*hz
           greenwake(kz,1) = coef1*(1.0-(1.+sqrt(zz/zz00))*exp(-sqrt(zz/zz00)))
         enddo
       endif

       do kz = Nz+2, 2*Nz
         greenwake(kz,1) = greenwake(twonz-kz+2,1)
       enddo
       do kz = 1, Nz
         greenwake(kz,1) = 0.0
       enddo

       ksign = 1
       !scale = 1./(twonz)
       scale = 1.0d0
       call fftrclocal2_FFT(ksign,scale,greenwake,twonz,one,&
                greenwakeout)

       do kz = 1, 2
          greenwake(kz,1) = densz2nout(kz,1)*greenwakeout(kz,1)
       enddo
       do kz = 2, twonz/2
          greenwake(2*kz-1,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz-1,1)-&
                            densz2nout(2*kz,1)*greenwakeout(2*kz,1)
          greenwake(2*kz,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz,1)+&
                            densz2nout(2*kz,1)*greenwakeout(2*kz-1,1)
       enddo

       !scale = 1.0
       scale = 1.0d0/twonz
       ksign = -1
       call fftcrlocal2_FFT(ksign,scale,greenwake,twonz,one,&
                   densz2nout)

       do kz = 1, Nz
         exwake(kz) = densz2nout(kz,1)*hz
       enddo

       !The following is for comparison purpose
       !calculate the short range transverse wakefield by direct summation
       !do kz = 1, Nz
       !  densz2n(kz,1) = 0.0
       !  do kkzz = kz+1, Nz
       !    zz = (kkzz-kz)*hz
       !    !densz2n(kz,1) = densz2n(kz,1) + 4*Z0*Clight*zz00*(1.0-(1.+sqrt(zz/zz00))*&
       !    !                 exp(-sqrt(zz/zz00)))/(pilc*aa*aa*aa*aa)*recvdensz(kkzz,1)*&
       !    !                 recvdensz(kkzz,2)*hz
       !    !densz2n(kz,1) = densz2n(kz,1) + coef1*(zz/zz00)*recvdensz(kkzz,1)* &
       !    !                                       recvdensz(kkzz,2)*hz
       !    densz2n(kz,1) = densz2n(kz,1) + coef1*(1.0-(1.+sqrt(zz/zz00))*&
       !                     exp(-sqrt(zz/zz00)))*recvdensz(kkzz,1)*&
       !                     recvdensz(kkzz,2)*hz
       !
       !  enddo
       !enddo
     
       !write(17,*)erwake
       !write(17,*)densz2n(1:Nz,1)
       !call flush_(17)
       !do i = 1, Nz
       !  zz = (i-1)*hz
       !  write(19,*)zz,densz2n(i,1),coef2*(0.5*zz**2-zzmax*zz+0.5*zzmax**2),exwake(i)
       !enddo
       !exwake(:) = densz2n(1:Nz,1)

       !vertical "Y" wakefields
       do kz = 1, Nz
          densz2n(kz,1) = recvdensz(kz,1)*ywakez(kz)
       enddo
       do kz = Nz + 1, 2*Nz
          densz2n(kz,1) = 0.0d0
       enddo
 
       twonz = 2*Nz
       one = 1
       ksign = 1
       !scale = 1./(twonz)
       scale = 1.0d0
       call fftrclocal2_FFT(ksign,scale,densz2n,twonz,one,densz2nout)
 
       do kz = 1, 2
          greenwake(kz,1) = densz2nout(kz,1)*greenwakeout(kz,1)
       enddo
       do kz = 2, twonz/2
          greenwake(2*kz-1,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz-1,1)-&
                            densz2nout(2*kz,1)*greenwakeout(2*kz,1)
          greenwake(2*kz,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz,1)+&
                            densz2nout(2*kz,1)*greenwakeout(2*kz-1,1)
       enddo

       !scale = 1.0
       scale = 1.0d0/twonz
       ksign = -1
       call fftcrlocal2_FFT(ksign,scale,greenwake,twonz,one,&
                   greenwakeout)

       do kz = 1, Nz
         eywake(kz) = greenwakeout(kz,1)*hz
       enddo

100    continue

       end subroutine wakefieldana_FieldQuant

!-------------------------------------------------------------------------
       !longitudinal and transverse wakefield using read in discrete data
       subroutine wakefield_FieldQuant(Nz,xwakez,ywakez,recvdensz,exwake,eywake,ezwake,&
                                       hz,aa,gg,leng,flagbtw)
       implicit none
       include 'mpif.h'
       integer, intent(in) :: Nz,flagbtw
       double precision, intent(in) :: hz, aa, gg, leng
       double precision, dimension(Nz,2), intent(in) :: recvdensz
       double precision, dimension(Nz), intent(in) :: xwakez,ywakez
       double precision, dimension(Nz), intent(out) :: exwake,ezwake,eywake
       double precision, dimension(2*Nz,1) :: densz2n,densz2nout,&
                         greenwake,greenwakeout
       integer :: kz,twonz,one,ksign,kkzz,i,iwk,iwk1
       double precision :: scale,zz00,pilc,Z0,alpha1,alpha,zz
       double precision :: coef1,coef2,densconst,offset,zzmax
       real*8 :: tmptmp,pbtw1,pbtw2,pbtw3,pbtw4
       real*8 :: leng1,leng2,hwk,zwk,ab

       pilc = 2*dasin(1.0d0)
  

       if(flagbtw.eq.4) then !skip longitudinal wake
         ezwake = 0.0d0
         goto 50
       endif

       do kz = 1, Nz
          densz2n(kz,1) = recvdensz(kz,1) 
       enddo
       do kz = Nz + 1, 2*Nz
          densz2n(kz,1) = 0.0d0
       enddo

       twonz = 2*Nz
       one = 1
       ksign = 1
       scale = 1.0d0
       call fftrclocal2_FFT(ksign,scale,densz2n,twonz,one,densz2nout)

       hwk = zdatwk(2)-zdatwk(1) !assumed uniform grid 
       do kz = 1, Nz+1
         zz = (kz-1)*hz
         zwk = zz - zdatwk(1)
         iwk = zwk/hwk + 1
         iwk1 = iwk+1
         ab = iwk*hwk-zwk
         if(kz.eq.1) then
           greenwake(kz,1) = 0.5d0*(edatwk(iwk)*ab + edatwk(iwk1)*(1.0d0-ab))
         else
           greenwake(kz,1) = edatwk(iwk)*ab + edatwk(iwk1)*(1.0d0-ab)
         endif
       enddo
       
       do kz = Nz+2, twonz
         greenwake(kz,1) = greenwake(twonz-kz+2,1)
       enddo
       do kz = 2, Nz
         greenwake(kz,1) = 0.0
       enddo

       call fftrclocal2_FFT(ksign,scale,greenwake,twonz,one,&
             greenwakeout)

       !do kz = 1, twonz
       !  greenwake(kz,1) = densz2nout(kz,1)*greenwakeout(kz,1)
       !enddo
       do kz = 1, 2
          greenwake(kz,1) = densz2nout(kz,1)*greenwakeout(kz,1)
       enddo
       do kz = 2, twonz/2
          greenwake(2*kz-1,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz-1,1)-&
                            densz2nout(2*kz,1)*greenwakeout(2*kz,1)
          greenwake(2*kz,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz,1)+&
                            densz2nout(2*kz,1)*greenwakeout(2*kz-1,1)
       enddo

       scale = 1.0d0/twonz
       ksign = -1
       call fftcrlocal2_FFT(ksign,scale,greenwake,twonz,one,&
            greenwakeout)

       !the "-" sign is from definition of longitudinal wake function, see Wangler's book
       do kz = 1, Nz
         ezwake(kz) = -greenwakeout(kz,1)*hz
       enddo

50     continue
       
       if((flagbtw.eq.2) .or. (flagbtw.eq.3)) then !no transverse wake for standing wave cavity so far
         exwake = 0.0d0
         eywake = 0.0d0
         goto 100
       endif

!------------------------------------------------------------------------------
       !calculate the transverse wakefield effects
       do kz = 1, Nz
          densz2n(kz,1) = recvdensz(kz,1)*xwakez(kz)
       enddo
       do kz = Nz + 1, 2*Nz
          densz2n(kz,1) = 0.0d0
       enddo
 
       twonz = 2*Nz
       one = 1
       ksign = 1
       !scale = 1./(twonz)
       scale = 1.0d0
       call fftrclocal2_FFT(ksign,scale,densz2n,twonz,one,densz2nout)

       hwk = zdatwk(2)-zdatwk(1) !assumed uniform grid
       do kz = 1, Nz+1
         zz = (kz-1)*hz
         zwk = zz - zdatwk(1)
         iwk = zwk/hwk + 1
         iwk1 = iwk+1
         ab = iwk*hwk-zwk
         if(kz.eq.1) then
           greenwake(kz,1) = 0.5d0*(epdatwk(iwk)*ab + epdatwk(iwk1)*(1.0d0-ab))
         else
           greenwake(kz,1) = epdatwk(iwk)*ab + epdatwk(iwk1)*(1.0d0-ab)
         endif
       enddo
       do kz = Nz+2, 2*Nz
         greenwake(kz,1) = greenwake(twonz-kz+2,1)
       enddo
       do kz = 1, Nz
         greenwake(kz,1) = 0.0
       enddo

       ksign = 1
       !scale = 1./(twonz)
       scale = 1.0d0
       call fftrclocal2_FFT(ksign,scale,greenwake,twonz,one,&
                greenwakeout)

       do kz = 1, 2
          greenwake(kz,1) = densz2nout(kz,1)*greenwakeout(kz,1)
       enddo
       do kz = 2, twonz/2
          greenwake(2*kz-1,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz-1,1)-&
                            densz2nout(2*kz,1)*greenwakeout(2*kz,1)
          greenwake(2*kz,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz,1)+&
                            densz2nout(2*kz,1)*greenwakeout(2*kz-1,1)
       enddo

       !scale = 1.0
       scale = 1.0d0/twonz
       ksign = -1
       call fftcrlocal2_FFT(ksign,scale,greenwake,twonz,one,&
                   densz2nout)

       do kz = 1, Nz
         exwake(kz) = densz2nout(kz,1)*hz
       enddo

       !vertical "Y" wakefields
       do kz = 1, Nz
          densz2n(kz,1) = recvdensz(kz,1)*ywakez(kz)
       enddo
       do kz = Nz + 1, 2*Nz
          densz2n(kz,1) = 0.0d0
       enddo

       twonz = 2*Nz
       one = 1
       ksign = 1
       !scale = 1./(twonz)
       scale = 1.0d0
       call fftrclocal2_FFT(ksign,scale,densz2n,twonz,one,densz2nout)

       hwk = zdatwk(2)-zdatwk(1) !assumed uniform grid
       do kz = 1, Nz+1
         zz = (kz-1)*hz
         zwk = zz - zdatwk(1)
         iwk = zwk/hwk + 1
         iwk1 = iwk+1
         ab = iwk*hwk-zwk
         if(kz.eq.1) then
           greenwake(kz,1) = 0.5d0*(eppdatwk(iwk)*ab + eppdatwk(iwk1)*(1.0d0-ab))
         else
           greenwake(kz,1) = eppdatwk(iwk)*ab + eppdatwk(iwk1)*(1.0d0-ab)
         endif
       enddo
       do kz = Nz+2, 2*Nz
         greenwake(kz,1) = greenwake(twonz-kz+2,1)
       enddo
       do kz = 1, Nz
         greenwake(kz,1) = 0.0
       enddo

       ksign = 1
       !scale = 1./(twonz)
       scale = 1.0d0
       call fftrclocal2_FFT(ksign,scale,greenwake,twonz,one,&
                greenwakeout)

       do kz = 1, 2
          greenwake(kz,1) = densz2nout(kz,1)*greenwakeout(kz,1)
       enddo
       do kz = 2, twonz/2
          greenwake(2*kz-1,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz-1,1)-&
                            densz2nout(2*kz,1)*greenwakeout(2*kz,1)
          greenwake(2*kz,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz,1)+&
                            densz2nout(2*kz,1)*greenwakeout(2*kz-1,1)
       enddo

       !scale = 1.0
       scale = 1.0d0/twonz
       ksign = -1
       call fftcrlocal2_FFT(ksign,scale,greenwake,twonz,one,&
                   greenwakeout)

       do kz = 1, Nz
         eywake(kz) = greenwakeout(kz,1)*hz
       enddo

100    continue

       end subroutine wakefield_FieldQuant

!    Longitudinal CSR wakefields (1D, steady state) by I. Pogorelov
!    Two-level denoising, since we'll have Nz=32-128 in simulations 
!    where this will be tested. "Ndetail" to be used in ver.2 to specify 
!    on how many levels denoising should be performed. 
!    Nz has to be a power of 2.
        subroutine csrwakefield_FieldQuant(Nz, & 
                                  recvdensz,ezwakecsr,hzi,Rbend, ssll)
        implicit none 
        include 'mpif.h' 
        integer, intent(in):: Nz !, Ndetail 
        integer:: i, k 
        double precision, intent(in):: hzi, Rbend, ssll 
        double precision, dimension(1:Nz), intent(in):: recvdensz 
        double precision, dimension(1:Nz), intent(out):: ezwakecsr 
        double precision, dimension(:), allocatable:: uodd, ueven 
        double precision, dimension(-1:Nz+2):: u 
        double precision, dimension(1:Nz):: AA, BB 
        double precision, dimension(1:Nz-1):: kernel13  
        double precision :: hz
        integer :: Npts
        real*8 :: eps0,pilc

!------------------------------------------------------------
        !hz = hzi*Scxl 
        hz = hzi
        Npts = Nz
!        print*,"hz: ",hz
!
        u(:) = 0.0d0 
        u(1:Nz) = recvdensz(1:Nz) 
        allocate(uodd(-1:Nz/2+2), ueven(-1:Nz/2+2)) 
        ueven(:) = 0.0d0 
        uodd(:) = 0.0d0 
! 
        u(-1) = u(3) 
        u(0) = u(2) 
        u(Npts+1) = u(Npts-1) 
        u(Npts+2) = u(Npts-2) 
!
        uodd(1:Npts/2) = u(1:Npts-1:2) 
        ueven(1:Npts/2) = u(2:Npts:2) 
! 
        uodd(1:Npts/2) = (7.*(uodd(-1:Npts/2-2)+ uodd(3:Npts/2+2))  & 
                         +24.*(uodd(0:Npts/2-1) + uodd(2:Npts/2+1)) & 
                         +34.*uodd(1:Npts/2) )/96. 
        ueven(1:Npts/2) = (7.*(ueven(-1:Npts/2-2) +ueven(3:Npts/2+2))  &  
                          +24.*(ueven(0:Npts/2-1) +ueven(2:Npts/2+1))  & 
                          +34.*ueven(1:Npts/2) )/96. 
!
        u(1:Npts-1:2) = uodd(1:Npts/2) 
        u(2:Npts:2) = ueven(1:Npts/2) 
!
        deallocate(uodd, ueven) 
! 
        u(-1) = u(3) 
        u(0) = u(2) 
        u(Npts+1) = u(Npts-1) 
        u(Npts+2) = u(Npts-2) 
! 
! ("u" can be used in place of both AA and BB, if necessary) 
        AA(1:Npts) = (7.*(u(-1:Npts-2) +u(3:Npts+2)) & 
                     +24.*(u(0:Npts-1) +u(2:Npts+1)) & 
                     +34.*u(1:Npts) )/96. 
!
        u(1:Npts) = AA(1:Npts) 
        BB(1:Npts) = (u(-1:Npts-2) -8.*u(0:Npts-1)  &  
                   +8.*u(2:Npts+1) -u(3:Npts+2))/(12.*hz) 
! 
        do i = 1, Npts-1
          kernel13(i) = (hz*dble(Npts-i))**(-1./3.) 
        end do 
!
        ezwakecsr(:) = 0.0d0 
        do k = 3, Npts 
          ezwakecsr(k) = hz*dot_product(BB(1:k-1), kernel13(Npts-k+1:Npts-1))   & 
               -0.5*hz*(BB(1)*kernel13(Npts-k+1) +BB(k-1)*kernel13(Npts-1))  & 
               +1.5*BB(k)*hz**(2./3.) 
        end do 
!
        eps0 = 8.854187817d-12
        pilc = 2*asin(1.0d0)
        ezwakecsr(:) = -(2.d0*1.0d0/(3.d0*Rbend**2)**(1.d0/3.d0))*&
                        ezwakecsr(:)/(4*pilc*eps0) ! e=-1.0
!
!        do i = 1, Npts
!          write(15,*)(i-1)*hz,ezwakecsr(i),recvdensz(i)
!        enddo

        end subroutine csrwakefield_FieldQuant 

!------------------------------------------------------------------------
! calculate the steady state csr wake fields from the Fourier coeficients 
! of density data. Here, a threshold of the power spectral of the density
! has to be input. The 2nd derivatives of the density distribution is
! reconstructed based on the Fourier coefficients.
! J. Qiang, LBNL, 08/21/06
!------------------------------------------------------------------------
      subroutine csrFcoef_FieldQuant(zmin,zmax,edata,ndatareal,ncoefreal,&
                                     coeftol,r0,ezwake3) 
      implicit none
      integer :: ncoefreal,ndatareal
      real*8 :: zmin,zmax,r0,coeftol
      double precision, dimension(ndatareal) :: edata,ezwake3
      double precision, dimension(ncoefreal) :: Fcoef,Fcoef2
      double precision, dimension(ndatareal) :: rhopp
      integer :: i,j
      real*8 :: zlen,zhalf,zmid,h,pilc,zz,coefmax,tmpsum,tmpsump,tmpsumpp,tmpf1
      real*8 :: rk,rkzz,eps0

      zlen = zmax - zmin
      zhalf = zlen/2.0
      zmid = (zmax + zmin)/2
      h = zlen/(ndatareal-1)
      pilc = 2*asin(1.0d0)
!      print*,"The RF data number is: ",ndatareal,zlen,zmid,h

      rk = 2*pilc/zlen
      do j = 1, ncoefreal
        zz = zmin - zmid
        rkzz = (j-1)*zz*rk
        Fcoef(j) = (-0.5*edata(1)*cos(rkzz)*h)/zhalf
        Fcoef2(j) = (-0.5*edata(1)*sin(rkzz)*h)/zhalf
        zz = zmax - zmid
        rkzz = (j-1)*zz*rk
        Fcoef(j) = Fcoef(j)-(0.5*edata(ndatareal)*cos(rkzz)*h)&
                            /zhalf
        Fcoef2(j) = Fcoef2(j)-(0.5*edata(ndatareal)*sin(rkzz)*h)&
                            /zhalf
      enddo

      do i = 1, ndatareal
        zz = zmin+(i-1)*h - zmid
        do j = 1, ncoefreal
          rkzz = (j-1)*zz*rk
          Fcoef(j) = Fcoef(j) + (edata(i)*cos(rkzz)*h)/zhalf
          Fcoef2(j) = Fcoef2(j) + (edata(i)*sin(rkzz)*h)/zhalf
        enddo
      enddo

!      open(7,file="rfcoef.out",status="unknown")
!      do j = 1, ncoefreal
!        write(7,*)j,Fcoef(j),Fcoef2(j)
!      enddo
!      close(7)

      coefmax = 0.0
      do i = 1, ncoefreal
        if(coefmax.le.abs(Fcoef(i))) then
          coefmax = Fcoef(i)
        endif
        if(coefmax.le.abs(Fcoef2(i))) then
          coefmax = Fcoef2(i)
        endif
      enddo 
!      print*,"coefmax: ",coefmax

!      print*,"coeftol:",coeftol
      do i = 1, ncoefreal
        tmpf1 = abs(Fcoef(i)/coefmax)
        if(tmpf1.lt.coeftol) then
          Fcoef(i) = 0.0
        endif
        tmpf1 = abs(Fcoef2(i)/coefmax)
        if(tmpf1.lt.coeftol) then
          Fcoef2(i) = 0.0
        endif
      enddo

!      open(8,file="rhodata.out",status="unknown")
      do i = 1, ndatareal
        zz = zmin+(i-1)*h - zmid
!        tmpsum = 0.5*Fcoef(1)
!        tmpsump = 0.0
        tmpsumpp = 0.0
        do j = 2,ncoefreal
         rkzz = (j-1)*zz*rk
!         tmpsum = tmpsum + Fcoef(j)*cos(rkzz) + &
!                  Fcoef2(j)*sin(rkzz)
!         tmpsump = tmpsump-(j-1)*rk*(Fcoef(j)*sin(rkzz)-&
!                           Fcoef2(j)*cos(rkzz))
         tmpsumpp = tmpsumpp-((j-1)*rk)**2*(Fcoef(j)*cos(rkzz)+&
                                   Fcoef2(j)*sin(rkzz))
        enddo
        rhopp(i) = tmpsumpp
!        write(8,*)zmin+(i-1)*h,tmpsum,tmpsump,tmpsumpp,edata(i)
      enddo
!      close(8)

      do i = 1, ndatareal
          ezwake3(i)= 0.0d0
          do j = 1,i-1
            zz = (i-j)*h
            ezwake3(i) = ezwake3(i) + zz**(2.0d0/3.0d0)*rhopp(j)*h
          enddo
      enddo
      eps0 = 8.854187817d-12
      do i = 1, ndatareal
          ezwake3(i) = -1.5*ezwake3(i)*&
                    (2.d0*1.0d0/(3.d0*r0**2)**(1.d0/3.d0))/(4*pilc*eps0)
      enddo

        !print*,"zmin: ",zmin,zmax,h
        !do i = 1, Ndatareal
        !  write(14,*)zmin+(i-1)*h,ezwake3(i),edata(i)
        !enddo

      end subroutine csrFcoef_FieldQuant

        !this subroutine calculate the 1d csr wakefield including
        !entrance, stead-state, and transitions effects
        !here, hx, rho,...are in real unit. the return ezwake is V/m
        subroutine csrwakeTr_FieldQuant(Nx,r0,ptmin,hx,blength,rhonew,rhonewp,&
                             rhonewpp,ezwake)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: Nx
        real*8 :: r0,ptmin,hx,blength
        real*8, dimension(Nx) :: rhonew,rhonewp,rhonewpp,ezwake
        real*8 :: xx,xxl,xx2,dx,tmprho,epstol,deltas,deltasmax,xxbar,psi,&
                  phim,pilc,yy,hx2,xconst,xk,psitmp
        integer :: Ni,il,i,Nmax,j
        integer :: myidlc,ierr

        call MPI_COMM_RANK(MPI_COMM_WORLD, myidlc, ierr)

        epstol = 1.0d-10 !tolerance for root finding
        Nmax = 100 !maximum # of iteration in root finding
        pilc = 2*asin(1.0d0)
        phim = blength/r0
        xconst = 1./(4*pilc*8.854187817e-12)
        xk = -2.0d0/(3*r0**2)**(1.0d0/3.0d0)*xconst

        psitmp = 0.0
        ezwake(1) = 0.0 
        do i = 2, Nx
           xx = ptmin + (i-1)*hx 
           xxl = (xx/r0)**3*r0/24
           Ni = i
           ezwake(i) = 0.0
           if(xx.le.0) then
           else if(xx.le.blength) then
             if(xxl.ge.hx) then
               !analytical integration from si-sni-1 to si 
               !add the contribution from si-sni-1 to si
               !ezwake(i) = ezwake(i)+0.5*(xx-xx2)**(-1.0d0/3.0d0)*rhonewp(Ni-1)*hx&
               !          + 0.5*(2*hx**(-1.0d0/3.0d0)-(2*hx)**(-1.0d0/3.0d0))*&
               !            rhonewp(Ni)*hx
               ezwake(i) = ezwake(i) + 1.5d0*hx**(2.0d0/3.0d0)*(rhonewp(Ni-1)+&
                         0.5*rhonewpp(Ni-1)*hx)
               if(Ni.gt.2) then
                 il = 0
                 !trapzoid rule for integral from si-sl to si-sni-1
                 do j = 1, Ni-1
                   xx2 = ptmin + (j-1)*hx
                   if((xx-xx2).le.xxl) then
                     ezwake(i) = ezwake(i)+(xx-xx2)**(-1.0d0/3.0d0)*rhonewp(j)*hx
                   else
                     il = j
                   endif
                 enddo
                 il = il + 1
                 xx2 = ptmin + (il-1)*hx
                 ezwake(i) = ezwake(i)-0.5*(xx-xx2)**(-1.0d0/3.0d0)*rhonewp(il)*hx
                 xx2 = ptmin + (Ni-1-1)*hx
                 ezwake(i) = ezwake(i)-0.5*(xx-xx2)**(-1.0d0/3.0d0)*rhonewp(Ni-1)*hx
                 !add the contribtion from sl to xil within one hx
                 if(xxl.lt.(xx-ptmin)) then
                   xx2 = ptmin + (il-1)*hx
                   dx = xx2-(xx-xxl)
                   tmprho = rhonewp(il)+dx/hx*(rhonewp(il-1)-rhonewp(il))
                   ezwake(i) = ezwake(i)+dx*(rhonewp(il)*(xx-xx2)**(-1.0d0/3.0d0)+&
                               tmprho*xxl**(-1.0d0/3.0d0))/2
!                   tmprho = rhonew(il)+dx/hx*(rhonew(il-1)-rhonew(il))
!                   ezwake(i) = ezwake(i)+rhonew(il)/(xx-xx2)**(1.0d0/3.0d0) &
!                               -tmprho/xxl**(1.0d0/3.0d0)-0.5*dx*( &
!                                rhonew(il)/(xx-xx2)**(4.0d0/3.0d0)+ &
!                                tmprho/xxl**(4.0d0/3.0d0) )/3
                 endif
               endif
             else
             endif
             !add the entrance transient contribution
             if(xxl.ge.(xx-ptmin)) then
             else
               il = (xx-xxl-ptmin)/hx + 1
               dx = xx-xxl-(ptmin+(il-1)*hx)
               tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
               ezwake(i) = ezwake(i) + tmprho/xxl**(1.0d0/3.0d0)
             endif
             if(4*xxl.ge.(xx-ptmin)) then
             else
               il = (xx-4*xxl-ptmin)/hx + 1
               dx = xx-4*xxl-(ptmin+(il-1)*hx)
               tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
!               ezwake(i) = ezwake(i) - tmprho/xxl**(1.0d0/3.0d0)
             endif
             ezwake(i) = ezwake(i)*xk
             !if(myidlc.eq.0) print*,"xxl: ",xxl,xx,ezwake(i),r0,ptmin,phim

           else
             xxbar = (xx - blength)/r0
             deltasmax = (r0*phim**3/24.0d0)*(phim+4*xxbar)/(phim+xxbar)
             if(deltasmax.ge.hx) then
               !add the contribution from Ni-1 to Ni
               xx2 = ptmin + (Ni-1-1)*hx
               deltas = xx - xx2
               call psiroot(r0,xxbar,deltas,psi,epstol,Nmax)
               yy = psi + 2*xxbar
               ezwake(i) = ezwake(i)+0.5*yy*(3*yy**2-12*xxbar*yy+12*xxbar**2)&
                         /(yy-xxbar)**2*rhonewp(Ni-1)*psi*r0/24
               if(Ni.gt.2) then
                 il = 0
                 do j = 1, Ni-1
                   xx2 = ptmin + (j-1)*hx
                   deltas = xx - xx2
                   if(deltas .le. deltasmax) then
                     call psiroot(r0,xxbar,deltas,psi,epstol,Nmax)
                     ezwake(i) = ezwake(i) + rhonewp(j)/(psi+2*xxbar)*hx
                   else
                     il = j
                   endif
                 enddo
                 il = il + 1
                 !minus the half end points of integral to Ni-1
                 xx2 = ptmin + (il-1)*hx
                 deltas = xx - xx2
                 call psiroot(r0,xxbar,deltas,psi,epstol,Nmax)
                 ezwake(i) = ezwake(i)-0.5*rhonewp(il)/(psi+2*xxbar)*hx
                 psitmp = psi
                 xx2 = ptmin + (Ni-1-1)*hx
                 deltas = xx - xx2
                 call psiroot(r0,xxbar,deltas,psi,epstol,Nmax)
                 ezwake(i) = ezwake(i)-0.5*rhonewp(Ni-1)/(psi+2*xxbar)*hx
                 if(deltasmax.lt.xx-ptmin) then
                   xx2 = ptmin + (il-1)*hx
                   dx = xx2-(xx-deltasmax)
                   tmprho = rhonewp(il)+dx/hx*(rhonewp(il-1)-rhonewp(il))
                   ezwake(i) = ezwake(i)+dx*(rhonewp(il)/(psitmp+2*xxbar)+ &
                                tmprho/(phim+2*xxbar))/2
                 endif
               endif
             endif

             !add the transition effects
             if(deltasmax.ge.xx-ptmin) then
             else
               il = (xx-deltasmax-ptmin)/hx + 1
               dx = xx-deltasmax-(ptmin+(il-1)*hx)
               tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
               ezwake(i) = ezwake(i) + tmprho/(phim+2*xxbar)
             endif
             deltasmax = (r0*phim**2/6.0d0)*(phim+3*xxbar)
             if(deltasmax.ge.xx-ptmin) then
             else
               il = (xx-deltasmax-ptmin)/hx + 1
               dx = xx-deltasmax-(ptmin+(il-1)*hx)
               tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
               ezwake(i) = ezwake(i) - tmprho/(phim+2*xxbar)
             endif
             ezwake(i) = -ezwake(i)*4/r0*xconst
           endif
        enddo

        end subroutine csrwakeTr_FieldQuant

        !this subroutine calculate the 1d csr wakefield including
        !entrance, stead-state, and transitions effects
        !here, hx, rho,...are in real unit. the return ezwake is V/m
        subroutine csrwakeTr2_FieldQuant(Nx,r0,ptmin,hx,blength,rhonew,rhonewp,&
                             rhonewpp,gam,ezwake)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: Nx
        real*8 :: r0,ptmin,hx,blength,gam
        real*8, dimension(Nx) :: rhonew,rhonewp,rhonewpp,ezwake
        real*8 :: xx,xxl,xx2,dx,tmprho,epstol,deltas,deltasmax,xxbar,psi,&
                  phim,pilc,yy,hx2,xconst,xk,psitmp
        integer :: Ni,il,i,Nmax,j
        integer :: myidlc,ierr,islp
        real*8 :: aa,uuh,uuh24,xslp,csrss,xxp,ssh,xk2,tcoef
        real*8 :: bb,cc,yh,phih,phpy2,csrtr,xblg,xxbarh,phimh
        real*8 :: psip2x2,psipx2,csrdrm,csrdr,phpxy2

        call MPI_COMM_RANK(MPI_COMM_WORLD, myidlc, ierr)

        epstol = 1.0d-10 !tolerance for root finding
        Nmax = 100 !maximum # of iteration in root finding
        pilc = 2*asin(1.0d0)
        phim = blength/r0
        xconst = 1./(4*pilc*8.854187817e-12)
        xk = -2.0d0/(3*r0**2)**(1.0d0/3.0d0)*xconst

        xk2 = 4.0d0*gam**4/r0**2

        psitmp = 0.0
        ezwake(1) = 0.0 
        do i = 2, Nx
           xx = ptmin + (i-1)*hx 
           xxl = (xx/r0)**3*r0/24

           Ni = i
           ezwake(i) = 0.0

           xblg = (i-1)*hx

           if(xx.le.0) then
           else if(xx.le.blength) then

             xslp = (xx/r0)*r0/2/gam/gam + xxl
             ! case B to Si for all S' contributions
             ! steady state regime
             if(xslp >= xblg ) then 
               do j = 1, Ni
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Ni) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xxp = ptmin + (j-1)*hx
                 ssh = (xx-xxp)*gam**3/r0
                 !Eq. (44)
                 aa = sqrt(64.0+144.0*ssh**2)
                 uuh = (aa+12*ssh)**(1.0d0/3.0d0)-(aa-12*ssh)**(1.0d0/3.0d0)
                 uuh24 = uuh*uuh/4
                 !Eq. (32)
                 csrss = xk2*( (uuh24-1)/2/(1+uuh24)**3 + &
                  (1.0d0/6.0d0-2.0d0/9.0d0*uuh24-1.0d0/6.0d0*uuh24**2)/&
                  (1+uuh24)**3/(1+1.0d0/3.0d0*uuh24)**2 )
                 ezwake(i) = ezwake(i) + tcoef*rhonew(j)*csrss*hx
                 !if(myidlc.eq.0)print*,"csr1: ",i,j,csrss,ezwake(i)
               enddo

             ! case A to Si for all S' contributions
             else if(xslp < hx ) then 
               !excluding the contribution of S(Ni) to Si in 
               !renormalized Coulomb term
               do j = 1, Ni-1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Ni) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xxp = ptmin + (j-1)*hx
                 ssh = (xx-xxp)*gam**3/r0
                 phih = xx/r0*gam
                 bb = 2*phih+phih**3/3-2*ssh
                 cc = phih**2+phih**4/12-2*ssh*phih
                 yh = (-bb+sqrt(bb*bb-4*cc))/2
                 phpy2 = (phih+yh)**2
                 if(yh.gt.0.0) then
                   !Eq.30
                   csrtr = xk2*phpy2*((phpy2+phih**3*(3*phih/4+yh))/&
                         (phpy2+phih**4/4)**3 - &
                         1.0/(phpy2+phih**3/12*(phih+4*yh))**2 )
                 else
                   !csrtr = xk2*1.0d0/6.0d0
                   csrtr = 0.0d0
                 endif
                 ezwake(i) = ezwake(i) + tcoef*rhonew(j)*csrtr*hx
                 !if(myidlc.eq.0) print*,"cstr1: ",i,j,csrtr,ezwake(i)
               enddo

             ! case A to Si from S'min to Si-xslp then (+)
             ! case B to Si from Si-xslp to Si
             else
               islp = (xx-xslp-ptmin)/hx + 1
               !case A to Si, needs to add S' from 1 to Si-xslp
               do j = 1, islp
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Ni) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xxp = ptmin + (j-1)*hx
                 ssh = (xx-xxp)*gam**3/r0
                 phih = xx/r0*gam
                 bb = 2*phih+phih**3/3-2*ssh
                 cc = phih**2+phih**4/12-2*ssh*phih
                 yh = (-bb+sqrt(bb*bb-4*cc))/2
                 phpy2 = (phih+yh)**2
                 if(yh.gt.0.0) then
                   !Eq.30
                   csrtr = xk2*phpy2*((phpy2+phih**3*(3*phih/4+yh))/&
                         (phpy2+phih**4/4)**3 - &
                         1.0/(phpy2+phih**3/12*(phih+4*yh))**2 )
                 else
                   csrtr = 0.0d0
                 endif
                 ezwake(i) = ezwake(i) + tcoef*rhonew(j)*csrtr*hx
                 !if(myidlc.eq.0) print*,"cstr2: ",i,j,csrtr,ezwake(i)
               enddo

               !case B to Si, S' from Si-xslp to Si
               do j = islp+1, Ni
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Ni) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xxp = ptmin + (j-1)*hx
                 ssh = (xx-xxp)*gam**3/r0
                 !Eq. (44)
                 aa = sqrt(64.0+144.0*ssh**2)
                 uuh = (aa+12*ssh)**(1.0d0/3.0d0)-(aa-12*ssh)**(1.0d0/3.0d0)
                 uuh24 = uuh*uuh/4
                 !Eq. (32)
                 csrss = xk2*( (uuh24-1)/2/(1+uuh24)**3 + &
                  (1.0d0/6.0d0-2.0d0/9.0d0*uuh24-1.0d0/6.0d0*uuh24**2)/&
                  (1+uuh24)**3/(1+1.0d0/3.0d0*uuh24)**2 )
                 ezwake(i) = ezwake(i) + tcoef*rhonew(j)*csrss*hx
                 !if(myidlc.eq.0)print*,"csr2: ",i,j,csrss,ezwake(i)
               enddo
             endif

             !ezwake(i) = ezwake(i)*xk2
             ezwake(i) = ezwake(i)*xconst

             !if(myidlc.eq.0) print*,"xxl: ",xxl,xx,ezwake(i),r0,ptmin,phim
             !if(myidlc.eq.0) print*,"xxl: ",xx,ezwake(i),xxl,xx-ptmin,xslp,hx

           else

             xxbar = (xx - blength)
             xslp= (r0*phim+xxbar)/2/gam/gam + r0*phim**3/24*&
                   (r0*phim+4*xxbar)/(r0*phim+xxbar)
             if(xslp.ge.xblg) then
               !all points S' in case D
               xxbarh = xxbar*gam/r0
               phimh = phim*gam
               do j = 1,Ni-1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Ni) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xxp = ptmin + (j-1)*hx
                 ssh = (xx-xxp)*gam**3/r0
                 call psiroot2(r0,xxbar,ssh,psi,epstol,Nmax)
                 psipx2 = (psi+xxbarh)**2
                 psip2x2 = (psi+2*xxbarh)**2
                 if(psi.gt.0.0 .and. psi.le.phimh) then
                   !Eq. 36
                   csrdr = xk2*psipx2*(psi**2*(psi**2/4*psip2x2-psipx2)/&
                         (2*(psipx2+psi**2/4*psip2x2)**3)+&
                         (psipx2+psi**2*(3*psi**2/4+xxbarh**2+2*xxbarh*psi))/&
                         (psipx2+psi**2/4*psip2x2)**3-1.0d0/&
                         (psipx2+psi**2/12*(psi**2+4*xxbarh*psi))**2)
                 else
                   csrdr = 0.0d0
                 endif
 
                 ezwake(i) = ezwake(i) + tcoef*rhonew(j)*csrdr*hx
               enddo
             else if(xslp.lt.hx) then
               !all points S' in case C
               xxbarh = xxbar*gam/r0
               phimh = phim*gam
               !excluding the contribution of S(Ni) to Si in 
               !renormalized Coulomb term
               do j = 1,Ni-1 
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Ni) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xxp = ptmin + (j-1)*hx
                 ssh = (xx-xxp)*gam**3/r0
                 bb = 2*(phimh+xxbarh)-2*ssh+phimh**3/3+phimh**2*xxbarh
                 cc = (phimh+xxbarh)**2+phimh**2*(phimh**2+4*phimh*xxbarh)/12-&
                      2*ssh*(phimh+xxbarh)
                 yh = (-bb+sqrt(bb*bb-4*cc))/2
                 phpxy2 = (phimh+xxbarh+yh)**2
                 !Eq. 34
                 csrdrm = xk2*phpxy2*((phpxy2+phimh**2*(3*phimh**2/4+&
                          xxbarh**2+2*xxbarh*phimh+yh*phimh+2*xxbarh*yh))/&
                          (phpxy2+phimh**2/4*(phimh+2*xxbarh)**2)**3-&
                          1.0/(phpxy2+phimh**2/12*(phimh**2+4*xxbarh*phimh+&
                          4*yh*phimh+12*xxbarh*yh))**2)
                 ezwake(i) = ezwake(i) + tcoef*rhonew(j)*csrdrm*hx
               enddo

             else 
               !points S' from 1 to S-xslp in case C
               islp = (xx-xslp-ptmin)/hx + 1
               xxbarh = xxbar*gam/r0
               phimh = phim*gam
               do j = 1,islp
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Ni) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xxp = ptmin + (j-1)*hx
                 ssh = (xx-xxp)*gam**3/r0
                 bb = 2*(phimh+xxbarh)-2*ssh+phimh**3/3+phimh**2*xxbarh
                 cc = (phimh+xxbarh)**2+phimh**2*(phimh**2+4*phimh*xxbarh)/12-&
                      2*ssh*(phimh+xxbarh)
                 yh = (-bb+sqrt(bb*bb-4*cc))/2
                 phpxy2 = (phimh+xxbarh+yh)**2
                 !Eq. 34
                 csrdrm = xk2*phpxy2*((phpxy2+phimh**2*(3*phimh**2/4+&
                          xxbarh**2+2*xxbarh*phimh+yh*phimh+2*xxbarh*yh))/&
                          (phpxy2+phimh**2/4*(phimh+2*xxbarh)**2)**3-&
                          1.0/(phpxy2+phimh**2/12*(phimh**2+4*xxbarh*phimh+&
                          4*yh*phimh+12*xxbarh*yh))**2)
                 ezwake(i) = ezwake(i) + tcoef*rhonew(j)*csrdrm*hx
               enddo

               !points S' from S-xslp to S in case D
               do j = islp+1,Ni-1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Ni) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xxp = ptmin + (j-1)*hx
                 ssh = (xx-xxp)*gam**3/r0
                 call psiroot2(r0,xxbar,ssh,psi,epstol,Nmax)
                 psipx2 = (psi+xxbarh)**2
                 psip2x2 = (psi+2*xxbarh)**2
                 if(psi.gt.0.0 .and. psi.le.phimh) then
                   !Eq. 36
                   csrdr = xk2*psipx2*(psi**2*(psi**2/4*psip2x2-psipx2)/&
                         (2*(psipx2+psi**2/4*psip2x2)**3)+&
                         (psipx2+psi**2*(3*psi**2/4+xxbarh**2+2*xxbarh*psi))/&
                         (psipx2+psi**2/4*psip2x2)**3-1.0d0/&
                         (psipx2+psi**2/12*(psi**2+4*xxbarh*psi))**2)
                 else
                   csrdr = 0.0d0
                 endif
                 ezwake(i) = ezwake(i) + tcoef*rhonew(j)*csrdr*hx
               enddo
             endif

           endif
           ezwake(i) = ezwake(i)*xconst
        enddo

        end subroutine csrwakeTr2_FieldQuant

      !find the psi in equation 12 of Stupakov and Emma's paper
      subroutine psiroot(r0,xx,deltas,psi,eps,Nmax)
      implicit none
      integer :: Nmax
      real*8 :: r0, xx, deltas,eps,psi
      integer :: i
      real*8 :: ps0,ps1,fps0,dfps0
 
 
      ps0 = (24*deltas/r0)**(1.0d0/3.0d0)
 
      do i = 1, Nmax
        fps0 = r0*ps0**4+4*xx*r0*ps0**3-24*deltas*ps0-24*deltas*xx
        !print*,"ps0: ",i,ps0,fps0
        if(abs(fps0).le.eps) then
          psi = ps0
          return
        else
          dfps0 = 4*r0*ps0**3+12*xx*r0*ps0**2-24*deltas
          ps1 = ps0 - fps0/dfps0
          ps0 = ps1
        endif
      enddo

      if(i.ge.Nmax) then 
        print*,"Not converged in psiroot" 
        stop
      endif 

      end subroutine psiroot

      !find the psi in equation 37 of Saldin et al. paper
      subroutine psiroot2(r0,xx,deltas,psi,eps,Nmax)
      implicit none
      integer :: Nmax
      real*8 :: r0, xx, deltas,eps,psi
      integer :: i
      real*8 :: ps0,ps1,fps0,dfps0
 
 
      ps0 = (6*deltas-3*xx)**(1.0d0/3.0d0)
      if(ps0.lt.0.0) ps0 = 0.0d0
 
      do i = 1, Nmax
        fps0 = ps0**4/12+xx*ps0**3/3+ps0**2+(2*xx-2*deltas)*ps0-&
               2*deltas*xx+xx*xx

        !print*,"ps0: ",i,ps0,fps0
        if(abs(fps0).le.eps) then
          psi = ps0
          return
        else
          dfps0 = ps0**3/3+xx*ps0**2+2*ps0+2*xx-2*deltas
          ps1 = ps0 - fps0/dfps0
          ps0 = ps1
        endif
      enddo

      if(i.ge.Nmax) then
        print*,"Not converged in psiroot2"
        stop
      endif

      end subroutine psiroot2

        !this subroutine calculate the 1d csr wakefield including
        !entrance, stead-state, and transitions effects
        !here, hx, rho,...are in real unit. the return ezwake is V/m
        subroutine csrwakeTr3old_FieldQuant(Nx,r0,ptmin,hx,blength,rhonew,rhonewp,&
                             rhonewpp,ezwake)
        implicit none
        integer, intent(in) :: Nx
        real*8 :: r0,ptmin,hx,blength
        real*8, dimension(Nx) :: rhonew,rhonewp,rhonewpp,ezwake
        real*8 :: xx,xxl,xx2,dx,tmprho,epstol,deltas,deltasmax,xxbar,psi,&
                  phim,pilc,yy,hx2,xconst,xk,psitmp
        integer :: Ni,il,i,Nmax,j,Nsup,islp,Nisup,jj0,jj,jj1,isup
        real*8 :: tcoef,hxsup,csrss,xxp0,xxsup,rhoxxp,xblg

        epstol = 1.0d-10 !tolerance for root finding
        Nmax = 100 !maximum # of iteration in root finding
        pilc = 2*asin(1.0d0)
        phim = blength/r0
        xconst = 1./(4*pilc*8.854187817e-12)
        xk = -2.0d0/(3*r0**2)**(1.0d0/3.0d0)*xconst

        Nsup = 10

        psitmp = 0.0
        ezwake(1) = 0.0 
        do i = 2, Nx
           xx = ptmin + (i-1)*hx 
           xxl = (xx/r0)**3*r0/24
           Ni = i
           ezwake(i) = 0.0
           xblg = (i-1)*hx

           print*,"xxl: ",i,xxl,xblg
           if(xx.le.0) then
           else if(xx.le.blength) then
             if(xxl.ge.xblg) then !S-S
               do j = 1, Ni-1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Ni-1) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xx2 = ptmin + (j-1)*hx
                 !csrss = (xx-xx2)**(-1.0d0/3.0d0)
                 ezwake(i) = ezwake(i)+tcoef*&
                            (xx-xx2)**(-1.0d0/3.0d0)*rhonewp(j)*hx
               enddo

               !subgrid integration from Ni-1 to Ni
               hxsup = hx/Nsup
               isup = 1
               jj0 = Ni-isup
               xxp0 = ptmin+(jj0-1)*hx
 
               Nisup = isup*Nsup
 
               do j = 1, Nisup+1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Nisup+1) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
 
                 xx2 = xxp0 + (j-1)*hxsup
                 if(j.ne.Nisup+1) then
                   csrss = (xx-xx2)**(-1.0d0/3.0d0)
                 else
                   csrss = 2*(hxsup)**(-1.0d0/3.0d0)-&
                           (2*hxsup)**(-1.0d0/3.0d0)
                 endif
                 jj = jj0 
                 jj1 = jj0+1
                 xxsup = xx2 - (ptmin+(jj-1)*hx)
                 rhoxxp = rhonewp(jj)+(rhonewp(jj1)-rhonewp(jj))*xxsup/hx
                 ezwake(i) = ezwake(i)+tcoef*csrss*rhoxxp*hxsup
               enddo
             else if(xxl.lt.hx) then !Transient
               if(4*xxl.ge.(xx-ptmin)) then
               else 
                 il = (xx-4*xxl-ptmin)/hx + 1  
                 dx = xx-4*xxl-(ptmin+(il-1)*hx)
                 tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
                 ezwake(i) = ezwake(i) - tmprho/xxl**(1.0d0/3.0d0)   
               endif
             else
               islp = (xx-xxl-ptmin)/hx + 1

               !transient contribution from 0 to xxl 
               if(4*xxl.ge.(xx-ptmin)) then
               else
                 il = (xx-4*xxl-ptmin)/hx + 1                 
                 dx = xx-4*xxl-(ptmin+(il-1)*hx) 
                 tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
                 ezwake(i) = ezwake(i) - tmprho/xxl**(1.0d0/3.0d0)
               endif 

!               print*,"islp: ",islp,xxl,xx,hx

               ezwake(i) = ezwake(i) + rhonew(islp)/xxl**(1.0d0/3.0d0)
               do j = islp, Ni-1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Ni-1) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xx2 = ptmin + (j-1)*hx
                 !csrss = (xx-xx2)**(-1.0d0/3.0d0)
                 ezwake(i) = ezwake(i)+tcoef*&
                            (xx-xx2)**(-1.0d0/3.0d0)*rhonewp(j)*hx
               enddo
               !subgrid integration from Ni-1 to Ni
               hxsup = hx/Nsup
               isup = 1
               jj0 = Ni-isup
               xxp0 = ptmin+(jj0-1)*hx
               Nisup = isup*Nsup
               do j = 1, Nisup+1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Nisup+1) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xx2 = xxp0 + (j-1)*hxsup
                 if(j.ne.Nisup+1) then
                   csrss = (xx-xx2)**(-1.0d0/3.0d0)
                 else
                   csrss = 2*(hxsup)**(-1.0d0/3.0d0)-&
                           (2*hxsup)**(-1.0d0/3.0d0)
                 endif
                 jj = jj0
                 jj1 = jj0+1
                 xxsup = xx2 - (ptmin+(jj-1)*hx)
                 rhoxxp = rhonewp(jj)+(rhonewp(jj1)-rhonewp(jj))*xxsup/hx
                 ezwake(i) = ezwake(i)+tcoef*csrss*rhoxxp*hxsup
               enddo
             endif
             ezwake(i) = ezwake(i)*xk
           else
             xxbar = (xx - blength)/r0
             deltasmax = (r0*phim**3/24.0d0)*(phim+4*xxbar)/(phim+xxbar)
             if(deltasmax.ge.hx) then
               !add the contribution from Ni-1 to Ni
               xx2 = ptmin + (Ni-1-1)*hx
               deltas = xx - xx2
               call psiroot(r0,xxbar,deltas,psi,epstol,Nmax)
               yy = psi + 2*xxbar
               ezwake(i) = ezwake(i)+0.5*yy*(3*yy**2-12*xxbar*yy+12*xxbar**2)&
                         /(yy-xxbar)**2*rhonewp(Ni-1)*psi*r0/24
               if(Ni.gt.2) then
                 il = 0
                 do j = 1, Ni-1
                   xx2 = ptmin + (j-1)*hx
                   deltas = xx - xx2
                   if(deltas .le. deltasmax) then
                     call psiroot(r0,xxbar,deltas,psi,epstol,Nmax)
                     ezwake(i) = ezwake(i) + rhonewp(j)/(psi+2*xxbar)*hx
                   else
                     il = j
                   endif
                 enddo
                 il = il + 1
                 !minus the half end points of integral to Ni-1
                 xx2 = ptmin + (il-1)*hx
                 deltas = xx - xx2
                 call psiroot(r0,xxbar,deltas,psi,epstol,Nmax)
                 ezwake(i) = ezwake(i)-0.5*rhonewp(il)/(psi+2*xxbar)*hx
                 psitmp = psi
                 xx2 = ptmin + (Ni-1-1)*hx
                 deltas = xx - xx2
                 call psiroot(r0,xxbar,deltas,psi,epstol,Nmax)
                 ezwake(i) = ezwake(i)-0.5*rhonewp(Ni-1)/(psi+2*xxbar)*hx
                 if(deltasmax.lt.xx-ptmin) then
                   xx2 = ptmin + (il-1)*hx
                   dx = xx2-(xx-deltasmax)
                   tmprho = rhonewp(il)+dx/hx*(rhonewp(il-1)-rhonewp(il))
                   ezwake(i) = ezwake(i)+dx*(rhonewp(il)/(psitmp+2*xxbar)+ &
                                tmprho/(phim+2*xxbar))/2
                 endif
               endif
             endif

             !add the transition effects
             if(deltasmax.ge.xx-ptmin) then
             else
               il = (xx-deltasmax-ptmin)/hx + 1
               dx = xx-deltasmax-(ptmin+(il-1)*hx)
               tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
               ezwake(i) = ezwake(i) + tmprho/(phim+2*xxbar)
             endif
             deltasmax = (r0*phim**2/6.0d0)*(phim+3*xxbar)
             if(deltasmax.ge.xx-ptmin) then
             else
               il = (xx-deltasmax-ptmin)/hx + 1
               dx = xx-deltasmax-(ptmin+(il-1)*hx)
               tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
               ezwake(i) = ezwake(i) - tmprho/(phim+2*xxbar)
             endif
             ezwake(i) = -ezwake(i)*4/r0*xconst
           endif
        enddo

        end subroutine csrwakeTr3old_FieldQuant

        !this subroutine calculate the 1d csr wakefield including
        !entrance, stead-state, and transitions effects
        !here, hx, rho,...are in real unit. the return ezwake is V/m
        subroutine csrwakeTr3old2_FieldQuant(Nx,r0,ptmin,hx,blength,rhonew,rhonewp,&
                             rhonewpp,ezwake)
        implicit none
        integer, intent(in) :: Nx
        real*8 :: r0,ptmin,hx,blength
        real*8, dimension(Nx) :: rhonew,rhonewp,rhonewpp,ezwake
        real*8 :: xx,xxl,xx2,dx,tmprho,epstol,deltas,deltasmax,xxbar,psi,&
                  phim,pilc,yy,hx2,xconst,xk,psitmp
        integer :: Ni,il,i,Nmax,j,Nsup,islp,Nisup,jj0,jj,jj1,isup
        real*8 :: tcoef,hxsup,csrss,xxp0,xxsup,rhoxxp,xblg,deltasmax2

        epstol = 1.0d-10 !tolerance for root finding
        Nmax = 100 !maximum # of iteration in root finding
        pilc = 2*asin(1.0d0)
        phim = blength/r0
        xconst = 1./(4*pilc*8.854187817e-12)
        xk = -2.0d0/(3*r0**2)**(1.0d0/3.0d0)*xconst

        Nsup = 10

        psitmp = 0.0
        ezwake(1) = 0.0 
        do i = 2, Nx
           xx = ptmin + (i-1)*hx 
           xxl = (xx/r0)**3*r0/24
           Ni = i
           ezwake(i) = 0.0
           xblg = (i-1)*hx

           !print*,"xxl: ",i,xx,blength,xxl,xblg
           if(xx.le.0) then
           else if(xx.le.blength) then
             if(xxl.ge.xblg) then !S-S
               do j = 1, Ni-1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Ni-1) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xx2 = ptmin + (j-1)*hx
                 !csrss = (xx-xx2)**(-1.0d0/3.0d0)
                 ezwake(i) = ezwake(i)+tcoef*&
                            (xx-xx2)**(-1.0d0/3.0d0)*rhonewp(j)*hx
               enddo

               !subgrid integration from Ni-1 to Ni
               hxsup = hx/Nsup
               isup = 1
               jj0 = Ni-isup
               xxp0 = ptmin+(jj0-1)*hx
 
               Nisup = isup*Nsup
 
               do j = 1, Nisup+1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Nisup+1) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
 
                 xx2 = xxp0 + (j-1)*hxsup
                 if(j.ne.Nisup+1) then
                   csrss = (xx-xx2)**(-1.0d0/3.0d0)
                 else
                   csrss = 2*(hxsup)**(-1.0d0/3.0d0)-&
                           (2*hxsup)**(-1.0d0/3.0d0)
                 endif
                 jj = jj0 
                 jj1 = jj0+1
                 xxsup = xx2 - (ptmin+(jj-1)*hx)
                 rhoxxp = rhonewp(jj)+(rhonewp(jj1)-rhonewp(jj))*xxsup/hx
                 ezwake(i) = ezwake(i)+tcoef*csrss*rhoxxp*hxsup
               enddo
             else if(xxl.lt.hx) then !Transient
               if(4*xxl.ge.(xx-ptmin)) then
               else 
                 il = (xx-4*xxl-ptmin)/hx + 1  
                 dx = xx-4*xxl-(ptmin+(il-1)*hx)
                 tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
                 ezwake(i) = ezwake(i) - tmprho/xxl**(1.0d0/3.0d0)   
               endif
             else
               islp = (xx-xxl-ptmin)/hx + 1

               !transient contribution from 0 to xbl-xxl 
               if(4*xxl.ge.(xx-ptmin)) then
               else
                 il = (xx-4*xxl-ptmin)/hx + 1                 
                 dx = xx-4*xxl-(ptmin+(il-1)*hx) 
                 tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
                 ezwake(i) = ezwake(i) - tmprho/xxl**(1.0d0/3.0d0)
               endif 

!               print*,"islp: ",islp,xxl,xx,hx

               !case B
               ezwake(i) = ezwake(i) + rhonew(islp)/xxl**(1.0d0/3.0d0)
               do j = islp, Ni-1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Ni-1) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xx2 = ptmin + (j-1)*hx
                 !csrss = (xx-xx2)**(-1.0d0/3.0d0)
                 ezwake(i) = ezwake(i)+tcoef*&
                            (xx-xx2)**(-1.0d0/3.0d0)*rhonewp(j)*hx
               enddo
               !subgrid integration from Ni-1 to Ni
               hxsup = hx/Nsup
               isup = 1
               jj0 = Ni-isup
               xxp0 = ptmin+(jj0-1)*hx
               Nisup = isup*Nsup
               do j = 1, Nisup+1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Nisup+1) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xx2 = xxp0 + (j-1)*hxsup
                 if(j.ne.Nisup+1) then
                   csrss = (xx-xx2)**(-1.0d0/3.0d0)
                 else
                   csrss = 2*(hxsup)**(-1.0d0/3.0d0)-&
                           (2*hxsup)**(-1.0d0/3.0d0)
                 endif
                 jj = jj0
                 jj1 = jj0+1
                 xxsup = xx2 - (ptmin+(jj-1)*hx)
                 rhoxxp = rhonewp(jj)+(rhonewp(jj1)-rhonewp(jj))*xxsup/hx
                 ezwake(i) = ezwake(i)+tcoef*csrss*rhoxxp*hxsup
               enddo
             endif
             ezwake(i) = ezwake(i)*xk
           else
             xxbar = (xx - blength)/r0
             deltasmax = (r0*phim**3/24.0d0)*(phim+4*xxbar)/(phim+xxbar)
             !print*,"xxbar: ",i,xx,xxbar,deltasmax,xblg
             if(deltasmax.ge.xblg) then !case D
               !neglect the contributions from s to s
               do j = 1, Ni-1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Ni) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xx2 = ptmin + (j-1)*hx
                 deltas = xx - xx2
                 call psiroot(r0,xxbar,deltas,psi,epstol,Nmax)
                 ezwake(i) = ezwake(i) + tcoef*rhonewp(j)/(psi+2*xxbar)*hx
               enddo
             else if(deltasmax.lt.hx) then !case C
               deltasmax2 = (r0*phim**2/6.0d0)*(phim+3*xxbar)
               if(deltasmax2.ge.xx-ptmin) then
               else
                 il = (xx-deltasmax2-ptmin)/hx + 1
                 dx = xx-deltasmax2-(ptmin+(il-1)*hx)
                 tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
                 ezwake(i) = ezwake(i) - tmprho/(phim+2*xxbar)
               endif
             else
               islp = (xx-deltasmax-ptmin)/hx + 1

               !case C from 0 to xblg-deltasmax
               deltasmax2 = (r0*phim**2/6.0d0)*(phim+3*xxbar)
               if(deltasmax2.ge.xx-ptmin) then
               else
                 il = (xx-deltasmax2-ptmin)/hx + 1
                 dx = xx-deltasmax2-(ptmin+(il-1)*hx)
                 tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
                 ezwake(i) = ezwake(i) - tmprho/(phim+2*xxbar)
               endif 

               !case D
               ezwake(i) = ezwake(i) + rhonew(islp)/(phim+2*xxbar)
               do j = islp, Ni-1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Ni) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xx2 = ptmin + (j-1)*hx
                 deltas = xx - xx2
                 call psiroot(r0,xxbar,deltas,psi,epstol,Nmax)
                 ezwake(i) = ezwake(i) + tcoef*rhonewp(j)/(psi+2*xxbar)*hx
               enddo
             endif
             ezwake(i) = -ezwake(i)*4/r0*xconst
           endif
        enddo

        end subroutine csrwakeTr3old2_FieldQuant

        !this subroutine calculate the 1d csr wakefield including
        !entrance, stead-state, and transitions effects
        !here, hx, rho,...are in real unit. the return ezwake is V/m
        subroutine csrwakeTr3_FieldQuant(Nx,r0,ptmin,hx,blength,rhonew,rhonewp,&
                             rhonewpp,ezwake)
        implicit none
        integer, intent(in) :: Nx
        real*8 :: r0,ptmin,hx,blength
        real*8, dimension(Nx) :: rhonew,rhonewp,rhonewpp,ezwake
        real*8 :: xx,xxl,xx2,dx,tmprho,epstol,deltas,deltasmax,xxbar,psi,&
                  phim,pilc,yy,hx2,xconst,xk,psitmp
        integer :: Ni,il,i,Nmax,j,Nsup,islp,Nisup,jj0,jj,jj1,isup
        real*8 :: tcoef,hxsup,csrss,xxp0,xxsup,rhoxxp,xblg,deltasmax2

        epstol = 1.0d-10 !tolerance for root finding
        Nmax = 100 !maximum # of iteration in root finding
        pilc = 2*asin(1.0d0)
        phim = blength/r0
        xconst = 1./(4*pilc*8.854187817e-12)
        xk = -2.0d0/(3*r0**2)**(1.0d0/3.0d0)*xconst

        Nsup = 10

        psitmp = 0.0
        ezwake(1) = 0.0 
        do i = 2, Nx
           xx = ptmin + (i-1)*hx 
           xxl = (xx/r0)**3*r0/24
           Ni = i
           ezwake(i) = 0.0
           xblg = (i-1)*hx

           !print*,"xxl: ",i,xx,blength,xxl,xblg
           if(xx.le.0) then
           else if(xx.le.blength) then
             if(xxl.ge.xblg) then !S-S
               do j = 1, Ni-1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Ni-1) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xx2 = ptmin + (j-1)*hx
                 !csrss = (xx-xx2)**(-1.0d0/3.0d0)
                 ezwake(i) = ezwake(i)+tcoef*&
                            (xx-xx2)**(-1.0d0/3.0d0)*rhonewp(j)*hx
               enddo

               !subgrid integration from Ni-1 to Ni
               hxsup = hx/Nsup
               isup = 1
               jj0 = Ni-isup
               xxp0 = ptmin+(jj0-1)*hx
 
               Nisup = isup*Nsup
 
               do j = 1, Nisup+1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Nisup+1) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
 
                 xx2 = xxp0 + (j-1)*hxsup
                 if(j.ne.Nisup+1) then
                   csrss = (xx-xx2)**(-1.0d0/3.0d0)
                 else
                   csrss = 2*(hxsup)**(-1.0d0/3.0d0)-&
                           (2*hxsup)**(-1.0d0/3.0d0)
                 endif
                 jj = jj0 
                 jj1 = jj0+1
                 xxsup = xx2 - (ptmin+(jj-1)*hx)
                 rhoxxp = rhonewp(jj)+(rhonewp(jj1)-rhonewp(jj))*xxsup/hx
                 ezwake(i) = ezwake(i)+tcoef*csrss*rhoxxp*hxsup
               enddo
             else if(xxl.lt.hx) then !Transient
               if(4*xxl.ge.(xx-ptmin)) then
               else 
                 il = (xx-4*xxl-ptmin)/hx + 1  
                 dx = xx-4*xxl-(ptmin+(il-1)*hx)
                 tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
                 ezwake(i) = ezwake(i) - tmprho/xxl**(1.0d0/3.0d0)   
               endif

               !contribution from case B
               il = Ni-1
               dx = xx-xxl-(ptmin+(il-1)*hx)
               tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
               ezwake(i) = ezwake(i) + tmprho/xxl**(1.0d0/3.0d0)

               hxsup = xxl/Nsup
               jj0 = Ni-1
               xxp0 = xx - xxl
               Nisup = Nsup
               jj = jj0
               jj1 = jj0+1
               do j = 1, Nisup+1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Nisup+1) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif

                 xx2 = xxp0 + (j-1)*hxsup
                 if(j.ne.Nisup+1) then
                   csrss = (xx-xx2)**(-1.0d0/3.0d0)
                 else
                   csrss = 2*(hxsup)**(-1.0d0/3.0d0)-&
                           (2*hxsup)**(-1.0d0/3.0d0)
                 endif

                 xxsup = xx2 - (ptmin+(jj-1)*hx)
                 rhoxxp = rhonewp(jj)+(rhonewp(jj1)-rhonewp(jj))*xxsup/hx
                 ezwake(i) = ezwake(i)+tcoef*csrss*rhoxxp*hxsup
               enddo

             else
               islp = (xx-xxl-ptmin)/hx + 1

               !transient contribution from 0 to xbl-xxl 
               if(4*xxl.ge.(xx-ptmin)) then
               else
                 il = (xx-4*xxl-ptmin)/hx + 1                 
                 dx = xx-4*xxl-(ptmin+(il-1)*hx) 
                 tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
                 ezwake(i) = ezwake(i) - tmprho/xxl**(1.0d0/3.0d0)
               endif 

!               print*,"islp: ",islp,xxl,xx,hx

               !case B

               !from islp to islp + 1 
               isup = xxl/hx
               hxsup = (xxl-isup*hx)/Nsup
               jj0 = islp
               xxp0 = xx - xxl
               Nisup = Nsup
               dx = xxp0 - (ptmin+(jj0-1)*hx)
               tmprho = rhonew(jj0) + dx/hx*(rhonew(jj0+1)-rhonew(jj0))
               ezwake(i) = ezwake(i) + tmprho/xxl**(1.0d0/3.0d0)
               jj = jj0
               jj1 = jj0+1
               do j = 1, Nisup+1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Nisup+1) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
 
                 xx2 = xxp0 + (j-1)*hxsup
                 csrss = (xx-xx2)**(-1.0d0/3.0d0)
 
                 xxsup = xx2 - (ptmin+(jj-1)*hx)
                 rhoxxp = rhonewp(jj)+(rhonewp(jj1)-rhonewp(jj))*xxsup/hx
                 ezwake(i) = ezwake(i)+tcoef*csrss*rhoxxp*hxsup
               enddo

               !ezwake(i) = ezwake(i) + rhonew(islp)/xxl**(1.0d0/3.0d0)
               !do j = islp, Ni-1
               do j = islp+1, Ni-1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Ni-1) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xx2 = ptmin + (j-1)*hx
                 !csrss = (xx-xx2)**(-1.0d0/3.0d0)
                 ezwake(i) = ezwake(i)+tcoef*&
                            (xx-xx2)**(-1.0d0/3.0d0)*rhonewp(j)*hx
               enddo
               !subgrid integration from Ni-1 to Ni
               hxsup = hx/Nsup
               isup = 1
               jj0 = Ni-isup
               xxp0 = ptmin+(jj0-1)*hx
               Nisup = isup*Nsup
               do j = 1, Nisup+1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Nisup+1) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xx2 = xxp0 + (j-1)*hxsup
                 if(j.ne.Nisup+1) then
                   csrss = (xx-xx2)**(-1.0d0/3.0d0)
                 else
                   csrss = 2*(hxsup)**(-1.0d0/3.0d0)-&
                           (2*hxsup)**(-1.0d0/3.0d0)
                 endif
                 jj = jj0
                 jj1 = jj0+1
                 xxsup = xx2 - (ptmin+(jj-1)*hx)
                 rhoxxp = rhonewp(jj)+(rhonewp(jj1)-rhonewp(jj))*xxsup/hx
                 ezwake(i) = ezwake(i)+tcoef*csrss*rhoxxp*hxsup
               enddo
             endif
             ezwake(i) = ezwake(i)*xk
           else
             xxbar = (xx - blength)/r0
             deltasmax = (r0*phim**3/24.0d0)*(phim+4*xxbar)/(phim+xxbar)
             !print*,"xxbar: ",i,xx,xxbar,deltasmax,xblg
             if(deltasmax.ge.xblg) then !case D
               !neglect the contributions from s to s
               !do j = 1, Ni-1
               do j = 1, Ni
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Ni) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xx2 = ptmin + (j-1)*hx
                 deltas = xx - xx2
                 call psiroot(r0,xxbar,deltas,psi,epstol,Nmax)
                 ezwake(i) = ezwake(i) + tcoef*rhonewp(j)/(psi+2*xxbar)*hx
               enddo
             else if(deltasmax.lt.hx) then !case C
               !print*,"deltasmaxC: ",deltasmax,hx
               deltasmax2 = (r0*phim**2/6.0d0)*(phim+3*xxbar)
               if(deltasmax2.ge.xx-ptmin) then
               else
                 il = (xx-deltasmax2-ptmin)/hx + 1
                 dx = xx-deltasmax2-(ptmin+(il-1)*hx)
                 tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
                 ezwake(i) = ezwake(i) - tmprho/(phim+2*xxbar)
               endif

               deltasmax2 = deltasmax
               il = (xx-deltasmax2-ptmin)/hx + 1
               dx = xx-deltasmax2-(ptmin+(il-1)*hx)
               tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
               ezwake(i) = ezwake(i) + tmprho/(phim+2*xxbar)

               !neglect the integral from xx-deltasmax to xx within hx
             else
               islp = (xx-deltasmax-ptmin)/hx + 1

               !case C from 0 to xblg-deltasmax
               deltasmax2 = (r0*phim**2/6.0d0)*(phim+3*xxbar)
               if(deltasmax2.ge.xx-ptmin) then
               else
                 il = (xx-deltasmax2-ptmin)/hx + 1
                 dx = xx-deltasmax2-(ptmin+(il-1)*hx)
                 tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
                 ezwake(i) = ezwake(i) - tmprho/(phim+2*xxbar)
               endif 

               !case D
               ezwake(i) = ezwake(i) + rhonew(islp)/(phim+2*xxbar)
               !do j = islp, Ni-1
               do j = islp, Ni
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Ni) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xx2 = ptmin + (j-1)*hx
                 deltas = xx - xx2
                 call psiroot(r0,xxbar,deltas,psi,epstol,Nmax)
                 ezwake(i) = ezwake(i) + tcoef*rhonewp(j)/(psi+2*xxbar)*hx
               enddo
             endif
             ezwake(i) = -ezwake(i)*4/r0*xconst
           endif
        enddo

        end subroutine csrwakeTr3_FieldQuant

        !this subroutine calculate the 1d csr wakefield including
        !stead-state and transitions drift effects
        !here the entrance effects are commented out for testing purpose
        !here, hx, rho,...are in real unit. the return ezwake is V/m
        subroutine csrwaketmp_FieldQuant(Nx,r0,ptmin,hx,blength,rhonew,rhonewp,&
                             rhonewpp,ezwake)
        implicit none
        integer, intent(in) :: Nx
        real*8 :: r0,ptmin,hx,blength
        real*8, dimension(Nx) :: rhonew,rhonewp,rhonewpp,ezwake
        real*8 :: xx,xxl,xx2,dx,tmprho,epstol,deltas,deltasmax,xxbar,psi,&
                  phim,pilc,yy,hx2,xconst,xk,psitmp
        integer :: Ni,il,i,Nmax,j,Nsup,islp,Nisup,jj0,jj,jj1,isup
        real*8 :: tcoef,hxsup,csrss,xxp0,xxsup,rhoxxp,xblg,deltasmax2

        epstol = 1.0d-10 !tolerance for root finding
        Nmax = 100 !maximum # of iteration in root finding
        pilc = 2*asin(1.0d0)
        phim = blength/r0
        xconst = 1./(4*pilc*8.854187817e-12)
        xk = -2.0d0/(3*r0**2)**(1.0d0/3.0d0)*xconst

        Nsup = 10

        psitmp = 0.0
        ezwake(1) = 0.0 
        do i = 2, Nx
           xx = ptmin + (i-1)*hx 
           xxl = (xx/r0)**3*r0/24
           Ni = i
           ezwake(i) = 0.0
           xblg = (i-1)*hx

           !print*,"xxl: ",i,xx,blength,xxl,xblg
           if(xx.le.0) then
           else if(xx.le.blength) then
!             if(xxl.ge.xblg) then !S-S
               do j = 1, Ni-1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Ni-1) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xx2 = ptmin + (j-1)*hx
                 !csrss = (xx-xx2)**(-1.0d0/3.0d0)
                 ezwake(i) = ezwake(i)+tcoef*&
                            (xx-xx2)**(-1.0d0/3.0d0)*rhonewp(j)*hx
               enddo

               !subgrid integration from Ni-1 to Ni
               hxsup = hx/Nsup
               isup = 1
               jj0 = Ni-isup
               xxp0 = ptmin+(jj0-1)*hx
 
               Nisup = isup*Nsup
 
               do j = 1, Nisup+1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Nisup+1) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
 
                 xx2 = xxp0 + (j-1)*hxsup
                 if(j.ne.Nisup+1) then
                   csrss = (xx-xx2)**(-1.0d0/3.0d0)
                 else
                   csrss = 2*(hxsup)**(-1.0d0/3.0d0)-&
                           (2*hxsup)**(-1.0d0/3.0d0)
                 endif
                 jj = jj0 
                 jj1 = jj0+1
                 xxsup = xx2 - (ptmin+(jj-1)*hx)
                 rhoxxp = rhonewp(jj)+(rhonewp(jj1)-rhonewp(jj))*xxsup/hx
                 ezwake(i) = ezwake(i)+tcoef*csrss*rhoxxp*hxsup
               enddo

             ezwake(i) = ezwake(i)*xk
           else
             xxbar = (xx - blength)/r0
             deltasmax = (r0*phim**3/24.0d0)*(phim+4*xxbar)/(phim+xxbar)
             !print*,"xxbar: ",i,xx,xxbar,deltasmax,xblg
             if(deltasmax.ge.xblg) then !case D
               !neglect the contributions from s to s
               !do j = 1, Ni-1
               do j = 1, Ni
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Ni) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xx2 = ptmin + (j-1)*hx
                 deltas = xx - xx2
                 call psiroot(r0,xxbar,deltas,psi,epstol,Nmax)
                 ezwake(i) = ezwake(i) + tcoef*rhonewp(j)/(psi+2*xxbar)*hx
               enddo
             else if(deltasmax.lt.hx) then !case C
               !print*,"deltasmaxC: ",deltasmax,hx
               deltasmax2 = (r0*phim**2/6.0d0)*(phim+3*xxbar)
               if(deltasmax2.ge.xx-ptmin) then
               else
                 il = (xx-deltasmax2-ptmin)/hx + 1
                 dx = xx-deltasmax2-(ptmin+(il-1)*hx)
                 tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
                 ezwake(i) = ezwake(i) - tmprho/(phim+2*xxbar)
               endif

               deltasmax2 = deltasmax
               il = (xx-deltasmax2-ptmin)/hx + 1
               dx = xx-deltasmax2-(ptmin+(il-1)*hx)
               tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
               ezwake(i) = ezwake(i) + tmprho/(phim+2*xxbar)

               !neglect the integral from xx-deltasmax to xx within hx
             else
               islp = (xx-deltasmax-ptmin)/hx + 1

               !case C from 0 to xblg-deltasmax
               deltasmax2 = (r0*phim**2/6.0d0)*(phim+3*xxbar)
               if(deltasmax2.ge.xx-ptmin) then
               else
                 il = (xx-deltasmax2-ptmin)/hx + 1
                 dx = xx-deltasmax2-(ptmin+(il-1)*hx)
                 tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
                 ezwake(i) = ezwake(i) - tmprho/(phim+2*xxbar)
               endif 

               !case D
               ezwake(i) = ezwake(i) + rhonew(islp)/(phim+2*xxbar)
               !do j = islp, Ni-1
               do j = islp, Ni
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Ni) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xx2 = ptmin + (j-1)*hx
                 deltas = xx - xx2
                 call psiroot(r0,xxbar,deltas,psi,epstol,Nmax)
                 ezwake(i) = ezwake(i) + tcoef*rhonewp(j)/(psi+2*xxbar)*hx
               enddo
             endif
             ezwake(i) = -ezwake(i)*4/r0*xconst
           endif
        enddo

        end subroutine csrwaketmp_FieldQuant
      end module FieldQuantclass
