!----------------------------------------------------------------
! (c) Copyright, 2006 by the Regents of the University of California.
! Filterclass: numerical filter function for a given noisy density distribution
!                of FUNCTION layer.
! Version: 1.0
! Author: Ji Qiang, LBNL, 12/13/06
! Description: This class defines functions to filter out the sampling
!              noise in a density distribution function.
! Comments:
!         The two-level filter is developed by Ilya Pogorelov.
!----------------------------------------------------------------
        module Filterclass
          use FFTclass
        contains

          SUBROUTINE convlv(data,n,respns,m,tmp1)
          INTEGER m,n,NMAX
          REAL*8 respns(n)
          REAL*8 data(n),tmp1(n)
          real*8 fft1(n),fft2(n)
          INTEGER i,no2,isg

          do 11 i=1,(m-1)/2
            respns(n+1-i)=respns(m+1-i)
11        continue
          do 12 i=(m+3)/2,n-(m-1)/2
            respns(i)=0.0
12        continue

          isg = 1
          fft1 = data
          call realft(fft1,n,isg)

          fft2 = respns
          call realft(fft2,n,isg)

          tmp1(1) = fft1(1)*fft2(1)
          tmp1(2) = fft1(2)*fft2(2)
          do i = 2, n/2
            tmp1(2*i-1) = fft1(2*i-1)*fft2(2*i-1)-fft1(2*i)*fft2(2*i)
            tmp1(2*i) = fft1(2*i-1)*fft2(2*i)+fft1(2*i)*fft2(2*i-1)
          enddo

          isg = -1
          call realft(tmp1,n,isg)

          do i = 1, n
            tmp1(i) = tmp1(i)*2.0d0/n
          enddo

          return
          END subroutine convlv

      !filtering using Savizky-Golay method from NR.
      SUBROUTINE savgol(c,np,nl,nr,ld,m)
      INTEGER ld,m,nl,np,nr,MMAX
      REAL*8 c(np)
      PARAMETER (MMAX=6)
!CU    USES lubksb,ludcmp
      INTEGER imj,ipj,j,k,kk,mm,indx(MMAX+1)
      REAL*8 d,fac,sum,a(MMAX+1,MMAX+1),b(MMAX+1)
      if(np.lt.nl+nr+ &
      1.or.nl.lt.0.or.nr.lt.0.or.ld.gt.m.or.m.gt.MMAX.or.nl+nr.lt.m) &
      pause 'bad args in savgol'
      do 14 ipj=0,2*m
        sum=0.
        if(ipj.eq.0)sum=1.
        do 11 k=1,nr
          sum=sum+float(k)**ipj
11      continue
        do 12 k=1,nl
          sum=sum+float(-k)**ipj
12      continue
        mm=min(ipj,2*m-ipj)
        do 13 imj=-mm,mm,2
          a(1+(ipj+imj)/2,1+(ipj-imj)/2)=sum
13      continue
14    continue
      call ludcmp(a,m+1,MMAX+1,indx,d)
      do 15 j=1,m+1
        b(j)=0.
15    continue
      b(ld+1)=1.
      call lubksb(a,m+1,MMAX+1,indx,b)
      do 16 kk=1,np
        c(kk)=0.
16    continue
      do 18 k=-nl,nr
        sum=b(1)
        fac=1.
        do 17 mm=1,m
          fac=fac*k
          sum=sum+b(mm+1)*fac
17      continue
        kk=mod(np-k,np)+1
        c(kk)=sum
18    continue
      return
      END subroutine savgol


      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL*8 sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END subroutine lubksb

      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      REAL*8 d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0d-20)
      INTEGER i,imax,j,k
      REAL*8 aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END subroutine ludcmp

        !The following is a two-level filter developed by Ilya Pogorelov
        subroutine filterIP_FieldQuant(Nz,recvdensz,&
                   rhoout,rhopout,hzi)
        implicit none 
        integer, intent(in):: Nz !, Ndetail 
        integer:: i, k 
        double precision, intent(in):: hzi
        double precision, dimension(1:Nz), intent(in):: recvdensz 
        double precision, dimension(1:Nz), intent(out):: rhopout,rhoout 
        double precision, dimension(:), allocatable:: uodd, ueven 
        double precision, dimension(-1:Nz+2):: u 
        double precision, dimension(1:Nz):: AA, BB 
        double precision :: hz
        integer :: Npts
!------------------------------------------------------------
        hz = hzi 
        Npts = Nz
!
        u(:) = 0.0d0 
        !u(1:Nz) = recvdensz(1:Nz,1) 
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
!
        rhoout(1:Npts) = AA(1:Npts)
        rhopout(1:Npts) = BB(1:Npts)

        end subroutine filterIP_FieldQuant 

      !filtering the density distribution using a relative power spectral
      !threshold in frequency domain
      subroutine Efourier(zmin,zmax,edata,ndatareal,ncoefreal,coeftol,&
                         edatanew,edatanewp,edatanewpp)
      implicit none
      integer :: ncoefreal,ndatareal
      real*8 :: zmin,zmax,r0,coeftol
      double precision, dimension(ndatareal) :: edata,edatanew,&
                                                edatanewp,edatanewpp
      double precision, dimension(ncoefreal) :: Fcoef,Fcoef2
      double precision, dimension(ndatareal) :: rhopp
      integer :: i,j
      real*8 :: zlen,zhalf,zmid,h,pilc,zz,coefmax,tmpsum,tmpsump,tmpsumpp,tmpf1
      real*8 :: rk,rkzz

      zlen = zmax - zmin
      zhalf = zlen/2.0
      zmid = (zmax + zmin)/2
      h = zlen/(ndatareal-1)
      pilc = 2*asin(1.0d0)
      print*,"The RF data number is: ",ndatareal,zlen,zmid,h

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
      print*,"coefmax: ",coefmax

      print*,"coeftol:",coeftol
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
        edatanew(i) = 0.5*Fcoef(1)
        edatanewp(i) = 0.0
        edatanewpp(i) = 0.0
        do j = 2,ncoefreal
         rkzz = (j-1)*zz*rk
         edatanew(i) = edatanew(i) + Fcoef(j)*cos(rkzz) + &
                  Fcoef2(j)*sin(rkzz)
         edatanewp(i) = edatanewp(i)-(j-1)*rk*(Fcoef(j)*sin(rkzz)-&
                           Fcoef2(j)*cos(rkzz))
         edatanewpp(i) = edatanewpp(i)-((j-1)*rk)**2*(Fcoef(j)*cos(rkzz)+&
                                   Fcoef2(j)*sin(rkzz))
!         print*,"i,j: ",i,j,edatanew(i),edatanewp(i)
        enddo
!        write(8,*)zmin+(i-1)*h,edata(i),edatanew(i),edatanewp(i),&
!                  edatanewpp(i)
      enddo
!      close(8)

      end subroutine Efourier

        end module Filterclass
