!----------------------------------------------------------------
! (c) Copyright, 2001 by the Regents of the University of California.
! Besselclass: Bessel function class in Math Function module of 
!              FUNCTION layer.
! Version: 1.0
! Author: Ji Qiang, LANL, 3/7/01
! Description: This class finds the coefficients of Bessel expansion
!              and its inverse expansion. It also calculates the roots 
!              and norms of Bessel function.
! Comments:
!----------------------------------------------------------------
      module Besselclass

      contains
        !Given the discrete input function "rho", find its Bessel
        !function expansion coefficients along radial "r" and store
        !back in "rho".
        subroutine Bessr_Bessel(nx,nsizey,nsizez,rho,bess,&
                                bessnorm,hx,modth,myidy,nmod)
        implicit none
        include 'mpif.h'
        integer,intent(in) :: nx,nsizez,nsizey,myidy,nmod
        double precision, dimension(nx,nsizey,nsizez), intent(inout) &
                            :: rho
        double precision, dimension(nx,nx,nmod), intent(in) :: bess
        double precision, dimension(nx,nmod), intent(in) :: bessnorm
        double precision, intent(in) :: hx
        integer, dimension(nsizey), intent(in) :: modth
        double precision, dimension(nx) :: right
        double precision :: simp
        integer :: i,j,k,ii,jj,jm
        
        do k = 1, nsizez
          do j = 1, nsizey
            !find the azmuthal mode number on each local PE.
            if(myidy.ne.0) then
              jm = modth(j)-modth(1) + 1
            else
!              jm = modth(j)-modth(1) + 1
!              if(j.eq.2) jm = j
              if(j.le.2) then
                 jm = j
              else
                jm = modth(j)-modth(1) + 2
              endif
            endif
            do i = 1, nx
              right(i) = rho(i,j,k)
              rho(i,j,k) = 0.0
            enddo
            !Compute the Bessel expansion coeficients from integral.
            !"2 and nx-1" is due to the integrand vanish at r=0 and
            !r = a, ie. boundary. 
            !Trapezoid rule.
            !do jj = 2, nx-1
            !do jj = 1, nx
            !  do ii = 1, nx
            !    rho(ii,j,k) = rho(ii,j,k)+bess(ii,jj,jm)*right(jj)*&
            !                  (jj-1)*hx*hx/bessnorm(ii,jm)
            !  enddo
            !enddo
            !Simpson's rule integration.
            do jj = 1, nx
              if(jj.eq.1) then
                simp = 1.0d0/3.0d0
              else if(jj.eq.nx) then
                simp = 1.0d0/3.0d0
              else if(mod(jj,2).eq.0) then
                simp = 4.0d0/3.0d0
              else
                simp = 2.0d0/3.0d0
              endif
              do ii = 1, nx
                rho(ii,j,k) = rho(ii,j,k)+simp*bess(ii,jj,jm)*right(jj)*&
                              (jj-1)*hx*hx/bessnorm(ii,jm)
              enddo
            enddo

            !do jj = 1, nx
            !  print*,"rhocoef: ",jj,j,k,rho(jj,j,k)
            !enddo

          enddo
        enddo

        return
        end subroutine Bessr_Bessel

        !Given the Bessel expansion coefficients stored in "rho, find
        !the 3d discrete function "rho(r,theta,z)".
        subroutine BessrI_Bessel(nx,nsizey,nsizez,rho,bess,modth,&
                                     myidy,nmod)
        implicit none
        include 'mpif.h'
        integer,intent(in) :: nx,nsizez,nsizey,myidy,nmod
        double precision, dimension(nx,nsizey,nsizez), intent(inout) &
                            :: rho
        double precision, dimension(nx,nx,nmod), intent(in) :: bess
        integer, dimension(nsizey), intent(in) :: modth
        double precision, dimension(nx) :: right
        integer :: i,j,k,ii,jj,jm
        
        do k = 1, nsizez
          do j = 1, nsizey
            !find the azmuthal mode number on each local PE.
            if(myidy.ne.0) then
              jm = modth(j)-modth(1) + 1
            else
!              jm = modth(j)-modth(1) + 1
!              if(j.eq.2) jm = j
              if(j.le.2) then
                jm = j
              else
                jm = modth(j)-modth(1) + 2
              endif
            endif
            do i = 1, nx
              right(i) = rho(i,j,k)
            enddo
            !rho(1,j,k) = right(1)
            !Compute function value from its Bessel expansion coeficients.
            do jj = 1, nx-1
              rho(jj,j,k) = 0.0
              !do ii = 1, nx
              !using a cut-off Bessel expansion to avoid the noise from higher order mode.
              !do ii = 1, 20
              do ii = 1, nx/5
                rho(jj,j,k) = rho(jj,j,k)+bess(ii,jj,jm)*right(ii)
                !print*,"ii,jj: ",ii,jj,bess(ii,jj,jm),right(ii),rho(jj,j,k)
              enddo
            enddo
            rho(nx,j,k) = 0.0
          enddo
        enddo

        return
        end subroutine BessrI_Bessel

      !preparation for finding the first derivative of the function using
      !the Bessel expansion.
      subroutine Bessprep2_Bessel(hr,innx,inny,nmod,modth,besscoef,&
                          besscoefP,bessnorm,gml)
      implicit none
      double precision, intent(in) :: hr
      integer, intent(in) :: innx,inny,nmod
      integer, dimension(inny), intent(in) :: modth
      double precision, dimension(innx,innx,nmod), intent(out) :: besscoef,&
                                                   besscoefP
      double precision, dimension(innx,nmod), intent(out) :: bessnorm,gml
      double precision :: hx,f1,f2,rad,tmp1
      real*8, dimension(2,innx) :: range
      real*8 :: x,zeros,rr,tol
      integer iroot,nroot,i,norder,imod,j,jj
      !double precision :: bessj, zbrent

      rad = (innx-1)*hr

      nroot = innx
      hx = 0.25
      norder = -10
      imod = 0
      do i = 1, inny
        if(norder.ne.modth(i)) then
          imod  = imod + 1
          norder = modth(i)
          iroot = 0
          x = 0.01
          do
            f1 = bessj(norder,x)
            x = x + hx
            f2 = bessj(norder,x)
            if(f1*f2.lt.0.0) then
              iroot = iroot + 1
              range(1,iroot) = x - hx
              range(2,iroot) = x
            endif
            if(iroot.ge.nroot) exit
          enddo
      
          !tol = 1.0e-4
          !tol = 1.0e-6
          tol = 1.0e-8
          do j = 1, nroot
            !gml(j,imod) = zbrent(norder,range(1,j),range(2,j),tol)
            gml(j,imod) = zbrent(norder,range(1,j),range(2,j),tol)/rad
            !print*,"roots: ",j,imod,gml(j,imod)
            zeros = rad*gml(j,imod)
            !calculate J_m(gamma_lm r)
            if(norder.ne.0) then
              tmp1 = (bessj(norder-1,zeros)-bessj(norder+1,zeros))/2
              besscoef(j,1,imod) = 0.0
            else
              tmp1 = -bessj(1,zeros)
              besscoef(j,1,imod) = 1.0
            endif
            do jj = 2, nroot-1
              rr = (jj-1)*hr*gml(j,imod)
              !calculate J_m(gamma_lm r)
              besscoef(j,jj,imod) = bessj(norder,rr)
            enddo
            !boundary condition
            besscoef(j,nroot,imod) = 0.0
            bessnorm(j,imod) = rad*rad*tmp1*tmp1/2
            !calculate J'_m(gamma_lm r)
            if(norder.ne.1) then
              besscoefP(j,1,imod) = 0.0
            else
              besscoefP(j,1,imod) = 0.5
            endif
            !print*,"bes1: ",j,1,imod,besscoefP(j,1,imod)
            do jj = 2, nroot
              rr = (jj-1)*hr*gml(j,imod)
              !calculate J'_m(gamma_lm r)
              if(norder.ne.0) then
                besscoefP(j,jj,imod) = gml(j,imod)*&
                                    (bessj(norder-1,rr)-bessj(norder+1,rr))/2
              else
                besscoefP(j,jj,imod) = -gml(j,imod)*&
                                    bessj(norder+1,rr)
              endif
              !print*,"bes2: ",j,jj,imod,norder,besscoefP(j,jj,imod),gml(j,imod)
            enddo
          enddo
        else
        endif
      enddo

      return

      end subroutine BessPrep2_Bessel

      subroutine Bessprep_Bessel(hr,innx,inny,nmod,modth,besscoef,&
                          bessnorm,gml)
      implicit none
      double precision, intent(in) :: hr
      integer, intent(in) :: innx,inny,nmod
      integer, dimension(inny), intent(in) :: modth
      double precision, dimension(innx,innx,nmod), intent(out) :: besscoef
      double precision, dimension(innx,nmod), intent(out) :: bessnorm,gml
      double precision :: hx,f1,f2,rad,tmp1
      real*8, dimension(2,innx) :: range
      real*8 :: x,zeros,rr,tol
      integer iroot,nroot,i,norder,imod,j,jj
      !double precision :: bessj, zbrent

      rad = (innx-1)*hr

      nroot = innx
      hx = 0.25
      norder = -10
      imod = 0
      do i = 1, inny
        if(norder.ne.modth(i)) then
          imod  = imod + 1
          norder = modth(i)
          iroot = 0
          x = 0.01
          do
            f1 = bessj(norder,x)
            x = x + hx
            f2 = bessj(norder,x)
            if(f1*f2.lt.0.0) then
              iroot = iroot + 1
              range(1,iroot) = x - hx
              range(2,iroot) = x
            endif
            if(iroot.ge.nroot) exit
          enddo
      
          !tol = 1.0e-4
          !tol = 1.0e-6
          tol = 1.0e-8
          do j = 1, nroot
            !gml(j,imod) = zbrent(norder,range(1,j),range(2,j),tol)
            gml(j,imod) = zbrent(norder,range(1,j),range(2,j),tol)/rad
            !print*,"roots: ",j,imod,gml(j,imod)
            zeros = rad*gml(j,imod)
            if(norder.ne.0) then
              tmp1 = (bessj(norder-1,zeros)-bessj(norder+1,zeros))/2
              !tmp1 = bessj(norder+1,zeros)
              besscoef(j,1,imod) = 0.0
            else
              tmp1 = -bessj(1,zeros)
              besscoef(j,1,imod) = 1.0
            endif
            bessnorm(j,imod) = rad*rad*tmp1*tmp1/2
            do jj = 2, nroot-1
              rr = (jj-1)*hr*gml(j,imod)
              besscoef(j,jj,imod) = bessj(norder,rr)
            enddo
            besscoef(j,nroot,imod) = 0.0
          enddo
        else
        endif
      enddo

      return

      end subroutine BessPrep_Bessel
     
      FUNCTION zbrent(nord,x1,x2,tol)
      INTEGER ITMAX, nord
      REAL*8 zbrent,tol,x1,x2,EPS
      PARAMETER (ITMAX=100,EPS=3.e-8)
      INTEGER iter
      REAL*8 a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      a=x1
      b=x2
      fa=bessj(nord,a)
      fb=bessj(nord,b)
      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))pause &
      'root must be bracketed for zbrent'
      c=b
      fc=fb
      do 11 iter=1,ITMAX
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.)then
          zbrent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q
          p=abs(p)
          if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=bessj(nord,b)
11    continue
      pause 'zbrent exceeding maximum iterations'
      zbrent=b
      return
      END function zbrent

      FUNCTION bessj(n,x)
      INTEGER n,IACC
      REAL*8 bessj,x,BIGNO,BIGNI
      PARAMETER (IACC=40,BIGNO=1.e10,BIGNI=1.e-10)
      INTEGER j,jsum,m
      REAL*8 ax,bj,bjm,bjp,sum,tox
      if(n.eq.0) then
        bessj = bessj0(x)
        return
      else if(n.eq.1) then
        bessj = bessj1(x)
        return
      else

      ax=abs(x)
      if(ax.eq.0.)then
        bessj=0.
      else if(ax.gt.float(n))then
        tox=2./ax
        bjm=bessj0(ax)
        bj=bessj1(ax)
        do 11 j=1,n-1
          bjp=j*tox*bj-bjm
          bjm=bj
          bj=bjp
11      continue
        bessj=bj
      else
        tox=2./ax
        m=2*((n+int(sqrt(float(IACC*n))))/2)
        bessj=0.
        jsum=0
        sum=0.
        bjp=0.
        bj=1.
        do 12 j=m,1,-1
          bjm=j*tox*bj-bjp
          bjp=bj
          bj=bjm
          if(abs(bj).gt.BIGNO)then
            bj=bj*BIGNI
            bjp=bjp*BIGNI
            bessj=bessj*BIGNI
            sum=sum*BIGNI
          endif
          if(jsum.ne.0)sum=sum+bj
          jsum=1-jsum
          if(j.eq.n)bessj=bjp
12      continue
        sum=2.*sum-bj
        bessj=bessj/sum
      endif
      if(x.lt.0..and.mod(n,2).eq.1)bessj=-bessj
      return

      endif

      END function bessj

      FUNCTION bessj1(x)
      REAL*8 bessj1,x
      REAL*8 ax,xx,z
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,&
      s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,&
      s5,s6
      DATA r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0, &
      242396853.1d0,-2972611.439d0,15704.48260d0,-30.16036606d0/,s1,s2,&
      s3,s4,s5,s6/144725228442.d0,2300535178.d0,18583304.74d0,&
      99447.43394d0,376.9991397d0,1.d0/
      DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4, &
      .2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0, &
      -.2002690873d-3,.8449199096d-5,-.88228987d-6,.105787412d-6/
      if(abs(x).lt.8.)then
        y=x**2
        bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+&
        y*(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8./ax
        y=z**2
        xx=ax-2.356194491
        bessj1=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y* &
      p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*sign(1.0d0,x)
      endif
      return
      END function bessj1

      FUNCTION bessj0(x)
      REAL*8 bessj0,x
      REAL*8 ax,xx,z
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,&
      s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,&
      s5,s6
      DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4, &
      -.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1,&
      .1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0, &
      651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/,s1,s2,&
      s3,s4,s5,s6/57568490411.d0,1029532985.d0,9494680.718d0, &
      59272.64853d0,267.8532712d0,1.d0/
      if(abs(x).lt.8.)then
        y=x**2
        bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*&
      (s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8./ax
        y=z**2
        xx=ax-.785398164
        bessj0=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*&
      p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      END function bessj0

      end module Besselclass
