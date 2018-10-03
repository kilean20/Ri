************************************************************************
* header              GENM (General GENMAP Routines)                   *
*  All routines common to the various GENMAP programs                  *
************************************************************************
      subroutine putmap(y,fa,fm)
c  Written by Rob Ryne, Spring 1986
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension y(monoms+15)
      dimension fa(monoms),fm(6,6)
      call clear(fa,fm)
      do 10 i=28,monoms
   10 fa(i)=y(i+15)
      do 20 i=1,6
      do 20 j=1,6
   20 fm(i,j)=y(i+6*j)
      call revf(1,fa,fm)
      return
      end
c
************************************************************************
c
      subroutine feval(t,yy,f,ne)
c subroutine called by numerical map integration routines
c  Written by Rob Ryne, Spring 1986
c  Modified by F. Neri and A. Dragt 9/29/89
c  Generic Vector Potential Version by F. Neri, 3/28/90
c  Merged generic version with existing ML/I version June 2004 RDR/FN/PW
c
c The equations are expressed as dy/dt=f(y;t), where:
c t denotes the independent variable (usually s, the arc length),
c y denotes the vector representation of the transfer map (elements 7 - 938)
c along with the reference trajectory (elements 1-6)
c This subroutine calls subroutine hlmtnN, which should return the
c coefficients of the Hamiltonian, which will depend on the particular
c beamline element (vector potential).
c This routine returns the RHS of the equation dy/dt=f, i.e. it returns f.
c 
cryneneriwalstrom===
      use lieaparam, only : monoms
      use beamdata, only : brho,gamma,beta
cryneneriwalstrom===
c
      include 'impli.inc' 
c
cryne fix this parameter statement later to use monoms instead of hardwire:
      parameter (itop=923,iplus=938)
      include 'hmflag.inc'
cryneneriwalstrom===
cryneneriwalstrom:param.inc replaced w/ "use beamdata..."
cryneneriwalstrom:parm.inc replaced w/ "use lieaparam..."(needed by vecpot.inc)
      include 'bfield.inc'
      include 'gronax.inc'
      include 'vecpot.inc'
cryneneriwalstrom===
      dimension h(itop),hint(itop),temp(itop),t5a(itop),t5b(itop)
      dimension yy(iplus),f(iplus)
      dimension ajsm(6,6)
cryne 7/11/2002 this missing "save" cost me many hours of lost sleep!
cryne 7/11/2002 would be nice to deal w/ the issue of save statements
cryne 7/11/2002 once and for all.
      save ajsm,h
cryneneriwalstrom: check to make sure there are not other variables that
cryneneriwalstrom  need to be saved in the merged version
c
c  y(1-6) = given (design) trajectory
c  y(7-42) = matrix
c  y(43-98) = f3
c  y(99-224) = f4
c  y(225-476) = f5
c  y(477-938) = f6
c
c  design trajectory:
c  NB We have assumed the coordinates have been selected so that
c  the design trajectory is yy(i)=0.d0 for i=1 to 6.
      f(1)=0.d0
      f(2)=0.d0
      f(3)=0.d0
      f(4)=0.d0
      f(5)=0.d0
      f(6)=0.d0
c
c  matrix and polynomials:
c  the expressions for f below are from p. 2742 of
c  J. Math.Phys.,Vol 24,No.12,Dec 1983 (A.J.Dragt and E.Forest)
c
c  select and compute hamiltonian, s matrix, and ajsm=-j*s
c  skip these steps if iflag=0
      if (iflag.eq.0) go to 100
c  otherwise carry out the selection and required computations
cryneneriwalstrom      go to (10,20),iflag
      go to (10,20,30,40,50),iflag
      write(6,*) 'trouble with iflag in subroutine feval'
      call myexit
   10 call hmltn1(h)
      call matify(ajsm,h)
      goto 100
   20 call hmltn2(t,yy,h)
      call matify(ajsm,h)
      go to 100
cryneneriwalstrom eventually, need to check that this works for
cryneneriwalstrom arbitrary choice of units.
c
cryneneriwalstrom===== code from Neri generic version follows:
   30 continue
      call hmltn3(t,yy,h)
      call matify(ajsm,h)
c Change in design trajectory in the general case: no midplane symmetry
      q = 1./brho
      x = yy(1)
      y = yy(3)
c      t = yy(5) ! No time dependence; Use other symbol otherwise (t is taken!)
      Px = yy(2)
      Py = yy(4)
      Pt = yy(6)
c
      k1 = 1
      k2 = 2
      k3 = 3
      k4 = 4
      k5 = 5
      k6 = 6
c
c Hamilton equations , given vector potentials ( Ax, Ay, Az ) :
c  Ax, Ay, Az come thru common from hmltn3
c
      Pix = yy(2) - Ax(0)
      Piy = yy(4) - Ay(0)
      xm2 = 1.d0/(gamma*beta)**2
c
      root = Dsqrt(-xm2 + yy(6)**2 - Pix**2 - Piy**2)
c
      f(k1) = Pix/root
      f(k2) = (Ax(1)*Pix+Ay(1)*Piy)/root + Az(1)
      f(k3) = Piy/root
      f(k4) = (Ax(3)*Pix+Ay(3)*Piy)/root + Az(3)
      f(k5) = (-yy(6))/root
      f(k6) = 0.0d0
c      write(6,*) t, Pix, yy(1)
c      write(36,*) t, Pix, yy(1)
      goto 100
c
c  Insert your code here!
c
   40 call hmltn4(t,yy,h)
      call matify(ajsm,h)
      go to 100
cfneri
cfneri The pbend rutines are here
cfneri
cfneri hmltn5 routine from gendip5:
   50 call hmltn5(t,yy,h)
      call matify(ajsm,h)
c
c Hamilton equations , given vector potentials ( Ax, Ay, Az ) :
c Ax, Ay, Az come thru common from hmltn5.
c
      Pix = yy(2) - Ax(0)
      Piy = yy(4) - Ay(0)
      xm2 = 1.d0/(gamma*beta)**2
c
      root = Dsqrt(-xm2 + yy(6)**2 - Pix**2 - Piy**2)
c
      f(1) = Pix/root
      f(2) = (Ax(1)*Pix+Ay(1)*Piy)/root + Az(1)
      f(3) = Piy/root
      f(4) = (Ax(3)*Pix+Ay(3)*Piy)/root + Az(3)
      f(5) = (-yy(6))/root
      f(6) = 0.0d0
cfneri
c      write(6,*) 'FEVAL: z, Pix, x', t, Pix, yy(1)
c      write(36,*)'FEVAL: z, Pix, x', t, Pix, yy(1)
cfneri
      go to 100
cryneneriwalstrom===== code from Neri generic version above
  100 continue
c
c  continue to evaluate f's
c
c  matrix part:  dm/dt = j*s*m = -ajsm*m
      do 110 i=7,42
  110 f(i)=0.d0
      do 120 i=1,6
      do 120 j=1,6
      do 120 k=1,6
c compute  xm2dot(i,j)=-ajsm(i,k)*xm2(k,j)
      ij=i+6*j
  120 f(ij)=f(ij) - ajsm(i,k)*yy(k+6*j)
      if(ne.eq.42)return
c
c  compute f3dot **********************************
      call xform5(h,3,yy(7),hint)
      do 130 i=43,98
  130 f(i)=-hint(i-15)
      if(ne.eq.98)return
c  compute f4dot **********************************
      call xform5(h,4,yy(7),hint)
      call pbkt1(yy(16),3,f(16),3,temp)
      do 140 i=99,224
  140 f(i)=-hint(i-15)+0.5*temp(i-15)
      if(ne.eq.224)return
c  compute f5dot **********************************
      call xform5(h,5,yy(7),hint)
      call pbkt1(yy(16),3,temp,4,t5a)
      call pbkt1(yy(16),3,f(16),4,t5b)
      do 150 i=225,476
  150 f(i)=-hint(i-15)-t5a(i-15)/6.d0+t5b(i-15)
      if(ne.eq.476)return
c  compute f6dot **********************************
      call xform5(h,6,yy(7),hint)
      do 160 i=210,461
  160 temp(i)=t5a(i)/24.d0-t5b(i)/2.d0+f(15+i)
      call pbkt1(yy(16),3,temp,5,t5a)
      call pbkt1(yy(16),4,f(16),4,t5b)
      do 170 i=477,938
  170 f(i)=-hint(i-15)+t5a(i-15)+t5b(i-15)/2.d0
cfneri
c      write(6,*) ' End of FEVAL!'
cfneri
      return
      end
c
************************************************************************
c
      subroutine adam11(h,ns,nf,t,y,ne)
c  Written by Rob Ryne, Spring 1986, based on a routine of Alex Dragt
c  This integration routine makes local truncation errors at each
c  step of order h**11.  That is, it is locally correct through
c  order h**10.  Due to round off errors, its true precision is
c  realized only when more than 64 bits are used.
      use lieaparam, only : monoms
      include 'impli.inc'
      character*6 nf
cryneneriwalstrom fix later to use monoms instead of hardwire:
      dimension y(938),yp(938),yc(938),f1(938),f2(938),f3(938),f4(938),
     & f5(938),f6(938),f7(938),f8(938),f9(938),f10(938),f11(938)
      dimension a(10),am(10),b(10),bm(10)
c
      data (a(i),i=1,10)/57281.d0,-583435.d0,2687864.d0,
     & -7394032.d0,13510082.d0,-17283646.d0,16002320.d0,
     & -11271304.d0,9449717.d0,2082753.d0/
      data (b(i),i=1,10)/-2082753.d0,20884811.d0,-94307320.d0,
     & 252618224.d0,-444772162.d0,538363838.d0,-454661776.d0,
     & 265932680.d0,-104995189.d0,30277247.d0/
cryne 7/23/2002
      save a,b
cryne 1 August 2004      ne=monoms+15
c
      nsa=ns
      if (nf.eq.'cont') go to 20
c  rk start
      iqt=5
      qt=float(iqt)
      hqt=h/qt
      call feval(t,y,f1,ne)
      call rk78ii(hqt,iqt,t,y,ne)
      call feval(t,y,f2,ne)
      call rk78ii(hqt,iqt,t,y,ne)
      call feval(t,y,f3,ne)
      call rk78ii(hqt,iqt,t,y,ne)
      call feval(t,y,f4,ne)
      call rk78ii(hqt,iqt,t,y,ne)
      call feval(t,y,f5,ne)
      call rk78ii(hqt,iqt,t,y,ne)
      call feval(t,y,f6,ne)
      call rk78ii(hqt,iqt,t,y,ne)
      call feval(t,y,f7,ne)
      call rk78ii(hqt,iqt,t,y,ne)
      call feval(t,y,f8,ne)
      call rk78ii(hqt,iqt,t,y,ne)
      call feval(t,y,f9,ne)
      call rk78ii(hqt,iqt,t,y,ne)
      call feval(t,y,f10,ne)
      nsa=ns-9
      hdiv=h/7257600.0d+00
      do 10 i=1,10
      am(i)=hdiv*a(i)
  10  bm(i)=hdiv*b(i)
  20  tint=t
      do 100 i=1,nsa
      do 30 j=1,ne
      yp(j)=y(j)+bm(1)*f1(j)+bm(2)*f2(j)+bm(3)*f3(j)
     &+bm(4)*f4(j)+bm(5)*f5(j)+bm(6)*f6(j)+bm(7)*f7(j)
     & +bm(8)*f8(j)+bm(9)*f9(j)+bm(10)*f10(j)
   30 continue
      call feval(t+h,yp,f11,ne)
      do 40 j=1,ne
      yp(j)=y(j)+am(1)*f2(j)+am(2)*f3(j)+am(3)*f4(j)+am(4)*f5(j)
     & +am(5)*f6(j)+am(6)*f7(j)+am(7)*f8(j)+am(8)*f9(j)+am(9)*f10(j)
  40  yc(j)=yp(j)+am(10)*f11(j)
  41  call feval(t+h,yc,f11,ne)
      do 50 j=1,ne
   50 y(j)=yp(j)+am(10)*f11(j)
      do 60 j=1,ne
      f1(j)=f2(j)
      f2(j)=f3(j)
      f3(j)=f4(j)
      f4(j)=f5(j)
      f5(j)=f6(j)
      f6(j)=f7(j)
      f7(j)=f8(j)
      f8(j)=f9(j)
      f9(j)=f10(j)
   60 f10(j)=f11(j)
      t=tint+i*h
  100 continue
      return
      end
c
***********************************************************************
c
      subroutine rk78ii(h,ns,t,y,ne)
c  Written by Rob Ryne, Spring 1986, based on a routine of
c  J. Milutinovic.
c  For a reference, see page 76 of F. Ceschino and J Kuntzmann,
c  Numerical Solution of Initial Value Problems, Prentice Hall 1966.
c  This integration routine makes local truncation errors at each
c  step of order h**7.
c  That is, it is locally correct through terms of order h**6.
c  Each step requires 8 function evaluations.

      use lieaparam, only : monoms
      include 'impli.inc'
cryneneriwalstrom fix later to use monoms instead of hardwire:
      dimension y(938),yt(938),f(938),a(938),b(938),c(938),d(938),
     &e(938),g(938),o(938),p(938)
cryne 1 August 2004      ne=monoms+15
c
      tint=t
      do 200 i=1,ns
      call feval(t,y,f,ne)
      do 10 j=1,ne
   10 a(j)=h*f(j)
      do 20 j=1,ne
   20 yt(j)=y(j)+a(j)/9.d+0
      tt=t+h/9.d+0
      call feval(tt,yt,f,ne)
      do 30 j=1,ne
   30 b(j)=h*f(j)
      do 40 j=1,ne
   40 yt(j)=y(j) + (a(j) + 3.d+0*b(j))/24.d+0
      tt=t+h/6.d+0
      call feval(tt,yt,f,ne)
      do 50 j=1,ne
   50 c(j)=h*f(j)
      do 60 j=1,ne
  60  yt(j)=y(j)+(a(j)-3.d+0*b(j)+4.d+0*c(j))/6.d+0
      tt=t+h/3.d+0
      call feval(tt,yt,f,ne)
      do 70 j=1,ne
  70  d(j)=h*f(j)
      do 80 j=1,ne
   80 yt(j)=y(j) + (-5.d+0*a(j) + 27.d+0*b(j) -
     &  24.d+0*c(j) + 6.d+0*d(j))/8.d+0
      tt=t+.5d+0*h
      call feval(tt,yt,f,ne)
      do 90 j=1,ne
   90 e(j)=h*f(j)
      do 100 j=1,ne
  100 yt(j)=y(j) + (221.d+0*a(j) - 981.d+0*b(j) +
     &  867.d+0*c(j)- 102.d+0*d(j) + e(j))/9.d+0
      tt = t+2.d+0*h/3.d+0
      call feval(tt,yt,f,ne)
      do 110 j=1,ne
  110 g(j)=h*f(j)
      do 120 j=1,ne
  120 yt(j) = y(j)+(-183.d+0*a(j)+678.d+0*b(j)-472.d+0*c(j)-
     &  66.d+0*d(j)+80.d+0*e(j) + 3.d+0*g(j))/48.d+0
      tt = t + 5.d+0*h/6.d+0
      call feval(tt,yt,f,ne)
      do 130 j=1,ne
  130 o(j)=h*f(j)
      do 140 j=1,ne
  140 yt(j) = y(j)+(716.d+0*a(j)-2079.d+0*b(j)+1002.d+0*c(j)+
     & 834.d+0*d(j)-454.d+0*e(j)-9.d+0*g(j)+72.d+0*o(j))/82.d+0
      tt = t + h
      call feval(tt,yt,f,ne)
      do 150 j=1,ne
  150 p(j)=h*f(j)
      do 160 j=1,ne
  160 y(j) = y(j)+(41.d+0*a(j)+216.d+0*c(j)+27.d+0*d(j)+
     &  272.d+0*e(j)+27.d+0*g(j)+216.d+0*o(j)+41.d+0*p(j))/840.d+0
      t=tint+i*h
  200 continue
      return
      end
cryneneriwalstrom===
c removed subroutine poly1 and pmult, which are now in integ.f 
c also, left in place below a commented out version of subroutine drift
c as a reminder to ourselves
cryneneriwalstrom===
c ******************************************************************
c
c  Drift and aux poly routines.
c  They don't really belong here.
c
c *****************************************************************
c
c      subroutine drift(l,h,mh)
c
c     generates linear matrix mh and
c     array h containing nonlinearities
c     for the transfer map describing
c     a drift section of length l metres
c
c      implicit double precision (a-h,o-z)
c      double precision l,lsc,mh
c      dimension j(6)
c      dimension h(923)
c      dimension mh(6,6)
c      dimension X(0:923), Y(0:923), A(0:6)
c      common/parm/brho,c,gamma,gamm1,beta,achg,sl,ts
c
c      call clear(h,mh)
c      lsc=l/sl
c
c     add drift terms to mh
c
c      do 40 k=1,6
c      mh(k,k)=+1.0d0
c   40 continue
c      mh(1,2)=+lsc
c      mh(3,4)=+lsc
c      mh(5,6)=+(lsc/((gamma**2)*(beta**2)))
c
c     add drift terms to h
c
c     do degree 3-6
c
c   F. Neri    Aug. 28 1987
c
c      A(0) = 0.d0
c      A(1) = 1.d0/2.d0
c      do 50 i=2, 6
c        A(i) =  A(i-1) * (1.d0/2.d0 - (i-1.d0))/i
c   50 continue
c 
c      do 60 i=0,923
c        X(i) = 0.d0
c   60 continue
c      X(6) = -2.d0 / beta
c      X(13) = -1.d0
c      X(22) = -1.d0
c      X(27) =  1.d0
c 
c      maxord = 6
c      nn = 6
c      call poly1(nn,A,X,Y,maxord)
c      do 70 i=28, 923
c        h(i) = Y(i) * lsc
c   70 continue
c      return
c      end
c *******************************************************************
