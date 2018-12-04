!----------------------------------------------------------------
! (c) Copyright, 2001 by the Regents of the University of California.
! NonlinearLensclass: Nonlinear lens beam line element class
!             in Lattice module of APPLICATION layer.
! Version: 1.0
! Author: Ji Qiang, Robert Ryne, LANL, 7/13/01
! Description: This class defines the linear transfer map and field
!              for the nonlinear lens beam line elment.
! Comments:
!----------------------------------------------------------------
      module NonlinearLensclass
        use PhysConstclass
        use Dataclass
        use Utilityclass
!        use Multipoleclass, only: NonlinearLensPropagatorCmplx, &
!             NonlinearLensPropagator
        use Multipoleclass
        integer, private, parameter :: Nparam = 10
        type NonlinearLens
          !Itype = 6
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
          ! Param(1) : zedge
          !      (2) : tn dimensionless strength of NL lens 
          !      (3) : cn dimensional parameter of NL lens
          !      (4) : phase advance across NL lens
          !      (5) : radius
          !      (6) : x misalignment error
          !      (7) : y misalignment error
          !      (8) : rotation error x
          !      (9) : rotation error y
          !      (10) : rotation error z
        end type NonlinearLens
        interface getparam_NonlinearLens
          module procedure getparam1_NonlinearLens,  &
                          getparam2_NonlinearLens,   &
                          getparam3_NonlinearLens
        end interface
        interface setparam_NonlinearLens
          module procedure setparam1_NonlinearLens,  &
                          setparam2_NonlinearLens, setparam3_NonlinearLens
        end interface
      contains
        subroutine construct_NonlinearLens(this,numseg,nmpstp,type,blength)
        implicit none
        type (NonlinearLens), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_NonlinearLens
   
        subroutine setparam1_NonlinearLens(this,i,value)
        implicit none
        type (NonlinearLens), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_NonlinearLens

        subroutine setparam2_NonlinearLens(this,values)
        implicit none
        type (NonlinearLens), intent(inout) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param = values

        end subroutine setparam2_NonlinearLens

        subroutine setparam3_NonlinearLens(this,numseg,nmpstp,type,blength)
        implicit none
        type (NonlinearLens), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_NonlinearLens
   
        subroutine getparam1_NonlinearLens(this,i,blparam) 
        implicit none 
        type (NonlinearLens), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_NonlinearLens
  
        subroutine getparam2_NonlinearLens(this,blparams)
        implicit none
        type (NonlinearLens), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams = this%Param

        end subroutine getparam2_NonlinearLens

        subroutine getparam3_NonlinearLens(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (NonlinearLens), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_NonlinearLens
       
        subroutine maplinear_NonlinearLens(t,tau,xm,this,refpt,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,tau,Bchg,Bmass
        double precision, dimension(6,6), intent(out) :: xm
        double precision, dimension(6), intent(inout) :: refpt
        type (NonlinearLens), intent(in) :: this
        double precision, dimension(18) :: y
        double precision :: h,t0,gi,betai,ui,squi,sqi3
        double precision :: dlti,thli,gf,betaf,uf,squf,sqf3
        double precision :: dltf,thlf,zedge,tfin
        integer :: mpstp

        xm = 0.0

        zedge = this%Param(1)
        mpstp = this%Mapstp

        y(1)=0.0
        y(2)=0.0
        y(3)=0.0
        y(4)=0.0
        y(5)=refpt(5)
        y(6)=refpt(6)
        y(7)=1.0
        y(8)=0.0
        y(9)=0.0
        y(10)=1.0
        y(11)=1.0
        y(12)=0.0
        y(13)=0.0
        y(14)=1.0
        y(15)=1.0
        y(16)=0.0
        y(17)=0.0
        y(18)=1.0
        h=tau/mpstp
        t0=t
        gi=-refpt(6)
        betai=sqrt(gi**2-1.)/gi
        ui=gi*betai
        squi=sqrt(ui)
        sqi3=sqrt(ui)**3

        dlti=0.0
        thli=0.0

        call rk6i_NonlinearLens(h,mpstp,t0,y,18,this,Bchg,Bmass)

        refpt(5)=y(5)
        refpt(6)=y(6)
        gf=-refpt(6)
        betaf=sqrt(gf**2-1.)/gf
        uf=gf*betaf
        squf=sqrt(uf)
        sqf3=sqrt(uf)**3

        tfin = t + tau

        dltf = 0.0
        thlf = 0.0

        xm(1,1)= y(7)*squi/squf+y(9)*squi/squf*dlti
        xm(2,1)=(y(8)-y(7)*dltf)*squi*squf+(y(10)-y(9)*dltf)*squi*squf*&
                 dlti
        xm(1,2)= y(9)/(squi*squf)
        xm(2,2)=(y(10)-y(9)*dltf)*squf/squi
        xm(3,3)= y(11)*squi/squf+y(13)*squi/squf*dlti
        xm(4,3)= &
        (y(12)-y(11)*dltf)*squi*squf+(y(14)-y(13)*dltf)*squi*squf*dlti
        xm(3,4)= y(13)/(squi*squf)
        xm(4,4)=(y(14)-y(13)*dltf)*squf/squi
        xm(5,5)= y(15)*sqi3/sqf3+y(17)*sqi3/sqf3*thli
        xm(6,5)= &
        (y(16)-y(15)*thlf)*sqi3*sqf3+(y(18)-y(17)*thlf)*sqi3*sqf3*thli
        xm(5,6)= y(17)/(sqi3*sqf3)
        xm(6,6)=(y(18)-y(17)*thlf)*sqf3/sqi3

        end subroutine maplinear_NonlinearLens

        subroutine rk6i_NonlinearLens(h,ns,t,y,nvar,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ns,nvar
        double precision, intent(inout) :: t
        double precision, intent(in) :: h,Bchg,Bmass
        double precision, dimension(nvar), intent(inout) :: y
        type (NonlinearLens), intent(in) :: this
        double precision, dimension(nvar) :: yt,a,b,c,d,e,f,g,o,p 
        double precision:: tint,tt
        integer :: i
        tint=t
        do i=1,ns
          call intfunc1_NonlinearLens(t,y,f,this,Bchg,Bmass)
          a=h*f
          yt=y+a/9.d0
          tt=t+h/9.d0
          call intfunc1_NonlinearLens(tt,yt,f,this,Bchg,Bmass)
          b=h*f
          yt=y + (a + 3.d0*b)/24.d0
          tt=t+h/6.d0
          call intfunc1_NonlinearLens(tt,yt,f,this,Bchg,Bmass)
          c=h*f
          yt=y+(a-3.d0*b+4.d0*c)/6.d0
          tt=t+h/3.d0
          call intfunc1_NonlinearLens(tt,yt,f,this,Bchg,Bmass)
          d=h*f
          yt=y + (278.d0*a - 945.d0*b + 840.d0*c + 99.d0*d)/544.d0
          tt=t+.5d0*h
          call intfunc1_NonlinearLens(tt,yt,f,this,Bchg,Bmass)
          e=h*f
          yt=y + (-106.d0*a+273.d0*b-104.d0*c-107.d0*d+48.d0*e)/6.d0
          tt = t+2.d0*h/3.d0
          call intfunc1_NonlinearLens(tt,yt,f,this,Bchg,Bmass)
          g=h*f
          yt = y+(110974.d0*a-236799.d0*b+68376.d0*c+ &
                  103803.d0*d-10240.d0*e + 1926.d0*g)/45648.d0
          tt = t + 5.d0*h/6.d0
          call intfunc1_NonlinearLens(tt,yt,f,this,Bchg,Bmass)
          o=h*f
          yt = y+(-101195.d0*a+222534.d0*b-71988.d0*c-  &
                  26109.d0*d-2.d4*e-72.d0*g+22824.d0*o)/25994.d0
          tt = t + h
          call intfunc1_NonlinearLens(tt,yt,f,this,Bchg,Bmass)
          p=h*f
          y = y+(41.d0*a+216.d0*c+27.d0*d+ &
                 272.d0*e+27.d0*g+216.d0*o+41.d0*p)/840.d0
          t=tint+i*h
        enddo

        end subroutine rk6i_NonlinearLens

        subroutine intfunc1_NonlinearLens(t,y,f,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,Bchg,Bmass
        double precision, dimension(:), intent(in) :: y
        double precision, dimension(:), intent(out) :: f
        type (NonlinearLens), intent(in) :: this
        double precision :: gamma0,beta0,gbet
        double precision :: zedge,s11,s33,s55

        zedge = this%Param(1)

        f(1) = 0.0
        f(2) = 0.0
        f(3) = 0.0
        f(4) = 0.0

        ! synchronous particle:
        gamma0=-y(6)
        beta0=sqrt(gamma0**2-1.)/gamma0
        gbet=sqrt((gamma0-1.0)*(gamma0+1.0))

        f(5)=1.0/(beta0*Scxl)
        f(6)=0.0

        ! matrix elements
        s11=0.0
        s33=0.0 
        s55=  0.0

        f(7)=y(8)/Scxl
        f(8)=-s11*y(7)
        f(9)=y(10)/Scxl
        f(10)=-s11*y(9)
        f(11)=y(12)/Scxl
        f(12)=-s33*y(11)
        f(13)=y(14)/Scxl
        f(14)=-s33*y(13)
        f(15)=y(16)/Scxl
        f(16)=-s55*y(15)
        f(17)=y(18)/Scxl
        f(18)=-s55*y(17)

        end subroutine intfunc1_NonlinearLens

        subroutine  getfld_NonlinearLens(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (NonlinearLens), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,bgrad,zedge

        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0
        extfld(4) = 0.0
        extfld(5) = 0.0
        extfld(6) = 0.0

        end subroutine getfld_NonlinearLens

        !get external field with displacement and rotation errors.
        subroutine  getflderr_NonlinearLens(pos,extfld,this,dx,dy,anglex,&
                                         angley,anglez)
        implicit none  
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (NonlinearLens), intent(in) :: this
        double precision, intent(in) :: dx,dy,anglex,angley,anglez
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,bgrad,zedge
        double precision, dimension(3) :: temp,tmp

        zedge = this%Param(1)
        if(this%Param(3).gt.1.0e-5) then
          zz = pos(3)-zedge
!          call getfldfrg_Quadrupole(zz,this,bgrad)
        else
          bgrad = this%Param(2)
        endif

!        dx = this%Param(5)
!        dy = this%Param(6)
!        anglex = this%Param(7)
!        angley = this%Param(8)
!        anglez = this%Param(9)

        temp(1) = pos(1) - dx
        temp(2) = pos(2) - dy
        tmp(1) = temp(1)*cos(anglez) + temp(2)*sin(anglez)
        tmp(2) = -temp(1)*sin(anglez) + temp(2)*cos(anglez)
        tmp(3) = pos(3) - zedge
        temp(1) = tmp(1)*cos(angley)+tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = -tmp(1)*sin(angley)+tmp(3)*cos(angley)
        tmp(1) = temp(1)
        tmp(2) = temp(2)*cos(anglex)+temp(3)*sin(anglex)
        tmp(3) = -temp(2)*sin(anglex)+temp(3)*cos(anglex)
        zz = tmp(3)

        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0
        extfld(4) = bgrad*tmp(2)
        extfld(5) = bgrad*tmp(1)
        extfld(6) = 0.0

        tmp(1) = extfld(4)
        tmp(2) = extfld(5)*cos(anglex)-extfld(6)*sin(anglex)
        tmp(3) = extfld(5)*sin(anglex)+extfld(6)*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(4) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(5) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(6) = temp(3)

        end subroutine getflderr_NonlinearLens

        subroutine propagator_NonlinearLens(t,tau,this,refpt,Nplc,pts,qmass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: Nplc
        double precision, intent(inout) :: t
        double precision, intent(in) :: tau,qmass
        double precision, dimension(6), intent(inout) :: refpt
        double precision, dimension(6) :: coord
        double precision, dimension(6,Nplc) :: invariants
        type (NonlinearLens), intent(in) :: this
        double precision, pointer, dimension(:,:) :: pts
        real*8 :: gambet0,zedge,tn,cn,bn,dn,mu0,f0,l0,sn, &
                  knll,cnll,b,d,beta0,ds,sedge,ds2,smid, &
                  an,xn,yn,pxn,pyn,u,v,Hinv,Iinv,f,g,mismatch,test

        real*8 :: snf,bnf,anf,cnllf,knllf
        integer  :: i,j,nsg,mpstp,nn
        
        zedge = this%Param(1)    !Location of NLL element entry
        tn = this%Param(2)       !Dimensionless strength of NLL
        cn = this%Param(3)       !Dimensional parameter of NLL
        mu0 = this%Param(4)      !Phase advance across NLL
        l0 = this%Length         !Length of NLL element
        nsg = this%Nseg          !Number of SC kicks
        mpstp = this%Mapstp      !Number of map steps/slice
        nn = (nsg+1)*mpstp       !Total number of NLL segments
        b = -pi/2.0d0                !Dimensionless parameter b/t
        d = 0.0d0            !Dimensionless parameter d/t

        l0 = dabs(l0) !To allow reversible tracking
        f0 = l0/4.0*(1.0+1.0/tan(pi*mu0)**2)  !Focal length of optics section
        gambet0 = sqrt(refpt(6)**2-1.0d0)     !Relativistic beta*gamma
        beta0 = sqrt(1.0d0-1.0d0/(refpt(6)**2)) !Relativistic beta

        !print*, 'zedge,tn,cn,mu0,l0:',zedge,tn,cn,mu0,l0
        !print*, 'nsg,mpstp,nn,f0:',nsg,mpstp,nn,f0
        !print*, 'Scxl,gambet0 = ',Scxl,gambet0

        sedge = t-zedge          !Location of SC slice entry
        if(tau<0.0d0) sedge = l0+sedge  !Begin at end of slice (reverse tracking)
        ds = tau/dble(mpstp)     !Size of each map step
        ds2 = ds/2.0d0
        !print*, 'Step thickness:',ds
        !print*, 'Parameters sn,bn,knll,cnll,an:'
        do j = 1, mpstp
          sn = sedge + ds*(j-0.5d0)  !Location of jth map kick
          bn=l0*(1-sn*(l0-sn)/l0/f0)/sqrt(1.0-(1.0-l0/2.0/f0)**2)
          an=((l0-2.0d0*sn)/f0)/(2.0d0*sqrt(1.0-(1.0-l0/2.0/f0)**2))
          knll=ds*tn*cn**2/bn
          cnll=cn*sqrt(bn)
          !write(*,20) sn,bn,knll,cnll,an

          do i = 1, Nplc
            coord(1) = pts(1,i)*Scxl          !Convert to normalized units for update
            coord(2) = pts(2,i)/gambet0
            coord(3) = pts(3,i)*Scxl
            coord(4) = pts(4,i)/gambet0
            coord(5) = pts(5,i)*Scxl
            coord(6) = pts(6,i)/gambet0
            call  DriftPropagator(ds2,beta0,gambet0,coord)    !Half step for drift
!   Old tracking algorithm should be commented here:
      !      call NonlinearLensPropagator(knll,cnll,b,d,coord,u,v)   !Full step in (px,py)
!   Complex potential tracking algorithm:
            call NonlinearLensPropagatorCmplx(knll,cnll,coord(1:4))  !Complex version of step in (px,py)
            call DriftPropagator(ds2,beta0,gambet0,coord)    !Half step for drift
            pts(1,i) = coord(1)/Scxl                   !Convert back to internal units
            pts(2,i) = coord(2)*gambet0
            pts(3,i) = coord(3)/Scxl 
            pts(4,i) = coord(4)*gambet0
            pts(5,i) = coord(5)/Scxl
            pts(6,i) = coord(6)*gambet0
!  Move diagnostic computation to coincide with particle output location
            snf = sedge + ds*j  !Location of jth step output
            bnf=l0*(1-snf*(l0-snf)/l0/f0)/sqrt(1.0-(1.0-l0/2.0/f0)**2)
            anf=((l0-2.0d0*snf)/f0)/(2.0d0*sqrt(1.0-(1.0-l0/2.0/f0)**2))
            knllf=ds*tn*cn**2/bnf
            cnllf=cn*sqrt(bnf)
!   Compute diagnostic quantities (temporarily added the "if" statement for speed)
          smid = l0/2.0d0
        if(abs(snf-smid)<1.0d-4) then
            xn = coord(1)/cnllf
            yn = coord(3)/cnllf
            pxn = coord(2)*sqrt(bnf)/cn + anf*xn
            pyn = coord(4)*sqrt(bnf)/cn + anf*yn
            call InvariantPotentials(xn,yn,Hinv,Iinv)
            Hinv = (xn**2+yn**2+pxn**2+pyn**2)/2.d0+tn*Hinv
            Iinv = (xn*pyn-yn*pxn)**2+pxn**2+xn**2+tn*Iinv
            
            !<<<<<<<<<<<<<<<<<< Kilean <<<<<<<<<<<<<<<<<<<<<
            if(Hinv.ne.Hinv) print*,'Hinv NaN, xn,yn = ',xn,yn
            if(Iinv.ne.Iinv) print*,'Iinv NaN, xn,yn = ',xn,yn
            if(xn.ne.xn) print*,'xn NaN, cnllf=',cnllf
            if(yn.ne.yn) print*,'yn NaN, cnllf=',cnllf
            if(pxn.ne.pxn) print*,'pxn NaN, anf, bnf=',anf,bnf
            if(pyn.ne.pyn) print*,'pyn NaN, anf, bnf=',anf,bnf
            if(Iinv<0) Iinv=0d0
            !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            
!   Old diagnostic algorithm should be commented here:
      !      Hinv = u*sqrt(u**2-1.0d0)*(d+acosh(u))
      !      Hinv = Hinv + v*sqrt(1.0d0-v**2)*(b+acos(v))
      !      Hinv = Hinv/(u**2-v**2)
      !      Hinv = (xn**2+yn**2+pxn**2+pyn**2)/2.d0-tn*Hinv
      !      f = u**2*(u**2-1.0d0)/2.0d0 - tn*u*sqrt(u**2-1.0d0)*(d+acosh(u))
      !      g = v**2*(1.0d0-v**2)/2.0d0 - tn*v*sqrt(1.0d0-v**2)*(b+acos(v))
      !      Iinv = 2.0d0*(f*v**2+g*u**2)/(u**2-v**2)
      !      Iinv = (xn*pyn-yn*pxn)**2+pxn**2+Iinv
!            write(71,20) snf+zedge,xn,pxn,yn,pyn,Hinv,Iinv
!            write(71,20) snf+zedge,an,bn,Hinv,Iinv
!   End old diagnostic algorithm
            invariants(1,i) = Hinv
!            invariants(2,i) = Iinv
            invariants(2,i) = dsqrt(Iinv)  !Use this for benchmark with CH.
            invariants(3,i) = xn
            invariants(4,i) = pxn
            invariants(5,i) = yn
            invariants(6,i) = pyn
!   Test for occurrence of NaN:
            test = Hinv*dsqrt(Iinv)*xn*pxn*yn*pyn
            if(test.ne.test) then
              write(*,*) 'NaN encountered (particle,s):'
              write(*,*) i,snf
              stop
            endif
         endif
!   End computation of diagnostic quantities

          enddo

!        t = t + tau
20      format(7(1x,g20.12))
        enddo

!  Use the following "if" statement if we want to write out diagnostics
!  only near the midpoint of the insert.
        smid = l0/2.0d0
         if(abs(sedge-smid+tau)<1.0d-4) then
          call diagnostics_NonlinearLens(t+tau,invariants,Nplc)
!          call mismatch_NonlinearLens(t+tau,invariants,Nplc)  !Commented for speed
         endif

        refpt(5) = refpt(5) + tau/beta0/Scxl

        end subroutine propagator_NonlinearLens

        subroutine diagnostics_NonlinearLens(t,invariants,Nplc)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: Nplc
        double precision, intent(in) :: t
        double precision, intent(in), dimension(6,Nplc) :: invariants
        real*8 :: Hinv1,Iinv1,Hinv2,Iinv2,HIinv,sigH,sigI,covHI
        real*8 :: x1,x2,y1,y2,px1,px2,py1,py2,xpx,ypy
        real*8 :: sigx,sigy,sigpx,sigpy,covxpx,covypy
        double precision, dimension(15) :: tmplc,tmpgl
        integer  :: j,my_rank,ierr,Nptot

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

        tmplc = 0.0d0
        tmpgl = 0.0d0
        Hinv1 = 0.0d0
        Iinv1 = 0.0d0
        HIinv = 0.0d0
        Hinv2 = 0.0d0
        Iinv2 = 0.0d0
        x1 = 0.0d0
        x2 = 0.0d0
        y1 = 0.0d0
        y2 = 0.0d0
        px1 = 0.0d0
        px2 = 0.0d0
        py1 = 0.0d0
        py2 = 0.0d0
        xpx = 0.0d0
        ypy = 0.0d0
        do j = 1,Nplc
          Hinv1 = Hinv1 + invariants(1,j)
          Iinv1 = Iinv1 + invariants(2,j)
          HIinv = HIinv + invariants(1,j)*invariants(2,j)
          Hinv2 = Hinv2 + invariants(1,j)**2
          Iinv2 = Iinv2 + invariants(2,j)**2
          x1 = x1 + invariants(3,j)
          x2 = x2 + invariants(3,j)**2
          px1 = px1 + invariants(4,j)
          px2 = px2 + invariants(4,j)**2
          y1 = y1 + invariants(5,j)
          y2 = y2 + invariants(5,j)**2
          py1 = py1 + invariants(6,j)
          py2 = py2 + invariants(6,j)**2
          xpx = xpx + invariants(3,j)*invariants(4,j)
          ypy = ypy + invariants(5,j)*invariants(6,j)
        enddo
        tmplc(1) = Hinv1
        tmplc(2) = Iinv1
        tmplc(3) = HIinv
        tmplc(4) = Hinv2
        tmplc(5) = Iinv2
        tmplc(6) = x1
        tmplc(7) = x2
        tmplc(8) = px1
        tmplc(9) = px2
        tmplc(10) = y1
        tmplc(11) = y2
        tmplc(12) = py1
        tmplc(13) = py2
        tmplc(14) = xpx
        tmplc(15) = ypy
        call MPI_REDUCE(tmplc,tmpgl,15,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(Nplc,Nptot,1,MPI_INTEGER,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if(my_rank.eq.0) then
        Hinv1 = tmpgl(1)/dble(Nptot)
        Iinv1 = tmpgl(2)/dble(Nptot)
        HIinv = tmpgl(3)/dble(Nptot)
        Hinv2 = tmpgl(4)/dble(Nptot)
        Iinv2 = tmpgl(5)/dble(Nptot)
        sigH = dsqrt(Hinv2-Hinv1**2)
        sigI = dsqrt(Iinv2-Iinv1**2)
        covHI = HIinv - Hinv1*Iinv1
        x1 = tmpgl(6)/dble(Nptot)
        x2 = tmpgl(7)/dble(Nptot)
        px1 = tmpgl(8)/dble(Nptot)
        px2 = tmpgl(9)/dble(Nptot)
        y1 = tmpgl(10)/dble(Nptot)
        y2 = tmpgl(11)/dble(Nptot)
        py1 = tmpgl(12)/dble(Nptot)
        py2 = tmpgl(13)/dble(Nptot)
        xpx = tmpgl(14)/dble(Nptot)
        ypy = tmpgl(15)/dble(Nptot)
        sigx = dsqrt(x2-x1**2)
        sigpx = dsqrt(px2-px1**2)
        sigy = dsqrt(y2-y1**2)
        sigpy = dsqrt(py2-py1**2)
        covxpx = xpx - x1*px1
        covypy = ypy - y1*py1
        write(81,20) t,Hinv1,Iinv1,sigH,sigI,covHI
        write(82,20) t,x1,px1,y1,py1,sigx,sigpx,sigy,sigpy,covxpx,covypy
        flush(81)
        flush(82)
      endif

20      format(11(1x,g20.12))

        end subroutine diagnostics_NonlinearLens

        subroutine mismatch_NonlinearLens(t,invariants,Nplc)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: Nplc
        double precision, intent(in) :: t
        double precision, intent(in), dimension(6,Nplc) :: invariants
        double precision, dimension(4,4) :: sigma_d,sigma_p,sigma
        double precision, dimension(4,4) :: sigmalc,sigmagl
        double precision, dimension(2) :: det
        double precision:: determinant,tr,mismatch
        integer  :: j,k,l,lda,n,info,job,my_rank,ierr,Nptot

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

        job = 11
        sigmalc = 0.0d0
        sigmagl = 0.0d0
        sigma_d = 0.0d0
        sigma = 0.0d0
        do j = 1,Nplc
          do l = 1,4
            do k = 1,l
              sigmalc(k,l)=sigmalc(k,l)+invariants(2+l,j)*invariants(2+k,j)
            enddo
          enddo
        enddo
        call MPI_REDUCE(sigmalc,sigmagl,16,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(Nplc,Nptot,1,MPI_INTEGER,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
        sigma = sigmagl/dble(Nptot)

!  Test matrix
!        sigma = 0.0d0
!        sigma_d = 0.0d0
!        sigma(1,1) = 0.5d0
!        sigma(1,2) = 0.2d0
!        sigma(2,2) = 4.0d0
!        sigma(3,3) = 0.5d0
!        sigma(3,4) = 0.2d0
!        sigma(4,4) = 4.0d0
!        sigma(2,1) = sigma(1,2)
!        sigma(4,3) = sigma(3,4)
!        sigma_d = sigma
!  Nominal 4D Twiss matrix for nonlinear KV
      if(my_rank.eq.0) then
!        sigma_d(1,1) = 0.607210558112741
!        sigma_d(1,2) = -7.189922950934375E-004
!        sigma_d(2,2) = 0.854230738529633
!        sigma_d(1,3) = -2.580942629625784E-003
!        sigma_d(2,3) = -3.487490180048264E-003
!        sigma_d(3,3) = 2.25445980775678
!        sigma_d(1,4) = -2.490184832359349E-004
!        sigma_d(2,4) = 1.205146145443822E-003
!        sigma_d(3,4) = -1.702149308806317E-003
!        sigma_d(4,4) = 0.855164869938001
!  Nominal 4D Twiss matrix for KV
!        sigma_d(1,1) = 1.0d0
!        sigma_d(2,2) = 1.0d0
!        sigma_d(3,3) = 1.0d0
!        sigma_d(4,4) = 1.0d0
!  Nominal 4D Twiss matrix for nonlinear KV test
!        sigma_d(1,1) = 0.478522840872335
!        sigma_d(1,2) = 2.037961684476454E-003
!        sigma_d(2,2) = 0.810448309500901
!        sigma_d(1,3) = -4.336749860086196E-004
!        sigma_d(2,3) = -6.879253530927047E-003
!        sigma_d(3,3) = 3.21482559940535
!        sigma_d(1,4) = -4.332846363060930E-004
!        sigma_d(2,4) = -1.251432543531919E-004
!        sigma_d(3,4) = -8.333613637319727E-003
!        sigma_d(4,4) = 0.802119548766217
!  New nominal 4D Twiss matrix for nonlinear KV w SC
        sigma_d(1,1) = 0.504540711130766
        sigma_d(1,2) = 8.136364265543093E-003
        sigma_d(2,2) = 0.831334243885491
        sigma_d(1,3) = -1.132190353023146E-002
        sigma_d(2,3) = 2.577385793436604E-003
        sigma_d(3,3) = 2.80040438958216
        sigma_d(1,4) = -2.512397796421949E-003
        sigma_d(2,4) = 3.941467993038517E-003
        sigma_d(3,4) = 8.713593277802075E-003
        sigma_d(4,4) = 0.851621184371757
!  
        call dpofa(sigma,4,4,info)
        call dpodi(sigma,4,4,det,job)
        do l=1,4
         do k=l,4
           sigma(k,l) = sigma(l,k)
           sigma_d(k,l) = sigma_d(l,k)
         enddo
        enddo
!        write(*,*) 'det = ',det       
!        write(*,*) 'sigma, sigma_d = '
!        do l=1,4
!          do k=1,4
!           write(*,*) k,l,sigma(k,l),sigma_d(k,l)
!          enddo
!        enddo
        sigma_p = matmul(dble(sigma),dble(sigma_d))
!        write(*,*) 'product = '
!        do l=1,4
!          do k=1,l
!           write(*,*) k,l,sigma_p(k,l)
!          enddo
!        enddo
        determinant = det(1)*10.0**(det(2))
        tr = 0.0d0
        do k=1,4
           tr = tr + sigma_p(k,k)
        enddo
        mismatch = determinant**(1.0d0/4.0d0)*tr/4.0d0
        write(83,*) t,mismatch
      endif
        
        call flush(83)

        end subroutine mismatch_NonlinearLens
        
      end module NonlinearLensclass
