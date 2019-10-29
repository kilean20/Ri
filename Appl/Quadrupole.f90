!----------------------------------------------------------------
! (c) Copyright, 2007 by the Regents of the University of California.
! Quadrupoleclass: Quadrupole beam line element class
!             in Lattice module of APPLICATION layer.
! Version: 1.0
! Author: Ji Qiang, Robert Ryne, LANL, 7/13/01
! Description: This class defines the linear transfer map and field
!              for the quadrupole beam line elment.
! Comments:
!----------------------------------------------------------------
      module Quadrupoleclass
        use PhysConstclass
        use Dataclass
        integer, private, parameter :: Nparam = 9
        type Quadrupole
          !Itype = 1
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
          ! Param(1) : zedge
          !      (2) : quad gradient or scaling const. of fringe field
          !      (3) : file ID
          !      (4) : radius
          !      (5) : x misalignment error
          !      (6) : y misalignment error
          !      (7) : rotation error x
          !      (8) : rotation error y
          !      (9) : rotation error z
        end type Quadrupole
        interface getparam_Quadrupole
          module procedure getparam1_Quadrupole,  &
                          getparam2_Quadrupole,   &
                          getparam3_Quadrupole
        end interface
        interface setparam_Quadrupole
          module procedure setparam1_Quadrupole,  &
                          setparam2_Quadrupole, setparam3_Quadrupole
        end interface
      contains
        subroutine construct_Quadrupole(this,numseg,nmpstp,type,blength)
        implicit none
        type (Quadrupole), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_Quadrupole
   
        subroutine setparam1_Quadrupole(this,i,value)
        implicit none
        type (Quadrupole), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_Quadrupole

        subroutine setparam2_Quadrupole(this,values)
        implicit none
        type (Quadrupole), intent(inout) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param = values

        end subroutine setparam2_Quadrupole

        subroutine setparam3_Quadrupole(this,numseg,nmpstp,type,blength)
        implicit none
        type (Quadrupole), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_Quadrupole
   
        subroutine getparam1_Quadrupole(this,i,blparam) 
        implicit none 
        type (Quadrupole), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_Quadrupole
  
        subroutine getparam2_Quadrupole(this,blparams)
        implicit none
        type (Quadrupole), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams = this%Param

        end subroutine getparam2_Quadrupole

        subroutine getparam3_Quadrupole(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (Quadrupole), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_Quadrupole
       
        subroutine maplinear_Quadrupole(t,tau,xm,this,refpt,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,tau,Bchg,Bmass
        double precision, dimension(6,6), intent(out) :: xm
        double precision, dimension(6), intent(inout) :: refpt
        type (Quadrupole), intent(in) :: this
        double precision, dimension(18) :: y
        double precision :: h,t0,gi,betai,ui,squi,sqi3
        double precision :: dlti,thli,gf,betaf,uf,squf,sqf3
        double precision :: dltf,thlf,zedge,escale,tfin
        integer  :: mpstp

        xm = 0.0

        zedge = this%Param(1)
        escale = this%Param(2)
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

        call rk6i_Quadrupole(h,mpstp,t0,y,18,this,Bchg,Bmass)

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
        !print*,"xm: ",xm(1,1),xm(2,1),xm(1,2),xm(2,2),xm(3,3),xm(4,3),xm(4,4)

        end subroutine maplinear_Quadrupole

        subroutine rk6i_Quadrupole(h,ns,t,y,nvar,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ns,nvar
        double precision, intent(inout) :: t
        double precision, intent(in) :: h,Bchg,Bmass
        double precision, dimension(nvar), intent(inout) :: y
        type (Quadrupole), intent(in) :: this
        double precision, dimension(nvar) :: yt,a,b,c,d,e,f,g,o,p 
        double precision:: tint,tt
        integer :: i
        tint=t
        do i=1,ns
          call intfunc1_Quadrupole(t,y,f,this,Bchg,Bmass)
          a=h*f
          yt=y+a/9.d0
          tt=t+h/9.d0
          call intfunc1_Quadrupole(tt,yt,f,this,Bchg,Bmass)
          b=h*f
          yt=y + (a + 3.d0*b)/24.d0
          tt=t+h/6.d0
          call intfunc1_Quadrupole(tt,yt,f,this,Bchg,Bmass)
          c=h*f
          yt=y+(a-3.d0*b+4.d0*c)/6.d0
          tt=t+h/3.d0
          call intfunc1_Quadrupole(tt,yt,f,this,Bchg,Bmass)
          d=h*f
          yt=y + (278.d0*a - 945.d0*b + 840.d0*c + 99.d0*d)/544.d0
          tt=t+.5d0*h
          call intfunc1_Quadrupole(tt,yt,f,this,Bchg,Bmass)
          e=h*f
          yt=y + (-106.d0*a+273.d0*b-104.d0*c-107.d0*d+48.d0*e)/6.d0
          tt = t+2.d0*h/3.d0
          call intfunc1_Quadrupole(tt,yt,f,this,Bchg,Bmass)
          g=h*f
          yt = y+(110974.d0*a-236799.d0*b+68376.d0*c+ &
                  103803.d0*d-10240.d0*e + 1926.d0*g)/45648.d0
          tt = t + 5.d0*h/6.d0
          call intfunc1_Quadrupole(tt,yt,f,this,Bchg,Bmass)
          o=h*f
          yt = y+(-101195.d0*a+222534.d0*b-71988.d0*c-  &
                  26109.d0*d-2.d4*e-72.d0*g+22824.d0*o)/25994.d0
          tt = t + h
          call intfunc1_Quadrupole(tt,yt,f,this,Bchg,Bmass)
          p=h*f
          y = y+(41.d0*a+216.d0*c+27.d0*d+ &
                 272.d0*e+27.d0*g+216.d0*o+41.d0*p)/840.d0
          t=tint+i*h
        enddo

        end subroutine rk6i_Quadrupole

        subroutine intfunc1_Quadrupole(t,y,f,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,Bchg,Bmass
        double precision, dimension(:), intent(in) :: y
        double precision, dimension(:), intent(out) :: f
        type (Quadrupole), intent(in) :: this
        double precision :: gamma0,beta0,gbet,s11,s33,s55
        double precision :: zedge,qmcc,brho,bgrad,zz

        zedge = this%Param(1)
        if(this%Param(3).gt.1.0e-5) then
          zz = t - zedge
          call getfldfrg_Quadrupole(zz,this,bgrad)
        else
          bgrad = this%Param(2)
        endif
        qmcc = Bchg/Bmass

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
        brho=gbet/Clight/qmcc
        s11=bgrad/brho*Scxl
        s33=-bgrad/brho*Scxl
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

        end subroutine intfunc1_Quadrupole

        !get external field with displacement and rotation errors.
        subroutine  getflderr_Quadrupole(pos,extfld,this,dx,dy,anglex,&
                                         angley,anglez)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Quadrupole), intent(in) :: this
        double precision, intent(in) :: dx,dy,anglex,angley,anglez
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,bgrad,zedge
        double precision, dimension(3) :: temp,tmp

        zedge = this%Param(1)
        if(this%Param(3).gt.1.0e-5) then
          zz = pos(3)-zedge
          call getfldfrg_Quadrupole(zz,this,bgrad)
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

        end subroutine getflderr_Quadrupole
        
        !get external field without displacement and rotation errors.
        subroutine  getfld_Quadrupole(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Quadrupole), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,bgrad,zedge

        zedge = this%Param(1)
        zz=pos(3)-zedge
        if(this%Param(3).gt.0.0) then
          call getfldfrg_Quadrupole(zz,this,bgrad)
        else
          bgrad = this%Param(2)
        endif

        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0
        extfld(4) = bgrad*pos(2)
        extfld(5) = bgrad*pos(1)
        extfld(6) = 0.0

        end subroutine getfld_Quadrupole

        !interpolate the field from the SC rf cavity onto bunch location.
        subroutine getfldfrg_Quadrupole(zz,this,bgrad)
        implicit none
        include 'mpif.h'
        type (Quadrupole), intent(in) :: this
        double precision, intent(in) :: zz
        double precision, intent(out) :: bgrad
        double precision:: hstep,slope
        integer :: klo,khi,k
        integer :: my_rank,ierr

        klo=1
        khi=Ndata
1       if(khi-klo.gt.1) then
          k=(khi+klo)/2
          if(zdat(k).gt.zz)then
             khi=k
          else
             klo=k
          endif
          goto 1
        endif
        hstep=zdat(khi)-zdat(klo)
        slope=(edat(khi)-edat(klo))/hstep
        bgrad = this%Param(2)*(edat(klo)+slope*(zz-zdat(klo)))

        end subroutine getfldfrg_Quadrupole

        subroutine transfmap_Quadrupole(t,tau,this,refpt,Nplc,pts)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: Nplc
        double precision, intent(inout) :: t
        double precision, intent(in) :: tau
        double precision, dimension(6), intent(inout) :: refpt
        double precision, dimension(6) :: tmp
        type (Quadrupole), intent(in) :: this
        double precision, pointer, dimension(:,:) :: pts
        real*8 :: xm11,xm12,xm21,xm22,xm33,xm34,xm43,xm44,gam,gambetz,&
                  betaz,beta0,rtkstrzz,rtkstr,kstr
        integer  :: i

        do i = 1, Nplc
          gam = -refpt(6) - pts(6,i)
          gambetz = sqrt(gam**2-1.0-pts(2,i)**2-pts(4,i)**2)
          betaz = gambetz/gam
          beta0 = sqrt(1.-1./(refpt(6)**2))
          kstr = pts(7,i)*2.997928e8/gambetz*this%Param(2)
          rtkstr = sqrt(abs(kstr))
          rtkstrzz = rtkstr*tau
          if(kstr.gt.0.0) then
            xm11 = cos(rtkstrzz)
            xm12 = sin(rtkstrzz)/rtkstr
            xm21 = -rtkstr*sin(rtkstrzz)
            xm22 = cos(rtkstrzz)
            xm33 = cosh(rtkstrzz)
            xm34 = sinh(rtkstrzz)/rtkstr
            xm43 = sinh(rtkstrzz)*rtkstr
            xm44 = cosh(rtkstrzz)
          else if(kstr.lt.0.0) then
            xm11 = cosh(rtkstrzz)
            xm12 = sinh(rtkstrzz)/rtkstr
            xm21 = rtkstr*sinh(rtkstrzz)
            xm22 = cosh(rtkstrzz)
            xm33 = cos(rtkstrzz)
            xm34 = sin(rtkstrzz)/rtkstr
            xm43 = -sin(rtkstrzz)*rtkstr
            xm44 = cos(rtkstrzz)
          else
            xm11 = 1.0d0
            xm12 = tau
            xm21 = 0.0d0
            xm22 = 1.0d0
            xm33 = 1.0d0
            xm34 = tau
            xm43 = 0.0
            xm44 = 1.0d0
          endif
          tmp(1) = xm11*pts(1,i)+xm12*pts(2,i)/gambetz/Scxl
          tmp(2) = gambetz*Scxl*xm21*pts(1,i)+xm22*pts(2,i)
          tmp(3) = xm33*pts(3,i)+xm34*pts(4,i)/gambetz/Scxl
          tmp(4) = gambetz*Scxl*xm43*pts(3,i)+xm44*pts(4,i)
          tmp(5) = pts(5,i) + (1.0/betaz-1.0/beta0)*tau/Scxl 
          tmp(6) = pts(6,i)
          pts(1,i) = tmp(1)
          pts(2,i) = tmp(2)
          pts(3,i) = tmp(3)
          pts(4,i) = tmp(4)
          pts(5,i) = tmp(5)
          pts(6,i) = tmp(6)
        enddo
        refpt(5) = refpt(5) + tau/beta0/Scxl
!        t = t + tau

        end subroutine transfmap_Quadrupole

        subroutine transfmapK_Quadrupole(tt,tau,this,refpt,Nplc,pts,qmass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: Nplc
        double precision, intent(inout) :: tt
        double precision, intent(in) :: tau,qmass
        double precision, dimension(6), intent(inout) :: refpt
        double precision, dimension(6) :: tmp
        type (Quadrupole), intent(in) :: this
        double precision, pointer, dimension(:,:) :: pts
        real*8 :: xm11,xm12,xm21,xm22,xm33,xm34,xm43,xm44,gam,gambetz,&
                  betaz,beta0,rtkstrzz,rtkstr,kstr,gambet0
        integer  :: i
        real*8 :: t,cs,ss

        gambet0 = sqrt(refpt(6)**2-1.0d0)
        beta0 = sqrt(1.0d0-1.0d0/(refpt(6)**2))
        do i = 1, Nplc
          gam = -refpt(6) - pts(6,i)
          gambetz = sqrt(gam**2-1.0d0-pts(2,i)**2-pts(4,i)**2)
          betaz = gambetz/gam
          !each reference particle momentum
          gambetz = sqrt(gam**2-1.0d0) 
!          kstr = pts(7,i)*2.997928e8/gambetz*this%Param(2)
!Param(2) is the K defined in MAD, i.e. G/Brho
          kstr = pts(7,i)/qmass*gambet0/gambetz*this%Param(2)
          rtkstr = sqrt(abs(kstr))
          rtkstrzz = rtkstr*tau
          if(kstr.gt.0.0) then
            cs = cos(rtkstrzz)
            ss = sin(rtkstrzz)
            xm11 = cs
            xm12 = ss/rtkstr
            xm21 = -rtkstr*ss
            xm22 = cs
            !xm11 = cos(rtkstrzz)
            !xm12 = sin(rtkstrzz)/rtkstr
            !xm21 = -rtkstr*sin(rtkstrzz)
            !xm22 = cos(rtkstrzz)
            t = exp(rtkstrzz)
            cs= (t*t +1.0d0) / (t + t)
            ss= (t*t -1.0d0) / (t + t) 
            xm33 = cs
            xm34 = ss/rtkstr
            xm43 = ss*rtkstr
            xm44 = cs
            !xm33 = cosh(rtkstrzz)
            !xm34 = sinh(rtkstrzz)/rtkstr
            !xm43 = sinh(rtkstrzz)*rtkstr
            !xm44 = cosh(rtkstrzz)
          else if(kstr.lt.0.0) then
            t = exp(rtkstrzz)
            cs= (t*t +1.0d0) / (t + t)
            ss= (t*t -1.0d0) / (t + t) 
            xm11 = cs
            xm12 = ss/rtkstr
            xm21 = rtkstr*ss
            xm22 = cs
            !xm11 = cosh(rtkstrzz)
            !xm12 = sinh(rtkstrzz)/rtkstr
            !xm21 = rtkstr*sinh(rtkstrzz)
            !xm22 = cosh(rtkstrzz)
            cs = cos(rtkstrzz)
            ss = sin(rtkstrzz)
            xm33 = cs
            xm34 = ss/rtkstr
            xm43 = -ss*rtkstr
            xm44 = cs
            !xm33 = cos(rtkstrzz)
            !xm34 = sin(rtkstrzz)/rtkstr
            !xm43 = -sin(rtkstrzz)*rtkstr
            !xm44 = cos(rtkstrzz)
          else
            xm11 = 1.0d0
            xm12 = tau
            xm21 = 0.0d0
            xm22 = 1.0d0
            xm33 = 1.0d0
            xm34 = tau
            xm43 = 0.0
            xm44 = 1.0d0
          endif
          tmp(1) = xm11*pts(1,i)+xm12*pts(2,i)/gambetz/Scxl
          tmp(2) = gambetz*Scxl*xm21*pts(1,i)+xm22*pts(2,i)
          tmp(3) = xm33*pts(3,i)+xm34*pts(4,i)/gambetz/Scxl
          tmp(4) = gambetz*Scxl*xm43*pts(3,i)+xm44*pts(4,i)
          tmp(5) = pts(5,i) + (1.0/betaz-1.0/beta0)*tau/Scxl 
          tmp(6) = pts(6,i)
          pts(1,i) = tmp(1)
          pts(2,i) = tmp(2)
          pts(3,i) = tmp(3)
          pts(4,i) = tmp(4)
          pts(5,i) = tmp(5)
          pts(6,i) = tmp(6)
        enddo
        refpt(5) = refpt(5) + tau/beta0/Scxl
!        tt = tt + tau

        end subroutine transfmapK_Quadrupole
        
!<<<<<<<<<<<<<<<<<<<<< Kilean, Quad fringe field map <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        subroutine monomial_maps2D(npt,pdata,C,nx,npx,ny,npy)
          implicit none
          integer, intent(in) :: npt,nx,npx,ny,npy
          double precision, pointer, dimension(:,:) :: pdata
          double precision, intent(in) :: C(npt)
          double precision :: xdata(2,npt), monoFactor, Ceff
          integer :: i
          
          
          xdata = pdata(1:2,1:npt)
          
          ! === x,px kicks ====
          if(nx/=0) then
            if(nx==npx) then
              if(nx==1) then     ! nx=1, npx=1
                do i=1,npt
                  Ceff = C(i)*pdata(3,i)**ny*pdata(4,i)**npy
                  monoFactor = exp(Ceff)
                  pdata(1,i) = pdata(1,i)/monoFactor
                  pdata(2,i) = pdata(2,i)*monoFactor
                end do
              else               ! nx = npx,   nx,npx > 1
                do i=1,npt
                  Ceff = C(i)*pdata(3,i)**ny*pdata(4,i)**npy
                  monoFactor = exp(Ceff*nx*(pdata(1,i)*pdata(2,i))**(nx-1))
                  pdata(1,i) = pdata(1,i)/monoFactor
                  pdata(2,i) = pdata(2,i)*monoFactor
                end do
              endif 
            else if(npx==0) then ! nx>0 npx=0
              do i=1,npt
                Ceff = C(i)*pdata(3,i)**ny*pdata(4,i)**npy
                pdata(2,i) = pdata(2,i) + Ceff*nx*pdata(1,i)**(nx-1)
              end do
            else                 ! nx != npx,  nx>=1 npx>=1
              do i=1,npt
                Ceff = C(i)*pdata(3,i)**ny*pdata(4,i)**npy
                monoFactor = 1d0 +Ceff*(nx-npx)*pdata(1,i)**(nx-1)*pdata(2,i)**(npx-1)
                pdata(1,i) = pdata(1,i)*(monoFactor**(npx*1d0/(npx-nx)))
                pdata(2,i) = pdata(2,i)*(monoFactor**(nx *1d0/(nx-npx)))
              end do
            endif
          else if(npx/=0) then  ! nx=0, npx>0
            do i=1,npt
              Ceff = C(i)*pdata(3,i)**ny*pdata(4,i)**npy
              pdata(1,i) = pdata(1,i) - Ceff*npx*pdata(2,i)**(npx-1)
            end do
          endif
          
          ! === y,py kicks ====
          if(ny/=0) then
            if(ny==npy) then
              if(ny==1) then     ! ny=1, npy=1
                do i=1,npt
                  Ceff = C(i)*xdata(1,i)**nx*xdata(2,i)**npx
                  monoFactor = exp(Ceff)
                  pdata(3,i) = pdata(3,i)/monoFactor
                  pdata(4,i) = pdata(4,i)*monoFactor
                end do
              else               ! ny = npy,   ny,npx > 1
                do i=1,npt
                  Ceff = C(i)*xdata(1,i)**nx*xdata(2,i)**npx
                  monoFactor = exp(Ceff*ny*(pdata(3,i)*pdata(4,i))**(ny-1))
                  pdata(3,i) = pdata(3,i)/monoFactor
                  pdata(4,i) = pdata(4,i)*monoFactor
                end do
              endif 
            else if(npy==0) then ! ny>0 npy=0
              do i=1,npt
                Ceff = C(i)*pdata(1,i)**nx*pdata(2,i)**npx
                pdata(4,i) = pdata(4,i) + Ceff*ny*pdata(3,i)**(ny-1)
              end do
            else                 ! ny != npy,  ny>1 npy>1
              do i=1,npt
                Ceff = C(i)*xdata(1,i)**nx*xdata(2,i)**npx
                monoFactor = 1d0 +Ceff*(ny-npy)*pdata(3,i)**(ny-1)*pdata(4,i)**(npy-1)
                pdata(3,i) = pdata(3,i)*monoFactor**(npy*1d0/(npy-ny))
                pdata(4,i) = pdata(4,i)*monoFactor**(ny *1d0/(ny-npy))
              end do
            endif
          else if(npy/=0) then  ! ny=0, npy>0
            do i=1,npt
              Ceff = C(i)*xdata(1,i)**nx*xdata(2,i)**npx
              pdata(3,i) = pdata(3,i) - Ceff*npy*pdata(4,i)**(npy-1)
            end do
          endif
          
        end subroutine monomial_maps2D
        
        subroutine HardEdgeQuadFringeMap(flagEntrance,Kx,nSteps,refpt,npt,pts,qmass)
        implicit none
        integer, intent(in) :: npt,nSteps
        double precision, intent(in) :: Kx,qmass
        double precision, dimension(6), intent(in) :: refpt
        logical, intent(in) :: flagEntrance
        double precision, pointer, dimension(:,:) :: pts
        integer :: i
        double precision :: gam,gambet(npt),gambet0,kstr(npt)
        
        
        gambet0 = sqrt(refpt(6)**2-1.0d0)
        print*, 'flagEntrance=',flagEntrance
        
        do i = 1, npt
          gam = -refpt(6) - pts(6,i)
          gambet(i) = sqrt(gam**2-1.0d0) 
          kstr(i) = pts(7,i)/qmass*gambet0/gambet(i)*Kx
        end do
        do i = 1, npt
          pts(1,i) = pts(1,i)*Scxl
          pts(3,i) = pts(3,i)*Scxl
          pts(2,i) = pts(2,i)/gambet(i)
          pts(4,i) = pts(4,i)/gambet(i)
        end do
   
        if (flagEntrance) then
          kstr=-kstr/(24d0*nSteps)
        else
          kstr = kstr/(24d0*nSteps)
        endif
        do i=1, nSteps
          call  monomial_maps2D(npt,pts,     kstr,3,1,0,0)
          call  monomial_maps2D(npt,pts,    -kstr,0,0,3,1)
          call  monomial_maps2D(npt,pts, 3d0*kstr,1,1,2,0)
          call  monomial_maps2D(npt,pts,-6d0*kstr,2,0,1,1)
          call  monomial_maps2D(npt,pts, 3d0*kstr,1,1,2,0)
          call  monomial_maps2D(npt,pts,    -kstr,0,0,3,1)
          call  monomial_maps2D(npt,pts,     kstr,3,1,0,0)
        end do
        !kstr = kstr/(12d0*nSteps)
        !do i=1, nSteps
        !  call  monomial_maps2D(npt,pts, kstr,3,1,0,0)
        !  call  monomial_maps2D(npt,pts,-kstr,0,0,3,1)
        !  call  monomial_maps2D(npt,pts, 3*kstr,1,1,2,0)
        !  call  monomial_maps2D(npt,pts,-3*kstr,2,0,1,1)
        !end do
        
        
        do i = 1, npt
          pts(1,i) = pts(1,i)/Scxl
          pts(3,i) = pts(3,i)/Scxl
          pts(2,i) = pts(2,i)*gambet(i)
          pts(4,i) = pts(4,i)*gambet(i)
        end do
        
        end subroutine HardEdgeQuadFringeMap
        
      end module Quadrupoleclass
