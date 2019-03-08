
    module SpaceChargeSF
!*********************************************************************
!The following module contains subroutines required to track particles
!in the analytically-provided space charge potential on the rectangular
!domain [-xmax,xmax] by [-ymax,ymax], as provided by the 2D spectral
!mode coefficients provided in ''.
!*********************************************************************
    implicit none
!  Numerical parameters:
    double precision, parameter:: xmax = 1.5d0  !half x-domain
    double precision, parameter:: ymax = 1.5d0  !half y-domain
    integer, parameter:: lmax = 15  !number of horizontal modes (nom 15)
    integer, parameter:: mmax = 15  !number of vertical modes (nom 15)
    integer, parameter:: nmode = lmax*mmax  !number of modes (lmax*mmax) (nom 225)
!  Storage of the array of mode coefficients:
    double precision:: SCpotential(nmode)

    contains


    subroutine initializeSpaceChargeSF
!*****************************************************************
!The following subroutine initializes the array of spectral mode
!coefficients required to evaluate the analytically-provided space
!charge potential on a rectangular domain.
!*****************************************************************
    implicit none
    integer:: j,l,m
    open(unit=10,file='SFcoefficients.dat',status='old')
    do j=1,nmode
       read(10,*) l,m,SCpotential(j)
    enddo
    close(10)
    end subroutine initializeSpaceChargeSF


    subroutine SpaceChargeSmoothFocusingPropagator(ksc,cnll,coord)
!*****************************************************************
!The following subroutine computes the nonlinear momentum kick
!associated with a single step in the 2D space charge potential
!provided in analytical form by the set of 2D Fourier coefficients
!stored in 'arr'.  Analytical derivatives of the potential are
!returned by the subroutine 'duapprox'.
!The arguments are as follows:
!         ksc - integrated strength of the kick (m)            
!         cnll - dimensional parameter of the lens (m)
!         coord = (x [m], px/p0, y [m], py/p0)
!C. Mitchell, 2/12/2019.
!*****************************************************************
    implicit none
    include 'mpif.h'
    double precision, intent(in) :: ksc,cnll
    double precision, dimension(4), intent(inout) :: coord
    double precision :: x,y,kick,dudx,dudy
    x = coord(1)/cnll                       !Dimensionless horizontal coord
    y = coord(3)/cnll                       !Dimensionless vertical coord
    kick = ksc/cnll                         !Dimensionless kick strength
  !Evaluate derivatives of the SC potential at (x,y)
    call duapprox(nmode,SCpotential,xmax,ymax,x,y,dudx,dudy)
  !Momentum update
    coord(2)=coord(2)-kick*dudx
    coord(4)=coord(4)-kick*dudy
    end subroutine SpaceChargeSmoothFocusingPropagator


    function eigenvalue(a,b,n)
!***************************************
! This function returns the eigenvalue
! indexed by n of the Laplacian in 
! the 2D domain given by the rectangle
! [-a,a]x[-b,b].
!***************************************
    implicit none
    double precision:: a,b,eigenvalue,pi
    double precision:: lx,ly
    integer:: l,m,n
    pi = 4.0d0*atan(1.0d0)
    call index2D(n,l,m)
    lx = -(l*pi/(2.0d0*a))**2
    ly = -(m*pi/(2.0d0*b))**2
    eigenvalue = lx+ly
    end function eigenvalue


    function eigenmode(a,b,n,x,y)
!*********************************************
! This function returns the eigenmode indexed
! by n of the Laplacian in the 2D domain
! given by the rectangle [-a,a]x[-b,b].
!*********************************************
    implicit none
    double precision:: a,b,x,y,eigenmode,pi
    double precision:: ex,ey
    integer:: l,m,n
    pi = 4.0d0*atan(1.0d0)
    call index2D(n,l,m)
    ex = dsin(l*pi*(x+a)/(2.0d0*a))/dsqrt(a)
    ey = dsin(m*pi*(y+b)/(2.0d0*b))/dsqrt(b)
    eigenmode = ex*ey
    end function eigenmode    


    subroutine deigenmode(a,b,n,x,y,dedx,dedy)
!*********************************************
! This subroutine returns the two partial
! derivatives (dedx,dedy) of the eigenmode en 
! indexed by n of the Laplacian in the 2D 
! domain given by the rectangle [-a,a]x[-b,b].
!*********************************************
    implicit none
    double precision:: a,b,x,y,dedx,dedy,pi
    double precision:: ex,ey,dex,dey,kx,ky
    integer:: l,m,n
    pi = 4.0d0*atan(1.0d0)
    call index2D(n,l,m)
    kx = l*pi/(2.0d0*a)
    ky = m*pi/(2.0d0*b)
    ex = dsin(kx*(x+a))/dsqrt(a)
    ey = dsin(ky*(y+b))/dsqrt(b)
    dex = kx*dcos(kx*(x+a))/dsqrt(a)
    dey = ky*dcos(ky*(y+b))/dsqrt(b)
    dedx = dex*ey
    dedy = ex*dey
    end subroutine deigenmode


    function uapprox(nmode,arr,a,b,x,y)
!***********************************************
! This subroutine evaluates the approximate
! solution u using nmode spectral modes, whose
! coefficients are stored in the array arr,
! at the point (x,y).
!**********************************************
    implicit none
    integer:: nmode,j
    double precision:: x,y,arr(nmode),term
    double precision:: a,b,uapprox
    uapprox = 0.0d0
    do j=1,nmode
       term = arr(j)*eigenmode(a,b,j,x,y)
       uapprox = uapprox + term
    enddo
    end function uapprox



    subroutine duapprox(nmode,arr,a,b,x,y,dudx,dudy)
!***********************************************
! This subroutine evaluates the two partial
! derivatives (dudx,dudy) of the solution u using
! nmode spectral modes, whose coefficients are 
! stored in the array arr, at the point (x,y).
!**********************************************
    implicit none
    integer:: nmode,j
    double precision:: x,y,arr(nmode),term
    double precision:: a,b,dudx,dudy
    double precision:: dedx,dedy
    dudx = 0.0d0
    dudy = 0.0d0
    do j=1,nmode
       call deigenmode(a,b,j,x,y,dedx,dedy) 
       dudx = dudx + arr(j)*dedx
       dudy = dudy + arr(j)*dedy  
    enddo
    end subroutine duapprox


    subroutine index(l,m,n)
!************************************************
! Given the 2D index (l,m), this function returns
! the single index n.
!***********************************************
    implicit none
    integer:: l,m,n
    n = (l-1)*mmax + m
    end subroutine index


    subroutine index2D(n,l,m)
!************************************************
! Given the single index n, this function returns
! the 2D index (l,m).
!************************************************
    implicit none
    integer:: l,m,n
    l = idint(dble(n-1)/mmax) + 1
    m = n - (l-1)*mmax
    end subroutine index2D


    end module SpaceChargeSF

