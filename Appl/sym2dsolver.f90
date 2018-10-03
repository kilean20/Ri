!innx - mode number in x direction
!inny - mode number in y direction
!xaper - x aperture/2
!yaper - y aperture/2
        subroutine sym2dsolver_BeamBunch(rays,perv,tau,innp,Npt,innx,inny,xaper,yaper)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp, innx, inny, Npt
        double precision, intent (inout), dimension (9,innp) :: rays
        real*8, dimension(innx,inny) :: philm,xden
        real*8, dimension(innx,innp) :: tmpl
        real*8, dimension(inny,innp) :: tmpm
        real*8, intent(in) :: xaper,yaper,tau
        real*8, dimension(innx,innp) :: tmpcl
        real*8, dimension(inny,innp) :: tmpcm
        real*8 :: perv
        integer :: i, j, ip,Nl,Nm,ierr
        real*8 :: aa,bb,xmin,xmax,ymin,ymax,alphal,betam,pi,egx,egy,xx,yy
        real*8 :: gamlm2,ForceX,ForceY

        Nl = innx
        Nm = inny

        pi = 2*asin(1.0d0)

        xmin = -xaper
        xmax = xaper
        ymin = -yaper
        ymax = yaper

        aa = xmax-xmin
        bb = ymax-ymin

        do ip = 1, innp
          xx = rays(1,ip)-xmin 
          do i = 1, Nl
            alphal = pi*i/aa
            tmpl(i,ip) = sin(alphal*xx)
          enddo
        enddo
        do ip = 1, innp
          yy = rays(3,ip)-ymin
          do j =1, Nm
            betam = pi*j/bb
            tmpm(j,ip) = sin(betam*yy)
          enddo
        enddo

        do ip = 1, innp
          xx = rays(1,ip) - xmin
          do i = 1, innx
            alphal = pi*i/aa
            tmpcl(i,ip) = cos(alphal*xx)
          enddo
        enddo
        do ip = 1, innp
          yy = rays(3,ip) - ymin
          do j =1, inny
            betam = pi*j/bb
            tmpcm(j,ip) = cos(betam*yy)
          enddo
        enddo

        do j =1, Nm
          do i = 1, Nl
            philm(i,j) = 0.0d0
            do ip = 1, innp
              philm(i,j) = philm(i,j) + tmpm(j,ip)*tmpl(i,ip)
            enddo
          enddo
        enddo

        xden = 0.0d0
        call MPI_ALLREDUCE(philm,xden,Nl*Nm,MPI_DOUBLE_PRECISION,MPI_SUM,&
                                  MPI_COMM_WORLD,ierr)

        philm = xden*4/aa/bb/dble(Npt)

       !here, q/epsilon0 is missing
        !we need to check that.
        do j =1, Nm
          betam = pi*j/bb
          do i = 1, Nl
            alphal = pi*i/aa
            gamlm2 = alphal**2+betam**2
            philm(i,j) = philm(i,j)/gamlm2
          enddo
        enddo

        philm = philm*4*pi

        do ip = 1, innp
          egx = 0.0d0
          egy = 0.0d0
          do j = 1, inny
            betam = pi*j/bb
            do i = 1, innx
              alphal = pi*i/aa
              egx = egx - philm(i,j)*alphal*tmpcl(i,ip)*tmpm(j,ip)
              egy = egy - philm(i,j)*betam*tmpl(i,ip)*tmpcm(j,ip)
            enddo
          enddo
          ForceX = 0.5d0 * egx * perv
          ForceY = 0.5d0 * egy * perv
        !Here, we apply the momentum update:
          rays(2,ip) = rays(2,ip) + tau*ForceX
          rays(4,ip) = rays(4,ip) + tau*ForceY
        enddo

        end subroutine sym2dsolver_BeamBunch


!innx - mode number in x direction
!inny - mode number in y direction
!xaper - x aperture/2
!yaper - y aperture/2
        subroutine sym2dDiag_BeamBunch(rays,perv,tau,innp,Npt,innx,inny,xaper,yaper,Htotal)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp, innx, inny, Npt
        double precision, intent (in), dimension (9,innp) :: rays
        double precision, intent (out) :: Htotal
        real*8, dimension(innx,inny) :: philm,xden
        real*8, dimension(innx,innp) :: tmpl
        real*8, dimension(inny,innp) :: tmpm
        real*8, intent(in) :: xaper,yaper,tau
        real*8, dimension(innx,innp) :: tmpcl
        real*8, dimension(inny,innp) :: tmpcm
        real*8, dimension(innp) :: Hdiag
        real*8 :: perv
        integer :: i, j, ip,Nl,Nm,ierr,my_rank
        real*8 :: aa,bb,xmin,xmax,ymin,ymax,alphal,betam,pi,egx,egy,xx,yy
        real*8 :: gamlm2,ForceX,ForceY,H,Hlcl

        Nl = innx
        Nm = inny

        pi = 2*asin(1.0d0)

        xmin = -xaper
        xmax = xaper
        ymin = -yaper
        ymax = yaper

        aa = xmax-xmin
        bb = ymax-ymin

        do ip = 1, innp
          xx = rays(1,ip)-xmin 
          do i = 1, Nl
            alphal = pi*i/aa
            tmpl(i,ip) = sin(alphal*xx)
          enddo
        enddo
        do ip = 1, innp
          yy = rays(3,ip)-ymin
          do j =1, Nm
            betam = pi*j/bb
            tmpm(j,ip) = sin(betam*yy)
          enddo
        enddo

        do j =1, Nm
          do i = 1, Nl
            philm(i,j) = 0.0d0
            do ip = 1, innp
              philm(i,j) = philm(i,j) + tmpm(j,ip)*tmpl(i,ip)
            enddo
          enddo
        enddo

        xden = 0.0d0
        call MPI_ALLREDUCE(philm,xden,Nl*Nm,MPI_DOUBLE_PRECISION,MPI_SUM,&
                                  MPI_COMM_WORLD,ierr)

        philm = xden*4/aa/bb/dble(Npt)

       !here, q/epsilon0 is missing
        !we need to check that.
        do j =1, Nm
          betam = pi*j/bb
          do i = 1, Nl
            alphal = pi*i/aa
            gamlm2 = alphal**2+betam**2
            philm(i,j) = philm(i,j)/gamlm2
          enddo
        enddo

        philm = philm*4*pi

        Hdiag = 0.0d0
        do ip = 1, innp
          H = 0.0d0
          do j = 1, inny
            do i = 1, innx
                H = H + 0.5d0*philm(i,j)*tmpl(i,ip)*tmpm(j,ip)
            enddo
          enddo
          H = 0.5d0 * H * perv
        !Here, we apply the momentum update:
          Hdiag(ip)  = (rays(2,ip)**2+rays(4,ip)**2)/2.0d0 + H
        enddo
        Hlcl = sum(Hdiag)
        Htotal = 0.0d0
        call MPI_ALLREDUCE(Hlcl,Htotal,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
                                  MPI_COMM_WORLD,ierr)

        Htotal = Htotal/dble(Npt)

        end subroutine sym2dDiag_BeamBunch

