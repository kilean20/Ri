        subroutine Sector_Dipole(len,beta,h0,k1,ptarry1,Nplocal,qm0)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: Nplocal
        double precision, pointer, dimension(:,:) :: ptarry1
        double precision :: h0,len,beta,k1,qm0
        double precision, dimension(6) :: ptarry2
        double precision :: kx2,kx,cx,sx,dx,j1,ky2,ky,cy,sy,dy,&
                            gambet,gam2,qmrel
        integer :: i

        gambet = beta/sqrt(1.0-beta**2)
        gam2 = 1.0/(1.0-beta**2) 

        kx2 = h0**2 + k1
        if(kx2*len**2 .gt. 1.0d-14) then
          kx = dsqrt(kx2)
          cx = dcos(kx*len)
          dx = (1.0d0-cx)/kx2
          sx = dsin(kx*len)/kx 
          j1 = (len-sx)/kx2 
        else if(abs(kx2*len**2) .le. 1.0d-14) then
          kx = dsqrt(kx2)
          cx = dcos(kx*len)
          dx = len**2/2. 
          sx = len 
          j1 = len**3/6. 
        else
          kx = dsqrt(-kx2)
          cx = dcosh(kx*len)
          dx = (1.0d0-cx)/kx2
          sx = dsinh(kx*len)/kx 
          j1 = (len-sx)/kx2 
        endif

        ky2 = -k1
        if(ky2*len**2 .gt. 1.0d-14) then
          ky = dsqrt(ky2)
          cy = dcos(ky*len)
          dy = (1.0d0-cy)/ky2
          sy = dsin(ky*len)/ky
        else if(abs(ky2*len**2) .le. 1.0d-14) then
          ky = dsqrt(ky2)
          cy = dcos(ky*len)
          dy = len**2/2. 
          sy = len
        else
          ky = dsqrt(-ky2)
          cy = dcosh(ky*len)
          dy = (1.0d0-cy)/ky2
          sy = dsinh(ky*len)/ky
        endif
!        print*,"h0,k1: ",h0,k1,beta,len,kx,cx,dx,sx,ky,cy,dy,sy

        ptarry2 = 0.0d0
        do i = 1, Nplocal
          ptarry2(1) = cx*ptarry1(1,i)+sx*ptarry1(2,i)+&
                     h0*dx*ptarry1(6,i) 
          ptarry2(2) = -kx2*sx*ptarry1(1,i)+cx*ptarry1(2,i)+&
                     h0*sx*ptarry1(6,i) 
          ptarry2(3) = cy*ptarry1(3,i)+sy*ptarry1(4,i)
          ptarry2(4) = -ky2*sy*ptarry1(3,i)+cy*ptarry1(4,i)
! (5) is defined as -v dt = -c beta dt
! see p.147 of F. Iselin paper
          qmrel = (ptarry1(7,i)-qm0)/qm0
          ptarry2(5) = -h0*sx*ptarry1(1,i)-h0*dx*ptarry1(2,i)+&
                     ptarry1(5,i) - h0**2*j1*ptarry1(6,i) + &
                     len/gam2*(ptarry1(6,i)+qmrel)
! Transport defines (5) as path differnce which is v dt
!          ptarry2(5) = h0*sx*ptarry1(1,i)+h0*dx*ptarry1(2,i)+&
!                     ptarry1(5,i) + h0**2*j1*ptarry1(6,i)
          ptarry2(6) = ptarry1(6,i)
          ptarry1(1,i) = ptarry2(1)
          ptarry1(2,i) = ptarry2(2)
          ptarry1(3,i) = ptarry2(3)
          ptarry1(4,i) = ptarry2(4)
          ptarry1(5,i) = ptarry2(5)
          ptarry1(6,i) = ptarry2(6)
        enddo

        end subroutine Sector_Dipole

