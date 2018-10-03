      subroutine hmltn1(h)
c  this subroutine is used to specify a constant hamiltonian
c Written by Alex Dragt, Spring 1987
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'extalk.inc'
      include 'hmflag.inc'
      dimension h(monoms)
c 
c  begin computation
      iflag=0
      do 10 i=1,monoms
   10 h(i)=-fa(i)
c     
      return
      end
c
      subroutine hmltn2(t,yy,h)
      write(6,*)'inside dummy routine hmltn2'
      call myexit
      end
c
      subroutine hmltn3(t,yy,h)
      write(6,*)'inside dummy routine hmltn3'
      call myexit
      end
c
      subroutine hmltn4(t,yy,h)
      write(6,*)'inside dummy routine hmltn4'
      call myexit
      end
c
      subroutine hmltn5(t,yy,h)
      write(6,*)'inside dummy routine hmltn5'
      call myexit
      end
c
