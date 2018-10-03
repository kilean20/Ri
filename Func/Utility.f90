!----------------------------------------------------------------
! (c) Copyright, 2001 by the Regents of the University of California.
! Utilityclass: 2D array (matrix) manipulation utilities in Linear  
!                 Algebra module of FUNCTION layer.
! Version: 1.0
! Author: C. Mitchell, LBNL, 11/30/16
! Description: This class defines a utility class which contains
!              matrix manipulation tools including inversion of
!              a symmetric, positive definite matrix.
! Comments: The following are standard utility tools taken from
!           the package LINPACK, which is distributed under the
!           GNU LGPL license.
!----------------------------------------------------------------
      module Utilityclass

      contains

subroutine dpodi ( a, lda, n, det, job )

!*****************************************************************************80
!
!! DPODI computes the determinant and inverse of a certain matrix.
!
!  Discussion:
!
!    The matrix is real symmetric positive definite.
!    DPODI uses the factors computed by DPOCO, DPOFA or DQRDC.
!
!    A division by zero will occur if the input factor contains
!    a zero on the diagonal and the inverse is requested.
!    It will not occur if the subroutines are called correctly
!    and if DPOCO or DPOFA has set INFO == 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(LDA,N).  On input, the output A from
!    DPOCO or DPOFA, or the output X from DQRDC.  On output, if DPOCO or
!    DPOFA was used to factor A then DPODI produces the upper half of
!    inverse(A).  If DQRDC was used to decompose X then DPODI produces
!    the upper half of inverse(X'*X) where X' is the transpose.
!    Elements of A below the diagonal are unchanged.  If the units digit
!    of JOB is zero, A is unchanged.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input, integer ( kind = 4 ) JOB, specifies the task.
!    11, both determinant and inverse.
!    01, inverse only.
!    10, determinant only.
!
!    Output, real ( kind = 8 ) DET(2), the determinant of A or of X'*X
!    if requested.
!      determinant = DET(1) * 10.0^DET(2)
!    with 1.0D+00 <= DET(1) < 10.0D+00 or DET(1) == 0.0D+00.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  real ( kind = 8 ) t
!
!  Compute the determinant.
!
  if ( job / 10 /= 0 ) then

    det(1) = 1.0D+00
    det(2) = 0.0D+00

    do i = 1, n

      det(1) = det(1) * a(i,i) * a(i,i)

      if ( det(1) == 0.0D+00 ) then
        exit
      end if

      do while ( det(1) < 1.0D+00 )
        det(1) = det(1) * 10.0D+00
        det(2) = det(2) - 1.0D+00
      end do

      do while ( 10.0D+00 <= det(1) )
        det(1) = det(1) / 10.0D+00
        det(2) = det(2) + 1.0D+00
      end do

    end do

  end if
!
!  Compute inverse(R).
!
  if ( mod ( job, 10 ) /= 0 ) then

    do k = 1, n

      a(k,k) = 1.0D+00 / a(k,k)
      t = -a(k,k)
      call dscal ( k-1, t, a(1,k), 1 )

      do j = k + 1, n
        t = a(k,j)
        a(k,j) = 0.0D+00
        call daxpy ( k, t, a(1,k), 1, a(1,j), 1 )
      end do

    end do
!
!  Form inverse(R) * (inverse(R))'.
!
    do j = 1, n
      do k = 1, j - 1
        t = a(k,j)
        call daxpy ( k, t, a(1,j), 1, a(1,k), 1 )
      end do
      t = a(j,j)
      call dscal ( j, t, a(1,j), 1 )
    end do

  end if

  return
end

subroutine dpofa ( a, lda, n, info )

!*****************************************************************************80
!
!! DPOFA factors a real symmetric positive definite matrix.
!
!  Discussion:
!
!    DPOFA is usually called by DPOCO, but it can be called
!    directly with a saving in time if RCOND is not needed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 April 2013
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(LDA,N).  On input, the symmetric matrix
!    to be  factored.  Only the diagonal and upper triangle are used.
!    On output, an upper triangular matrix R so that A = R'*R
!    where R' is the transpose.  The strict lower triangle is unaltered.
!    If INFO /= 0, the factorization is not complete.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    0, for normal return.
!    K, signals an error condition.  The leading minor of order K is not
!    positive definite.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) s
  real ( kind = 8 ) t

  do j = 1, n

    s = 0.0D+00

    do k = 1, j - 1
      t = a(k,j) - dot_product ( a(1:k-1,k), a(1:k-1,j) )
      t = t / a(k,k)
      a(k,j) = t
      s = s + t * t
    end do

    s = a(j,j) - s

    if ( s <= 0.0D+00 ) then
      info = j
      return
    end if

    a(j,j) = sqrt ( s )

  end do

  info = 0

  return
end

      end module Utilityclass
