!c--- This file is from hypoDD by Felix Waldhauser ---------
!c-------------------------Modified by Haijiang Zhang-------
!c Multiply a matrix by a vector
!c  Version for use with sparse matrix specified by
!c  output of subroutine sparse for use with LSQR

	subroutine aprod(mode, m, n, x, y, leniw, lenrw, iw, rw)

	implicit none

!c	Parameters:
	integer	mode		! ==1: Compute  y = y + a*x
				! 	y is altered without changing x
				! ==2: Compute  x = x + a(transpose)*y
				!	x is altered without changing y
	integer	m, n		! Row and column dimensions of a
	real	x(n), y(m)	! Input vectors
	integer(kind=8) :: leniw
	integer(kind=8) :: lenrw
	integer	iw(leniw)	! Integer work vector containing:
				! iw[1]  Number of non-zero elements in a
				! iw[1] will be set seperately: Modified by Xin
				! lenrw*2 will replace iw[1]
				! iw[1:lenrw]  Row indices of non-zero elements
				! iw[lenrw+1:2*lenrw]  Column indices
	real	rw(lenrw)	! Non-zero elements of a

!c	Local variables:
	integer(kind=8) i1
	integer(kind=8) j1
	integer(kind=8) k
	integer(kind=8) kk

!c	set the ranges the indices in vector iw

	kk=lenrw
	! i1=1
	! j1=kk+1
	i1=0
	j1=kk

!c	main iteration loop

	do k = 1,kk
	
	if (mode.eq.1) then

!c	compute  y = y + a*x

	y(iw(i1+k)) = y(iw(i1+k)) + rw(k)*x(iw(j1+k))

	else

!c	compute  x = x + a(transpose)*y

	x(iw(j1+k)) = x(iw(j1+k)) + rw(k)*y(iw(i1+k))

	endif
        enddo

!  100	continue

	return
	end

	subroutine aprod_dp(mode, m, n, x, y, leniw, lenrw, iw, rw)
	use lsmrDataModule, only : dp
	implicit none

	integer	mode
	integer	m, n
	real(dp) x(n), y(m)
	integer(kind=8) :: leniw
	integer(kind=8) :: lenrw
	integer	iw(leniw)
	real	rw(lenrw)

	integer(kind=8) i1
	integer(kind=8) j1
	integer(kind=8) k
	integer(kind=8) kk

	kk=lenrw
	i1=0
	j1=kk

	do k = 1,kk
	if (mode.eq.1) then
	y(iw(i1+k)) = y(iw(i1+k)) + real(rw(k),dp)*x(iw(j1+k))
	else
	x(iw(j1+k)) = x(iw(j1+k)) + real(rw(k),dp)*y(iw(i1+k))
	endif
	enddo

	return
	end
