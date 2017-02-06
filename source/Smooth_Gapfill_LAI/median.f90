subroutine median(xin,n,xmed)

! Find the median of X(1), ... , X(N), using as much of the quicksort
! algorithm as is needed to isolate it.


!     Latest revision - 26 November 1996
implicit none

integer :: n
double precision :: xin(n),x(n),xmed

! Local variables

integer          :: hi, lo, nby2, nby2p1, mid, i, j, k
double precision :: temp, xhi, xlo, xmax, xmin
logical          :: odd


x = xin

nby2 = n / 2
nby2p1 = nby2 + 1
odd = .true.

!     HI & LO are position limits encompassing the median.

if (n == 2 * nby2) odd = .false.
lo = 1
hi = n
if (n < 3) then
  if (n < 1) then
    xmed = 0.d0
    return
  end if
  xmed = x(1)
  if (n == 1) return
  xmed = 0.d5*(xmed + x(2))
  return
end if

!     Find median of 1st, middle & last values.

10 mid = (lo + hi)/2
xmed = x(mid)
xlo = x(lo)
xhi = x(hi)
if (xhi < xlo) then          ! Swap xhi & xlo
  temp = xhi
  xhi = xlo
  xlo = temp
end if
if (xmed > xhi) then
  xmed = xhi
else if (xmed < xlo) then
  xmed = xlo
end if

! The basic quicksort algorithm to move all values <= the sort key (XMED)
! to the left-hand end, and all higher values to the other end.

i = lo
j = hi
50 do
  if (x(i) >= xmed) exit
  i = i + 1
end do
do
  if (x(j) <= xmed) exit
  j = j - 1
end do
if (i < j) then
  temp = x(i)
  x(i) = x(j)
  x(j) = temp
  i = i + 1
  j = j - 1

!     Decide which half the median is in.

  if (i <= j) go to 50
end if

if (.not. odd) then
  if (j == nby2 .and. i == nby2p1) go to 130
  if (j < nby2) lo = i
  if (i > nby2p1) hi = j
  if (i /= j) go to 100
  if (i == nby2) lo = nby2
  if (j == nby2p1) hi = nby2p1
else
  if (j < nby2p1) lo = i
  if (i > nby2p1) hi = j
  if (i /= j) go to 100

! Test whether median has been isolated.

  if (i == nby2p1) return
end if
100 if (lo < hi - 1) go to 10

if (.not. odd) then
  xmed = 0.5d0*(x(nby2) + x(nby2p1))
  return
end if
temp = x(lo)
if (temp > x(hi)) then
  x(lo) = x(hi)
  x(hi) = temp
end if
xmed = x(nby2p1)
return

! Special case, N even, J = N/2 & I = J + 1, so the median is
! between the two halves of the series.   Find max. of the first
! half & min. of the second half, then average.

130 xmax = x(1)
do k = lo, j
  xmax = max(xmax, x(k))
end do
xmin = x(n)
do k = i, hi
  xmin = min(xmin, x(k))
end do
xmed = 0.5d0*(xmin + xmax)

end subroutine median
