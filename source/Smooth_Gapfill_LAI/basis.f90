!-------------------------------------------------------------------------------!
!                                                                               !                                                                      

 subroutine basis(npt,nptperyear,f,g,yag,ydl)                        

!                                                                               !           
! This subroutine computes sine and cosine basis functions                      !
! to be used in the subroutine SEASON. It also computes                         !
! polynomial basis functions to be used in the subroutine SAVGOL                !
!                                                                               !
!  Authors:                                                                     !
!  Per J\"onsson, Malm\"o University, Sweden                                    !
!  e-mail per.jonsson@ts.mah.se                                                 !
!                                                                               !
!  Lars Eklundh, Lund University, Sweden                                        !
!  e-mail lars.eklundh@nateko.lu.se                                             !
!                                                                               !
!  Updated 1/11-2005                                                            !
!                                                                               !
!-------------------------------------------------------------------------------!

implicit none

integer :: i,j,npt,nptperyear
double precision :: pi,width,fact(7),x(5),t(npt),f(npt,5),g(npt,3)
double precision :: yag(2*nptperyear,7,2),ydl(2*nptperyear,5,2),arg(2*nptperyear)

!---- Generate sine basis -------------------------------------------------------

pi = 2.d0*asin(1.d0)
t = dble((/(i,i=1,npt)/))/dble(nptperyear)

f(:,1) = 1.d0
f(:,2) = sin(2.d0*pi*t)
f(:,3) = cos(2.d0*pi*t)
f(:,4) = sin(4.d0*pi*t)
f(:,5) = cos(4.d0*pi*t)

!---- Generate polynomial basis -------------------------------------------------

g(:,1) = 1.d0
g(:,2) = dble((/(i,i=1,npt)/))
g(:,3) = dble((/(i,i=1,npt)/))**2.d0

!---- Generate asymmetric gaussian basis ---------------------------------------

width = nptperyear/4.d0

t = dble((/(i,i=1,npt)/))

fact(1:7) = (/ 0.2d0, 0.4d0, 0.7d0, 1.d0, 1.3d0, 1.7d0, 2.1d0 /)
x(1) = nptperyear + 0.5d0
x(3) = 3.d0
x(5) = 3.d0

do i = 1,7
  do j = 1,2
    x(2) = width*fact(i)/(j*log(2.d0)**0.25d0)
	x(4) = x(2)
	arg = (/ ((x(1)-t(1:nptperyear))/x(2))**x(3), &
	              ((t(nptperyear+1:2*nptperyear)-x(1))/x(4))**x(5) /)
    yag(:,i,j) = exp(-arg) 
  end do
end do

!---- Generate double logistic basis -------------------------------------------

fact(1:5) = (/ 0.3d0, 0.6d0, 1.d0, 1.5d0, 1.9d0 /)	
x(2) = nptperyear/20.d0
x(4) = x(2)

do i = 1,5
  do j = 1,2
    x(1) = nptperyear - width*fact(i)/j
	x(3) = nptperyear + width*fact(i)/j
	ydl(:,i,j) = 1.d0/(1.d0 + exp((x(1)-t(1:2*nptperyear))/x(2))) &
	                - 1.d0/(1.d0 + exp((x(3)-t(1:2*nptperyear))/x(4))) 
  end do
end do

end subroutine basis
   