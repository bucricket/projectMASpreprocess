!-------------------------------------------------------------------------------!
!                                                                               !

subroutine spike(y,w,spikecutoff,error)

!  Determine single spikes and set the corresponding weight to zero             !      
!                                                                               !
!  Criteria for spike:                                                          !
!                                                                               !
!  |y(i)-ymed| > cutoff .AND.                                                   !
!                                                                               !
!  (y(i) < (y(i-1)+y(i+1))/2 - cutoff .OR. y(i) > max(y(i-1),y(i+1)) + cutoff)  !
!                                                                               !
!                                                                               !
!  where ymed is the median filtered value                                      !
!                                                                               !
!                                                                               !
!                                                                               !
!  Authors:                                                                     !
!  Per J\"onsson, Malm\"o University, Sweden                                    !
!  e-mail per.jonsson@ts.mah.se                                                 !
!                                                                               !
!  Lars Eklundh, Lund University, Sweden                                        !
!  e-mail lars.eklundh@nateko.lu.se                                             !
!                                                                               !
!  Updated 12/9-2006                                                            !
!                                                                               !
!-------------------------------------------------------------------------------!

implicit none

integer :: i,j,k,m1,m2,winmax,nyear,nptperyear,npt,nenvi,printflag,debugflag
integer :: win(3),error
double precision :: y(npt),w(npt),yext(npt+2*nptperyear/5),wext(npt+2*nptperyear/5)
double precision :: ymean,ystd,ymed,spikecutoff,dist,ynzw(2*nptperyear/5+1),wfact

common /parameters/ nyear,nptperyear,npt,nenvi,printflag,debugflag,wfact,win 

ymean = sum(y)/dble(npt)
ystd = sqrt(dot_product(y-ymean,y-ymean)/dble(npt-1))
dist = spikecutoff*ystd

winmax = nptperyear/7  ! 5/7/10????????????

wext = (/ w(npt-winmax+1:npt),w,w(1:winmax) /)
yext = (/ y(npt-winmax+1:npt),y,y(1:winmax) /)

!---- Find single spikes by comparing with median filtered values and with -----
!     closest neighbors. If spike set corresponding weight to zero                           

do i = 1+winmax,npt+winmax
  m1 = i - winmax
  m2 = i + winmax
  
  k = 0
  do j = m1,m2
    if (wext(j) > 0.d0) then
	  k = k + 1
      ynzw(k) = yext(j)
    end if
  end do
  if (k > 0) then
!    write(*,*) 'error in spike'
!	error = 1
!	return
!  end if
    call median(ynzw(1:k),k,ymed)
    if ( (abs(y(i-winmax) - ymed) >= dist) .and. &
      ( (y(i-winmax) < (yext(i-1)+yext(i+1))/2.d0 - dist) .or. &
        (y(i-winmax) > max( yext(i-1),yext(i+1) ) + dist) ) ) then
      w(i-winmax) = 0.d0
    end if
  end if
end do

end subroutine spike