!-------------------------------------------------------------------------------!
!                                                                               !                                      

subroutine savgol(y,w,f,y1,next)

!  Iterative and adaptive Sawitzky-Golay filering                               !
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

integer :: i,j,nyear,nptperyear,npt,nenvi,printflag,debugflag,win(3),next,m1,m2
integer :: failleft, failright
double precision :: wfact,y(npt),w(npt),y1(npt),yext(1:npt+2*next)  !! OBS next
double precision :: yfit(1:npt+2*next),wfit(1:npt+2*next),ymed
double precision :: yfitmean,yfitstd,ydiff,f(npt,3),c(3)

common /parameters/ nyear,nptperyear,npt,nenvi,printflag,debugflag,wfact,win 

!---- Extend data circularly to improve fitting near the boundary of the -------
!     original data                                                      

yext = (/y(npt-next+1:npt),y,y(1:next)/)
wfit = (/w(npt-next+1:npt),w,w(1:next)/)
yfit = yext

!---- Iterative and adaptive Savitzky-Golay filter -----------------------------

do j = 1,nenvi

!---- Compute mean and standard deviation -------------------------------------- 

  yfitmean = sum(yfit)/(npt+2.d0*next)
  yfitstd = sqrt(dot_product(yfit-yfitmean,yfit-yfitmean)/(npt+2.d0*next-1.d0))

  do i = 1+next,npt+next
      
!---- Set fitting window -------------------------------------------------------      
      
    m1 = i - win(j)   
    m2 = i + win(j) 
 
!---- Adapt fitting interval. Large variation use a smaller window -------------    
   
    if (maxval(yfit(m1:m2)) - minval(yfit(m1:m2)) > 1.2d0*2.d0*yfitstd) then 
      m1 = m1 + win(j)/3
      m2 = m2 - win(j)/3
    end if
    
!---- Check so that there are enough points, at least 3 at either side, with ---    
!     weights different from zero. If not, extend fitting window 

    failleft = 0
    do while ( (count(abs(wfit(m1:i)) > 1.d-10) < 3) .and. (failleft == 0) )    
      m1 = m1 - 1
	  if (m1 < 1) then 
        failleft = 1
        m1 = 1
      end if
    end do

    failright = 0
    do while ( (count(abs(wfit(i:m2)) > 1.d-10) < 3) .and. (failright == 0) )    
      m2 = m2 + 1
	  if (m2 > npt+2*next) then
	    failright = 1
		m2 = npt+2*next
	  end if
    end do
        
!---- Fit polynomial if enough data values with non-zero weight ----------------

    if ((failleft == 0) .and. (failright == 0)) then

!---- Fit polynomial of degree two to the datavalues in the interval [m1,m2] ---
!     Observe that the basis functions need to be shifted in order to get a  
!     better conditioning                                                                 

      call linlsq(m2-m1+1,3,f(1:m2-m1+1,:),yext(m1:m2),wfit(m1:m2),c)
     
 ! Compute the fitted function

      yfit(i) = c(1)*f(i-m1+1,1) + c(2)*f(i-m1+1,2) + c(3)*f(i-m1+1,3)

!---- If not enough values with non-zero weigth replace with median value ------

	else

	   call median(yext(m1:m2),m2-m1+1,ymed)
	   yfit(i) = ymed

    end if
    
  end do

!---- Modify weights -----------------------------------------------------------

  if (j < nenvi) then
    call modweight(j,npt+2*next,wfit,yext,yfit)
  end if
 
end do

!---- Transfer fitted points in the intervall [1,npt] to y1 --------------------

y1 = yfit(1+next:npt+next)


!---- Output for debug purposes ------------------------------------------------

! write(*,*) 1,sum(y1)

end subroutine savgol