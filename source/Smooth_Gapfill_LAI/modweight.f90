!-------------------------------------------------------------------------------!
!                                                                               !

subroutine modweight(j,n,wfit,y,yfit)

!                                                                               !
!  Authors:                                                                     !
!  Per J\"onsson, Malm\"o University, Sweden                                    !
!  e-mail per.jonsson@ts.mah.se                                                 !
!                                                                               !
!  Lars Eklundh, Lund University, Sweden                                        !
!  e-mail lars.eklundh@nateko.lu.se                                             !
!                                                                               !
!  Updated 26/3-2005                                                            !
!                                                                               !
!-------------------------------------------------------------------------------!

implicit none

integer :: i,j,m1,m2,nyear,nptperyear,npt,n,nenvi,printflag,debugflag,win(3)
double precision :: wfit(n),y(n),yfit(n),yfitmean,yfitstd,ydiff,wfact

common /parameters/ nyear,nptperyear,npt,nenvi,printflag,debugflag,wfact,win 

!---- Compute mean and standard deviation -------------------------------------- 

yfitmean = sum(yfit)/n
yfitstd = sqrt(dot_product(yfit-yfitmean,yfit-yfitmean)/(n-1.d0))

  
!---- Adjust the weights dependent on if the values are above or below the -----
!     fitted values                                                         

do i = 1,n
  m1 = max(1,i - nptperyear/7)  
  m2 = min(n,i + nptperyear/7)
  if (y(i) < yfit(i) - 1.d-8) then    
          
!---- Modify weights of all points in the first iteration. After that ----------
!     only weights of high points are modified                        
          
    if ((minval(yfit(m1:m2)) > yfitmean) .or. (j < 2)) then

!---- If there is a low variation in an interval, i.e. if the interval ---------
!     is at a peak or at a minima compute the normalized distance      
!     between the data point and the fitted point                      

      if (maxval(yfit(m1:m2))-minval(yfit(m1:m2)) < 0.4d0*2.d0*yfitstd) then
        ydiff = 2.d0*(yfit(i)-y(i))/yfitstd  
      else
        ydiff = 0.d0
      end if

!---- Use the computed distance to modify the weight. Large distance -----------
!     will give a small weight                                       
          
      wfit(i) = wfact*wfit(i)*exp(-ydiff**2.d0)
	  
    end if      
  end if
end do

end subroutine modweight