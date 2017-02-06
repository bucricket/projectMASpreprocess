!-------------------------------------------------------------------------------!
!                                                                               !

subroutine season(y,w,f,s,nlocal,seasonpar,amplitudecutoff,minimum)

!  Detrend the time-series by a running average                                 !
!                                                                               !
!  Fit sine and cosine functions to determine the number of annual seasons      !
!  The current version of the code only allows one or two annual seasons        !
!                                                                               !
!  Determine points that define the fitting intervals for local functions       !  
!                                                                               !
!  Authors:                                                                     !
!  Per J\"onsson, Malm\"o University, Sweden                                    !
!  e-mail per.jonsson@ts.mah.se                                                 !
!                                                                               !
!  Lars Eklundh, Lund University, Sweden                                        !
!  e-mail lars.eklundh@nateko.lu.se                                             !
!                                                                               !
!  Updated May 2006                                                             !
!                                                                               !
!-------------------------------------------------------------------------------!

implicit none

integer :: i,j,nyear,nptperyear,npt,nseason,nenvi,printflag,debugflag,win(3),m1,m2
integer :: s(4*nyear+2),nlocal,yextpos(4),posmin1(1),posmin2(1),posmax1(1)
integer :: posmax2(1),yminpos,minimum
double precision :: y(npt),w(npt),ydt(npt),yav(npt),yfit(npt),f(npt,5),c(5)
double precision :: wfact,amp,yext(4), seasonpar, amplitudecutoff

common /parameters/ nyear,nptperyear,npt,nenvi,printflag,debugflag,wfact,win 

!---- Use running average with large window to remove trend in time-series -----

do i = 1,npt
    
!---- Set window ---------------------------------------------------------------     
    
  m1 = max(1,i - nptperyear/2)
  m2 = min(npt,i + nptperyear/2)
  
!---- Compute average in the window --------------------------------------------  
  
  yav(i) = sum(w(m1:m2)*y(m1:m2))/sum(w(m1:m2))
  
end do

!---- Remove the trend ---------------------------------------------------------

ydt = y - yav

!---- Fit sine and cosine functions to detrended time-series -------------------

call linlsq(npt,5,f,ydt,w,c)

!---- Evaluate the fitted function ---------------------------------------------

yfit = c(1)*f(:,1)+c(2)*f(:,2)+c(3)*f(:,3)+c(4)*f(:,4)+c(5)*f(:,5) 
        
!---- Check if the amplitude is larger than the criteria -----------------------

if ( (maxval(yfit) - minval(yfit)) <= amplitudecutoff ) then
  nlocal = 0
  s = 0
  return
end if

!---- Determine points that define the fitting intervals for local functions ---  

j = 0
if ( (yfit(1)-yfit(nptperyear))*(yfit(2)-yfit(1)) < 0.d0 ) then
  j = j + 1
  yext(j) = yfit(1)
  yextpos(j) = 1
end if

do i = 2,nptperyear-1
  if ( (yfit(i)-yfit(i-1))*(yfit(i+1)-yfit(i)) < 0.d0 ) then
    j = j + 1
    yext(j) = yfit(i)
    yextpos(j) = i
  end if
end do

if ( (yfit(nptperyear)-yfit(nptperyear-1))*(yfit(1)-yfit(nptperyear)) < 0.d0 ) then
  j = j + 1
  yext(j) = yfit(nptperyear)
  yextpos(j) = nptperyear
end if

if (j == 2) then
  s(1:2*nyear) = (/ ((i-1)*nptperyear+yextpos(1:2), i = 1,nyear) /)
  s(1:2*nyear+2) = (/ 1,s(1:2*nyear),npt /)
  nlocal = 2*nyear 
else
  posmax1(1:1) = maxloc(yext)
  posmax2(1:1) = maxloc(yext,mask=yext /= yext(posmax1(1)))
  posmin1(1:1) = minloc(yext)
  posmin2(1:1) = minloc(yext,mask=yext /= yext(posmin1(1))) 
  
  if ( (yext(posmax2(1)) - yext(posmin2(1)))/(yext(posmax1(1)) - yext(posmin1(1))) > seasonpar ) then  

!---- Two seasons ------------------------------------------------------------

    s(1:4*nyear) = (/ ((i-1)*nptperyear+yextpos(1:4), i = 1,nyear) /)
    s(1:4*nyear+2) = (/ 1,s(1:4*nyear),npt /)
    nlocal = 4*nyear 

  else

!---- One season -------------------------------------------------------------

!---- Fit single sine and cosine functions to detrended time-series ----------

!   write(*,*) sum(ydt),sum(w)
    call linlsq(npt,3,f(:,1:3),ydt,w,c)

!---- Evaluate the fitted function -------------------------------------------

    yfit = c(1)*f(:,1)+c(2)*f(:,2)+c(3)*f(:,3) 

    j = 0
    if ( (yfit(1)-yfit(nptperyear))*(yfit(2)-yfit(1)) < 0.d0 ) then
      j = j + 1
      yext(j) = yfit(1)
      yextpos(j) = 1
    end if

    do i = 2,nptperyear-1
      if ( (yfit(i)-yfit(i-1))*(yfit(i+1)-yfit(i)) < 0.d0 ) then
        j = j + 1
        yext(j) = yfit(i)
        yextpos(j) = i
      end if
    end do

    if ( (yfit(nptperyear)-yfit(nptperyear-1))*(yfit(1)-yfit(nptperyear)) < 0.d0 ) then
      j = j + 1
      yext(j) = yfit(nptperyear)
      yextpos(j) = nptperyear
    end if
    
	s(1:2*nyear) = (/ ((i-1)*nptperyear+yextpos(1:2), i = 1,nyear) /)
    s(1:2*nyear+2) = (/ 1,s(1:2*nyear),npt /)
    nlocal = 2*nyear 

  end if
end if

!---- Find out if position corresponds to a minimum or a maximum. This information is ---
!     used in PHENOLOGY

if ( yfit(s(2)) < yfit(s(3)) ) then
  minimum = 1
else
  minimum = 0
end if

end subroutine season   