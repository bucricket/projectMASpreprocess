!-------------------------------------------------------------------------------!
!                                                                               !

subroutine phenology(nfile,row,col,s,y,yfit,nlocal,startcutoff,minimum) 

!  This routine extracts phenological parameters from processed data            !
!                                                                               !
!   1.  time for which the left edge has increased to the 20 % level            !
!   2.  time for which the right edge has decreased to the 20 % level           !
!   3.  length of the season                                                    !
!   4.  average of left and right minimum values (base level)                   !
!   5.  time for the peak of the season                                         !
!   6.  value for the peak                                                      !
!   7.  seasonal amplitude                                                      !
!   8.  approximate left slope                                                  !
!   9.  approximate right slope                                                 !
!  10.  large seasonal integral                                                 !
!  11.  small seasonal integral                                                 !
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

integer :: row,col,nyear,nptperyear,npt,nenvi,printflag,debugflag,win(3),nlocal
integer :: i,j,n,n1,n2,s(4*nyear+2),nfile,minposleft,minposright
integer :: maxposcenter,pos(1),maxminsum,minimum
real*4 :: dummy(11)
double precision :: wfact,y(npt),yfit(npt),minvalleft,minvalright,maxvalcenter,amplitude
double precision :: baseval,poscenter,posright20,posleft20,posright80,posleft80
double precision :: seasonlength,leftderivative,rightderivative,largeintegral,smallintegral
double precision :: valleft20,valleft80,valright20,valright80,startcutoff

common /parameters/ nyear,nptperyear,npt,nenvi,printflag,debugflag,wfact,win 

n = nlocal + 2  

if (minimum == 0) then       ! the position s(2) corresponds to a maximum
  n1 = 4
  n2 = n-2
else                         ! the position s(2) corresponds to a minimum
  n1 = 3
  n2 = n-3 
end if

!---- Dump row, column and number of seasons on file ---------------------------
!    (n-4)/2 = number of seasons for which we have phenological metrics

write(nfile)  row,col,(n-4)/2

do i = n1,n2,2
  minvalleft = minval(yfit(s(i-2):s(i)))
  pos(1:1) = minloc(yfit(s(i-2):s(i)))
  minposleft = pos(1) + s(i-2) - 1
  
  minvalright = minval(yfit(s(i):s(i+2)))
  pos(1:1) = minloc(yfit(s(i):s(i+2)))
  minposright = pos(1) + s(i) - 1
  
  maxvalcenter = maxval(yfit(s(i-1):s(i+1)))
  pos(1:1) = maxloc(yfit(s(i-1):s(i+1)))
  maxposcenter = pos(1) + s(i-1) - 1
  
  amplitude = maxvalcenter - minvalleft
  
!---- Find left 80% level. Step from left min to max -------------------------- 
  
  valleft80 = minvalleft + 0.8d0*amplitude;
  posleft80 = 0.d0                   
  do j = minposleft,maxposcenter
    if (yfit(j) > valleft80) then
      posleft80 = j - 1.d0 + (valleft80 - yfit(j-1))/(yfit(j) - yfit(j-1))
!      posleft80 = j - 0.5d0 
      exit
    end if
  end do

!---- Find left X % level. Step from max to left min --------------------------   
  
  valleft20 = minvalleft + startcutoff*amplitude
  posleft20 = 0.d0
  do j = maxposcenter,minposleft,-1
    if (yfit(j) < valleft20) then
      posleft20 = j + (valleft20 - yfit(j))/(yfit(j+1) - yfit(j))
!      posleft20 = j + 0.5d0         
	  exit
    end if
  end do
  
  amplitude = maxvalcenter - minvalright
  
!---- Find right 80% level. Step from right min to max ------------------------
    
  valright80 = minvalright + 0.8d0*amplitude
  posright80 = 0.d0               
  do j = minposright,maxposcenter,-1
    if (yfit(j) > valright80) then
	  posright80 = j + 1.d0 - (valright80 - yfit(j+1))/(yfit(j) - yfit(j+1))  
!      posright80 = j + 0.5d0
	  exit
    end if
  end do
  
!---- Find right X % level. Step from max to right min ------------------------  
    
  valright20 = minvalright + startcutoff*amplitude
  posright20 = 0.d0
  do j = maxposcenter,minposright
    if (yfit(j) < valright20) then
      posright20 = j - (valright20 - yfit(j))/(yfit(j-1) - yfit(j));        
!      posright20 = j - 0.5d0        
	  exit
    end if
  end do

!--- Comment: We have added 0.01 in the denominator for the derivatives to ensure that
!    we never divide by zero 
  
  if (posleft20*posleft80*posright20*posright80 > 0.d0) then
    seasonlength = posright20 - posleft20                                   ! 3) season length  
    baseval = (minvalleft + minvalright)/2.d0                               ! 4) basevalue
    poscenter = (posleft80 + posright80)/2.d0                               ! 5) time for peak
    amplitude = maxvalcenter - baseval                                      ! 7) seasonal amplitude 
    leftderivative = (valleft80 - valleft20)/(0.001 + posleft80 - posleft20) ! 8) left derivative
    rightderivative = -(valright80 - valright20)/(-0.001 + posright80 - posright20)  ! 9) right derivative
    largeintegral = sum(yfit(floor(posleft20):ceiling(posright20)))         ! 10) large integral
    smallintegral = sum(yfit(floor(posleft20):ceiling(posright20))-baseval) ! 11) small integral  
  else 
    posleft20 = 0.d0         ! 1)
    posright20 = 0.d0        ! 2)
    seasonlength = 0.d0      ! 3)
    baseval = 0.d0           ! 4)
    poscenter = 0.d0         ! 5)
    maxvalcenter = 0.d0      ! 6)
    amplitude = 0.d0         ! 7)          
    leftderivative = 0.d0    ! 8)
    rightderivative = 0.d0   ! 9)
    largeintegral = 0.d0     ! 10)
    smallintegral = 0.d0     ! 11)
  end if
  dummy(1) = posleft20 
  dummy(2) = posright20 
  dummy(3) =   seasonlength 
  dummy(4) =   baseval 
  dummy(5) =   poscenter 
  dummy(6) =   maxvalcenter 
  dummy(7) =   amplitude       
  dummy(8) =   leftderivative 
  dummy(9) =   rightderivative 
  dummy(10) =   largeintegral      
  dummy(11) =   smallintegral      
  write(nfile) dummy
end do

end subroutine phenology

