!-------------------------------------------------------------------------------!
!                                                                               !

subroutine fitlogistic(s,ydl,y,w,nlocal,yfit)

!  (1) Fit local double logistic functions                                      !
!                                                                               !
!  y(t) = c(1) + c(2)*(1/(1+exp((x(1)-t)/x(2))) - 1/(1+exp((x(3)-t)/x(4)))      !
!                                                                               !
!  in an iterative proceedure to approach the upper/lower envelope of           !
!  the data values. The linear, c(1) and c(2), and non-linear x(1),...,x(4)     !
!  parameters are determined using Marquard's method. Note that we separate     !
!  the two parameter sets in order to improve performance. The non-linear       !
!  parameters are box constrained bl(i) <= x(i) <= bu(i) to ensure that we      !
!  find a phenologically sound solution.                                        !  
!                                                                               !
!  (2) Merge the local functions to a global function                           !
!                                                                               !
!  Authors:                                                                     !
!  Per J\"onsson, Malm\"o University, Sweden                                    !
!  e-mail per.jonsson@ts.mah.se                                                 !
!                                                                               !
!  Lars Eklundh, Lund University, Sweden                                        !
!  e-mail lars.eklundh@nateko.lu.se                                             !
!                                                                               !
!  Improved estimates of starting parameters                                    !
!                                                                               !
!  Updated September 2006                                                       !
!                                                                               !
!-------------------------------------------------------------------------------!

implicit none

integer :: i,j,k,l,nmid,nyear,nptperyear,npt,nenvi,printflag,debugflag,win(3)
integer :: m0,m1,m2,n,n1,n2,nlocal,s(4*nyear+2),term,kmin,lmin,nseason
integer :: sumnonconstant
double precision :: wfact,y(npt),w(npt),yfit(npt),wfit(npt),ylocal(nlocal,npt)
double precision :: width,x(4),xn(4),ones(npt),t(npt),cut,ndiff,mid
double precision :: exp12(npt),ytilde(npt),Btilde(npt,2),c(2),bl(4),bu(4)
double precision :: ydl(2*nptperyear,5,2),fact(5),chi2,chi2min,ymed

common /parameters/ nyear,nptperyear,npt,nenvi,printflag,debugflag,wfact,win 

yfit = y
wfit = w
ylocal = 0.d0
ones = 1.d0

!---- Initiate matrix that will contain the local functions ---------------------  

ylocal = 0.d0

!---- Define widths for one and two seasons respectively ------------------------

if (nlocal/nyear == 2) then
  width = nptperyear/4.d0
  nseason = 1
else 
  width = nptperyear/8.d0
  nseason = 2
end if

fact = (/ 0.3d0, 0.6d0, 1.d0, 1.5d0, 1.9d0 /)

!---- Iterative fits to local double logistic functions in all intervalls -------

do j = 1,nenvi
  do i = 1,nlocal     
      
!---- Fitting intervall given by [m0,m2] ----------------------------------------      
      
    m0 = s(i)
    m1 = s(i+1)
    m2 = s(i+2)
	n = m2 - m0 + 1

!---- Check if the are enought non constant data values in the intervall ---------

    sumnonconstant = 0
	do k = m0,m2-1
	  if (abs(y(k) - y(k+1)) < 1.d-6) then
	    sumnonconstant = sumnonconstant + 1
	  end if
	end do

	if ((sumnonconstant < 3.d0*(m2-m0+1)/4.d0) .and. (count(w(m0:m2) == 0.d0) < 3.d0*(m2-m0+1)/4.d0)) then

!---- Initial estimates of fitting parameters -----------------------------------    
!     Loop over width and position in time     

      chi2min = 1.d60

	  do k = 1,5
	    do l = 1,5
          nmid = floor(m1 - m0 + (l-3.d0)*width/3.d0)
	      exp12(1:n) = ydl(nptperyear-nmid+1:nptperyear+n-nmid,k,nseason)
          ytilde(1:n) = wfit(m0:m2)*y(m0:m2)
          Btilde(1:n,1:2) = reshape(source = (/ wfit(m0:m2), wfit(m0:m2)*exp12(1:n) /), &
	                           shape = (/ n,2 /)) 
	      call gauss(2,matmul(transpose(Btilde(1:n,1:2)),Btilde(1:n,1:2)), &
	             matmul(ytilde(1:n),Btilde(1:n,1:2)),c)
    
!---- Evaluate the fitted function ----------------------------------------------    
   
          ylocal(i,m0:m2) = c(1) + c(2)*exp12(1:n)		  
		  chi2 = sum((wfit(m0:m2)*(ylocal(i,m0:m2)-y(m0:m2)))**2.d0)
		  if (chi2 < chi2min) then
		    chi2min = chi2
	        kmin = k
		    lmin = l
	      end if 

	    end do
	  end do 
		  
!---- Best starting values ------------------------------------------------------		  
		  
      x(1) = m1 - m0 - width*fact(kmin) + (lmin-3.d0)*width/3.d0
      x(3) = m1 - m0 + width*fact(kmin) + (lmin-3.d0)*width/3.d0
      x(2) = nptperyear/20.d0
      x(4) = nptperyear/20.d0

!---- Set constraints -----------------------------------------------------------

      bl = (/ x(1)-width/2.d0, 0.5d0*x(2), x(3)-width/2.d0, 0.5d0*x(4) /)
      bu = (/ x(1)+width/2.d0, 10.d0*x(2), x(3)+width/2.d0, 10.d0*x(4) /)
    
!---- Separable box constrained non-linear least squares fit -------------------- 

!     write(*,*) 'Log,x,chi2min',x,chi2min
	  call marquardt(2,term,x,xn,bl,bu,4,y(m0:m2),wfit(m0:m2),n) 
!     write(*,*) 'xn',xn
!	  pause

!---- Final linear fit to determine c(1) and c(2) -------------------------------

      t = dble((/(i,i=1,npt)/))
	  exp12(1:n) = 1.d0/(1.d0 + exp((xn(1)-t(1:n))/xn(2))) &
	                            - 1.d0/(1.d0 + exp((xn(3)-t(1:n))/xn(4)))
    
      ytilde(1:n) = wfit(m0:m2)*y(m0:m2)
      Btilde(1:n,1:2) = reshape(source = (/ wfit(m0:m2), wfit(m0:m2)*exp12(1:n) /), &
	                           shape = (/ n,2 /)) 
	  call gauss(2,matmul(transpose(Btilde(1:n,1:2)),Btilde(1:n,1:2)), &
	           matmul(ytilde(1:n),Btilde(1:n,1:2)),c)
    
!---- Evaluate the fitted function ----------------------------------------------    
   
      ylocal(i,m0:m2) = c(1) + c(2)*exp12(1:n)

    else 
	
!---- If many constant values in interval replace with median value -------------
	  
	  call median(y(m0:m2),m2-m0+1,ymed)
	  ylocal(i,m0:m2) = ymed
    
	end if

  end do ! end loop over local functions 
  
!---- Merge local functions to a global function -------------------------------- 
  
  do i = 2,nlocal
    n1 = s(i)
    n2 = s(i+1)
	ndiff = dble(n2 - n1)
    mid = (n1 + n2)/2.d0   
    do k = n1,n2   
      cut = (atan(10.d0*(k-mid)/ndiff) - atan(10.d0*(n1-mid)/ndiff))/ &
            (atan(10.d0*(n2-mid)/ndiff) - atan(10.d0*(n1-mid)/ndiff))
      yfit(k) = cut*ylocal(i,k) + (1.d0-cut)*ylocal(i-1,k)
    end do
  end do
  
!---- Add points in the beginning and at the end -------------------------------
  
  yfit(s(1):s(2)) = ylocal(1,s(1):s(2))
  yfit(s(nlocal+1):s(nlocal+2)) = ylocal(nlocal,s(nlocal+1):s(nlocal+2))  

!---- Modify weights -----------------------------------------------------------

  if (j < nenvi) then
    call modweight(j,npt,wfit,y,yfit)
  end if

end do ! end loop over envelope iterations

!---- Output for debug purposes ------------------------------------------------

! write(*,*) 3,sum(yfit) 

end subroutine fitlogistic