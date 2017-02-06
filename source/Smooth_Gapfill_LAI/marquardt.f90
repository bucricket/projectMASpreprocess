!-------------------------------------------------------------------------------!
!                                                                               !

subroutine marquardt(model,term,x,xn,bl,bu,n,y,w,m)

!  Marquardt's method for least squares                                         !
!  Find x that minimizes F(x) = 0.5*sum_{i=1}^m (f_i(x1,x2,..,xn)^2)            !
!  under box constraints bl(i) <= xi <= bu(i), i = 1,2,...,n                    !   
!                                                                               !
!  Refs:                                                                        !  
!  Hans Bruun Nielsen, IMM, DTU, Denmark                                        !   
!                                                                               !
!  Authors:                                                                     !
!  Per J\"onsson, Malm\"o University, Sweden                                    !
!  e-mail per.jonsson@ts.mah.se                                                 !
!                                                                               !
!  Lars Eklundh, Lund University, Sweden                                        !
!  e-mail lars.eklundh@nateko.lu.se                                             !
!                                                                               !
!  Updated 17/8-2005                                                            !
!                                                                               !
!  Variables                                                                    !
!  x(n)      row vector with parameter values                                   !    
!  f(m)      row vector with function values                                    !
!  J(m,n)    Jacobian matrix                                                    !
!                                                                               !
!  Termination criteria                                                         !
!  term = 1  terminated by small gradient                                       !
!  term = 2  terminated by small x-step                                         !
!  term = 3  singular matrix, restart from current x with increased             !
!            value for mu                                                       !
!  term = 4  terminated by kmax                                                 !
!                                                                               !
!-------------------------------------------------------------------------------!

implicit none

integer          :: m,n,i,jjj,k,kmax,term,model
double precision :: delta,x(n),xn(n),bl(n),bu(n),y(m),w(m),f(m),fn(m),J(m,n)
double precision :: Jn(m,n),eye(n,n),h(n),dL,g(n),F_,Fn_,dF,ng,mu,nu,nh,nx,A(n,n)

A = 0.d0
eye = 0.d0
do i = 1,n
  eye(i,i) = 1.d0
end do
kmax = 100                    ! 100
delta = 1.d-4                 ! 1.d-8           
if (model == 1) then 
  call fungauss(m,f,J,x,y,w)
else if (model == 2) then
  call funlogistic(m,f,J,x,y,w)
end if              
A(1:n,1:n) = matmul(transpose(J),J)                        
g = matmul(f,J)                          
F_ = dot_product(f,f)/2.d0                  
ng = maxval(abs(g))                      
mu = 1.d0*maxval((/(A(i,i),i=1,n)/))         
k = 1   
nu = 2.d0   
nh = 0.d0   
term = 0

do while (term == 0) 
  do i = 1,n
    do jjj = 1,n 
      !if (isnan(A(i,jjj)) == .true.) then
       if (isnan(A(i,jjj))) then
	    write(*,*) 'NaN in marquardt'
		term = 99
		return
	  end if
	end do
  end do                     
  if (ng <= delta) then
    term = 1  
  else 
    call gauss(n,A(1:n,1:n)+mu*eye,-g,h)          
    nh = sqrt(dot_product(h,h))                                                
    nx = delta + sqrt(dot_product(x,x))  
    if (nh <= delta*nx .and. k > 1) then  
      term = 2
    else if (nh >= nx/epsilon(1.d0)) then ! Almost singular, restart from current x with increased value for  mu .  
      term = 3 
	  write(*,*) 'stop in marquard'
	  stop
    end if   
  end if
  if (term == 0) then                     
    xn = x + h               
          
    where (xn < bl) xn = bl               ! Project onto the feasible parameter set
	where (xn > bu) xn = bu
    
    dL =  dot_product(h,mu*h-g)/2.d0
	if (model == 1) then       
      call fungauss(m,fn,Jn,xn,y,w)
	else if (model == 2) then
	  call funlogistic(m,fn,Jn,xn,y,w)
	end if       
    Fn_ = dot_product(fn,fn)/2.d0            
    dF = F_ - Fn_
    if  ((dL > 0) .and. (dF > 0)) then    ! Update x and modify mu
      x = xn   
      F_ = Fn_  
      J = Jn  
      f = fn 
      A(1:n,1:n) = matmul(transpose(J),J)              
	  g = matmul(f,J)                         
      ng = maxval(abs(g))                  
      mu = mu*max(1.d0/3.d0, 1 - (2.d0*dF/dL - 1.d0)**3.d0)   
	  nu = 2.d0
    else                                  ! Same  x, increase  mu
	  mu = mu*nu  
      nu = 2.d0*nu 
    end if
    k = k + 1
    if  (k > kmax) then
      term = 4 
    end if      
  end if  
end do

end subroutine marquardt
