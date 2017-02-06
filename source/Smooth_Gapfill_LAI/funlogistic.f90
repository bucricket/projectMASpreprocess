!-------------------------------------------------------------------------------!
!                                                                               !

subroutine funlogistic(n,f,J,x,y,w)
  
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
!  Local variables                                                              !
!  x(n)      vector with parameter values                                       !    
!  f(m)      vector with function values                                        !
!  J(m,n)    Jacobian matrix                                                    !
!                                                                               !
!                                                                               !
!-------------------------------------------------------------------------------!

implicit none

integer :: i,n
double precision :: zeros(n),t(n),exp1(n),exp2(n),exp3(n),exp4(n),x(4),y(n)
double precision :: ytilde(n),w(n),f(n),J(n,4),J1(n),J2(n),J3(n),J4(n),Jb(n,4)
double precision :: A(2,2),Btilde(n,2),Btildeprim(n,2),c(2),cprim(2)

zeros = 0.d0
t = dble((/ (i, i = 1,n) /)) 

exp1 = exp((x(1)-t)/x(2))
exp2 = exp((x(3)-t)/x(4))
exp3 = exp1/(1+exp1)**2.d0
exp4 = exp2/(1+exp2)**2.d0

!---- Compute the Jacobian for the logistic function ------------------------------

J1 = -exp3/x(2)
J2 = (x(1)-t)*exp3/x(2)**2.d0
J3 = exp4/x(4)
J4 = -(x(3)-t)*exp4/x(4)**2.d0
Jb = reshape(source = (/ J1,J2,J3,J4 /), shape = (/ n,4 /))

ytilde = w*y
Btilde = reshape(source = (/ w, w*(1/(1+exp1) - 1/(1+exp2)) /), shape = (/ n,2 /))    
A(1:2,1:2) = matmul(transpose(Btilde),Btilde)                        
call gauss(2,A,matmul(ytilde,Btilde),c)                     
f = ytilde - matmul(c,transpose(Btilde))               

do i = 1,4
  Btildeprim = reshape(source = (/ zeros, w*Jb(1:n,i)  /), shape = (/ n,2 /))
  call gauss(2,matmul(transpose(Btilde),Btilde),matmul(f,Btildeprim) &
                       - matmul(c,matmul(transpose(Btildeprim),Btilde)),cprim) 
  J(1:n,i) = -matmul(c,transpose(Btildeprim)) - matmul(cprim,transpose(Btilde)) 
end do

end subroutine funlogistic