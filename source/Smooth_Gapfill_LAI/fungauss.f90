!-------------------------------------------------------------------------------!
!                                                                               !

subroutine fungauss(n,f,J,x,y,w)
  
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

integer :: i,n,nmid
double precision :: zeros(n),t(n),arg(n),exparg(n),x(5),y(n),ytilde(n),w(n),f(n)
double precision :: J(n,5),J1(n),J2(n),J3(n),J4(n),J5(n),Jb(n,5),A(2,2)
double precision :: Btilde(n,2),Btildeprim(n,2),c(2),cprim(2)

nmid = floor(x(1)) 
zeros = 0.d0
t = dble((/ (i, i = 1,n) /)) 
arg = (/ ((x(1)-t(1:nmid))/x(2))**x(3),((t(nmid+1:n)-x(1))/x(4))**x(5) /)
exparg = exp(-arg)

!---- Compute the Jacobian for asymmetric Gaussian ------------------------------

J1 = (/ -(x(3)/(x(1)-t(1:nmid)))*arg(1:nmid)*exparg(1:nmid), & 
         (x(5)/(t(nmid+1:n)-x(1)))*arg(nmid+1:n)*exparg(nmid+1:n) /)
J2 = (/ (x(3)/x(2))*arg(1:nmid)*exparg(1:nmid), zeros(nmid+1:n) /)
J3 = (/ -log((x(1)-t(1:nmid))/x(2))*arg(1:nmid)*exparg(1:nmid), zeros(nmid+1:n) /)
J4 = (/ zeros(1:nmid) , (x(5)/x(4))*arg(nmid+1:n)*exparg(nmid+1:n) /)
J5 = (/ zeros(1:nmid) , -log((t(nmid+1:n)-x(1))/x(4))*arg(nmid+1:n)*exparg(nmid+1:n) /)
Jb = reshape(source = (/ J1,J2,J3,J4,J5 /), shape = (/ n,5 /))

ytilde = w*y
Btilde = reshape(source = (/ w, w*exparg /), shape = (/ n,2 /))    
A(1:2,1:2) = matmul(transpose(Btilde),Btilde)                        
call gauss(2,A,matmul(ytilde,Btilde),c)                     
f = ytilde - matmul(c,transpose(Btilde))               

do i = 1,5
  Btildeprim = reshape(source = (/ zeros, w*Jb(1:n,i)  /), shape = (/ n,2 /))
  call gauss(2,matmul(transpose(Btilde),Btilde),matmul(f,Btildeprim) &
                 - matmul(c,matmul(transpose(Btildeprim),Btilde)),cprim) 
  J(1:n,i) = -matmul(c,transpose(Btildeprim)) - matmul(cprim,transpose(Btilde)) 
end do

end subroutine fungauss