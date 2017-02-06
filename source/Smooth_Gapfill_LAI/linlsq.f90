!-------------------------------------------------------------------------------!
!                                                                               !

subroutine linlsq(n,m,f,y,w,c)

! This subroutine performs a fit of the datavalues (i,y(i)), i = 1,..,n         !
! to the linear modelfunction c(1)*f(:,1) + c(2)*f(:,2) + ... + c(m)*f(:,m)     !    
!                                                                               !
! see A. Garcia, Numerical methods for physics, ISBN 0-13-906744-2              !
! pages 142-146                                                                 !
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

integer :: i,j,n,m
double precision :: f(n,m),y(n),w(n),wsq(n),a(m,m),b(m),c(m)

!---- Construct the matrices a = A^T A  and b = A^T y --------------------------

wsq = w**2.d0         
do i = 1,m
  do j = i,m
    a(i,j) = sum(f(:,i)*f(:,j)*wsq)
  end do 
  b(i) = sum(f(:,i)*y*wsq)
end do
  
do i = 2,m
  do j = 1,i-1
    a(i,j) = a(j,i)   
  end do
end do
  
!---- Solve the m x m normal equation  A^T A c = A^T b -------------------------
      
call gauss(m,a,b,c)

end subroutine linlsq
