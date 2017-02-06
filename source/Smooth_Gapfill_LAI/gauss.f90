!-------------------------------------------------------------------------------!
!                                                                               !

 subroutine gauss(n,a,b,c)                                        

!                                                                               !
! This subroutine solves the n x n linear system ac = b                         !
! by Gauss-elimination with pivoting                                            !
!                                                                               !
! On input:                                                                     !
! n              dimension of the system                                        !
! a(n,n)         coefficient matrix                                             !
! b(n)           righthand side                                                 !
!                                                                               !
! On output:                                                                    !
! c(n)           solutions to the system                                        !
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
 
 integer :: i,j,k,n,pos(1)
 double precision :: a(n,n),b(n),c(n),f,temp1(n),temp2
      
!---- Elimination ---------------------------------------------------------------
           
 do i = 1,n-1
   pos(1:1) = maxloc(abs(a(i:n,i)))
   if (pos(1) > 1) then
     temp1  = a(i,:)
     temp2  = b(i)
     a(i,:) = a(i+pos(1)-1,:)
     b(i)   = b(i+pos(1)-1)
     a(i+pos(1)-1,:) = temp1
     b(i+pos(1)-1)   = temp2
   end if
   do j = i+1,n
     f = - a(j,i)/a(i,i)
     a(j,:) = a(j,:) + f*a(i,:)
     b(j) = b(j) + f*b(i)
   end do
 end do

!---- Back substitution ---------------------------------------------------------
 
 do i = 1,n
   if (abs(a(i,i)) < 1d-80) then
     c = 0.d0
!	 write(*,*) 'Bad condition in gauss'
     return
   end if
 end do

 c(n) = b(n)/a(n,n)
 do i = n-1,1,-1
   c(i) = (b(i) - sum(a(i,i+1:n)*c(i+1:n)) )/a(i,i) 	
 end do
    
 end subroutine gauss
