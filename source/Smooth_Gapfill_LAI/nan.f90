program nantest
!====NAN.F90 illustrates what works and what doesn't when
!    detecting a NaN
! Platforms: Windows 9x/Me/NT/2000, AIX 4.3, Linux
! Compiler Notes:
!  Compaq Visual Fortran 6.6a with default
!     fpe settings (/fpe:3 /nocheck) and /OPTIMIZE:0
!     (ISNAN is an Elemental Intrinsic Function)
!     Options /fpe:0 /traceback will cause this
!     program to stop with error message,
!     James Giles points out that in order to actually print
!     minus zero from CVF, you have to compile with the
!     /assume:minus0  option.
!
!   AIX XLF90 without optimization.
!   (ISNAN is part of a BOS Runtime C library;
!   thus ISNAN must be declared LOGICAL.)
!
!   Linux Intel IFORT, use O0 (no optimization).
!
! Author: hdkLESS at SPAM psu dot edu
! Date: March, 2002, January 2005.
!
! References:
! http://www.psc.edu/general/software/packages/ieee/ieee.html
! http://homepages.borland.com/efg2lab/Mathematics/NaN.htm
!
       logical :: ISNAN
       integer :: i
       real, Dimension (6) :: y
       real ::  PInf, MInf, MZero, DivNan
       Character (Len=10), Dimension(6) :: NType
       data NType/'+Infinity','-Infinity','-0', 'NaN','NaN','0/0'/
       data PInf/B'01111111100000000000000000000000'/    ! +Infinity
       data MInf/B'11111111100000000000000000000000'/    ! -Infinity
       data MZero/B'10000000000000000000000000000000'/   ! -0

       data y(1)/B'01111111100000000000000000000000'/       ! +Infinity
       data y(2)/B'11111111100000000000000000000000'/       ! -Infinity
       data y(3)/B'10000000000000000000000000000000'/       ! -0
       data y(4)/B'01111111100000100000000000000000'/       ! NaN
       data y(5)/B'11111111100100010001001010101010'/       ! NaN
       DivNan=0
       y(6)=DivNan/DivNan

       Do i=1,6
        print *, 'Test#',i,' ->',NType(i)
        if (y(i).eq.PInf) print *, 'Y = Plus Infinity'
        if (y(i).eq.MInf) print *, 'Y = Minus Infinity'
        if (y(i).eq.Mzero) print *, 'Y = Minus Zero'
        print *, 'y(Test#)=',y(i)
        print *, 'Testing each of three NaN detection methods:'
! EQV -> true iff both A and B are True or iff both A and B are False.
        if( (y(i) > 0.0) .EQV. (y(i) <= 0.0)) then
           print *, '1) (y(Test#) > 0.0) .EQV. (y(Test#) <= 0.0)'
        end if
        if (y(i)/=y(i)) then
           print *, '2) (y(Test#)/=(y(Test#))'
        end if
!        write(*,*) (IBIT(TRANSFER(y(i),1),23,8) == 2047)
! If ISNAN is not available for a specific compiler, comment the
! following line.
!        if (ISNAN(y(i))) then
!          print *, '3) ISNAN(y(Test#))'
!        end if
        print *, ' '
       End Do

 end program nantest
