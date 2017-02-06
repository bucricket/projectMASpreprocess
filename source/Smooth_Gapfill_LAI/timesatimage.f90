!-------------------------------------------------------------------------------!
!                                                                               !

 program timesatimage

!  TIMESAT image version 2.3                                                    !
!                                                                               !
!  A package for extracting phenological pararameters from                      !
!  time-series of satellite sensor data.                                        !
!                                                                               !
!  Authors:                                                                     !
!  Per J\"onsson, Malm\"o University, Sweden                                    !
!  e-mail per.jonsson@ts.mah.se                                                 !
!                                                                               !
!  Lars Eklundh, Lund University, Sweden                                        !
!  e-mail lars.eklundh@nateko.lu.se                                             !
!                                                                               !
!  Updated September 2006                                                       !
!                                                                               !
!  Global variables                                                             !    
!  nyear           number of years                                              !  
!  nptperyear      number of datapoints per year                                !
!  npt             total number of datapoints npt = nyear*nptperyear            ! 
!  nseason         number of annual seasons                                     !      
!  nenvi           number of envelop iterations                                 !
!  plotflag        plot flag                                                    !
!  wfact           factor determining the strength of envelope adaptation       !
!  win(3)          vector with Savitzky-Golay window sizes                      !    
!                                                                               !
!  Local variables                                                              !                                              !
!  nc              number of columns                                            !
!  nrstart         start row                                                    !
!  nrstop          stop row                                                     !
!  ncstart         start column                                                 !
!  ncstop          stop column                                                  !
!  y(npt)          vector containing time-series for a pixel                    !
!  w(npt)          vector containing time-series for a pixel                    !
!  ymatrixchar(npt,nc) matrix containing time-series for all pixels in a row    !
!  wmatrixchar(npt,nc) matrix containing time-series for all pixels in a row    !
!  ymatrixint2(npt,nc) matrix containing time-series for all pixels in a row    !
!  wmatrixint2(npt,nc) matrix containing time-series for all pixels in a row    !
!  ymatrixreal(npt,nc) matrix containing time-series for all pixels in a row    !
!  wmatrixreal(npt,nc) matrix containing time-series for all pixels in a row    !
!  f1(npt,7)       sine basis used in subroutine SEASON                         !
!  f2(npt,7)       polynomial basis used in subroutine SAVGOL                   !
!  yag(2*nptperyear,7,2) basis used in fitgauss                                 !
!  ydl(2*nptperyear,5,2) basis used in fitlogistic                              !
!  s(:)            vector containing approximate maximum and minimum positions  !
!                                                                               !
!-------------------------------------------------------------------------------!

implicit none

integer          :: i,j,k,nr,nc,nrstart,nrstop,ncstart,ncstop,nyear,nptperyear
integer          :: npt,nenvi,printflag,debugflag,filetype,win(3),nlocal
integer          :: mask,nbyte,missing,sg,ag,dl,date_time(8),minimum,error
double precision :: yl,yu,wfact,bl(3),bu(3),weight(3),ymean,ystd,seasonpar
double precision :: startcutoff,amplitudecutoff,spikecutoff,timestart,timeend
character*20     :: jobname
character*80     :: namedata,namemask,imagefile
character*12     :: real_clock(3)

integer,          allocatable, dimension(:)     :: s
double precision, allocatable, dimension(:)     :: y,w,y1,y2,y3,diff
real*4,           allocatable, dimension(:)     :: ysp
double precision, allocatable, dimension(:,:)   :: f1,f2
double precision, allocatable, dimension(:,:,:) :: yag,ydl
real*4,           allocatable, dimension(:,:)   :: ymatrixreal,wmatrixreal
character,        allocatable, dimension(:,:)   :: ymatrixchar,wmatrixchar 
integer*2,        allocatable, dimension(:,:)   :: ymatrixint2,wmatrixint2  

common /parameters/ nyear,nptperyear,npt,nenvi,printflag,debugflag,wfact,win 

call date_and_time (real_clock(1), real_clock(2), real_clock(3), date_time) 
call cpu_time(timestart) 

write(*,*)
write(*,*) '  Timesat image version 2.3 for Fortran'
write(*,*) '  Per Jonsson and Lars Eklundh'
write(*,*) '  per.jonsson@ts.mah.se, lars.eklundh@nateko.lu.se' 
write(*,*) '  September 2006'
write(*,*)
write(*,*)

!---- File name -----------------------------------------------------------------

write(*,*) '   Name of the sensor data file list'  
read(*,*)    namedata
write(*,*)

!---- Set filetype --------------------------------------------------------------

write(*,*) '   Specify the file type for the image files '
write(*,*) '   1 = 8-bit unsigned integer (byte)'
write(*,*) '   2 = 16-bit signed integer'
write(*,*) '   3 = 32-bit real'
read(*,*) filetype              
write(*,*)

select case (filetype)
case (1)
  nbyte = 1
case (2)
  nbyte = 2
case (3)
  nbyte = 4
end select      

!---- Define image position for processing --------------------------------------

write(*,*) '   Number of rows and columns in the image tile ' 
read(*,*)    nr, nc
write(*,*)                                            
write(*,*) '   Define processing area: start row, stop row, start column, stop column'                       
read(*,*)    nrstart, nrstop, ncstart, ncstop                  
write(*,*)

!---- Number of years, number of points per year and total number of points -----

write(*,*) '   Number of years and number of points per year in the time-series'
read(*,*)    nyear, nptperyear                                 
write(*,*)

npt = nyear*nptperyear

!---- Define range for sensor data ----------------------------------------------

write(*,*) '   Sensor data values in the range [min,max] are accepted. Data outside'
write(*,*) '   this range are assigned weight 0. Give min and max '
read(*,*)    yl,yu
write(*,*)

!---- Define mask data ----------------------------------------------------------

write(*,*) '   Name of the mask data file list: write NONE if mask data is unavailable ' 
read(*,*)    namemask
if ( (trim(namemask) == 'NONE') .or. (trim(namemask) == 'none') ) then
  mask = 0
else
  mask = 1
end if

!---- Set boundaries and weights for mask data conversion -----------------------

if (mask == 1) then
  write(*,*) '   Sensor data values are assigned weights depending on the '
  write(*,*) '   corresponding mask data values: '
  write(*,*) '   if mask data in the range [a1,b1] then set weight to w1 '
  write(*,*) '   if mask data in the range [a2,b2] then set weight to w2 '
  write(*,*) '   if mask data in the range [a3,b3] then set weight to w3 '
  write(*,*) '   for mask data values outside these ranges weight is set to 0'  
  write(*,*) '   Give a1,b1,w1,a2,b2,w2,a3,b3,w3 '
  read(*,*)    bl(1),bu(1),weight(1),bl(2),bu(2),weight(2),bl(3),bu(3),weight(3)
  write(*,*)
end if

!---- Cutoff for low variation in sensordata ------------------------------------

write(*,*) '   Time-series with an amplitude less than A are not processed  '
write(*,*) '   A = 0 forces all time-series to be processed. Give A '
read(*,*)    amplitudecutoff
write(*,*)

!---- Cutoff for spike ----------------------------------------------------------

write(*,*) '   Single spikes are detected by a comparison with median filtered values '
write(*,*) '   and with closest neighbors. If the distance is greater than S*ystd,   ' 
write(*,*) '   where ystd is the standard deviation for data values, we have a spike. '
write(*,*) '   S = 2 is the normal value. Give S.'
read(*,*)    spikecutoff
write(*,*)

!---- Parameter for determining the number of seasons ---------------------------

write(*,*) '   Give parameter that is used to determine the number of annual seasons. '
write(*,*) '   Parameter value should be in the range [0,1]. A value close to 0 '
write(*,*) '   will force the program to interpret a small depression in the main '
write(*,*) '   curve as a second annual season. A value close to 1 will force the '
write(*,*) '   program to always use only one annual season. Give season parameter'
read(*,*)    seasonpar
write(*,*) 

!---- Number of envelope iterations ---------------------------------------------

write(*,*) '   Fitted curves adapt to the upper envelope of the sensor data '
write(*,*) '   values in an iterative procedure. Give the number of fitting steps: 1, 2 or 3'
read(*,*)    nenvi
write(*,*)  
write(*,*) '   Give the strength of the adaptation. Value should be in the '
write(*,*) '   range [1,10] where 2 is the normal value. '
read(*,*)    wfact
write(*,*)

wfact = 1.d0/wfact         

!---- Specify the processing methods --------------------------------------------

write(*,*) '   Specify the processing methods'
write(*,*) '   Savitzky-Golay (0/1), Asymmetric-Gauss (0/1), Double-Logistic (0/1)'
read(*,*)     sg,ag,dl
write(*,*) 

!---- Window size in Savgol for each of the fitting steps -----------------

if (sg == 1) then
  write(*,*) '   Savitzky-Golay window size for each of the fitting steps'
  read(*,*)    win(1:nenvi)
  write(*,*)
end if 

!---- Phenology parameter -------------------------------------------------------

write(*,*) '   The time for the season start (end) is defined as the time for which the'
write(*,*) '   sensor data value ,measured from the base level, has increased (decreased) '
write(*,*) '   to X % of the seasonal amplitude. X = 20 is the normal value. Give X '
read(*,*)    startcutoff
write(*,*)
 
startcutoff = startcutoff/100.d0

!-----Job name------------------------------------------------------------------

write(*,*) '   Give an identification tag for the job (text string of max 20 chars): '
read(*,*) jobname

!---- Allocate memory -----------------------------------------------------------

allocate ( y(npt),w(npt),y1(npt),y2(npt),y3(npt),ysp(npt),s(4*nyear+2),diff(npt-1) )
allocate ( f1(npt,5),f2(npt,3) )
allocate ( yag(2*nptperyear,7,2),ydl(2*nptperyear,5,2) )
allocate ( ymatrixchar(npt,nc),wmatrixchar(npt,nc) )
allocate ( ymatrixint2(npt,nc),wmatrixint2(npt,nc) )
allocate ( ymatrixreal(npt,nc),wmatrixreal(npt,nc) )

!---- Open sensordata file list -------------------------------------------------   

open(1,file=namedata,status='old',err=91)

!---- Get number of sensordata imagefiles in the list ---------------------------

read(1,*) npt
if (npt /= nyear*nptperyear) stop 'Wrong number of sensordata imagefiles'

!---- Read name of sensordata imagefiles and open corresponding files -----------

do i = 1,npt
  read(1,'(a)') imagefile
! open(10000+i,file=imagefile,status='old',access='stream',err=92)                ! alt1 g95
! open(10000+i,file=imagefile,status='old',access='direct',recl=nbyte*nc,err=92)  ! alt2 g95
  open(10000+i,file=imagefile,status='old',form='binary',recordtype='fixed', &    ! f95
       access='direct',recl=nbyte*nc,err=92) 
  write(*,*) '  Sensordata imagefile ',i, ' opened'    
end do 

close(1)

!---- Open cloudmask file list --------------------------------------------------

if (mask == 1) then
  open(1,file=namemask,status='old',err=93) 

!---- Get number of cloudmask imagefiles in the list ----------------------------

  read(1,*) npt
  if (npt /= nyear*nptperyear) stop 'Wrong number of cloudmask imagefiles'

!---- Read name of cloudmask imagefiles and open corresponding files ------------

  do i = 1,npt
    read(1,'(a)') imagefile
!   open(20000+i,file=imagefile,status='old',access='stream',err=94)                ! alt1 g95 
!   open(20000+i,file=imagefile,status='old',access='direct',recl=nbyte*nc,err=94)  ! alt2 g95 
    open(20000+i,file=imagefile,status='old',form='binary',recordtype='fixed', &    ! f95
         access='direct',recl=nbyte*nc,err=94) 

    write(*,*) '  Cloudmask imagefile ',i, ' opened'    
  end do 

  close(1)
end if

!---- Open output files ---------------------------------------------------------

open(13,file='phenologySG_'//trim(adjustl(jobname)),form='binary',status='unknown') 
open(14,file='phenologyAG_'//trim(adjustl(jobname)),form='binary',status='unknown')
open(15,file='phenologyDL_'//trim(adjustl(jobname)),form='binary',status='unknown') 
open(16,file='sensordata_'//trim(adjustl(jobname)),form='binary',status='unknown') 
open(17,file='fitSG_'//trim(adjustl(jobname)),form='binary',status='unknown') 
open(18,file='fitAG_'//trim(adjustl(jobname)),form='binary',status='unknown')
open(19,file='fitDL_'//trim(adjustl(jobname)),form='binary',status='unknown') 
open(20,file='input_'//trim(adjustl(jobname))//'.txt',status='unknown')

!---- Write input parameters to file --------------------------------------------

write(20,'(a)') namedata
write(20,'(10i6)') filetype
write(20,'(10i6)') nr,nc
write(20,'(10i6)') nrstart,nrstop,ncstart,ncstop
write(20,'(10i6)') nyear,nptperyear
write(20,'(10f11.4)') yl,yu
write(20,'(a)') namemask
if (mask == 1) write(20,'(10f11.4)') bl(1),bu(1),weight(1),bl(2),bu(2),weight(2), &
                                                bl(3),bu(3),weight(3)
write(20,'(10f11.4)') amplitudecutoff
write(20,'(10f11.4)') spikecutoff
write(20,'(10f11.4)') seasonpar
write(20,'(10i6)') nenvi
write(20,'(10f11.4)') 1.d0/wfact
write(20,'(10i6)') sg,ag,dl
if (sg == 1) write(20,'(10i6)') win(1:nenvi)
write(20,'(10f11.4)') startcutoff*100.d0
write(20,'(a)') trim(adjustl(jobname))

!---- Print header information to outputfiles -----------------------------------

do i = 13,19
  write(i) nyear, nptperyear, nrstart, nrstop, ncstart, ncstop 
end do

!---- Generate basis functions needed for SEASON, SAVGOL, FITGAUSS, FITLOGISTIC -

call basis(npt,nptperyear,f1,f2,yag,ydl)

!---- Read data from sensordata and cloudmask imagefiles ------------------------
!     Loop over rows, for each row loop over all npt imagefiles                 

write(*,*)
write(*,*) '  Processing started'
write(*,*)

do i = nrstart,nrstop        ! f95
!do i = 1,nrstop             ! g95

!---- Read time-series for all pixels (columns) in the row ----------------------   

  select case (filetype)
  case (1)
    do j = 1,npt
      read(10000+j,rec=i) ymatrixchar(j,:)
      if (mask == 1) read(20000+j,rec=i) wmatrixchar(j,:)
    end do
  case (2)
    do j = 1,npt
      read(10000+j,rec=i) ymatrixint2(j,:)
      if (mask == 1) read(20000+j,rec=i) wmatrixint2(j,:)
    end do
  case (3)
    do j = 1,npt
      read(10000+j,rec=i) ymatrixreal(j,:)
      if (mask == 1) read(20000+j,rec=i) wmatrixreal(j,:)
    end do
  end select
  
!  if (i < nrstart) cycle    ! g95
  
  write(*,*) '  Row ',i

  do j = ncstart,ncstop

!    write(*,*) '  Col ',j
  
!---- Transfer time-series of sensor and cloudmask data to vektors y and w ------
 
    select case (filetype)
    case (1)
      y = dble( ichar(ymatrixchar(:,j)) )
      if (mask == 1) then
	    w = dble( ichar(wmatrixchar(:,j)) )
      else
	    w = 1.d0
	  end if
    case (2)
      y = dble( ymatrixint2(:,j) )
      if (mask == 1) then
	    w = dble( wmatrixint2(:,j) )
      else
	    w = 1.d0
	  end if
    case (3)
      y = dble( ymatrixreal(:,j) )
      if (mask == 1) then
	    w = dble( wmatrixreal(:,j) )
      else
	    w = 1.d0
	  end if
    end select
      
!---- Convert mask data to weights used in the fitting procedure --------------

    if (mask == 1) then
	  where ( (bl(1) <= w) .and. (w <= bu(1)) )  
	    w = weight(1)
	  elsewhere ( (bl(2) <= w) .and. (w <= bu(2)) ) 
	    w = weight(2)
	  elsewhere ( (bl(3) <= w) .and. (w <= bu(3)) ) 
	    w = weight(3)
	  elsewhere
	    w = 0.d0
	  end where
	end if   

!---- Check if sensordata is in the specified range, if not assign weight zero --

    where ( (y < yl) .or. (y > yu) ) w = 0.d0
      
!---- Time-series with too many data values with weight zero are not processed --

    missing = 0

    if (count(w == 0.d0) >= 3.d0*npt/4.d0) missing = 1
    do k = 1,npt-nptperyear/3
	  if ( abs(sum(w(k:k+nptperyear/3))) == 0.d0 ) missing = 1
	end do

	if (missing == 1) cycle  

!---- Only time-series with non-constant datavalues are processed ---------------

    do k = 1,npt-1
      diff(k) = abs(y(k) - y(k+1)) 
    end do
!    if (count(diff < 1.d-6) >= npt/4.d0 ) cycle      
!    if (count(diff < 1.d-6) >= 3*npt/4 ) cycle   ! PIETER      
    if (count(diff < 1.d-6) >= 3.d0*npt/4.d0 ) cycle   ! Duarte  

!---- Identify spikes in the time-series, set corresponding weigths to zero -----   
    
	error = 0  
    call spike(y,w,spikecutoff,error)
	if (error == 1) cycle

!---- Fit sine and cosine functions to determine the number of annual seasons ---
!     The current version of the code only allows one or two annual seasons     
!     The function returns a vector s(nlocal) that defines the intervals for subsequent 
!     fits to local functions. If the amplitude of the fitted curve is less than the
!     amplitude cutoff then nlocal is set to zero and we cycle                                                    

    call season(y,w,f1,s,nlocal,seasonpar,amplitudecutoff,minimum) 
	if (nlocal == 0) cycle
      
!    if (mod(j,100) == 0) write(*,*) 'Row, Col',i,j
    ysp = y
	write(16) i,j
	write(16) ysp

!---- Iterative Savitzky-Golay filtering, analyze and dump phenological data ----  

    if (sg == 1) then
	  call savgol(y,w,f2,y1,maxval(win))
	  call phenology(13,i,j,s,y,y1,nlocal,startcutoff,minimum)
      ysp = y1
      write(17) i,j
      write(17) ysp

	end if
     
!---- Iterative fits to asymmetric Gaussians, analyze and dump phenological data -      

    if (ag == 1) then      
      call fitgauss(s,yag,y,w,nlocal,y2)
	  call phenology(14,i,j,s,y,y2,nlocal,startcutoff,minimum)
	  ysp = y2
      write(18) i,j
	  write(18) ysp

	end if

!---- Iterative fits to double logistic functions, analyze and dump phenological data            
    
	if (dl == 1) then  
      call fitlogistic(s,ydl,y,w,nlocal,y3)       
	  call phenology(15,i,j,s,y,y3,nlocal,startcutoff,minimum)
      ysp = y3
      write(19) i,j
	  write(19) ysp

	end if	
	                                        
  end do ! end loop over columns 
end do ! end loop over rows

!---- Close all imagefiles ------------------------------------------------------

do i = 1,npt
  close(10000+i)
  close(20000+i)
end do
	
!---- Close remaining files -----------------------------------------------------

do i = 13,20
  close(i)
end do

!----- Print information to screen ----------------------------------------------

write(*,*)
write(*,*) '  Processing finished  ' 
write(*,*) 
write(*,*) '  Phenological parameters written to:'
if (sg == 1) then 
  write(*,*) '  phenologySG_'//trim(adjustl(jobname))
end if
if (ag == 1) then 
  write(*,*) '  phenologyAG_'//trim(adjustl(jobname))
end if
if (dl == 1) then 
  write(*,*) '  phenologyDL_'//trim(adjustl(jobname))
end if
write(*,*) 
write(*,*) '  Sensor data and fitted functions written to:'
write(*,*) '  sensordata_'//trim(adjustl(jobname))
if (sg == 1) then
  write(*,*)  '  fitSG_'//trim(adjustl(jobname))
end if
if (ag == 1) then
  write(*,*)  '  fitAG_'//trim(adjustl(jobname))
end if
if (dl == 1) then
  write(*,*)  '  fitDG_'//trim(adjustl(jobname))
end if
write(*,*) 
write(*,*) '  Input data written to:'
write(*,*) '  input_'//trim(adjustl(jobname))//'.txt'

call cpu_time(timeend)
write(*,*)
write(*,*) '  CPU time used ',timeend - timestart,' seconds' 

stop

91 stop ' Error opening sensordatalist'
92 stop ' Error opening sensordata imagefiles'  
93 stop ' Error opening maskdatalist'
94 stop ' Error opening cloudmask imagefiles'     
       
end program timesatimage
