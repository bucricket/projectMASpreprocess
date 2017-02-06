program func2img
implicit none

character        :: paralist*80,namedata*80,namemask*80,date*8,version*80,jobname*80
character        :: namefunc*80,nameout*80,nameoutf*80,num*4,stopcode*80
integer          :: i,j,k,nr,nc,nrstart,nrstop,ncstart,ncstop,nyear,nptperyear
integer          :: npt,nenvi,printflag,debugflag,filetype,win(3),row,col
integer          :: mask,nbyte,missing,sg,ag,dl,seasonpara,header
double precision :: yl,yu,wfact,bl(3),bu(3),weight(3),ymean,ystd,seasonpar
double precision :: startcutoff,amplitudecutoff,spikecutoff
integer          :: datemin,datemax,i1,i2,imno
real             :: misspix

integer*2,        allocatable, dimension(:,:) :: imageint2
real*4,           allocatable, dimension(:,:) :: imagereal
real*4,           allocatable, dimension(:)   :: imagedata
character,        allocatable, dimension(:,:) :: imagechar 
101 format(a,a)
write(*,*)
write(*,*)' Program func2img'
write(*,*)' Program for generating images from TIMESAT fitted function files'
write(*,*)
write(6,*)
write(6,*)' Name of input fitted function file'
read(5,*) namefunc
open(9,file=namefunc,status='old',form='binary',err=998)
write(*,*)' Code for missing function value for pixel '
read(5,*) misspix
write(*,*)' Name of output files (no extension!)'
read(5,*) nameout
write(*,*) ' Specify the file type for the output image files (1,2 or 3)'
write(*,*) ' 1 = 8-bit binary'
write(*,*) ' 2 = 16-bit signed integer'
write(*,*) ' 3 = 32-bit real'
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

read(9,end=779)nyear,nptperyear,nrstart, nrstop, ncstart, ncstop
write(6,103) '  Number of years: ',nyear, ',  Number of points per year: ',nptperyear
write(6,103) '  Image window (rows, cols): ',nrstart, ' - ',nrstop, ' ; ',ncstart, ' - ',ncstop
write(6,*)
103 format(a,i4,a,i4,a,i4,a,i4)

write(6,*)'  Image number to map (-1 = all)'
read(5,*)imno

!total no. of images
npt = nyear * nptperyear

!--- Open output files (Compaq / Digital Fortran 90 compiler)
if (imno == -1) then
  !Open all output files
  do i = 1,npt
    !write the number to internal character variable num
    write(num,110) i
110 format(i4.4) 
    nameoutf=trim(nameout)//'_'//trim(adjustl(num))
    write(6,101)trim(nameoutf),' opened'
    open(100+i,file=nameoutf,status='unknown',form='binary',err=999)
  end do
else
  !write the number to internal character variable num
  write(num,110) imno
  nameout=trim(nameout)//'_'//trim(adjustl(num))
  open(100+imno,file=nameout,status='unknown',form='binary',err=999)
  write(6,101)trim(nameout),' opened'
end if   
   
ftype: select case (filetype)

  !!!!!!!!!!!!!!!!!!!!!filetype = 8-bit integer !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
case (1) 
  ! -------- Allocate memory -----------------------------------------------------------
  allocate ( imagechar(ncstart:ncstop,npt) ) !image matrix for fitted data
  allocate ( imagedata(npt) ) ! variable for holding the function data
 
  read(9,end=779)row,col     ! Read first row, column 
  read(9)imagedata           ! Read data for each point in time
  do i = nrstart,nrstop      ! For each row in output images
    imagechar = char(nint(misspix))! Put missing data code in all pixels
    do while (row == i)    
      imagechar(col,:) = char(nint(imagedata))
      read(9,end=771)row,col
      read(9)imagedata
    end do
    !write one row of data to files, all columns in row
    
771 if (imno == -1) then
      do k = 1,npt
        write(100+k)imagechar(:,k)
      end do
    else
      write(100+imno)imagechar(:,imno) 
    end if
  end do !rows
! close all files
  close(9)
  if (imno == -1) then
    do k = 1,npt
      close(100+k)
    end do   
  else
    close(100+imno)
  end if   
  write(6,101)'  Seasonality data written'
  write(6,101)'  File format 8-bit integer'

  !!!!!!!!!!!!!!!!!!!!!filetype = 16-bit integer !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
case (2) 
  ! -------- Allocate memory -----------------------------------------------------------
  allocate ( imageint2(ncstart:ncstop,npt) ) !image matrix for fitted data
  allocate ( imagedata(npt) ) ! variable for holding the function data
 
  read(9,end=779)row,col     ! Read first row, column 
  read(9)imagedata           ! Read data for each point in time
  do i = nrstart,nrstop      ! For each row in output images
    imageint2 = nint(misspix)! Put missing data code in all pixels
    do while (row == i)    
      imageint2(col,:) = nint(imagedata)
      read(9,end=772)row,col
      read(9)imagedata
    end do
    !write one row of data to files, all columns in row
    
772 if (imno == -1) then
      do k = 1,npt
        write(100+k)imageint2(:,k)
      end do
    else
      write(100+imno)imageint2(:,imno) 
    end if
  end do !rows
! close all files
  close(9)
  if (imno == -1) then
    do k = 1,npt
      close(100+k)
    end do   
  else
    close(100+imno)
  end if   
  write(6,101)'  Seasonality data written'
  write(6,101)'  File format 16-bit integer'


 !!!!!!!!!!!!!!!!!!!!!filetype = 32-bit real !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
case (3) !filetype = 32-bit real
  ! -------- Allocate memory -----------------------------------------------------------
  allocate ( imagereal(ncstart:ncstop,npt) ) !image matrix for fitted data
  allocate ( imagedata(npt) ) ! variable for holding the function data
 
  read(9,end=779)row,col     ! Read first row, column 
  read(9)imagedata           ! Read data for each point in time
  do i = nrstart,nrstop      ! For each row in output images
    imagereal = misspix      ! Put missing data code in all pixels
    do while (row == i)    
      imagereal(col,:) = imagedata
      read(9,end=773)row,col
      read(9)imagedata
    end do
    !write one row of data to files, all columns in row
    
773 if (imno == -1) then
      do k = 1,npt
        write(100+k)imagereal(:,k)
      end do
    else
      write(100+imno)imagereal(:,imno) 
    end if
  end do !rows
! close all files
  close(9)
  if (imno == -1) then
    do k = 1,npt
      close(100+k)
    end do   
  else
    close(100+imno)
  end if   
  write(6,101)'  Seasonality data written'
  write(6,101)'  File format 32-bit real'

end select ftype

stop

779 write(6,*) '  File is empty'
stop
998 write(6,*) '  Error opening input file '//trim(namefunc)
stop
999 write(6,*) '  Error opening output file '//nameout
stop 

end
