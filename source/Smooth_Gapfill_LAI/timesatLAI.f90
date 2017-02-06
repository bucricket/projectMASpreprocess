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
  !  Updated January 2006                                                         !
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
  double precision :: startcutoff,amplitudecutoff,spikecutoff
  character*20     :: jobname
  character*180    :: namedata,namemask,imagefile
  character*12     :: real_clock(3)

  integer,          allocatable, dimension(:)   :: s
  double precision, allocatable, dimension(:)   :: y,w,y1,y2,y3,diff
  real*4,           allocatable, dimension(:)   :: ysp
  double precision, allocatable, dimension(:,:) :: f1,f2
  double precision, allocatable, dimension(:,:,:) :: yag,ydl
  real*4,           allocatable, dimension(:,:) :: ymatrixreal,wmatrixreal
  !character,        allocatable, dimension(:,:) :: ymatrixchar,wmatrixchar 
  !integer*2,        allocatable, dimension(:,:) :: ymatrixint2,wmatrixint2

  common /parameters/ nyear,nptperyear,npt,nenvi,printflag,debugflag,wfact,win 

  ! parameters for HDF file - Feng
  integer            status, gddetach, gdclose
  integer*4          gdfid, gdid, gdopen, gdattach
  integer*4          rank, ntype, size 
  integer*4          gdgridinfo
  integer            sdid, sfstart
  integer            index, sfn2index
  integer            sfselect
  integer            sffattr, sfrnatt, sfrdata, att_id 

  integer            lai_id(512), fpar_id(512), laiQC_id(512), laiQC_ext_id(512)
  integer            ncols, nrows 
  integer            start(2), length(2), stride(2)
  real*8             upleft(2), lowright(2)
  real*8             scale_factor
  double precision,  allocatable, dimension(:)   :: w1
  integer*1,         allocatable, dimension(:) :: buf, csp, qc, extqc
  integer*1,         allocatable, dimension(:,:) :: fpar_matrixchar, lai_matrixchar
  integer*1          fillvalue, range(2)
  integer            bits(0:7), tmpbits(0:2), SCF_QC, CLD_QC, DET_QC, ii 
  integer            num_valid, ind(512)
  real,              allocatable, dimension(:)   :: ty, sy
  real                sumx, sumx2, stdev, dy, ydiff
  integer             nt, flag 
  real                HIGH_LAT_ULY
  ! 50 degrees and higher 
  parameter (HIGH_LAT_ULY = 926.625433*1200*6-100)
  integer             isHighLatitudeTile
  integer             MIN_SAT_POINT
  parameter (MIN_SAT_POINT = 60)

  integer DFACC_READ
  parameter (DFACC_READ=1)
  ! end of HDF parameters

  ! debug a single pixel for a given location - Feng
  integer debug_icol, debug_irow
  parameter (debug_icol = 565)
  parameter (debug_irow = 293)


  call date_and_time (real_clock(1), real_clock(2), real_clock(3), date_time) 

  write(*,*)
  write(*,*) '  Timesat image version 2.2 for Fortran'
  write(*,*) '  Per Jonsson and Lars Eklundh'
  write(*,*) '  per.jonsson@ts.mah.se, lars.eklundh@nateko.lu.se' 
  write(*,*) '  January 2006'
  write(*,*) '  Adapted for MODIS LAI process'
  write(*,*) '  March 2006' 
  write(*,*)
  write(*,*)

  !---- File name -----------------------------------------------------------------

  write(*,*) '   Name of the sensor data file list'  
  read(*,*)    namedata
  write(*,*)

  !---- Open sensordata file list -------------------------------------

  open(1,file=namedata,status='old',err=91)

  !---- Get number of sensordata imagefiles in the list ---------------------------

  read(1,*) npt
  !if (npt /= nyear*nptperyear) stop 'Wrong number of sensordata imagefiles'

  !---- Read name of sensordata imagefiles and open corresponding files -----------

  do i = 1,npt
     read(1,'(a)') imagefile
     ! open(10000+i,file=imagefile,status='old',access='stream',err=92)                ! alt1 g95
     ! open(10000+i,file=imagefile,status='old',access='direct',recl=nbyte*nc,err=92)  ! alt2 g95
     !  open(10000+i,file=imagefile,status='old',form='binary',recordtype='fixed', &    ! f95
     !       access='direct',recl=nbyte*nc,err=92) 

     ! open image file as a HDF file and get metadata from it - Feng
     gdfid = gdopen(imagefile, DFACC_READ)
     write(*, '(a,I5,a)') "===============", i, "================"

     if (gdfid .NE. -1) then
        gdid = gdattach(gdfid, "MOD_Grid_MOD15A2")    
        if (gdid .NE. -1) then
           ! get grid info 
           status = gdgridinfo(gdid, ncols, nrows, upleft, lowright);
           write(*, '(2a)') "HDF file   = ", imagefile
           write(*, '(a,2I5)') "image size = ", nrows, ncols
           write(*, '(a,2f12.3)') "upperleft  = ", upleft(1), upleft(2)
           status = gddetach(gdid)
        endif
        status = gdclose(gdfid)    
     else
        write(*, '(10a)') "Can't Open File ", imagefile
        write(*, *)
        lai_id(i) = -1
        laiQC_id(i) = -1
     endif

     sdid = sfstart(imagefile, DFACC_READ)
     if (sdid .NE. -1) then

        ! open Lai_1km SDS layer for reading - Feng
        index=sfn2index(sdid, "Lai_1km")
        lai_id(i) = sfselect(sdid, index);

        att_id = sffattr(lai_id(i), "_FillValue")
        if(att_id .NE. -1) then
           status = sfrnatt(lai_id(i), att_id, fillValue)
!           write(*, '(a, I3)') "Fill Value = ", fillValue
        endif

        att_id = sffattr(lai_id(i), "scale_factor")
        if(att_id .NE. -1) then
           status = sfrnatt(lai_id(i), att_id, scale_factor)
!           write(*, '(a, f6.4)') "Scale Factor = ", scale_factor
        endif

        att_id = sffattr(lai_id(i), "valid_range")
        if(att_id .NE. -1) then
           status = sfrnatt(lai_id(i), att_id, range)
!           write(*, '(a, 2I4)') "Valid Range=", range(1), range(2)
        endif

        if(lai_id(i) .NE. -1) then
           write(*, '(2a)') 'LAI data  layer opened!'   
        endif

        ! open fpar layer
        index=sfn2index(sdid, "Fpar_1km")
        fpar_id(i) = sfselect(sdid, index);
        if(fpar_id(i) .NE. -1) then
           write(*, '(2a)') 'FPAR data layer opened!'   
        endif
        
        ! use HDF LAI QA layer for mask operation
        index=sfn2index(sdid, "FparLai_QC")
        laiQC_id(i) = sfselect(sdid, index);
        if(laiQC_id(i) .NE. -1) then
           write(*, '(2a)') 'quality   layer opened!'
        endif

        ! use HDF LAI Extra QA layer for mask operation
        index=sfn2index(sdid, "FparExtra_QC")
        laiQC_ext_id(i) = sfselect(sdid, index);
        if(laiQC_ext_id(i) .NE. -1) then
           write(*, '(2a)') 'extra quality layer opened!'
        endif

        write(*,*)

     else
        lai_id(i) = -1
        fpar_id(i) = -1
        laiQC_id(i) = -1
        laiQC_ext_id(i) = -1
     endif

  end do

  if(abs(upleft(2)) > HIGH_LAT_ULY) then
     isHighLatitudeTile = 1
  else
     isHighLatitudeTile = 0
  endif

  close(1)


  !---- Set filetype --------------------------------------------------------------

  !write(*,*) '   Specify the file type for the image files '
  !write(*,*) '   1 = 8-bit unsigned integer (byte)'
  !write(*,*) '   2 = 16-bit signed integer'
  !write(*,*) '   3 = 32-bit real'
  !read(*,*) filetype              
  !write(*,*)

  !select case (filetype)
  !case (1)
  !  nbyte = 1
  !case (2)
  !  nbyte = 2
  !case (3)
  !  nbyte = 4
  !end select      

  !---- Define image position for processing --------------------------------------

  ! process whole tile - Feng
  nr = nrows
  nc = ncols
  nrstart = 1
  nrstop = nrows
  ncstart = 1
  ncstop = ncols

  !write(*,*) '   Number of rows and columns in the image tile ' 
  !read(*,*)    nr, nc
  !write(*,*)                                            
  !write(*,*) '   Define processing area: start row, stop row, start column, stop column'                       
  !read(*,*)    nrstart, nrstop, ncstart, ncstop                  
  !write(*,*)

  !---- Number of years, number of points per year and total number of points -----

  write(*,*) '   Number of years and number of points per year in the time-series'
  read(*,*)    nyear, nptperyear                                 
  write(*,*)

  ! how about the data set only include part of year - Feng
  !npt = nyear*nptperyear

  !---- Define range for sensor data ----------------------------------------------

  ! use values retrieved from HDF file - Feng
  yl = range(1)
  yu = range(2)
  !write(*,*) '   Sensor data values in the range [min,max] are accepted. Data outside'
  !write(*,*) '   this range are assigned weight 0. Give min and max '
  !read(*,*)    yl,yu
  !write(*,*)

  !---- Define mask data ----------------------------------------------------------

  ! need include LAI QA info here and assign weight according to LAI quality - Feng
  !write(*,*) '   Name of the mask data file list: write NONE if mask data is unavailable ' 
  !read(*,*)    namemask
  !if ( (trim(namemask) == 'NONE') .or. (trim(namemask) == 'none') ) then
  !  mask = 0
  !else
  !  mask = 1
  !end if
  ! open LAI QA layer here - Feng

  !---- Open cloudmask file list --------------------------------------------------

  !if (mask == 1) then
  !  open(1,file=namemask,status='old',err=93) 

  !---- Get number of cloudmask imagefiles in the list ----------------------------

  !  read(1,*) npt
  !  if (npt /= nyear*nptperyear) stop 'Wrong number of cloudmask imagefiles'

  !---- Read name of cloudmask imagefiles and open corresponding files ------------

  !  do i = 1,npt

  !    read(1,'(a)') imagefile
  !   open(20000+i,file=imagefile,status='old',access='stream',err=94)                ! alt1 g95 
  !   open(20000+i,file=imagefile,status='old',access='direct',recl=nbyte*nc,err=94)  ! alt2 g95 
  !    open(20000+i,file=imagefile,status='old',form='binary',recordtype='fixed', &    ! f95
  !        access='direct',recl=nbyte*nc,err=94) 
  !  end do 

  !  close(1)
  !end if


  !---- Set boundaries and weights for mask data conversion -----------------------

  ! need new way to interpolate LAI QA - Feng
  !if (mask == 1) then
  !  write(*,*) '   Sensor data values are assigned weights depending on the '
  !  write(*,*) '   corresponding mask data values: '
  !  write(*,*) '   if mask data in the range [a1,b1] then set weight to w1 '
  !  write(*,*) '   if mask data in the range [a2,b2] then set weight to w2 '
  !  write(*,*) '   if mask data in the range [a3,b3] then set weight to w3 '
  !  write(*,*) '   for mask data values outside these ranges weight is set to 0'  
  !  write(*,*) '   Give a1,b1,w1,a2,b2,w2,a3,b3,w3 '
  !  read(*,*)    bl(1),bu(1),weight(1),bl(2),bu(2),weight(2),bl(3),bu(3),weight(3)
  !  write(*,*)
  !end if

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
  !allocate ( ymatrixchar(npt,nc),wmatrixchar(npt,nc) )
  !allocate ( ymatrixint2(npt,nc),wmatrixint2(npt,nc) )
  allocate ( ymatrixreal(npt,nc),wmatrixreal(npt,nc) )
  allocate ( fpar_matrixchar(npt,nc), lai_matrixchar(npt,nc) )
  allocate ( buf(nc), csp(npt), qc(nc), extqc(nc) )
  allocate ( ty(npt), sy(npt), w1(npt) )

  !---- Open output files ---------------------------------------------------------
!  open(13,file='phenologySG_'//trim(adjustl(jobname)),form='unformatted',status='unknown') 
!  open(14,file='phenologyAG_'//trim(adjustl(jobname)),form='unformatted',status='unknown')
!  open(15,file='phenologyDL_'//trim(adjustl(jobname)),form='unformatted',status='unknown') 
!  open(16,file='sensordata_'//trim(adjustl(jobname)),form='unformatted',access='stream',status='unknown') 
  open(17,file='fitSG_'//trim(adjustl(jobname)),form='unformatted',access='stream',status='unknown') 
  open(18,file='fitAG_'//trim(adjustl(jobname)),form='unformatted',access='stream',status='unknown')
  open(19,file='fitDL_'//trim(adjustl(jobname)),form='unformatted',access='stream',status='unknown') 
  open(20,file='input_'//trim(adjustl(jobname))//'.txt',status='unknown')

  ! open files to save data and weigth in BIP format - Feng  
  open(21,file='originalLAI_'//trim(adjustl(jobname)),form='unformatted',access='stream',status='unknown') 
  open(22,file='assignedweight_'//trim(adjustl(jobname)),form='unformatted',access='stream',status='unknown')
  open(23,file='originalFPAR_'//trim(adjustl(jobname)),form='unformatted',access='stream',status='unknown')  

  !---- Write input parameters to file --------------------------------------------

  ! only write input data for batch process - Feng
  write(20,'(a)') namedata
  !write(20,'(10i6)') filetype
  !write(20,'(10i6)') nr,nc
  !write(20,'(10i6)') nrstart,nrstop,ncstart,ncstop
  write(20,'(10i6)') nyear,nptperyear
  !write(20,'(10f11.4)') yl,yu
  !write(20,'(a)') namemask
  !if (mask == 1) write(20,'(10f11.4)') bl(1),bu(1),weight(1),bl(2),bu(2),weight(2), &
  !                                                bl(3),bu(3),weight(3)
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

  do i = 17,19
     write(i) nyear, nptperyear, nrstart, nrstop, ncstart, ncstop 
  end do

  !---- Generate basis functions needed for SEASON and SAVGOL ---------------------

  call basis(npt,nptperyear,f1,f2,yag,ydl)

  !---- Read data from sensordata and cloudmask imagefiles ------------------------
  !     Loop over rows, for each row loop over all npt imagefiles                 

  write(*,*)
  write(*,*) '  Processing started'
  write(*,*)

  do i = nrstart,nrstop        ! f95
     !do i = 1,nrstop             ! g95
     !if(i .NE. debug_irow) cycle
     !---- Read time-series for all pixels (columns) in the row ----------------------   
     ! read from HDF file - uses column * row in Fortran Feng
!if(i < 422) cycle
     start(1) = 0
     start(2) = i-1
     stride(1) = 1
     stride(2) = 1
     length(1) = ncols
     length(2) = 1

     do j = 1,npt

        ! read LAI 
        status = sfrdata(lai_id(j), start, stride, length, buf)
        do k = 1, ncols
           if(status .NE. -1)then
              ymatrixreal(j,k) = buf(k)
              lai_matrixchar(j,k) = buf(k)
           else
              ymatrixreal(j,k) = fillValue
              lai_matrixchar(j,k) = fillValue
              wmatrixreal(j,k) = 0.0
           endif
        end do

        ! read FPAR
        status = sfrdata(fpar_id(j), start, stride, length, buf)
        do k = 1, ncols
           if(status .NE. -1)then
              fpar_matrixchar(j,k) = buf(k)
           else
              fpar_matrixchar(j,k) = fillValue
           endif
        end do

        status = sfrdata(laiQC_id(j), start, stride, length, qc)
        status = sfrdata(laiQC_ext_id(j), start, stride, length, extqc)

        ! process each pixel
        do k = 1, ncols
           ! interpret QA bits and assign weight here (MOD15A)
           ! "MODLAND_QC START 0 END 1 VALIDS 4\n",
           ! "MODLAND_QC   00=0 Best Possible\n",
           ! "MODLAND_QC   01=1 OK, but not the best\n",
           ! "MODLAND_QC   10=2 Not produced,due to cloud\n",
           ! "MODLAND_QC   11=3 Not produced,due to other reasons\n",
           ! "SCF_QC START 5 END 7 VALIDS 4\n",
           ! "SCF_QC       000=0 Main (RT) method used with the best possible results\n",    
           ! "SCF_QC       001=1 Main (RT) method used with saturation\n",
           ! "SCF_QC       010=2 Main (RT) method failed due to geometry problems, empirical method used\n",
           ! "SCF_QC       011=3 Main (RT) method failed due to problems other than geometry, emprical method used\n",
           ! "SCF_QC       100=4 Couldn\'t retrieve pixel\n",
           ! "SCF_QC       111=7 NOT PRODUCED AT ALL (non-terrestrial biome)\n",

           !!! MCD15A3 QC flags
           ! "MODLAND_QC START 0 END 0 VALIDS 2\n",
           ! "MODLAND_QC   0 = Good Quality (main algorithm with or without saturation)\n",
           ! "MODLAND_QC   1 = Other Quality (back-up algorithm or fill value)\n",
           ! "SENSOR START 1 END 1 VALIDS 2\n",
           ! "SENSOR       0  = Terra\n",
           ! "SENSOR       1  = Aqua\n",
           ! "DEADDETECTOR START 2 END 2 VALIDS 2\n",
           ! "DEADDETECTOR 0 = Detectors apparently fine for up to 50% of channels 1,2\n",
           ! "DEADDETECTOR 1 = Dead detectors caused >50% adjacent detector retrieval\n",
           ! "CLOUDSTATE START 3 END 4 VALIDS 4 (this inherited from Aggregate_QC bits {0,1} cloud state)\n",
           ! "CLOUDSTATE   00 = 0 Significant clouds NOT present (clear)\n",
           ! "CLOUDSTATE   01 = 1 Significant clouds WERE present\n",
           ! "CLOUDSTATE   10 = 2 Mixed cloud present on pixel\n",
           ! "CLOUDSTATE   11 = 3 Cloud state not defined,assumed clear\n",
           ! "SCF_QC START 5 END 7 VALIDS 5\n",
           ! "SCF_QC       000=0 Main (RT) algorithm used, best result possible (no saturation)\n",
           ! "SCF_QC       001=1 Main (RT) algorithm used, saturation occured. Good, very usable.\n",
           ! "SCF_QC       010=2 Main algorithm failed due to bad geometry, empirical algorithm used\n",
           ! "SCF_QC       011=3 Main algorithm faild due to problems other than geometry, empirical algorithm used\n",
           ! "SCF_QC       100=4 Pixel not produced at all, value coudn\'t be retrieved 

           !wmatrixchar(j,k) = buf(k)

           if(status .NE. -1) then
              !bits = iand(ishft(buf(k),[-7:0]),1) ! unpack word
              !tmpbits(0) = bits(0)  ! bit 7
              !tmpbits(1) = bits(1)  ! bit 6
              !tmpbits(2) = bits(2)  ! bit 5
              !SCF_QC = sum(ishft(tmpbits,[2:0:-1])) ! pack bits
              !           if(k==debug_icol .AND. i==debug_irow) then
              !              write(*,'(I5, I5, a, 8z1)') k, buf(k)," = ", bits
              !              write(*,'(3z1, a, I5)') tmpbits, " = ", SCF_QC
              !           endif

              ! change for gfortran - Feng (09/12)
              !tmpbits(7) = iand(ishft(buf(k), -7), 1) ! bit 7
              !tmpbits(6) = iand(ishft(buf(k), -6), 1) ! bit 6
              !tmpbits(5) = iand(ishft(buf(k), -5), 1) ! bit 5
              !tmpbits(4) = iand(ishft(buf(k), -4), 1) ! bit 4

              do ii = 0, 7
                 bits(ii) = iand(ishft(qc(k), -ii), 1)
              end do
              SCF_QC = bits(5)+bits(6)*2+bits(7)*4
              CLD_QC = bits(3)+bits(4)*2
	      DET_QC = bits(2)

             ! if(k==debug_icol .AND. i==debug_irow) then
             !    write(*,'(I5, I5, 8z1, I5, I5, I5)') k, buf(k), bits, SCF_QC, CLD_QC, DET_QC
             ! endif

              select case (SCF_QC)
              case(0)      
                 wmatrixreal(j,k) = 1.0
              case(1)
                 wmatrixreal(j,k) = 1.0
              case(2)
                 wmatrixreal(j,k) = 0.25
              case(3)
                 wmatrixreal(j,k) = 0.25
              case(4)
                 ! for some high latitude tiles (50 degrees and higher) 
                 ! includes v00, v01, v02, v03, v14, v15, h16, h17 
                 ! if it is not produced due to large solar zenith angle in winter
                 ! we will use LAI=0.0 and w=0.1 to feed to timesat
                 ! so it will "eliminate" missing problem and force timesat to run 
                 if(isHighLatitudeTile .EQ. 1) then
                    wmatrixreal(j,k) = 0.25
                    ymatrixreal(j,k) = 0.0
                    ! use original MODIS value for writing
                    ! lai_matrixchar(j,k) = 0
                    ! fpar_matrixchar(j,k) = 0
                 else
                    wmatrixreal(j,k) = 0.0
                    ymatrixreal(j,k) = fillValue
                 end if
              case default
                 wmatrixreal(j,k) = 0.0
              end select

              ! need to check extra QC flags to ensure data quality
              do ii = 0, 7
                 bits(ii) = iand(ishft(extqc(k), -ii), 1)
              end do
              !    "FparExtra_QC 6 BITFIELDS IN 8 BITWORD\n",
              !    "LANDSEA PASS-THROUGH START 0 END 1 VALIDS 4\n",
              !    "LANDSEA   00 = 0 LAND       AggrQC(3,5)values{001}\n",
              !    "LANDSEA   01 = 1 SHORE      AggrQC(3,5)values{000,010,100}\n",
              !    "LANDSEA   10 = 2 FRESHWATER AggrQC(3,5)values{011,101}\n",
              !    "LANDSEA   11 = 3 OCEAN      AggrQC(3,5)values{110,111}\n",
              !    "SNOW_ICE (from Aggregate_QC bits) START 2 END 2 VALIDS 2\n",
              !    "SNOW_ICE  0 = No snow/ice detected\n",
              !    "SNOW_ICE  1 = Snow/ice were detected\n",
              !    "AEROSOL START 3 END 3 VALIDS 2\n",
              !    "AEROSOL   0 = No or low atmospheric aerosol levels detected\n",
              !    "AEROSOL   1 = Average or high aerosol levels detected\n",
              !    "CIRRUS (from Aggregate_QC bits {8,9} ) START 4 END 4 VALIDS 2\n",
              !    "CIRRUS    0 = No cirrus detected\n",
              !    "CIRRUS    1 = Cirrus was detected\n",
              !    "INTERNAL_CLOUD_MASK START 5 END 5 VALIDS 2\n",
              !    "INTERNAL_CLOUD_MASK 0 = No clouds\n",
              !    "INTERNAL_CLOUD_MASK 1 = Clouds were detected\n",
              !    "CLOUD_SHADOW START 6 END 6 VALIDS 2\n",
              !    "CLOUD_SHADOW        0 = No cloud shadow detected\n",
              !    "CLOUD_SHADOW        1 = Cloud shadow detected\n",
              !    "SCF_BIOME_MASK START 7 END 7 VALIDS 2\n",
              !    "SCF_BIOME_MASK  0 = Biome outside interval <1,4>\n",
              !    "SCF_BIOME_MASK  1 = Biome in interval <1,4>" ;

              SCF_QC = bits(2)+bits(3)+bits(4)+bits(5)+bits(6)
              if((SCF_QC .NE. 0).and.(wmatrixreal(j,k) .LT. 0.99)) then
                 wmatrixreal(j,k) = 0.25   ! use zero will leave too many missings? (10/12)
              endif

              ! set high weight for value larger than saturation point
              if (ymatrixreal(j,k) .GE. MIN_SAT_POINT) then
                 wmatrixreal(j,k) = 1.0
              endif
              ! However, if other bits say not good, we are better to be more conservative
              !if (bits(6)==1 .OR. bits(5)==1 .OR. bits(4)==1 .OR. bits(3)==1) then
              if(CLD_QC==2) then
                 wmatrixreal(j,k) = 0.25
              endif
	      if(CLD_QC==1 .or. DET_QC==1) then
                 wmatrixreal(j,k) = 0.0
              end if
!if(i==debug_irow .and. k==debug_icol) then
!   write(*, '(I6, I5, "=", 8Z1, I3, F5.2, F6.2)') k, buf(k), bits, SCF_QC,  wmatrixreal(j,k), ymatrixreal(j,k)
!endif

           else
              ymatrixreal(j,k) = fillValue
              lai_matrixchar(j,k) = fillValue
              fpar_matrixchar(j,k) = fillValue
              wmatrixreal(j,k) = 0.0
           endif

        end do
     end do

     !  select case (filetype)
     !  case (1)
     !    do j = 1,npt
     !      read(10000+j,rec=i) ymatrixchar(j,:)
     !      if (mask == 1) read(20000+j,rec=i) wmatrixchar(j,:)
     !    end do
     !  case (2)
     !    do j = 1,npt
     !      read(10000+j,rec=i) ymatrixint2(j,:)
     !      if (mask == 1) read(20000+j,rec=i) wmatrixint2(j,:)
     !    end do
     !  case (3)
     !    do j = 1,npt
     !      read(10000+j,rec=i) ymatrixreal(j,:)
     !      if (mask == 1) read(20000+j,rec=i) wmatrixreal(j,:)
     !    end do
     !  end select

     !  if (i < nrstart) cycle    ! g95

     write(*,*) '  Row ',i

     do j = ncstart,ncstop

        !---- Transfer time-series of sensor and cloudmask data to vektors y and w ------
        !if(j .NE. debug_icol) cycle

        !  convert for LAI only - Feng
        y = dble( (ymatrixreal(:,j)) )
        w = dble( wmatrixreal(:,j) )

        !---- Check if sensordata is in the specified range, if not assign weight zero --

        where ( (y < yl) .or. (y > yu) ) w = 0.d0

        ! write original time series data and weight in binary format
        csp = lai_matrixchar(:,j)
	write(21) csp
        csp = fpar_matrixchar(:,j)
        write(23) csp
        ysp = w	
        write(22) ysp

        !---- Time-series with too many data values with weight zero are not processed --

        missing = 0

        if (count(w == 0.d0) >= 2.d0*npt/4.d0) missing = 1
        do k = 1,npt-nptperyear/4
           if ( abs(sum(w(k:k+nptperyear/4))) == 0.d0 ) missing = 1
        end do

        if (missing == 1) cycle  

        !---- Only time-series with non-constant datavalues are processed ---------------

        do k = 1,npt-1
           diff(k) = abs(y(k) - y(k+1)) 
        end do
        ! remove constant value checking - Feng
        !if (count(diff < 1.d-6) >= npt/4 ) cycle      
        if (count(diff < 1.d-6) >= 3*npt/4 ) cycle   ! PIETER      
        if (count(diff < 1.d-6) >= 3.d0*npt/4.d0 ) cycle  ! Duarte 

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
!        ysp = y
!	write(16) i,j
!	write(16) ysp

        !---- Iterative Savitzky-Golay filtering, analyze and dump phenological data ----  

        if (sg == 1) then
           call savgol(y,w,f2,y1,maxval(win))
         !  call phenology(13,i,j,s,y,y1,nlocal,startcutoff,minimum)
           where (y1 < yl) y1 = yl
           where (y1 > yu) y1 = yu                         
           ysp = y1
           write(17) i,j
           write(17) ysp
        end if

        !---- Iterative fits to asymmetric Gaussians, analyze and dump phenological data -      

        if (ag == 1) then      

           call fitgauss(s,yag,y,w,nlocal,y2)
           ! check if fitted values are valid number 
           flag = 1
           do k = 1, npt 
              !    if(isnan(y2(k))) then 
              if(y2(k) .NE. y2(k)) then
                 flag = 0
                 exit
              endif
           end do

           if(flag == 1) then
              where (y2 < yl) y2 = yl
              where (y2 > yu) y2 = yu                         
              ysp = y2
              write(18) i,j
              write(18) ysp
              ! call phenology(14,i,j,s,y,y2,nlocal,startcutoff,minimum)
           end if
           
        end if

!if(j==debug_icol .AND. i==debug_irow) then
!do k = 1, npt
!   write(*,'(I5, I3, 3f10.2)') j, k, y(k)/10.0, w(k), ysp(k)/10.0
!end do
!endif
        !---- Iterative fits to double logistic functions, analyze and dump phenological data            

	if (dl == 1) then  
           call fitlogistic(s,ydl,y,w,nlocal,y3)
           where (y3 < yl) y3 = yl
           where (y3 > yu) y3 = yu
        !   call phenology(15,i,j,s,y,y3,nlocal,startcutoff,minimum)
           ysp = y3
           write(19) i,j
           write(19) ysp
	end if

     end do ! end loop over columns 
  end do ! end loop over rows

  !---- Close all imagefiles ------------------------------------------------------

  !do i = 1,npt
  !    close(10000+i)
  !    close(20000+i)
  !end do

  !---- Close remaining files -----------------------------------------------------

  do i = 17,23
     close(i)
  end do

  !----- Print information to screen ----------------------------------------------

  write(*,*)
  write(*,*) '  Processing finished  ' 
  write(*,*) 
!  write(*,*) '  Phenological parameters written to:'
!  if (sg == 1) then 
!     write(*,*) '  phenologySG_'//trim(adjustl(jobname))
!  end if
!  if (ag == 1) then 
!     write(*,*) '  phenologyAG_'//trim(adjustl(jobname))
!  end if
!  if (dl == 1) then 
!     write(*,*) '  phenologyDL_'//trim(adjustl(jobname))
!  end if
!  write(*,*) 
  write(*,*) '  Sensor data and fitted functions written to:'
! write(*,*) '  sensordata_'//trim(adjustl(jobname))
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

  stop

91 stop ' Error opening sensordatalist'
92 stop ' Error opening sensordata imagefiles'  
93 stop ' Error opening maskdatalist'
94 stop ' Error opening cloudmask imagefiles'     

end program timesatimage
