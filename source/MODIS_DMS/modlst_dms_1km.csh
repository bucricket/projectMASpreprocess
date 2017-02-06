#! /bin/csh -f

# sharpening off-nadir daily MODIS LST (1km) using MODIS NDVI composite (1km)
# energy conservation is applied at 2km  
  
# v1.0 original version for Landsat TIR imagery
# v1.1 modified for MODIS TIR sharpening using MODIS 8-day NDVI 
# v1.2 combines previous revisions from 201405 and 201408 (10/29/2015)
# v1.3 modified sharpening according to view zenith angle (11/11/2015)
# v1.4 modified script for processing multiple LST files (yyyydddtttt) from one day (12/9/2015) 

if $#argv != 1 then
    echo "Usage:   modlst_dms_1km.csh <yyyyddd> or modlst_dms_1km.csh <yyyydddtttt>"
    echo "Example: modlst_dms_1km.csh 2010171 or modlst_dms_1km.csh 2010171750"
    exit
else
    set doyt = $argv[1]
    set doy = `echo $doyt | cut -c1-7`
endif

# define program directory 
set exepath = $modis_sharpen

# define DMS executable programs
set tsamples = "$exepath/get_samples.exe"
set predict = "$exepath/predict_fineT.exe"
set combine = "$exepath/combine_models.exe"
set convert_lst = "$exepath/modlst_C2K.exe"
set convert_ndvi = "$exepath/modndvi_F2I.exe"

# define cubist regression tree program
set cubist = "/stor/array01/tools/Cubist/cubist"
set nrules = 1
set min_samples = 50

# define input and output directories
set ndvi_dir = $MODIS_NDVI_DIR
set lst_dir = $MODIS_LST_DIR
set out_dir = $MODIS_LST_SHARP_OUT
if(! $?MODIS_ANGLE_DIR) then
    set angle_dir = `echo $ndvi_dir | sed 's/NDVI/view_angle/'`
else    
    set angle_dir = $MODIS_ANGLE_DIR
endif

if ! -e $out_dir then
  mkdir $out_dir
endif
cd $out_dir

# define maximum numbers of samples (1M)
set MAX_NSAMS = 1000000

# find all lndsr file under Landsat directory
set modis_ndvi = `ls $ndvi_dir/NDVI_$doy*`
set modis_lst = `ls $lst_dir/LST_$doyt*`
set hdr = `ls $ndvi_dir/NDVI_$doy*.HDR`
set view_angle = `ls $angle_dir/theta_$doyt*`

echo "$ndvi_dir"
if ($#modis_ndvi == 0) then
   echo "MODIS NDVI on $doy does not exist!"
   exit
endif

if ($#modis_lst == 0) then
   echo "MODIS LST on $doyt does not exist!"
   exit
endif

if ($#view_angle == 0) then
   echo "MODIS view angle on $doyt does not exist, do sharpening anyway"
endif
 
# extract image info from lndsr file 
set nrows = `grep -m 1 "lines" $hdr | awk '{print $3}'`
set ncols = `grep -m 1 "samples" $hdr | awk '{print $3}'`
set res = "0.01"
set ul = `grep -m 1 "map info" $hdr | awk '{print $8, $9}' | sed 's/,//g'`
set zone = "-1"

set th_nrows = $nrows
set th_ncols = $ncols
set th_res = "0.02"

set outf = "sharpened_modlst_$doyt.bin"

echo $modis_ndvi[1] $modis_lst[1] $nrows $ncols $ul

# convert thermal data from Celsius degree (DN*0.01)to Kelvin for DMS processing
set modis_lst_k = "modlst_kelvin.$doyt.bin"
$convert_lst $modis_lst[1] $modis_lst_k

# convert NDVI from float to integer
set modis_ndvi_i = "modndvi_int.$doy.bin"
$convert_ndvi $modis_ndvi[1] $modis_ndvi_i

set inp = "dms_$doyt.inp"
set cubistf = "th_samples"

echo "# input file for Data Mining Sharpener" > $inp
echo "NFILES = 1" >> $inp
echo "SW_FILE_NAME = $modis_ndvi_i" >> $inp
echo "SW_CLOUD_MASK = none" >> $inp
echo "SW_FILE_TYPE = binary" >> $inp
echo "SW_CLOUD_TYPE = binary" >> $inp
echo "SW_NROWS = $nrows" >> $inp
echo "SW_NCOLS = $ncols" >> $inp
echo "SW_PIXEL_SIZE = $res" >> $inp
echo "SW_FILL_VALUE = -9999" >> $inp
echo "SW_CLOUD_CODE = 1" >> $inp
echo "SW_DATA_RANGE = 2000, 10000" >> $inp
echo "SW_UPPER_LEFT_CORNER = $ul" >> $inp
echo "SW_PROJECTION_CODE = 0" >> $inp
echo "SW_PROJECTION_PARAMETERS = 6371007.181000 0 0 0 0 0 0 0 0 0 0 0 0 0 0" >> $inp
echo "SW_PROJECTION_ZONE = $zone" >> $inp
echo "SW_PROJECTION_UNIT = 4" >> $inp
echo "SW_PROJECTION_DATUM = 12" >> $inp
 
echo "ORG_TH_FILE_NAME = $modis_lst_k" >> $inp
if ($#view_angle != 0) then
    echo "ORG_TH_VIEW_ANGLE = $view_angle[1]" >> $inp
else 
    echo "ORG_TH_VIEW_ANGLE = none" >> $inp
endif
echo "ORG_TH_FILE_TYPE = BINARY" >> $inp
echo "ORG_TH_DATA_RANGE = 200., 400." >> $inp
echo "ORG_TH_PIXEL_SIZE = $res" >> $inp
echo "ORG_NROWS = $nrows" >> $inp
echo "ORG_NCOLS = $ncols" >> $inp

echo "RES_TH_PIXEL_SIZE = $th_res " >> $inp

echo "PURE_CV_TH = 0.05" >> $inp
echo "SMOOTH_FLAG = 0" >> $inp
echo "CUBIST_FILE_STEM = $cubistf" >> $inp
echo "OUT_FILE = $outf" >> $inp
echo "end" >> $inp

echo "Sharpening $modis_lst[1] using $modis_ndvi[1]"

set log = "tree.log" 
set scale = `echo $th_res $res | awk '{print $1/$2}'`

@ cnrows = $nrows / $scale
@ cncols = $ncols / $scale

# skip if output file exists 
#if (-e $outf) then
#    echo "File $outf has been existed, skip sharpening"
#    rm $inp $modis_lst_k
#    exit
#endif

echo ""
echo "*** Doing Prediction with Global Model ***"
# global model prediction
$tsamples $inp
set nsams = `wc -l $cubistf.data`
set percent = `echo $MAX_NSAMS $nsams[1] | awk '{printf("%d", $1/$2*100)}'`
if($percent > 99) @ percent = 99
if($percent < 1) @ percent = 1
$cubist -f $cubistf -r $nrules -S $percent -u > $log
#$cubist -r 10 -f $cubistf > $log
echo "total # of samples =" $nsams[1]  "; used samples =" $percent "%"
if ($nsams[1] < $min_samples) then
   echo "not enough samples to build regression trees, skip"
   rm $inp $log th_samples.* $modis_lst_k $modis_ndvi_i
   exit
endif
$predict $inp 
set global_out=`echo $outf| sed 's/.bin/.global/'`
echo "global model results are saved in: " $global_out
mv $outf $global_out
cp $outf.hdr $global_out.hdr

echo ""
echo "*** Doing Prediction with Local Model ***"
# local model prediction
@ s_row = 0
@ s_col = 0
#@ wsize1 = 200
#@ overlap1 = 50
# use variant moving windoe size (50 to 100) decided by image size (or in 100-200km)
set wsize1 = `echo $nrows $ncols $scale | awk '{if($1>$2) printf("%d", $1/$3/10); else printf("%d", $2/$3/10);}'`
if($wsize1 < 50) @ wsize1 = 50
if($wsize1 > 100) @ wsize1 = 100 
set overlap1 = `echo $wsize1 | awk '{printf("%d", $1/4)}'`
# use same size of sub-area as TM 
#  set wsize = `echo $wsize1 $th_res | awk '{printf("%d", $1*120./$2)}'` 
#  set overlap = `echo $overlap1 $th_res | awk '{printf("%d", $1*120./$2)}'` 
set wsize = $wsize1
set overlap = $overlap1
echo $wsize $overlap   

while($s_row<$cnrows)
    while($s_col<$cncols) 
      @ e_row = $s_row + $wsize
      @ e_col = $s_col + $wsize
      echo ""
      echo "--- processing central  window" $s_row $s_col $e_row $e_col
      @ os_row = $s_row - $overlap
      @ os_col = $s_col - $overlap
      @ oe_row = $e_row + $overlap
      @ oe_col = $e_col + $overlap
      echo "    overlapped sampling window" $os_row $os_col $oe_row $oe_col
      $tsamples $inp $os_row $os_col $oe_row $oe_col
      set nsamples = `wc -l $cubistf.data | awk '{print $1}'`
      if($nsamples > $min_samples) then 
        # do cubist prediction
        $cubist -f $cubistf -r $nrules -u >> $log
        $predict $inp $s_row $s_col $e_row $e_col
        # cp $outf $outf.$s_row$s_col
      endif
      @ s_col += $wsize
    end
    @ s_col = 0
    @ s_row += $wsize
end

set local_out=`echo $outf| sed 's/.bin/.local/'`
echo "local model results are saved in:" $local_out
mv $outf $local_out
cp $outf.hdr $local_out.hdr

echo ""
echo "*** Combing Local and Gloal Model ***"
#  replace global prediction with local prediction if errors are smaller
echo "Combining" $global_out  "and"  $local_out  "to"  $outf
$combine $inp

grep "map info" $hdr >> $outf.hdr
grep "coordinate" $hdr >> $outf.hdr

rm $inp th_samples.* *.global* *.local* $modis_lst_k $modis_ndvi_i tree.log




