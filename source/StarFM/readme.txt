Description (v1.2.1)
====================

MODIS and Landsat surface reflectance products have complementary characteristics 
in terms of spatial and temporal resolution. To fully exploit these individual 
datasets, a Spatial and Temporal Adaptive Reflectance Fusion Model (STARFM) was 
developed (Gao et al., 2006). It combines the spatial resolution of Landsat with 
the temporal frequency of MODIS. STARFM uses comparisons of one or more pairs of 
observed Landsat/MODIS maps, collected on the same day, to predict maps at Landsat-scale 
on other MODIS observation dates. This program was revised from the original version
of STARFM package (v1.1.2) and allows multiple MODIS dates for the prediction. The 
program enables OpenMP parallel computing that allows program running on multiple cores. 

The program was tested in Linux using gcc compiler (v4.4 or newer to support OpenMP 3.0). 

Reference
  
Gao, F., J. Masek, M. Schwaller and H. Forrest, (2006). On the Blending of the 
Landsat and MODIS Surface Reflectance: Predict Daily Landsat Surface Reflectance, 
IEEE Transactions on Geoscience and Remote Sensing, Vol. 44, No. 8, pp. 2207-2218.


Installation
============

gzip -d StarFM.tar.gz
tar -xvf StarFM.tar

StarFM/source
StarFM.h
StarFM_main.c
StarFM_compute.c
StarFM_util.c
StarFM_alloc.c
Makefile

StarFM/simulation_test
StarFM/simulation_test/change_veg  
StarFM/simulation_test/road_veg  
StarFM/simulation_test/water_veg
(include simulated testing files)

StarFM/reflectance_test (full version only)
include MODIS and Landsat testing files)

then run 
> make 
under source directory and create executive program "StarFM.exe"


Usage 
=====

> StarFM.exe input.txt
where input.txt contains the input files and parameters (see details below)


Inputs
======

All Landsat and MODIS data should be first pre-processed and co-registered and saved in 
  - same resolution   (Landsat resolution)
  - same image size   
  - same map projection
  - same image extent (orthorectifed and precisely registered)
  - same scale factor (10000 for MODIS reflectance products)
  - "short int" type (2 byte/pixel) for the scaled reflectance data
  - "unsigned char" type (1 byte/pixel) for optional mask and classification file 

Major inputs include:
  - Landsat and MODIS surface reflectance image pair in binary format
  - MODIS surface reflectance for the prediction date
  - optional Landsat classification map for input pair 
  - optional Landsat and MODIS cloud and poor quality data mask
  - control parameters 

Detailed input parameters and explanations:

STARFM_PARAMETER_START
# number of input pairs, maximum = 2
        NUM_IN_PAIRS = 
# input MODIS files
        IN_PAIR_MODIS_FNAME = 
# optional mask file, same order MODIS files and use NONE if not exist
        IN_PAIR_MODIS_MASK = 
# input Landsat files in the same order as MODIS
        IN_PAIR_LANDSAT_FNAME = 
# optional mask file for Landsat
        IN_PAIR_LANDSAT_MASK = 
# optional classification map in the same order as input pairs
# classification map has to be generated from input Landsat data 
# with as many separable classes as possible
        IN_PAIR_CLASSIFICATION_MAP =
# MODIS input for the prediction date 
        IN_PDAY_MODIS_FNAME = (accept multiple MODIS dates)
# optional mask file for MODIS
        IN_PDAY_MODIS_MASK = 
# output Landsat prediction file
        OUT_PDAY_LANDSAT_FNAME = (include same number of MODIS prediction dates) 
# number of rows 
        NROWS =
# number of columns
        NCOLS =
# spatial resolution
        RESOLUTION =
# scale factor for input reflectance file 
        SCALE_FACTOR =
# fill value for Landsat surface reflectance
        LANDSAT_FILLV =
# Landsat data range
        LANDSAT_DATA_RANGE =
# uncertainty in the same unit as the scaled reflectance
        LANDSAT_UNCERTAINTY = 
# fill value for MODIS surface reflectance
        MODIS_FILLV =
# MODIS data range
        MODIS_DATA_RANGE =
# uncertainty for MODIS 
        MODIS_UNCERTAINTY =
# spatial information flag, "ON" is strongly suggested
        USE_SPATIAL_FLAG = ON
# maximum search distance for the spectral similar neighbor pixels 
        MAX_SEARCH_DISTANCE =
# number of slice for the spectral similar test (pure pixel)
        NUM_SLICE_PURE_TEST =
STARFM_PARAMETER_END

(use # for comments)

MODIS inputs can choose daily surface reflectance products (MOD09) or nadir BRDF-adjusted 
reflectance (NBAR) (MOD43). They can be reprojected and resampled to Landsat projection 
and resolution using MODIS Reprojection Tool (MRT). The MRT can be downloaded from 
http://edcdaac.usgs.gov/landdaac/tools/modis/index.asp 


Output
======

- Landsat surface reflectance for the prediction date in binary format. 
  It also comes with an associated ENVI header file.
- samples used for prediction within searching window in basis point (unit: 1/10000)
  valid data range is 0 - 10000 (high values normally represent better prediction)
  32767 means the replacement from the best prediction look-up-table


Testing Data
============

1) simulated data fusion

- StarFM/simulation_test/water_veg (see Fig. 2 in Gao et al.)

"input_use_spatial.txt" 
use t1 image pair to predict reflectance for t2 using spatial information. 

"input_no_spatial.txt"  
use t1 image pair to predict reflectance for t2 without using spatial information. 

This example shows the difference with and without using spatial information. 
The "USE_SPATIAL_FLAG" is set to "ON" for the rest of tests. We strongly recommend
to use spatial information in the prediction.  

- StarFM/simulation_test/change_veg (see Fig. 3 in Gao et al.)

"input_t2.txt" 
uses t1 and t4 pairs to predict t2

"input_t3.txt" 
uses t1 and t4 pairs to predict t3

This example shows how StarFM algorithm handles changing objects.

- StarFM/simulation_test/road_veg (see Fig. 5 in Gao et al.)

"input_t2.txt" 
uses t1 and t4 image pairs to predict t2 
(with replacement option on, for MIN_SAMPLE_PERCENT <= 2% )

"input_t3.txt" 
uses t1 and t4 image pairs to predict t3

This example shows how StarFM handles linear objects.

2) Actual Landsat data fusion 
(full version only, compact version need to download separately)

"input.06-04-01.*.txt"
use 05-24-01 image pair to predict 06-04-01 Landsat surface reflectance
(used default slice approach for spectral similar testing)

"input.07-11-01.*.txt"
use 05-24-01 image pair to predict 07-11-01 Landsat surface reflectance
(used defined classification map for spectral similar testing)

"input.08-12-01.*.txt"
use 07-11-01 image pair to predict 08-12-01 Landsat surface reflectance
(used mask file to exclude low quality data or clouds)
(with replacement option on, for MIN_SAMPLE_PERCENT <= 1%)

"input2.07-11-01.*.txt"
use 05-24-01 and 08-12-01 image pairs to predict 07-11-01
(with replacement option on, for MIN_SAMPLE_PERCENT <= 2%)


Appendix A
==========

1) MODIS data preparation
MODIS inputs can be either daily surface reflectance (MOD09) or 16-day nadir BRDF-adjusted
surface reflectance (NBAR) product (MOD43B4 or MCD43A4). They need to be pre-processed
for STARFM inputs.
a) order MODIS tiles that cover whole Landsat scene on the Landsat acquisition date and
prediction dates
b) use MODIS reprojection tool (MRT) to mosaic and reproject MODIS tiles to Landsat 
projection, resolution and extents. Save otuputs in binary format.

2) Landsat data preparation
Landsat DN values need to be calibrated and atmosphericly corrected and saved in the
same format as MODIS surface reflectance (2-byte integer with scale factor of 10000) 
with one file per band. 
Note that MODIS bands use different band number. You can relate them with 
Landsat  MODIS
   1       3
   2       4
   3       1
   4       2
   5       6
   7       7

 





