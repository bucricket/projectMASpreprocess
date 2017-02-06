================================================
Landsat TIR sharpening Standalone Version (v2.0)
================================================ 

This package was developed for sharpening Landsat thermal imagery using SWIR bands.
Current version supports Landsat 5, and 8 surface reflectance (SR) data products generated 
by the USGS EROS data center. LEDAPS processing has been removed from the package. The Cubist 
software is needed for building regression tree and should be downloaded and installed separately. 

[REF]
Gao, F., Kustas, W.P., Anderson, M.C. (2012). A Data Mining Approach for Sharpening Satellite Thermal 
Imagery over Land. Remote Sensing, 4(11), 3287-3319, doi:10.3390/rs4113287.

Version History
- v2.0 by Feng Gao on March 2015 
- v1.0 by Feng Gao on January 17, 2014

***
I. Data Preparation
***

- make a new directory for the test (e.g. test)

> mkdir test
> cd test
> mkdir landsat
> cd landsat
> mkdir org_tar

- order Landsat SR data product through the USGS Earth Explorer (earthexplorer.usgs.gov) 
or ESPA (espa.cr.usgs.gov) websites. Landsat SR data ordered through Earth Explorer are 
saved in GeoTiff format. Data ordered through ESPA can be saved in GeoTiff, ENVI and HDF format. 
ESPA provides bulk order option and includes additional processing for TOA reflectance, SR and 
Brightness Temperature. The list of Landsat scenes required by the ESPA website can be created 
using the USGS GloVis or Earth Explorer websites. 

- download Landsat SR data when they are ready. Bulk downloads are possible by using 
the DownThemAll FireFox plugin (needs to be first installed in the FireFox Web browser),
then save all Landsat SR data in /landsat/org_tar

- create a text file (e.g. env.sh) in working directory to save path environment as

> more env.sh
export BIN="/home/fgao/Tools/Landsat_DMS/bin"
export PATH=$BIN:$PATH

export DATA=`pwd`
export LANDSAT_DIR=$DATA/landsat
export LANDSAT_LST_DIR=$DATA/landsat_lst
export LANDSAT_TAR_DIR=$DATA/landsat/org_tar

===
II. Usage

> source env.sh
- setup program and data directories 

> lndsr2envi_tir.csh
- convert Landsat SR and TIR from GeoTIFF and HDF file to ENVI file and save in $LANDSAT_DIR by 
  calling GeoTiff2ENVI.exe and HDF commands
- copy Landsat ENVI file to $LANDSAT_DIR

> lndlst_dms_sa.csh
- run Landsat thermal sharpening program automatically 
