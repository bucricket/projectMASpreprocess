/**
 * header file for data mining sharpening (DMS) 
 * written by Feng Gao on June 2012
 * v1.0_HDF: uses LEDAPS SR in HDF format (in-house LEDAPS product)
 * v1.1_ENVI: accepts GeoTIFF, ENVI and HDF format (USGS Landsat SR product)
 * v1.2_ENVI: uses Landsat cloud mask (USGS Landsat fmask product)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>

/* header files for GeoTIFF */
#include "geotiffio.h"
#include "xtiffio.h"
#include "tiffio.h"

/* header file for GNU Scientific Libraries */
//#include <gsl/gsl_multifit.h>

#define MIN_NSAMPLES 20
#define MAX_STRLEN   1000
#define MAX_NBANDS   10
#define EC_VALID_TH  0.5

#define SUCCESS 0
#define FAILURE -1
#define FILLV   -9999   


/* define structure to hold all variables from input */
typedef struct {

  /* for original shortwave spectral bands */
  int    nbands;
  char   inFile[MAX_NBANDS][MAX_STRLEN];
  char   cloudFile[MAX_STRLEN];
  char   fileType[MAX_STRLEN];   
  int    nrows;
  int    ncols;
  int    fillv;
  int    range[2];
  float  res;

  /* geolocation information (same for shortwave spectral bands and thermal band */
  float  ulx, uly;
  long   insys;
  long   inzone;
  double inparm[15];
  long   inunit;
  long   indatum;

  /* original temperature (float for binary file or byte for GeoTIFF) */
  char  org_th_File[MAX_STRLEN];           
  char  org_th_fileType[MAX_STRLEN];   
  int   org_th_nrows;
  int   org_th_ncols;
  float org_th_res;
  float org_th_range[2];
   
  /* reduced and rescaled T in coarse resolution */
  int   res_th_nrows;
  int   res_th_ncols;
  float res_th_res;

  char  cubistFile[MAX_STRLEN]; /* filestem for cubist */
  char  outFile[MAX_STRLEN];    /* sharpened T (output) */

  float CV_TH;        /* coeficient of variation threshold */
  int   SMOOTH_FLAG;  /* 1=apply smooth operation to difference map in EC application; 0=do not apply */ 
 
  /* processing window */
  int s_row;  /* start row */
  int s_col;  /* start column */
  int e_row;  /* end row */
  int e_col;  /* end column */

} INPUT_PARS;  


/* define structure to store data (shortwave bands or thermal band ) */
typedef struct {
  int    nbands;
  char   fileName[MAX_NBANDS][MAX_STRLEN];
  char   fileType[MAX_STRLEN];   
  int    nrows;           
  int    ncols;
  int    fillValue;
  float  range[2];
  long   total_valid_data;
  double ulx;              /* x value of up-left corner */
  double uly;              /* y value of up-left corner */
  double res;              /* spatial resolution */

  long   insys;
  long   inzone;
  double inparm[15];
  long   inunit;
  long   indatum;
  float  homoTest[MAX_NBANDS];   /* homogeneous test (stdev/mean) */

  FILE   *fp[MAX_NBANDS];        /* input file pointer for binary file */
  TIFF   *fp_tiff[MAX_NBANDS];   /* input file pointer for GeoTIFF file */
  int16  **data;                 /* one row in [iband][icol] for shortwave spectral bands */
  float  **fdata;                /* whole image in [irow][icol] for thermal band in float point */
  int8   **qa;                   /* qa data for output (not used) */

  FILE   *fp_cloud;              /* for cloud mask (added on 10/2015) */
  TIFF   *fp_tiff_cloud;
  uint8  **cloud;
  int    CMASK_FLAG;

  float  min[MAX_NBANDS], max[MAX_NBANDS];
  int scale;     /* zoom in scale to the coarse resolution temperature (shortwave and thermal band may different) */

} SENSOR;   


/* structure to store samples */
typedef struct {
  float ref[MAX_NBANDS];  /* shortwave bands reflectance */
  float st;               /* skin/surface temperature */
  float w;                /* weight */
} SAMPLES;

void parseParameters(char *ifile, INPUT_PARS *ipar);
void printoutParameters(INPUT_PARS *ipars);
int  getSensorMetaInfo(SENSOR *sensor, SENSOR *st, SENSOR *th, SENSOR *sth, INPUT_PARS *pars);
int  getMetaDataFromGeoTIFF(SENSOR *sensor, int iband);
int  openForWrite(SENSOR *sth);
int  loadThermal(SENSOR *th);
int  applyEC(SENSOR *spec, SENSOR *st, SENSOR *th, SENSOR *sth, INPUT_PARS *ipar);
int  savePureSamples(SENSOR *spec, SENSOR *th, INPUT_PARS *ipar);
int  extractSamples(SENSOR *spec, SENSOR *th, INPUT_PARS *ipar, SAMPLES *sam);
int  resize(SENSOR *st, SENSOR *th);
int  loadSensorRow(SENSOR *sensor, int irow);
int  cleanUpSensor(SENSOR *sensor, SENSOR *st, SENSOR *th, SENSOR *sth);
void dec2bin(unsigned short int num, int *bitpattern, int limit);
void gslRegression(SENSOR *spec, SAMPLES *sam, int NPT, double *coeff);
void localPrediction(SENSOR *spec, SENSOR *sth, INPUT_PARS *pars, double *coeff);
int  openOutput(SENSOR *sth);

void alloc_1dim_contig (void **, int, int);
void alloc_2dim_contig (void ***, int, int, int);
void alloc_3dim_contig (void ****, int, int, int, int);
void alloc_4dim_contig (void *****, int, int, int, int, int);
void free_2dim_contig  (void **);
void free_3dim_contig  (void ***);
void free_4dim_contig  (void ****);
void writeENVIHeader(SENSOR *sensor, char *fname, int nbands, int dtype);

/*#define DEBUG*/
#ifdef DEBUG
#define DEBUG_icol 450
#define DEBUG_irow 450
#endif

/*#define SAVE_COARSE_T*/
