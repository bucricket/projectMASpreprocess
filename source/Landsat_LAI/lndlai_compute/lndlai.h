#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "proj.h"
#include "hdf.h"
#include "mfhdf.h"
#include "HdfEosDef.h"

#define SUCCESS 1
#define FAILURE -1
#define FILLV   -9999   
#define SCALEF  10000  
#define NROWSTEP   10
#define MAX_STRING_LENGTH 2000

#define MATHPI	3.1415926535897932384626433832795
#define D2R	0.017453292519943295769236907684886
#define R2D	57.295779513082320876798154814105

typedef enum 
{
  BLUE=0, GREEN, RED, NIR, SWIR1, SWIR2, CLOUD, NBANDS  
}SRBANDS;


typedef struct 
{
  int32 nrows;           
  int32 ncols;
  float64 ulx;              /* x value of up-left corner */
  float64 uly;              /* y value of up-left corner */
  float64 res;              /* spatial resolution */
  int16 path;
  int16 row;
  /*int16 fillValue;*/
  /* image metadata */
  char Satellite[MAX_STRING_LENGTH];
  char Instrument[MAX_STRING_LENGTH];
  float32 SolarZenith;
  float32 SolarAzimuth;
  float64 WestBoundingCoordinate;
  float64 EastBoundingCoordinate;
  float64 NorthBoundingCoordinate;
  float64 SouthBoundingCoordinate;
  float64 GD_upleft[2];
  float64 GD_lowright[2];
  int32   GD_projcode;
  int32   GD_zonecode;
  int32   GD_spherecode;
  float64 GD_projparm[16];
  int32   GD_origincode;

} META_LANDSAT_SR;

typedef struct 
{
  int32 nrows;           
  int32 ncols;
  float64 ulx;              /* x value of up-left corner */
  float64 uly;              /* y value of up-left corner */
  float64 res;              /* spatial resolution */
  /*int16 fillValue;*/
  /* image metadata */
  float64 GD_upleft[2];
  float64 GD_lowright[2];
  int32   GD_projcode;
  int32   GD_zonecode;
  int32   GD_spherecode;
  float64 GD_projparm[16];
} META_MODIS_LAI;

typedef struct
{
	char longname[MAX_STRING_LENGTH];
	char unit[MAX_STRING_LENGTH];
	int16 range[2];
	int16 fillvalue;
	float64 scale;
	int32 nrows, ncols;
} META_SDS_SR;

typedef struct {
  char fileName[NBANDS][MAX_STRING_LENGTH];     /* SR file name in binary*/
  char gdname[MAX_STRING_LENGTH];
  int32 SD_ID;            /* HDF scientific data grid ID */
  int32 GDfid;            /* HDF-EOS file ID */
  int32 GDid;             /* HDF-EOS grid ID */
  int datatype[NBANDS];
  int fillv[NBANDS];
  float scale[NBANDS];
  FILE *in[NBANDS];
} GRID_LANDSAT_SR;    /* structure to store SR info */

typedef struct {
  char fileName[MAX_STRING_LENGTH];     /* SR file name in binary*/
  char gdname[MAX_STRING_LENGTH];
  int32 SD_ID;            /* HDF scientific data grid ID */
  int32 GDfid;            /* HDF-EOS file ID */
  int32 GDid;             /* HDF-EOS grid ID */
} GRID_LAI;    /* structure to store LAI info */

typedef struct {
  char LandsatFile[NBANDS][MAX_STRING_LENGTH];   /* Landsat SR and QA file (directly from Ledaps)  */
  char ModisFile[MAX_STRING_LENGTH];       /* output composition  file         */
  char SampleFile[MAX_STRING_LENGTH];       /* output composition  file         */
  float pure_threshold;
  char LandsatAncFile[MAX_STRING_LENGTH];
  char LandsatBioPhysFile[MAX_STRING_LENGTH];       /* output composition  file         */
} INPUT_PARS;  /* variables from input file */

typedef struct
{
  double sum[NBANDS-1];
  double sum2[NBANDS-1];
  int land;
  int count;
  int irow;   /* approximate pixel location in Landsat */
  int icol;
} MODIS_CELL;

void InitModisCell(MODIS_CELL* mc);
unsigned char ExtractBit(unsigned char byt, int from, int n);
void dec2bin(unsigned short int num, int *bitpattern, int limit);

void GDSRInit(GRID_LAI* gdsr);
//int getLandsatMetaInfo(META_LANDSAT_SR *metasr, GRID_LANDSAT_SR* gdsr );
int WriteLandsatMetaInfo(META_LANDSAT_SR* metasr, GRID_LAI* gdsr);
int CreateLandsatField(GRID_LAI* gdsr, char* fieldname,int32 datatype);
int ReadnRow(GRID_LAI* gdsr, char* fieldname, int irow,int nrow,int ncol,int n, void* buffer);
int ReadENVInRow(GRID_LANDSAT_SR* gdsr, int iband, int irow, int nrow, int ncol, int n, void* buffer);
int WritenRow(GRID_LAI* gdsr, char* fieldname,int irow,int nrow,int ncol,int n, void* buffer);
int Cleanup(GRID_LAI* gdsr);
int WriteSDSAttr(GRID_LAI* gdsr, META_SDS_SR* metasr, char* fieldname);
int WriteQAAttr(GRID_LAI* gdsr, META_SDS_SR* metasr, char* fieldname);
int getMODISLAIMetaInfo(GRID_LAI* gdsr, META_MODIS_LAI* meta);

void usage(char *command);
void parseParameters(char *ifile, INPUT_PARS *ipar); /* parse input parameters from input parameter file */

void InitLandsatMODISProjInv(META_LANDSAT_SR* metasr, META_MODIS_LAI* metalai);
int Landsat2MODIS(META_LANDSAT_SR* metasr, int row,int col,  META_MODIS_LAI* metalai, int* mrow, int *mcol);

/* error handling and log file */
#define logFile "LogReport.txt"
char msg[MAX_STRING_LENGTH];
FILE *log_fp;
void Errorm(const char *message, const char *module,
           const char *source, long line);
#define ERRORMSG(message, module)   Errorm((message), (module), (__FILE__), (long)(__LINE__))
void WARNING(char *message, char *module);
int openLog(char process[]);
void closeLog();



