#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "geotiffio.h"
#include "xtiffio.h"
#include "tiffio.h"

#define FAILURE -1
#define SUCCESS 0

typedef struct {
  char   fileName[500];
  char   outName[500];
  char   instrument[500];
  char   bandname[500];
  TIFF   *fp_tiff;
  int    nrows;
  int    ncols;  
  double res;
  double ulx;
  double uly;
  int fillv;
  int type;
  int utm_zone;
  void  *buf;
} GEOTIFF;

int getGeoTiffInfo(GEOTIFF *sdem);
int writeENVIheader(GEOTIFF *mydata);

int main(int argc, char *argv[])
{
  int i;
  FILE *out;
  GEOTIFF *mydata;

  mydata = malloc(sizeof(GEOTIFF));

  if(argc!=3) {
    printf("Usage: %s <in_GEOTIFF> <out_binary>\n", argv[0]);
    return FAILURE;
  }

  strcpy(mydata->fileName, argv[1]);
  strcpy(mydata->outName, argv[2]);
  if(getGeoTiffInfo(mydata)==FAILURE) {
    printf("Retrieve GeoTiff %s metadata error!\n", mydata->fileName);
    return FAILURE;
  }

  if((out=fopen(mydata->outName, "wb"))==NULL) {
    printf("Open %s error\n", argv[2]);
    return FAILURE;
  }
  
  //printf("nrows=%d ncols=%d\n", mydata->nrows, mydata->ncols);
  //printf("ulx=%7.1f uly=%7.1f res=%f\n", mydata->ulx, mydata->uly, mydata->res);
  
  for(i=0; i<mydata->nrows; i++) {
    if (!TIFFReadScanline(mydata->fp_tiff, mydata->buf, i, 0)) {
      printf("Read line %d Error\n", i);
      return -1;
    } 
    fwrite(mydata->buf, mydata->type, mydata->ncols, out);
  }

  fclose(out);
  XTIFFClose(mydata->fp_tiff);
  writeENVIheader(mydata);

  return SUCCESS;
}



/* get GEOTIFF metadata */
int getGeoTiffInfo(GEOTIFF *mydata) 
{
  int i;
  char proj_info[100], str[100];
  uint16 count, coor_sys;
  double *tiePoint, *pixelScale;
  GTIF *gtif;

  /* get instrument name from filename */
  if(strstr(mydata->fileName, "LC8")) strcpy(mydata->instrument, "LANDSAT_8 OLI_TIRS");
  else if(strstr(mydata->fileName, "LE7")) strcpy(mydata->instrument, "LANDSAT_7 ETM");
  else if(strstr(mydata->fileName, "LT5")) strcpy(mydata->instrument, "LANDSAT_5 TM");
  else strcpy(mydata->instrument, "unkown");

  /* get bandname, data type and fill value based on input filename */
  for(i=1; i<=7; i++) {
    sprintf(str, "band%d", i);
    // if it's a surface reflectance file
    if(strstr(mydata->fileName, str)) {
      sprintf(mydata->bandname, "band %d surface reflectance", i);
      mydata->fillv = -9999;
      mydata->type = sizeof(short int);
    }
  }
  // if it's a cloud mask file
  if(strstr(mydata->fileName, "cfmask")) {
    sprintf(mydata->bandname, "cfmask_band");
    mydata->fillv = 255;
    mydata->type = sizeof(char);
  } 

  if((mydata->fp_tiff = XTIFFOpen(mydata->fileName, "r"))==NULL) {
    fprintf(stderr, "Can't open file %s\n", mydata->fileName);
    return FAILURE;
  }

  /* get  metadata from tiff file */
  if(TIFFGetField(mydata->fp_tiff, TIFFTAG_IMAGEWIDTH, &(mydata->ncols))==0) {
    printf("Retrieve  file %s error\n", mydata->fileName);
    return FAILURE;
  }
  if(TIFFGetField(mydata->fp_tiff, TIFFTAG_IMAGELENGTH, &(mydata->nrows))==0) {
    printf("Retrieve  file %s error\n", mydata->fileName);
    return FAILURE;
  }
  count=6;
  if(TIFFGetField(mydata->fp_tiff, TIFFTAG_GEOTIEPOINTS, &count, &tiePoint)==0) {
    printf("Retrieve  file %s error\n", mydata->fileName);
    return FAILURE;
  }
  count=3;
  if(TIFFGetField(mydata->fp_tiff, TIFFTAG_GEOPIXELSCALE, &count, &pixelScale)==0) {
    printf("Retrieve  file %s error\n", mydata->fileName);
    return FAILURE;
  }

  mydata->res = pixelScale[0];

  /* GeoKey 1025 (GTRasterTypeGeoKey) dictates whether the reference
     coordinate is the UL (*RasterPixelIsArea*, code 1) or center
     (*RasterPixelIsPoint*, code 2) of the UL pixel. If this key is missing,
     the default (as defined by the specification) is to be
     *RasterPixelIsArea*, which is the UL of the UL pixel. */
  gtif = GTIFNew(mydata->fp_tiff);
  if (GTIFKeyGet(gtif, GTRasterTypeGeoKey, &coor_sys, 0, 1) != 1) {
    printf("Coordinate system is not defined in %s\n", mydata->fileName);
    printf("assume used UL of the UL pixel\n");
  }
  if (coor_sys == RasterPixelIsPoint){
    mydata->ulx = tiePoint[3] - 0.5 * mydata->res;
    mydata->uly = tiePoint[4] + 0.5 * mydata->res;
  }
  else {  /* default use RasterPixelIsArea */
    mydata->ulx = tiePoint[3];
    mydata->uly = tiePoint[4];
  }

  if (GTIFKeyGet(gtif, GTCitationGeoKey, proj_info, 0, 33) == 0)
    printf("Projection was not defined\n");
  sscanf(proj_info, "%s %s %d", str, str, &(mydata->utm_zone));  
  if(strstr(proj_info, "Northern Hemisphere") == NULL)
    mydata->utm_zone = -1 * mydata->utm_zone;


  GTIFFree(gtif);
  
  mydata->buf = calloc(mydata->ncols, mydata->type);
  return SUCCESS;
}

int writeENVIheader(GEOTIFF *mydata)
{
  char fname[100];
  FILE *fp;

  sprintf(fname, "%s.hdr", mydata->outName);
  if((fp=fopen(fname, "w"))==NULL) {
    fprintf(stderr, "can't write header for output file %s\n", fname);
    return FAILURE;
  }

  fprintf(fp, "ENVI\n");
  fprintf(fp, "description = {converted from GeoTiff %s}\n", mydata->fileName);
  fprintf(fp, "samples = %d\n", mydata->ncols);
  fprintf(fp, "lines = %d\n", mydata->nrows);
  fprintf(fp, "bands = 1\n");
  fprintf(fp, "header offset = 0\n");
  fprintf(fp, "byte order = 0\n");
  fprintf(fp, "file type = ENVI Standard\n");
  fprintf(fp, "data type = %d\n", mydata->type);
  fprintf(fp, "data ignore value = %d\n", mydata->fillv);
  fprintf(fp, "interleave = BSQ\n");
  fprintf(fp, "sensor type = %s\n", mydata->instrument);
  fprintf(fp, "map info = {UTM, 1.000, 1.000, %f, %f, %f, %f, ", mydata->ulx, mydata->uly, mydata->res, mydata->res); 
  if(mydata->utm_zone > 0)
    fprintf(fp, "%d, North, WGS-84, units=Meters}\n", abs(mydata->utm_zone));
  else
    fprintf(fp, "%d, South, WGS-84, units=Meters}\n", abs(mydata->utm_zone));  
  fprintf(fp, "band names = {%s}\n", mydata->bandname);

  fclose(fp);

  return SUCCESS;
}
