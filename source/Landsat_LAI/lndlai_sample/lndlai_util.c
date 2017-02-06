/**
 * combined subroutines for lndlai programs
 * reads lndsr in ENVI, MODIS LAI in HDF and
 * writes lndlai in HDF format  
 * Feng Gao (2/27/2015) 
 */ 

#include "lndlai.h"
//#define DEBUG

/** 
 * extract input parameters from input text file 
 */
void parseParameters(char *ifile, INPUT_PARS *ipar)
{
  FILE *in;
  int i;
  char  buffer[MAX_STRING_LENGTH] = "\0";
  char  tmpstr[MAX_STRING_LENGTH] = "\0";
  char  *label = NULL;
  char  *tokenptr = NULL;
  char  *seperator = "= ,";

  if((in=fopen(ifile,"r"))==NULL) {
    sprintf(msg, "Can't open input %s", ifile);
    ERRORMSG(msg, "parseParameters");
    return;
  }

  /* process line by line */
  while(fgets(buffer, MAX_STRING_LENGTH, in) != NULL) {
 
    /* get string token */
    tokenptr = strtok(buffer, seperator);
    label=tokenptr;

    /* skip comment line */
    if(strcmp(label,"#") == 0) continue;

    while(tokenptr != NULL) {
  
      tokenptr = strtok(NULL, seperator);

      /* get file names */
      if(strcmp(label, "LANDSAT_BASE_BLUE") == 0) 
	sscanf(tokenptr, "%s", ipar->LandsatFile[BLUE]); 
      else if(strcmp(label, "LANDSAT_BASE_GREEN") == 0) 
	sscanf(tokenptr, "%s", ipar->LandsatFile[GREEN]); 
      else if(strcmp(label, "LANDSAT_BASE_RED") == 0) 
	sscanf(tokenptr, "%s", ipar->LandsatFile[RED]); 
      else if(strcmp(label, "LANDSAT_BASE_NIR") == 0) 
	sscanf(tokenptr, "%s", ipar->LandsatFile[NIR]); 
      else if(strcmp(label, "LANDSAT_BASE_SWIR1") == 0) 
	sscanf(tokenptr, "%s", ipar->LandsatFile[SWIR1]); 
      else if(strcmp(label, "LANDSAT_BASE_SWIR2") == 0) 
	sscanf(tokenptr, "%s", ipar->LandsatFile[SWIR2]); 
      else if(strcmp(label, "LANDSAT_BASE_CLOUD") == 0) 
	sscanf(tokenptr, "%s", ipar->LandsatFile[CLOUD]); 
      else if(strcmp(label, "MODIS_BASE_FILE") == 0)
	sscanf(tokenptr, "%s", ipar->ModisFile);
      else if(strcmp(label, "SAMPLE_FILE_OUT") == 0)
	sscanf(tokenptr, "%s", ipar->SampleFile);
      else if(strcmp(label, "PURE_SAMPLE_TH") == 0)
	sscanf(tokenptr, "%f", &(ipar->pure_threshold));
      else if(strcmp(label, "LANDSAT_ANC_FILE") == 0)
	sscanf(tokenptr, "%s", ipar->LandsatAncFile); 
      else if(strcmp(label, "BIOPHYSICS_PARA_FILE_OUT") == 0)
	sscanf(tokenptr, "%s", ipar->LandsatBioPhysFile);

      /* in case label (key words) is no the first word in a line */
      label = tokenptr;

    }  /* while token */
  } /* while line */

#ifdef DEBUG
  for(i=0; i<NBANDS; i++) printf("%s\n", ipar->LandsatFile[i]);
#endif

  fclose(in);

}

/* allocate memory, get metadata and open specific sds for landsat surface reflectance */
int getLandsatMetaInfo(META_LANDSAT_SR *metasr, GRID_LANDSAT_SR* gdsr ) {

  int  ret;
  char GD_gridlist[100];
  int32 ngrid=0;
  int32 bufsize=100;
  int32 att_id, hdfid, datatype, n_values;
  char s[100], hemi[100];
  float64 t;

  char hdr[MAX_STRING_LENGTH];
  FILE *fp_hdr;

  char  buffer[MAX_STRING_LENGTH] = "\0";
  char  *label = NULL;
  char  *tokenptr = NULL;
  char  *separator = "=\" \t,";
  int i;

  // get HAEDER info from first band
  sprintf(hdr, "%s.hdr", gdsr->fileName[0]);
  if((fp_hdr=fopen(hdr, "r"))==NULL) {
    sprintf(msg, "Failed to open file %s", hdr);
    ERRORMSG(msg, "getLandsatMetaInfo");
    return FAILURE;
  }

  /* process line by line */
  while(fgets(buffer, MAX_STRING_LENGTH, fp_hdr) != NULL) {
   /* get string token */
    tokenptr = strtok(buffer, separator);
    label=tokenptr;
 
    while(tokenptr != NULL) {
 
      tokenptr = strtok(NULL, separator);
      if(strcmp(label, "lines") == 0)
	sscanf(tokenptr, "%d", &(metasr->nrows));
      if(strcmp(label, "samples") == 0)
	sscanf(tokenptr, "%d", &(metasr->ncols));
      // reads projection info
      if(strcmp(label, "map") == 0) {
	for(i=0; i<4; i++) 
	  tokenptr = strtok(NULL, separator);
	sscanf(tokenptr, "%lf", &(metasr->ulx));
	tokenptr = strtok(NULL, separator);
	sscanf(tokenptr, "%lf", &(metasr->uly));
	tokenptr = strtok(NULL, separator);
	sscanf(tokenptr, "%lf", &(metasr->res));
	tokenptr = strtok(NULL, separator);
	tokenptr = strtok(NULL, separator);
	sscanf(tokenptr, "%d", &(metasr->GD_zonecode));
	tokenptr = strtok(NULL, separator);
	sscanf(tokenptr, "%s", hemi);
      }
      /* in case label (key words) is no the first word in a line */
      label = tokenptr;
    }
  }
  fclose(fp_hdr);

  metasr->GD_upleft[0] = metasr->ulx;
  metasr->GD_upleft[1] = metasr->uly;
  metasr->GD_lowright[0] =  metasr->ulx + metasr->ncols * metasr->res;
  metasr->GD_lowright[1] =  metasr->uly - metasr->nrows * metasr->res;

  metasr->GD_projcode = 1;
  // use negative zone number for southern hemisphere
  if(strcasecmp(hemi, "South")==0) metasr->GD_zonecode = -metasr->GD_zonecode;
  // use default Landsat setting: WGS-84
  metasr->GD_spherecode = 12; 
  metasr->GD_origincode = 0;

  for(i=0; i<16; i++) metasr->GD_projparm[i] = 0.0;

  for(i=0; i<NBANDS-1; i++) {
    gdsr->fillv[i] = -9999;
    gdsr->scale[i] = 0.0001;
    gdsr->datatype[i] = 2;
  }
  gdsr->fillv[CLOUD] = 255;
  gdsr->scale[CLOUD] = 1;
  gdsr->datatype[CLOUD] = 1;

  // open all bands and cfmask
  for(i=0; i<NBANDS; i++) {
    if((gdsr->in[i]=fopen(gdsr->fileName[i], "rb"))==NULL) {
      sprintf(msg, "Failed to open file %s", gdsr->fileName[i]);
      ERRORMSG(msg, "getLandsatMetaInfo");
      return FAILURE;
    }
  }
  
#ifdef DEBUG
  printf("\nFile: %s\n",gdsr->fileName[0]);
  printf("nrows:%d, ncols:%d\n",metasr->nrows, metasr->ncols);
  printf("ulx:%8.1f, uly:%8.1f\n", metasr->ulx, metasr->uly);
  printf("resolution:%5.1f\n", metasr->res);
  printf("projcode:%d zonecode:%d  spherecode:%d origincode:%d\n", 
	 metasr->GD_projcode, metasr->GD_zonecode, metasr->GD_spherecode, metasr->GD_origincode);
  for(i=0; i<16; i++)
    printf("%d: %f\n", i, metasr->GD_projparm[i]);   
#endif

  return SUCCESS;
}


int getMODISLAIMetaInfo(GRID_LAI* gdsr, META_MODIS_LAI* meta)
{
  int  ret;
  char GD_gridlist[100];
  int32 ngrid=0;
  int32 bufsize=100;


  /* open a hdf file */
  gdsr->GDfid = GDopen(gdsr->fileName, DFACC_READ);
  if(gdsr->GDfid==FAILURE){
    sprintf(msg, "Not successful in retrieving grid file %s ID/open", gdsr->fileName);
    ERRORMSG(msg, "getMODISLAIMetaInfo");
    return FAILURE;
  }

  /* find out about grid type */
  ngrid=GDinqgrid(gdsr->fileName, GD_gridlist, &bufsize);
  if(ngrid==FAILURE){
    sprintf(msg, "Not successful in retrieving %s grid name list", gdsr->fileName);
    ERRORMSG(msg, "getMODISLAIMetaInfo");
    return FAILURE;
  }
  strcpy(gdsr->gdname,GD_gridlist);
  /* attach grid */
  gdsr->GDid = GDattach(gdsr->GDfid, GD_gridlist);
  if(gdsr->GDid==FAILURE){
    sprintf(msg, "Not successful in attaching grid for %s.", gdsr->fileName);
    ERRORMSG(msg, "getMODISLAIMetaInfo");
    return FAILURE;
  }

  /* get grid info */
  ret = GDgridinfo(gdsr->GDid, &(meta ->ncols), &(meta->nrows), meta->GD_upleft, meta->GD_lowright);
  if(ret==FAILURE){
    sprintf(msg, "Failed to read grid info. from %s", gdsr->fileName);
    ERRORMSG(msg, "getMODISLAIMetaInfo");
    return FAILURE;
  }
  meta->ulx = meta->GD_upleft[0];
  meta->uly = meta->GD_upleft[1];
  meta->res = fabs(meta->GD_upleft[0] - meta->GD_lowright[0])/meta->ncols;
  
  /* get projection parameters */
  ret = GDprojinfo(gdsr->GDid, &meta->GD_projcode, &meta->GD_zonecode, &meta->GD_spherecode, meta->GD_projparm);
  if(ret==FAILURE){
    sprintf(msg, "Not successful in reading grid projection info from %s", gdsr->fileName);
    ERRORMSG(msg, "getMODISLAIMetaInfo");
    return FAILURE;
  }

  return SUCCESS;

}


int WriteLandsatMetaInfo(META_LANDSAT_SR* metasr, GRID_LAI* gdsr)
{
  int ret;


  /* create output use HDF-EOS functions */
  gdsr->GDfid = GDopen(gdsr->fileName, DFACC_CREATE);
  if(gdsr->GDfid == FAILURE){
    sprintf(msg, "Not successful in creating grid file %s", gdsr->fileName);
    ERRORMSG(msg, "writeMetaInfo");
    return FAILURE;
  }
 
  /* create a new grid in output */ 
  gdsr->GDid = GDcreate(gdsr->GDfid, gdsr->gdname, metasr->ncols, metasr->nrows, metasr->GD_upleft, metasr->GD_lowright);
  if(gdsr->GDid == FAILURE) {
    sprintf(msg, "Not successful in creating grid ID/ create for %s", gdsr->fileName);
    ERRORMSG(msg, "writeMetaInfo");
    return FAILURE;
  }

  /* define grid projection */
  ret = GDdefproj(gdsr->GDid, metasr->GD_projcode, metasr->GD_zonecode, metasr->GD_spherecode, metasr->GD_projparm);
  if(ret==FAILURE){
    sprintf(msg, "Not successful in defining grid projection");
    ERRORMSG(msg, "writeMetaInfo");
    return FAILURE;
  }

  /* define grid origin */
  ret = GDdeforigin(gdsr->GDid, metasr->GD_origincode);
  if(ret==FAILURE){
    sprintf (msg, "Not successful in defining grid origin");
    ERRORMSG(msg, "writeMetaInfo");
    return FAILURE;
  }

  ret = GDdetach(gdsr->GDid);
  if(ret==FAILURE){
    sprintf (msg, "Not successful in detach grid");
    ERRORMSG(msg, "writeMetaInfo");
    return FAILURE;
  }
  gdsr->GDid = -1;

  if(GDclose(gdsr->GDfid) == FAILURE)
    {
      sprintf (msg, "Not successful in closing grid");
      ERRORMSG(msg, "writeMetaInfo");
    }

  /*get SD_ID for this file*/
  if( (gdsr->SD_ID = SDstart(gdsr->fileName, DFACC_RDWR)) == FAILURE)
    {
      sprintf (msg, "Not successful in openning SD");
      ERRORMSG(msg, "writeMetaInfo");
    }
 
  /*if ((SDsetattr(gdsr->SD_ID, "PixelSize", DFNT_FLOAT64, 1, &(metasr->res))) == FAILURE)
    {
    ERRORMSG("Can't write PixelSize to output HDF file", "writeMetaInfo");
    return FAILURE;
    }*/
  
  if(SDend(gdsr->SD_ID)== FAILURE)
    {
      sprintf (msg, "Not successful in closing SD");
      ERRORMSG(msg, "writeMetaInfo");
    }
  gdsr->SD_ID = -1;

  if( (gdsr->GDfid = GDopen(gdsr->fileName, DFACC_RDWR)) == FAILURE)
    {
      sprintf (msg, "Not successful in openning grid");
      ERRORMSG(msg, "writeMetaInfo");

    }

  return SUCCESS;
}


int CreateLandsatField(GRID_LAI* gdsr, char* fieldname,int32 datatype)
{
  int32 ret;
	

  if(gdsr->GDid == -1)
    {
      gdsr->GDid = GDattach(gdsr->GDfid, gdsr->gdname);
      if(gdsr->GDid == FAILURE)
	{
	  sprintf (msg, "Not successful in attaching GRID %s in %s", gdsr->gdname, gdsr->fileName);
	  ERRORMSG(msg, "CreateLandsatField");
	  return FAILURE;
	}
    }

  ret = GDdeffield(gdsr->GDid, fieldname, "YDim,XDim", datatype, HDFE_NOMERGE);
  if(ret==FAILURE)
    {
      sprintf (msg, "Not successful in defining field %s in grid %s", fieldname,gdsr->gdname );
      ERRORMSG(msg, "CreateLandsatField");
      return FAILURE;
    }
 
  return SUCCESS;
}

int WriteSDSAttr(GRID_LAI* gdsr, META_SDS_SR* metasr, char* fieldname)
{
	
  int32 sr_id, index, ret;
  /*get SD_ID for this file*/
  if( (gdsr->SD_ID = SDstart(gdsr->fileName, DFACC_RDWR)) == FAILURE)
    {
      sprintf (msg, "Not successful in openning SD");
      ERRORMSG(msg, "writeSDSAttr");
    } 

  /*write SDS level attribute */
  if ((index = SDnametoindex(gdsr->SD_ID,fieldname))<0) {
    sprintf(msg, "Not successful in convert SR sdsName %s to index", fieldname);
    ERRORMSG(msg, "WriteSDSAttr");
    return FAILURE;
  }
  sr_id = SDselect(gdsr->SD_ID, index);
 
  /*	ret = SDsetattr(sr_id, "long_name", DFNT_CHAR8, (int32)strlen(metasr->longname), metasr->longname);
	if (ret == FAILURE) {
	sprintf(msg, "Can't write SR long_name for SDS %s\n", fieldname);
	ERRORMSG(msg, "WriteSDSAttr");
	return FAILURE;
	}

  */  
  ret = SDsetattr(sr_id, "units", DFNT_CHAR8, (int32)strlen(metasr->unit), metasr->unit);
  if (ret == FAILURE) {
    sprintf(msg, "Can't write SR units for SDS %s", fieldname);
    ERRORMSG(msg, "WriteSDSAttr");
    return FAILURE;
  } 

  ret = SDsetattr(sr_id, "valid_range", DFNT_INT16, 2, metasr->range);
  if (ret == FAILURE) {
    sprintf (msg, "Can't write SR valid_range for SDS %s", fieldname);
    ERRORMSG(msg, "WriteSDSAttr");
    return FAILURE;
  } 

  ret = SDsetattr(sr_id, "_FillValue", DFNT_INT16, 1, &(metasr->fillvalue));
  if (ret == FAILURE) {
    sprintf (msg, "Can't write SR _FillValue for SDS %s",fieldname);
    ERRORMSG(msg, "WriteSDSAttr");
    return FAILURE;
  } 

  ret = SDsetattr(sr_id, "scale_factor", DFNT_FLOAT64, 1, &(metasr->scale));
  if (ret == FAILURE) {
    sprintf (msg, "Can't write SR scale_factor for SDS %s", fieldname);
    ERRORMSG(msg, "WriteSDSAttr");
    return FAILURE;
  } 
    
  if((SDendaccess(sr_id)) == FAILURE) {
    sprintf(msg, "SDendaccess for %s reflectance error!", gdsr->fileName);
    ERRORMSG(msg, "WriteSDSAttr");
    return FAILURE;
  }
	
  if( SDend(gdsr->SD_ID) == FAILURE)
    {
      sprintf(msg, "SDend for %s reflectance error!", gdsr->fileName);
      ERRORMSG(msg, "WriteSDSAttr");
      return FAILURE;
    }

  return SUCCESS;

}

int WriteQAAttr(GRID_LAI* gdsr, META_SDS_SR* metasr, char* fieldname)
{
  int32 sr_id, index, ret;
  uint8 range[2];
  char *bitmapindex = "0=clear_land; 1=water; 2=cloud_shadow; 3=snow; 4=cloud; 5=adjacent_cloud; 255=fill";

  /*get SD_ID for this file*/
  if( (gdsr->SD_ID = SDstart(gdsr->fileName, DFACC_RDWR)) == FAILURE)
    {
      sprintf (msg, "Not successful in openning SD");
      ERRORMSG(msg, "WriteQAAttr");
    } 

  /*write SDS level attribute */
  if ((index = SDnametoindex(gdsr->SD_ID,fieldname))<0) {
    sprintf(msg, "Not successful in convert SR sdsName %s to index", fieldname);
    ERRORMSG(msg, "WriteSDSAttr");
    return FAILURE;
  }
  sr_id = SDselect(gdsr->SD_ID, index);
 
  /*ret = SDsetattr(sr_id, "valid_range", DFNT_UINT8, 2, metasr->range);
  if (ret == FAILURE) {
    sprintf (msg, "Can't write SR valid_range for SDS %s", fieldname);
    ERRORMSG(msg, "WriteQAAttr");
    return FAILURE;
  } 

  ret = SDsetattr(sr_id, "_FillValue", DFNT_UINT8, 1, &(metasr->fillvalue));
  if (ret == FAILURE) {
    sprintf (msg, "Can't write SR _FillValue for SDS %s",fieldname);
    ERRORMSG(msg, "WriteQAAttr");
    return FAILURE;
    } */

  ret = SDsetattr(sr_id, "Legend", DFNT_CHAR8, strlen(bitmapindex), bitmapindex);
  if (ret == FAILURE) {
    sprintf (msg, "Can't write QAbitmap index for SDS %s", fieldname);
    ERRORMSG(msg, "WriteQAAttr");
    return FAILURE;
  } 

  if((SDendaccess(sr_id)) == FAILURE) {
    sprintf(msg, "SDendaccess for %s reflectance error!", gdsr->fileName);
    ERRORMSG(msg, "WriteQAAttr");
    return FAILURE;
  }
	
  if( SDend(gdsr->SD_ID) == FAILURE)
    {
      sprintf(msg, "SDend for %s reflectance error!", gdsr->fileName);
      ERRORMSG(msg, "WriteQAAttr");
      return FAILURE;
    }

  return SUCCESS;

}


/* read Landsat SR from ENVI file starting from irow to irow+n for iband */
int ReadENVInRow(GRID_LANDSAT_SR* gdsr, int iband, int irow, int nrow, int ncol, int n, void* buffer)
{
  int32 start[2];
  int32 stride[2];
  int32 length[2];
  long offset;
  
  if(irow >= nrow) return 0;
  start[0] = irow;
  start[1] = 0;
  stride[0] = 1;
  stride[1] = 1;
  if(irow +n <nrow)  	length[0] = n;
  else length[0] =nrow-irow;
  length[1] = ncol;

  offset = (long)(start[0] * ncol * gdsr->datatype[iband]);
  fseek(gdsr->in[iband], offset, 0);
  fread(buffer, gdsr->datatype[iband], length[0]*length[1], gdsr->in[iband]);

  return length[0];
}


/* read MODIS data from irow to irow+n */
int ReadnRow(GRID_LAI* gdsr, char* fieldname, int irow,int nrow,int ncol,int n, void* buffer)
{
  int32 start[2];
  int32 stride[2];
  int32 length[2];
  
  if(irow >= nrow) return 0;
  start[0] = irow;
  start[1] = 0;
  stride[0] = 1;
  stride[1] = 1;
  if(irow +n <nrow)  	length[0] = n;
  else length[0] =nrow-irow;
  length[1] = ncol;

  if((GDreadfield(gdsr->GDid, fieldname, start, stride, length, buffer)) == FAILURE) 
    {
      sprintf(msg, "reading %s reflectance error irow: %d", gdsr->fileName, irow);
      ERRORMSG(msg, "ReadnRow");
      return FAILURE;
    }
  return length[0];
}


/* write n rows data to HDF file (use HDF-EOS function)*/
int WritenRow(GRID_LAI* gdsr, char* fieldname,int irow,int nrow,int ncol,int n, void* buffer)
{
  int32 start[2];
  int32 stride[2];
  int32 length[2];
  
  start[0] = irow;
  start[1] = 0;
  stride[0] = 1;
  stride[1] = 1;
  if(irow +n <nrow)  	length[0] = n;
  else length[0] =nrow-irow;
  if(length[0]<0) return 0;
  length[1] = ncol;

  if((GDwritefield(gdsr->GDid, fieldname, start, stride, length, buffer)) == FAILURE) 
    {
      sprintf(msg, "writing %s reflectance error irow: %d", gdsr->fileName, irow);
      ERRORMSG(msg, "WritenRow");
      return FAILURE;
    }
  return SUCCESS;
}


void GDSRInit(GRID_LAI* gdsr)
{
  gdsr->GDfid = -1;
  gdsr->GDid = -1;
  gdsr->SD_ID = -1;
}


void InitLandsatMODISProjInv(META_LANDSAT_SR* metasr, META_MODIS_LAI* metalai)
{
  long iflg;
  long (*inv_trans[MAXPROJ+1])();

	
  /*UTM projection inverse initialize */
  inv_init((long)metasr->GD_projcode, (long)metasr->GD_zonecode, metasr->GD_projparm, (long)metasr->GD_spherecode, 
	   NULL, NULL, &iflg, inv_trans);

  /*SIN projection forward initialize */
  for_init((long)metalai->GD_projcode, (long)metalai->GD_zonecode, metalai->GD_projparm, (long)metalai->GD_spherecode, 
	   NULL, NULL, &iflg, inv_trans);

}

void InitModisCell(MODIS_CELL* mc)
{
  int i;
  for(i=0; i<NBANDS-1; i++)
    {
      mc->sum[i] = mc->sum2[i] =0;	
    }
  mc->count =0;
  mc->land  =0;
}


int Cleanup(GRID_LAI* gdsr)
{
  int flag =SUCCESS;

  if(gdsr->GDid != -1)
    {
      if(GDdetach(gdsr->GDid)== FAILURE) 
	{
	  sprintf(msg, "GDdetach for %s error!", gdsr->fileName);
	  ERRORMSG(msg, "Cleanup_Landsat_SR");
	  flag =FAILURE; 
	}
      gdsr->GDid = -1;
    }
  if(gdsr->GDfid != -1)
    {
      if(GDclose(gdsr->GDfid)== FAILURE) 
	{
	  sprintf(msg, "GDclose for %s error!", gdsr->fileName);
	  ERRORMSG(msg, "Cleanup_Landsat_SR");
	  flag =FAILURE; 
	}
      gdsr->GDfid = -1;
    }

  return flag;

}


int Landsat2MODIS(META_LANDSAT_SR* metasr, int row,int col,  META_MODIS_LAI* metalai, int* mrow, int *mcol)
{
  double LndEast, LndNorth, ModisEast, ModisNorth, lat, lon;
	
  LndEast = col*metasr->res + metasr->ulx;
  LndNorth = metasr->uly - row*metasr->res;

  utminv(LndEast, LndNorth, &lon, &lat);
  sinfor(lon, lat, &ModisEast, &ModisNorth);
  *mrow = (int)((metalai->uly - ModisNorth)/metalai->res);
  *mcol = (int)((ModisEast - metalai->ulx)/metalai->res);
  return SUCCESS;
}


unsigned char ExtractBit(unsigned char byt, int from, int n)
{
  static unsigned char mask[8] = {0x01, 0x03, 0x07, 0x0f, 0x1f, 0x3f, 0x7f, 0xff};	
	
  if(n==8) return byt;
  if(n==0 || from >7) return 0;
  if(from>0)
    {
      byt = byt>>from;
    }
  return (byt&mask[n-1]);
}


/* convert decimal number to binary number */
void dec2bin(unsigned short int num, int *bitpattern, int limit)
{
 
  register int i=0;       
 
  for(i=0; i<limit; i++)
    bitpattern[i]=0;
 
  for(i=0; num>0; i++)
    {
      bitpattern[i] = (int)num & 1;
      num >>= 1;
    }
}


void Errorm(const char *message, const char *module,
           const char *source, long line)
{
  fprintf(stderr, " error [%s, %s:%ld] : %s\n", module, source, line, message);
  fprintf(log_fp, " error [%s, %s:%ld] : %s\n", module, source, line, message);
}

void WARNING(char *message, char *module) {
  fprintf(stderr, "WARNING [%s, %s:%ld] : %s\n", module, __FILE__, __LINE__, message);
  fprintf(log_fp, "WARNING [%s, %s:%ld] : %s\n", module, __FILE__, __LINE__, message);
} 

/* open log file for append and add time stamp on it */
int openLog(char process[])
{
  struct tm *currtime;
  time_t t; 
  char str[MAX_STRING_LENGTH];

  if((log_fp=fopen(logFile, "a"))==NULL) {
    ERRORMSG("open log file as input", "openLog");
    return FAILURE;
  }
  t = time(NULL);
  currtime = (struct tm *)gmtime(&t);
  strftime(str, 100, "%FT%H:%M:%SZ", currtime);
  fprintf(log_fp, "\n\n##################################\n");
  fprintf(log_fp, "Start Process <%s> on %s\n", process, str);
  return SUCCESS;
}

void closeLog()
{
  struct tm *currtime;
  time_t t; 
  char str[MAX_STRING_LENGTH];

  t = time(NULL);
  currtime = (struct tm *)gmtime(&t);
  strftime(str, 100, "%FT%H:%M:%SZ", currtime);
  fprintf(log_fp, "End Process on %s\n", str);
  fclose(log_fp);
}

