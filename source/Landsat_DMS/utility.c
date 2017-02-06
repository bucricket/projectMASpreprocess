/**
 * this program includes sub-routines that support main program 
 * written and revised by Feng Gao on July 2012
 **/

#include "landsat.h"

/* parse and extract input parameters from text file */
void parseParameters(char *ifile, INPUT_PARS *ipar)
{
  FILE *in;
  int i;

  char  buffer[MAX_STRLEN] = "\0";
  char  *label = NULL;
  char  *tokenptr = NULL;
  char  *seperator = "= ,";

  if((in=fopen(ifile,"r"))==NULL) {
    fprintf(stderr, "\nCan't open input %s: parseParameters\n", ifile);
    return;
  }

  ipar->org_th_nrows = -1;
  ipar->org_th_ncols = -1;

  /* process line by line */
  while(fgets(buffer, MAX_STRLEN, in) != NULL) {
 
    /* get string token */
    tokenptr = strtok(buffer, seperator);
    label=tokenptr;

    /* skip comment line */
    if(strcmp(label,"#") == 0) continue;

    while(tokenptr != NULL) {
  
      tokenptr = strtok(NULL, seperator);

      /* get file names */

      if(strcasecmp(label, "NFILES") == 0) 
	sscanf(tokenptr, "%d", &(ipar->nbands)); 
      if(strcasecmp(label, "SW_FILE_NAME") == 0)
	for(i=0; i<ipar->nbands; i++) {  
	  sscanf(tokenptr, "%s ", ipar->inFile[i]);
	  tokenptr = strtok(NULL, seperator);
	} 
      else if(strcasecmp(label, "SW_CLOUD_MASK") == 0)
	sscanf(tokenptr, "%s", ipar->cloudFile);
      else if(strcasecmp(label, "SW_FILE_TYPE") == 0) 
	sscanf(tokenptr, "%s", ipar->fileType); 
      else if(strcasecmp(label, "SW_NCOLS") == 0)
	sscanf(tokenptr, "%d", &(ipar->ncols));
      else if(strcasecmp(label, "SW_NROWS") == 0)
	sscanf(tokenptr, "%d", &(ipar->nrows));
      else if(strcasecmp(label, "SW_PIXEL_SIZE") == 0)
	sscanf(tokenptr, "%f", &(ipar->res));
      else if(strcasecmp(label, "SW_FILL_VALUE") == 0)
	sscanf(tokenptr, "%d", &(ipar->fillv));
      else if(strcasecmp(label, "SW_DATA_RANGE") == 0) {
	sscanf(tokenptr, "%d", &(ipar->range[0]));
	tokenptr = strtok(NULL, seperator);
	sscanf(tokenptr, "%d", &(ipar->range[1]));
      }
      else if(strcasecmp(label, "SW_UPPER_LEFT_CORNER") == 0) {
	sscanf(tokenptr, "%f", &(ipar->ulx));
 	tokenptr = strtok(NULL, seperator);
	sscanf(tokenptr, "%f", &(ipar->uly));
      }
      else if(strcasecmp(label, "SW_PROJECTION_CODE") == 0)
	sscanf(tokenptr, "%ld", &(ipar->insys));
      else if(strcasecmp(label, "SW_PROJECTION_PARAMETERS") == 0)
	for(i=0; i<15; i++) {
	  sscanf(tokenptr, "%lf ", &(ipar->inparm[i]));
	  tokenptr = strtok(NULL, seperator);
	}
      else if(strcasecmp(label, "SW_PROJECTION_ZONE") == 0)
	sscanf(tokenptr, "%ld", &(ipar->inzone));
      else if(strcasecmp(label, "SW_PROJECTION_UNIT") == 0)
	sscanf(tokenptr, "%ld", &(ipar->inunit));
      else if(strcasecmp(label, "SW_PROJECTION_DATUM") == 0)
	sscanf(tokenptr, "%ld", &(ipar->indatum));

      else if(strcasecmp(label, "ORG_TH_FILE_NAME") == 0)
	sscanf(tokenptr, "%s", ipar->org_th_File); 
      else if(strcasecmp(label, "ORG_TH_FILE_TYPE") == 0) 
	sscanf(tokenptr, "%s", ipar->org_th_fileType); 
      else if(strcasecmp(label, "ORG_TH_NROWS") == 0) 
	sscanf(tokenptr, "%d", &(ipar->org_th_nrows)); 
      else if(strcasecmp(label, "ORG_TH_NCOLS") == 0) 
	sscanf(tokenptr, "%d", &(ipar->org_th_ncols)); 
      else if(strcasecmp(label, "ORG_TH_PIXEL_SIZE") == 0)
	sscanf(tokenptr, "%f", &(ipar->org_th_res));
      else if(strcasecmp(label, "ORG_TH_DATA_RANGE") == 0) {
	sscanf(tokenptr, "%f", &(ipar->org_th_range[0]));
	tokenptr = strtok(NULL, seperator);
	sscanf(tokenptr, "%f", &(ipar->org_th_range[1]));
      }

      else if(strcasecmp(label, "RES_TH_PIXEL_SIZE") == 0)
	sscanf(tokenptr, "%f", &(ipar->res_th_res));

      else if(strcasecmp(label, "PURE_CV_TH") == 0)
	sscanf(tokenptr, "%f", &(ipar->CV_TH));
      else if(strcasecmp(label, "SMOOTH_FLAG") == 0)
	sscanf(tokenptr, "%d", &(ipar->SMOOTH_FLAG));      

      else if(strcasecmp(label, "CUBIST_FILE_STEM") == 0) 
	sscanf(tokenptr, "%s", ipar->cubistFile); 
      else if(strcasecmp(label, "OUT_FILE") == 0) 
	sscanf(tokenptr, "%s", ipar->outFile); 

      
      /* in case label (key words) is not the first word in a line */
      label = tokenptr;

    }  /* while token */
  } /* while line */
}



/* print out input parameters */
void printoutParameters(INPUT_PARS *ipars)
{
  int i;
  printf("\nPARAMETER_FILE");
  printf("\n\n# define shortwave band files");
  printf("\nNFILES = %d", ipars->nbands);
  printf("\nSW_FILE_NAME = ");
  for(i=0; i<ipars->nbands; i++)
    printf("%s ", ipars->inFile[i]);
  printf("\nSW_FILE_TYPE = %s", ipars->fileType); 
  printf("\nSW_NCOLS = %d", ipars->ncols);
  printf("\nSW_NROWS = %d", ipars->nrows);
  printf("\nSW_PIXEL_SIZE = %f", ipars->res);
  printf("\nSW_FILL_VALUE = %d", ipars->fillv);
  printf("\nSW_DATA_RANGE = %d %d", ipars->range[0], ipars->range[1]);
  printf("\nSW_UPPER_LEFT_CORNER = %f %f", ipars->ulx, ipars->uly);
  printf("\nSW_PROJECTION_CODE = %ld", ipars->insys);
  printf("\nSW_PROJECTION_PARAMETERS = ");
  for(i=0; i<16; i++)
    printf("%lf ", ipars->inparm[i]);
  printf("\nSW_PROJECTION_ZONE = %ld", ipars->inzone);
  printf("\nSW_PROJECTION_UNIT = %ld", ipars->inunit);
  printf("\nSW_PROJECTION_DATUM = %ld", ipars->indatum);

  printf("\n\n# define original thermal band file");
  printf("\nORG_TH_FILE_NAME = %s", ipars->org_th_File);
  printf("\nORG_TH_FILE_TYPE = %s", ipars->org_th_fileType); 
  printf("\nORG_TH_NROWS = %d", ipars->org_th_nrows);
  printf("\nORG_TH_NCOLS = %d", ipars->org_th_ncols);
  printf("\nORG_TH_PIXEL_SIZE = %f", ipars->org_th_res);
  printf("\nORG_TH_DATA_RANGE = %f %f", ipars->org_th_range[0], ipars->org_th_range[1]);

  printf("\n\n# define resized thermal band file");
  printf("\nRES_TH_PIXEL_SIZE = %f", ipars->res_th_res);

  printf("\n\n# define pure pixel threshold and residuals smooth option");
  printf("\nPURE_CV_TH = %f", ipars->CV_TH);
  printf("\nSMOOTH_FLAG = %d", ipars->SMOOTH_FLAG);

  printf("\n\n# define cubist (decision tree) filestem and final output file");
  printf("\nCUBIST_FILE_STEM = %s", ipars->cubistFile);
  printf("\nOUT_FILE = %s", ipars->outFile);
  printf("\n\nEND\n");
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


/* allocate memory */
void alloc_1dim_contig (void **ptr, int d1, int elsize)
{
   void *p = NULL;

   p = calloc (d1, elsize);
   if (!p) {
     fprintf(stderr, "\nMemory allocation error in alloc_1dim_contig");
     return;
   }
   *ptr = p;
   return;
}


void alloc_2dim_contig (void ***ptr, int d1, int d2, int elsize)
{
   void *p = NULL;
   void **pp = NULL;
   int i = 0;

   /* alloc array for data */
   alloc_1dim_contig ((void **) (&p), d1 * d2, elsize);

   /* alloc array for pointers */
   alloc_1dim_contig ((void **) (&pp), d1, sizeof (void *));

   /* Set the pointers to indicate the appropriate elements of the data array. */
   for (i = 0; i < d1; i++) {
      pp[i] = (char *) p + (i * d2 * elsize);
   }

   *ptr = pp;

   return;
}


void alloc_3dim_contig (void ****ptr, int d1, int d2, int d3, int elsize)
{
   void *p = NULL;
   void **pp = NULL;
   void ***ppp = NULL;
   int i = 0;

   /* allocate the data array */
   alloc_1dim_contig ((void **) (&p), d1 * d2 * d3, elsize);

   /* alloc the double pointers */
   alloc_1dim_contig ((void **) (&pp), d1 * d2, sizeof (void *));

   /* and again for the first dimensional pointers */
   alloc_1dim_contig ((void **) (&ppp), d1, sizeof (void **));

   /* first set the d1 pointers */
   for (i = 0; i < d1; i++) {
      ppp[i] = pp + (i * d2);
   }

   /* next set all of the d2 pointers */
   for (i = 0; i < d1 * d2; i++) {
      pp[i] = (char *) p + (i * d3 * elsize);
   }

   *ptr = ppp;

   return;
}


void alloc_4dim_contig (void *****ptr, int d1, int d2, int d3, int d4, int elsize)
{
   void *p = NULL;
   void **pp = NULL;
   void ***ppp = NULL;
   void ****pppp = NULL;
   int i = 0;

   /* allocate the data array */
   alloc_1dim_contig ((void **) (&p), d1 * d2 * d3 * d4, elsize);

   /* alloc the double pointers */
   alloc_1dim_contig ((void **) (&pp), d1 * d2 * d3, sizeof (void *));

   /* and again for the triple pointers */
   alloc_1dim_contig ((void **) (&ppp), d1 * d2, sizeof (void **));

   alloc_1dim_contig ((void **) (&pppp), d1, sizeof (void ***));

   for (i = 0; i < d1; i++) {
      pppp[i] = ppp + (i * d2);
   }

   for (i = 0; i < d1 * d2; i++) {
      ppp[i] = pp + (i * d3);
   }

   for (i = 0; i < d1 * d2 * d3; i++) {
      pp[i] = (char *) p + (i * d4 * elsize);
   }

   *ptr = pppp;

   return;
}


void free_2dim_contig (void **a)
{
  free (a[0]);
  free (a);
  return;
} 

void free_3dim_contig (void ***a)
{
   free (a[0][0]);
   free (a[0]);
   free (a);
   return;
}

void free_4dim_contig (void ****a)
{
   free (a[0][0][0]);
   free (a[0][0]);
   free (a[0]);
   free (a);
   return;
}


/* write header file in ENVI format */ 
void writeENVIHeader(SENSOR *sensor, char *fname, int nbands, int dtype) 
{
  char str[MAX_STRLEN];
   FILE *out;
  
  sprintf(str, "%s.hdr", fname);
  if((out=fopen(str, "w"))==NULL) exit(1);
  fprintf(out, "ENVI\n");
  strcpy(str, "description = {normalized reflectance} \n");
  fprintf(out, "%s", str);

  fprintf(out, "samples = %d\n", sensor->ncols);
  fprintf(out, "lines = %d\n", sensor->nrows);
  fprintf(out, "bands = %d\n", nbands);
  fprintf(out, "file type = ENVI Standard\n");
  fprintf(out, "data type = %d\n", dtype);
  fprintf(out, "interleave = bsq\n");
  fprintf(out, "byte order = 0\n");
  if(sensor->insys == 1) {
    if(sensor->inzone > 0) 
      fprintf(out, "map info = {UTM, 1.0, 1.0, %f, %f, %f, %f, %d, North, WGS-84, units=Meters}\n", sensor->ulx, sensor->uly, sensor->res, sensor->res, abs(sensor->inzone));
    else
      fprintf(out, "map info = {UTM, 1.0, 1.0, %f, %f, %f, %f, %d, South, WGS-84, units=Meters}\n", sensor->ulx, sensor->uly, sensor->res, sensor->res, abs(sensor->inzone));
  }
  fclose(out);
}

