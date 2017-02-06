/**
   This program combines prediction results from global and local models 
   written by Feng Gao (feng.gao@ars.usda.gov) on June 2012
   v1.0_HDF: uses LEDAPS SR in HDF format (in-house LEDAPS product)
   v1.1_ENVI: accepts GeoTIFF, ENVI and HDF format (USGS Landsat SR product)
   v1.2_ENVI: uses Landsat cloud mask (USGS Landsat fmask product)
**/

#include "landsat.h"

int main(int argc, char *argv[])
{
  char  name[MAX_STRLEN], buffer[MAX_STRLEN], choice;
  int   i, j, m, n, scale;
  float p1, p2, d1, d2, w1, w2, tw;
  double gnum, lnum;
  FILE  *fp1, *fp2, *fps, *tp;
 
  INPUT_PARS *pars;
  SENSOR     *spec, *fine_th, *coarse_th, *global_th, *local_th, *combined_th;

  /* allocate memory for variables */
  pars = malloc(sizeof(INPUT_PARS)); /* input parameters */
  spec = malloc(sizeof(SENSOR));     /* input fine resolution spectral data */
  
  /* NOTE: define original (fine) thermal data and intermediate (coarse) thermal data for Landsat
     since USGS EROS data center only provide resampled 30 m resolution data and need to 
     aggregtae to coarse resolution before using DMS approach */
  fine_th   = malloc(sizeof(SENSOR));     /* input original fine resolution thermal data */
  coarse_th = malloc(sizeof(SENSOR));     /* intermediate coarse resolution thermal data */

  global_th   = malloc(sizeof(SENSOR));    /* input global prediction */
  local_th   = malloc(sizeof(SENSOR));     /* input local prediction  */
  combined_th  = malloc(sizeof(SENSOR));   /* output sharpened thermal data (combined) */  

  printf("- parsing parameters and extracting metadata ");
  if(argc==2) {
    /* get input parameters from input file */
    parseParameters(argv[1], pars);
  }
  else {
    fprintf(stderr, "\nUsage: %s <input_pars_file>\n", argv[0]);
    exit(1);
  }

#ifdef DEBUG
  fps=fopen("weight.bin", "wb");
  tp=fopen("temp.bin", "wb"); 
#endif

  /* get metadata from inputs */
  if(getSensorMetaInfo(spec, fine_th, coarse_th, combined_th, pars)==FAILURE) {
    fprintf(stderr, "\nmain: Retrieve Sensor metadata error\n");
    exit(1);
  }

  openForWrite(combined_th);

  printf("\n- loading original fine resolution temperature file ...");
  loadThermal(fine_th);

  printf("\n- reducing to coarse resolution ...");
  printf("\n\tfor original TIR");
  resize(fine_th, coarse_th);

#ifdef DEBUG
  for(i=0; i<coarse_th->nrows; i++)
    for(j=0; j<coarse_th->ncols; j++)
      fwrite(&(coarse_th->fdata[i][j]), sizeof(float), 1, tp);
#endif

  strncpy(buffer, combined_th->fileName[0], strlen(combined_th->fileName[0])-4);
  buffer[strlen(combined_th->fileName[0])-4] = '\0';
  
  printf("\n\tfor global model prediction");
  sprintf(name, "%s.global", buffer);
  printf("\n\tfilename: %s", name);
  if((fp1=fopen(name, "rb"))==NULL) {
    printf("Open file %s error!", name);
    return -1;
  }

  global_th->scale = coarse_th->scale;
  global_th->res = coarse_th->res;
  global_th->nrows = coarse_th->nrows;
  global_th->ncols = coarse_th->ncols;
  global_th->range[0] = coarse_th->range[0];
  global_th->range[1] = coarse_th->range[1];
  alloc_2dim_contig((void ***) (&global_th->fdata), global_th->nrows, global_th->ncols, sizeof(float));
  for(i=0; i<combined_th->nrows; i++)
    for(j=0; j<combined_th->ncols; j++) {
      fread(&p1, sizeof(float), 1, fp1);
      fine_th->fdata[i][j] = p1;
    }
  resize(fine_th, global_th);

#ifdef DEBUG
  for(i=0; i<global_th->nrows; i++)
    for(j=0; j<global_th->ncols; j++)
      fwrite(&(global_th->fdata[i][j]), sizeof(float), 1, tp);
#endif

  printf("\n\tfor local model prediction");
  sprintf(name, "%s.local", buffer);
  printf("\n\tfilename: %s", name);
  if((fp2=fopen(name, "rb"))==NULL) {
    printf("Open file %s error!", name);
    return -1;
  }
  alloc_2dim_contig((void ***) (&local_th->fdata), coarse_th->nrows, coarse_th->ncols, sizeof(float));
  local_th->scale = coarse_th->scale;
  local_th->res = coarse_th->res;
  local_th->nrows = coarse_th->nrows;
  local_th->ncols = coarse_th->ncols;
  local_th->range[0] = coarse_th->range[0];
  local_th->range[1] = coarse_th->range[1];
  for(i=0; i<combined_th->nrows; i++)
    for(j=0; j<combined_th->ncols; j++) {
      fread(&p2, sizeof(float), 1, fp2);
      fine_th->fdata[i][j] = p2;
    }
  resize(fine_th, local_th);

#ifdef DEBUG
  for(i=0; i<local_th->nrows; i++)
    for(j=0; j<local_th->ncols; j++)
      fwrite(&(local_th->fdata[i][j]), sizeof(float), 1, tp);
#endif

  rewind(fp1);
  rewind(fp2);
  
  printf("\n- combining predictions from global and local models ....");
  // fine_th->CMASK_FLAG = 1;  /* 0=include cloud 1=exclude cloud */
  loadThermal(fine_th);  
  gnum = 0;
  lnum = 0;
  scale = fine_th->scale;

  for(i=0; i<combined_th->nrows; i++)
    for(j=0; j<combined_th->ncols; j++) {

      fread(&p1, sizeof(float), 1, fp1); /* global model result */
      fread(&p2, sizeof(float), 1, fp2); /* local model result */
      m = i/scale;
      n = j/scale;

      /*if(coarse_th->fdata[m][n] < coarse_th->rang[0] || coarse_th->fdata[m][n] > coarse_th->rang[1]) {
	combined_th->fdata[i][j] = coarse_th->fillValue;
	continue;
	}*/

      d1 = global_th->fdata[m][n] - coarse_th->fdata[m][n]; /* error from global model */
      d2 = local_th->fdata[m][n] - coarse_th->fdata[m][n];  /* error from local model */

      /*if(d1<d2-1.0) {
	combined_th->fdata[i][j] = p1;
	gnum++;
	choice = 1;
      }
      else {
	combined_th->fdata[i][j] = p2;
	lnum++;
	choice = 2;
	}*/

      /* weight global model */
      if(fabs(d1)<0.1)
	w1 = 100.0;
      else 
	w1 = 1.0/(d1*d1);
 
      /* weight local model */
      if(fabs(d2)<0.1) 
	w2 = 100.0;
      else
	w2 = 1.0/(d2*d2);
    
      tw = w1+w2;
      w1 = w1/tw;
      w2 = w2/tw;
      combined_th->fdata[i][j] = p1*w1+p2*w2;

      /* printf("(%d, %d) o=%5.1f g=%5.1f l=%5.1f d1=%5.2f d2=%5.2f w1=%5.3f w2=%5.3f c=%5.1f\n", n+1, m+1, coarse_th->fdata[m][n], global_th->fdata[m][n], local_th->fdata[m][n], d1, d2, w1, w2, combined_th->fdata[i][j]);  */

      gnum += w1;
      lnum += w2;

      if(combined_th->fdata[i][j]<combined_th->range[0]||combined_th->fdata[i][j]>combined_th->range[1]) {
	combined_th->fdata[i][j] = FILLV;
	choice = 0;
      }

#ifdef DEBUG
      fwrite(&w1, sizeof(float), 1, fps);
      fwrite(&w2, sizeof(float), 1, fps);
#endif
    }

  printf("\n\tglobal model used %4.1f%%; local model used %4.1f%%", gnum/(gnum+lnum)*100.0, lnum/(gnum+lnum)*100.0);
  printf("\n- applying energy conservation to the combined model ...");
  applyEC(spec, fine_th, coarse_th, combined_th, pars);


#ifdef DEBUG
  /* save combined results for debuging */
  resize(combined_th, coarse_th);
  for(i=0; i<coarse_th->nrows; i++)
    for(j=0; j<coarse_th->ncols; j++)
      fwrite(&(coarse_th->fdata[i][j]), sizeof(float), 1, tp);
#endif      
  
  printf("\n- cleaning up memory and closing files ...");
  cleanUpSensor(spec, fine_th, coarse_th, combined_th);

  free(coarse_th);
  free(fine_th);
  free(combined_th);
  free(global_th);
  free(local_th);
  free(spec);
  free(pars);

  fclose(fp1);
  fclose(fp2); 

#ifdef DEBUG
  fclose(fps); 
  fclose(tp); 
#endif

  printf("\n- program ends successfully\n");

  return 0;
}
