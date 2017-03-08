/**
 * ! This program extracts "pure" and homogeneous temperature (T) and reflectance (R) samples for training T-R 
 *   relationships using regression tree program (cubist)
 *   
 * !Usage
 *  for the entire image use:  cubist_sample.exe <in_parameter_file>
 *  for the partial image use: cubist_sample.exe <in_parameter_file> <s_row> <s_col> <e_row> <e_col>
 *
 * !Input
 *  spectral bands and thermal band data
 *
 * !Output
 *  T-R samples and name that are acceptable by cubist program
 * 
 * !Developer
 *  Feng Gao (Feng.Gao@ars.usda.gov)
 *
 * !Revision
 *
 * Original version - 06/2012 by Feng Gao (USDA-ARS Hydrology and Remote Sensing Lab)
 * v1.0_HDF: uses LEDAPS SR in HDF format (in-house LEDAPS product)
 * v1.1_ENVI: accepts GeoTIFF, ENVI and HDF format (USGS Landsat SR product)
 * v1.2_ENVI: uses Landsat cloud mask (USGS Landsat fmask product)
 */


#include "landsat.h"
#define SAVE_COARSE_T

int main(int argc, char *argv[])
{
  int  npures;
  INPUT_PARS *pars;
  SENSOR     *spec, *th, *sth, *st;
 
  /* allocate memory for variables */
  pars = malloc(sizeof(INPUT_PARS));
  spec = malloc(sizeof(SENSOR));  /* input spectral */
  th   = malloc(sizeof(SENSOR));  /* original input thermal band */
  st   = malloc(sizeof(SENSOR));  /* resized thermal band */
  sth  = malloc(sizeof(SENSOR));  /* output sharpened thermal */
  
  switch(argc) {
  case 2: 
    parseParameters(argv[1], pars); 
    break;
  case 6: 
    parseParameters(argv[1], pars); 
    pars->s_row = atoi(argv[2]);
    pars->s_col = atoi(argv[3]);
    pars->e_row = atoi(argv[4]);
    pars->e_col = atoi(argv[5]);
    break;
  default: 
    fprintf(stderr, "\nUsage: %s <in_parameter_file>", argv[0]);
    fprintf(stderr, "\nor:    %s <in_parameter_file> <s_row> <s_col> <e_row> <e_col> (in coarse T res)\n", argv[0]);
    return FAILURE;
  }
 

#ifdef DEBUG
    printoutParameters(pars);
#endif
  
  printf("- parsing parameters and extracting metadata ");
  /* get metadata from inputs */
  if(getSensorMetaInfo(spec, st, th, sth, pars)==FAILURE) {
    fprintf(stderr, "\nmain: Retrieve Sensor metadata error\n");
    return FAILURE;
  }

  /* use the entire image if not defined */
  if(argc == 2) {
    pars->s_row = 0;
    pars->s_col = 0;
    pars->e_row = th->nrows;
    pars->e_col = th->ncols;
  }

  printf("\n- loading thermal data from original resolution file ...");
  /* read all thermal data to memory */
  loadThermal(st);

  printf("\n- aggregating original temperature to corase (actual) resolution ...");
  resize(st, th);
  
  printf("\n- extracting and saving pure pixels in coarse resolution ... ");
  npures=savePureSamples(spec, th, pars);
  printf("\n\ttotal number of samples = %d", npures);

  printf("\n- cleaning up memory and closing files ...");
  cleanUpSensor(spec, st, th, sth);
 
  free(spec);
  free(st);
  free(th);
  free(sth);
  free(pars);

  printf("\n- program ends successfully\n");
  return SUCCESS;
}


int savePureSamples(SENSOR *spec, SENSOR *th, INPUT_PARS *pars) 
{
  char name[MAX_STRLEN];
  int16 ***dn, sdn, save[MAX_NBANDS];
  float sd, sd2;
  int   i, j, k, m, n, scale, bn, num;
  float sum, sum2, count, total, ave[MAX_NBANDS], ave_cv, weight;
  float ndvi, fc;
  FILE *th_sample, *th_names;

#ifdef DEBUG
  FILE *out[MAX_NBANDS];
  for(k=0; k<spec->nbands; k++) {
    sprintf(name, "%s_%dm.bin", spec->fileName[k], (int)th->res);
    out[k] = fopen(name, "wb");
  }
#endif
  
  sprintf(name, "%s.data", pars->cubistFile);
  th_sample = fopen(name, "w");
  sprintf(name, "%s.names", pars->cubistFile);
  th_names = fopen(name, "w");

  scale = spec->scale;
  alloc_3dim_contig((void ****) (&dn), scale, spec->ncols, spec->nbands, sizeof(int16));

  total = 0.0;  
  count = 0.0;

  for(i=0; i<th->nrows; i++) {

    if(i < pars->s_row || i > pars->e_row) continue;
    /* load lines to memory to compute mean values */
    for(m=i*scale; m<(i+1)*scale; m++) {
      if(m<spec->nrows)
	loadSensorRow(spec, m);
      for(n=0; n<spec->ncols; n++)
	for(k=0; k<spec->nbands; k++)
	  dn[m-i*scale][n][k] = spec->data[k][n]; 
    }
    
    for(j=0; j<th->ncols; j++) {
      
      if(j < pars->s_col || j > pars->e_col) continue;
      if(th->fdata[i][j]>=th->range[0] && th->fdata[i][j]<=th->range[1]) {
	sd2 = 0.0;
	bn = 0;	
	for(k=0; k<spec->nbands; k++) {
	  sum  = 0.0;
	  sum2 = 0.0;
	  num  = 0;
	  for(m=0; m<scale; m++)
	    for(n=0; n<scale; n++) 
	      if(j*scale+n < spec->ncols) {
		sdn = dn[m][j*scale+n][k];
		
		if(sdn != spec->fillValue && sdn>=spec->range[0] && sdn<=spec->range[1]) {
		  sum += sdn;
		  sum2 += sdn * sdn;
		  num++;  
		}
	      }
	  if(num==scale*scale) {
	    ave[k] = sum/num;
	    sd = sqrt(sum2/num - ave[k]*ave[k]); 
	    sd2 += sd/ave[k];
	    bn ++;
	    save[k] = (int)(ave[k]+0.5);
	  }
	  else
	    save[k] = spec->fillValue;
	}

	if(bn==spec->nbands) {
	  total++;
	  ave_cv = sd2/bn;
	  if(ave_cv < pars->CV_TH) {	
	    /* print out spectral band and thermal band values for pure pixels */
	    /* in irow, icol, icls, six spectral bands, ndvi, thermal_dn */
	    fprintf(th_sample, "%4d, %4d, ", i*scale, j*scale);

	    for(k=0; k<spec->nbands; k++)
	      fprintf(th_sample, "%6.1f, ", ave[k]);
	    if(ave[3]+ave[2] != 0) 
	      ndvi = (ave[3]-ave[2])/(ave[3]+ave[2]);
            else
	      ndvi = -1.0;
	    /* for airborne data test only */
	    /*ndvi = (ave[0]-ave[1])/(ave[0]+ave[1]);*/
	    fprintf(th_sample, "%7.4f, ", ndvi);
            if(ndvi<1.0)
	      fc = pow((1.0-ndvi), 0.625);
            else
	      fc = 0.0;
	    fprintf(th_sample, "%7.4f, ", fc);
	    fprintf(th_sample, "%7.3f, ", th->fdata[i][j]);
	    /* use inverse of coefficient of variance as weight */
	    if(ave_cv<0.01) weight = pars->CV_TH/0.01;
	    else weight = pars->CV_TH/ave_cv;	    
	    fprintf(th_sample, "%6.4f\n", weight);
	    count++;
	  }	
	}
      } /* end of if TIR is valid */
      else 
	for(k=0; k<spec->nbands; k++)
	  save[k] = spec->fillValue;
#ifdef DEBUG
      for(k=0; k<spec->nbands; k++)
	fwrite(&(save[k]), sizeof(int16), 1, out[k]);
#endif
    } /* end of icol */
  } /* end of irow */
  
  
  fprintf(th_names, "th.\n\n");
  fprintf(th_names, "irow: continuous\n");
  fprintf(th_names, "icol: continuous\n");

  for(k=0; k<spec->nbands; k++)
    fprintf(th_names, "B%d: continuous\n", k+1);
  
  fprintf(th_names, "ndvi: continuous\n");
  fprintf(th_names, "fc: continuous\n");
  fprintf(th_names, "th: continuous\n");
  fprintf(th_names, "case weight: continuous\n");
  
  fprintf(th_names, "attributes excluded: irow, icol, B1, ndvi, fc\n\n"); 
  

  fclose(th_sample);
  fclose(th_names);

#ifdef DEBUG
  for(k=0; k<spec->nbands; k++) {
    fclose(out[k]);
    sprintf(name, "%s_%dm.bin", spec->fileName[k], (int)th->res);
    writeENVIHeader(th, name, 1, sizeof(int16));
  }
#endif
  
  
  free_3dim_contig((void ***) dn);

  printf("\n\toriginal: nrows=%4d; ncols=%4d; res=%6.2f", spec->nrows, spec->ncols, spec->res);
  printf("\n\tcoarse:   nrows=%4d; ncols=%4d; res=%6.2f (resize factor =%d)", th->nrows, th->ncols, th->res, scale);
  printf("\n\tselected_nsamples = %d from %d (~%4.1f%%)", (int)count, (int)total, count*100.0/total);
  return count;
}
  

