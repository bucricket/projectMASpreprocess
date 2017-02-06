/**
 * this program includes sub-routines that read, write  
 * written and revised by Feng Gao on July 2012
 * v1.0_HDF: uses LEDAPS SR in HDF format (in-house LEDAPS product)
 * v1.1_ENVI: accepts GeoTIFF, ENVI and HDF format (USGS Landsat SR product)
 * v1.2_ENVI: uses Landsat cloud mask (USGS Landsat fmask product)
 **/

/* variables
   shortw: shortwave spectral bands 
   oth:  input temperature (could be at same resolution as shortwave bands
   cth:  coarse resolution T (could be same as st) aggregated from st
   sth : sharpened T
*/

#include "landsat.h"

/* read one row for all shortwave spectral bands */
int loadSensorRow(SENSOR *shortw, int irow)
{
  unsigned char *one_row_tiff;
  int i, j;
  long offset;

  alloc_1dim_contig((void **) (&one_row_tiff), shortw->ncols, sizeof(char));

  if(strcasecmp(shortw->fileType, "GEOTIFF")==0) {
    for(i=0; i<shortw->nbands; i++) {
      /* read data array from GeoTiff file (1 byte/pixel) */
      if (!TIFFReadScanline(shortw->fp_tiff[i], one_row_tiff, irow, 0)) {
	fprintf(stderr, "loadSensorRow: Read Sensor DN GeoTiff error\n");
	return FAILURE;
      }
      for(j=0; j<shortw->ncols; j++)
	shortw->data[i][j] = one_row_tiff[j];
    }
  }
  else {
    offset = (long) irow * shortw->ncols * sizeof(int16);
    for(i=0; i<shortw->nbands; i++) {
      fseek(shortw->fp[i], offset, 0); 
      /* read data array from binary file (2 byte short integer) */
      fread(shortw->data[i], sizeof(int16), shortw->ncols, shortw->fp[i]);
    }
  }

  free(one_row_tiff);
  return SUCCESS;
}


/* get metadata, allocate memory, and open files */
int getSensorMetaInfo(SENSOR *shortw, SENSOR *oth, SENSOR *cth, SENSOR *sth, INPUT_PARS *pars) 
{
  int i;
  
  /* set parameters from input */
  shortw->nbands = pars->nbands;
  shortw->nrows  = pars->nrows;
  shortw->ncols  = pars->ncols;
  shortw->fillValue = pars->fillv;
  for(i=0; i<2; i++) 
    shortw->range[i] = pars->range[i];
  shortw->res = pars->res;
  shortw->ulx = pars->ulx;
  shortw->uly = pars->uly;
  shortw->insys  = pars->insys;
  shortw->inzone = pars->inzone;
  shortw->inunit = pars->inunit;
  for(i=0; i<16; i++)
    shortw->inparm[i] = pars->inparm[i];
  shortw->indatum = pars->indatum;
  strcpy(shortw->fileType, pars->fileType);
  /* open shortw files (one file per band) and update metadata if it is a GeoTIFF file */
  for(i=0; i<shortw->nbands; i++) {
    strcpy(shortw->fileName[i], pars->inFile[i]);
    if(getMetaDataFromGeoTIFF(shortw, i)==FAILURE) {
      fprintf(stderr, "\ngetMetaDataFromGeoTIFF error\n");
      return FAILURE;
    }
  }

  /* set original surface temperature (oth) */
  oth->nbands = 1;
  if(pars->org_th_nrows != -1) { 
    oth->nrows = pars->org_th_nrows;
    oth->ncols = pars->org_th_ncols;
  }
  else {
    /* if not define, use nrows and ncols from shortwave bands */ 
    oth->nrows = shortw->nrows;
    oth->ncols = shortw->ncols;
  }
  oth->res = pars->org_th_res;
  oth->fillValue = 0;
  oth->range[0] = pars->org_th_range[0];
  oth->range[1] = pars->org_th_range[1];
  strcpy(oth->fileType, pars->org_th_fileType);
  /* open original ST file and update info if it is a GeoTIFF file */
  strcpy(oth->fileName[0], pars->org_th_File);
  if(getMetaDataFromGeoTIFF(oth, 0) == FAILURE) {
      fprintf(stderr, "\ngetMetaDataFromGeoTIFF error\n");
      return FAILURE;
  }

  /* open cloud mask file */
  oth->CMASK_FLAG = 0;
  if(strcasecmp(pars->fileType, "binary") == 0) {
    if((oth->fp_cloud = fopen(pars->cloudFile, "rb")) != NULL) 
      oth->CMASK_FLAG = 1;
  }
  else if(strcasecmp(pars->fileType, "GeoTiff") == 0) {
    if((oth->fp_tiff_cloud = XTIFFOpen(pars->cloudFile, "r")) != NULL) 
      oth->CMASK_FLAG = 1;
  } 
  
  /* set reduced and rescaled coarse resolution T (cth) */
  cth->insys = pars->insys;
  cth->res = pars->res_th_res;
  cth->ulx = shortw->ulx;
  cth->uly = shortw->uly;
  cth->inzone = pars->inzone;
  cth->nbands = 1;
  if( fabs(cth->res - oth->res) < 0.00001) {
    cth->nrows = oth->nrows;
    cth->ncols = oth->ncols;
  }
  else {
    /* if original and coarse T has different spatial resolution */
    cth->nrows = (int)(oth->nrows * oth->res / cth->res) + 1;
    cth->ncols = (int)(oth->ncols * oth->res / cth->res) + 1;
  }
  cth->fillValue = 0;
  cth->range[0] = oth->range[0];
  cth->range[1] = oth->range[1];
  /* this will be a temp file in binary format */
  strcpy(cth->fileType, "binary");

  /* initialize sharpened T (output) in same dimension as shortwave bands */ 
  sth->nbands = 1;
  sth->nrows = shortw->nrows;   
  sth->ncols = shortw->ncols;
  sth->res = shortw->res;
  sth->fillValue = 0;
  sth->ulx = shortw->ulx;
  sth->uly = shortw->uly;
  sth->insys = shortw->insys;
  sth->inzone = shortw->inzone;
  sth->inunit = shortw->inunit;
  sth->range[0] = oth->range[0];
  sth->range[1] = oth->range[1];
  for(i=0; i<16; i++)
    sth->inparm[i] = shortw->inparm[i];
  sth->indatum = shortw->indatum;
  /* only write in binary format */
  strcpy(sth->fileType, "Binary");

  /* compute zoom-in scale to coarse T */
  /* original T and sharpened T could be different or same. For Landsat, thermal bands have been
     resampled to shortwave band resolution, so they are same. However, if coarse T (60m for ETM+ 
     or 120 m for TM) is provided, the two scales are different. 
     scale for shortwave R always equals to the sharpened T */
  shortw->scale = (int)(cth->res/shortw->res+0.0001);
  oth->scale = (int)(cth->res/oth->res+0.0001);
  sth->scale = (int)(cth->res/sth->res+0.0001);
  cth->scale = 1;
  
  i = 0;   /* only accepts one thermal band */
  strcpy(sth->fileName[0], pars->outFile);
 
  /* allocate memory to load one row for shortw data (shortwave bands) */
  alloc_2dim_contig((void ***) (&shortw->data), shortw->nbands, shortw->ncols, sizeof(int16));

  /* allocate memory to store the entire coarse thermal image) */
  alloc_2dim_contig((void ***) (&cth->fdata), cth->nrows, cth->ncols, sizeof(float));

  /* allocate memory for the entire original thermal image */
  alloc_2dim_contig((void ***) (&oth->fdata), oth->nrows, oth->ncols, sizeof(float));
  alloc_2dim_contig((void ***) (&oth->cloud), oth->nrows, oth->ncols, sizeof(uint8));
 
  /* allocate memory for the entire sharpened thermal image (T and qa) */
  alloc_2dim_contig((void ***) (&sth->fdata), sth->nrows, sth->ncols, sizeof(float));
  alloc_2dim_contig((void ***) (&sth->qa), sth->nrows, sth->ncols, sizeof(uint8));

#ifdef SAVE_COARSE_T
  /* open files for writing */
  sprintf(cth->fileName[0], "%s_%dm.bin", oth->fileName[0], (int)cth->res);
  if((cth->fp[0] = fopen(cth->fileName[0], "wb"))==NULL) {
    fprintf(stderr, "\ngetSensorMetaInfo: Error in open thermal file as binary\n");
    return FAILURE;
  }
#endif
 
  return SUCCESS;
}


/* open sharpened T for write */
int openForWrite(SENSOR *sth)
{
  int i,j;

  /* if it is exist then read it so the result from next subwindow will be overwrite to it */
  if((sth->fp[0] = fopen(sth->fileName[0], "rb")) != NULL) {
    for(i=0; i<sth->nrows; i++)
      for(j=0; j<sth->ncols; j++)
	fread(&(sth->fdata[i][j]),sizeof(float), 1, sth->fp[0]); 
    /*fprintf(stderr, "\n\tread from previous prediction");*/
    fclose(sth->fp[0]);
  }
  else {
    for(i=0; i<sth->nrows; i++)
      for(j=0; j<sth->ncols; j++) {
	sth->fdata[i][j] = sth->fillValue;
	sth->qa[i][j] = -1;
      }
  }

  if((sth->fp[0] = fopen(sth->fileName[0], "wb"))==NULL) {
    fprintf(stderr, "\ngetSensorMetaInfo: Error in open thermal file as binary\n");
    return FAILURE;
  }

  return SUCCESS;
}


/* load the entire thermal image and compute min and max values */
int loadThermal(SENSOR *th)
{
  unsigned char *one_row_tiff;
  int i, j;

  alloc_1dim_contig((void **) (&one_row_tiff), th->ncols, sizeof(char));
  if(strcasecmp(th->fileType, "binary")==0) rewind(th->fp[0]);

  for(i=0; i<th->nrows; i++) {

    /* initialize cloud mask (0 or 1 isclear) */
    for(j=0; j<th->ncols; j++)
      th->cloud[i][j] = 0;
    
    if(strcasecmp(th->fileType, "GEOTIFF")==0) { 
      /* read data array from GeoTiff file */
      if (!TIFFReadScanline(th->fp_tiff[0], one_row_tiff, i, 0)) {
        fprintf(stderr, "loadSensorRow: Read Sensor DN GeoTiff error\n");
        return FAILURE;
      } 
      for(j=0; j<th->ncols; j++)
	th->fdata[i][j] = one_row_tiff[j];

      /* read cloud mask from GeoTiff file */
      if(th->CMASK_FLAG == 1) {
	if (!TIFFReadScanline(th->fp_tiff_cloud, one_row_tiff, i, 0)) {
	  fprintf(stderr, "loadSensorRow: Read Cloud Mask from GeoTiff error\n");
	  return FAILURE;
	} 
	for(j=0; j<th->ncols; j++)
	  th->cloud[i][j] = one_row_tiff[j];
      }
      
    }
    else {
      /* read data array from binary file */
      fread(th->fdata[i], sizeof(float), th->ncols, th->fp[0]);
      if(th->CMASK_FLAG == 1)
	fread(th->cloud[i], sizeof(uint8), th->ncols, th->fp_cloud);
    }
  }

  th->min[0] = th->range[1];
  th->max[0] = th->range[0];
  for(i=0; i<th->nrows; i++)
    for(j=0; j<th->ncols; j++) {

      /* Landsat fmask (0=clear_land; 1=clear_water; 2=cloud_shadow; 3=snow; 4=cloud; 255=no_obs */
      if(th->cloud[i][j] > 1) th->fdata[i][j] = FILLV; 
      if(th->fdata[i][j]>=th->range[0]&&th->fdata[i][j]<=th->range[1]) {
	if(th->fdata[i][j]<th->min[0]) th->min[0] = th->fdata[i][j];
	if(th->fdata[i][j]>th->max[0]) th->max[0] = th->fdata[i][j];
      }      
    }

  printf("\n\tmin=%5.1f; max=%5.1f", th->min[0], th->max[0]);
  if(strcasecmp(th->fileType, "binary")==0) {
    rewind(th->fp[0]);
    if(th->CMASK_FLAG == 1)
      rewind(th->fp_cloud);
  }
  free(one_row_tiff);
  return SUCCESS;
}


/* rescale original T to coarse T */
int resize(SENSOR *oth, SENSOR *cth)
{
  int i, j, ti, tj, m, n, scale, num;
  float sum, min, max, f;
 
  cth->range[0] = oth->range[0];
  cth->range[1] = oth->range[1];
  cth->fillValue = oth->fillValue;
  scale = oth->scale;

  printf("\n\toriginal: nrows=%4d; ncols=%4d; res=%6.2f", oth->nrows, oth->ncols, oth->res);
  printf("\n\tcoarse:   nrows=%4d; ncols=%4d; res=%6.2f (resize factor =%d)", cth->nrows, cth->ncols, cth->res, scale);
 
  /* find min and max from aggregated data */
  min = oth->range[1];
  max = oth->range[0];

  for(i=0; i<cth->nrows; i++) {
    for(j=0; j<cth->ncols; j++) {

      sum = 0.0;
      num = 0;
      f = oth->fillValue;

      for(m=0; m<scale; m++)
	for(n=0; n<scale; n++) {
	  ti = i*scale + m;
	  tj = j*scale + n;
	  if(ti>=oth->nrows || tj>=oth->ncols) continue;
	  // both shortwave and thermal band data should be valid  
	  // otherwise, energy conservation will be broken  
	  if(oth->fdata[ti][tj]>=oth->range[0]&&oth->fdata[ti][tj]<=oth->range[1]) {
	    sum += pow(oth->fdata[ti][tj],4);
	    num ++;
	  }
	}
      if(num>=scale*scale*EC_VALID_TH) {
	f = pow(sum/num, 0.25);
	if(f < min) min = f;
	if(f > max) max = f;
      }     
      /*else
	if(i*scale<oth->nrows && j*scale<oth->ncols)
	  f = oth->fdata[i*scale][j*scale];*/
      
      cth->fdata[i][j] = f;
#ifdef SAVE_COARSE_T
      if(cth->fp[0]!=NULL)
	fwrite(&f, sizeof(float), 1, cth->fp[0]);
#endif
    }
  }
  
#ifdef SAVE_COARSE_T
  printf("\nagg_min=%f; agg_max=%f", min, max);
  if(cth->fp[0]!=NULL)
    writeENVIHeader(cth, cth->fileName[0], 1, 4); 
#endif

  return SUCCESS;
}



/* close opened Landsat sds and file, free memory */
int cleanUpSensor(SENSOR *shortw, SENSOR *oth, SENSOR *cth, SENSOR *sth)
{
  int i;

  /* close shortwave bands */
  for(i=0; i<shortw->nbands; i++) {
    if(strcasecmp(shortw->fileType, "binary") == 0) 
      fclose(shortw->fp[i]);
    else if(strcasecmp(shortw->fileType, "GeoTiff") == 0)
      XTIFFClose(shortw->fp_tiff[i]);
  }

  /* close original thermal band */
  if(strcasecmp(oth->fileType, "binary") == 0) 
    fclose(oth->fp[0]);
  else if(strcasecmp(oth->fileType, "GeoTiff") == 0) 
    XTIFFClose(oth->fp_tiff[0]);

#ifdef SAVE_COARSE_T
  /* clean up for coarse thermal band */
  if(strcasecmp(cth->fileType, "binary") == 0) 
    fclose(cth->fp[0]);
  else if(strcasecmp(cth->fileType, "GeoTiff") == 0) 
    XTIFFClose(cth->fp_tiff[0]);
#endif

  if(sth->fp[0]!=NULL) 
    fclose(sth->fp[0]);

  free_2dim_contig((void **) shortw->data);

  free_2dim_contig((void **) cth->fdata);
  free_2dim_contig((void **) oth->fdata);
  free_2dim_contig((void **) oth->cloud);

  free_2dim_contig((void **) sth->fdata);
  free_2dim_contig((void **) sth->qa);

  return SUCCESS;
}



/* apply energy conservation to sharpened T */
int applyEC(SENSOR *shortw, SENSOR *oth, SENSOR *cth, SENSOR *sth, INPUT_PARS *pars)
{  
  int i, j, m, n, num, scale;
  int ii, jj, si, sj, flag, oflag, sflag;
  float x, y, dx, dy, v[2][2], s_diff;
  float  **diff;
  double sum;

#ifdef DEBUG
  float org_diff, pdata, odata;
  double  asum, sum1, sum2, sumx1, sumx2, sumy1, sumy2, sumxy, r, stdev1, stdev2;
  double npixs;
#endif

  /* write cubist prediction in BSQ format */
  for(i=0; i<sth->nrows; i++)
    for(j=0; j<sth->ncols; j++) {
      /* set cloud pixel to fill value (0 or 1 is clear) */
      if(oth->cloud[i][j]>1) 
	sth->fdata[i][j] = sth->fillValue;
      fwrite(&(sth->fdata[i][j]), sizeof(float), 1, sth->fp[0]);
    }

#ifdef DEBUG
  /* write the difference between original TIR and cubist prediction at original resolution */
  asum = 0.0;  sum1 = 0.0;  sum2 = 0.0;
  sumx1 = 0.0; sumx2 = 0.0;
  sumy1 = 0.0; sumy2 = 0.0; sumxy = 0.0;
  npixs = 0;
  for(i=0; i<sth->nrows; i++)
    for(j=0; j<sth->ncols; j++) {
      pdata = sth->fdata[i][j];
      odata = oth->fdata[i][j];
      if(odata<oth->range[0]||odata>oth->range[1]||pdata<oth->range[0]||pdata>oth->range[1]||isnan(odata)||isnan(pdata))
	org_diff = FILLV;
      else {
	org_diff = pdata-odata;
	asum += fabs(org_diff);
	sum1 += org_diff;
	sum2 += org_diff*org_diff;
	sumx1 += oth->fdata[i][j];
	sumx2 += oth->fdata[i][j] * oth->fdata[i][j];
	sumy1 += sth->fdata[i][j];
	sumy2 += sth->fdata[i][j] * sth->fdata[i][j];
	sumxy += oth->fdata[i][j] * sth->fdata[i][j]; 
	npixs++;
      }
      fwrite(&(org_diff), sizeof(float), 1, sth->fp[0]);
    }
  stdev1 = sqrt(sumx2/npixs - sumx1/npixs * sumx1/npixs);
  stdev2 = sqrt(sum2/npixs - sum1/npixs * sum1/npixs);
  r = (npixs*sumxy-sumx1*sumy1)/sqrt((npixs*sumx2-sumx1*sumx1)*(npixs*sumy2-sumy1*sumy1));
  printf("\n\nModel Prediction: Mean=%6.3f STDEV=%6.3f MAE=%6.3f MBE=%6.3f MBE_STDEV=%6.3f rmse=%6.3f r=%6.3f\n", (float)sumx1/npixs, stdev1, (float)asum/npixs, (float)sum1/npixs, stdev2, sqrt(sum2/npixs), r);

#endif

  /* sharpened image has same resolution as shortwave spectral bands */
  scale = sth->scale;

  alloc_2dim_contig((void ***) (&diff), cth->nrows, cth->ncols, sizeof(float));

  /* compute difference from aggregated temperature */
  for(i=0; i<cth->nrows; i++)
    for(j=0; j<cth->ncols; j++) {

      sum = 0.0;
      num = 0;

      flag = 1;
      for(m=i*scale; m<(i+1)*scale; m++)
	for(n=j*scale; n<(j+1)*scale; n++) {
	  if(m>=sth->nrows || n>=sth->ncols) continue;

	  // predicted data (so shortwave bands) and thermal data (TIR band) must be consistent
	  // found ETM+ gaps are not consistent between shortwave and TIR band 
	  if(sth->fdata[m][n]>=sth->range[0]&&sth->fdata[m][n]<=sth->range[1]) sflag = 1;
	  else sflag = 0;
	  if(oth->fdata[m][n]>=oth->range[0]&&oth->fdata[m][n]<=oth->range[1]) oflag = 1;
	  else oflag = 0;
	  if(sflag*oflag == 0) {flag = 0; break;}
	  if(sflag == 1 && oflag == 1) {
	    sum += pow(sth->fdata[m][n], 4);
	    num++;
	  }
	}
      
      if(num < scale*scale*EC_VALID_TH || flag==0 || cth->fdata[i][j]==cth->fillValue) 
	diff[i][j] = 0;
      else 
	/* mean difference of energy */
	diff[i][j] = cth->fdata[i][j] - pow(sum/num, 0.25);
      /* if((i>=100/2-1&&i<=100/2+1)  && (j>=1445/2-1&& j<=1445/2+1)) 
	printf("\n(%d %d) %6.4f %6.4f %6.4f", i, j, cth->fdata[i][j], pow(sum/num, 0.25), diff[i][j]);
      */      
    } /* end of image loop */
  
  /* apply residuals */
  for(i=0; i<sth->nrows; i++) {
    for(j=0; j<sth->ncols; j++) {
      /* convert to coarse space (x,y)*/
      y = (float)(i+0.5)/scale;
      x = (float)(j+0.5)/scale;
      /* actual pixel location (m,n) */
      m = (int) y;
      n = (int) x; 
      /*if(m>=cth->nrows||n>=cth->ncols) continue; */    
      if(pars->SMOOTH_FLAG) {
	/* apply bilinear interoplation on difference */
	dy = y-(m+0.5);
	dx = x-(n+0.5);
	/* decide upper-left pixel */
	if(dy<0) si=m-1; else si=m;
	if(dx<0) sj=n-1; else sj=n;
	flag = 1;
	for(ii=0; ii<2; ii++)
	  for(jj=0; jj<2; jj++)
	    if(si+ii>0&&si+ii<cth->nrows&&sj+jj>0&&sj+jj<cth->ncols)
	      v[ii][jj] = diff[si+ii][sj+jj];    
	    else
	      flag = 0;
	if(flag != 0) {
	  /* dx and sy should be always positive here */
	  dy = y-(si+0.5);
	  dx = x-(sj+0.5);
	  s_diff = v[0][0]*(1-dx)*(1-dy)+v[0][1]*dx*(1-dy)+v[1][0]*(1-dx)*dy+v[1][1]*dx*dy;
	}
	else
	  s_diff = diff[m][n];
	/*if(i==100 && j==1445) {
	  printf("\n(%d %d) %6.4f\n", i, j, s_diff);
	  for(ii=0; ii<2; ii++)
	    for(jj=0; jj<2; jj++)
	      printf("%6.4f ", v[ii][jj]);
	      }*/
	
      }
      else
	s_diff = diff[m][n];

#ifdef DEBUG
      /* write the difference (either smoothed or not) between TIR and cubist prediction at coarse resolution */
      /* exclude invalid value in difference image */
      /* if(cth->fdata[m][n]<cth->range[0]||cth->fdata[m][n]>cth->range[1])
	s_diff = FILLV;
	fwrite(&(s_diff), sizeof(float), 1, sth->fp[0]);      */
#endif
      if(sth->fdata[i][j]<oth->range[0]||sth->fdata[i][j]>oth->range[1]||oth->cloud[i][j]>1)
	/*sth->fdata[i][j] = FILLV;*/
	sth->fdata[i][j] = oth->fdata[i][j];
      else
	sth->fdata[i][j] += s_diff;
    }
  }

  /* final prediction with residual applied */
  for(i=0; i<sth->nrows; i++)
    for(j=0; j<sth->ncols; j++) 
      fwrite(&(sth->fdata[i][j]), sizeof(float), 1, sth->fp[0]);
  
#ifdef DEBUG
  /* write difference between original TIR and final prediction */
  asum = 0.0;  sum1 = 0.0;  sum2 = 0.0;
  sumx1 = 0.0; sumx2 = 0.0;
  sumy1 = 0.0; sumy2 = 0.0; sumxy = 0.0;
  npixs = 0;
  for(i=0; i<sth->nrows; i++)
    for(j=0; j<sth->ncols; j++) {
      pdata = sth->fdata[i][j];
      odata = oth->fdata[i][j];
      if(odata<oth->range[0]||odata>oth->range[1]||pdata<oth->range[0]||pdata>oth->range[1]||isnan(odata)||isnan(pdata))
	org_diff = FILLV;
      else {
	org_diff = pdata-odata;
	asum += fabs(org_diff);
	sum1 += org_diff;
	sum2 += org_diff*org_diff;
	sumx1 += oth->fdata[i][j];
	sumx2 += oth->fdata[i][j] * oth->fdata[i][j];
	sumy1 += sth->fdata[i][j];
	sumy2 += sth->fdata[i][j] * sth->fdata[i][j];
	sumxy += oth->fdata[i][j] * sth->fdata[i][j]; 
	npixs++;
      }
      fwrite(&(org_diff), sizeof(float), 1, sth->fp[0]);
    }
  stdev1 = sqrt(sumx2/npixs - sumx1/npixs * sumx1/npixs);
  stdev2 = sqrt(sum2/npixs - sum1/npixs * sum1/npixs);
  r = (npixs*sumxy-sumx1*sumy1)/sqrt((npixs*sumx2-sumx1*sumx1)*(npixs*sumy2-sumy1*sumy1));
  printf("Residual Applied: Mean=%6.3f STDEV=%6.3f MAE=%6.3f MBE=%6.3f MBE_STDEV=%6.3f rmse=%6.3f r=%6.3f\n", (float)sumx1/npixs, stdev1, (float)asum/npixs, (float)sum1/npixs, stdev2, sqrt(sum2/npixs), r);

  /* write out: 
     1) cubist prediction
     2) difference between orginal fine T and cubist model prediction, 
     3) residuals applied sharpened T
     4) difference between original fine T and final prediction */
  writeENVIHeader(sth, pars->outFile, 4, 4); 

#else
  writeENVIHeader(sth, pars->outFile, 2, 4); 
#endif
  
  free_2dim_contig((void **) diff);
  return SUCCESS;
}


/* retrieve metadata from GeoTIFF file and open it */
int getMetaDataFromGeoTIFF(SENSOR *sensor, int iband)
{
  int i;
  uint16 count, coor_sys;
  double *tiePoint, *pixelScale;
  GTIF *gtif;

  i = iband;
  if(strcasecmp(sensor->fileType, "binary") == 0) {
    if((sensor->fp[i] = fopen(sensor->fileName[i], "rb"))==NULL) {
      fprintf(stderr, "\ngetSensorMetaInfo: Error in open input file as binary\n");
      return FAILURE;
    }
  }
  else if(strcasecmp(sensor->fileType, "GeoTiff") == 0) {
    if((sensor->fp_tiff[i] = XTIFFOpen(sensor->fileName[i], "r"))==NULL) {
      fprintf(stderr, "\ngetSensorMetaInfo: Error in open input file as GeoTiff\n");
      return FAILURE;
    }
    
    /* update metadata if they can be found in GeoTiff file */
    if(TIFFGetField(sensor->fp_tiff[i], TIFFTAG_IMAGEWIDTH, &(sensor->ncols))==0) {
      fprintf(stderr, "\nRetrieve image width error: getSensorMetaInfo\n");
      return FAILURE;
    }
    if(TIFFGetField(sensor->fp_tiff[i], TIFFTAG_IMAGELENGTH, &(sensor->nrows))==0) {
      fprintf(stderr, "\nRetrieve image length error: getSensorMetaInfo\n");
      return FAILURE;
    }
    count=6;
    if(TIFFGetField(sensor->fp_tiff[i], TIFFTAG_GEOTIEPOINTS, &count, &tiePoint)==0) {
      fprintf(stderr, "\nRetrieve tiePoint error: getSensorMetaInfo\n");
      return FAILURE;
    }
    count=3;
    if(TIFFGetField(sensor->fp_tiff[i], TIFFTAG_GEOPIXELSCALE, &count, &pixelScale)==0) {
      fprintf(stderr, "\nRetrieve pixelScale file error: getSensorMetaInfo\n");
      return FAILURE;
    }
      
      sensor->res = pixelScale[0];

      /* GeoKey 1025 (GTRasterTypeGeoKey) dictates whether the reference
	 coordinate is the UL (*RasterPixelIsArea*, code 1) or center
	 (*RasterPixelIsPoint*, code 2) of the UL pixel. If this key is missing,
	 the default (as defined by the specification) is to be
	 *RasterPixelIsArea*, which is the UL of the UL pixel. */
      gtif = GTIFNew(sensor->fp_tiff[i]);
      if (GTIFKeyGet(gtif, GTRasterTypeGeoKey, &coor_sys, 0, 1) != 1) {
	printf("Coordinate system is not defined in %s\n", sensor->fileName[i]);
	printf("assume used UL of the UL pixel\n");
      }
      if (coor_sys == RasterPixelIsPoint){
	sensor->ulx = tiePoint[3] - 0.5 * sensor->res;
	sensor->uly = tiePoint[4] + 0.5 * sensor->res;
      }
      else {  /* default use RasterPixelIsArea */
	sensor->ulx = tiePoint[3];
	sensor->uly = tiePoint[4];
      }
      GTIFFree(gtif);

#ifdef DEBUG
      if(i==0) {
	printf("\nExtracted information from %s\n", sensor->fileName[i]);
	printf("\nNCOLS = %d", sensor->ncols);
	printf("\nNROWS = %d", sensor->nrows);
	printf("\nPIXEL_SIZE = %f", sensor->res);
	printf("\nUPPER_LEFT_CORNER = %f %f\n", sensor->ulx, sensor->uly);
      }
#endif
    }
  else {
    fprintf(stderr, "\nUnknown file type %s: in getSensorMetaInfo \n", sensor->fileType);
    return FAILURE;
  }

  return SUCCESS;
}


