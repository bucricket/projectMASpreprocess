/**
 * This program reads TIMESAT smoothed data and fills data gap
 * using temporal curve from neighbor pixels with same land cover 
 * type and then adjusts curve to current pixel values.
 *
 * Inputs for this program includes
 * - timesat smoothed data set (e.g. fitAG_*)
 * - HDF input lists used for timesat program
 * - MODIS land cover data
 *
 * For processing on LAI data by Feng Gao on 11/30/2006
 * Revision 2/26/2007: compress HDF output 
 * Revision 8/7/2012: use MCD12Q1 (500m, 2400*2400) instead of MOD12Q1 (1km, 1200*1200)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include "hdf.h"
#include "mfhdf.h"
#include "HdfEosDef.h"
typedef long long off64_t;

#define Max_StrLen 1000
#define FAILURE -1
#define SUCCESS 1
#define MAX_OUTS 512
#define LAI_INVALID 249 /* 249 and up */
#define MAX_STEP 60
#define ACCEPTABLE_ERR 20
#define ACCEPTABLE_DIFF 40
uint8 laiRange[2] = {0, 100};
#define NUM_CLASSES 20
#define NUM_LAYERS  9
#define TIMESAT_LAI_FV -100

#define COM_ORG_HQ       1
#define COM_FIT          2
#define COM_FILL         3

#define FIT_TIMESAT      1
#define FIT_GAPFILLED    2
#define FIT_ROUNDED      3
#define FIT_FILL         4

#define ORG_HQ_TS_HQ     1
#define ORG_HQ_TS_LQ     2
#define ORG_BACKUP       3
#define ORG_FILL         4

#define WEIGHT_TS_HQ     0.667
#define WEIGHT_TS_LQ     0.251
#define WEIGHT_BACKUP    0.249
#define WEIGHT_VALID     0.001

/*#define ShortName "MOD15COM"
#define GridName  "MODIS_NACP_LAI"*/
#define ShortName "MOD15A2"
#define GridName  "MOD_Grid_MOD15A2"
#define FitFileInBIP "temp_fitting_bip_file.bin"

/*#define DEBUG*/
#ifdef DEBUG
  #define DEBUG_icol 564
  #define DEBUG_irow 292
#endif

#define ERROR(message, module) \
          Error((message), (module), (__FILE__), (long)(__LINE__))

typedef struct {
  char job_name[Max_StrLen];
  char lai_fn[Max_StrLen];
  char fpar_fn[Max_StrLen];
  char weight_fn[Max_StrLen];
  char input_fn[Max_StrLen];
  char lai_list[Max_StrLen];
  char fit_fn[Max_StrLen];
  char lc_fn[Max_StrLen];
  char biome_fn[Max_StrLen];
  char lut_fn[Max_StrLen];
  char temp_fn[Max_StrLen];
  int map_start;
  int map_end;
  int nyears, nptperyear, tnpt;
  int rowstart, rowstop, colstart, colstop;
  int nrows, ncols;
  uint8 lut[10][100];
} PARAMETERS;

typedef struct {
  uint8   fillValue;
  uint8   range[2];
  float64 scaleFactor;
  char  name[Max_StrLen];
  char  longName[Max_StrLen];
  char  units[Max_StrLen];
  char  explain[Max_StrLen];
} SDS;

typedef struct {
  FILE  *lai;    /* original MODIS LAI */
  FILE  *fpar;   /* original MODIS FPAR */
  FILE  *qc;     
  FILE  *slai;   /* smoothed LAI */
  FILE  *composed_lai[MAX_OUTS];

  int32 nrows;           
  int32 ncols;

  /* image metadata */
  float64 GD_upleft[2];
  float64 GD_lowright[2];
  int32   GD_projcode;
  int32   GD_zonecode;
  int32   GD_spherecode;
  float64 GD_projparm[16];
  int32   GD_origincode;

  int32 lc_GDfid;
  int32 lc_GDid;
  int32 biome_GDfid;
  int32 biome_GDid;
  int32 out_GDfid[MAX_OUTS];
  int32 out_GDid[MAX_OUTS];

  SDS modis;
  SDS timesat;
  SDS compose;
  SDS modis_fpar;
  SDS lut_fpar;
  SDS compose_fpar;
  SDS modis_qc;
  SDS timesat_qc;
  SDS compose_qc;
  
  char allNames[NUM_LAYERS][Max_StrLen];

} GRID; 


int openFiles(PARAMETERS *par, GRID *grid);
int fillGaps(PARAMETERS *par, GRID *grid);
int computeFillfromLocal(PARAMETERS *par, float cdata[], float weight[], float ndata[], float pdata[]); 
int computeFillfromAll(PARAMETERS *par, float cdata[], float weight[], float ndata[], float pdata[]);
int lai2fpar(int lai, int lc, PARAMETERS *par);
int openLUT(PARAMETERS *par, GRID *grid);
int closeFiles(PARAMETERS *par, GRID *grid);
int convertFittingFile(PARAMETERS *par); 
int parseParameters(PARAMETERS *par, int argc, char *argv[]);
int getModisMetaInfo(char *fname, GRID *g);
int writeMetaInfo(char *fname, GRID *g, int index);
void usage(char *);
void dec2bin(uint32 num, int *bitpattern, int limit);
void bin2dec (uint32 *num, int bitpattern[], int limit);
void getAffine(float cp[][3], int num_cps, int n, float a[]);
void solveEquation(float s[][10], float ss[], int m, float a[]);

void alloc_1dim_contig (void **, int, int);
void alloc_2dim_contig (void ***, int, int, int);
void alloc_3dim_contig (void ****, int, int, int, int);
void free_2dim_contig (void **);
void free_3dim_contig (void ***);
void Error(const char *message, const char *module,
           const char *source, long line);

int compare(const void *f1, const void *f2);

main(int argc, char *argv[])
{
  PARAMETERS par;
  GRID grid;

  /* retrieve input parameters */
  if(parseParameters(&par, argc, argv)==FAILURE) {
    ERROR("parse parameters error", "main");
    exit(1);
  }

  /* convert timesat fitting file to BIP image format for easy operate */
  printf("Converting timesat fitting file to BIP format ...\n");
  if(convertFittingFile(&par)==FAILURE) {
    ERROR("convert fitting function data error", "main");
    exit(1);
    }

  /* open input files and create output files */
  printf("Opening files for read and write\n");
  if(openFiles(&par, &grid)==FAILURE) {
    ERROR("open files for read and write error", "main");
    exit(1);
  }

  /* fill spatial gaps for each invalid value */
  printf("Doing gap filling\n");
  if(fillGaps(&par, &grid) == FAILURE) {
    ERROR("gap fill error", "main");
    exit(1);
  }

  /* close all inputs and outputs */
  printf("\nClosing all files\n");
  if(closeFiles(&par, &grid) == FAILURE) {
    ERROR("close files error", "main");
    exit(1);
  }

}


int fillGaps(PARAMETERS *par, GRID *grid)
{
  int32 start[2];
  int32 stride[2];
  int32 length[2];
  uint8 bip[NUM_LAYERS], *sdata, **can_data, ***oneline, *buffer; 
  uint8 **igbp, **biome, **igbp_rows, **biome_rows, igbp_lc, lc;
  uint8 *lai_char, *lai_smooth, *lai_composed;
  uint8 *fpar_char, *fpar_lut, *fpar_composed;
  uint8 *org_qc, *fit_qc, *com_qc;
  uint8 **smooth_flag; 
  int8  *data, **data_lc, temp; 
  float *cdata, *weight, *cweight, *ndata, *pdata, **can_weight;
  float diff, ave_diff;
  float min_top, max_bottom;
  int num, hq_num, cnum, max_hq_num, best;
  int i, j, k, m, n, ri, rj, step, valid_flag;
  int num_valid, max_freq;
  int hq_num_lc[NUM_CLASSES], tot_ave[NUM_CLASSES], freq[NUM_CLASSES];
  double **ave_data_lc;

  off64_t offset;
  /* FILE *out;
     out=fopen("lc_agg_1km.bin", "wb");*/

  printf("\tAllocating memory");
  alloc_1dim_contig((void **)(&org_qc), par->tnpt, sizeof(uint8));
  alloc_1dim_contig((void **)(&fit_qc), par->tnpt, sizeof(uint8));
  alloc_1dim_contig((void **)(&com_qc), par->tnpt, sizeof(uint8));
  alloc_1dim_contig((void **)(&data), par->tnpt, sizeof(int8));
  alloc_1dim_contig((void **)(&sdata), par->tnpt, sizeof(uint8));  
   alloc_1dim_contig((void **)(&buffer), par->ncols, sizeof(uint8)); 
  alloc_2dim_contig((void ***)(&data_lc), NUM_CLASSES, par->tnpt, sizeof(int8));
  alloc_2dim_contig((void ***)(&ave_data_lc), NUM_CLASSES, par->tnpt, sizeof(double));
  alloc_2dim_contig((void ***)(&igbp_rows), 2, 2400, sizeof(uint8)); 
  alloc_2dim_contig((void ***)(&biome_rows), 2, 2400, sizeof(uint8)); 
  alloc_2dim_contig((void ***)(&igbp), 1200, 1200, sizeof(uint8));
  alloc_2dim_contig((void ***)(&biome), 1200, 1200, sizeof(uint8));
  alloc_2dim_contig((void ***)(&smooth_flag), par->nrows, par->ncols, sizeof(uint8));    
  alloc_2dim_contig((void ***)(&can_data), MAX_STEP*MAX_STEP, par->tnpt, sizeof(uint8));  
  alloc_1dim_contig((void **)(&lai_char),  par->tnpt, sizeof(uint8)); 
  alloc_1dim_contig((void **)(&lai_smooth), par->tnpt, sizeof(uint8));  
  alloc_1dim_contig((void **)(&lai_composed), par->tnpt, sizeof(uint8));  
  alloc_1dim_contig((void **)(&fpar_char), par->tnpt, sizeof(uint8)); 
  alloc_1dim_contig((void **)(&fpar_lut), par->tnpt, sizeof(uint8));  
  alloc_1dim_contig((void **)(&fpar_composed), par->tnpt, sizeof(uint8));  
  alloc_1dim_contig((void **)(&cdata), par->tnpt, sizeof(float));
  alloc_1dim_contig((void **)(&weight), par->tnpt, sizeof(float));
  alloc_1dim_contig((void **)(&cweight), par->tnpt, sizeof(float));
  alloc_1dim_contig((void **)(&ndata), par->tnpt, sizeof(float));
  alloc_1dim_contig((void **)(&pdata), par->tnpt, sizeof(float));
  alloc_2dim_contig((void ***)(&can_weight), MAX_STEP*MAX_STEP, par->tnpt, sizeof(float));
  alloc_3dim_contig((void ****)(&oneline), par->tnpt, NUM_LAYERS, par->ncols, sizeof(uint8));

  printf("\n\tRetrieving IGBP & biome land cover");
  /* load IGBP land cover type for this tile */
  /* change input from 1km MOD12Q1 to 500m MCD12Q1 */
  for(i=0; i<1200; i++) {

    for(k=0; k<2; k++) {
      start[0] = i*2+k;
      start[1] = 0;
      stride[0] = 1;
      stride[1] = 1;
      length[0] = 1;
      length[1] = 2400;   
      if((GDreadfield(grid->lc_GDid, "Land_Cover_Type_1", start, stride, length, igbp_rows[k])) == FAILURE) {
	ERROR("read MODIS IGBP land cover error", "fillGaps");
	return FAILURE;
      }
      if((GDreadfield(grid->biome_GDid, "Land_Cover_Type_3", start, stride, length, biome_rows[k])) == FAILURE) {
	ERROR("read MODIS Biome land cover error", "fillGaps");
	return FAILURE;
      }
    }

    for(j=0; j<1200; j++) {

      /* process IGBP class */
      for(k=0; k<17; k++) freq[k] = 0;
      for(m=0; m<2; m++)			
	for(n=0; n<2; n++) {
	  k = igbp_rows[m][j*2+n];
	  if(k>=0&&k<17) freq[k]++;
	}
      
      /* convert igbp from 500m to 1km using majority class */ 
      max_freq = 0;
      lc = 255;
      for(k=0; k<17; k++) {
	if(freq[k] >= max_freq) {
	  max_freq = freq[k];
	  lc = k;
	}
      }
      igbp[i][j] = lc;   
      /* fwrite(&lc, 1, 1, out);
	 fwrite(&max_freq, 1, 1, out);*/

      /* process biome class */
      for(k=0; k<17; k++) freq[k] = 0;
      for(m=0; m<2; m++)			
	for(n=0; n<2; n++) {
	  k = biome_rows[m][j*2+n];
	  if(k>=0&&k<17) freq[k]++;
	}

      /* convert biome from 500m to 1km using majority class */ 
      max_freq = 0;
      lc = 255;
      for(k=0; k<17; k++) {
	if(freq[k] >= max_freq) {
	  max_freq = freq[k];
	  lc = k;
	}
      }
      biome[i][j] = lc;    /* convert 500m lc to 1km using majority class */       
    }    

  }
  /* fclose(out);*/

  printf("\n\tChecking smooth data and building up temporal data foreach IGBP type");
  /* TIMESAT validation check and 
     prepare TIMESAT values (best one from entire tile ) for backup */
  for(i=0; i<17; i++) {
    hq_num_lc[i] = 0;
    tot_ave[i] = 0;
    for(k=0; k<par->tnpt; k++) {
      data_lc[i][k] = TIMESAT_LAI_FV;
      ave_data_lc[i][k] = 0.0;
    }
  }

  for(i=par->rowstart; i<=par->rowstop; i++)
    for(j=par->colstart; j<=par->colstop; j++) {
      ri = i - par->rowstart;
      rj = j - par->colstart;
      /* get data from TIMESAT */
      fread(data, sizeof(int8), par->tnpt, grid->slai);
      
      /* get weight */
      fread(weight, sizeof(float), par->tnpt, grid->qc);

      for(k=0; k<par->tnpt; k++) { 
	if(data[k]<laiRange[0] || data[k]>laiRange[1]) {
	  smooth_flag[ri][rj] = 0;
	  break;
	}
	else
	  smooth_flag[ri][rj] = 1;
      }
      if(smooth_flag[ri][rj] == 0) continue;
      igbp_lc = igbp[i-1][j-1];
      if(igbp_lc >= 17) continue;

      hq_num = 0;
      for(n=0; n<par->tnpt; n++) {
	/* count the number of high quality LAI for this curve */
	if(weight[n] > WEIGHT_TS_HQ) 
	  hq_num++;
      }
      /* pick the one with best quality */
      if(hq_num > hq_num_lc[igbp_lc] && smooth_flag[ri][rj] == 1) {
	hq_num_lc[igbp_lc] = hq_num;
	for(k=0; k<par->tnpt; k++)
	  data_lc[igbp_lc][k] = data[k];
      }
      /* do average of all high quality TIMESAT results (from 50% HQ MODIS data) */
      if(hq_num > 0.50 * par->tnpt) {
	for(k=0; k<par->tnpt; k++)
	  ave_data_lc[igbp_lc][k] += data[k];
	tot_ave[igbp_lc]++;
      }
    }

  /* if there are more 2 HQ TIMESAT curve then use average value 
     otherwise use the best quality one */
  for(k=0; k<par->tnpt; k++)
    for(i=0; i<17; i++) 
      if(tot_ave[i] >= 2)
	data_lc[i][k] = (int8) (ave_data_lc[i][k]/tot_ave[i] + 0.5);

#ifdef DEBUG 
  printf("\nGeneric Temporal Curve for Each IGBP Class");
  printf("\nNumIn");
  for(i=0; i<17; i++)
    printf("%4d ", i);
  printf("\nPerHQ");
  for(i=0; i<17; i++)
    printf("%4.1f ",hq_num_lc[i]*100.0/par->tnpt);
  for(k=0; k<par->tnpt; k++) {
    printf("\n%3d: ", k);
    for(i=0; i<17; i++) {
      printf("%4d ", data_lc[i][k]);
    }
    /*printf("\n     ");
    for(i=0; i<17; i++)
    printf("%4.1f ", ave_data_lc[i][k]/tot_ave[i]);*/
  }
#endif

  printf("\n\tProcessing row ...");
  /* start processing for each pixel */
  for(i=par->rowstart-1; i<par->rowstop; i++) {
    printf("%4d\b\b\b\b", i+1);

    for(j=par->colstart-1; j<par->colstop; j++) {
      
      ri = (i-par->rowstart+1);
      rj = (j-par->colstart+1);

      /* get smoothed LAI data */
      offset = (off64_t)(ri*par->ncols+rj) * par->tnpt * sizeof(int8);
      fseeko64(grid->slai, offset, 0);
      fread(data, sizeof(int8), par->tnpt, grid->slai); 
      /* get MODIS LAI for current pixel */
      fseeko64(grid->lai, offset, 0);
      fread(lai_char, sizeof(int8), par->tnpt, grid->lai); 
      /* get MODIS FPAR for current pixel */
      fseeko64(grid->fpar, offset, 0);
      fread(fpar_char, sizeof(int8), par->tnpt, grid->fpar);  
      /* get QC weight for current pixel */
      offset = (off64_t)(ri*par->ncols+rj) * par->tnpt * sizeof(float);
      fseeko64(grid->qc, offset, 0);
      fread(cweight, sizeof(float), par->tnpt, grid->qc);

      /* decide QA for original LAI product */
      for(k=par->map_start; k<=par->map_end; k++) {
	/* all MODIS HQ data in 1 sigma range of first fitting */
	/* w = 1.0 / (1 + sigma / (2*sigma)) = 2/3 */
	if(cweight[k] > WEIGHT_TS_HQ)
	  org_qc[k] = ORG_HQ_TS_HQ;
	else if(cweight[k] > WEIGHT_TS_LQ)
	  org_qc[k] = ORG_HQ_TS_LQ;						    
	else if(cweight[k] > WEIGHT_BACKUP)
	  org_qc[k] = ORG_BACKUP;
	else
	  /* use original LAI fill values */
	  /*org_qc[k] = cdata[k];*/
	  org_qc[k] = ORG_FILL;
      }

      /* check number of valid LAI values of current pixel */
      num_valid = 0;
      for(k=0; k<par->tnpt; k++) {
	cdata[k] = lai_char[k];
	if(cweight[k]>WEIGHT_VALID && cdata[k]>laiRange[0] && cdata[k]<laiRange[1]) {
	  num_valid++;
	}
      }

      /* use timesat data if it's valid or IGBP class is water */
      valid_flag = SUCCESS;
      for(k=par->map_start; k<=par->map_end; k++) {

	if(data[k] >= laiRange[1] + 10 || data[k] <= laiRange[0] - 10) {
	  /* exclude water, urban, snow and ice, unclassified and fill value */
	  /* ALSO can't do anything if num_valid = 0 for current pixel */
	  if(igbp[i][j] == 0 || igbp[i][j] == 13 || igbp[i][j] ==15 || igbp[i][j] == 17 
	     || igbp[i][j] == 254 || igbp[i][j] == 255 || num_valid == 0) { 
	    sdata[k] = grid->modis.fillValue;
	    fit_qc[k] = FIT_FILL;
	  }
	  else {
	    valid_flag = FAILURE;
	    break;
	  }
	}
	else if(data[k] >= laiRange[1]) {
	  sdata[k] = laiRange[1] - 1;
	  fit_qc[k] = FIT_ROUNDED;
	}
	else if(data[k] <= laiRange[0]) {
	  sdata[k] = laiRange[0] + 1;
	  fit_qc[k] = FIT_ROUNDED;
	}
	else {
	  sdata[k] = data[k];
	  fit_qc[k] = FIT_TIMESAT;
	}	  
      }

#ifdef DEBUG
      if(i==DEBUG_irow && j==DEBUG_icol) {
	printf("\nOriginal MODIS LAI and TIMESAT LAI for (%d, %d) IGBP=%d\n", j, i, igbp[i][j]);
	for(k=par->map_start; k<=par->map_end; k++) 
	  printf("Modis=%4.1f  Weight=%4.1f  Timesat=%4d  Adj_timesat=%4d\n", 
		 cdata[k], cweight[k], data[k], sdata[k]); 
      }
#endif	
	  
      if(valid_flag == FAILURE)  {
	step = 5;  /* first search around +-5 pixels */
	num = 0;   /* total number of replacement candidates */
	do {
	  /* search small window and find pixels with same IGBP class */
	  for(m=i-step; m<=i+step; m++)
	    for(n=j-step; n<=j+step; n++) {

	      /* search must within image range */
	      if(m>=0 && m<par->nrows && n>=0 && n<par->ncols)
		/* process neighbor pixel with same IGBP class */
		if(!(m==i && n==j) && igbp[m][n]==igbp[i][j]) {
		  /*#ifdef DEBUG
		  printf("(%d, %d) = (%d, %d) lc=%d\n", m, n, i, j, igbp[m][n]); 	  
		  #endif*/
		  
		  ri = (m-par->rowstart+1);
		  rj = (n-par->colstart+1);
		  
		  /* check if neighbor pixel is valid */		  
		  if(smooth_flag[ri][rj] == 1) {

		    /* get timesat fit LAI of neighbor pixel */
		    offset = (off64_t)(ri*par->ncols+rj) * par->tnpt * sizeof(int8);
		    fseeko64(grid->slai, offset, 0);
		    fread(data, sizeof(int8), par->tnpt, grid->slai);

		    /* get weight for neigbor pixel */  
		    offset = (off64_t)(ri*par->ncols+rj) * par->tnpt * sizeof(float);
		    fseeko64(grid->qc, offset, 0);
		    fread(weight, sizeof(float), par->tnpt, grid->qc);

		    /* save as a candiadte */
		    cnum = 0;
		    diff = 0.0;
		    for(k=0; k<par->tnpt; k++) {
		      can_data[num][k] = data[k];
		      can_weight[num][k] = weight[k];
		      /* see if this curve is close the current obs 
			 MODIS IGBP classes have been filtered and 
			 may not represent the same class as it shows */
		      if(weight[k] >= WEIGHT_TS_HQ && cweight[k] >= WEIGHT_TS_HQ) {		
			diff += fabs(data[k]-cdata[k])/(data[k]+cdata[k])*100.0*2;
			cnum++;
		      }
		    } /* endfor - save to stack */
		    ave_diff = diff/cnum;
		    /* only accept curve which has close value to current pixel */
		    if(ave_diff < ACCEPTABLE_DIFF) 
		      num++; 
		  } /* endif - candidate is valid */
		  
		} /* endif - same IGBP class */		
	    } /* endfor - window search */
	  step *= 1.5;  /* increase searching area */
	} while(num==0 && step<=MAX_STEP); /* loop and increase window size if no class is found */
	
	if(num != 0) {
	  /* select best fitting data from candiadtes with highest quality */
	  max_hq_num = 0;
	  for(m=0; m<num; m++) {
	    hq_num = 0;
	    for(n=0; n<par->tnpt; n++) {
	      /* count the number of high quality LAI for this curve */
	      if(can_weight[m][n] >= WEIGHT_TS_HQ) 
		hq_num++;
	    }

	    /*#ifdef DEBUG
	    printf("%d %d %d\n", m, hq_num, max_hq_num);
	    #endif*/
	    if(hq_num>max_hq_num) {
	      max_hq_num = hq_num;
	      best = m;
	    }
	  } /* endfor - select highest quality candidate */

	  for(k=0; k<par->tnpt; k++) 
	    ndata[k] = can_data[best][k];
	}
	else {
	  /* use backup one from entire tile */
	  best = igbp[i][j];
	  for(k=0; k<par->tnpt; k++) 
	    ndata[k] = data_lc[best][k];
	}

	/* first try local fitting using smoothed LAI from neighbor pixel */
	if(computeFillfromLocal(par, cdata, cweight, ndata, pdata) != FAILURE) {
	  for(k=par->map_start; k<=par->map_end; k++) {
	    
	    temp = (int8)(pdata[k]+0.5);
	    sdata[k] = (uint8)(pdata[k]+0.5);

	    if(sdata[k] != LAI_INVALID) {
	      fit_qc[k] = FIT_GAPFILLED;
		
	      if(temp >= laiRange[1]) {
		sdata[k] = laiRange[1] - 1;
		fit_qc[k] = FIT_ROUNDED;
	      }
	      if(temp <= laiRange[0]) {
		sdata[k] = laiRange[0] + 1;
		fit_qc[k] = FIT_ROUNDED;
	      }
	    }
	    else {
	      fit_qc[k] = FIT_FILL;
	      sdata[k] = grid->modis.fillValue;
	    }
	  }

#ifdef DEBUG
	  if(i==DEBUG_irow && j==DEBUG_icol) {
	    /* if(igbp[i][j] == 4) {*/
	    printf("\ncompute from local fitting for (%d, %d) num_can=%d IGBP=%d", j, i, num, igbp[i][j]);
	    printf("\nMODIS, Weight, Neighbor, Computed, Saved\n");
	    for(k=par->map_start; k<=par->map_end; k++) 
	      if(cweight[k] >= WEIGHT_TS_HQ) 
		printf("%3d, %5.1f, %5.2f, %5.1f, %5.1f, %3d\n", k, 
		       cdata[k], cweight[k], ndata[k], pdata[k], sdata[k]);
	      else
		printf("%3d,      , %5.2f, %5.1f, %5.1f, %3d\n", k, 
		       cweight[k], ndata[k], pdata[k], sdata[k]);
	    getchar();
	  }
#endif
	  
	}
	else {
	  /* if local try fails then use global fitting */
	  if(computeFillfromAll(par, cdata, cweight, ndata, pdata) != FAILURE) {
	    for(k=par->map_start; k<=par->map_end; k++) {
	      
	      temp = (int8)(pdata[k]+0.5); 
	      sdata[k] = (uint8) (pdata[k]+0.5);
	      
	      if(sdata[k] != grid->modis.fillValue) {
		
		fit_qc[k] = FIT_GAPFILLED;
		
		if(temp >= laiRange[1]) {
		  sdata[k] = laiRange[1] - 1;
		  fit_qc[k] = FIT_ROUNDED;
		}
		
		if(temp <= laiRange[0]) {
		  sdata[k] = laiRange[0] + 1;
		  fit_qc[k] = FIT_ROUNDED;
		}
	      }
	      else
		fit_qc[k] = FIT_FILL;
	    }
#ifdef DEBUG
	    if(i==DEBUG_irow && j==DEBUG_icol) {
	      printf("\ncompute from global fitting for (%d, %d) IGBP=%d", j, i, igbp[i][j]);
	      printf("\nMODIS, Weight, Neighbor, Computed, Saved\n");
	      for(k=par->map_start; k<=par->map_end; k++) 
		if(cweight[k] >= WEIGHT_TS_HQ) 
		  printf("%3d, %5.1f, %5.2f, %5.1f, %5.1f, %3d\n", k, 
			 cdata[k], cweight[k], ndata[k], pdata[k], sdata[k]);
		else
		  printf("%3d,      , %5.2f, %5.1f, %5.1f, %3d\n", k, 
			 cweight[k], ndata[k], pdata[k], sdata[k]);
	    }
#endif
	    
	  }
	  else { 
	    /* if global fitting fails, then set fill value */
	    for(k=par->map_start; k<=par->map_end; k++) {
	      sdata[k] = grid->modis.fillValue;
	      fit_qc[k] = FIT_FILL;
	    }
	  }
	} 
	
      } /* endif - else do gap filling */

      /* do look up table conversion from LAI to FPAR */
      for(k=par->map_start; k<=par->map_end; k++) {
	if(fit_qc[k] != FIT_FILL) {
	  lai_smooth[k] = sdata[k];
	  fpar_lut[k] = lai2fpar(lai_smooth[k], biome[i][j], par);
	}
	else {
	  fpar_lut[k] = grid->modis_fpar.fillValue;
	  lai_smooth[k] = grid->modis.fillValue;
	}
      }

      /* do composed data -- keep high quality LAI but replace low quality data with smoothed LAI */
      for(k=par->map_start; k<=par->map_end; k++) {
	/* Feng - used MODIS HQ LAI regardless of timesat fit (9/2012) */ 
	if(org_qc[k] == ORG_HQ_TS_HQ || org_qc[k] == ORG_HQ_TS_LQ ) {
	  lai_composed[k] = cdata[k];
	  fpar_composed[k] = fpar_char[k];
	  com_qc[k] = COM_ORG_HQ;
	} else if (fit_qc[k] != FIT_FILL) {
	  lai_composed[k] = sdata[k];  
	  fpar_composed[k] = fpar_lut[k];
	  com_qc[k] = COM_FIT;
	} else {
	  lai_composed[k] = cdata[k];  
	  fpar_composed[k] = fpar_char[k];
	  com_qc[k] = COM_FILL;
	}
      }
 
      /* write composited LAI and QA to binary file in BIP format */
      for(k=par->map_start; k<=par->map_end; k++) {
	bip[0] = lai_char[k];
	bip[1] = lai_smooth[k];
	bip[2] = lai_composed[k];
	bip[3] = fpar_char[k];
	bip[4] = fpar_lut[k];
	bip[5] = fpar_composed[k];
	bip[6] = org_qc[k];
	bip[7] = fit_qc[k];
	bip[8] = com_qc[k];
	fwrite(bip, sizeof(uint8), NUM_LAYERS, grid->composed_lai[k]);
	for(m=0; m<NUM_LAYERS; m++) 
	  oneline[k][m][j-par->colstart+1] = bip[m];
      }

#ifdef DEBUG
      if(i==DEBUG_irow && j==DEBUG_icol) {
	printf("\nMODIS, TIMESAT & Gap-filled, and Composed LAI for (%d, %d)\n", j, i);
	for(k=par->map_start; k<=par->map_end; k++)
	  printf("M_LAI=%3d T_LAI=%3d C_LAI=%3d M_FPAR=%3d LUT_FPAR=%3d C_FPAR=%3d Modis_QC=%3d  T&G_QC=%1d  Com_QC=%1d\n", 
		 lai_char[k], lai_smooth[k], lai_composed[k], fpar_char[k], fpar_lut[k], fpar_composed[k], 
		 org_qc[k], fit_qc[k], com_qc[k]);
      }
#endif	
      
    } /* end of icol */        

    /* write gap-filled smoothed LAI to HDF files */
    for(k=par->map_start; k<=par->map_end; k++) {
      start[0] = i;
      start[1] = par->colstart-1;
      length[0] = 1;
      length[1] = par->ncols;
      for(m=0; m<NUM_LAYERS; m++) {
	for(n=0; n<par->ncols; n++)
	  buffer[n] = oneline[k][m][n];
	if((GDwritefield(grid->out_GDid[k], grid->allNames[m], start, stride, length, buffer)) == FAILURE) {
	  ERROR("write MODIS_NACP LAI or QC error", "fillGaps");
	  return FAILURE;
	}
      }
    } /* end of file */

  } /* end of irow */

  free(sdata);
  free(buffer);
  free(lai_char);
  free(lai_smooth);
  free(lai_composed);
  free(fpar_char);
  free(fpar_lut);
  free(fpar_composed);
  free(org_qc);
  free(fit_qc);
  free(com_qc);
  free(data);
  free(cdata);
  free(weight);
  free(cweight);
  free(ndata);
  free(pdata);

  free_2dim_contig((void **)can_data);
  free_2dim_contig((void **)igbp_rows);
  free_2dim_contig((void **)biome_rows);
  free_2dim_contig((void **)igbp);
  free_2dim_contig((void **)biome);
  free_2dim_contig((void **)smooth_flag);
  free_2dim_contig((void **)data_lc);
  free_2dim_contig((void **)ave_data_lc);
  free_2dim_contig((void **)can_weight);

  free_3dim_contig((void ***)oneline);

  return SUCCESS;
}


int lai2fpar(int lai, int lc, PARAMETERS *par)
{
  if(lc==0 || lc>6) lc = 9;
  if( par->lut[lc][lai] == LAI_INVALID ) {
    ERROR("convert lai to fpar error", "lai2fpar");
    return FAILURE;
  }
  else
    return par->lut[lc][lai];
}


int computeFillfromLocal(PARAMETERS *par, float cdata[], float weight[], float ndata[], float pdata[])
{
  int i, j, k, sp, ep, num;
  /*double sum_x, sum_y, sum_x2, sum_xy;
    double a, b;*/
  float cp[MAX_OUTS][3], a[3], tmp, worst;

  if(ndata[0] == TIMESAT_LAI_FV) return FAILURE;

  for(j=par->map_start; j<=par->map_end; j++) {
    
    sp = j - par->nptperyear/2;
    ep = j + par->nptperyear/2;
    if(sp<0) {
      sp = 0;
      ep = par->nptperyear;
    }
    if(ep>par->tnpt) {
      sp = par->tnpt - par->nptperyear; 
      ep = par->tnpt;
    }

    /* remove one maximum outlier */
    worst = 0.0;
    for(i=sp; i<ep; i++) {
      if(weight[i]>=WEIGHT_TS_HQ && cdata[i]>=laiRange[0] && cdata[i]<=laiRange[1]) {
	tmp = fabs(ndata[i]-cdata[i])/ndata[i];
	if(tmp>worst) {
	  worst = tmp;
	  k = i;
	}
      }
    }
#ifdef DEBUG
    /*    printf("worst=%d\n", k);*/
#endif
    /* use least square fit approach to build relationship y=a+b*x+c*x*x */
    num = 0;
    for(i=sp; i<ep; i++) {
      if(weight[i]>=WEIGHT_TS_HQ && cdata[i]>laiRange[0] && cdata[i]<laiRange[1] && i!=k) {
	cp[num][0] = ndata[i];
	cp[num][1] = ndata[i]*ndata[i];
	cp[num][2] = cdata[i];
	num++;
      }
    } 
    if(num >=5 ) {
      getAffine(cp, num, 1, a);
      pdata[j] = a[0] + a[1]*ndata[j] + a[2]*ndata[j]*ndata[j];
    }
    else
      pdata[j] = LAI_INVALID;

    /* use least square fit approach to build relationship y=a*x+b */
    /*num = 0;
    sum_x = 0.0;
    sum_y =0.0;
    sum_x2 = 0.0;
    sum_xy = 0.0;
    for(i=sp; i<ep; i++) {
      if(weight[i]>=0.5 && cdata[i]>=laiRange[0] && cdata[i]<=laiRange[1]) {
	sum_x += ndata[i];
	sum_y += cdata[i];
	sum_x2 += ndata[i]*ndata[i];
	sum_xy += ndata[i]*cdata[i];
	num++;
      }
    } 

    if(num >= 5) {
      a = (num*sum_xy - sum_x*sum_y)/(num*sum_x2 - sum_x*sum_x);
      b = (sum_y - a*sum_x)/num;
      pdata[j] =  a*ndata[j] + b;*/
      /*if(pdata[j] < laiRange[0]) pdata[j] = laiRange[0];
	if(pdata[j] > laiRange[1]) pdata[j] = laiRange[1];
	printf("%d %d %d %f %f\n", j, sp, ep, a, b);*/
    /*}
      else
      pdata[j] = grid->modis.fillValue;*/
  }

  return SUCCESS;
}


int computeFillfromAll(PARAMETERS *par, float cdata[], float weight[], float ndata[], float pdata[])
{
  int i, num;
  /*double sum_x, sum_y, sum_x2, sum_xy;
    double a, b;*/
  float cp[MAX_OUTS][3], a[3];

  if(ndata[0] == TIMESAT_LAI_FV) return FAILURE;

  /* use least square fit approach to build relationship y=a+b*x+c*x*x */
  num =0;
  for(i=0; i<par->tnpt; i++) {
    if(weight[i]>=WEIGHT_TS_HQ && cdata[i]>=laiRange[0] && cdata[i]<=laiRange[1]) {
      cp[num][0] = ndata[i];
      cp[num][1] = ndata[i]*ndata[i];
      cp[num][2] = cdata[i];
      num++;
    }
  }
  
  if(num >= 5) {
    getAffine(cp, num, 1, a);
    for(i=0; i<par->tnpt; i++) 
      pdata[i] = a[0] + a[1]*ndata[i] + a[2]*ndata[i]*ndata[i];
    return SUCCESS;
  }
  else
    return FAILURE;

  /* use least square fit approach to build relationship */
  /*
  num = 0;
  sum_x = 0.0;
  sum_y =0.0;
  sum_x2 = 0.0;
  sum_xy = 0.0;
    
  for(i=0; i<par->tnpt; i++) {
    if(weight[i]>=0.5 && cdata[i]>=laiRange[0] && cdata[i]<=laiRange[1]) {
      sum_x += ndata[i];
      sum_y += cdata[i];
      sum_x2 += ndata[i]*ndata[i];
      sum_xy += ndata[i]*cdata[i];
      num++;
    }
  } 
  
  if(num >= 5) {
    a = (num*sum_xy - sum_x*sum_y)/(num*sum_x2 - sum_x*sum_x);
    b = (sum_y - a*sum_x)/num;
    for(i=0; i<par->tnpt; i++) {
      pdata[i] = a*ndata[i] + b;
    //  if(pdata[i] < laiRange[0]) pdata[i] = laiRange[0];
    //  if(pdata[i] > laiRange[1]) pdata[i] = laiRange[1];
    }
    return SUCCESS;
  }
  else
    return FAILURE; */
}


/* convert timsat fitting output file to standard BIP format */
int convertFittingFile(PARAMETERS *par) 
{
  int i, j;
  int irow, icol;
  float data[MAX_OUTS];
  int8 lai[MAX_OUTS];
  off64_t  offset;

  FILE *in, *out;  
  
  if((in=fopen(par->fit_fn, "rb"))==NULL) {
    printf("Can't open file %s\n", par->fit_fn);
    return FAILURE;
  }

  if((out=fopen(par->temp_fn, "wb"))==NULL) {
    printf("Can't open file for write\n");
    return FAILURE; 
  }

  /* re-organize sensor data for each operation */
  fread(&(par->nyears), sizeof(int), 1, in);
  fread(&(par->nptperyear), sizeof(int), 1, in);
  fread(&(par->rowstart), sizeof(int), 1, in);
  fread(&(par->rowstop), sizeof(int), 1, in);
  fread(&(par->colstart), sizeof(int), 1, in);
  fread(&(par->colstop), sizeof(int), 1, in);
  par->tnpt = par->nyears * par->nptperyear;
  
#ifdef DEBUG 
  printf("nyears=%d nptperyear=%d\n", par->nyears, par->nptperyear);
  printf("rowstart=%d rowstop=%d colstart=%d colstop=%d\n", 
	 par->rowstart, par->rowstop, par->colstart, par->colstop);
#endif

  par->ncols = (par->colstop - par->colstart + 1);
  par->nrows = (par->rowstop - par->rowstart + 1);

  for(i=0; i<par->tnpt; i++)
    lai[i] = TIMESAT_LAI_FV;
  for(i=par->rowstart; i<=par->rowstop; i++)
    for(j=par->colstart; j<=par->colstop; j++)
      fwrite(lai, sizeof(int8), par->tnpt, out);
  
  do {
    fread(&irow, sizeof(int), 1, in);
    fread(&icol, sizeof(int), 1, in);
    fread(data, sizeof(float), par->tnpt, in);
    /*printf("\n%d %d\n", irow, icol);
    for(i=0; i<tnpt; i++)
      printf("%3.0f ", data[i]);
      getchar();*/
    for(i=0; i<par->tnpt; i++)
      lai[i] = (int8) data[i];

    offset = (off64_t)
      ((irow - par->rowstart) * par->ncols + icol - par->colstart) 
      * par->tnpt * sizeof(int8);
    fseeko64(out, offset, 0);
    fwrite(lai, sizeof(int8), par->tnpt, out);
    
  } while(!feof(in));


  fclose(in);
  fclose(out);
  return SUCCESS;
  
}


/* open files for input and output */
int openFiles(PARAMETERS *par, GRID *grid)
{
  int   i, j, index, nfiles, ret, meta_flag;
  int   lc, lai, fpar;
  char  name[Max_StrLen];
  char  fname[MAX_OUTS][Max_StrLen];
  char  oname[Max_StrLen];
  char  cname[Max_StrLen];
  char  tmpname[Max_StrLen];
  char  str[Max_StrLen];
  char  all_meta[53][Max_StrLen];
  char  *p;
  float f1, f2;
  int32 SD_ID, att_id, GDfid, GDid, tsat_id;
  time_t tp;
  struct tm *tm;
  FILE *fp, *hdr;

  /* open observation file */
  if((grid->lai=fopen(par->lai_fn, "rb"))==NULL) {
    ERROR("open data file error", "openFiles");
    return FAILURE;
  }

  if((grid->fpar=fopen(par->fpar_fn, "rb"))==NULL) {
    ERROR("open data file error", "openFiles");
    return FAILURE;
  }

  /* open weighting file */
  if((grid->qc=fopen(par->weight_fn, "rb"))==NULL) {
    ERROR("open data QC/weight file error", "openFiles");
    return FAILURE;
  }

  /* open timesat fitting data */
  if((grid->slai=fopen(par->temp_fn, "rb"))==NULL) {
    ERROR("open temporary timesat fitting file error", "openFiles");
    return FAILURE; 
  }
  
  /* open MODIS land cover file */
  if ((grid->lc_GDfid = GDopen(par->lc_fn, DFACC_READ))<0) {
    ERROR("open MODIS land cover file error", "openFiles");
    grid->lc_GDid = -1;
  }
  else {
    /* open land cover layer */
    if((grid->lc_GDid = GDattach(grid->lc_GDfid, "MOD12Q1")) < 0) {
      grid->lc_GDid = -1;
      ERROR("attach MODIS land cover data error", "openFiles");
      return FAILURE;
    }
  }

  /* open MODIS biome file used for this collection */
  if ((grid->biome_GDfid = GDopen(par->biome_fn, DFACC_READ))<0) {
    ERROR("open MODIS biome file error", "openFiles");
    grid->biome_GDid = -1;
  }
  else {
    /* open land cover layer */
    if((grid->biome_GDid = GDattach(grid->biome_GDfid, "MOD12Q1")) < 0) {
      grid->biome_GDid = -1;
      ERROR("attach MODIS biome data error", "openFiles");
      return FAILURE;
    }
  }
  
  /* get time series file list from timesat input file */
  if((fp=fopen(par->input_fn, "r"))==NULL) {
    ERROR("open timesat input file error", "openFiles");
    return FAILURE;
  }
  fscanf(fp, "%s", par->lai_list);
  fclose(fp);

  if((fp=fopen(par->lai_list, "r"))==NULL) {
    ERROR("open LAI file list error", "openFiles");
    return FAILURE;
  }

  fscanf(fp, "%d", &nfiles);

  meta_flag = 0;
  /* read metadata from first valid file */
  for(i=0; i<nfiles; i++) {
    fscanf(fp, "%s", fname[i]);
    /* metedata only needs retrieve once */
    if(meta_flag == 0) {      
      if(getModisMetaInfo(fname[i], grid) == FAILURE)
	ERROR("get MODIS metadata error, try next file", "openFiles");
      else
	meta_flag = 1;
    }
  }
  

  for(i=0; i<nfiles; i++) {
    
    /* create output file for write */
    if(i>=par->map_start && i<=par->map_end) {
      
      /* get current time */
      time(&tp);
      tm = gmtime(&tp);
    
      /* compose output filename */
      strftime(str, 14, "%Y%j%H%M%S", tm);
  
      /* remove previous ShortName from string and use new prefix */
      strcpy(name, fname[i]);
      p = strtok(name, ".");
      strcpy(oname, ShortName);
      strcpy(cname, ShortName);
      /* compose new LocalGranuleID */
      for(j=0; j<3; j++) {
	p = strtok(NULL, ".");
	strcat(oname,".");
	strcat(oname, p);
	strcat(cname,".");
	strcat(cname, p);
      }
      /*strcat(oname, ".");
	strcat(oname, str);*/

      strcat(oname, ".hdf");
      strcat(cname, ".bip");

      if(writeMetaInfo(oname, grid, i) == FAILURE) {
	ERROR("write output metadata error", "openFiles");
	return FAILURE;
      }

      /* open files for LAI binary output */
      if((grid->composed_lai[i] = fopen(cname, "wb"))==NULL) {
	ERROR("opening composed output file error", "openFiles");
	return FAILURE;
	}	

      /* write ENVI header file for binary output */
      sprintf(tmpname, "%s.hdr", cname);
      if((hdr=fopen(tmpname, "w"))==NULL) {
	ERROR("opening header file for writing error", "openFiles");
	return FAILURE;
      }
      
      fprintf(hdr, "ENVI\n");
      fprintf(hdr, "description = {MODIS_NACP LAI}\n");
      fprintf(hdr, "samples = %d\n", par->ncols);
      fprintf(hdr, "lines   = %d\n", par->nrows);
      fprintf(hdr, "bands   = %d\n", NUM_LAYERS);
      fprintf(hdr, "header offset = 0\n");
      fprintf(hdr, "file type = ENVI Standard\n");
      fprintf(hdr, "data type = 1\n");
      fprintf(hdr, "interleave = bip\n");
      fprintf(hdr, "sensor type = MODIS\n");
      fprintf(hdr, "byte order = 0\n");
      fprintf(hdr, "map info = {Sinusoidal, 1.0000, 1.0000, %lf, %lf,", grid->GD_upleft[0], grid->GD_upleft[1]);
      fprintf(hdr, "9.2739826494e+02, 9.2739826494e+02, , units=Meters}\n");
      fprintf(hdr, "projection info = {16, 6371007.2, 0.000000, 0.0, 0.0, Sinusoidal, units=Meters}\n");
      fprintf(hdr, "band names = {MODIS_LAI, Smoothed_LAI, Composed_LAI, ");
      fprintf(hdr, "MODIS_FPAR, Smoothed_FPAR, Composed_FPAR, ");
      fprintf(hdr, "MODIS Quality, Smoothed Quality, Composed Quality}\n");
      fclose(hdr);

      
    } /* end of creating output files */   
    
  } /* end of opening input files */

  fclose(fp);

  /* load look up table for lai-fpar */
  if(openLUT(par, grid)==FAILURE) {
    ERROR("load lai-fpar LUT error", "openFiles");
    return FAILURE;
  }
  
}


/**
   build lai to fpar LUT per biome type 
   BIOME TYPES
   Land_Cover_Type_3:water = '\0' ;
   Land_Cover_Type_3:grasses_cereal = '\1' ;
   Land_Cover_Type_3:shrubs = '\2' ;
   Land_Cover_Type_3:broadleaf_crops = '\3' ;
   Land_Cover_Type_3:savannah = '\4' ;
   Land_Cover_Type_3:broadleaf_forest = '\5' ;
   Land_Cover_Type_3:needleleaf_forest = '\6' ;
   Land_Cover_Type_3:unvegetated = '\7' ;
   Land_Cover_Type_3:urban = '\8' ;

   int par->lut[10][100] 
   LAI range: 0-10; scale=0.1
   FPAR range: 0-1.0; scale=0.01
   type 0-8 same as above, type 9 represents average case when lc is not valid/available
*/
int openLUT(PARAMETERS *par, GRID *grid) 
{
  char  str[Max_StrLen];
  int   i, j, lc, lai, fpar;
  float f1, f2, a, b, x1, y1, x2, y2;
  FILE  *fp; 

  if((fp=fopen(par->lut_fn, "r"))==NULL) {
    ERROR("open lai-fpar LUT error", "openLUT");
    return FAILURE;
  }
  fscanf(fp, "%s %s %s\n", str,str,str);
  for(i=0; i<10; i++)
    for(j=0; j<100; j++)
      par->lut[i][j] = LAI_INVALID;
  do {
    fscanf(fp, "%d %f %f\n", &lc, &f1, &f2);
    lai = (int)(f1 / grid->modis.scaleFactor + 0.5);
    fpar = (int) (f2 / grid->modis_fpar.scaleFactor + 0.5);
    par->lut[lc][lai] = fpar;
    /*printf("%d %d %d\n", lc, lai, fpar);*/
  } while(!feof(fp));

  fclose(fp);

  /* do interpolation to fill all LUT */
  for(j=0; j<10; j++)
    for(lai=0; lai<100; lai++) {

      if(j==0 || j>=7) lc = 9;  /* use general case */
      else lc = j;

      if(par->lut[lc][lai] ==  LAI_INVALID) {

	/* find two closest valid values */
	x1 = 0;
	for(i=lai; i>=0; i--)
	  if(par->lut[lc][i] != LAI_INVALID) {
	    x1 = i;
	    y1 = par->lut[lc][i];
	    break;
	  }  
	x2 = 99;
	for(i=lai; i<100; i++)
	  if(par->lut[lc][i] != LAI_INVALID) {
	    x2 = i;
	    y2 = par->lut[lc][i];
	    break;
	  }  

	if(x1 == 0) fpar = y2;
	else if(x2 == 99) fpar = y1;
	else {
	  /* do interpolate */
	  a = (float)(y2-y1)/(x2-x1);
	  b = y1-a*x1;	  
	  fpar = (int)(a*lai+b+0.5);
	}
      }
      else
	fpar = par->lut[lc][lai];      

      par->lut[j][lai] = fpar;
    }
  
  /*for(lc=0; lc<10; lc++) {
    for(lai=0; lai<100; lai++)
      printf("%3d %3d %3d\n", lc, lai, par->lut[lc][lai]);
    getchar();
    }*/

}


int closeFiles(PARAMETERS *par, GRID *grid)
{
  int i, ret;

  fclose(grid->lai);
  fclose(grid->fpar);
  fclose(grid->qc);
  
  for(i=par->map_start; i<=par->map_end; i++) {

    if((GDdetach(grid->out_GDid[i])) == FAILURE) {
      ERROR("Failed to detach grid.", "closeFiles");
      return FAILURE;
    }

    /* close for grid access */
    if((GDclose(grid->out_GDfid[i])) == FAILURE){
      ERROR("GD-file close failed.", "closeFiles");
      return FAILURE;
    }

    fclose(grid->composed_lai[i]);
  }
  
  if(grid->lc_GDid != -1) {
    if(GDdetach(grid->lc_GDid) == FAILURE)  {
      ERROR("GDdetach for land cover error", "closeFiles");
      return FAILURE;
    }
    if((GDclose(grid->lc_GDfid)) == FAILURE){
      ERROR("GD-file close failed.", "closeFiles");
      return FAILURE;
    }
  }

  if(grid->biome_GDid != -1) {
    if(GDdetach(grid->biome_GDid) == FAILURE)  {
      ERROR("GDdetach for MODIS biome error", "closeFiles");
      return FAILURE;
    }
    if((GDclose(grid->biome_GDfid)) == FAILURE){
      ERROR("GD-file close failed.", "closeFiles");
      return FAILURE;
    }
  }
  
}  


/* extract input parameters */
int parseParameters(PARAMETERS *par, int argc, char *argv[]) 
{ 
  int i;

  /* command line should include 13 args */ 
  if(argc != 13) {
    usage(argv[0]);
    return FAILURE;
  }

  /* parse command line */
  for(i=1; i<argc; i++){
    if(strcasecmp(argv[i],"-j")==0)
      strcpy(par->job_name, argv[++i]);
    else
      if(strcasecmp(argv[i], "-f")==0)
	strcpy(par->fit_fn, argv[++i]);
      else
	if(strcasecmp(argv[i],"-si")==0)
	  par->map_start = atoi(argv[++i]);
	else
	  if(strcasecmp(argv[i], "-ei")==0)
	    par->map_end = atoi(argv[++i]);
	  else
	    if(strcasecmp(argv[i], "-l")==0)
	      strcpy(par->lc_fn, argv[++i]);
	    else
	      if(strcasecmp(argv[i], "-b")==0)
		strcpy(par->biome_fn, argv[++i]);    
	      else{
		printf("\nWrong option:%s\n",argv[i]);
		usage(argv[0]);
		return FAILURE;
	      }    
  }

  strcpy(par->temp_fn, FitFileInBIP);
  sprintf(par->lai_fn, "originalLAI_%s", par->job_name);
  sprintf(par->fpar_fn, "originalFPAR_%s", par->job_name);
  sprintf(par->weight_fn, "assignedweight_%s", par->job_name);
  sprintf(par->input_fn, "lai_input.%s.txt", par->job_name);
  strcpy(par->lut_fn, "lai_fpar.lut");

  return SUCCESS;
}  


/* display usage */
void usage(char *command)
{
  printf("\nUsage: %s [-j][-f][-si][-ei][-l][-b] \n\n", command);
  printf("   -j  <timesat_job_name>  job name used in TIMESAT\n");
  printf("   -f  <fitting_file>      selected TIMESAT fitted function file\n");
  printf("   -si <map_index>         start image number to create (start from 0)\n");
  printf("   -ei <map_index>         end image number to create (start from 0)\n");
  printf("   -l  <MODIS_LC_file>     MODIS land cover file for production period\n");
  printf("   -b  <MODIS_BIOME_file>  MODIS biome type file used for this collection\n");
}


/* allocate memory, get metadata and open specific sds */
int getModisMetaInfo(char *fname, GRID *g) {

  char name[100], msg[100];
  int i, index, ret;
  char GD_gridlist[100];
  int32 gfid=0, ngrid=0, gid=0;
  int32 SD_ID, lai_id, fpar_id; 
  int32 bufsize=100;
  
  int32 rank,data_type,attributes,dim_sizes[2];
  int32 att_id;

  /* open a hdf file */
  gfid = GDopen(fname, DFACC_READ);
  if(gfid==FAILURE){
    sprintf(msg, "Not successful in retrieving MODIS grid file ID/open from %f", fname);
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }
  
  /* find out about grid type */
  ngrid=GDinqgrid(fname, GD_gridlist, &bufsize);
  if(ngrid==FAILURE){
    sprintf(msg, "Not successful in retrieving MODIS grid name list from %s", fname);
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }

  /* attach grid */
  gid = GDattach(gfid, GD_gridlist);
  if(gid==FAILURE){
    sprintf(msg, "Not successful in attaching MODIS grid from %s", fname);
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }

  /* get grid info */
  ret = GDgridinfo(gid, &(g->ncols), &(g->nrows), g->GD_upleft, g->GD_lowright);
  if(ret==FAILURE){
      sprintf(msg, "Failed to read grid info from %s", fname);
      ERROR(msg, "getModisMetaInfo");
      return FAILURE;
  }

  /* get projection parameters */
  ret = GDprojinfo(gid, &g->GD_projcode, &g->GD_zonecode, &g->GD_spherecode, g->GD_projparm);
  if(ret==FAILURE){
    sprintf(msg, "Not successful in reading grid projection info. from %s", fname);
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }

  /* get grid origin */
  ret = GDorigininfo(gid, &g->GD_origincode);
  if(ret==FAILURE){
    sprintf(msg, "Failed to read grid origin info from %s", fname);
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }

  /* detach grid */
  ret = GDdetach(gid);
  if(ret==FAILURE){
      sprintf(msg, "Failed to detach grid for %s", fname);
      ERROR(msg, "getModisMetaInfo");
      return FAILURE;
  }

  /* close for grid access */
  ret = GDclose(gfid);
  if(ret==FAILURE){
      sprintf(msg, "GD-file close failed for %s", fname);
      ERROR(msg, "getModisMetaInfo");
      return FAILURE;
  }

  /* open hdf file and get sds_id from given sds_name */  
  if ((SD_ID = SDstart(fname, DFACC_READ))<0) {
    sprintf(msg, "Can't open file %s", fname);
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }

  /* open SDS layer */
  if ((index=SDnametoindex(SD_ID, "Lai_1km"))<0) {
    sprintf(msg, "Not successful in convert SR sdsName Lai_1km to index");
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }
  lai_id = SDselect(SD_ID, index);
  
  /* retrieve fill value */
  if ((att_id = SDfindattr(lai_id, "_FillValue")) == FAILURE) {
    sprintf(msg, "Can't retrieve fill value from SDS attr");
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }
  if (SDreadattr(lai_id, att_id, &(g->modis.fillValue)) == FAILURE) {
    sprintf(msg, "Can't retrieve fill value from SDS attr.");
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }

  /* retrieve scale factor */
  if ((att_id = SDfindattr(lai_id, "scale_factor")) == FAILURE) {
    sprintf(msg, "Can't retrieve scale factor from SDS attr");
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }
  if (SDreadattr(lai_id, att_id, &(g->modis.scaleFactor)) == FAILURE) {
    sprintf(msg, "Can't retrieve scale factor from SDS attr.");
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }

  /* retrieve data range */
  if ((att_id = SDfindattr(lai_id, "valid_range")) == FAILURE) {
    sprintf(msg, "Can't retrieve valid range from SDS attr");
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }
  if (SDreadattr(lai_id, att_id, g->modis.range) == FAILURE) {
    sprintf(msg, "Can't retrieve valid range from SDS attr.");
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }

  /* retrieve fill value description */
  if ((att_id = SDfindattr(lai_id, "MOD15A2_FILLVALUE_DOC")) == FAILURE) {
    sprintf(msg, "Can't retrieve valid range from SDS attr");
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }
  if (SDreadattr(lai_id, att_id, g->modis.explain) == FAILURE) {
    sprintf(msg, "Can't retrieve valid range from SDS attr.");
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }

  if((SDendaccess(lai_id)) == FAILURE) {
    sprintf(msg, "SDendaccess for %s QA error!", fname);
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }

  /* open SDS layer */
  if ((index=SDnametoindex(SD_ID, "Fpar_1km"))<0) {
    sprintf(msg, "Not successful in convert SR sdsName Fpar_1km to index");
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }
  fpar_id = SDselect(SD_ID, index);
  
  /* retrieve fill value */
  if ((att_id = SDfindattr(fpar_id, "_FillValue")) == FAILURE) {
    sprintf(msg, "Can't retrieve fill value from SDS attr");
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }
  if (SDreadattr(fpar_id, att_id, &(g->modis_fpar.fillValue)) == FAILURE) {
    sprintf(msg, "Can't retrieve fill value from SDS attr.");
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }

  /* retrieve scale factor */
  if ((att_id = SDfindattr(fpar_id, "scale_factor")) == FAILURE) {
    sprintf(msg, "Can't retrieve scale factor from SDS attr");
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }
  if (SDreadattr(fpar_id, att_id, &(g->modis_fpar.scaleFactor)) == FAILURE) {
    sprintf(msg, "Can't retrieve scale factor from SDS attr.");
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }

  /* retrieve data range */
  if ((att_id = SDfindattr(fpar_id, "valid_range")) == FAILURE) {
    sprintf(msg, "Can't retrieve valid range from SDS attr");
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }
  if (SDreadattr(fpar_id, att_id, g->modis_fpar.range) == FAILURE) {
    sprintf(msg, "Can't retrieve valid range from SDS attr.");
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }

  /* retrieve fill value description */
  if ((att_id = SDfindattr(fpar_id, "MOD15A2_FILLVALUE_DOC")) == FAILURE) {
    sprintf(msg, "Can't retrieve MOD15A2_FILLVALUE_DOC from SDS attr");
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }
  if (SDreadattr(fpar_id, att_id, g->modis_fpar.explain) == FAILURE) {
    sprintf(msg, "Can't retrieve MOD15A2_FILLVALUE_DOC from SDS attr.");
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }

  if((SDendaccess(fpar_id)) == FAILURE) {
    sprintf(msg, "SDendaccess for %s QA error!", fname);
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }

  if((SDend(SD_ID)) == FAILURE) {
    sprintf(msg, "SDend for %s error!", fname);
    ERROR(msg, "getModisMetaInfo");
    return FAILURE;
  }

  /* set parameters for LAI */
  strcpy(g->modis.name, "MODIS_LAI");
  strcpy(g->modis.longName, "Original MODIS LAI");

  strcpy(g->timesat.name, "Smoothed_LAI");
  strcpy(g->timesat.longName, "Gap-Filled and Timesat-Smoothed LAI");
  g->timesat.fillValue = g->modis.fillValue;
  g->timesat.range[0] = g->modis.range[0];
  g->timesat.range[1] = g->modis.range[1];
  g->timesat.scaleFactor = g->modis.scaleFactor;

  strcpy(g->compose.name, "Composed_LAI");
  strcpy(g->compose.longName, "Composed LAI from MODIS and Timesat");
  g->compose.fillValue = g->modis.fillValue;
  g->compose.range[0] = g->modis.range[0];
  g->compose.range[1] = g->modis.range[1];
  g->compose.scaleFactor = g->modis.scaleFactor;

  /* set parameters for FPAR */
  strcpy(g->modis_fpar.name, "MODIS_FPAR");
  strcpy(g->modis_fpar.longName, "Original MODIS FPAR");

  strcpy(g->lut_fpar.name, "Converted_FPAR");
  strcpy(g->lut_fpar.longName, "Converted FPAR from LAI-FPAR LUT");
  g->lut_fpar.fillValue = g->modis_fpar.fillValue;
  g->lut_fpar.range[0] = g->modis_fpar.range[0];
  g->lut_fpar.range[1] = g->modis_fpar.range[1];
  g->lut_fpar.scaleFactor = g->modis_fpar.scaleFactor;

  strcpy(g->compose_fpar.name, "Composed_FPAR");
  strcpy(g->compose_fpar.longName, "Composed FPAR from MODIS and Converted FPAR");
  g->compose_fpar.fillValue = g->modis_fpar.fillValue;
  g->compose_fpar.range[0] = g->modis_fpar.range[0];
  g->compose_fpar.range[1] = g->modis_fpar.range[1];
  g->compose_fpar.scaleFactor = g->modis_fpar.scaleFactor;

  /* set parameters for QC */
  strcpy(g->modis_qc.name, "MODIS_LAI_FPAR_QC");
  strcpy(g->modis_qc.longName, "MODIS LAI/FPAR Retrieval Quality");
  strcpy(g->timesat_qc.name, "Smoothed_LAI_FPAR_QC");
  strcpy(g->timesat_qc.longName, "Gap-Filled and Timesat-Smoothed LAI/FPAR Quality");
  strcpy(g->compose_qc.name, "Composed_LAI_FPAR_QC");
  strcpy(g->compose_qc.longName, "Composed LAI/FPAR Quality");

  strcpy(g->allNames[0], g->modis.name);
  strcpy(g->allNames[1], g->timesat.name);
  strcpy(g->allNames[2], g->compose.name);
  strcpy(g->allNames[3], g->modis_fpar.name);
  strcpy(g->allNames[4], g->lut_fpar.name);
  strcpy(g->allNames[5], g->compose_fpar.name);
  strcpy(g->allNames[6], g->modis_qc.name);
  strcpy(g->allNames[7], g->timesat_qc.name);
  strcpy(g->allNames[8], g->compose_qc.name);

#ifdef DEBUG
  printf("\nFile:  %s",fname);
  printf("\nnrows: %d,  ncols: %d",g->nrows, g->ncols);
  printf("\nfillV: %d", g->modis.fillValue); 
  printf("\nrange: (%d, %d)", g->modis.range[0], g->modis.range[1]);
  printf("\nscale: %f\n", g->modis.scaleFactor);
  printf("%s\n", g->modis.explain);
#endif

  return SUCCESS;
}


/* write grid metadata to output HDF file */
int writeMetaInfo(char *fname, GRID *g, int i)
{
  char msg[100], str[Max_StrLen];
  int32 SD_ID, GDfid, GDid, lai_id, fpar_id;
  int32 tile_rank;
  int32 tile_dims[2];
  intn comp_param[4];
  int ret, index;
 
  /* set compression parameters */
  tile_rank = 2;
  tile_dims[0] = g->nrows;
  tile_dims[1] = g->ncols;
  comp_param[0] = 8;

  /* create output use HDF-EOS functions */
  GDfid = GDopen(fname, DFACC_CREATE);
  if(GDfid == FAILURE){
    sprintf(msg, "Not successful in creating grid file %s", fname);
    ERROR(msg, "writeMetaInfo");
    return FAILURE;
  }
 
  /* create a new grid in output */ 
  GDid = GDcreate(GDfid, GridName, g->ncols, g->nrows, g->GD_upleft, g->GD_lowright);
  if(GDid == FAILURE) {
    sprintf(msg, "Not successful in creating grid ID/ create for %s", fname);
    ERROR(msg, "writeMetaInfo");
    return FAILURE;
  }

  /* define grid projection */
  ret = GDdefproj(GDid, g->GD_projcode, g->GD_zonecode, g->GD_spherecode, g->GD_projparm);
  if(ret==FAILURE){
    sprintf(msg, "Not successful in defining grid projection");
    ERROR(msg, "writeMetaInfo");
    return FAILURE;
  }

  /* define grid origin */
  ret = GDdeforigin(GDid, g->GD_origincode);
  if(ret==FAILURE){
    sprintf (msg, "Not successful in defining grid origin");
    ERROR(msg, "writeMetaInfo");
    return FAILURE;
  }

  /* define MODIS LAI layer */
  ret = GDdeffield(GDid, g->modis.name, "YDim,XDim", DFNT_UINT8, HDFE_NOMERGE);
  if(ret==FAILURE){
    ERROR("define MODIS LAI sds layer errro", "writeMetaInfo");
    return FAILURE;
  }

  ret=GDsettilecomp(GDid, g->modis.name, tile_rank, tile_dims, HDFE_COMP_DEFLATE, comp_param);
  if(ret==FAILURE)
    {
      sprintf (msg, "Error in put compression info for sds %s", g->modis.name);
      ERROR(msg, "WriteMetaDataInfo");
      return FAILURE;
    }
  
  /* define smoothed LAI layer */
  ret = GDdeffield(GDid, g->timesat.name, "YDim,XDim", DFNT_UINT8, HDFE_NOMERGE);
  if(ret==FAILURE){
    ERROR("define smoothed LAI sds layer errro", "writeMetaInfo");
    return FAILURE;
  }

  ret=GDsettilecomp(GDid, g->timesat.name, tile_rank, tile_dims, HDFE_COMP_DEFLATE, comp_param);
  if(ret==FAILURE)
    {
      sprintf (msg, "Error in put compression info for sds %s", g->timesat.name);
      ERROR(msg, "WriteMetaDataInfo");
      return FAILURE;
    }
  
  /* define composite LAI layer */
  ret = GDdeffield(GDid, g->compose.name, "YDim,XDim", DFNT_UINT8, HDFE_NOMERGE);
  if(ret==FAILURE){
    ERROR("define composed LAI sds layer errro", "writeMetaInfo");
    return FAILURE;
  }

  ret=GDsettilecomp(GDid, g->compose.name, tile_rank, tile_dims, HDFE_COMP_DEFLATE, comp_param);
  if(ret==FAILURE)
    {
      sprintf (msg, "Error in put compression info for sds %s", g->compose.name);
      ERROR(msg, "WriteMetaDataInfo");
      return FAILURE;
    }
  

  /* define FPAR layers */
  ret = GDdeffield(GDid, g->modis_fpar.name, "YDim,XDim", DFNT_UINT8, HDFE_NOMERGE);
  if(ret==FAILURE){
    ERROR("define MODIS LAI sds layer errro", "writeMetaInfo");
    return FAILURE;
  }

  ret=GDsettilecomp(GDid, g->modis_fpar.name, tile_rank, tile_dims, HDFE_COMP_DEFLATE, comp_param);
  if(ret==FAILURE)
    {
      sprintf (msg, "Error in put compression info for sds %s", g->modis_fpar.name);
      ERROR(msg, "WriteMetaDataInfo");
      return FAILURE;
    }
  

  ret = GDdeffield(GDid, g->lut_fpar.name, "YDim,XDim", DFNT_UINT8, HDFE_NOMERGE);
  if(ret==FAILURE){
    ERROR("define smoothed LAI sds layer errro", "writeMetaInfo");
    return FAILURE;
  }

  ret=GDsettilecomp(GDid, g->lut_fpar.name, tile_rank, tile_dims, HDFE_COMP_DEFLATE, comp_param);
  if(ret==FAILURE)
    {
      sprintf (msg, "Error in put compression info for sds %s", g->lut_fpar.name);
      ERROR(msg, "WriteMetaDataInfo");
      return FAILURE;
    }
  
  ret = GDdeffield(GDid, g->compose_fpar.name, "YDim,XDim", DFNT_UINT8, HDFE_NOMERGE);
  if(ret==FAILURE){
    ERROR("define composed LAI sds layer errro", "writeMetaInfo");
    return FAILURE;
  }

  ret=GDsettilecomp(GDid, g->compose_fpar.name, tile_rank, tile_dims, HDFE_COMP_DEFLATE, comp_param);
  if(ret==FAILURE)
    {
      sprintf (msg, "Error in put compression info for sds %s", g->compose_fpar.name);
      ERROR(msg, "WriteMetaDataInfo");
      return FAILURE;
    }
  

  /* define MODIS QC layer */
  ret = GDdeffield(GDid, g->modis_qc.name, "YDim,XDim", DFNT_UINT8, HDFE_NOMERGE);
  if(ret==FAILURE){
    ERROR("define MODIS LAI QC sds layer errro", "writeMetaInfo");
    return FAILURE;
  }

  ret=GDsettilecomp(GDid, g->modis_qc.name, tile_rank, tile_dims, HDFE_COMP_DEFLATE, comp_param);
  if(ret==FAILURE)
    {
      sprintf (msg, "Error in put compression info for sds %s", g->modis_qc.name);
      ERROR(msg, "WriteMetaDataInfo");
      return FAILURE;
    }  

  /* define TIMESAT QC layer */
  ret = GDdeffield(GDid, g->timesat_qc.name, "YDim,XDim", DFNT_UINT8, HDFE_NOMERGE);
  if(ret==FAILURE){
    ERROR("define smoothed LAI QC sds layer errro", "writeMetaInfo");
    return FAILURE;
  }

  ret=GDsettilecomp(GDid, g->timesat_qc.name, tile_rank, tile_dims, HDFE_COMP_DEFLATE, comp_param);
  if(ret==FAILURE)
    {
      sprintf (msg, "Error in put compression info for sds %s", g->timesat_qc.name);
      ERROR(msg, "WriteMetaDataInfo");
      return FAILURE;
    }

  /* define composite QC layer */
  ret = GDdeffield(GDid, g->compose_qc.name, "YDim,XDim", DFNT_UINT8, HDFE_NOMERGE);
  if(ret==FAILURE){
    ERROR("define composed LAI QC sds layer errro", "writeMetaInfo");
    return FAILURE;
  }

  ret=GDsettilecomp(GDid, g->compose_qc.name, tile_rank, tile_dims, HDFE_COMP_DEFLATE, comp_param);
  if(ret==FAILURE)
    {
      sprintf (msg, "Error in put compression info for sds %s", g->compose_qc.name);
      ERROR(msg, "WriteMetaDataInfo");
      return FAILURE;
    } 

  /* detach grid */
  ret = GDdetach(GDid);
  if(ret==FAILURE){
    sprintf (msg, "Failed to detach grid.");
    ERROR(msg, "writeMetaInfo");
    return FAILURE;
  }

  /* close for grid access */
  ret = GDclose(GDfid);
  if(ret==FAILURE){
    sprintf (msg, "GD-file close failed.");
    ERROR(msg, "writeMetaInfo");
    return FAILURE;
  }

  /* reopen use HDF function to write SDS level metadata (don't know if there is related function in HDF-EOS?) */
  /* open hdf file and get sds_id from given sds_name and then write metedata */  
  if ((SD_ID = SDstart(fname, DFACC_RDWR))<0) {
    sprintf(msg, "Can't open file %s", fname);
    ERROR(msg, "writeMetaInfo");
    return FAILURE;
  }
  
  /* write metadata for MODIS LAI */
  if ((index = SDnametoindex(SD_ID, g->modis.name))<0) {
    ERROR("convert MODIS LAI sds layer error", "writeMetaInfo");
    return FAILURE;
  }
  lai_id = SDselect(SD_ID, index);           

  ret = SDsetattr(lai_id, "long_name", DFNT_CHAR8, strlen(g->modis.longName), g->modis.longName);
  if (ret == FAILURE) {
    ERROR("write long name for MODIS LAI error", "writeMetaInfo");
    return FAILURE;
  } 

  ret = SDsetattr(lai_id, "scale_factor", DFNT_FLOAT64, 1, &(g->modis.scaleFactor));
  if (ret == FAILURE) {
    ERROR("write scale factor error", "writeMetaInfo");
    return FAILURE;
  } 

  ret = SDsetattr(lai_id, "valid_range", DFNT_UINT8, 2, g->modis.range);
  if (ret == FAILURE) {
    ERROR("write valid range error", "writeMetaInfo");
    return FAILURE;
  } 
      
  ret = SDsetattr(lai_id, "_FillValue", DFNT_UINT8, 1, &(g->modis.fillValue));
  if (ret == FAILURE) {
    ERROR("write _FillValue error", "writeMetaInfo");
    return FAILURE;
  } 

  ret = SDsetattr(lai_id, "fill_value_legend", DFNT_CHAR8, strlen(g->modis.explain), g->modis.explain);
  if (ret == FAILURE) {
    ERROR("write description for MODIS LAI error", "writeMetaInfo");
    return FAILURE;
  } 

  if((SDendaccess(lai_id)) == FAILURE) {
    ERROR("SDendaccess error", "writeMetaInfo");
    return FAILURE;
  }

  /* write metadata for SMOOTH LAI */
  if ((index = SDnametoindex(SD_ID, g->timesat.name))<0) {
    ERROR("convert smoothed LAI sds layer error", "writeMetaInfo");
    return FAILURE;
  }
  lai_id = SDselect(SD_ID, index);           

  ret = SDsetattr(lai_id, "long_name", DFNT_CHAR8, strlen(g->timesat.longName), g->timesat.longName);
  if (ret == FAILURE) {
    ERROR("write long name for smoothed LAI error", "writeMetaInfo");
    return FAILURE;
  } 

  ret = SDsetattr(lai_id, "scale_factor", DFNT_FLOAT64, 1, &(g->timesat.scaleFactor));
  if (ret == FAILURE) {
    ERROR("write scale factor error", "writeMetaInfo");
    return FAILURE;
  } 

  ret = SDsetattr(lai_id, "valid_range", DFNT_UINT8, 2, g->timesat.range);
  if (ret == FAILURE) {
    ERROR("write valid range error", "writeMetaInfo");
    return FAILURE;
  } 
      
  ret = SDsetattr(lai_id, "_FillValue", DFNT_UINT8, 1, &(g->timesat.fillValue));
  if (ret == FAILURE) {
    ERROR("write _FillValue error", "writeMetaInfo");
    return FAILURE;
  } 
   
  if((SDendaccess(lai_id)) == FAILURE) {
    ERROR("SDendaccess error", "writeMetaInfo");
    return FAILURE;
  }

  /* write metadata for COMPOSE LAI */
  if ((index = SDnametoindex(SD_ID, g->compose.name))<0) {
    ERROR("convert composed LAI sds layer error", "writeMetaInfo");
    return FAILURE;
  }
  lai_id = SDselect(SD_ID, index);           

  ret = SDsetattr(lai_id, "long_name", DFNT_CHAR8, strlen(g->compose.longName), g->compose.longName);
  if (ret == FAILURE) {
    ERROR("write long name for composed LAI error", "writeMetaInfo");
    return FAILURE;
  } 

  ret = SDsetattr(lai_id, "scale_factor", DFNT_FLOAT64, 1, &(g->compose.scaleFactor));
  if (ret == FAILURE) {
    ERROR("write scale factor error", "writeMetaInfo");
    return FAILURE;
  } 

  ret = SDsetattr(lai_id, "valid_range", DFNT_UINT8, 2, g->compose.range);
  if (ret == FAILURE) {
    ERROR("write valid range error", "writeMetaInfo");
    return FAILURE;
  } 
      
  ret = SDsetattr(lai_id, "_FillValue", DFNT_UINT8, 1, &(g->compose.fillValue));
  if (ret == FAILURE) {
    ERROR("write _FillValue error", "writeMetaInfo");
    return FAILURE;
  } 

  ret = SDsetattr(lai_id, "fill_value_legend", DFNT_CHAR8, strlen(g->modis.explain), g->modis.explain);
  if (ret == FAILURE) {
    ERROR("write description for MODIS LAI error", "writeMetaInfo");
    return FAILURE;
  } 

  if((SDendaccess(lai_id)) == FAILURE) {
    ERROR("SDendaccess error", "writeMetaInfo");
    return FAILURE;
  }

  /* write metadata for MODIS FPAR */
  if ((index = SDnametoindex(SD_ID, g->modis_fpar.name))<0) {
    ERROR("convert MODIS FPAR sds layer error", "writeMetaInfo");
    return FAILURE;
  }
  fpar_id = SDselect(SD_ID, index);           

  ret = SDsetattr(fpar_id, "long_name", DFNT_CHAR8, strlen(g->modis_fpar.longName), g->modis_fpar.longName);
  if (ret == FAILURE) {
    ERROR("write long name for MODIS FPAR error", "writeMetaInfo");
    return FAILURE;
  } 

  ret = SDsetattr(fpar_id, "scale_factor", DFNT_FLOAT64, 1, &(g->modis_fpar.scaleFactor));
  if (ret == FAILURE) {
    ERROR("write scale factor error", "writeMetaInfo");
    return FAILURE;
  } 

  ret = SDsetattr(fpar_id, "valid_range", DFNT_UINT8, 2, g->modis_fpar.range);
  if (ret == FAILURE) {
    ERROR("write valid range error", "writeMetaInfo");
    return FAILURE;
  } 
      
  ret = SDsetattr(fpar_id, "_FillValue", DFNT_UINT8, 1, &(g->modis_fpar.fillValue));
  if (ret == FAILURE) {
    ERROR("write _FillValue error", "writeMetaInfo");
    return FAILURE;
  } 

  ret = SDsetattr(fpar_id, "fill_value_legend", DFNT_CHAR8, strlen(g->modis_fpar.explain), g->modis_fpar.explain);
  if (ret == FAILURE) {
    ERROR("write description for MODIS FPAR error", "writeMetaInfo");
    return FAILURE;
  } 

  if((SDendaccess(fpar_id)) == FAILURE) {
    ERROR("SDendaccess error", "writeMetaInfo");
    return FAILURE;
  }

  /* write metadata for LUT FPAR */
  if ((index = SDnametoindex(SD_ID, g->lut_fpar.name))<0) {
    ERROR("convert smoothed FPAR sds layer error", "writeMetaInfo");
    return FAILURE;
  }
  fpar_id = SDselect(SD_ID, index);           

  ret = SDsetattr(fpar_id, "long_name", DFNT_CHAR8, strlen(g->lut_fpar.longName), g->lut_fpar.longName);
  if (ret == FAILURE) {
    ERROR("write long name for smoothed FPAR error", "writeMetaInfo");
    return FAILURE;
  } 

  ret = SDsetattr(fpar_id, "scale_factor", DFNT_FLOAT64, 1, &(g->lut_fpar.scaleFactor));
  if (ret == FAILURE) {
    ERROR("write scale factor error", "writeMetaInfo");
    return FAILURE;
  } 

  ret = SDsetattr(fpar_id, "valid_range", DFNT_UINT8, 2, g->lut_fpar.range);
  if (ret == FAILURE) {
    ERROR("write valid range error", "writeMetaInfo");
    return FAILURE;
  } 
      
  ret = SDsetattr(fpar_id, "_FillValue", DFNT_UINT8, 1, &(g->lut_fpar.fillValue));
  if (ret == FAILURE) {
    ERROR("write _FillValue error", "writeMetaInfo");
    return FAILURE;
  } 
   
  if((SDendaccess(fpar_id)) == FAILURE) {
    ERROR("SDendaccess error", "writeMetaInfo");
    return FAILURE;
  }

  /* write metadata for COMPOSE FPAR */
  if ((index = SDnametoindex(SD_ID, g->compose_fpar.name))<0) {
    ERROR("convert composed FPAR sds layer error", "writeMetaInfo");
    return FAILURE;
  }
  fpar_id = SDselect(SD_ID, index);           

  ret = SDsetattr(fpar_id, "long_name", DFNT_CHAR8, strlen(g->compose_fpar.longName), g->compose_fpar.longName);
  if (ret == FAILURE) {
    ERROR("write long name for composed FPAR error", "writeMetaInfo");
    return FAILURE;
  } 

  ret = SDsetattr(fpar_id, "scale_factor", DFNT_FLOAT64, 1, &(g->compose_fpar.scaleFactor));
  if (ret == FAILURE) {
    ERROR("write scale factor error", "writeMetaInfo");
    return FAILURE;
  } 

  ret = SDsetattr(fpar_id, "valid_range", DFNT_UINT8, 2, g->compose_fpar.range);
  if (ret == FAILURE) {
    ERROR("write valid range error", "writeMetaInfo");
    return FAILURE;
  } 
      
  ret = SDsetattr(fpar_id, "_FillValue", DFNT_UINT8, 1, &(g->compose_fpar.fillValue));
  if (ret == FAILURE) {
    ERROR("write _FillValue error", "writeMetaInfo");
    return FAILURE;
  } 

  ret = SDsetattr(fpar_id, "fill_value_legend", DFNT_CHAR8, strlen(g->modis_fpar.explain), g->modis_fpar.explain);
  if (ret == FAILURE) {
    ERROR("write description for MODIS FPAR error", "writeMetaInfo");
    return FAILURE;
  } 

  if((SDendaccess(fpar_id)) == FAILURE) {
    ERROR("SDendaccess error", "writeMetaInfo");
    return FAILURE;
  }

  /* write metadata for MODIS LAI QC */
  if ((index = SDnametoindex(SD_ID, g->modis_qc.name))<0) {
    ERROR("convert modis_qc sds layer error", "writeMetaInfo");
    return FAILURE;
  }
  lai_id = SDselect(SD_ID, index);           

  ret = SDsetattr(lai_id, "long_name", DFNT_CHAR8, strlen(g->modis_qc.longName), g->modis_qc.longName);
  if (ret == FAILURE) {
    ERROR("write long name for MODIS LAI QC error", "writeMetaInfo");
    return FAILURE;
  } 

  strcpy(str, "\n1 = MODIS high quality, timesat good fit\n2 = MODIS high quality, timesat moderate fit\n3 = MODIS low quality empirical model\n4 = not produced due to cloud etc.");
  ret = SDsetattr(lai_id, "description", DFNT_CHAR8, strlen(str), str);
  if (ret == FAILURE) {
    ERROR("Write MODIS LAI QC description error", "writeMetaInfo");
    return FAILURE;
  } 

  ret = SDsetattr(lai_id, "fill_value_legend", DFNT_CHAR8, strlen(g->modis.explain), g->modis.explain);
  if (ret == FAILURE) {
    ERROR("write description for MODIS LAI error", "writeMetaInfo");
    return FAILURE;
  } 

  if((SDendaccess(lai_id)) == FAILURE) {
    ERROR("SDendaccess error", "writeMetaInfo");
    return FAILURE;
  }

  /* write metadata for SMOOTH LAI QC */
  if ((index = SDnametoindex(SD_ID, g->timesat_qc.name))<0) {
    ERROR("convert smoothed LAI qc sds layer error", "writeMetaInfo");
    return FAILURE;
  }
  lai_id = SDselect(SD_ID, index);           

  ret = SDsetattr(lai_id, "long_name", DFNT_CHAR8, strlen(g->timesat_qc.longName), g->timesat_qc.longName);
  if (ret == FAILURE) {
    ERROR("write long name for smoothed LAI QC error", "writeMetaInfo");
    return FAILURE;
  } 

  strcpy(str, "\n1 = timesat\n2 = gap-filled\n3 = rounded\n4 = fill value");
  ret = SDsetattr(lai_id, "description", DFNT_CHAR8, strlen(str), str);
  if (ret == FAILURE) {
    ERROR("Write smoothed LAI QC description error", "writeMetaInfo");
    return FAILURE;
  } 

  if((SDendaccess(lai_id)) == FAILURE) {
    ERROR("SDendaccess error", "writeMetaInfo");
    return FAILURE;
  }

  /* write metadata for COMPOSE LAI QC */
  if ((index = SDnametoindex(SD_ID, g->compose_qc.name))<0) {
    ERROR("convert composed_qc sds layer error", "writeMetaInfo");
    return FAILURE;
  }
  lai_id = SDselect(SD_ID, index);           

  ret = SDsetattr(lai_id, "long_name", DFNT_CHAR8, strlen(g->compose_qc.longName), g->compose_qc.longName);
  if (ret == FAILURE) {
    ERROR("write long name for composed_qc error", "writeMetaInfo");
    return FAILURE;
  } 

  strcpy(str, "\n1 = high quality MODIS data\n2 = smoothed data\n3 = fill value from MODIS");
  ret = SDsetattr(lai_id, "description", DFNT_CHAR8, strlen(str), str);
  if (ret == FAILURE) {
    ERROR("Write composed LAI QC description error", "writeMetaInfo");
    return FAILURE;
  } 
  
  if((SDendaccess(lai_id)) == FAILURE) {
    ERROR("SDendaccess error", "writeMetaInfo");
    return FAILURE;
  }

  /* end */      
  if((SDend(SD_ID)) == FAILURE) {
    ERROR("SDend error", "writeMetaInfo");
    return FAILURE;
  }

  /* open again for field array write */      
  g->out_GDfid[i] = GDopen(fname, DFACC_RDWR);
  if( g->out_GDfid[i] == FAILURE){
    ERROR("opening grid file error", "writeMetaInfo");
    return FAILURE;
  } 
  g->out_GDid[i] = GDattach(g->out_GDfid[i], GridName);
  if( g->out_GDid[i] == FAILURE){
    ERROR("opening grid file error", "writeMetaInfo");
    return FAILURE;
  } 

  return SUCCESS;
}



void dec2bin(uint32 num, int *bitpattern, int limit)
{

  register int i=0;	    /*Counter*/

  for(i=0; i<limit; i++)
    bitpattern[i]=0;

  for(i=0; num>0; i++)
    {
      bitpattern[i] = (int)num & 1;
      num >>= 1;
    }
}


void bin2dec (uint32 *num, int bitpattern[], int limit)
{

  register int i=0;	    /*Counter*/
  unsigned long int x;

  *num=0;
  for(i=0; i<limit; i++)
    {
      x=bitpattern[i];
      x <<= i;
      *num = *num+x;
     }
}

void alloc_1dim_contig (void **ptr, int d1, int elsize)
{
   void *p = NULL;

   p = calloc (d1, elsize);
   if (!p) {
     ERROR ("Memory allocation error in alloc_1dim_contig", "alloc_1dim_contig");
     exit(1);
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


void free_3dim_contig (void ***a)
{
   free (a[0][0]);
   free (a[0]);
   free (a);
   return;
}

void free_2dim_contig (void **a)
{
  free (a[0]);
  free (a);
  return;
} 


void Error(const char *message, const char *module,
           const char *source, long line)
{
  fprintf(stderr, " error [%s, %s:%ld] : %s\n", module, source, line, message);
}


/**
 * create coefficient matrix and retrieve affine parameters
 * ! Input
 *     float cp[][3]  sample array 0=x, 1=y(x^2), 2 = z
 *     int num_cps: number of samples
 *     int n: order of equation used (1 in here)
 * ! Output
 *     float a[]: equation solution 
 */
void getAffine(float cp[][3], int num_cps, int n, float a[])       
{
  int i0, i1, i2, i3, j, j1, m, jf;
  float s[10][10], ss[10];

  m=(n+1)*(n+2)/2;
  for(i1=0;i1<m;i1++)
    for(i2=0;i2<m;i2++)
      s[i1][i2]=0.0;
 
  for(i0=0;i0<=n;i0++)
    for(i1=0;i1<=i0;i1++) {
      jf=i0*(i0+1)/2+i1;
      for(i2=0;i2<=n;i2++)
	for(i3=0;i3<=i2;i3++)
	  {
	    j=(i2+1)*i2/2+i3;
	    if(jf>j)
	      continue;
	    s[jf][j]=0.0;
	    for(j1=0;j1<num_cps;j1++)
	      s[jf][j]+=pow(cp[j1][0],i2-i3)*pow(cp[j1][1],i3)*pow(cp[j1][0],i0-i1)*pow(cp[j1][1],i1);
	  }
    }

  for(i1=1;i1<m;i1++)
    for(i2=0;i2<i1;i2++)
      s[i1][i2]=s[i2][i1];
  for(i1=0;i1<=n;i1++)
    for(i2=0;i2<=i1;i2++) {
      j=i1*(i1+1)/2+i2;   
      ss[j]=0.0;
      for(j1=0;j1<num_cps;j1++)
	ss[j]+=pow(cp[j1][0],i1-i2)*pow(cp[j1][1],i2)*cp[j1][2];
    }

  solveEquation(s, ss, m, a);
}

 
/**
 * find solution of a equation 
 */
void solveEquation(float s[][10], float ss[], int m, float a[])  
{
  int i1, i2, i3, ii;
  float sum, t[10][10], c;
 
  for(i1=0;i1<m;i1++) {
    for(i2=0;i2<m;i2++)
      t[i1][i2]=s[i1][i2];
    t[i1][m]=ss[i1];
  }
  
  for(i1=0;i1<m-1;i1++) {   
    c=t[i1][i1];
    ii=i1;
    for(i2=i1+1;i2<m;i2++) {
      if(t[i2][i1]>t[i1][i1])
	ii=i2; 
    }
    
    for(i3=i1;i3<=m;i3++) {
      c=t[i1][i3];
      t[i1][i3]=t[ii][i3];
      t[ii][i3]=c;
    }
    
    for(i2=i1+1;i2<m;i2++)
      for(i3=i1+1;i3<=m;i3++)
	t[i2][i3]=t[i2][i3]-t[i1][i3]*t[i2][i1]/t[i1][i1];
  }
 
  a[m-1]=t[m-1][m]/t[m-1][m-1];
  for(i1=m-2;i1>=0;i1--) {
    sum=t[i1][m];
    for(i2=m-1;i2>=i1+1;i2--)
      sum-=t[i1][i2]*a[i2];
    a[i1]=sum/t[i1][i1];
  } 
}


int compare(const void *f1, const void *f2)
{ return (*(float*) f1 > *(float*)f2) ? 1 : -1; }






