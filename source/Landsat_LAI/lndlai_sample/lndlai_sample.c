/**
 * !Description
 * Source Filename: laisample.c ( Landsat & MODIS LAI matching sample data )
 * This program match Landsat and MODIS pixels, select relatively homogenious 
 * MODIS pixels, output the relationship table of mean Landsat surface reflectance 
 * and MODIS LAI.
 *
 * !Inputs: 
 *   Landsat surface reflectance
 *	 MODIS LAI product (MOD15A2)
 *   
 *
 PARAMETER_FILE
 
 # define input Landsat and MODIS surface reflectance
 LANDSAT_BASE_FILE = lndsr.05242001.hdf 
 MODIS_BASE_FILE = MOD15A2.A2001241.h12v04.004.2003148062201.hdf
 SAMPLE_FILE_OUT = sample.dat

 END
 *
 * !Output: L[2]
 *  Landsat surface reflectance and MODIS LAI relations table
 * 
 * !Developer
 *Yujie Wang (ywang@pop900.gsfc.nasa.gov)
 *
 * !Revision
 * revision - 02/15 by Feng Gao (reads Landsat SR in ENVI format)
 * revision - 09/12 by Feng Gao (reads in extra QC bits to constrain HQ LAI samples) 
 * Revision - 08/12 by Feng Gao (adds sample location in Landsat for combing with in-situ measurements
 * Original version - 10/06 by Yujie Wang (NASA/GSFC through ERT)
 */

#include "lndlai.h"
extern char BandName[NBANDS][20];

int main(int argc, char *argv[])
{
  INPUT_PARS        *pars;
  META_LANDSAT_SR   *metaGD;
  META_MODIS_LAI    metaLai;
  GRID_LANDSAT_SR   srLandsat;
  GRID_LAI          gdLai;
  int16 *pbfBand[NBANDS];
  uint8 *cfmask;
  int i, j, k, buffsize, irow, n, Modiscol, Modisrow, kk, jj;
  int flag, retrieval_QA, cloudstate;
  double var, var_av, sigma, weight, ndvi, ndwi, mean_coeff_var;
  double  Threshold[NBANDS-1] = {200, 200, 200, 500, 400, 400};
  FILE *fp, *tmp;
  MODIS_CELL* mc;
  uint8* pbfLai, *pbfLaiQA, *pbfLaiExtraQA;
  char str[1000], sat_flag;
	
  int half_pixel;

  /*tmp=fopen("mask_test.bin", "wb");*/

  /* allocate memory for variables */
  pars = malloc(sizeof(INPUT_PARS));
  metaGD = malloc(sizeof(META_LANDSAT_SR));
  openLog("Landsat surface reflectance-MODIS LAI match up");
  if(!pars || !metaGD)
    {
      sprintf (msg, "Not successful in allocating memory for grid metadata");
      ERRORMSG(msg, "main");
      closeLog();
      return FAILURE;
    }

  pars->pure_threshold = 0.15;  /* set default value - Feng (10/12) */
  if(argc==2) 
    {
      parseParameters(argv[1], pars);
    }
  else 
    {
      printf("Usage: %s <input_parameter_file>\n", argv[0]);
      exit(1);
    }

  GDSRInit(&gdLai);
  for(i=0; i<NBANDS; i++) 
    strcpy(srLandsat.fileName[i], pars->LandsatFile[i]);
  strcpy(gdLai.fileName, pars->ModisFile);
  if( (fp = fopen(pars->SampleFile, "a+")) == NULL)
    {
      sprintf (msg, "Not successful in opening file %s", pars->SampleFile);
      ERRORMSG(msg, "main");
      closeLog();
      return FAILURE;
    }

  if(getLandsatMetaInfo(metaGD, &srLandsat) == FAILURE)
    {
      sprintf(msg, "Retrieve Landsat %s metadata error", srLandsat.fileName);
      ERRORMSG(msg, "main");
      return FAILURE;
    }

  if(getMODISLAIMetaInfo(&gdLai, &metaLai) == FAILURE)
    {
      sprintf(msg, "Retrieve MODIS LAI %s metadata error", gdLai.fileName);
      ERRORMSG(msg, "main");
      return FAILURE;
    }

  buffsize = metaGD->ncols* NROWSTEP;
  cfmask = (uint8*)malloc(sizeof(uint8)*buffsize);
  for(i=0; i<NBANDS; i++)
    {
      pbfBand[i] = (int16*)malloc(sizeof(int16)*buffsize);
      if(!pbfBand[i])
	{
	  sprintf(msg, "Not successful in allocating buffer memory");
	  ERRORMSG(msg, "main");
	  return FAILURE;
	}
    }
  pbfLai = (uint8*)malloc(sizeof(uint8)*metaLai.ncols*metaLai.nrows);
  if(!pbfLai)
    {
      sprintf(msg, "Not successful in allocating buffer memory for MODIS LAI data");
      ERRORMSG(msg, "main");
      return FAILURE;
    }

  pbfLaiQA = (uint8*)malloc(sizeof(uint8)*metaLai.ncols*metaLai.nrows);
  if(!pbfLaiQA)
    {
      sprintf(msg, "Not successful in allocating buffer memory for MODIS LAI QA data");
      ERRORMSG(msg, "main");
      return FAILURE;
    }
  pbfLaiExtraQA = (uint8*)malloc(sizeof(uint8)*metaLai.ncols*metaLai.nrows);
  if(!pbfLaiExtraQA)
    {
      sprintf(msg, "Not successful in allocating buffer memory for MODIS LAI Extra QA data");
      ERRORMSG(msg, "main");
      return FAILURE;
    }

  mc = (MODIS_CELL*)malloc(sizeof(MODIS_CELL)*metaLai.ncols*metaLai.nrows);
  if(!mc)
    {
      sprintf(msg, "Not successful in allocating buffer memory for MODIS cells");
      ERRORMSG(msg, "main");
      return FAILURE;
    }
  for(i=0; i<metaLai.ncols*metaLai.nrows; i++)
    InitModisCell(&(mc[i]));

  InitLandsatMODISProjInv(metaGD, &metaLai);
  n=ReadnRow(&gdLai, "Lai_1km", 0, metaLai.nrows, metaLai.ncols, metaLai.nrows, pbfLai);
  n=ReadnRow(&gdLai, "FparLai_QC", 0, metaLai.nrows, metaLai.ncols, metaLai.nrows, pbfLaiQA);
  n=ReadnRow(&gdLai, "FparExtra_QC", 0, metaLai.nrows, metaLai.ncols, metaLai.nrows, pbfLaiExtraQA);
  irow =0;
  while(1)
    {
      flag =0;
      for(i=0; i<metaGD->ncols*NROWSTEP; i++)  cfmask[i] = srLandsat.fillv[CLOUD];
      for(i=0; i<NBANDS; i++)
	{
	  if(i==CLOUD)
	    n=ReadENVInRow(&srLandsat, CLOUD, irow, metaGD->nrows, metaGD->ncols, NROWSTEP, cfmask);
	  else
	    n=ReadENVInRow(&srLandsat, i, irow, metaGD->nrows, metaGD->ncols, NROWSTEP, pbfBand[i]);
	  if(n<=0)
	    {
	      flag = 1;
	      break;
	    }
	}

      if(flag == 1) break;
      for(i=0; i<metaGD->ncols*NROWSTEP; i++)  pbfBand[CLOUD][i] = cfmask[i];
      for(i=0; i<n; i++)
	{
	  for(j=0; j<metaGD->ncols; j++)
	    {
	      jj = i*metaGD->ncols + j;

	      flag = 0;
	      for(k=0; k<NBANDS-1; k++)
		if(pbfBand[k][jj] == srLandsat.fillv[k] || pbfBand[k][jj] < 0) 
		  {
		    flag = 1;
		    break;
		  }

	      /*fwrite(&flag, 1, 1, tmp);*/
	      if(flag == 1 || pbfBand[CLOUD][jj] > 1) continue;
	      Landsat2MODIS(metaGD, irow+i, j, &metaLai, &Modisrow, &Modiscol);
	      //printf("%d %d\n", Modisrow, Modiscol);
	      if(Modisrow<0 || Modisrow>=metaLai.nrows || Modiscol<0 || Modiscol>=metaLai.ncols) continue;
	      kk = Modisrow*metaLai.ncols + Modiscol;
	      for(k=0; k<NBANDS-1; k++)
		{
		  mc[kk].sum[k] += pbfBand[k][jj];
		  mc[kk].sum2[k] += pbfBand[k][jj]*pbfBand[k][jj];
		}

	      mc[kk].count++;
	      /* count land pixel */
	      if(cfmask[jj] !=1 ) mc[kk].land ++;
	      /* only last pixel in mc[kk] will be saved, suppose it's in the lower right corner
		 then I need to move half pixel to the center of MODIS cell approxiamtely - Feng (08/12) */  
	      mc[kk].irow = irow+i;
	      mc[kk].icol = j;
	      /*if(mc[kk].irow==3603&&mc[kk].icol==4086) printf("%d %d %d\n", kk, mc[kk].count, mc[kk].land);*/ 

	    }
	}
      printf("%5d\b\b\b\b\b",irow);
      irow += n;
    }

  for(i=0; i<metaLai.ncols*metaLai.nrows; i++)
    {
      if(mc[i].count<900) continue;

      /* check MODIS LAI QA flag 
	 "MODLAND_QC START 0 END 0 VALIDS 2\n",
	 "MODLAND_QC   0 = Good Quality (main algorithm with or without saturation)\n",
	 "MODLAND_QC   1 = Other Quality (back-up algorithm or fill value)\n",
	 "SENSOR START 1 END 1 VALIDS 2\n",
	 "SENSOR       0  = Terra\n",
	 "SENSOR       1  = Aqua\n",
	 "DEADDETECTOR START 2 END 2 VALIDS 2\n",
	 "DEADDETECTOR 0 = Detectors apparently fine for up to 50% of channels 1,2\n",
	 "DEADDETECTOR 1 = Dead detectors caused >50% adjacent detector retrieval\n",
	 "CLOUDSTATE START 3 END 4 VALIDS 4 (this inherited from Aggregate_QC bits {0,1} cloud state)\n",
	 "CLOUDSTATE   00 = 0 Significant clouds NOT present (clear)\n",
	 "CLOUDSTATE   01 = 1 Significant clouds WERE present\n",
	 "CLOUDSTATE   10 = 2 Mixed cloud present on pixel\n",
	 "CLOUDSTATE   11 = 3 Cloud state not defined,assumed clear\n",
	 "SCF_QC START 5 END 7 VALIDS 5\n",
	 "SCF_QC       000=0 Main (RT) algorithm used, best result possible (no saturation)\n",
	 "SCF_QC       001=1 Main (RT) algorithm used, saturation occured. Good, very usable.\n",
	 "SCF_QC       010=2 Main algorithm failed due to bad geometry, empirical algorithm used\n",
	 "SCF_QC       011=3 Main algorithm faild due to problems other than geometry, empirical algorithm used\n",
	 "SCF_QC       100=4 Pixel not produced at all, value coudn\'t be retrieved (possible reasons: bad L1B data, unusabl
	 e MODAGAGG data)" */

      /* include homogeneous water pixels to help building a balanced regression */
      /*if(mc[i].land ==0) {
	pbfLai[i] = 0;
	pbfLaiQA[i] = 0;
	pbfLaiExtraQA[i] = 0;
	pars->pure_threshold = 0.01;
	}
	else */
      if(ExtractBit(pbfLaiQA[i], 0, 1) !=0 || mc[i].land<900) continue;  /* main algorithm and land only */
	  
      if(pbfLai[i]>100) continue; /* valid MODIS LAI */

      if(ExtractBit(pbfLaiQA[i], 2, 1) !=0) continue;  /* detectors must be good */
      cloudstate = ExtractBit(pbfLaiQA[i], 3, 2);
      if(cloudstate==1 || cloudstate==2 ) continue;    /* must be clear */

      retrieval_QA = ExtractBit(pbfLaiQA[i], 5, 3);
      /* main alg. no saturation */
      if(retrieval_QA == 0) {
	weight = 1.0; 
	sat_flag = 'N';       
      } 
      /* main alg. with saturation */
      else {
	sat_flag = 'S';	    
	if(pbfLai[i]<55) weight = 0.50;
	else if(pbfLai[i]<60) weight = 0.25; 
	else if(pbfLai[i]<65) weight = 0.10; 
	else continue;
      }
      /* else continue;*/
		
      /* "FparExtra_QC 6 BITFIELDS IN 8 BITWORD\n",
	 "LANDSEA PASS-THROUGH START 0 END 1 VALIDS 4\n",
	 "LANDSEA   00 = 0 LAND       AggrQC(3,5)values{001}\n",
	 "LANDSEA   01 = 1 SHORE      AggrQC(3,5)values{000,010,100}\n",
	 "LANDSEA   10 = 2 FRESHWATER AggrQC(3,5)values{011,101}\n",
	 "LANDSEA   11 = 3 OCEAN      AggrQC(3,5)values{110,111}\n",
	 "SNOW_ICE (from Aggregate_QC bits) START 2 END 2 VALIDS 2\n",
	 "SNOW_ICE  0 = No snow/ice detected\n",
	 "SNOW_ICE  1 = Snow/ice were detected\n",
	 "AEROSOL START 3 END 3 VALIDS 2\n",
	 "AEROSOL   0 = No or low atmospheric aerosol levels detected\n",
	 "AEROSOL   1 = Average or high aerosol levels detected\n",
	 "CIRRUS (from Aggregate_QC bits {8,9} ) START 4 END 4 VALIDS 2\n",
	 "CIRRUS    0 = No cirrus detected\n",
	 "CIRRUS    1 = Cirrus was detected\n",
	 "INTERNAL_CLOUD_MASK START 5 END 5 VALIDS 2\n",
	 "INTERNAL_CLOUD_MASK 0 = No clouds\n",
	 "INTERNAL_CLOUD_MASK 1 = Clouds were detected\n",
	 "CLOUD_SHADOW START 6 END 6 VALIDS 2\n",
	 "CLOUD_SHADOW        0 = No cloud shadow detected\n",
	 "CLOUD_SHADOW        1 = Cloud shadow detected\n",
	 "SCF_BIOME_MASK START 7 END 7 VALIDS 2\n",
	 "SCF_BIOME_MASK  0 = Biome outside interval <1,4>\n",
	 "SCF_BIOME_MASK  1 = Biome in interval <1,4>" ; */

      /* LAI extra QC bits from 2-6 must be all zeros - Feng (9/28/2012) */
      if(ExtractBit(pbfLaiExtraQA[i], 2, 5) !=0) continue; 

      flag =0;
      var_av=0;
      for(k=0; k<NBANDS-1; k++)
	{
	      
	  mc[i].sum[k] /=mc[i].count;
	  mc[i].sum2[k] /=mc[i].count;
	  sigma = mc[i].sum2[k] - mc[i].sum[k]*mc[i].sum[k];
	  var = sqrt(sigma);
	  /* The coefficient of variation (CV) is defined as the ratio of the standard deviation to the mean */
	  var_av += var/mc[i].sum[k];
	  /*if(var> Threshold[k]) 
	    {
	    flag =1; 
	    break;
	    }*/
	}

      /* keep some higher LAI values - Feng (10/11) */
      /*if(pbfLai[i]>30 && var_av/QA<300)
	flag = 0;

	if(flag ==0 && pbfLai[i]/10.0<6.0) */
      mean_coeff_var = var_av/(NBANDS-1);
      /* keep some higher LAI values (> 3.0)  - Feng (10/11) */
      /*if(mean_coeff_var < PURE_TH || (pbfLai[i]>30 && mean_coeff_var<PURE_TH*((pbfLai[i]-30.0)/80.0+1.0) ))*/
      if(mean_coeff_var < pars->pure_threshold) 
	{
	      
	  /* convert to signal-to-noise ratio and use it as a weight */
	  if(mean_coeff_var > 0.01)
	    weight *= 0.01/mean_coeff_var;
	  else
	    weight *= 1.0;
	  /* compute approxiamte location in Landsat image - Feng (08/12) */
	  half_pixel = (int)(0.5 * metaLai.res / metaGD->res + 0.5); 
	  fprintf(fp, "%7.1f, %7.1f, ", metaGD->ulx+(mc[i].irow-half_pixel)*metaGD->res, 
		  metaGD->uly-(mc[i].icol-half_pixel)*metaGD->res);
	      
	  for(j=0; j<NBANDS-1; j++)
	    {
	      fprintf(fp, "%6.1f, ", mc[i].sum[j]);
	    }
	  ndvi = (float)(mc[i].sum[3]-mc[i].sum[2])/(mc[i].sum[3]+mc[i].sum[2]);
	  ndwi = (float)(mc[i].sum[3]-mc[i].sum[4])/(mc[i].sum[3]+mc[i].sum[4]);
	  fprintf(fp, "%7.4f, %7.4f, %5.3f, %6.4f %c\n", ndvi, ndwi, pbfLai[i]/10.0, weight, sat_flag);
	  /*fprintf(fp, "%7.4f, %7.4f, %5.3f, %6.4f, %d, %d\n", 
	    ndvi, ndwi, pbfLai[i]/10.0, weight, pbfLaiQA[i], pbfLaiExtraQA[i]);*/
	}
		
    }
	
  Cleanup(&gdLai);
  //	Cleanup_Landsat_SR(&gdLai);

  free(pars);
  free(metaGD);
  for(i=0; i<NBANDS; i++) free(pbfBand[i]);
  free(cfmask);
  free(pbfLai);
  free(pbfLaiQA);
  free(pbfLaiExtraQA);
  fclose(fp);
  closeLog();
  return SUCCESS;
	
}

