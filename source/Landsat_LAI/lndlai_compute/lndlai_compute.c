/**
 * !Description
 * Source Filename: Indvi_biosphysics.c ( Landsat bio-physics parameter product)
 * This program calculate NDVI and EVI for one landsat image
 *
 * !Inputs: 
 *   Landsat surface reflectance
 * 
 *   A sample input parameter file looks like
 *
 PARAMETER_FILE
 
 # define input Landsat and MODIS surface reflectance
 LANDSAT_BASE_FILE = lndsr.05242001.hdf 
 BIOPHYSICS_PARA_FILE_OUT = biopara.A2001185.hdf

 END
 *
 * !Output: L[2]
 *  Landsat NDVI and EVI
 * 
 * !Developer
 *Yujie Wang (ywang@pop900.gsfc.nasa.gov)
 *
 * !Revision
 * Original version - 10/06 by Yujie Wang (NASA/GSFC through ERT)
 */

#include "lndlai.h"
#include "defns.h"

extern Description Case;

void GetBioPhyPara(float64 red, float64 nir, float64 blue, float64* ndvi,float64* evi);
int InitLAIModel(char *namestem);
float Predict_LAI();
void CleanUpLAIModel();
void DefaultMeta(META_SDS_SR* meta);

int main(int argc, char *argv[])
{
	INPUT_PARS        *pars;
	META_LANDSAT_SR   *metaGD;
	GRID_LANDSAT_SR   srLandsat;
	GRID_LAI          srBio;
	META_SDS_SR	  metaSDSNDVI, metaSDSEVI, metaSDSLai, metaSDSQA;
	int16 *pbfBands[NBANDS];
	int16 *pbfNDVI, *pbfEVI, *pbfLai;
	uint8 *cfmask;
	int i, buffsize, irow, n, flag, iband;
	float64 blue, red, nir, ndvi, evi;
	float temp;
	char qaattr[MAX_STRING_LENGTH];

	/* allocate memory for variables */
	pars = malloc(sizeof(INPUT_PARS));
	metaGD = malloc(sizeof(META_LANDSAT_SR));
	openLog("Landsat Bio-physics parameters");
	if(!pars || !metaGD)
	  {
	    sprintf (msg, "Not successful in allocating memory for grid metadata");
	    ERRORMSG(msg, "main");
	    closeLog();
	    return FAILURE;
	  }
	if(argc==2) 
	  {
	    /* get input parameters from input file */
	    parseParameters(argv[1], pars);
	  }
	else 
	  {
	    printf("Usage: %s <input_parameter_file>\n", argv[0]);
	    exit(1);
	  }

	GDSRInit(&srBio);
	for(i=0; i<NBANDS; i++) 
	  strcpy(srLandsat.fileName[i], pars->LandsatFile[i]);
	strcpy(srBio.fileName, pars->LandsatBioPhysFile);
	InitLAIModel(pars->LandsatAncFile);
	free(pars);

	if(getLandsatMetaInfo(metaGD, &srLandsat) == FAILURE)
	  {
	    sprintf(msg, "Retrieve Landsat %s metadata error", srLandsat.fileName);
	    ERRORMSG(msg, "main");
	    return FAILURE;
	  }
	strcpy(srBio.gdname, "LANDSAT");
	if(WriteLandsatMetaInfo(metaGD, &srBio) == FAILURE)
	  {
	    sprintf(msg, "create output file %s metadata error", srBio.fileName);
	    ERRORMSG(msg, "main");
	    return FAILURE;
	  }

        metaSDSQA.ncols= metaSDSNDVI.ncols = metaSDSEVI.ncols =metaSDSLai.ncols = metaGD->ncols;
        metaSDSQA.nrows= metaSDSNDVI.nrows = metaSDSEVI.nrows =metaSDSLai.nrows = metaGD->nrows;

	DefaultMeta(&metaSDSNDVI);
	DefaultMeta(&metaSDSEVI);
	DefaultMeta(&metaSDSLai);
	metaSDSLai.range[0] = (int)(0.01/metaSDSLai.scale+0.5); 
	metaSDSLai.range[1] = (int)(10.0/metaSDSLai.scale+0.5);

	if(CreateLandsatField(&srBio, "NDVI", DFNT_INT16) == FAILURE) 
	  {
	    sprintf(msg, "create grid field NDVI error");
	    ERRORMSG(msg, "main");
	    return FAILURE;
	  }
	if(CreateLandsatField(&srBio, "EVI", DFNT_INT16) == FAILURE) 
	  {
	    sprintf(msg, "create grid field EVI error");
	    ERRORMSG(msg, "main");
	    return FAILURE;
	  }
	if(CreateLandsatField(&srBio, "LAI", DFNT_INT16) == FAILURE) 
	  {
	    sprintf(msg, "create grid field LAI error");
	    ERRORMSG(msg, "main");
	    return FAILURE;
	  }

	if(CreateLandsatField(&srBio, "cfmask", DFNT_INT8) == FAILURE) 
	  {
	    sprintf(msg, "create grid field indsr_QA error");
	    ERRORMSG(msg, "main");
	    return FAILURE;
	  }

	
	buffsize = metaGD->ncols* NROWSTEP;
	for(i=0; i<NBANDS; i++)
	  {
	    pbfBands[i] = (int16*)malloc(sizeof(int16)*buffsize);
	    if(!pbfBands[i])
	      {
		sprintf(msg, "Not successful in allocating buffer memory");
		ERRORMSG(msg, "main");
		return FAILURE;
	      }
	  }
	cfmask = (uint8*)malloc(sizeof(uint8)*buffsize);
	pbfLai = (int16*)malloc(sizeof(int16)*buffsize);
	pbfNDVI = (int16*)malloc(sizeof(int16)*buffsize);
	pbfEVI = (int16*)malloc(sizeof(int16)*buffsize);
	if(!pbfEVI)
	  {
	    sprintf(msg, "Not successful in allocating buffer memory");
	    ERRORMSG(msg, "main");
	    return FAILURE;
	  }

	irow =0;
	while(1) 
	  {
	    flag =0;
	    for(i=0; i<NBANDS; i++)
	      {
		if(i==CLOUD)
		  n=ReadENVInRow(&srLandsat, CLOUD, irow, metaGD->nrows, metaGD->ncols, NROWSTEP, cfmask);
		else
		  n=ReadENVInRow(&srLandsat, i, irow, metaGD->nrows, metaGD->ncols, NROWSTEP, pbfBands[i]);
		if(n<=0)
		  {
		    flag = 1;
		    break;
		  }
	      }
	    if(flag == 1) break;

	    for(i=0; i<buffsize; i++)
	      {
			
		if(cfmask[i] != 255)
		  {
			
		    if(pbfBands[RED][i] == srLandsat.fillv[RED] || pbfBands[NIR][i] == srLandsat.fillv[NIR])
		      {
			pbfNDVI[i] = metaSDSNDVI.fillvalue;
			pbfEVI[i] = metaSDSEVI.fillvalue;
		      }
		    else
		      {
			blue = pbfBands[BLUE][i] * srLandsat.scale[BLUE];
			red = pbfBands[RED][i] * srLandsat.scale[RED];
			nir = pbfBands[NIR][i] * srLandsat.scale[NIR];
			GetBioPhyPara(red, nir, blue, &ndvi, &evi);
			if(ndvi >-9990)
			  pbfNDVI[i] = (int16) (ndvi/metaSDSNDVI.scale);
			else
			  pbfNDVI[i] = metaSDSNDVI.fillvalue;
			if(evi >-9990)
			  pbfEVI[i] = (int16) (evi/metaSDSEVI.scale);
			else
			  pbfEVI[i] = metaSDSEVI.fillvalue;
		      }
		  
		    /* exclude water and clouds */
		    if(cfmask[i] != 0)
		      pbfLai[i] = metaSDSLai.fillvalue;
		    else
		      {
			flag =0;
			for(iband=0; iband<NBANDS-1; iband++)
			  {
			    if(pbfBands[iband][i]<=0) 
			      {
				flag =1;
				break;
			      }
			    /* assign value for this Case */
			    CVal(Case, iband+1) = pbfBands[iband][i];
					
			  }

			if(flag == 1) 
			  {
			    pbfLai[i] = metaSDSLai.fillvalue;
			    continue;
			  }

			/* include NDVI and NDWI */
			CVal(Case, NBANDS) = (float)(pbfBands[NIR][i]-pbfBands[RED][i])/(pbfBands[NIR][i]+pbfBands[RED][i]);
			CVal(Case, NBANDS+1) = (float)(pbfBands[NIR][i]-pbfBands[SWIR1][i])/(pbfBands[NIR][i]+pbfBands[SWIR1][i]);
			
			temp = Predict_LAI();
				
			if(temp >0.01)
			  /*pbfLai[i] = (int16)(temp/metaSDSLai.scale);*/
			  pbfLai[i] = (int16)(temp/metaSDSLai.scale+0.5); 
			else
			  /*pbfLai[i] = metaSDSLai.fillvalue;*/
			  // if fails then uses a small value (0.01) - Feng (2/29/16) 
			  pbfLai[i] = metaSDSLai.range[0];
		      }
		  }
		else
		  {
		    pbfNDVI[i] = metaSDSNDVI.fillvalue;
		    pbfEVI[i] = metaSDSEVI.fillvalue;
		    pbfLai[i] = metaSDSLai.fillvalue;
		  }

	      }
		
	    WritenRow(&srBio, "NDVI", irow, metaSDSNDVI.nrows, metaSDSNDVI.ncols, NROWSTEP, pbfNDVI);
	    WritenRow(&srBio, "EVI", irow, metaSDSEVI.nrows, metaSDSEVI.ncols, NROWSTEP, pbfEVI);
	    WritenRow(&srBio, "cfmask", irow, metaSDSQA.nrows, metaSDSQA.ncols, NROWSTEP, cfmask);
	    WritenRow(&srBio, "LAI", irow, metaSDSLai.nrows, metaSDSLai.ncols, NROWSTEP, pbfLai);

	    printf("%4d\b\b\b\b",irow);
	    irow += n;
	  }

	Cleanup(&srBio);

	WriteSDSAttr(&srBio, &metaSDSNDVI, "NDVI");
	WriteSDSAttr(&srBio, &metaSDSEVI, "EVI");
	WriteSDSAttr(&srBio, &metaSDSLai, "LAI");
	WriteQAAttr(&srBio, &metaSDSQA, "cfmask");
	//	WriteQAAttr(&srBio, &metaSDSQA, "cfmask", qaattr);

	free(metaGD);
	for(i=0; i<NBANDS; i++) free(pbfBands[i]);
	free(cfmask);
	free(pbfNDVI);
	free(pbfEVI);
	/*CleanUpLAIModel();*/
	return SUCCESS;
	
}


void GetBioPhyPara(float64 red, float64 nir, float64 blue, float64* ndvi,float64* evi)
{
  double G=2.5, L=1, C1=6, C2=7.5; 
	
  if(red>0 && nir>0) 
    {
      *ndvi = (nir-red)/(nir+red);
      if(blue >0)
	*evi = G*(nir-red)/(nir+C1*red - C2*blue +L);
      else *evi = -9999;
    }
  else
    {
      *ndvi = -9999;
      *evi = -9999;
    }
}


void DefaultMeta(META_SDS_SR* meta)
{
  strcpy(meta->unit, "none");
  meta->longname[0] = '\0';
  meta->range[0] = -1000;
  meta->range[1] = 1000;
  meta->fillvalue = -9999;
  meta->scale = 0.001; 
}


































