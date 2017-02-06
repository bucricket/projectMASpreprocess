Generating Smoothed and Gap-filled MODIS LAI
(by Feng Gao on 5/15/2014)
=============================================

A. Software Package 
-------------------

A complete MODIS LAI processing chain includes smoothing and gap-filling processes. All processes
have been included in a C shell script (do_smooth_gapfill.csh). The directories for the programs, 
input and output data were pre-defined in the script. You will need to change them for your system. 

Usage example: do_smooth_gapfill.csh h11v04 2003

1. Smoothing LAI (timesatLAI.exe)

This program was revised based on the TIMESAT software package (in Fortran, v2.3). The original 
TIMESAT program has been modified to accept MODIS LAI input (in HDF format). The MODIS LAI quality 
flags were checked in the program. We uses Savitzky-Golay approach (window size =11) to fit MODIS LAI. 

This program can be compiled under 32-bit linux system only using static libraries under the same 
directory. The package includes all necessary libraries for 32-bit linux system. 
(make -f timesatLAI.mk)

2. Gap-filling LAI (gapfillLAI.exe)

The LAI smoothing program can only produce fitted values if there are enough high quality data in 
the time-series. No result will be produced if there are too many missing values because the fitting 
function becomes unreliable if it is forced to do so. Also, some fits may be unrealistic 
(e.g., out of the data range) due to noise or limitations of the fitting function. In these cases, 
instead of trying to fit a curve to the smaller number of high quality data points, a separate gap 
filling process was applied.  The gap-filling algorithm (in C program) searches an appropriate 
seasonal variation curve (same IGBP type) for the gap and then adjusts the seasonal variation curve 
to the sparsely available high-quality observations of the gap.

This program was written in C and can be compiled in either 32-bit or 64-bit linux system.
(make -f gapfillLAI.mk)

B. Required Inputs
------------------

- MCD15A3 (Terra and Aqua combined MODIS LAI, 4-day composite, 1km)
Note that MCD15A2, MOD15A2 or MYD15A2 (1km, 8-day composite) is also possible but need change 
the number of files (46 instead of 92 production periods per year)
 
- MCD12Q1 (annual land cover map in 500m)
The gap-filling program reads 500m IGBP data and then aggregates (majority) to 1km resolution 
automatically before it's been used for processing 1km MODIS LAI. 

C. LAI Outputs 
--------------

Output includes original, smoothed and composite MODIS LAI (and Fpar) for each production period 
as well as QA layers. Both binary and HDF format files are produced. The binary file comes with 
an ENVI header in BIP format. The HDF format contains necessary metadata information and can be 
re-projected and mosaiced using MRT software. The header from HDF file includes details of the LAI 
value and its quality flags (see an example below).

netcdf MOD15COM.A2009125.h12v04.005 {
dimensions:
	YDim:MODIS_NACP_LAI = 1200 ;
	XDim:MODIS_NACP_LAI = 1200 ;

variables:
	byte MODIS_LAI(YDim:MODIS_NACP_LAI, XDim:MODIS_NACP_LAI) ;
		MODIS_LAI:long_name = "Original MODIS LAI" ;
		MODIS_LAI:scale_factor = 0.1 ;
		MODIS_LAI:valid_range = '\0', 'd' ;
		MODIS_LAI:_FillValue = '\377' ;
		MODIS_LAI:fill_value_legend = "MOD15A2 FILL VALUE LEGEND\n",
    "255 = _Fillvalue, assigned when:\n",
    "    * the MODAGAGG suf. reflectance for channel VIS, NIR was assigned its _Fillvalue, or\n",
    "    * land cover pixel itself was assigned _Fillvalus 255 or 254.\n",
    "254 = land cover assigned as perennial salt or inland fresh water.\n",
    "253 = land cover assigned as barren, sparse vegetation (rock, tundra, desert.)\n",
    "252 = land cover assigned as perennial snow, ice.\n",
    "251 = land cover assigned as \"permanent\" wetlands/inundated marshlands.\n",
    "250 = land cover assigned as urban/built-up.\n",
    "249 = land cover assigned as \"unclassified\" or not able to determine." ;
	byte Smoothed_LAI(YDim:MODIS_NACP_LAI, XDim:MODIS_NACP_LAI) ;
		Smoothed_LAI:long_name = "Gap-Filled and Timesat-Smoothed LAI" ;
		Smoothed_LAI:scale_factor = 0.1 ;
		Smoothed_LAI:valid_range = '\0', 'd' ;
		Smoothed_LAI:_FillValue = '\377' ;
	byte Composed_LAI(YDim:MODIS_NACP_LAI, XDim:MODIS_NACP_LAI) ;
		Composed_LAI:long_name = "Composed LAI from MODIS and Timesat" ;
		Composed_LAI:scale_factor = 0.1 ;
		Composed_LAI:valid_range = '\0', 'd' ;
		Composed_LAI:_FillValue = '\377' ;
		Composed_LAI:fill_value_legend = "MOD15A2 FILL VALUE LEGEND\n",
    "255 = _Fillvalue, assigned when:\n",
    "    * the MODAGAGG suf. reflectance for channel VIS, NIR was assigned its _Fillvalue, or\n",
    "    * land cover pixel itself was assigned _Fillvalus 255 or 254.\n",
    "254 = land cover assigned as perennial salt or inland fresh water.\n",
    "253 = land cover assigned as barren, sparse vegetation (rock, tundra, desert.)\n",
    "252 = land cover assigned as perennial snow, ice.\n",
    "251 = land cover assigned as \"permanent\" wetlands/inundated marshlands.\n",
    "250 = land cover assigned as urban/built-up.\n",
    "249 = land cover assigned as \"unclassified\" or not able to determine." ;
	byte MODIS_FPAR(YDim:MODIS_NACP_LAI, XDim:MODIS_NACP_LAI) ;
		MODIS_FPAR:long_name = "Original MODIS FPAR" ;
		MODIS_FPAR:scale_factor = 0.01 ;
		MODIS_FPAR:valid_range = '\0', 'd' ;
		MODIS_FPAR:_FillValue = '\377' ;
		MODIS_FPAR:fill_value_legend = "MOD15A2 FILL VALUE LEGEND\n",
    "255 = _Fillvalue, assigned when:\n",
    "    * the MODAGAGG suf. reflectance for channel VIS, NIR was assigned its _Fillvalue, or\n",
    "    * land cover pixel itself was assigned _Fillvalus 255 or 254.\n",
    "254 = land cover assigned as perennial salt or inland fresh water.\n",
    "253 = land cover assigned as barren, sparse vegetation (rock, tundra, desert.)\n",
    "252 = land cover assigned as perennial snow, ice.\n",
    "251 = land cover assigned as \"permanent\" wetlands/inundated marshlands.\n",
    "250 = land cover assigned as urban/built-up.\n",
    "249 = land cover assigned as \"unclassified\" or not able to determine." ;
	byte Converted_FPAR(YDim:MODIS_NACP_LAI, XDim:MODIS_NACP_LAI) ;
		Converted_FPAR:long_name = "Converted FPAR from LAI-FPAR LUT" ;
		Converted_FPAR:scale_factor = 0.01 ;
		Converted_FPAR:valid_range = '\0', 'd' ;
		Converted_FPAR:_FillValue = '\377' ;
	byte Composed_FPAR(YDim:MODIS_NACP_LAI, XDim:MODIS_NACP_LAI) ;
		Composed_FPAR:long_name = "Composed FPAR from MODIS and Converted FPAR" ;
		Composed_FPAR:scale_factor = 0.01 ;
		Composed_FPAR:valid_range = '\0', 'd' ;
		Composed_FPAR:_FillValue = '\377' ;
		Composed_FPAR:fill_value_legend = "MOD15A2 FILL VALUE LEGEND\n",
    "255 = _Fillvalue, assigned when:\n",
    "    * the MODAGAGG suf. reflectance for channel VIS, NIR was assigned its _Fillvalue, or\n",
    "    * land cover pixel itself was assigned _Fillvalus 255 or 254.\n",
    "254 = land cover assigned as perennial salt or inland fresh water.\n",
    "253 = land cover assigned as barren, sparse vegetation (rock, tundra, desert.)\n",
    "252 = land cover assigned as perennial snow, ice.\n",
    "251 = land cover assigned as \"permanent\" wetlands/inundated marshlands.\n",
    "250 = land cover assigned as urban/built-up.\n",
    "249 = land cover assigned as \"unclassified\" or not able to determine." ;
	byte MODIS_LAI_FPAR_QC(YDim:MODIS_NACP_LAI, XDim:MODIS_NACP_LAI) ;
		MODIS_LAI_FPAR_QC:long_name = "MODIS LAI/FPAR Retrieval Quality" ;
		MODIS_LAI_FPAR_QC:description = "\n",
    "1 = MODIS high quality, timesat good fit\n",
    "2 = MODIS high quality, timesat moderate fit\n",
    "3 = MODIS low quality empirical model\n",
    "4 = not produced due to cloud etc." ;
		MODIS_LAI_FPAR_QC:fill_value_legend = "MOD15A2 FILL VALUE LEGEND\n",
    "255 = _Fillvalue, assigned when:\n",
    "    * the MODAGAGG suf. reflectance for channel VIS, NIR was assigned its _Fillvalue, or\n",
    "    * land cover pixel itself was assigned _Fillvalus 255 or 254.\n",
    "254 = land cover assigned as perennial salt or inland fresh water.\n",
    "253 = land cover assigned as barren, sparse vegetation (rock, tundra, desert.)\n",
    "252 = land cover assigned as perennial snow, ice.\n",
    "251 = land cover assigned as \"permanent\" wetlands/inundated marshlands.\n",
    "250 = land cover assigned as urban/built-up.\n",
    "249 = land cover assigned as \"unclassified\" or not able to determine." ;
	byte Smoothed_LAI_FPAR_QC(YDim:MODIS_NACP_LAI, XDim:MODIS_NACP_LAI) ;
		Smoothed_LAI_FPAR_QC:long_name = "Gap-Filled and Timesat-Smoothed LAI/FPAR Quality" ;
		Smoothed_LAI_FPAR_QC:description = "\n",
    "1 = timesat\n",
    "2 = gap-filled\n",
    "3 = rounded\n",
    "4 = fill value" ;
	byte Composed_LAI_FPAR_QC(YDim:MODIS_NACP_LAI, XDim:MODIS_NACP_LAI) ;
		Composed_LAI_FPAR_QC:long_name = "Composed LAI/FPAR Quality" ;
		Composed_LAI_FPAR_QC:description = "\n",
    "1 = high quality MODIS data\n",
    "2 = smoothed data\n",
    "3 = fill value from MODIS" ;
    
