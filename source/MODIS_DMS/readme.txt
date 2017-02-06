Sharpening MODIS LST (1km) using MODIS NDVI (1km) 
(by Feng Gao on 5/21/2014)
==========================

A. Software Package 
-------------------
The MODIS surface temperature sharpening package was modified based on the Data Mining Sharpener (DMS) 
approach. NDVI (converted to fc in the package) is used for sharpening, which makes this approach 
similar to the TsHARP approach except that cubist regression may provide multiple trees (rules). 
The entire package includes one shell command (modlst_dms_1km.csh). The directories for programs, 
input and output data were pre-defined in the script. You will need to revise them in your system. 

Usage example: modlst_dms_1km.csh 2003227

B. Required Inputs
------------------
This program uses disALEXI preprocessed LST (daily MODIS swath) and NDVI (8-day composite MODIS tile) data.
There could be some inconsistent data gaps or qualities (e.g. clouds etc.) between two products.

C. LST Outputs 
--------------
A sharpened MODIS LST file and a ENVI header file will be created if there are enough samples. 
The sharpened LST file includes two layers. The first layer stores the predictions from 
regression tree model. The second layer stores energy conserved MODIS LST. 

NOTE
----
1. If there are missings/gaps from NDVI map, model prediction will have gaps in the first layer. However, 
if LST map still shows valid values, these values will be populated in the final energy conserved LST (second layer)

2. If correlation between T and NDVI is bad, the cubist will not be able to generate a regression tree. 
Instead it will use mean value for the model prediction. In such case, the first layer will only have an 
uniform value. The second layer will use original LST (no sharpening).

3. If total number of high quality T-NDVI samples is less than 50 (adjustable in shell script, in 2x2 pixels), 
no result will be generated

4. If one of NDVI or LST file is missing, program will stop

5. T-NDVI sampling and energy conservation application are applied at a 2x2 pixel size (adjustable) 
 

