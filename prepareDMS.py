#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 30 09:31:45 2016

@author: mschull
"""

import landsatTools
import os
import glob
import datetime
import subprocess

# set project base directory structure
p = os.getcwd()
base = os.path.abspath(os.path.join(p,os.pardir))
#base = os.path.join(os.sep,'Users','mschull','umdGD','pyDisALEXI')
modisBase = os.path.join(base,'data','MODIS')
dataBase = os.path.join(base,'data')
landsatDataBase = os.path.join(dataBase,'Landsat-8')
if not os.path.exists(modisBase):
    os.mkdir(modisBase)
landsatSR = os.path.join(landsatDataBase,'SR')
if not os.path.exists(landsatSR):
    os.mkdir(landsatSR)
lstBase = os.path.join(landsatDataBase,'LST')
landsatTemp = os.path.join(landsatSR,'temp')
if not os.path.exists(landsatTemp):
    os.mkdir(landsatTemp)
processLST = os.path.join(base,'processData','LST')
if not os.path.exists(processLST):
    os.mkdir(processLST)


def geotiff2envi():   
    geotiffConvert = os.path.join(base,'processData','bin','GeoTiff2ENVI')
    bands = ["blue","green","red","nir","swir1","swir2","cloud"]
    l8bands = ["sr_band2","sr_band3","sr_band4","sr_band5","sr_band6","sr_band7","cfmask"] 
    
    landsatFiles = glob.glob(os.path.join(landsatTemp,"*.xml"))
    for i in range(len(landsatFiles)):
        fstem = landsatFiles[i][:-4]
        for i in range(len(bands)):
            tifFile = fstem+"_%s.tif" % l8bands[i]
            datFile = fstem+"_%s.%s.dat" % (l8bands[i],bands[i])            
            subprocess.call(["%s" % geotiffConvert ,"%s" % tifFile, "%s" % datFile])
        LSTtifFile = fstem+"_lst.tiff"
        LSTdatFile = fstem+"_lst.dat"
        subprocess.call(["%s" % geotiffConvert ,"%s" % LSTtifFile, "%s" % LSTdatFile])
          
bands = ["blue","green","red","nir","swir1","swir2","cloud"]
l8bands = ["sr_band2","sr_band3","sr_band4","sr_band5","sr_band6","sr_band7","cfmask"] 
fn = os.path.join(landsatSR,'028032','LC80280322015224LGN00_sr_band1.tif')
ls = landsatTools.GeoTIFF(fn)
sceneID = fn.split(os.sep)[-1][:-13]
# extract the Landsat doy and year
ldoy = sceneID[13:16]
year = int(sceneID[9:13])
# convert to date    
dd = datetime.datetime(year, 1, 1) + datetime.timedelta(int(ldoy) - 1)
date = '%d-%02d-%02d' % (dd.year,dd.month,dd.day)
if int(sceneID[2])==5:
    th_res="120.0"
elif int(sceneID[2])==7:
    th_res="60.0"
elif int(sceneID[2])==8:
    th_res="90.0"

nrows = ls.nrow
ncols = ls.ncol

res = ls.delx
ul = "%f,%f" %(ls.ulx,ls.uly)
#ul = "%f,%f" %(ls.uly,ls.ulx )
zone = ls.proj4.split(' ')[1].split('=')[-1]
sbands = ""
for i in range(len(bands)-1):
#    sbands = sbands + os.path.join(landsatTemp,"%s_%s.%s.tif  " % (sceneID,l8bands[i],bands[i]))
#cmask = os.path.join(landsatTemp,"%s_%s.%s.tif" % (sceneID,l8bands[-1],bands[-1]))
    sbands = sbands + os.path.join(landsatTemp,"%s_%s.tif " % (sceneID,l8bands[i]))
cmask = os.path.join(landsatTemp,"%s_%s.tif" % (sceneID,l8bands[-1]))
inf = os.path.join(lstBase,"%s_lst.tiff" % sceneID)
outf = os.path.join(landsatTemp,"lndsr.%s.sharpened_band6.bin" % sceneID)

#====create dms.inp ====================================
fn = os.path.join(processLST,"dms.inp")
file = open(fn, "w")
file.write("# input file for Data Mining Sharpener\n")
file.write("NFILES = 6\n")
file.write("SW_FILE_NAME = %s\n" % sbands)
file.write("SW_CLOUD_MASK = %s\n" % cmask)
file.write("SW_FILE_TYPE = GeoTiff\n")
file.write("SW_CLOUD_TYPE = GeoTiff\n")
file.write("SW_NROWS = %d\n" % nrows)
file.write("SW_NCOLS = %d\n" % ncols)
file.write("SW_PIXEL_SIZE = %f\n" % res)
file.write("SW_FILL_VALUE = -9999\n")
file.write("SW_CLOUD_CODE = 1\n")
file.write("SW_DATA_RANGE = -2000, 16000\n")
file.write("SW_UPPER_LEFT_CORNER = %s\n" % ul)
file.write("SW_PROJECTION_CODE = 1\n")
file.write("SW_PROJECTION_PARAMETERS = 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n")
file.write("SW_PROJECTION_ZONE = %s\n" % zone)
file.write("SW_PROJECTION_UNIT = 1\n")
file.write("SW_PROJECTION_DATUM = 12\n")
 
file.write("ORG_TH_FILE_NAME = %s\n" % inf)
file.write("ORG_TH_FILE_TYPE = GeoTiff\n")
file.write("ORG_TH_DATA_RANGE = 230., 370.\n")
file.write("ORG_TH_PIXEL_SIZE = %f\n" % res)
file.write("ORG_NROWS = %d\n" % nrows)
file.write("ORG_NCOLS = %d\n" % ncols)

file.write("RES_TH_PIXEL_SIZE = %s\n" % th_res)

file.write("PURE_CV_TH = 0.1\n")
file.write("ZONE_SIZE = 240\n")
file.write("SMOOTH_FLAG = 1\n")
file.write("CUBIST_FILE_STEM = th_samples\n")
file.write("OUT_FILE = %s\n" % outf)
file.write("end")

file.close()