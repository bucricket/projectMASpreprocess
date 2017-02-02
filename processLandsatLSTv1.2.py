# -*- coding: utf-8 -*-
"""
Created on Wed May 18 10:08:01 2016

@author: schull
"""

import os
import sys
import pyrttov
import numpy as np
from osgeo import gdal
from osgeo import osr
import glob
import subprocess
import landsatTools
import pyDisALEXI_utils as dA
import h5py
import urllib
from pydap.client import open_url
from pydap.client import open_dods
from landsat.downloader import Downloader
from landsat.search import Search
import datetime
import shutil
import tarfile
import pycurl
from pydap.cas.urs import setup_session
session = setup_session('mschull', 'sushmaMITCH12')


base = os.path.join(os.sep,'data','smcd4','mschull')
base = os.path.join(os.sep,'Users','mschull','umdGD')


def setFolders(base):
    dataBase = os.path.join(base,'pyDisALEXI','data')
    landsatDataBase = os.path.join(dataBase,'Landsat-8')
    lstBase = os.path.join(landsatDataBase,'LST')
    asterDataBase = os.path.join(dataBase,'ASTER')
    if not os.path.exists(lstBase):
        os.makedirs(lstBase)
    landsatDN = os.path.join(landsatDataBase,'DN')
    if not os.path.exists(landsatDN):
        os.makedirs(landsatDN)
    landsatProcessed = os.path.join(landsatDataBase,'processed')
    if not os.path.exists(landsatProcessed):
        os.makedirs(landsatProcessed)
    mosaicTemp = os.path.join(landsatProcessed,'mosaic')
    if not os.path.exists(mosaicTemp):
        os.makedirs(mosaicTemp)
    asterEmissivityBase = os.path.join(asterDataBase,'asterEmissivity')
    if not os.path.exists(asterEmissivityBase):
        os.makedirs(asterEmissivityBase)
    landsatEmissivityBase = os.path.join(asterDataBase,'landsatEmissivity')
    if not os.path.exists(landsatEmissivityBase):
        os.makedirs(landsatEmissivityBase)
    ASTERmosaicTemp = os.path.join(asterDataBase,'mosaicTemp')    
    if not os.path.exists(ASTERmosaicTemp):
        os.makedirs(ASTERmosaicTemp)
    out = {'dataBase':dataBase,'landsatDataBase':landsatDataBase,\
    'landsatDN':landsatDN,'landsatProcessed':landsatProcessed,\
    'mosaicTemp':mosaicTemp,'lstBase':lstBase,'ASTERmosaicTemp':ASTERmosaicTemp,\
    'asterEmissivityBase':asterEmissivityBase,'landsatEmissivityBase':landsatEmissivityBase}
    return out
    
Folders = setFolders(base)    
landsatDataBase = Folders['landsatDataBase']
landsatDN = Folders['landsatDN']
landsatProcessed=Folders['landsatProcessed']
mosaicTemp = Folders['mosaicTemp']
lstBase = Folders['lstBase']
ASTERmosaicTemp = Folders['ASTERmosaicTemp']
asterEmissivityBase = Folders['asterEmissivityBase']
landsatEmissivityBase = Folders['landsatEmissivityBase']
rttov_installdir = os.path.join(base,'pyDisALEXI','processData','rttov113')

def writeArray2Tiff(data,lats,lons,outfile):
    Projection = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    xres = lons[1] - lons[0]
    yres = lats[1] - lats[0]

    ysize = len(lats)
    xsize = len(lons)

    ulx = lons[0] #- (xres / 2.)
    uly = lats[0]# - (yres / 2.)
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(outfile, xsize, ysize, 1, gdal.GDT_Float32)
    
    srs = osr.SpatialReference()
    if isinstance(Projection, basestring):        
        srs.ImportFromProj4(Projection)
    else:
        srs.ImportFromEPSG(Projection)        
    ds.SetProjection(srs.ExportToWkt())
    
    gt = [ulx, xres, 0, uly, 0, yres ]
    ds.SetGeoTransform(gt)
    
    outband = ds.GetRasterBand(1)
    outband.WriteArray(data)    
    ds.FlushCache()  
    
    ds = None
def writeImageData(dataset,geo,proj,shape,fileFormat,fName,ddtype):    
    driver = gdal.GetDriverByName(fileFormat)
    dst_ds = driver.Create(fName, shape[1], shape [0],1, ddtype)
    dst_ds.SetGeoTransform(geo)
    dst_ds.SetProjection(proj)
    dst_ds.GetRasterBand(1).WriteArray(dataset)
    dst_ds = None
    return None 
    
def expand2nprofiles(n, nprof):
    # Transform 1D array to a [nprof, nlevels] array
    outp = np.empty((nprof, len(n)), dtype=n.dtype)
    for i in range(nprof):
        outp[i, :] = n[:]
    return outp
    
def untar(fname, fpath):
    if (fname.endswith('tar.gz') or fname.endswith('tar.bz')):
        tar = tarfile.open(fname)
        tar.extractall(path = fpath)
        tar.close()
        os.remove(fname)

def downloadEarthdata(fileUrl,outfile):
    status = 0
    c = pycurl.Curl()
    with open(outfile, 'w') as f:
        c.setopt(c.WRITEFUNCTION, f.write)
        c.setopt(c.NETRC,1)
        
        c.setopt(c.FOLLOWLOCATION,1)
        c.setopt(c.COOKIEJAR,'cookiefile')
        c.setopt(c.COOKIE,'cookiefile')
        c.setopt(c.FAILONERROR,True)
        try:

            c.setopt(c.URL,fileUrl)
            c.perform()
        except pycurl.error, error:
            errno, errstr = error
            print 'An error occurred: ', errstr
            os.remove(outfile)
            status = 1
    return status
                
        
def processASTERemis(rawLandsatFolder,landsatGTiff):  
    ls = landsatTools.GeoTIFF(landsatGTiff)
    #sc = raster.Landsatscene(rawLandsatFolder)
    meta = landsatTools.landsat_metadata(landsatGTiff)
    ullat = meta.CORNER_UL_LAT_PRODUCT
    ullon = meta.CORNER_UL_LON_PRODUCT
    inProj4 = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

    ASTERurlBase = 'http://e4ftl01.cr.usgs.gov/ASTT/AG100.003/2000.01.01'
    # use Landsat scene area
    UL = [int(np.ceil(ullat)),int(np.floor(ullon))]
    
    for i in xrange(4):
        for j in xrange(5):
            status = 0
            if UL[1]>0:
                ulLon=str((UL[1]+j)).zfill(3)
            else:
                ulLon=str((UL[1]+j)).zfill(4)
            asterFN = 'AG100.v003.%d.%s.0001.h5' % (UL[0]-i,ulLon)
            # ASTER Emissivity product AG100 comes in 1 x 1 degree tiles where UL is given in the filename.
            ASTERurl = os.path.join(ASTERurlBase,asterFN)
            print ASTERurl
            localAsterFN = os.path.join(asterEmissivityBase,asterFN)
            if not os.path.isfile(localAsterFN):
#                #first download HDF file
#                req = Request(ASTERurl)
#                try:
#                    response = urlopen(req)
#                except URLError as e:
#                    if hasattr(e, 'reason'):
#                        print('We failed to reach a server.')
#                        print('Reason: ', e.reason)
#                    elif hasattr(e, 'code'):
#                        print('The server couldn\'t fulfill the request.')
#                        print('Error code: ', e.code)
#                else:
                print "downloading ASTER..."
                #urllib.urlretrieve(ASTERurl,localAsterFN) # THIS NO LONGER WORKS BECAUSE NASA REQUIRES LOGIN
                status = downloadEarthdata(ASTERurl,localAsterFN)
       
        
        
        #open HDF file, extract the desired dataset and save to GTiff
            
            if status == 0:
                fh5  = h5py.File(localAsterFN , 'r')
                EmisBand4 = np.array(fh5["/Emissivity/Mean/"])[3]/1000.
                lats = np.array(fh5["/Geolocation/Latitude/"])
                lons = np.array(fh5["/Geolocation/Longitude/"])
                inUL = [lats[0,0],lons[0,0]]
                inRes = [90.0,90.0]
                tempName = os.path.join(ASTERmosaicTemp,'emis%d%s.tiff'% (UL[0]-i,ulLon))
                writeArray2Tiff(EmisBand4,lats[:,0],lons[0,:],tempName)
                outFormat = gdal.GDT_Float32
                #dA.writeArray2Tiff(EmisBand4,inRes,inUL,inProj4,tempName,outFormat)
    
    
    #mosaic ,reproject and save as geotiff
    mosaicTempFN = '%s/mosaic.vrt' % ASTERmosaicTemp
    mosaicVRTcommand = 'gdalbuildvrt -srcnodata 0 %s %s/*.tiff' % (mosaicTempFN,ASTERmosaicTemp)
    out = subprocess.check_output(mosaicVRTcommand, shell=True)
    sceneID = meta.LANDSAT_SCENE_ID
    resampName = os.path.join(landsatEmissivityBase,'%s_EMIS.tiff' % sceneID)
    #command = "gdalwarp -overwrite -s_srs EPSG:%d -t_srs '%s' -r bilinear -tr 30 30 -te %f %f %f %f -te_srs EPSG:%d -of GTiff %s %s" % (4326,ls.proj4,ullon,lrlat,lrlon,ullat,4326, mosaicTempFN,resampName)
    command = "gdalwarp -overwrite -s_srs '%s' -t_srs '%s' -r bilinear -tr 90 90 -te %f %f %f %f -of GTiff %s %s" % (inProj4,ls.proj4,ls.ulx,ls.lry,ls.lrx,ls.uly, mosaicTempFN,resampName)
    out = subprocess.check_output(command, shell=True)
    print 'done processing ASTER'
    shutil.rmtree(ASTERmosaicTemp)
    os.makedirs(ASTERmosaicTemp)
    return resampName
        
def km2deg(x,y,lat):
    
    # how many KM in 1 deg
    degLat = 110.54 # KM
    degLon = 111.320*np.cos(np.deg2rad(lat))  #KM
    
    degOut = []
    degOut.append(x/degLat)
    degOut.append(y/degLon) 
    return degOut
    
def downloadLandsat(startDate,endDate,GCP):
    #find all path/rows associates with GCP
    s = Search()
    scenes = s.search(lat=GCP[0],lon=GCP[1],limit = 100, start_date = startDate,end_date=endDate, cloud_max=5)
    print '%d scenes found' % scenes['total_returned']
    if len(scenes)==3:
        print 'no cloud free scenes found'
        return
    sceneIDlist =[]
    for i in xrange(len(scenes['results'])):        
        sceneID = np.str(scenes['results'][i]['sceneID'])
        if sceneID.startswith('LO'): #OLI only
            continue
        #sceneIDlist.append(sceneID)
        filelist = glob.glob(os.path.join(landsatDN,sceneID))
        if not filelist:
            d = Downloader(download_dir = landsatDN)
            d.download([sceneID],['10'])
            path =  os.path.join(landsatDN,sceneID)
            if not os.path.exists(path):
                tarFN =  os.path.join(landsatDN,'%s.tar.bz' % sceneID)
                if os.path.exists(tarFN):
                    try:            
                        if not os.path.exists(path):
                            os.makedirs(path)
                        untar(tarFN, path)
                        sceneIDlist.append(sceneID)
                        continue
                    except Exception:
                        print 'something wrong with tar file'
                        os.remove(tarFN)
                        shutil.rmtree(path)
                        continue
                elif os.path.exists(os.path.join(landsatDN,'%s.tar.gz' % sceneID)):
                    tarFN = os.path.join(landsatDN,'%s.tar.gz' % sceneID)
                    try:
                        if not os.path.exists(path):
                            os.makedirs(path)
                        untar(tarFN, path)
                        sceneIDlist.append(sceneID)
                        continue
                    except Exception:
                        print 'something wrong with tar file'
                        os.remove(tarFN)
                        shutil.rmtree(path)
                        continue
                else:
                    print 'no files exist on AWS or google'
                    continue
            sceneIDlist.append(sceneID)
        else:
            sceneIDlist.append(sceneID)
    return sceneIDlist
    
def prepareMERRA2data(landsatSceneID):
    print "using landsat scene: %s" % landsatSceneID
    mtlFile = os.path.join(landsatDN,landsatSceneID,'%s_MTL.txt' % landsatSceneID)
    meta = landsatTools.landsat_metadata(mtlFile)
    ulLat = meta.CORNER_UL_LAT_PRODUCT
    ulLon = meta.CORNER_UL_LON_PRODUCT
    lrLat = meta.CORNER_LR_LAT_PRODUCT
    lrLon = meta.CORNER_LR_LON_PRODUCT
    solZen = meta.SUN_ELEVATION
    solAzi = meta.SUN_AZIMUTH
    landsatDate = meta.DATE_ACQUIRED
    landsatTime = meta.SCENE_CENTER_TIME
    
    d = datetime.datetime.strptime('%s%s' % (landsatDate,landsatTime),'%Y-%m-%d%H:%M:%S.%f')
    
    year = d.year
    month = d.month
    day = d.day
    hr = d.hour
    
    ul = [ulLon-1.5,ulLat+1.5]
    lr = [lrLon+1.5,lrLat-1.5]
    # The data is lat/lon and upside down so [0,0] = [-90.0,-180.0]
    maxX = int((lr[0]-(-180))/0.625)
    minX = int((ul[0]-(-180))/0.625)
    minY = int((lr[1]-(-90))/0.5)
    maxY = int((ul[1]-(-90))/0.5)
    
    
    if year <1992:
        fileType = 100
    elif year >1991 and year < 2001:
        fileType=200
    elif year > 2000 and year<2011:
        fileType = 300
    else:
        fileType = 400
    
    #Instantaneous Two-Dimensional Collections
    #inst1_2d_asm_Nx (M2I1NXASM): Single-Level Diagnostics
    #=============================================================================
    opendap_url = 'http://goldsmr4.sci.gsfc.nasa.gov:80/opendap/MERRA2/'
    product = 'M2I1NXASM.5.12.4'
    filename = 'MERRA2_%d.inst1_2d_asm_Nx.%04d%02d%02d.nc4' % (fileType,year,month,day)
    fullUrl =os.path.join(opendap_url,product,'%04d'% year,'%02d'% month,filename)
    #d=open_dods(fullUrl+'?PS[1:1:23][0:1:360][0:1:575]')
    d = open_url(fullUrl, session=session)
#    d.keys()
    #surface presure [Pa]
    surfacePressure=d.PS
    sp = np.squeeze(surfacePressure[hr,minY:maxY,minX:maxX]/100) # Pa to kPa
    sprshp =np.reshape(sp,sp.shape[0]*sp.shape[1])
    
    #2m air Temp (K)
    Temp2 = d.T2M
    #Temp2=open_dods(fullUrl+'?T2M[1:1:23][0:1:360][0:1:575]')
    t2 = np.squeeze(Temp2[hr,minY:maxY,minX:maxX])
    t2rshp =np.reshape(t2,t2.shape[0]*t2.shape[1])
    
    #2m specific humidity [kg kg -1] -> 2 m water vapor [ppmv]
    spcHum = d.QV2M
    #spcHum=open_dods(fullUrl+'?QV2M[1:1:23][0:1:360][0:1:575]')
    # wv_mmr = 1.e-6 * wv_ppmv_layer * (Rair / Rwater)
    # wv_mmr in kg/kg, Rair = 287.0, Rwater = 461.5
    q = np.squeeze(spcHum[hr,minY:maxY,minX:maxX])
    q2 = q/(1e-6*(287.0/461.5))
    q2rshp =np.reshape(q2,q2.shape[0]*q2.shape[1])
    
    # skin temp [K]
    sktIn = d.TS
    #sktIn=open_dods(fullUrl+'?TS[1:1:23][0:1:360][0:1:575]')
    skt = np.squeeze(sktIn[hr,minY:maxY,minX:maxX])
    sktrshp =np.reshape(skt,skt.shape[0]*skt.shape[1])
    
    # U10M 10-meter_eastward_wind [m s-1]
    u10In = d.U10M
    #u10In=open_dods(fullUrl+'?U10[1:1:23][0:1:360][0:1:575]')
    u10 = np.squeeze(u10In[hr,minY:maxY,minX:maxX])
    u10rshp =np.reshape(u10,u10.shape[0]*u10.shape[1])
    
    # V10M 10-meter_northward_wind [m s-1]
    v10In = d.V10M
    #v10In=open_dods(fullUrl+'?V10M[1:1:23][0:1:360][0:1:575]')
    v10 = np.squeeze(v10In[hr,minY:maxY,minX:maxX])
    v10rshp =np.reshape(v10,v10.shape[0]*v10.shape[1])
    
    opendap_url = 'http://goldsmr5.sci.gsfc.nasa.gov:80/opendap/MERRA2/'
    product = 'M2I3NVASM.5.12.4'
    filename = 'MERRA2_%d.inst3_3d_asm_Nv.%04d%02d%02d.nc4' % (fileType,year,month,day)
    fullUrl =os.path.join(opendap_url,product,'%04d'% year,'%02d'% month,filename)
    d = open_url(fullUrl,session=session)
    hr = int(np.round(hr/3.)) # convert from 1 hr to 3 hr dataset
    
    #layers specific humidity [kg kg -1] -> 2 m water vapor [ppmv]
    qvIn = d.QV
    #qvIn=open_dods(fullUrl+'?QV[0:1:7][0,:1:71][0:1:360][0:1:575]')
    # wv_mmr = 1.e-6 * wv_ppmv_layer * (Rair / Rwater)
    # wv_mmr in kg/kg, Rair = 287.0, Rwater = 461.5
    qv = np.squeeze(qvIn[hr,:,minY:maxY,minX:maxX])
    qv = qv/(1e-6*(287.0/461.5))
    qvrshp =np.reshape(qv,[qv.shape[0],qv.shape[1]*qv.shape[2]]).T
    
    
    #layers air temperature [K]
    tIn = d.T
    #tIn=open_dods(fullUrl+'?T[0:1:7][0,:1:71][0:1:360][0:1:575]')
    # wv_mmr = 1.e-6 * wv_ppmv_layer * (Rair / Rwater)
    # wv_mmr in kg/kg, Rair = 287.0, Rwater = 461.5
    t = np.squeeze(tIn[hr,:,minY:maxY,minX:maxX])
    trshp =np.reshape(t,[t.shape[0],t.shape[1]*t.shape[2]]).T
    
    #mid_level_pressure [Pa]
    
    plIn=d.PL
    #plIn=open_dods(fullUrl+'?PL[0:1:7][0,:1:71][0:1:360][0:1:575]')
    pl = np.squeeze(plIn[hr,:,minY:maxY,minX:maxX]/100) # Pa to kPa
    plrshp =np.reshape(pl,[pl.shape[0],pl.shape[1]*pl.shape[2]]).T
    #qrshp =np.reshape(q,q.shape[0]*q.shape[1])
    
    
    LAT = d.lat
    LON = d.lon
    lats = LAT[:]
    lons = LON[:]
    lat = np.tile(lats,(len(lons),1)).T
    latIn = np.squeeze(lat[minY:maxY,minX:maxX])
    latrshp =np.reshape(latIn,latIn.shape[0]*latIn.shape[1])
    lon = np.tile(lons,(len(lats),1))
    lonIn = np.squeeze(lon[minY:maxY,minX:maxX])
    lonrshp =np.reshape(lonIn,lonIn.shape[0]*lonIn.shape[1])
    el = np.repeat(0.0,v10.shape[0]*v10.shape[1]) #NEED DEM
    #check surface pressure
    
    
    sunzen = np.repeat(solZen,v10.shape[0]*v10.shape[1])
    sunazi = np.repeat(solAzi,v10.shape[0]*v10.shape[1])
    fetch = np.repeat(100000,v10.shape[0]*v10.shape[1])
    satzen = np.repeat(0.0,v10.shape[0]*v10.shape[1])
    satazi = np.repeat(0.0,v10.shape[0]*v10.shape[1])
    
    # Units for gas profiles
    gas_units = 2  # ppmv over moist air
    
    # datetimes[6][nprofiles]: yy, mm, dd, hh, mm, ss
    datetimes = np.tile([year, month, day, hr, 0, 0],(v10.shape[0]*v10.shape[1],1))
    
    # angles[4][nprofiles]: satzen, satazi, sunzen, sunazi
    #get from landsat MTL
    angles = np.vstack((satzen,satazi,sunzen,sunazi)).T
    
    # surftype[2][nprofiles]: surftype, watertype
    surftype = np.zeros([angles.shape[0],2]) #NEED LAND/WATER mask
    
    # surfgeom[3][nprofiles]: lat, lon, elev
    surfgeom = np.vstack((latrshp,lonrshp,el)).T
    
    # s2m[6][nprofiles]: 2m p, 2m t, 2m q, 10m wind u, v, wind fetch
    s2m = np.vstack((sprshp,t2rshp,q2rshp,u10rshp,v10rshp,fetch)).T
    
    # skin[9][nprofiles]: skin T, salinity, snow_frac, foam_frac, fastem_coefsx5
    sal = np.repeat(35.0,v10.shape[0]*v10.shape[1])
    snow_frac = np.repeat(0.0,v10.shape[0]*v10.shape[1])
    foam_frac= np.repeat(0.0,v10.shape[0]*v10.shape[1])
    fastem_coef1 = np.repeat(3.0,v10.shape[0]*v10.shape[1])
    fastem_coef2 = np.repeat(5.0,v10.shape[0]*v10.shape[1])
    fastem_coef3 = np.repeat(15.0,v10.shape[0]*v10.shape[1])
    fastem_coef4 = np.repeat(0.1,v10.shape[0]*v10.shape[1])
    fastem_coef5 = np.repeat(0.3,v10.shape[0]*v10.shape[1])
    
    skin= np.vstack((sktrshp,sal,snow_frac,foam_frac,fastem_coef1,fastem_coef2,fastem_coef3,fastem_coef4,fastem_coef5)).T

    outDict = {'P':plrshp,'T':trshp,'Q':qvrshp,'Angles':angles,'S2m':s2m,\
    'Skin': skin,'SurfType':surftype,'SurfGeom':surfgeom,'Datetimes':datetimes,\
    'origShape':t2.shape}
    
    return outDict
    
def runRTTOV(profileDict):
    nlevels = profileDict['P'].shape[1]
    nprofiles = profileDict['P'].shape[0]
    myProfiles = pyrttov.Profiles(nprofiles, nlevels)
    myProfiles.GasUnits = 2
    myProfiles.P = profileDict['P']
    myProfiles.T = profileDict['T']
    myProfiles.Q = profileDict['Q']
    myProfiles.Angles = profileDict['Angles']
    myProfiles.S2m = profileDict['S2m']
    myProfiles.Skin = profileDict['Skin']
    myProfiles.SurfType = profileDict['SurfType']
    myProfiles.SurfGeom =profileDict['SurfGeom']
    myProfiles.DateTimes = profileDict['Datetimes']
    month = profileDict['Datetimes'][0,1]

    # ------------------------------------------------------------------------
    # Set up Rttov instance
    # ------------------------------------------------------------------------

    # Create Rttov object for the TIRS instrument

    tirsRttov = pyrttov.Rttov()
    nchan_tirs = 1

    # Set the options for each Rttov instance:
    # - the path to the coefficient file must always be specified
    # - specify paths to the emissivity and BRDF atlas data in order to use
    #   the atlases (the BRDF atlas is only used for VIS/NIR channels so here
    #   it is unnecessary for HIRS or MHS)
    # - turn RTTOV interpolation on (because input pressure levels differ from
    #   coefficient file levels)
    # - set the verbose_wrapper flag to true so the wrapper provides more
    #   information
    # - enable solar simulations for SEVIRI
    # - enable CO2 simulations for HIRS (the CO2 profiles are ignored for
    #   the SEVIRI and MHS simulations)
    # - enable the store_trans wrapper option for MHS to provide access to
    #   RTTOV transmission structure

    
    tirsRttov.FileCoef = '{}/{}'.format(rttov_installdir,
                                       "rtcoef_rttov11/rttov7pred54L/rtcoef_landsat_8_tirs.dat")
    
    #tirsRttov.EmisAtlasPath = os.path.join(base,'ALEXIdisALEXIfusion','rttov113','emis_data')
    tirsRttov.EmisAtlasPath = '{}/{}'.format(rttov_installdir, "emis_data")
    print "%s" % tirsRttov.EmisAtlasPath
    tirsRttov.BrdfAtlasPath = '{}/{}'.format(rttov_installdir, "brdf_data")
    #tirsRttov.BrdfAtlasPath = os.path.join(base,'ALEXIdisALEXIfusion','rttov113','brdf_data')

    tirsRttov.Options.AddInterp = True
    tirsRttov.Options.StoreTrans = True
    tirsRttov.Options.StoreRad2 = True
    tirsRttov.Options.VerboseWrapper = True


    # Load the instruments:

    try:
        tirsRttov.loadInst()
    except pyrttov.RttovError as e:
        sys.stderr.write("Error loading instrument(s): {!s}".format(e))
        sys.exit(1)

    # Associate the profiles with each Rttov instance
    tirsRttov.Profiles = myProfiles
    # ------------------------------------------------------------------------
    # Load the emissivity and BRDF atlases
    # ------------------------------------------------------------------------

    # Load the emissivity and BRDF atlases:
    # - load data for August (month=8)
    # - note that we only need to load the IR emissivity once and it is
    #   available for both SEVIRI and HIRS: we could use either the seviriRttov
    #   or hirsRttov object to do this
    # - for the BRDF atlas, since SEVIRI is the only VIS/NIR instrument we can
    #   use the single-instrument initialisation

    tirsRttov.irEmisAtlasSetup(month)
    # ------------------------------------------------------------------------
    # Call RTTOV
    # ------------------------------------------------------------------------

    # Since we want the emissivity/reflectance to be calculated, the
    # SurfEmisRefl attribute of the Rttov objects are left uninitialised:
    # That way they will be automatically initialise to -1 by the wrapper

    # Call the RTTOV direct model for each instrument:
    # no arguments are supplied to runDirect so all loaded channels are
    # simulated
    try:
        tirsRttov.runDirect()
    except pyrttov.RttovError as e:
        sys.stderr.write("Error running RTTOV direct model: {!s}".format(e))
        sys.exit(1)
        
    return tirsRttov
def processLandsatLST(tirsRttov,landsatscene,merraDict):
    
    origShap = merraDict['origShape']
    surfgeom=merraDict['SurfGeom']
    nlevels = merraDict['P'].shape[1]
    MERRA2_ulLat = 90.0
    MERRA2_ulLon = -180.0
    MERRA2LatRes = 0.5
    MERRA2LonRes = 0.625
    #reshape and resize image to fit landsat
    lats = np.flipud(np.resize(surfgeom[:,0],origShap))
    lons = np.flipud(np.resize(surfgeom[:,1],origShap))
    inRes = [MERRA2LonRes,MERRA2LatRes]
    inUL = [MERRA2_ulLat,MERRA2_ulLon]
    channel=1
    rawLandsatFolder =os.path.join(landsatDN, landsatscene)
    mtlFile = os.path.join(landsatDN,landsatscene,'%s_MTL.txt' % landsatscene)
    meta = landsatTools.landsat_metadata(mtlFile)
    ulx = meta.CORNER_UL_PROJECTION_X_PRODUCT
    uly = meta.CORNER_UL_PROJECTION_Y_PRODUCT
    lrx = meta.CORNER_LR_PROJECTION_X_PRODUCT
    lry = meta.CORNER_LR_PROJECTION_Y_PRODUCT
    dn2Radslope = meta.RADIANCE_MULT_BAND_10
    dn2Radintcpt = meta.RADIANCE_ADD_BAND_10
    #chipName = pdap.landsatscene[3:9]
    #landsat = glob.glob(os.path.join(pdap.landsatProcessed,chipName,'*%s*RAD.tif' % chipName))[0]
    landsat = glob.glob(os.path.join(rawLandsatFolder,'*B10.TIF'))[0]
    #convert from 30 to 90 m
    resampName = os.path.join('%sResample.vrt' % landsat[:-4])
    command = "gdalwarp -overwrite -r average -tr 90 90 -of VRT %s %s" % (landsat,resampName)
    out = subprocess.check_output(command, shell=True)
    Lg = gdal.Open(resampName)
    Ls = landsatTools.GeoTIFF(resampName)
    ulx = Ls.ulx
    uly = Ls.uly
    lrx = Ls.lrx
    lry = Ls.lry
    L8therm = Lg.ReadAsArray()
    ThermalRad = L8therm*dn2Radslope+dn2Radintcpt
    print "L8 Size: %f,%f" % (ThermalRad.shape[0],ThermalRad.shape[1])
    geo = Lg.GetGeoTransform()
    proj = Lg.GetProjection()
    Lsrs=osr.SpatialReference()
    Lsrs.ImportFromWkt(proj)
    landsatProj4 = str(Lsrs.ExportToProj4())
    landsatProj4 = Ls.coordtrans.srs
    inProj4 = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    Lg = None
    

    nu4 = 1/(10.895*0.0001) # convert to cm
    
    #Process downwelling radiance
    RadDown = np.flipud(np.resize(tirsRttov.Rad2Down[:,channel,nlevels-2],origShap))
    tempName = os.path.join(landsatDataBase,'RadDown.tiff')
    resampName = os.path.join('%sReproj.tiff' % tempName[:-4])
    writeArray2Tiff(RadDown,lats[:,0],lons[0,:],tempName)
    #outFormat = gdal.GDT_Float32
    #dA.writeArray2Tiff(RadDown,inRes,inUL,inProj4,tempName,outFormat)
    command = "gdalwarp -overwrite -s_srs '%s' -t_srs '%s' -r bilinear -tr 90 90 -te %d %d %d %d -of GTiff %s %s" % (inProj4,landsatProj4,ulx,lry,lrx,uly,tempName,resampName)
    out = subprocess.check_output(command, shell=True)
    Lg = gdal.Open(resampName)
    RadDown = Lg.ReadAsArray()
    RadDown = (RadDown*(nu4**2/10**7))#*.001
    print "RadDown Size: %f,%f" % (RadDown.shape[0],RadDown.shape[1])
    Lg = None
    
    #Process upwelling radiance
    RadUp = np.flipud(np.resize(tirsRttov.Rad2Up[:,channel,nlevels-2],origShap))
    tempName = os.path.join(landsatDataBase,'RadUp.tiff')
    resampName = os.path.join('%sReproj.tiff' % tempName[:-4])
    #outFormat = gdal.GDT_Float32
    #dA.writeArray2Tiff(RadUp,inRes,inUL,inProj4,tempName,outFormat)
    writeArray2Tiff(RadUp,lats[:,0],lons[0,:],tempName)
    command = "gdalwarp -overwrite -s_srs '%s' -t_srs '%s' -r bilinear -tr 90 90 -te %d %d %d %d -of GTiff %s %s" % (inProj4,landsatProj4,ulx,lry,lrx,uly,tempName,resampName)
    out = subprocess.check_output(command, shell=True)
    Lg = gdal.Open(resampName)
    RadUp = Lg.ReadAsArray()
    RadUp = (RadUp*(nu4**2/10**7))#*.001
    Lg = None
    
    #Process transmission
    trans = np.flipud(np.resize(tirsRttov.TauTotal[:,channel],origShap))
    tempName = os.path.join(landsatDataBase,'trans.tiff')
    resampName = os.path.join('%sReproj.tiff' % tempName[:-4])
    writeArray2Tiff(trans,lats[:,0],lons[0,:],tempName)
    #outFormat = gdal.GDT_Float32
    #dA.writeArray2Tiff(trans,inRes,inUL,inProj4,tempName,outFormat)
    command = "gdalwarp -overwrite -s_srs '%s' -t_srs '%s' -r bilinear -tr 90 90 -te %d %d %d %d -of GTiff %s %s" % (inProj4,landsatProj4,ulx,lry,lrx,uly,tempName,resampName)
    out = subprocess.check_output(command, shell=True)
    Lg = gdal.Open(resampName)
    trans = Lg.ReadAsArray()
    Lg = None
      
    #get emissivity from ASTER
    path_row = landsatscene
    if not os.path.exists(os.path.join(landsatEmissivityBase,'%s_EMIS.tiff' % path_row)):    
        ASTERemisFN = processASTERemis(rawLandsatFolder,landsat)
    else:
        ASTERemisFN = os.path.join(landsatEmissivityBase,'%s_EMIS.tiff' % path_row)
    aster = gdal.Open(ASTERemisFN)
    emis = aster.ReadAsArray()
    print "emis Size: %f,%f" % (emis.shape[0],emis.shape[1])
    #emis = emis[:-1,:-1]
    aster = None
    # calcualte LST
    emis[emis<0.000001] = np.nan
    surfRad =(((ThermalRad-RadUp)/trans)-(1-emis)*RadDown)/emis
    #get Kappa constants from Landsat
    Kappa1 = meta.K1_CONSTANT_BAND_10
    Kappa2 = meta.K2_CONSTANT_BAND_10
    LST = Kappa2*(1/np.log(Kappa1/surfRad))
    lstName = os.path.join(lstBase,'%s_lst.tiff'% landsatscene)
    #write LST to a geoTiff
    writeImageData(LST,geo,proj,LST.shape,'GTiff',lstName,gdal.GDT_Float32)
    
    print 'done processing LST'
    
if __name__ == '__main__':
        # This example program simulates two profiles for the Landsat 8 TIRS instrument
    # The example profile data are defined in example_data
    startDate = '2015-01-01'
    endDate = '2015-12-31'
    GCP = [41.18,-96.00] #mead,NE
    #GCP = [30.8,32.2] #Nile Delta
    
    sceneIDlist = downloadLandsat(startDate,endDate,GCP)
    #save list to text file for ordering SR data for LAI calculation (THIS MAY CHANGE WHEN THE USGS RESTful SERVICE COMES ONLINE)
    outTextFn = os.path.join(landsatDataBase,'lansatScenes_%s_%d_%d.txt' % (startDate+'_'+endDate,GCP[0],GCP[1]))
    with open(outTextFn,'w') as file:
        for item in sceneIDlist:
            file.write("%s\n" % item)
            #print>>file, item

    # ------------------------------------------------------------------------
    # Set up the profile data
    # ------------------------------------------------------------------------
    for i in xrange(len(sceneIDlist)):
        landsatSceneID = sceneIDlist[i]
        tifFile = os.path.join(lstBase,'%s_lst.tiff' % landsatSceneID)
        binFile = os.path.join(lstBase,"lndsr."+landsatSceneID+".cband6.bin")
        
        if not os.path.exists(tifFile):
            print 'processing %s...' %landsatSceneID
            profileDict = prepareMERRA2data(landsatSceneID)
            tiirsRttov = runRTTOV(profileDict)
            processLandsatLST(tiirsRttov,landsatSceneID,profileDict)
        else:
            print '%s: already processed' %landsatSceneID
        
            
        subprocess.call(["gdal_translate","-of", "ENVI", "%s" % tifFile, "%s" % binFile])

        

    