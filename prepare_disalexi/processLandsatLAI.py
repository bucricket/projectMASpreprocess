#!/usr/bin/env python2
#!/Applications/anaconda/envs/root3 python
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 13:50:04 2016

@author: mschull
"""
#python

from landsat.search import Search
import os
import subprocess
import glob
import shutil
import pandas as pd
import datetime
import argparse
import getpass
import keyring
from pyproj import Proj
from utils import folders
from Clients import Client
from Order import Order
from OrderTemplate import OrderTemplate
import pycurl

base = os.getcwd()
Folders = folders(base)   
modisBase = Folders['modisBase']
landsatSR = Folders['landsatSR']
landsatLAI = Folders['landsatLAI']
landsatTemp = os.path.join(landsatSR,'temp')
if not os.path.exists(landsatTemp):
    os.mkdir(landsatTemp)

processData = os.path.join(base,'processData')
if not os.path.exists(processData ):
    os.mkdir(processData)

def getLandsatData(loc,startDate,endDate,auth):
    
    # build the various handlers to spec
    template = OrderTemplate('template')
    template.load(path='./L8SR_UTMorder.json' )
    order = Order(template, note="Lat%dLon%d-%s_%s" %(int(loc[0]),int(loc[1]),startDate,endDate))
    client = Client(auth) # will prompt user for username and password if auth argument not supplied
    #downloader = EspaLandsatLocalDownloader('USGS_downloads')
    
    # find cloud free landsat scenes
    s = Search()
    scenes = s.search(lat=loc[0],lon=loc[1],limit = 100, start_date = startDate,end_date=endDate, cloud_max=5)
    
    l8_tiles=[]
    for i in range(len(scenes['results'])):
        l8_tiles.append(scenes['results'][i]['sceneID'])
        
    # order the data
    order.add_tiles("olitirs8", l8_tiles)
    #order.add_tiles("etm7", l7_tiles)
    response = order.submit(client)
    
    # view the servers whole response. which might indicate an ordering error!
    print(response) 
    
    # assuming there were no order submission errors
    orderid = response['orderid']
    
    # now start the downloader!
    for download in client.download_order_gen(orderid):
        print(download)
    return l8_tiles
    # download is a tuple with the filepath, and True if the file
    # is a fresh download.
    
    # this is where data pipeline scripts go that can operate
    # on files as they are downloaded (generator),
    
    # See the Client class for further documentation.

def getMODISlai(tiles,product,version,startDate,endDate,auth):    

    subprocess.call(["modis_download.py", "-r", "-U", "%s" % auth[0], "-P", 
                    "%s" % auth[1],"-p", "%s.%s" % (product,version), "-t", 
                    "%s" % tiles,"-s","MOTA", "-f", "%s" % startDate,"-e", "%s" % endDate, 
                     "%s" % modisBase])
                     
def latlon2MODtile(lat,lon):
    # reference: https://code.env.duke.edu/projects/mget/wiki/SinusoidalMODIS
    p_modis_grid = Proj('+proj=sinu +R=6371007.181 +nadgrids=@null +wktext')
    x, y = p_modis_grid(lon, lat)
    # or the inverse, from x, y to lon, lat
    lon, lat = p_modis_grid(x, y, inverse=True)
    tileWidth = 1111950.5196666666
    ulx = -20015109.354
    uly = -10007554.677
    H = (x-ulx)/tileWidth
    V = 18-((y-uly)/tileWidth)
    return int(V),int(H)
    
def geotiff2envi():   
    #geotiffConvert = os.path.join(base,'processData','bin','GeoTiff2ENVI')
    geotiffConvert = 'GeoTiff2ENVI'
    bands = ["blue","green","red","nir","swir1","swir2","cloud"]
    l8bands = ["sr_band2","sr_band3","sr_band4","sr_band5","sr_band6","sr_band7","cfmask"] 
    
    landsatFiles = glob.glob(os.path.join(landsatTemp,"*.xml"))
    for i in range(len(landsatFiles)):
        fstem = landsatFiles[i][:-4]
        for i in range(len(bands)):
            tifFile = fstem+"_%s.tif" % l8bands[i]
            datFile = fstem+"_%s.%s.dat" % (l8bands[i],bands[i])
            subprocess.call(["%s" % geotiffConvert ,"%s" % tifFile, "%s" % datFile])

def sample():    
    #sample = os.path.join(base,'processData','bin','lndlai_sample')
    sample = 'lndlai_sample'
    bands = ["blue","green","red","nir","swir1","swir2","cloud"]
    l8bands = ["sr_band2","sr_band3","sr_band4","sr_band5","sr_band6","sr_band7","cfmask"] 
    
    landsatFiles = glob.glob(os.path.join(landsatTemp,"*.xml"))
    for i in range(len(landsatFiles)):
        sceneID = landsatFiles[i].split(os.sep)[-1][:-4]
        # extract the Landsat doy and year
        ldoy = sceneID[13:16]
        year = int(sceneID[9:13])
        # convert to date    
        dd = datetime.datetime(year, 1, 1) + datetime.timedelta(int(ldoy) - 1)
        date = '%d-%02d-%02d' % (dd.year,dd.month,dd.day)
        # find the 4 day MODIS doy prior to the Landsat doy
        mdoy = int((int((float(ldoy)-1.)/4.)*4.)+1)
        
        modFiles = glob.glob(os.path.join(modisBase,"MCD15A3.A%s%s.*.hdf" % (year,mdoy)))

        fstem = landsatFiles[i][:-4]
        laiPath = landsatLAI
        if not os.path.exists(laiPath):
            os.mkdir(laiPath)
        sam_file = os.path.join(laiPath,"SR_LAI.%s.%s.MCD15A3_A%s%s.txt" %(date,sceneID,year,mdoy))
        
        for i in range(len(modFiles)):  
            fn = os.path.join(laiPath,"slai%s.inp" % i)
            file = open(fn, "w")
            file.write("LANDSAT_BASE_BLUE = %s_%s.%s.dat\n" % (fstem,l8bands[0],bands[0]))
            file.write("LANDSAT_BASE_GREEN = %s_%s.%s.dat\n" % (fstem,l8bands[1],bands[1]))
            file.write("LANDSAT_BASE_RED = %s_%s.%s.dat\n" % (fstem,l8bands[2],bands[2]))
            file.write("LANDSAT_BASE_NIR = %s_%s.%s.dat\n" % (fstem,l8bands[3],bands[3]))
            file.write("LANDSAT_BASE_SWIR1 = %s_%s.%s.dat\n" % (fstem,l8bands[4],bands[4]))
            file.write("LANDSAT_BASE_SWIR2 = %s_%s.%s.dat\n" % (fstem,l8bands[5],bands[5]))
            file.write("LANDSAT_BASE_CLOUD = %s_%s.%s.dat\n" % (fstem,l8bands[6],bands[6]))
            file.write("MODIS_BASE_FILE = %s\n" % modFiles[i])
            file.write("SAMPLE_FILE_OUT = %s\n" % sam_file)
            file.write("PURE_SAMPLE_TH = 0.2\n")
            file.close()
        
            subprocess.call(["%s" % sample , "%s" % fn])
            os.remove(os.path.join(laiPath,"slai%s.inp" % i))
            
def train():    
    cubist = os.path.join(base,'processData','bin','cubist')
    
    landsatFiles = glob.glob(os.path.join(landsatLAI,"*.txt"))
    #======combine input data======================================
    df = pd.DataFrame(columns=['ulx','uly','blue',
        'green','red','nir','swir1','swir2','ndvi','ndwi','lai','weight','satFlag'])
    for i in range(len(landsatFiles)):
        sam_file = landsatFiles[i]
    
        df = df.append(pd.read_csv(sam_file,delim_whitespace=True,names=['ulx','uly','blue',
        'green','red','nir','swir1','swir2','ndvi','ndwi','lai','weight','satFlag']),ignore_index=True)
    
    #=====create filestem.data====================================
    df = df[(df.satFlag=='N')]
    df = df.sort_values(by='weight')
    startDate='200'
    endDate = '300'
    filestem = os.path.join(landsatLAI,"lndsr_modlai_samples.combined_%s-%s" %(startDate,endDate))
    df.to_csv(os.path.join(landsatLAI,filestem+".data"), columns = ['blue','green','red',
    'nir','swir1','swir2','ndvi','ndwi','lai','weight'],header=None, 
    index=None, mode='w',  sep="\t", encoding='utf-8')
    
    #====create filestem.names====================================
    fn = os.path.join(landsatLAI,"%s.names" % filestem)
    file = open(fn, "w")
    file.write("lai.\n")
    file.write("B1: continuous\n")
    file.write("B2: continuous\n")
    file.write("B3: continuous\n")
    file.write("B4: continuous\n")
    file.write("B5: continuous\n")
    file.write("B7: continuous\n")
    file.write("ndvi: continuous\n")
    file.write("ndwi: continuous\n")
    file.write("lai: continuous\n")
    file.write("case weight: continuous\n")
    file.write("attributes excluded: B1, B2, B7, ndvi, ndwi\n")
    file.close()
    
    nrules = 5
    subprocess.call(["%s" % cubist , "-f" ,"%s" % filestem, "-r", "%d" % nrules, "-u"])
    
def compute():    
    #lndbio = os.path.join(base,'processData','bin','lndlai_compute')
    lndbio ='lndlai_compute'
    bands = ["blue","green","red","nir","swir1","swir2","cloud"]
    l8bands = ["sr_band2","sr_band3","sr_band4","sr_band5","sr_band6","sr_band7","cfmask"] 
    
    landsatFiles = glob.glob(os.path.join(landsatTemp,"*.xml"))
    for i in range(len(landsatFiles)):
        sceneID = landsatFiles[i].split(os.sep)[-1][:-4]        
        fstem = landsatFiles[i][:-4]       
        # create a folder for lai if it does not exist
        #laiPath = os.path.join(landsatLAI,'%s' % sceneID[9:16])
        laiPath = os.path.join(landsatLAI,'%s' % sceneID[3:9])
        if not os.path.exists(laiPath):
            os.mkdir(laiPath)
        
        startDate='200'
        endDate = '300'
        filestem = os.path.join(landsatLAI,"lndsr_modlai_samples.combined_%s-%s" %(startDate,endDate))
        laiFN = os.path.join(landsatLAI,"lndlai.%s.hdf" % sceneID)
        fn = os.path.join(landsatLAI,"compute_lai%s.inp")
        file = open(fn, "w")
        file.write("LANDSAT_BASE_BLUE = %s_%s.%s.dat\n" % (fstem,l8bands[0],bands[0]))
        file.write("LANDSAT_BASE_GREEN = %s_%s.%s.dat\n" % (fstem,l8bands[1],bands[1]))
        file.write("LANDSAT_BASE_RED = %s_%s.%s.dat\n" % (fstem,l8bands[2],bands[2]))
        file.write("LANDSAT_BASE_NIR = %s_%s.%s.dat\n" % (fstem,l8bands[3],bands[3]))
        file.write("LANDSAT_BASE_SWIR1 = %s_%s.%s.dat\n" % (fstem,l8bands[4],bands[4]))
        file.write("LANDSAT_BASE_SWIR2 = %s_%s.%s.dat\n" % (fstem,l8bands[5],bands[5]))
        file.write("LANDSAT_BASE_CLOUD = %s_%s.%s.dat\n" % (fstem,l8bands[6],bands[6]))
        file.write("LANDSAT_ANC_FILE = %s\n" % filestem)
        file.write("BIOPHYSICS_PARA_FILE_OUT = %s\n" % laiFN)
        file.close()
        
        subprocess.call(["%s" % lndbio , "%s" % fn])
        shutil.move(laiFN,os.path.join(laiPath,"lndlai.%s.hdf" % sceneID))
        os.remove(fn)
    #=====CLEANING UP========
    filelist = [ f for f in os.listdir(landsatLAI) if f.startswith("lndsr_modlai_samples") ]
    for f in filelist:
        os.remove(os.path.join(landsatLAI,f))
    
def getLAI():    
    # Convert Landsat SR downloads to ENVI format
    # Note:  May be some warnings about unknown field - ignore.
    print("Converting Landsat SR to ENVI format...")
    geotiff2envi()
    
    # Generate MODIS-Landsat samples for LAI computation
    print("Generating MODIS-Landsat samples...")
    sample()    
    
    # Compute Landsat LAI
    print("Computing Landsat LAI...")
    train()
    compute()    

def main():
    # Get time and location from user
    parser = argparse.ArgumentParser()
    parser.add_argument("lat", type=float, help="latitude")
    parser.add_argument("lon", type=float, help="longitude")
    parser.add_argument("startDate", type=str, help="Start date yyyy-mm-dd")
    parser.add_argument("endDate", type=str, help="Start date yyyy-mm-dd")
    args = parser.parse_args()
      
    loc = [args.lat,args.lon] 
    startDate = args.startDate
    endDate = args.endDate

    # set project base directory structure
    #41.18,-96.43

    
    # =====USGS credentials===============
     # need to get this from pop up
    usgsUser = str(getpass.getpass(prompt="usgs username:"))
    if keyring.get_password("usgs",usgsUser)==None:
        usgsPass = str(getpass.getpass(prompt="usgs password:"))
        keyring.set_password("usgs",usgsUser,usgsPass)
    else:
        usgsPass = str(keyring.get_password("usgs",usgsUser)) 
    
    
     # =====earthData credentials===============
    earthLoginUser = str(getpass.getpass(prompt="earth login username:"))
    if keyring.get_password("nasa",earthLoginUser)==None:
        earthLoginPass = str(getpass.getpass(prompt="earth login password:"))
        keyring.set_password("nasa",earthLoginUser,earthLoginPass)
    else:
        earthLoginPass = str(keyring.get_password("nasa",earthLoginUser)) 
        
    
    #start Landsat order process
    l8_tiles = getLandsatData(loc,startDate,endDate,("%s"% usgsUser,"%s"% usgsPass))
    
    # find MODIS tiles that cover landsat scene
    # MODIS products   
    product = 'MCD15A3'
    version = '005'
    [v,h]= latlon2MODtile(args.lat,args.lon)
    tiles = "h%02dv%02d" %(h,v)
    #tiles = 'h10v04,h10v05'
    
    # download MODIS LAI over the same area and time
    print("Downloading MODIS data...")
    getMODISlai(tiles,product,version,startDate,endDate,("%s"% earthLoginUser,"%s"% earthLoginPass))
    
    
    # move surface relectance files and estimate get LAI
    downloadFolder = os.path.join(base,'processData','espa_downloads')
    folders2move = glob.glob(os.path.join(downloadFolder ,'*'))
    L8bands = ["sr_band2","sr_band3","sr_band4","sr_band5","sr_band6","sr_band7","cfmask"]
    bands = ["blue","green","red","nir","swir1","swir2","cloud"]
    for i in range(len(folders2move)):
        inputFN = folders2move[i]
        sceneID = (inputFN).split(os.sep)[-1].split('-')[0]
        scene = sceneID[3:9]
        folder = os.path.join(landsatSR,scene)
        if not os.path.exists(folder):
            os.mkdir(folder)
    
        for filename in glob.glob(os.path.join(inputFN, '*.*')):
            shutil.copy(filename, folder)  
        for filename in glob.glob(os.path.join(folder, '*.*')):
            shutil.copy(filename, landsatTemp)   
    getLAI()
    #======Clean up folder===============================
    shutil.rmtree(downloadFolder)
    #shutil.rmtree(landsatTemp)
    print("All done!")

if __name__ == "__main__":
    try:
        main()
    except (KeyboardInterrupt, pycurl.error):
        exit('Received Ctrl + C... Exiting! Bye.', 1)   
    