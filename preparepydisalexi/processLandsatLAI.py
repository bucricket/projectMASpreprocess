#!/usr/bin/env python2
#!/Applications/anaconda/envs/root3 python
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 13:50:04 2016

@author: mschull
"""
#python

#from .search import Search
import os
import subprocess
import glob
import shutil
import pandas as pd
from datetime import datetime
import argparse
import getpass
import keyring
import json
from pyproj import Proj
from .utils import folders,search,checkOrderCache
from .Clients import Client
from Downloaders import BaseDownloader
#from .Order import Order
#from .OrderTemplate import OrderTemplate
import pycurl
from .landsatTools import landsat_metadata
import requests
from time import sleep
import logging
logging.getLogger("urllib3").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.DEBUG)
#import json
#import getpass

# The current URL hosting the ESPA interfaces has reached a stable version 1.0
host = 'https://espa.cr.usgs.gov/api/v1/'
timeout=86400

base = os.getcwd()
Folders = folders(base)   
modisBase = Folders['modisBase']
landsatSR = Folders['landsatSR']
landsatLAI = Folders['landsatLAI']
landsatTemp = os.path.join(landsatSR,'temp')
if not os.path.exists(landsatTemp):
    os.mkdir(landsatTemp)

#def getLandsatData(collection,loc,startDate,endDate,auth):
#    
#    
#    data = {'olitirs8':{"inputs":[],"products": ["sr", "bt","cloud"]},"format":"gtiff",
#        "plot_statistics":False,"note":""}    
#    template = OrderTemplate('order',template_dir='./')
#    template.define(data)
#    template.save()
#    with open('order.json', 'w') as outfile:  
#        json.dump(data, outfile)
#    
#    # build the various handlers to spec
#    template = OrderTemplate('template')
#    template.load(path='./order.json' )
#    template.load('./order.json' )
#    order = Order('order')
#    order = Order('order', note="Lat%dLon%d-%s_%s" %(int(loc[0]),int(loc[1]),startDate,endDate))
#    client = Client(auth) # will prompt user for username and password if auth argument not supplied
#    #downloader = EspaLandsatLocalDownloader('USGS_downloads')
#    
#    # find cloud free landsat scenes
#    try:
#        s = Search()
#        scenes = s.search(lat=loc[0],lon=loc[1],limit = 100, start_date = startDate,end_date=endDate, cloud_max=5)
#        l8_tiles=[]
#        for i in range(len(scenes['results'])):
#            path = scenes['results'][i]['path']
#            row = scenes['results'][i]['row']
#            sceneID = scenes['results'][i]['sceneID'][:-1]+'%s' % collection
#            if sceneID.startswith('LC'):
#                dataFN = os.path.join(landsatSR,"%s%s" %(path,row),"%s%s.xml" % sceneID)
#                if not os.path.exists(dataFN):
#                    l8_tiles.append(sceneID)
#                else:
#                    files = glob.glob("%s*" % dataFN[:-4])
#                    for file in files:
#                        os.symlink(file,os.path.join(landsatTemp,file.split(os.sep)[-1]))
#                        #shutil.copy(file,landsatTemp)
#    except:
#        print("crappy version")
#        sceneIDs = search(loc[0],loc[1],startDate, endDate)
#
#        l8_tiles=[]
#        for i in range(len(sceneIDs)):
#            l8_tiles.append(sceneIDs[i])
#    print l8_tiles
#    if l8_tiles:    
#        # order the data
#        order.add_tiles("olitirs8", l8_tiles)
#        #order.add_tiles("etm7", l7_tiles)
#        response = order.submit(client)
#        
#        # view the servers whole response. which might indicate an ordering error!
#        print(response) 
#        
#        # assuming there were no order submission errors
#        orderid = response['orderid']
#        
#        # now start the downloader!
#        for download in client.download_order_gen(orderid):
#            print(download)
#        # download is a tuple with the filepath, and True if the file
#        # is a fresh download.
#        
#        # this is where data pipeline scripts go that can operate
#        # on files as they are downloaded (generator),
#        
#        # See the Client class for further documentation.

def getLandsatData(collection,loc,startDate,endDate,auth,cloud):
    username = auth[0]
    password = auth[1]
    client = Client(auth)
    def api_request(endpoint, verb='get', json=None, uauth=None):
        """
        Here we can see how easy it is to handle calls to a REST API that uses JSON
        """
        auth_tup = uauth if uauth else (username, password)
        response = getattr(requests, verb)(host + endpoint, auth=auth_tup, json=json)
        return response.json()

    #=====set products=======
    l8_prods = ['sr','bt','cloud']
    #=====search for data=======
    print("Searching...")
    sceneIDs = search(collection,loc[0],loc[1],startDate,endDate,cloud)
    orderedData = checkOrderCache(auth)
    l8_tiles =[]
    completedOrderedIDs = []
    notCompletedOrderedIDs = []
    completedSceneIDs = []
    notCompletedSceneIDs =[]
    for sceneID in sceneIDs:
        if sceneID.startswith('LC'):
            dataFN = os.path.join(landsatSR,"%s" % sceneID.split('_')[2],"%s_MTL.txt" % sceneID)
            if not os.path.exists(dataFN):
                if len(orderedData[(orderedData.productID==sceneID)])>0:
                    if len(orderedData[(orderedData.productID==sceneID) & (orderedData.status=='complete')]['orderid'])>0:
                        completedOrderedIDs.append(list(orderedData[(orderedData.productID==sceneID) & (orderedData.status=='complete')]['orderid'])[0])
                        completedSceneIDs.append(sceneID)
                    else:
                        notCompletedOrderedIDs.append(list(orderedData[(orderedData.productID==sceneID) & (orderedData.status=='oncache')]['orderid'])[-1])
                        notCompletedSceneIDs.append(sceneID)
                else:
                    l8_tiles.append(sceneID)
            else:
                files = glob.glob("%s*" % dataFN[:-4])
                for file in files:
                    os.symlink(file,os.path.join(landsatTemp,file.split(os.sep)[-1]))


    
    if l8_tiles:
        print("Ordering new data...")
        #========setup order=========
        order = api_request('available-products', verb='post', json=dict(inputs=l8_tiles))
        for sensor in order.keys():
            if isinstance(order[sensor], dict) and order[sensor].get('inputs'):
                order[sensor]['products'] = l8_prods
        
        order['format'] = 'gtiff'
        # =======order the data============
        resp = api_request('order', verb='post', json=order)
        print(json.dumps(resp, indent=4))
        orderidNew = resp['orderid']
                    
    if completedOrderedIDs:
        print("downloading completed existing orders...")
        i = -1
        for orderid in completedOrderedIDs:
            i+=1
            sceneID = completedSceneIDs[i]
            complete = False
            reached_timeout = False
            starttime = datetime.now()
            while not complete and not reached_timeout:
                resp = api_request('item-status/{0}'.format(orderid))
                for item in resp['orderid'][orderid]:
                    if item.get('name')==sceneID:
                        url = item.get('product_dload_url')                      
                        elapsed_time = (datetime.now() - starttime).seconds
                        reached_timeout = elapsed_time > timeout
                        print("Elapsed time is {0}m".format(elapsed_time / 60.0))
                        if len(url)>0:
                            downloader = BaseDownloader('espa_downloads')
                            downloader.download(url)
                        #if os.path.exists(os.path.join(os.getcwd,'espa_downloads',url.split(os.sep)[-1][:-7])):
                            complete = True
                        
                        if not complete:
                            sleep(300)
                        
                        
    if notCompletedOrderedIDs:
        print("waiting for cached existing orders...")
        i = -1
        for orderid in notCompletedOrderedIDs:
            i+=1
            complete = False
            reached_timeout = False
            starttime = datetime.now()
            sceneID = notCompletedSceneIDs[i]
            while not complete and not reached_timeout:
                resp = api_request('item-status/{0}'.format(orderid))
                for item in resp['orderid'][orderid]:
                    if item.get('name')==sceneID:
                        url = item.get('product_dload_url')
                        elapsed_time = (datetime.now() - starttime).seconds
                        reached_timeout = elapsed_time > timeout
                        print("Elapsed time is {0}m".format(elapsed_time / 60.0))
                        if len(url)>0:
                            downloader = BaseDownloader('espa_downloads')                            
                            downloader.download(url)
                        #if os.path.exists(os.path.join(os.getcwd,'espa_downloads',url.split(os.sep)[-1][:-7])):
                            complete = True
                        
                        if not complete:
                            sleep(300)

    if l8_tiles:
        print("Download new data...")       
        #======Download data=========    
        for download in client.download_order_gen(orderidNew):
            print(download)
    
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
    geotiffConvert = 'GeoTiff2ENVI'
    bands = ["blue","green","red","nir","swir1","swir2","cloud"]
    l8bands = ["sr_band2","sr_band3","sr_band4","sr_band5","sr_band6","sr_band7","cfmask"] 
    
    landsatFiles = glob.glob(os.path.join(landsatTemp,"*_MTL.txt"))
    for i in range(len(landsatFiles)):
        
        fn = landsatFiles[i][:-8]
        meta = landsat_metadata(landsatFiles[i])
        fstem = os.sep.join((fn.split(os.sep)[:-1]))+meta.LANDSAT_SCENE_ID

        for i in range(len(bands)):
            tifFile = fn+"_%s.tif" % l8bands[i]
            datFile = fstem+"_%s.%s.dat" % (l8bands[i],bands[i])
            subprocess.call(["%s" % geotiffConvert ,"%s" % tifFile, "%s" % datFile])

def sample():    
    sample = 'lndlai_sample'
    bands = ["blue","green","red","nir","swir1","swir2","cloud"]
    l8bands = ["sr_band2","sr_band3","sr_band4","sr_band5","sr_band6","sr_band7","cfmask"] 
    
    landsatFiles = glob.glob(os.path.join(landsatTemp,"*_MTL.txt"))
    
    for i in range(len(landsatFiles)):
        #sceneID = landsatFiles[i].split(os.sep)[-1][:-4]
        meta = landsat_metadata(landsatFiles[i])
        sceneID = meta.LANDSAT_SCENE_ID
        
        # extract the Landsat doy and year
        d = meta.DATETIME_OBJ
        year = d.year
        ldoy = sceneID[13:16]
#        year = int(sceneID[9:13])
        # convert to date    
#        dd = datetime.datetime(year, 1, 1) + datetime.timedelta(int(ldoy) - 1)
#        date = '%d-%02d-%02d' % (dd.year,dd.month,dd.day)
        date = meta.DATE_ACQUIRED
        # find the 4 day MODIS doy prior to the Landsat doy
        mdoy = int((int((float(ldoy)-1.)/4.)*4.)+1)
        
        modFiles = glob.glob(os.path.join(modisBase,"MCD15A3.A%s%s.*.hdf" % (year,mdoy)))

        #fstem = landsatFiles[i][:-4]
        fn = landsatFiles[i][:-8]
        fstem = os.sep.join((fn.split(os.sep)[:-1]))+meta.LANDSAT_SCENE_ID
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
    cubist = 'cubist'
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
    lndbio ='lndlai_compute'
    bands = ["blue","green","red","nir","swir1","swir2","cloud"]
    l8bands = ["sr_band2","sr_band3","sr_band4","sr_band5","sr_band6","sr_band7","cfmask"] 
    
    landsatFiles = glob.glob(os.path.join(landsatTemp,"*_MTL.txt"))
    for i in range(len(landsatFiles)):
#        sceneID = landsatFiles[i].split(os.sep)[-1][:-4]  
        meta = landsat_metadata(landsatFiles[i])
        sceneID = meta.LANDSAT_SCENE_ID
        #fstem = landsatFiles[i][:-4]       
        fn = landsatFiles[i][:-8]
        fstem = os.sep.join((fn.split(os.sep)[:-1]))+meta.LANDSAT_SCENE_ID
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
    parser.add_argument("cloud", type=str, help="cloud coverage")
    args = parser.parse_args()
      
    loc = [args.lat,args.lon] 
    startDate = args.startDate
    endDate = args.endDate
    cloud = args.cloud
    collection = 1
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
    getLandsatData(collection,loc,startDate,endDate,("%s"% usgsUser,"%s"% usgsPass),cloud)
    
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
    downloadFolder = os.path.join(base,'espa_downloads')
    folders2move = glob.glob(os.path.join(downloadFolder ,'*'))

    for i in range(len(folders2move)):
        inputFN = folders2move[i]
        metFN = glob.glob(os.path.join(inputFN,'*MTL.txt'))[0]
        meta = landsat_metadata(metFN)
        #sceneID = (inputFN).split(os.sep)[-1].split('-')[0]
        sceneID = meta.LANDSAT_SCENE_ID
        scene = sceneID[3:9]
        folder = os.path.join(landsatSR,scene)
        if not os.path.exists(folder):
            os.mkdir(folder)
    
        for filename in glob.glob(os.path.join(inputFN, '*.*')):
            shutil.copy(filename, folder)  
            os.symlink(os.path.join(folder,filename.split(os.sep)[-1]),
            os.path.join(landsatTemp,filename.split(os.sep)[-1]))
 
    if len(folders2move)>0:
            #======Clean up folder===============================
            shutil.rmtree(downloadFolder)
        
    getLAI()

    print("All done with LAI")
    print("========================================")
    print("==============process LST===============")
    subprocess.call(["processlst","%s" % earthLoginUser,"%s" % earthLoginPass])
    #shutil.rmtree(landsatTemp)

if __name__ == "__main__":
    try:
        main()
    except (KeyboardInterrupt, pycurl.error):
        exit('Received Ctrl + C... Exiting! Bye.', 1)   
    