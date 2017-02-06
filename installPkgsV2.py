#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 15:35:20 2016

@author: mschull
"""
from __future__ import print_function
import subprocess
import os
import glob
import platform

# set project base directory structure
p = os.getcwd()
base = os.path.abspath(os.path.join(p,os.pardir))

# Create folder structure
binDir = os.path.join(base,'processData','bin')
if not os.path.exists(binDir):
    os.makedirs(binDir)


# SET Channels
subprocess.call(["conda", "config", "--add", "channels", "conda-forge"])
subprocess.call(["conda", "config", "--add", "channels", "dbrown"]) #szip

# SET ENVIRONMENT
envName = "pyDMS34"
subprocess.call(["conda", "create", "--name", "%s" % envName, "python=3.4"])

# get Anaconda root location
p = subprocess.Popen(["conda", "info", "--root"],stdout=subprocess.PIPE)
out = p.communicate()
condaPath = out[0][:-1]
libEnv = os.path.join(condaPath,'lib')    

processDir = os.path.join(base,'processData')
libDir = os.path.join(processDir,'source','lib')


#print('installing libraries...')
##====INSTALL Libraries==================
#subprocess.call(["conda", "install", "scons"])
#subprocess.call(["conda", "install", "gcc"])
#subprocess.call(["conda", "install", "jpeg"])
#subprocess.call(["conda", "install", "geotiff"])
#subprocess.call(["conda", "install", "zlib"])
#subprocess.call(["conda", "install", "hdfeos2"])
#subprocess.call(["conda", "install", "szip"])
#subprocess.call(["conda", "install", "libtiff"])
#
#
#
##====Creating SYMBOLIC L:INKS===========
#subprocess.call(["ln","-s", "%s" % os.path.join(libEnv,'libtiff.a'), 
#"%s" % os.path.join(libDir,'libtiff.a')])
#subprocess.call(["ln","-s", "%s" % os.path.join(libEnv,'liblzma.a'), 
#"%s" % os.path.join(libDir,'liblzma.a')])
#subprocess.call(["ln","-s", "%s" % os.path.join(libEnv,'libjpeg.a'), 
#"%s" % os.path.join(libDir,'libjpeg.a')])
#subprocess.call(["ln","-s", "%s" % os.path.join(libEnv,'libgeotiff.a'), 
#"%s" % os.path.join(libDir,'libgeotiff.a')])
#subprocess.call(["ln","-s", "%s" % os.path.join(libEnv,'libsz.a'), 
#"%s" % os.path.join(libDir,'libsz.a')])
#subprocess.call(["ln","-s", "%s" % os.path.join(libEnv,'libGctp.a'), 
#"%s" % os.path.join(libDir,'libGctp.a')])
#subprocess.call(["ln","-s", "%s" % os.path.join(libEnv,'libdf.a'), 
#"%s" % os.path.join(libDir,'libdf.a')])
#subprocess.call(["ln","-s", "%s" % os.path.join(libEnv,'libz.a'), 
#"%s" % os.path.join(libDir,'libz.a')])
#subprocess.call(["ln","-s", "%s" % os.path.join(libEnv,'libmfhdf.a'), 
#"%s" % os.path.join(libDir,'libmfhdf.a')])
#subprocess.call(["ln","-s", "%s" % os.path.join(libEnv,'libhdfeos.a'), 
#"%s" % os.path.join(libDir,'libhdfeos.a')])

mkPath = os.path.join(processDir,'source','Landsat_LAI','lndlai_compute')
os.chdir(mkPath)
subprocess.call(["scons","-Q","install"])
subprocess.call(["scons","-c"])

mkPath = os.path.join(processDir,'source','Landsat_LAI','lndlai_sample')
os.chdir(mkPath)
subprocess.call(["scons","-Q","install"])
subprocess.call(["scons","-c"])

mkPath = os.path.join(processDir,'source','Landsat_LAI','GeoTiff2ENVI')
os.chdir(mkPath)
subprocess.call(["scons","-Q","install"])
subprocess.call(["scons","-c"])

mkPath = os.path.join(processDir,'source','Landsat_DMS')
os.chdir(mkPath)
subprocess.call(["scons","-Q","install"])
subprocess.call(["scons","-c"])

mkPath = os.path.join(processDir,'source','MODIS_DMS')
os.chdir(mkPath)
subprocess.call(["scons","-Q","install"])
subprocess.call(["scons","-c"])

mkPath = os.path.join(processDir,'source','Cubist')
os.chdir(mkPath)
subprocess.call(["scons","-Q","install"])
subprocess.call(["scons","-c"])


#if make_process.wait() != 0:
#     print("somethings wrong")
##=====INSTALLING FENG'S CODES==========
#
#print("installing Gao's packages...")
#source = os.path.join(base,'processData','source')
#os.chdir(source)
#subprocess.call(["bash","build_all.sh"])

##======install PYRTTOV==================
#
##!. edit Makefil.local
#
#
##2. interactive install
#
##3  make symbolic links to needed files
#
#subprocess.call(["ln","-s", "%s" % os.path.join(processDir,'rttov113','lib',
#'rttov_wrapper_f2py.so'), "%s" % os.path.join(processDir,'LST','rttov_wrapper_f2py.so')])
#
#subprocess.call(["ln","-s", "%s" % os.path.join(processDir,'rttov113','rttov_wrapper_pyrttov_20160127',
#'wrapper','pyrttov'), "%s" % os.path.join(processDir,'LST','pyrttov')])
