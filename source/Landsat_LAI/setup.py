# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#from distutils.core import setup, Extension
#
#module1 = Extension('lndlai_compute',
#                    define_macros = [('MAJOR_VERSION', '1'),
#                                     ('MINOR_VERSION', '0')],
#                    include_dirs = ['/Applications/anaconda/envs/pyDMS3/include'],
#                    libraries = ['hdfeos','Gctp','mfhdf','df','z','sz','jpeg','m'],
#                    library_dirs = ['/Applications/anaconda/envs/pyDMS3/lib'],
#                    sources = ['lndlai_compute/lndlai_compute.c'])
#
#setup (name = 'lndlai',
#       version = '1.0',
#       description = 'This is a demo package',
#       author = 'Mitch',
#       author_email = 'mitch.schull@noaa.gov',
#       url = 'https://docs.python.org/extending/building',
#       long_description = '''
#This is really just a demo package.
#''',
#       ext_modules = [module1],
#       packages=['lndlai_compute'])




import os
from setuptools import setup, Extension
 
setup(
    name='GeoTiff2ENVI',
    version='0.0.1',
    author='Mitch',
    author_email='bucricket@gmail.com',
    packages=['GeoTiff2ENVI'],
    license='Use as you wish. No guarantees whatsoever.',
    install_requires=[''],
    classifiers=['Development Status :: 3 - Alpha'],
    entry_points={'console_scripts': [
                                      'GeoTiff2ENVI=GeoTiff2ENVI.GeoTiff2ENVI.main']},
    ext_modules=[
        Extension('GeoTiff2ENVI.main',
                  ['GeoTiff2ENVI/GeoTiff2ENVI.c'],
                  include_dirs=['/Applications/anaconda/envs/pyDMS3/include'],
                  library_dirs=['/Applications/anaconda/envs/pyDMS3/lib'],
                  libraries=['tiff','geotiff','z','jpeg','m'],
                  extra_compile_args=['-g']
                 )
     ],
)
