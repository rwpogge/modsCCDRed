#! /usr/bin/env python
#
# modsMedian - Median combine images
#
# Usage: modsMedian [options] inFile1.fits [inFile2.fits ...] medFile.fits
# 
# where: inFile?.fits = input FITS images to be median combined
#        medFile.fits = output median image to create
# 
# Options:
#   -f      force overwrite of an existing median file ('clobber')
#   -v      print verbose debugging info during execution
#   -V      print version info and exit
#
# Description:
#
#   Median combine a list of images.  Because of the way python
#   allocates memory, this could be limited to 10 full-size MODS
#   images before users will encounter memory problems.  When this
#   happens depends on your computer hardware.
#
# Authors:
#   R. Pogge, OSU Astronomy Dept
#   pogge.1@osu.edu
#
# Module Dependencies:
#   numpy
#   astropy.io fits
#
# Modification History
#   2012 March 28 - fixed bug in file counter [rwp/osu]
#   2012 March 30 - added getopt features [rwp/osu]
#   2014 Aug 24 - using astropy.io.fits to replace unsupported pyfits 
#                 and suppressing nuisance FITS header warnings [rwp/osu]
#   2017 May 18: patch for astropy.io fits .update() deprecation [rwp/osu]
#
#-----------------------------------------------------------------------------

import string as str
import os 
import sys
import getopt
import numpy as np 
from astropy.io import fits

# Version Number and Date

versNum = '2.0.3'
versDate = '2017-05-18'

# Runtime flags

Verbose = False    # Default: verbose debugging output off
NoClobber = True   # Default: no-clobber on output

# Suppress nuisance warning from astropy.io.fits, which can
# be unusually strict about FITS standard, but the gunk it 
# complains about is otherwise harmless (blank keywords)

import warnings
warnings.filterwarnings('ignore',category=UserWarning, append=True)
warnings.filterwarnings('ignore',category=RuntimeWarning, append=True)
warnings.filterwarnings('ignore',category=FutureWarning, append=True)
from astropy.utils.exceptions import AstropyWarning
warnings.filterwarnings('ignore',category=AstropyWarning, append=True)

#----------------------------------------------------------------
#
# printUsage - print a simple usage message
#

def printUsage(): 
  print('\nUsage: modsMedian [options] inFile1.fits [inFile2.fits ...] medFile.fits')
  print('\nwhere: inFile?.fits = input FITS files to be median combined.')
  print('       medFile.fits = output median image to create')
  print('\nOptions:')
  print('  -f      force overwrite of an existing median file (\'clobber\')')
  print('  -v      print verbose debugging info during execution')
  print('  -V      print version info and exit\n')
  
#----------------------------------------------------------------
#
# Main Program starts here...
#

# Parse the command-line arguments (GNU-style getopt)
  
try:
  opts, files = getopt.gnu_getopt(sys.argv[1:],'fvV',['clobber','verbose','version',])
except getopt.GetoptError as err:
  print('\n** ERROR: %s' % (err))
  printUsage()
  sys.exit(2)

if len(opts)==0 and len(files)==0:
  printUsage()
  sys.exit(1)

for opt, arg in opts:
  if opt in ('-f','--clobber'):
    NoClobber = False
  elif opt in ('-v','--verbose'):
    Verbose = True
  elif opt in ('-V','--version'):
    print('modsMedian.py v%s [%s]' % (versNum, versDate))
    sys.exit(0)

numFiles = len(files)

if numFiles == 0:
  printUsage()
  sys.exit(0)
elif numFiles < 4:
  print('\n** ERROR: Only %d files to median combine.')
  print('          You need at least 3 files for a robust median' % (numFiles-1))
  print('          modsMedian aborting.\n')
  sys.exit(1)

# Create the output filename, make sure we won't clobber it unless -f is used

outFile = files[numFiles-1]
if NoClobber and os.path.isfile(outFile): 
  print('\n** ERROR: Operation would overwrite existing FITS file %s' % (outFile))
  print('          Use -f to force output (\'clobber\').')
  print('          modsMedian aborting.\n')
  sys.exit(1)

# Do it!

if Verbose:
  print('\nMedian combining %d images...' % (numFiles-1))

kount=0
for rawFile in files:
  if kount==(numFiles-1):
    if Verbose:
      print('Done, computing the median...')
  else:
    if os.path.isfile(rawFile):
      if kount==0:
        fitsFile = fits.open(rawFile,uint=False)
        if Verbose:
          print('  Reading image %d of %d: %s...' % (kount+1,numFiles-1,rawFile))
        d0 = fitsFile[0].data
        kount += 1
      else:
        stackFile = fits.open(rawFile,uint=False)
        if Verbose:
          print('  Reading image %d of %d: %s...' % (kount+1,numFiles-1,rawFile))
        d1 = stackFile[0].data
        if kount==1:
          cube = np.dstack((d0,d1))
        else:
          cube = np.dstack((cube,d1))
        kount += 1
        stackFile.close()

# Now that we have the data cube, median it

d0 = np.median(cube,axis=2)
fitsFile[0].data = d0
fitsFile[0].header.add_history('modsMedian v%s %s nfiles=%d' % (versNum,versDate,numFiles-1))

if Verbose:
  print('Writing median to output file %s...' % (outFile))
  
if os.path.isfile(outFile): 
  if Verbose:
    print('** WARNING: Overwriting existing FITS file %s' % (outFile) )
  os.remove(outFile) 

fitsFile.writeto(outFile,output_verify='ignore')
    
if Verbose:
  print('modsMedian Done')

sys.exit(0)
