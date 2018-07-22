#! /usr/bin/env python
#
# modsSub - Subtract an image from another image
#
# Usage: modsSub [options] inFile1 inFile2 outFile
# 
# Where:
#   inFile1,2 = FITS files to difference.
#   outFile = inFile1-inFile2
# 
# Options:
#   -f  = force overwrite of the output file ('clobber')
#   -v  = print verbose debugging info during execution
#   -V  = print version information and exit
#
# Description:
#
#   Subtracts inFile2 from inFile1, creating a new FITS file
#   with the difference (outFile = inFile1 - inFile2).
#
#   Provided as a convenience function.
#
# Authors:
#   R. Pogge, OSU Astronomy Dept
#   pogge.1@osu.edu
#
# Module Dependencies:
#   numpy
#   astropy.io.fits
#
# Modification History
#   2012 Mar 28 - fixed bug in file counter [rwp/osu]
#   2012 Apr 01 - added notes and getopt functions [rwp/osu]
#   2014 Aug 24 - using astropy.io.fits instead of unsupported pyfits and
#                 suppressing nuisance FITS header warnings [rwp/osu] 
#   2016 May 08 - first release version [rwp/osu]
#   2017 May 18 - patch for astropy .update() deprecation [rwp/osu]
#
#-----------------------------------------------------------------------------

import string as str
import os 
import sys
import getopt
import numpy as np 
#import pyfits as fits
from astropy.io import fits

# Version Number and Date

versNum = '2.1.3'
versDate = '2017-05-18'

# Runtime flags

Verbose=False      # Verbose debugging output on/off
NoClobber=True     # default: no clobber on output

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
  print('\nUsage: modsSub [options] inFile1 inFile2 outFile')
  print('\nWhere:')
  print('  inFile1,2 = FITS files to difference.')
  print('  outFile = inFile1-inFile2')
  print('\nOptions:')
  print('  -f  = force overwrite of the output file (\'clobber\')')
  print('  -v  = print verbose debugging info during execution')
  print('  -V  = print version information and exit\n')

#----------------------------------------------------------------
#
# Main Program starts here...
#

# Parse the command-line arguments (GNU-style getopt)
  
try:
  opts, files = getopt.gnu_getopt(sys.argv[1:],'fvV',
                                  ['clobber','verbose','version',])
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
    print('modsSub.py v%s [%s]' % (versNum, versDate))
    sys.exit(0)

numFiles = len(files)
if numFiles != 3:
  printUsage()
  sys.exit(1)

# Create the output filename, make sure we won't clobber it...

inFile1 = files[0]
inFile2 = files[1]
outFile = files[2]

# Make sure we don't clobber the output file unless -f is used

if NoClobber and os.path.isfile(outFile): 
  print('\n** ERROR: Operation would overwrite existing FITS file %s' % (outFile))
  print('          Use -f to force output (\'clobber\').')
  print('          modsSub aborting.\n')
  sys.exit(1)

# start!

if Verbose:
  print('\nSubtracting %s from %s...' % (inFile1,inFile2))

if os.path.isfile(inFile1):
  if Verbose:
    print('  Reading File 1: %s...' % (inFile1))
  fitsFile = fits.open(inFile1,uint=False)
  d0 = fitsFile[0].data
else:
  print('\n** ERROR: Input FITS file %s not found.')
  print('          modsSub aborting.\n')
  sys.exit(1)

if os.path.isfile(inFile2):
  if Verbose:
    print('  Reading File 2: %s...' % (inFile2))
  diffFile = fits.open(inFile2,uint=False)
  if Verbose:
    print('  Computing file1-file2...')
  d0 -= diffFile[0].data
  diffFile.close()
else:
  print('\n** ERROR: Input FITS file %s not found.')
  print('          modsSub aborting.\n')
  sys.exit(1)

fitsFile[0].data = d0
fitsFile[0].header.add_history('modsSub v%s %s file1=%s file2=%s' % (versNum,versDate,inFile1,inFile2))

if Verbose:
  print('Writing output to %s...' % (outFile))

if os.path.isfile(outFile): 
  if Verbose:
    print('** WARNING: Overwriting existing FITS file %s' % (outFile) )
  os.remove(outFile) 
  
fitsFile.writeto(outFile,output_verify='ignore')
    
if Verbose:
  print('modsSub Done')

sys.exit(0)
