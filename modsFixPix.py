#!/usr/bin/env python
#
# modsFixPix - Fix bad columns on raw MODS images using a bad pixel list
#
# Usage: modsFixPix.py [options] inFile.fits outFile.fits
#
# Where:
#   inFiles.fits = FITS image file to fix
#   outFile.fits = output FITS file to create
#
# Options:
#   -V      print version info and exit
#   -v      print verbose debugging info during execution
#   -w nPix width of the averaging region either side of the
#             bad pixel region.  Default: 5 pixels
#   -l file use the custom bad pixel list in file
#             Default: the standard instrument bad pixel file
#                      in MODS_DBDIR (e.g., mods1r.bpl)
#
# Description:
#
#   Fixes bad columns on a MODS CCD image by replacing bad pixels
#   within a region with the median of 'good' pixels flanking that
#   region by +/-numAvg.  The default averaging region size is
#   +/-avgSize pixels (defined below), but may be overridden with the
#   -w flag.  The median in the flanking regions is computed and the
#   pixels replaced by 0.5*(med_left+med_right).
#
#   The default bad pixel list is created from the INSTRUME keyword in
#   the FITS header, thus mods1r.bpl or mods2b.bpl.  If the MODS_DBDIR
#   environment variable is undefined, it looks in the current working
#   directory, otherwise it looks in MODS_DBDIR if no path
#   specification is given with the filename.  A custom bad pixel list
#   file may be used with the -l option.
#
# Bad Pixel List Format:
#
#   The bad pixel list format is the same as used by iraf: a 4-column
#   ASCII text file with a space character delimiter.
#
#   For example:
#     #
#     # Bad pixel list for my CCD
#     #
#     xstart1 xend1 ystart1 yend1
#     xstart2 xend2 ystart2 yend2
#     ...
#
#   Units are unbinned pixels starting with pixel (1,1)
#
# Environment Variables:
#   MODS_DBDIR - directory path to MODS database files.  Defaults to ./
#                if undefined.
#
# Standard Modules:
#   string, os, sys, getopt
#
# External Modules:
#   numpy
#   astropy.io fits
#
# Author:
#   Rick Pogge, OSU Astronomy Department
#   pogge.1@osu.edu
#
# Modification History
#   2012 Mar 28: alpha-testing [rwp/osu]
#   2013 Dec 01: v2 modifications
#   2014 Aug 24: using astropy.io.fits instead of unsupported pyfits
#                and suppressing nuisance FITS header warnings [rwp/osu]
#   2016 Mar 08: release version [rwp]
#   2017 May 18: patch for astropy.io fits .update() deprecation [rwp/osu]
#
#-----------------------------------------------------------------------------

import string as str
import os 
import sys
import getopt
import numpy as np 
#import pyfits as fits
from astropy.io import fits

# Version and Date

versNum = '2.0.4'
versDate = '2018-10-24'

# Globals

avgSize = 5       # Default size of the averaging region for bad columns
modsDir = os.getenv('MODS_DBDIR','./')  # Default database directory


# Suppress nuisance warning from astropy.io.fits, which can be
# unusually strict about FITS standard, but the gunk it complains
# about is otherwise harmless (blank keywords)

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
  print('\nUsage: modsFixPix [options] inFile outFile')
  print('\nWhere:')
  print('  inFiles = FITS image file to fix')
  print('  outFile = output FITS file to create')
  print('\nOptions:')
  print('  -w nPix width of the averaging region either side of the')
  print('            bad pixel region.  Default: %d pixels' % (avgSize))
  print('  -l file use the custom bad pixel list in file')
  print('            Default: the standard instrument bad pixel file')
  print('            in MODS_DBDIR (e.g., mods1r.bpl)')
  print('  -f      force overwrite of an existing outFile (\'clobber\')')
  print('  -v      print verbose debugging info during execution')
  print('  -V      print version info and exit')
  print('\nEnvironment:')
  print('  MODS_DBDIR = directory path to MODS database files')
  print('               currently MODS_DBDIR=%s\n' % (modsDir))


#----------------------------------------------------------------
#
# readBadPixList(file) - read contents of a bad pixel list
#
# IRAF-format bad pixel coordinate list:
#    xbegin xend ybegin yend
#
#

def readBadPixList(file):
  F = open(file,'r')
  M = F.readlines()
  F.close()

  xs=[]
  xe=[]
  ys=[]
  ye=[]

  for i in range(len(M)): 
    if str.find(M[i], '#') <0: 
      xs.append(float(M[i].split()[0])) 
      xe.append(float(M[i].split()[1])) 
      ys.append(float(M[i].split()[2])) 
      ye.append(float(M[i].split()[3])) 

  if len(xs) == 0:
    print('\n** ERROR: %s has no bad regions listed' % (file))
    print('          modsFixPix aborting.\n')
    sys.exit(1)

  return xs, xe, ys, ye

#----------------------------------------------------------------
#
# Main Program starts here...
#

#def main():

# Defaults for various flags

userBPL = False  # if True, use the user-supplied bad pixel list file
Verbose = False  # if True, print detailed debugging output
numAvg  = avgSize # default averaging size
NoClobber = True  # do not allow overwrite on output (override with -f)
  
# Parse the command-line arguments (GNU-style getopt)
  
try:
  opts, files = getopt.gnu_getopt(sys.argv[1:],'l:w:vVf',
                                  ['bpl=','list=','width=','verbose',
                                   'version','clobber'])
except getopt.GetoptError as err:
  print('\n** ERROR: %s' % (err))
  printUsage()
  sys.exit(2)

if len(opts)==0 and len(files)==0:
  printUsage()
  sys.exit(1)

for opt, arg in opts:
  if opt in ('-l','--bpl','--list'):
    bplFile = arg
    userBPL = True

  elif opt in ('-w','--width'):
    numAvg = int(arg)
    if numAvg < 3:
      print('** ERROR: width of the averaging region must be 3 pixel or larger')
      sys.exit(1)

  elif opt in ('-f','--clobber'):
    NoClobber = False

  elif opt in ('-v','--verbose'):
    Verbose = True

  elif opt in ('-V','--version'):
    print('modsFixPix.py v%s [%s]' % (versNum, versDate))
    sys.exit(0)
      
numFiles = len(files)

if numFiles < 2:
  printUsage()
  sys.exit(0)
else:
  inFile = files[0]
  outFile = files[1]

# Make sure the input FITS file exists before proceeding

if os.path.isfile(inFile) == 0:
  print('\n** ERROR: Input FITS image file %s not found.' % (inFile))
  print('          modsFixPix aborting.\n')
  sys.exit(1)

# Make sure we won't clobber the output file unless -f is used

if NoClobber and os.path.isfile(outFile): 
  print('\n** ERROR: Operation would overwrite existing FITS file %s' % (outFile))
  print('          Use -f to force output (\'clobber\').')
  print('          modsFixPix aborting.\n')
  sys.exit(1)

# Open the input FITS image file, and get header info as needed

fitsFile = fits.open(inFile,uint=False)
data = fitsFile[0].data

# Extract useful FITS info from the header

naxis1 = fitsFile[0].header['NAXIS1']
naxis2 = fitsFile[0].header['NAXIS2']
instID = fitsFile[0].header['INSTRUME']
channel = fitsFile[0].header['CHANNEL']

if Verbose:
  print('Processing %s [%d x %d]' % (inFile, naxis1, naxis2))
  print('%s %s Channel' % (instID,channel))

# If we were not given a bad pixel list on the command line, use the
# INSTRUME keyword in the header to find the default

if userBPL:
  pathBits = os.path.split(bplFile)
  if len(pathBits[0]) == 0:
    bplFile = os.path.join(modsDir,bplFile)
else:
  bplFile = os.path.join(modsDir,str.lower(instID)+'.bpl')

if os.path.isfile(bplFile) == 0:
  print('\n** ERROR: MODS bad pixel list file %s not found.' % (bplFile))
  print('          modsFixPix aborting.\n')
  sys.exit(1)

if Verbose:
  print('Bad Pixel List: %s' % (bplFile))

# Read and parse the bad pixel list

xs, xe, ys, ye = readBadPixList(bplFile)

if Verbose:
  print('%s contains %d bad pixel regions to fix' % (bplFile,len(xs)+1))

# Page through the bad pixel list and fix the pixel

for i in range(len(xs)):
  xstart = int(xs[i]-1)
  xend = int(xe[i])
  ystart = int(ys[i]-1)
  yend = int(ye[i])
  nx=xe[i]-xs[i]+1

  if Verbose:
    print('  Fixing pixels in region [%d:%d,%d:%d]' % (xs[i],xe[i],ys[i],ye[i]))
    
  # The 'fix' is to replace bad-column pixels by the median of
  # pixels +/-numAvg either side of the bad column(s).

  xsl=xstart-numAvg
  xel=xstart-1
  xsr=xend+1
  xer=xend+numAvg

  for y in range(ystart,yend,1):
    medLeft = np.median(data[y,xsl:xel])
    medRight= np.median(data[y,xsr:xer])
    corrVal = (medLeft+medRight)/2.0
    if (nx>1):
      data[y,xstart:xend] = corrVal
    else:
      data[y,xstart] = corrVal

if Verbose:
  print('Done fixing pixels, writing image %s...' % (outFile))

# Write the output file 

fitsFile[0].data = data
fitsFile[0].header['BPLFILE'] = (bplFile)
fitsFile[0].header['BPWIDTH'] = (numAvg,'Bad pixel averaging region width [pixels]')
fitsFile[0].header.add_history('modsFixPix v%s %s' % (versNum,versDate))

if os.path.isfile(outFile): 
  if Verbose:
    print('** WARNING: Overwriting existing FITS file %s' % (outFile) )
  os.remove(outFile) 

fitsFile.writeto(outFile,output_verify='ignore')

if Verbose:
  print('modsFixPix Done')

sys.exit(0)

#---------------------------------------------------------------------------
#
#if __name__ == '__main__':
#  main()

