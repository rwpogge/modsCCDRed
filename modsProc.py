#! /usr/bin/env python
#
# modsProc - bias-correct and flat field MODS raw CCD images.
#
# Usage: modsProc [options] rawFile.fits [rawFile2.fits ...] flatField.fits
# 
# Where:
#   rawFile*.fits  raw MODS 8x3K FITS image(s) to process
#   pixFlat.fits   is a normalized pixel flat to apply to each rawFile.fits
# Output file will be named rawFile_otf.fits
#
# If N images are listed on the command line, (N-1) are raw images and N is
# the flat field.  Beware!
# 
# Options:
#   -b       fix bad pixels (see -l)
#   -l file  bad pixel list to use (requires -b)
#   -w nPix  width of the averaging region either side of the
#             bad pixel region (requires -b).  Default: 5 pixels
#   -f       force overwrite of an existing output file ('clobber')
#   --noflip do not flip a Red channel image X-axis (ignored if blue)
#   -v       print verbose debugging info during execution
#   -V       print version info and exit
# 
# Description:
#
#   Processes a raw MODS image as follows:
#
#    1. Remove even/odd overscan and trim (see modsBias.py)
#    2. Divide by the flat field.
#    3. If -b is given, fix known bad columns
#    4. If a Red channel image, flip about X to give the same
#       blue-to-red/left-to-right orientation as the blue channel.
#       can be suppressed with --noflip
#
#   Creates rawFile_otf.fits, otf = overscan, trimmed, and flattened
#
# Quadrant Layout of the MODS Science CCDs:
#
#      ++--------+--------++(NAXIS1,NAXSI2)
#      ||   Q3   |   Q4   ||
#      ++--------+--------++
#      ||   Q1   |   Q2   ||
# (1,1)++--------+--------++
# 
#   Each MODS quadrant has even 'e' and odd 'o' readout chains,
#   for a total of 8 readout channels.
#
# NOTE: This version hardwires the extraction of overscan (really
#       'prescan') from full-frame 8x3K images.  A future revision
#       will use header info on the overscan columns, once the CCD
#       control software permits use of overscan (this works using the
#       relevant FITS header cards, but is currently disabled by a
#       stubborn software bug in the CCD readout program that is
#       messing up the overscan delacing on readout).
#
# Authors:
#   Rick Pogge (pogge.1@osu.edu), based on the MDM4K version by Paul Martini
#
# Module Dependencies:
#   numpy
#   astropy.io.fits
#
# Modification History
#  2011 Sep 08: initial version for just bias subtraction 
#  2011 Sep 12: tested, adapted to run on MDM computers, work on MDM4K data
#  2011 Oct 04: adapted for MODS from original proc4k.py by Paul Martini,
#               with significant changes for MODS
#  2011 Oct 21: beta-release version [rwp/osu]
#  2012 Jan 24: minor modifications [rwp/osu]
#  2012 Mar 20: added getopt command-line parsing and various updates [rwp/osu]
#  2012 Apr 01: added bad pixel fixing and other options [rwp/osu]
#  2013 Dec 01: allow multiple raw images (version 2 start) [rwp/osu]
#  2014 Aug 24: using astropy.io.fits instead of unsupported pyfits and
#               suppressing nuisance FITS header warnings [rwp/osu]
#  2016 Jul 05: version that can handle binned images, but note that binned
#               images disables the -b option [rwp/osu]
#  2017 May 18: patch for astropy deprecation of .update() [rwp/osu]
#  2017 May 21: fixed bug Barry Rothberg found out of the box [rwp/osu]
#-----------------------------------------------------------------------------

import string as str
import os 
import sys
import getopt
import numpy as np 
from astropy.io import fits

# Version number and date

versNum = '2.1.4'
versDate = '2017-05-21'

# Global Defaults

biasSingle = 0
biasRow = 1
biasFit = 2
biasType = biasSingle
outSuffix = '_otf'

avgSize = 5       # Default size of the averaging region for bad columns
modsDir = os.getenv('MODS_DBDIR','./')  # Default database directory

# Runtime flags

Verbose = False      # Verbose debugging output off
NoClobber = True     # default: no clobber on output
fixBadPix = False    # default: do not fix bad pixels
userBPL = False      # default: use the default bad pixel list
redFlip = True       # default: flip red images about X

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
  print('\nUsage: modsProc [options] rawFile(s).fits flatField.fits')
  print('\nWhere:')
  print('  rawFile(s).fits  raw MODS 8x3K FITS image(s)')
  print('  pixFlat.fits     normalized pixel flat to apply to all rawFiles')
  print('Output file will be named rawFile%s.fits' % (outSuffix))
  print('For N files listed on the command line:')
  print('  1..N-1 are the raw image files')
  print('  N is the flat field for *all* raw images')
  print('\nOptions:')
  print('  -b       fix bad pixels (see -l)')
  print('  -l file  bad pixel list to use (requires -b)')
  print('  -w nPix  width of the averaging region either side of the')
  print('            bad pixel region (requires -b).  Default: %d pixels' % (avgSize))
  print('  -f       force overwrite of an existing output file (\'clobber\')')
  print('  --noflip do not flip a Red channel image X-axis (ignored if blue)')
  print('  -v       print verbose debugging info during execution')
  print('  -V       print version info and exit\n'  )
  print('\nEnvironment:')
  print('  MODS_DBDIR = directory path to MODS database files')
  print('               currently MODS_DBDIR=%s\n' % (modsDir))

#----------------------------------------------------------------
#
# readBadPixList(file) - read contents of a bad pixel list
#
# IRAF-format bad pixel coordinate list:
#    xbegin xend ybegin yend
# comments are delimited with #
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
    print('          modsProc aborting.\n')
    sys.exit(1)

  return xs, xe, ys, ye

#----------------------------------------------------------------
#
# Main Program starts here...
#

# Parse the command-line arguments (GNU-style getopt)
  
try:
  opts, files = getopt.gnu_getopt(sys.argv[1:],'fvVbl:w:',['clobber','verbose','version',
                                                           'badpix','bpl','list','width',
                                                           'noflip',])

except getopt.GetoptError as err:
  print('\n** ERROR: %s' % (err))
  printUsage()
  sys.exit(2)

if len(opts)==0 and len(files)==0:
  printUsage()
  sys.exit(1)

# Quick Defaults as needed

numAvg = avgSize

for opt, arg in opts: 
  if opt in ('-l','--bpl','--list'):
    bplFile = arg
    userBPL = True
  elif opt in ('-b','--badpix'):
    fixBadPix = True
  elif opt in ('-w','--width'):
    numAvg = int(arg)
    if numAvg < 3:
      print('** ERROR: width of the averaging region must be 3 pixel or larger')
      sys.exit(1)
  elif opt in ('-f','--clobber'):
    NoClobber = False
  elif opt in ('--noflip'):
    redFlip = False
  elif opt in ('-v','--verbose'):
    Verbose = True
  elif opt in ('-V','--version'):
    print('modsBias.py v%s [%s]' % (versNum, versDate))
    sys.exit(0)

if len(files) < 2:
  printUsage()
  sys.exit(0)

# Number of files on the command line

numFiles = len(files)
numRaw = numFiles-1

# The flat file is the last in the line

flatFile = files[numFiles-1]

if not os.path.isfile(flatFile):
  print('Cannot find flat field FITS image %s' % (flatFile))
  printUsage()
  sys.exit(0)

# Open and load the flat field image

flatData = fits.open(flatFile,uint=False)

# Process all the raw files on the command line

for i in range(numRaw):
  rawFile = files[i]

  # Create the output filename and test for clobber/no-clobber 

  outFile = rawFile[:str.find(rawFile, '.fits')]+outSuffix+'.fits' 
  if os.path.isfile(outFile):
    if NoClobber:
      print('\n** ERROR: Operation would overwrite existing FITS file %s' % (outFile))
      print('          Use -f to force output (\'clobber\').')
      print('          modsPixFlat aborting.\n')
      sys.exit(1)
    else:
      if Verbose:
        print('** WARNING: Overwriting existing FITS file %s' % (outFile) )
      os.remove(outFile)

  print('Processing image %d of %d: %s --> %s' % (i+1,numRaw,rawFile,outFile))

  # Open the raw FITS file and extract useful info from the FITS header

  rawData = fits.open(rawFile,uint=False) 
  naxis1 = rawData[0].header['NAXIS1']
  naxis2 = rawData[0].header['NAXIS2']
  channel = rawData[0].header['CHANNEL']
  overscanx = rawData[0].header['OVERSCNX']
  overscany = rawData[0].header['OVERSCNY']	# should be 0
  detector = rawData[0].header['DETECTOR']	
  telescope = rawData[0].header['TELESCOP']
  instID = rawData[0].header['INSTRUME']
  ccdxbin = rawData[0].header['CCDXBIN']
  ccdybin = rawData[0].header['CCDYBIN']

  # If we are binned, we cannot do fixBadPix at present

  if fixBadPix and ccdxbin != 1:
    print('*** Warning: -b (bad pixel fixing) is not yet enabled for')
    print('             binned raw images.  Disabling.')
    fixBadPix = False

  # If fixing bad pixels (-b), open and read the bad pixel list

  if fixBadPix:
    if userBPL:
      pathBits = os.path.split(bplFile)
      if len(pathBits[0]) == 0:
        bplFile = os.path.join(modsDir,bplFile)
      else:
        bplFile = os.path.join(modsDir,str.lower(instID)+'.bpl')
    else:
      bplFile = os.path.join(modsDir,str.lower(instID)+'.bpl')

    if os.path.isfile(bplFile) == 0:
      print('\n** ERROR: MODS bad pixel list file %s not found.' % (bplFile))
      print('          modsProc aborting.\n')
      sys.exit(1)

    # read and parse the contents of bplFile
  
    xs, xe, ys, ye = readBadPixList(bplFile)

    if Verbose:
      print('Bad pixel list %s contains %d bad pixel regions to fix' % (bplFile,len(xs)+1))

  # Stupid bit of hardwiring for now, until the DOS IC code can read
  # X and Y overscan pixels for any subframe ROI without crashing
  # and burning... [rwp]
    
  if naxis1==8288:
    overscanx = 48
    overscany = 0
  elif naxis1==4144 and ccdxbin==2:
    overscanx = 24
    overscany = 0
  else:
    print('** ERROR: NAXIS1=%d, modsProc only implemented for full-frame readout' % (naxis1))
    print('          at the present time, sorry.')
    sys.exit(1)
      
  if Verbose:
    print('Processing %s[%d:%d] OVERSCANX=%d OVERSCANY=%d' % (rawFile, naxis1, naxis2, overscanx, overscany))
    print('Channel: %s' % (channel))
    print('Detector: %s' % (detector))
    print('Telescope: %s' % (telescope))
    
  if overscanx < (32/ccdxbin): 
    print('** ERROR: NAXIS1=%d OVERSCNX=%d less than %d in %s' % ((32/ccdxbin),naxis1,overscanx,rawFile))
    print('          modsProc only works with full-frame images for now')
    sys.exit(1)
  if overscany > 0: 
    print('Error: code not tested with OVERSCNY > 0!' )
    sys.exit(1)

  # Channel-specific setup, if needed, would go here...

  # CCD image 'metric' parameters
    
  c1 = overscanx 		# first image column counting from *zero*
  c2 = int(0.5*naxis1)-1	# last image column on first half 
  c3 = c2+1			# first image column on second half 
  c4 = naxis1-overscanx-1 	# last image column 
  r1 = overscany 		# first image row 
  r2 = int(0.5*naxis2)-1	# last image row on first half 
  r3 = r2+1			# first image row on second half 
  r4 = naxis2-overscany-1  	# last image row 
  outnaxis1 = c4-c1+1		# columns in output, trimmed image 
  outnaxis2 = r4-r1+1		# rows in output, trimmed image
  collen = int(0.5*outnaxis1)	# number of rows in an image quadrant
  rowlen = int(0.5*outnaxis2)	# number of rows in an image quadrant
  
  if Verbose:
    print('Quadrants:')
    print('  Q1: [%d:%d,%d:%d]' % (c1, c2, r1, r2) )
    print('  Q2: [%d:%d,%d:%d]' % (c3, c4, r1, r2) )
    print('  Q3: [%d:%d,%d:%d]' % (c1, c2, r3, r4) )
    print('  Q4: [%d:%d,%d:%d]' % (c3, c4, r3, r4) )

  # Calculate the bias level for each amplifier

  data = rawData[0].data

  # Identify the columns to use to calculate the bias level 
  # Skip the first and last columns of the overscan 

  # Identify the even and odd columns in the overscan

  cols_over_q1e = np.arange(4, overscanx-4, 2)
  cols_over_q1o = np.arange(5, overscanx-4, 2)
  cols_over_q2e = np.arange(naxis1-overscanx+4, naxis1-4, 2)
  cols_over_q2o = np.arange(naxis1-overscanx+5, naxis1-4, 2)
  cols_over_q3e = cols_over_q1e
  cols_over_q3o = cols_over_q1o 
  cols_over_q4e = cols_over_q2e
  cols_over_q4o = cols_over_q2o 

  # Identify the even and odd columns in each quadrant 

  cols_q1e = np.arange(c1,c2,2)
  cols_q1o = np.arange(c1+1,c2+2,2)
  cols_q2e = np.arange(c3,c4,2)
  cols_q2o = np.arange(c3+1,c4+2,2)
  cols_q3e = cols_q1e
  cols_q3o = cols_q1o
  cols_q4e = cols_q2e
  cols_q4o = cols_q2o

  # Create arrays with the median overscan vs. row for each amplifier
    
  bias_q1e = np.zeros(rowlen,dtype=float)
  bias_q1o = np.zeros(rowlen,dtype=float)
  bias_q2e = np.zeros(rowlen,dtype=float)
  bias_q2o = np.zeros(rowlen,dtype=float)
  bias_q3e = np.zeros(rowlen,dtype=float)
  bias_q3o = np.zeros(rowlen,dtype=float)
  bias_q4e = np.zeros(rowlen,dtype=float)
  bias_q4o = np.zeros(rowlen,dtype=float)

  # Calculate 1-D bias arrays for each amplifier
    
  for i in range(r1, r2+1, 1): 
    bias_q1e[i] = np.median(data[i,cols_over_q1e]) 	# data[rows, columns]
    bias_q1o[i] = np.median(data[i,cols_over_q1o]) 	
    bias_q2e[i] = np.median(data[i,cols_over_q2e])
    bias_q2o[i] = np.median(data[i,cols_over_q2o])
    bias_q3e[i] = np.median(data[i+rowlen,cols_over_q3e])
    bias_q3o[i] = np.median(data[i+rowlen,cols_over_q3o])
    bias_q4e[i] = np.median(data[i+rowlen,cols_over_q4e])
    bias_q4o[i] = np.median(data[i+rowlen,cols_over_q4o])

  # Subtract the bias from the output

  OverscanKeyValue = 'biasSingle' 
    
  # Subtract a single bias value for each amplifier
    
  bq1e = np.median(bias_q1e) 
  bq1o = np.median(bias_q1o) 
  bq2e = np.median(bias_q2e) 
  bq2o = np.median(bias_q2o) 
  bq3e = np.median(bias_q3e) 
  bq3o = np.median(bias_q3o) 
  bq4e = np.median(bias_q4e) 
  bq4o = np.median(bias_q4o) 

  data[r1:r2+1,cols_q1e] -= bq1e
  data[r1:r2+1,cols_q1o] -= bq1o
  data[r1:r2+1,cols_q2e] -= bq2e
  data[r1:r2+1,cols_q2o] -= bq2o
  data[r3:r4+1,cols_q3e] -= bq3e
  data[r3:r4+1,cols_q3o] -= bq3o
  data[r3:r4+1,cols_q4e] -= bq4e
  data[r3:r4+1,cols_q4o] -= bq4o

  # Write the output file to the FITS array

  rawData[0].data = data[r1:r4+1,c1:c4+1] 
  OverscanKeyComment = 'Overscan subtracted by modsProc.py' 

  rawData[0].header.add_history('modsProc v%s %s' % (versNum,versDate))
  rawData[0].header['BIASPROC'] = (OverscanKeyValue,OverscanKeyComment) 
  rawData[0].header['BIASQ1E'] = (bq1e,'Bias subtracted from Q1E') 
  rawData[0].header['BIASQ1O'] = (bq1o,'Bias subtracted from Q1O') 
  rawData[0].header['BIASQ2E'] = (bq2e,'Bias subtracted from Q2E') 
  rawData[0].header['BIASQ2O'] = (bq2o,'Bias subtracted from Q2O') 
  rawData[0].header['BIASQ3E'] = (bq3e,'Bias subtracted from Q3E') 
  rawData[0].header['BIASQ3O'] = (bq3o,'Bias subtracted from Q3O') 
  rawData[0].header['BIASQ4E'] = (bq4e,'Bias subtracted from Q4E') 
  rawData[0].header['BIASQ4O'] = (bq4o,'Bias subtracted from Q4O') 

  # We need a working copy of the data from here out

  workData = rawData[0].data

  # Divide the bias-corrected image by the flat field

  if Verbose:
    print('Dividing bias-corrected %s by flat field %s...' % (rawFile,flatFile))

  workData /= flatData[0].data
  rawData[0].header.add_history('modsProc - divided by flat field')
  rawData[0].header['FLATPROC'] = (flatFile, 'Flat field image') 

  if fixBadPix:
    # Page through the bad pixel list and fix the pixel
    rawData[0].header.add_history('modsProc - fixing bad pixels')
    rawData[0].header['BPLFILE'] = (bplFile)
    rawData[0].header['BPWIDTH'] = (numAvg,'Bad pixel averaging region width [pixels]')

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
        medLeft = np.median(workData[y,xsl:xel])
        medRight= np.median(workData[y,xsr:xer])
        corrVal = (medLeft+medRight)/2.0
        if (nx>1):
          workData[y,xstart:xend] = corrVal
        else:
          workData[y,xstart] = corrVal

  # If this is a red-channel image, we need to flip the red channel into
  # blue-to-red wavelength order

  if str.upper(channel)=='RED':
    if redFlip:
      if Verbose:
        print('Flipping Red-channel CCD into blue-to-red wavelength order...')
      workData = np.fliplr(workData)
      rawData[0].header.add_history('modsProc - flipped image about the X-axis')
    else:
      if Verbose:
        print('Red-channel CCD images were not flipped in X (--noflip)...')

  # Write the processed image data to disk

  rawData[0].data = workData

  if Verbose:
    print('Writing bias- and flat-corrected %s --> %s' % (rawFile, outFile))

  rawData.writeto(outFile,output_verify='ignore')

  # done close the fits file

  rawData.close()

# Bottom of the loop, all done

if Verbose:
  print('modsProc Done')

sys.exit(0)
