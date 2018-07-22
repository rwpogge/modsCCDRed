#! /usr/bin/env python
#
# modsBias - correct MODS CCD images for bias and even/odd
#            readout artifacts
#
# Usage: modsBias [options] rawFile1.fits [rawFile2.fits...]
# 
# Where:
#   rawFile?.fits  raw MODS 8x3K FITS files
#
# Bias-corrected output files will be named rawFile1_ot.fits, etc.
# 
# Options:
#   -f      force overwrite of an existing output file ('clobber')
#   -v      print verbose debugging info during execution
#   -V      print version info and exit
#
# Description:
#
# Perform the overscan subtraction and remove even-odd bias 
# differences for a single full-frame MODS 8x3K image.
#
# Steps: 
#  1. Determine if input is a file or list of files
#  2. Identify binning, size of overscan
#  3. Remove overscan and trim
#
# Experiments to date show that attempting to empirically correct for
# relative gain variations between the 8 readout chains is not
# indicated: they come out as hoped in the pixel-to-pixel flats formed
# from the slitless spectral flats.
#
# Quadrant Layout of the MODS Science CCDs:
#
#      ++--------+--------++(NAXIS1,NAXSI2)
#      ||   Q3   |   Q4   ||
#      ++--------+--------++
#      ||   Q1   |   Q2   ||
# (1,1)++--------+--------++
# 
# Each MODS quadrant has even 'e' and odd 'o' readout chains, for a
# total of 8 readout channels.
#
# NOTE: This version hardwires the extraction of overscan (really
#       'prescan') from full-frame 8x3K images.  A future revision
#       will use header info on the overscan columns, once the CCD
#       control software permits use of overscan (this is currently
#       disabled by a stubborn software bug).
#
# Authors:
#   Rick Pogge (pogge.1@osu.edu), based on an MDM4K CCD version by
#   Paul Martini
#
# Module Dependencies:
#   numpy
#   astropy.io.fits
#
# Modification History
#  2011 Sep 08: Initial version for just bias subtraction.
#  2011 Sep 12: Tested, adapted to run on MDM computers, work on MDM4K data.
#  2011 Oct 04: Adapted for MODS from original proc4k.py by Paul Martini,
#               with significant changes for MODS.
#  2011 Oct 21: Beta-release version. [rwp/osu]
#  2012 Mar 30: Added getopt command-line handling and new options. [rwp/osu]
#  2013 Dec 01: Made quadrant numbering consistent with the MODS Instrument
#               manual throughout the program not just in the BIASQnE/O 
#               header cards at the end. [rwp/osu]
#  2014 Aug 24: Using anaconda, and astropy.io to replace the unsupported
#               pyfits and suppressing nuisance FITS header warnings [rwp/osu]
#  2016 Mar 08: Release of version that can handled binned images [rwp/osu]
#  2016 Jul 05: Minor tweak to error message (inconsistent info) [rwp/osu]
#  2017 May 17: Replaces all instances of .update() in FITS header
#               handling with the [] syntax.  Recent astropy.io.fits update
#               eliminated the .update() method, but [] is back compatible
#               for all versions [rwp/osu]
#-----------------------------------------------------------------------------

import string as str
import os 
import sys
import getopt
import numpy as np 
from astropy.io import fits

# Version number and date

versNum = '2.1.3'
versDate = '2017-05-17'

# Global Defaults

biasSingle = 0
biasRow = 1
biasFit = 2
biasType = biasSingle
outSuffix = '_ot'

# Runtime flags

Verbose = False   # Verbose debugging output off
NoClobber = True  # default: no clobber on output

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
  print('\nUsage: modsBias [options] rawFile1.fits [rawFile2.fits...]')
  print('\nWhere:')
  print('  rawFile?.fits  raw MODS 8x3K FITS files')
  print('                 output files will be rawFile1%s.fits, etc.' % (outSuffix))
  print('\nOptions:')
  print('  -f      force overwrite of an existing output file (\'clobber\')')
  print('  -v      print verbose debugging info during execution')
  print('  -V      print version info and exit\n')

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
    print('modsBias.py v%s [%s]' % (versNum, versDate))
    sys.exit(0)

if len(files) < 1:
  printUsage()
  sys.exit(0)

# Loop over all the files and process one by one

for rawFile in files: 
  if os.path.isfile(rawFile):

    # Check the output file for clobber/noclobber
    
    outFile = rawFile[:str.find(rawFile, '.fits')]+outSuffix+'.fits' 
    if os.path.isfile(outFile):
      if NoClobber:
        print('\n** ERROR: Operation would overwrite existing FITS file %s' % (outFile))
        print('          Use -f to force output (\'clobber\').')
        print('          modsPixFlat aborting.\n')
        sys.exit(1)
      else:
        if Verbose:
          print('** WARNING: Overwriting existing FITS file %s' % (outFile))
        os.remove(outFile)

    # Open the FITS file and extract useful info from the header
    
    fitsFile = fits.open(rawFile,uint=False) 
    naxis1 = fitsFile[0].header['NAXIS1']
    naxis2 = fitsFile[0].header['NAXIS2']
    channel = fitsFile[0].header['CHANNEL']
    overscanx = fitsFile[0].header['OVERSCNX']
    overscany = fitsFile[0].header['OVERSCNY']	# should be 0
    detector = fitsFile[0].header['DETECTOR']	
    telescope = fitsFile[0].header['TELESCOP'] 
    ccdxbin = fitsFile[0].header['CCDXBIN']
    ccdybin = fitsFile[0].header['CCDYBIN']

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
      print('\n** ERROR: NAXIS1=%d, modsBias only implemented for full-frame readout' % (naxis1))
      print('          at the present time, sorry.')
      sys.exit(1)
      
    if Verbose:
      print('Processing %s[%d:%d] OVERSCANX=%d OVERSCANY=%d' % (rawFile, naxis1, naxis2, overscanx, overscany))
      print('Channel: %s' % (channel))
      print('Detector: %s' % (detector))
      print('Telescope: %s' % (telescope))
    
    if overscanx < (32/ccdxbin): 
      print('\n** ERROR: NAXIS1=%d OVERSCNX=%d less than %d in %s' % ((32/ccdxbin),naxis1,overscanx,rawFile))
      print('          modsBias only works with full-frame images for now.')
      sys.exit(1)
    if overscany > 0: 
      print('\n** Error: code not tested with OVERSCNY > 0!' )
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

    data = fitsFile[0].data

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
    
    bias_q1e = np.zeros(rowlen, dtype=float)
    bias_q1o = np.zeros(rowlen, dtype=float)
    bias_q2e = np.zeros(rowlen, dtype=float)
    bias_q2o = np.zeros(rowlen, dtype=float)
    bias_q3e = np.zeros(rowlen, dtype=float)
    bias_q3o = np.zeros(rowlen, dtype=float)
    bias_q4e = np.zeros(rowlen, dtype=float)
    bias_q4o = np.zeros(rowlen, dtype=float)

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

    overscanKeyValue = 'biasSingle' 
    
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

    # Transfer the data and update the FITS headers

    fitsFile[0].data = data[r1:r4+1,c1:c4+1]

    overscanKeyComment = 'Overscan subtracted by modsBias.py' 
    fitsFile[0].header.add_history('modsBias v%s %s' % (versNum,versDate))
    fitsFile[0].header['BIASPROC'] = (overscanKeyValue, overscanKeyComment) 
    fitsFile[0].header['BIASQ1E'] = (bq1e, 'Bias subtracted from Q1E') 
    fitsFile[0].header['BIASQ1O'] = (bq1o, 'Bias subtracted from Q1O') 
    fitsFile[0].header['BIASQ2E'] = (bq2e, 'Bias subtracted from Q2E')
    fitsFile[0].header['BIASQ2O'] = (bq2o, 'Bias subtracted from Q2O') 
    fitsFile[0].header['BIASQ3E'] = (bq3e, 'Bias subtracted from Q3E')
    fitsFile[0].header['BIASQ3O'] = (bq3o, 'Bias subtracted from Q3O') 
    fitsFile[0].header['BIASQ4E'] = (bq4e, 'Bias subtracted from Q4E') 
    fitsFile[0].header['BIASQ4O'] = (bq4o, 'Bias subtracted from Q4O') 

    # Write the output file and close the original FITS file
    
    fitsFile.writeto(outFile,output_verify='ignore')
    print('Processed image %s --> %s' % (rawFile, outFile))

    fitsFile.close()

if Verbose:
  print('modsBias Done')

sys.sys.exit(0)
