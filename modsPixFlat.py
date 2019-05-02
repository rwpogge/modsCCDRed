#! /usr/bin/env python
#
# modsPixFlat - create a normalized pixel flat from bias-corrected
#               MODS spectral flats.
#
# Usage: modsPixFlat inFlat pixFlat
# 
# Where:
#   inFlat  = input bias-corrected flat field image
#   pixFlat = output normalized pixel flat to create
# 
# Options:
#   -f      force overwrite of an existing pixFlat ('clobber')
#   -v      print verbose debugging info during execution
#   -V      print version info and exit
#
# Description:
#
#   Creates a normalized, color-free pixel flat from a bias-corrected
#   MODS spectral flat field image.
# 
#   In each quadrant, it averages along columns to create a mean color
#   vector.  It then smooths that vector by a 2-pixel boxcar filter to
#   avoid amplifying the even/odd effect in the MODS readout, and then
#   divides each line of the quadrant by the mean, smoothed color
#   vector.
#
#   To preserve the inter-quadrant normalization, it measures the
#   relative quadrant gain along the boundaries between them from the
#   original 2D spectrum, and then uses these data to rescale the
#   color-free pixel flat to restore the relative gains between the
#   quadrants.
#
# Coordinate System:
#      +---------------+---------------+
#      |               |               |
#      |        3      |       4       |
#      |               |               |
#      +---------------+---------------+
#      |               |               |
#      |       1       |       2       |
#      |               |               |
#      +---------------+---------------+
#
# Author:
#   R. Pogge, OSU Astronomy Dept
#   pogge.1@osu.edu
#
# Module Dependencies:
#   numpy
#   astropy.io fits
#
# Modification History
#   2012 Mar 26 - added boxcar filter of the color term vector [rwp/osu]
#   2012 Mar 30 - added getopt features [rwp/osu]
#   2013 Feb 25 - new normalization scheme [rwp/osu]
#   2014 Aug 24 - using astropy.io.fits instead of unsupported pyfits and
#                 suppressing nuisance FITS header warnings [rwp/osu]
#   2017 May 18 - patch for .update() deprecation in astropy [rwp/osu]
#
#-----------------------------------------------------------------------------

import string as str
import os 
import sys
import getopt
import numpy as np 
from astropy.io import fits

# Version number and date

versNum  = '2.0.10'
versDate = '2019-05-02'

# Global Defaults

stripeWid = 10   # normalization stripe width in pixels
stripeLen = 100  # normalization stripe length in pixels
stripeOff = 1    # stripe standoff from border in pixels
colorOff  = 150  # offset from top/bottom for color term extraction

# Runtime flags

Verbose   = False    # Verbose debugging output on/off
NoClobber = True     # default: no clobber on output

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
  print('\nUsage: modsPixFlat inFlat pixFlat')
  print('\nWhere:')
  print('  inFlat  = input bias-corrected flat field image')
  print('  pixFlat = output normalized pixel flat to create')
  print('\nOptions:')
  print('  -o ###  offset from top/bottom edges in pixels (default: %d)' % (colorOff))
  print('  -f      force overwrite of an existing pixFlat (\'clobber\')')
  print('  -v      print verbose debugging info during execution')
  print('  -V      print version info and exit\n')
  
#----------------------------------------------------------------
#
# Main Program starts here...
#

# Parse the command-line arguments (GNU-style getopt)
  
try:
  opts, files = getopt.gnu_getopt(sys.argv[1:],'fvVo:',
                                  ['clobber','verbose','version','offset='])
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
  elif opt in ('-o','--offset'):
    colorOff = int(arg)
  elif opt in ('-V','--version'):
    print('modsPixFlat.py v%s [%s]' % (versNum, versDate))
    sys.exit(0)

numFiles = len(files)

if numFiles < 2:
  printUsage()
  sys.exit(1)
  
inFile = files[0]
outFile = files[1]

# Verify files exist, abort if input file not found.

if not os.path.isfile(inFile):
  print('** ERROR: Input FITS file %s not found' % (inFile))
  print('          modsPixFlat aborting.')
  sys.exit(1)

# Do not clobber an existing output file unless -f is used

if NoClobber and os.path.isfile(outFile):
  print('\n** ERROR: Operation would overwrite existing FITS file %s' % (outFile))
  print('          Use -f to force output (\'clobber\').')
  print('          modsPixFlat aborting.\n')
  sys.exit(1)
  
# Read in the un-normalized flat field.

fitsFile = fits.open(inFile,uint=False)

# Get possibly useful FITS header info...

channel = fitsFile[0].header['CHANNEL']
naxis1  = fitsFile[0].header['NAXIS1']
naxis2  = fitsFile[0].header['NAXIS2']

# Read in the FITS data array
  
if Verbose:
  print('reading %s...' % (inFile))

data = fitsFile[0].data

# These are the starting pixels of each quadrant

nxquad = int(naxis1/2)
nyquad = int(naxis2/2)
(xs1,ys1) = (0,0)
(xs2,ys2) = (nxquad,0)
(xs3,ys3) = (0,nyquad)
(xs4,ys4) = (nxquad,nyquad)

# These are the coordinates of the 4-corners pixels for
# each quadrant used to compute the stripe vectors for
# measuring the relative inter-quadrant gain levels.

(x1,y1)=(nxquad-1,nyquad-1)
(x2,y2)=(nxquad,  nyquad-1)
(x3,y3)=(nxquad-1,nyquad)
(x4,y4)=(nxquad,  nyquad)

# Vertical Normalization Stripes

v13 = np.mean(data[...,x1-stripeWid:x1],axis=1)
v24 = np.mean(data[...,x2:x2+stripeWid],axis=1)
vrat = v24/v13  # vertical stripe ratio 

r21 = np.median(vrat[y1-stripeOff-stripeLen:y1-stripeOff])
r43 = np.median(vrat[y3+stripeOff:y3+stripeOff+stripeLen])

if Verbose:
  print('r21=%.5f r43=%.5f' % (r21,r43))

# Horizontal Normalization Stripes

h12 = np.mean(data[y1-stripeWid:y1,...],axis=0)
h34 = np.mean(data[y3:y3+stripeWid,...],axis=0)
hrat = h34/h12  # horizontal stripe ratio

r31 = np.median(hrat[x1-stripeOff-stripeLen:x1-stripeOff])
r42 = np.median(hrat[x2+stripeOff:x2+stripeOff+stripeLen])

if Verbose:
  print('r31=%.5f r42=%.5f' % (r31,r42))

# Quadrant scaling factors

s1 = 1.0
s2 = r21
s3 = r31
s4 = r43*r31
norm=(s1+s2+s3+s4)/4

s1 /= norm
s2 /= norm
s3 /= norm
s4 /= norm

if Verbose:
  print('   Quadrant Scaling Factors: s1=%.5f s2=%.5f s3=%.5f s4=%.5f' % (s1,s2,s3,s4))

# For each quadrant, create a color term vector by computing the
# mean along columns.

if Verbose:
  print('Removing the color terms...')
  print('   Quadrant 1: [%d:%d,%d:%d]' % (xs1,nxquad-1,ys1,nyquad-1))

cTerm = np.mean(data[ys1+colorOff:nyquad,xs1:nxquad],axis=0)
cTerm += np.roll(cTerm,1)
cTerm /= 2.0
data[ys1:nyquad,xs1:nxquad] /= cTerm
if Verbose:
  print('   Quadrant 2: [%d:%d,%d:%d]' % (xs2,naxis1-1,ys2,nyquad-1))

cTerm = np.mean(data[ys2+colorOff:nyquad,xs2:naxis1],axis=0)
cTerm += np.roll(cTerm,-1)
cTerm /= 2.0
data[ys2:nyquad,xs2:naxis1] /= cTerm

if Verbose:
  print('   Quadrant 3: [%d:%d,%d:%d]' % (xs3,nxquad-1,ys3,naxis2-1))

cTerm = np.mean(data[ys3:naxis2-colorOff,xs3:nxquad],axis=0)
cTerm += np.roll(cTerm,1)
cTerm /= 2.0
data[ys3:naxis2,xs3:nxquad] /= cTerm

if Verbose:
  print('   Quadrant 4: [%d:%d,%d:%d]' % (xs4,naxis1-1,ys4,naxis2-1))

cTerm = np.mean(data[ys4:naxis2-colorOff,xs4:naxis1],axis=0)
cTerm += np.roll(cTerm,-1)
cTerm /= 2.0
data[ys4:naxis2,xs4:naxis1] /= cTerm

# If any color term vector pixels were 0, this produces an infinity in
# corresponding pixels in the flat.  Set those pixels to 1.

if np.any(data==float('inf')):
  data[data==float('inf')] = 1.0

# Also take out 0 pixels in the flat to avoid divergent pixels when used

if np.any(data==0.0):
  data[data==0.0] = 1.0;

# Recompute the inter-quadrant scalings after removing the color term

vc13 = np.mean(data[...,x1-stripeWid:x1],axis=1)
vc24 = np.mean(data[...,x2:x2+stripeWid],axis=1)
vcrat = vc24/vc13  # vertical stripe ratio 

rc21 = np.median(vcrat[y1-stripeOff-stripeLen:y1-stripeOff])
rc43 = np.median(vcrat[y3+stripeOff:y3+stripeOff+stripeLen])

if Verbose:
  print('rc21=%.5f rc43=%.5f' % (rc21,rc43))

hc12 = np.mean(data[y1-stripeWid:y1,...],axis=0)
hc34 = np.mean(data[y3:y3+stripeWid,...],axis=0)
hcrat = hc34/hc12  # horizontal stripe ratio

rc31 = np.median(hcrat[x1-stripeOff-stripeLen:x1-stripeOff])
rc42 = np.median(hcrat[x2+stripeOff:x2+stripeOff+stripeLen])

if Verbose:
  print('rc31=%.5f rc42=%.5f' % (rc31,rc42))

# Post-color term quadrant scaling factors

sc1 = 1.0
sc2 = rc21
sc3 = rc31
sc4 = rc43*rc31

norm=(sc1+sc2+sc3+sc4)/4

sc1 /= norm
sc2 /= norm
sc3 /= norm
sc4 /= norm

if Verbose:
  print('sc1=%.5f sc2=%.5f sc3=%.5f sc4=%.5f' % (sc1,sc2,sc3,sc4))

# Now reconcile the inter-quadrant scalings

c1=s1/sc1
c2=s2/sc2
c3=s3/sc3
c4=s4/sc4

if Verbose:
  print('c1=%.5f c2=%.5f c3=%.5f c4=%.5f' % (c1,c2,c3,c4))

# Now reconcile the inter-quadrant scalings

data[ys1:nyquad,xs1:nxquad] *= c1
data[ys2:nyquad,xs2:naxis1] *= c2
data[ys3:naxis2,xs3:nxquad] *= c3
data[ys4:naxis2,xs4:naxis1] *= c4

# Write the output file

if Verbose:
  print('Writing normalized pixel flat to %s...' % (outFile))

fitsFile[0].header.add_history('modsPixFlat v%s %s orig=%s' % (versNum,versDate,inFile))
fitsFile[0].data = data

if os.path.isfile(outFile): 
  if Verbose:
    print('** WARNING: Overwriting existing FITS file %s' % (outFile) )
  os.remove(outFile) 

fitsFile.writeto(outFile,output_verify='ignore')
    
if Verbose:
  print('modsPixFlat Done')

sys.exit(0)
