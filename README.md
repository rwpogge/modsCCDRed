# modsCCDRed
MODS Basic 2D CCD Reduction Programs

## Overview

modsCCDRed is a collection of Python programs to perform basic 2D
reduction (bias, trim, flat field, and bad column fixing) on raw MODS
CCD data.

## System Requirements

The modsCCDRed Python programs have been tested on Linux (CentOS 5 to
7) and Mac OS X (v10.7 Lion thru 10.13 HighSierra) operating systems.

Version 2.0.5 of these scripts were ported to Python 3 in a way that
should still be back-compatible with people still running 
Python 2.7. We have adopted the free Anaconda Python for development of these program as it has
proven to be the most trouble free distribution across Linux and Mac
platforms, and we are becoming more invested in astropy across the
profession for utilities like this.

If you do not use anaconda python, you will need to load and install
the numpy and astropy modules before you can use modsCCDRed.

Your computer should have at least 4 Gb of memory in order to be able
to perform the image median stacking steps, but otherwise the programs
work with single images, but more is better.  Runtime performance will
depend on memory and processor speed, and is somewhat I/O limited
given the large (24 Mpixel) images involved.  On an older 2-processor
2.8GHz Pentium 4 computer with 8Gb of RAM running CentOS 5.8, the time
required to bias, flat-field, bad column fix and axis flip a single
full-frame unbinned 8x3K MODS red-channel spectrum is about 13
seconds.

## Auxiliary files

The database/ subdirectory includes a set of default bad pixel list
(bpl) files for the MODS CCDs.

The last version of the manual is in modsCCDRed.pdf, but it predates
the upload to github, so beware.

## Installation

You have two options: Personal or Public installation

For a personal installation, you have two ways to install the files to use:

1. Keep the modsCCDRed Python programs in place, and put the
modsCCDRed directory in your default execution path.

2. Copy the executable Python programs into your ~/bin/ or ~/programs/
directory where you put executables in your default execution path.

For a public installation, e.g., on a central disk to share one copy
of the package among many users, we suggest logging in as root and
unpacking the tarball in your usual place for public add-on programs,
e.g., /usr/local/, and then installing the executables in, e.g.,
/usr/local/bin/.

## Common Bad Pixel Lists

To use the [default bad pixel lists](database), users need to define this
environment variable

 > setenv MODS_DBDIR /usr/local/LBT/modsCCDRed/database/

for a csh/tcsh shells with a public installation where the source code
lives in /usr/local/LBT/.  Other shells (e.g., bash) use a different
syntax.

Users can always override the default bad pixel list used with one of
these methods:

1. Redefining the MODS_DBDIR environment variable to point to the
   directory with their personal copies of bpl files with the same
   names (e.g., mods1r.bpl).

2. Undefining MODS_DBDIR, at which point the modsCCDRed programs
   default to the current working directory (./), and look for the
   standard names.

3. Using the -l flag with no path to use a custom file in
   MODS_DBDIR (or ./ if undefined), e.g., "modsFixPixl -l mym1r.bpl"

4. Using the -l flag to give the full path and name of a particular
   bad pixel list file to use, e.g., "modsFixPix -l
   /my/path/to/mym1r.bpl"

There is no requirement of using the MODS_DBDIR environment variable,
but care should be taken to avoid having multiple, possibly
conflicting versions of bad pixel lists spread through many
directories.

## Python 3 Compatibility

We have a verified Python 3 compatible version starting with v2.0.5. This version
is also verified to run under Python 2.7 for those who still haven't made the
leap, but beware that external modules like numpy et al. may not be so forgiving
moving forward, so we recommend you switch to Python 3 if you haven't already.

## Acknowledging modsCCDRed

If your research used modsCCDRed, we ask that you follow emerging
[software citation principles](https://doi.org/10.7717/peerj-cs.86) being adopted by the astronomical community
to ensure the proper citation of software in scientific publications. 

 > [DOI 10.5281/zenodo.2550740](https://doi.org/10.5281/zenodo.2550740)

If you publish MODS data reduced with the modsCCDRed, we ask that you add this line to the acknowledgments section of
your paper:

    This paper made use of the modsCCDRed data reduction code developed in part with funds provided by NSF Grants 
    AST-9987045 and AST-1108693. 
    
modsCCDRed was developed for the MODS1 and MODS2 instruments at the Large Binocular Telescope Observatory, which were
built with major support provided by grants from the U.S. National Science Foundation's Division of Astronomical Sciences
Advanced Technologies and Instrumentation (AST-9987045), the NSF/NOAO TSIP Program, and matching funds provided by the Ohio
State University Office of Research and the Ohio Board of Regents.  Additional support for modsCCDRed was provided by NSF
Grant AST-1108693.

## Release Notes

<dl>
<dt>2019 Apr 19 - v2.0.5 - bug fixes
<dd>Fixed problems found in Python 2 to Python 3 conversion that will cause
problems with some Python 2.7 distros
 
 <dt>2018 Oct 24 - v2.0.4 - bug fixes
<dd>Fixed stupid errors in modsBias.py and modsFixPix.py
 
<dt>2017 May 21 - v2.0.3 - astropy patch
<dd>A recent update to the astropy.io FITS module removed support for
the .update() method for FITS headers, but a backwards-compatible
syntax is available that works regardless of the astropy.io update.
Too many cooks in astropy.io FITS of late ...

<dt>2016 Mar 09 - v2.0.1 - First Public Release
<dd>Release of changes in modsPixFlat in general

<dd>Adoption of astropy.io.fits for FITS handling with the
announcement by STScI of end of support for PyFITS (with the wholesale
move of PyFITS development to AstroPy, so we're just moving with it).

<dd>Suppression of nuisance FITS warnings, mostly because the AstroPy
FITS implementation went a more orthodox route regarding FITS
compliance than PyFITS.

<dd>Many small changes fixing bugs and responding to user requests.

<dd>Fixes an issues with astropy v1.1.0 and later in that in reading
integer FITS the fits.open() function now defaults to retaining
integer format instead of the desired (old) default behavior of
converting to floats.

</dl>

## Author & History

modsCCDRed was written by Rick Pogge, The Ohio State University
Department of Astronomy (pogge.1@osu.edu), including bias correction
programs developed for OSMOS by Paul Martini (OSU Astronomy).
