# MODS Bad Pixel Tables

These files are used by [modsFixPix.py](../modsFixPix.py) to repair
known CCD artifacts, mostly clusters of blocked columns.

The files lists the bad pixel regions as starting/ending x and y
pixels of the region corners in full-frame (8192x3088) unbinned pixel
coordinates:
<pre>
    xs xe ys ye
</pre>
one region per line.

Files may be updated if new features appear on the CCDs (or we find
regions that were overlooked).
