#!/usr/bin/env python

# Add the ZPN projection WCS headers to an all sky image

# Example usage:
# ./zpn_wcs.py image.fits new.fits
#
# This page was helpful
# https://docs.astropy.org/en/stable/wcs/example_create_imaging.html
#
# New version 6/3/2024

import sys
#import os
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
#from astropy.wcs import WCS # changed this to
#from astropy.wcs import WCS
from astropy import wcs
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy import units as u
from astropy.time import Time
import math

# Check input parameters
if ( len(sys.argv) != 3):
    print("There is a problem with your arguments.")
    print("Argument 1 should be an existing fits file.")
    print("Argument 2 should be a new fits file.")
    sys.exit('Parameters are wrong')

# If all OK, then first argument is the file to open
imagefile = sys.argv[1]
print("Opening ",imagefile)

# Open image
hdu = fits.open(imagefile)

# Assume we want the 0th header unit
hdr = hdu[0].header

# Get data array, again fron fits extension 0 
data_array = fits.getdata(imagefile,ext=0)

# Close input image
hdu.close()

# See if there was a WCS header in the original image and print it.
# There probably shouldn't be a WCS header, so all values will be zero.
wcs_orig = wcs.WCS(hdr)
print("Original WCS is \n",wcs_orig)

# Get some header values
obs_date = hdu[0].header['DATE']
obsgeob = hdu[0].header['OBSGEO-B']
obsgeol = hdu[0].header['OBSGEO-L']

# Create time object from the OBS_DATE header
obs_time = Time(obs_date, format='isot', scale='utc')

# Calculate a modified Juliad Date (MJD)
mjd = obs_time.mjd
print("MJD is ",mjd)

# Get camera latitude and longitude and print them
lat_degrees = float(obsgeob)
lon_degrees = float(obsgeol)

print("latitude is ",lat_degrees)
print("longitude is ",lon_degrees)

# Create a new WCS object.  The number of axes must be set
# from the start
print("\nCreating new ZPN WCS\n")

#new_wcs = WCS(hdu)
new_wcs = wcs.WCS(naxis=2)

# Calculate plate scale using distance from polaris to Zenith
#x_polaris =1456.0
#y_polaris = 1728.0

# These are the x,y values of the centre of the image, the zenith point
# 1599.55, 1060.60 from George project.
# Values of 1598.05, 1058.79 from Queenette and Jamie
xcentre = 1598.05
ycentre = 1058.79

# Set CRPIX WCS header using centre of the image as reference pixels
new_wcs.wcs.crpix = [xcentre,ycentre]

print("Plate scale is 3.3 arcmin per pixel or 0.055 degrees per pixel")
xplate_scale = 0.055
yplate_scale = 0.055

# Set CRDELT WCS header
new_wcs.wcs.cdelt = [xplate_scale, yplate_scale]

# Need to calculate RA of zenith, which is the same as the LST.
# Local Siderial Time which is the RA of an object that is due South.
LST = obs_time.sidereal_time('apparent', lon_degrees)
print("LST is ",LST.degree)

# Set CRVAL and CTYPE WCS header values
new_wcs.wcs.crval = [LST.degree,lat_degrees]
new_wcs.wcs.ctype = ["RA---ZPN", "DEC--ZPN"]

# Set CD matrix here for a clockwise rotation of 9 degrees
rot_angle = -12
a = math.sin(math.radians(rot_angle))
b = math.cos(math.radians(rot_angle))
new_wcs.wcs.pc=[(b,-a,), (a ,b )]

# Need the line below or geterror and no file created.
# Set at least one of the polynomial distortion coefficients
#new_wcs.wcs.set_pv([(2, 1, 1.0)])

# These coefficients do not work well, ra and dec only defined for a 
# small area around the zenith..
#new_wcs.wcs.set_pv([(2, 1, 0.05648),
#                    (2, 2, 8.227e-06),
#                    (2, 3, -1.089e-08)
#                ])

# Divide above coefficients with plate scale as PV terms are meant
# to be dimensionless. Using 3.3/60 = 0.055 degrees as plate scale.
# 3.3 pixels from the camera spec not from measurements
new_wcs.wcs.set_pv([(2, 1, 1.027),
                    (2, 2, 1.496e-04),
                    (2, 3, -1.980e-07)
                ])


# Set the LONPOLE and LATPOLE keywords in ZPN transform
new_wcs.wcs.lonpole = 180.0
new_wcs.wcs.latpole = lat_degrees

# That's it, we have created a new ZPN WCS.
print("New WCS is \n",new_wcs)

outfile = sys.argv[2]
print("\nWriting new file", outfile)

# Create a new header object from WCS and add the original header
# information from the input image.
new_header = new_wcs.to_header() + hdr

# Add Modified Julian Date to the header before writing it
new_header.set('MJD-OBS',mjd,'MJD of Exposure Start' )

# Create new image
new_hdu = fits.PrimaryHDU(data_array,header=new_header)

# Write header to new image
new_hdu.writeto(outfile,overwrite=True)

# END - finished.


