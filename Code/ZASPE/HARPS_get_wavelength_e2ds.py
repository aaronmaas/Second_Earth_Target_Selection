import numpy as np
from numpy.polynomial import Polynomial
from astropy.io import fits
import os

# Directory containing the FITS files
fits_directory = '/home/aaron/Desktop/ZASPE/Spectra/GJ277'

# Get a list of FITS files in the directory
fits_files = [file for file in os.listdir(fits_directory) if file.endswith('.fits')]
order_numes = fits.getdata("HARPS...")
# Number of pixels
num_pixels = 4096

# Initialize the array to store wavelengths
wavelengths = np.zeros((len(fits_files), order_number, num_pixels))

# Loop over FITS files
for i, fits_file in enumerate(fits_files):
    fits_path = os.path.join(fits_directory, fits_file)
    
    hdr = fits.getheader(fits_path)
    degree_polynomial = hdr['ESO DRS CAL TH DEG' ]
    coeffecients_orders = hdr['HIERARCH ESO DRS CAL TH COEFF LL*']
    
    order_number = int(np.round(len(coeffecients_orders) / degree_polynomial))
    
    xx = np.arange(num_pixels)
    
    # Loop over orders
    for order in range(order_number):
        coeffecients_order = coeffecients_orders[order * degree_polynomial:(order + 1) * degree_polynomial]
        p = Polynomial(coeffecients_order)
        wavelengths[i, order] = p(xx)

# Check the shape of the wavelengths array
print("Shape of wavelengths array:", wavelengths.shape)

