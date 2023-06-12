
import numpy as np
from numpy.polynomial import Polynomial
from astropy.io import fits

num_pixels = 4096

hdr = fits.getheader('HARPS.2004-09-13T00:16:57.806_e2ds_A.fits')
degree_polynomial = hdr('HIERACH ESO DRS CAL TH DEG')
coeffecients_orders = hdr('HIERACH ESO DRS CAL TH COEFF LL*') 
xx = np.arange(num_pixels)

order_number = np.round(len(coeffecients_order)/degree_polynomial)

wavelength = []
for order in range(len(order_number)): 
	coeffients_order = coeffecients_orders[order:order+4] 
	p = Polynomial([coeffecient_order[0],coeffecient_order[1],coeffecient_order[2],coeffecient_order[3]])
	wavelength.appen(np.array(p(xx))
	
#Check sanity with ceres