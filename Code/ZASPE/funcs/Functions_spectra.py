#import relevant libraries

import matplotlib.pyplot as plt
import os
import re
import numpy as np
from astropy.stats import sigma_clip
from astropy.io import fits
from astropy.timeseries import LombScargle
from scipy.signal import correlate, correlate2d
from scipy.interpolate import interp1d
from scipy import interpolate 
import scipy.constants as c
from numpy.polynomial import Polynomial
from tqdm.notebook import tqdm
from time import sleep 


#/home/aaron/Desktop/Second_Earth_Target_Selection/Code/ZASPE
from pyutil.pyutil.blazeFit import blazeFit


c_km = c.c/1000 #speed of light in km/s 


#computing wavelengths from ESO headers 

def calculate_wavelengths(fits_files, fits_directory):
    
    """
    Calculates the wavelengths corresponding to the pixels and orders in a list of FITS files.

    Args:
        fits_files (list): A list of FITS file names.
        fits_directory (str): The directory path where the FITS files are located.

    Returns:
    wavelengths
        numpy.ndarray: An array of shape (num_nights, num_order, num_pixels) that stores the wavelengths.

    This function reads the FITS files, extracts the necessary information from the FITS headers, and performs
    polynomial fitting to obtain the wavelengths for each pixel and order. The wavelengths are then stored in an array
    and returned.
    """
    

    # Initialize the array to store wavelengths
    d = fits.getdata(fits_files[0])

    num_pixels = np.shape(d)[1]
    num_order = np.shape(d)[0]
    num_nights = len(fits_files)
    wavelengths = np.zeros((num_nights, num_order, num_pixels))

    # Loop over nights
    for i, fits_file in enumerate(fits_files):

        fits_path = os.path.join(fits_directory, fits_file)
        hdr = fits.getheader(fits_path)
        degree_polynomial = hdr['ESO DRS CAL TH DEG LL'] + 1  # plus one because there is zeroth order
        coeffecients_orders = sorted(hdr['ESO DRS CAL TH COEFF LL*'], key=natural_sort_key)
        order_number = int(np.round(len(coeffecients_orders) / degree_polynomial))

        xx = np.arange(num_pixels)

        # Loop over orders
        for order in range(order_number):
            # Get coefficients of the fits headers
            if order == 0:
                idx1 = order
                idx2 = idx1 + degree_polynomial  # might have varying polynomial degrees per night
            else:
                idx1 = idx2
                idx2 = idx2 + degree_polynomial
            if order == 0:
                coeffecients_order = coeffecients_orders[idx1:idx2]
            else:
                coeffecients_order = coeffecients_orders[idx1:idx2]
            coeff = []
            for c in coeffecients_order:
                coeff.append(hdr[c])
            # Obtain wavelength array for night and order
            p = Polynomial(np.asarray(coeff))
            wavelengths[i, order] = p(xx)

    # Check the shape of the wavelengths array
    print("Shape of wavelengths array:", wavelengths.shape)

    return wavelengths

#Appliying Blaze fitting 

def apply_blaze_correction(observed_wavelengths,fits_files, ccf_files, output_directory, progress = True):
    
    """
    Applies the blaze correction to the input FITS files and generates blaze-corrected output files.

    Args:
        observed_wavelengths (numpy.ndarray): Array of observed wavelengths.
        fits_files (list): List of FITS file names.
        ccf_files (list): List of CCF file names.
        output_directory (str): Directory path to save the output files.
        progress (bool, optional): Flag to indicate whether to display progress. Defaults to True.

    Returns:
        numpy.ndarray: Array of RV values.
        numpy.ndarray: Array of BERV values.
        numpy.ndarray: Array of BJD values.

    This function reads the FITS files and their corresponding CCF files, applies the blaze correction, and generates
    blaze-corrected output FITS files. It also returns the arrays of RV (Radial Velocity), BERV (Barycentric Earth
    Radial Velocity), and BJD (Barycentric Julian Date) values extracted from the FITS headers.
    
    #you should seperate the RV,BJD and BERV collection from this function and build another one
    #also plot the blaze fit and make it available
    
    """
    
    
    
    c_km = c.c / 1000
    wav = observed_wavelengths
    # Initialize the array to store blaze functions and blaze corrected flux
    d = fits.getdata(fits_files[0])

    num_orders, num_pixels = np.shape(d)
    num_nights = len(fits_files)
    
    discarded_nights = []

    blaze_functions, flux_blaze_corr = np.zeros((num_nights, num_orders, num_pixels)), np.zeros((num_nights, num_orders, num_pixels))

    for night in tqdm(range(len(fits_files))):
        hdu = fits.PrimaryHDU()
        old_header = fits.getheader(fits_files[night])
        ccf_header = fits.getheader(ccf_files[night])
        
        #check if extension is available: 
        if 'ESO DRS CCF RVC' in ccf_header and 'ESO DRS DVRMS' in ccf_header and 'ESO DRS BJD' in old_header:
        
            BERV = ccf_header['ESO DRS BERV']
            RV = ccf_header['ESO DRS CCF RVC']
            RV_error = ccf_header['ESO DRS DVRMS']
            BJD = old_header["ESO DRS BJD"]

            data = fits.getdata(fits_files[night])  # get flux
            v_obs = RV + BERV
            barycentric_wavelength = observed_wavelengths[night] * (1 + v_obs / c_km)  # apply RV + BERV correction

            for order in range(num_orders):
                spec = data[order]
                norm_spec = spec / np.median(spec)
                maxrms = 0.01
                z = blazeFit(wav[night, order], norm_spec, maxrms, numcalls=5)

                wave = wav[night, order]
                wavspread = max(wave) - min(wave)
                wavcent = wave - min(wave) - wavspread / 2.

                cfit = np.poly1d(z)
                blaze_functions[night, order] = cfit(wavcent)
                flux_blaze_corr[night, order] = spec / cfit(wavcent)

            spec_name = old_header["INSTRUME"] + old_header["DATE-OBS"]
            os.makedirs(output_directory, exist_ok=True)

            hdu.header['SPEC_NAME'] = spec_name
            hdu.header['NUM_ORDERS'] = num_orders
            hdu.header['BERV'] = old_header["ESO DRS BERV"]
            hdu.header['BJD'] = old_header["ESO DRS BJD"]
            hdu.header['OBJECT'] = old_header["OBJECT"]
            hdu.header['DATE-OBS'] = old_header["DATE-OBS"]
            hdu.header['RV'] = ccf_header["ESO DRS CCF RVC"]

            barycentric_wavelength = barycentric_wavelength.astype(np.float64)
            data = data.astype(np.float64)
            blaze_functions = blaze_functions.astype(np.float64)
            flux_blaze_corr = flux_blaze_corr.astype(np.float64)

            hdu.data = [barycentric_wavelength, data, blaze_functions[night], flux_blaze_corr[night]]

            output_file = os.path.join(output_directory, f'{spec_name}_blaze_correction.fits')
            hdu.writeto(output_file, overwrite=True)
        else:
            print(f"Skipping {fits_files[night]}: Required header(s) not found.")
            discarded_nights.append(night)
        
        if progress == True:
            sleep(3)
        
    return np.asarray(discarded_nights)


def get_RV_BJD(ccf_files):
     
    """ Function to get out of CCF files in HARPS acillary data the RV, RV uncertainty and BJD. 
    
    Arguments:
    ccf_files (list):
        list of the Cross Correlation Fits files (1 per night)
    
    Returns:
    RV (ndarray)
        1-d array of Radial velocity corrected for BERV and instrumental shift
    RV_error (ndarray)
        1-d array of Radial Velocity uncertainty 
    BJD (ndarray)
        1-d array of Barycentric Julian Date 
    
    
    """
    
    #obtain header for CCFs and data from observation
    num_nights = len(ccf_files)
    
    #initialize arrays
    RV, BJD, RV_error = np.zeros(num_nights),np.zeros(num_nights),np.zeros(num_nights)
    
    #loop over nights 
    for night in range(num_nights): 
        
        ccf_header = fits.getheader(ccf_files[night]) 
        RV[night] = ccf_header['ESO DRS CCF RVC']
        RV_error[night] = ccf_header['ESO DRS DVRMS']
        BJD[night] = ccf_header["ESO DRS BJD"]
    
    return RV, RV_error, BJD





#Combining and aligning Spectra with RV correction 

def get_wavelength_flux_RV(filelist, sigma=10.0, ceres = False):
    """
    Read wavelength, flux, and RV from FITS files and apply sigma clipping.

    Args:
        filelist (list): List of FITS file paths.
        sigma (float): Sigma value for sigma clipping. If None, no sigma clipping will be applied.

    Returns:
        wave (ndarray): Array of wavelength values.
        flux (ndarray): Array of flux values after sigma clipping.
        eflux (ndarray): Array of error in flux values.
        RV (ndarray): Array of RV values.
        wave_corr (ndarray): Array of corrected wavelength values.
        median (ndarray): Array of median flux values.
    """
    dat, hdr = fits.getdata(filelist[0], header=True)
    wave_ref = dat[0]
    if ceres == True:
        wave_ref_corr = wave_ref / (1 + hdr['RV'] / c_km)
    else: 
        wave_ref_corr = wave_ref / (1 + hdr['RV'] / c_km)
    num_orders = np.shape(dat)[1]

    wave = np.zeros((len(filelist), np.shape(dat)[1], np.shape(dat)[2]), 'd')
    wave_corr = np.zeros((len(filelist), np.shape(dat)[1], np.shape(dat)[2]), 'd')
    #np full rather than zeros otherwise after shift first value is zero!
    flux = np.full((len(filelist), num_orders, np.shape(dat)[2]), np.nan, dtype='d')
    eflux = np.zeros((len(filelist), np.shape(dat)[1], np.shape(dat)[2]), 'd')
    RV = np.zeros(len(filelist), 'd')
    BERV = np.zeros(len(filelist), 'd')
    combined_orders = np.zeros((np.shape(dat)[1], np.shape(dat)[2]), 'd')

    for night, filename in enumerate(filelist):
        dat_aux, h_aux = fits.getdata(filename, header=True)
        wave[night] = dat_aux[0]
        flux[night] = dat_aux[3]#/np.median(dat_aux[3]) # normalized for quickcheck
        if ceres == True:
            eflux[night] = dat_aux[4]
            RV[night] = h_aux['RV']   #in CERES it is RV but use in general 
        else:
            RV[night] = h_aux['RV'] 
            BERV[night] = h_aux['BERV']
        wave_corr[night] = wave[night] / (1 + RV[night] / c_km)
        

        # Apply sigma clipping to flux data
        if sigma is not None:
            clip_mask = sigma_clip(flux[night], sigma_lower=sigma, sigma_upper=sigma/4, masked=True)
            clipped_flux = clip_mask.data
            median_night = np.median(clipped_flux)
            flux[night][clip_mask.mask] = median_night
            
    # Exclude nights with non-contributing flux for each order
    for order in range(num_orders):
        order_median = np.median(flux[:, order])

        for night in range(len(filelist)):
            night_median = np.median(flux[night, order])

            dev_threshold_low = 0.2 * order_median
            dev_threshold_high = 5.0 * order_median

            if night_median < dev_threshold_low or night_median > dev_threshold_high:
                print(f"Night {night} in Order {order} has been discarded due to non-contributing flux.")
                flux[night, order] = order_median
    
    for i in range(len(combined_orders)): 
        combined_orders[i] = np.ones(np.shape(dat)[2]) * i 

    return wave, flux, eflux, RV, BERV, wave_corr, combined_orders



def save_spectrum_as_text(order, wave_ref_corr, flux_median, save_dir, output_filename):
    """
    Save the spectrum as a text file.

    Args:
        order (ndarray): Array of order values.
        wave_ref_corr (ndarray): Corrected reference wavelength values.
        flux_median (ndarray): Median flux values.
        save_dir (str): Output directory name.
        output_filename (str): Output text file name.
    
    """
    os.makedirs(os.path.dirname(save_dir), exist_ok=True)
    np.savetxt(output_filename, np.column_stack((np.ravel(order), np.ravel(wave_ref_corr), np.ravel(flux_median))), delimiter=' ')



def interpolate_flux(wave, flux, wave_ref, RV, BERV, BERVs = True ):
    """
    Interpolate flux to a common wavelength reference.

    Args:
        wave (ndarray): Array of wavelength values.
        flux (ndarray): Array of flux values.
        wave_ref (ndarray): Reference wavelength values.
        RV (ndarray): Array of RV values.

    Returns:
        flux_interp (ndarray): Interpolated flux values.
    """
    flux_interp = np.zeros_like(flux)

    for ind in range(len(wave)):
        #if BERVs == False:
        wave_corr = wave[ind] / (1 + (RV[ind]) / c_km)  
        #else: 
        #    wave_corr = wave[ind] / (1 + (BERV[ind] + RV[ind]) / c_km) #wave is barycentric, put it into restframe
        for order in range(np.shape(wave[ind])[0]):
            #interpolated for the corrected wavelengths
            interp_func = interpolate.interp1d(wave_corr[order], flux[ind][order], fill_value='extrapolate')
            #the interpolation function takes the reference wavelength from night 0
            flux_interp[ind][order] = interp_func(wave_ref[order])

    return flux_interp


def plot_fluxes(wave_corr, flux, order):
    """
    Plot the fluxes for each order.

    Args:
        wave_corr (ndarray): Array of corrected wavelength values.
        flux (ndarray): Array of flux values.
        order to plot
    """
    num_nights = flux.shape[0]
    num_orders = flux.shape[1]
    num_pixels = flux.shape[2]

    plt.figure()
    for i in range(num_nights):
        x = wave_corr[i][order]
        y = flux[i][order]
        plt.plot(x, y) 
        plt.xlabel("Wavelength in [A]") 
        plt.ylabel("Flux")
    plt.show()



def plot_median_spectrum(wave_ref_corr, flux_median, order):
    """
    Plot the median spectrum.

    Args:
        wave_ref_corr (ndarray): Corrected reference wavelength values.
        flux_median (ndarray): Median flux values.
    """
    
    plt.plot(wave_ref_corr[order], flux_median[order], '--',label = "Median")
    plt.legend()
    plt.show()



#Combining and aligning Spectra with cross correlation


def combine_spectra(nights, order, template, saturation_threshold=None, sigma=10.0, CERES_spectra=True, blaze_corrected=False):
    '''
    Combine spectra from multiple observation nights for a specific order with crosscorrelation. 
        
    Params:
        nights (list): 
                List of FITS file paths corresponding to the observation nights.
        order (int): 
                Order number of the spectra to extract.
        CERES_spectra (bool, optional): 
                Whether the spectra are CERES spectra (default=True).
        blaze_corrected (bool, optional): 
                Whether the spectra are blaze-corrected in CERES (default=False).

    Returns:
            combined_spectrum (ndarray): 
                Combined spectrum for the specified order.

    This function takes a list of FITS file paths corresponding to different observation nights 
    and combines the spectra for a specific order. The function supports CERES spectra and provides 
    the option to specify whether the spectra are blaze-corrected.

    The function loads the FITS files for each night, extracts the specific order and flux based 
    on the provided options, and appends the order spectrum to a list. Then, it aligns the spectra
    using the `align_spectra` function and calculates the median 
    of the aligned spectra to obtain the combined spectrum.

    The combined spectrum for the specified order is returned as a NumPy array.
    '''
    spectra = []
    discarded_nights = []  # Track discarded nights
    

    for night in nights:
        # Load the FITS file for the night
        hdulist = fits.open(night)
        data = hdulist[0].data
        wavelength_extension = 0
        # Extract the specific order and flux
        if CERES_spectra:
            if blaze_corrected:
                extension = 3  # blaze corrected flux in CERES
            else:
                extension = 1  # flux in CERES
        else:
            print("Currently, only CERES spectra are supported, or same data structure.")
            return None

        order_data = data[extension][order]
        wavelength_spectrum = data[wavelength_extension][order]

        # Append the order spectrum to the list
        spectra.append(order_data)
        
        target_name, observed_night = extract_target_info_astropy(night)
        if target_name and observed_night == "Unknown":
            target_name, observed_night = extract_target_info(night)
        # Close the FITS file
        hdulist.close()

    # Convert the list of spectra into a NumPy array
    spectra = np.array(spectra)
    
    # Calculate the median of the spectra for the specific order
    order_median = np.median(np.median(spectra, axis=0))
    #print(np.median(order_median), order)

    # Apply sigma clipping; atm there is a problem here the sigma clipping does only work for one night and not all of them per order!
    if sigma is not None:
        clipped_spectra = []
        for spectrum, night in zip(spectra,nights):
            clip_mask = sigma_clip(spectrum, sigma_lower=sigma, sigma_upper=sigma, masked=True)
            clipped_mask = clip_mask.data
            median = np.median(clipped_mask)
            mask = clip_mask.mask
            #print(median,order)
            masked_spectrum = []
            for i in range(len(mask)): 
                if mask[i] == True: 
                    masked_spectrum.append(median)
                else: 
                    masked_spectrum.append(clipped_mask[i])
            clipped_spectrum = masked_spectrum
            median_ = np.median(masked_spectrum)
            dev_threshold_low = 0.2 * order_median
            dev_threshold_high = 5. * order_median
            
            #print(diff_median, dev_threshold)
            
            if median_ < dev_threshold_low or median_ > dev_threshold_high:
                print(f"Spectrum from night {night} in order {order} has been discarded due to significant deviation. Median value {median_}")
                discarded_nights.append(night)
                continue
            
        
            clipped_spectra.append(clipped_spectrum)
        clipped_spectra = np.array(clipped_spectra)
    else:
        clipped_spectra = spectra
        
    if len(clipped_spectra) == 0:
        print(f"No valid spectra available for order {order}. Skipping alignment.")
        return 0

    # Find suitable ref spectrum with most flux and deepest absorption lines; good Idea but does not make sense you want to shift to one reference frame!! If different frames for each order get chosen, it does not work!!
    #median_ref = 0
    #for i in range(len(clipped_spectra)):
    #    median = np.median(clipped_spectra[i])
    #    if median > median_ref:
    #        median_ref = median 
    #        order_ref = i 
            #print(i)
    
    reference_spectrum = template[order]
    aligned_spectra = [reference_spectrum]

    for spectrum in clipped_spectra:
        aligned_spectrum = align_spectra(reference_spectrum, spectrum)
        aligned_spectra.append(aligned_spectrum)
    
    

    # Convert the list of aligned spectra into a NumPy array
    aligned_spectra = np.array(aligned_spectra)

    # Calculate the median of the aligned spectra along the first axis (night axis)
    combined_spectrum = np.median(aligned_spectra, axis=0)
    
    if target_name.startswith("TIC") == False:
        target_name = query_tic_name(target_name)
        
    target_name = target_name + 'order' + str(order) +'.fits'
    
    combined_hdulist = fits.PrimaryHDU(combined_spectrum)
    combined_hdulist.writeto("Spectra/" + target_name, overwrite=True)
    
    spectra_order = np.ones(len(combined_spectrum))*order
    
    return spectra_order, wavelength_spectrum, combined_spectrum



def align_spectra(reference_spectrum, spectrum2):
    '''
    Align two spectra based on their flux values using cross-correlation.
    If the observation nights correspond to different time frames, using cross-correlation
    to align the spectra can account for any time shifts and ensure that the spectra are properly
    aligned in the time domain.

    Args:
        reference_spectrum (ndarray):
            Reference spectrum to align to.
        spectrum2 (ndarray):
            Spectrum to align.

    Returns:
        spectrum2_shifted (ndarray):
            Aligned spectrum.

    '''

    # Normalize the spectra by their continuum
    reference_spectrum_normalized = reference_spectrum / np.median(reference_spectrum)
    spectrum2_normalized = spectrum2 / np.median(spectrum2)

    # Find the index of the deepest absorption line in the reference spectrum (excluding zero flux values)
    min_flux_index = np.argmin(reference_spectrum_normalized[np.nonzero(reference_spectrum_normalized)])

    # Define the correlation range
    correlation_half_range = 15  # Number of indices on each side of the deepest absorption line for correlation
    correlation_start = max(min_flux_index - correlation_half_range, 0)
    correlation_end = min(min_flux_index + correlation_half_range + 1, len(reference_spectrum))

    # Perform cross-correlation within the defined range
    correlation = correlate(reference_spectrum_normalized[correlation_start:correlation_end],
                            spectrum2_normalized[correlation_start:correlation_end],
                            mode='same')

    # Find the index of the maximum correlation value
    shift_index = np.argmax(correlation)

    # Calculate the shift value
    shift = shift_index - correlation_half_range
    
    #print(shift)

    # Shift the second spectrum
    spectrum2_shifted = np.roll(spectrum2, shift)
    

    return spectrum2_shifted


# Combine all nights and orders to one spectrum with interpolation 

#Bearbeitung gerade # Combining to one spectrum in the crosscorrelation method

from scipy.interpolate import interp1d 
def combine_spectra_interpolation(wavelengths,spectra):
    num_orders, num_wavelengths = spectra.shape
    # Determine the common wavelength grid
    common_wavelength = np.unique(wavelengths.flatten())

    combined_flux = np.zeros_like(common_wavelength)
    combined_orders = []
    order_array = []
    mean_flux = []
    overlap_indices = []
    
    for order in range(num_orders):
    # Get the wavelength and flux for the current order
        order_wavelengths = wavelengths[order]
        order_flux = spectra[order]

        # Determine the valid range of wavelengths for the current order
        min_wavelength = np.min(order_wavelengths)
        max_wavelength = np.max(order_wavelengths)

        # Exclude regions where the common wavelength grid exceeds the range of the order
        mask = (common_wavelength >= min_wavelength) & (common_wavelength <= max_wavelength)
        #print(mask)
        common_wavelength_masked = common_wavelength[mask]

        # Interpolate the flux onto the common wavelength grid
        interpolated_flux = interp1d(order_wavelengths, order_flux,  fill_value="extrapolate")(common_wavelength_masked)
        
        
        #Combine the interpolated flux
        order_next = order + 1
        
     
        if order_next < num_orders:
            order_wavelengths_next = wavelengths[order_next]
            order_flux_next = spectra[order_next] 
            min_wavelength_next = np.min(order_wavelengths_next)
            max_wavelength_next = np.max(order_wavelengths_next)
            #print(max_wavelength_next,min_wavelength)
            mask_overlap = (common_wavelength <= max_wavelength_next) & (common_wavelength >= min_wavelength)
            next_wavelength_masked = common_wavelength[mask_overlap]
            #print(next_wavelength_masked)
            #print(next_wavelength_masked)
            interpolated_flux_1 = interp1d(order_wavelengths, order_flux, fill_value="extrapolate")(next_wavelength_masked)
            interpolated_flux_2 = interp1d(order_wavelengths_next, order_flux_next, fill_value="extrapolate")(next_wavelength_masked)
            mean_flux_between_order = (interpolated_flux_1 + interpolated_flux_2)/2 
            #here i could save mean flux from previous order and overwrite it again
            mean_flux.append(mean_flux_between_order)
            overlap_indices.append(mask_overlap)
        
        # the 
        if order == 0: 
            combined_flux[mask] = interpolated_flux #this I dont want I want only the overlapping regions to be
            combined_flux[mask_overlap] = mean_flux_between_order # extrapolated not the whole thing.
            combined_orders.extend([order] * np.sum(mask)) #that is not true at the moment
          
        else:
            prev_order = order - 1
            mask_overlap_prev = np.asarray(overlap_indices[prev_order]) 
            # if the order is not zero i need to pay attention to the overlaps. I cannot assign the mask
            # because the mask will include overlaping regions from the previous order
            #indices_where_zero = np.where(combined_flux[mask] == 0)
            #combined_flux[indices_where_zero] = interpolated_flux
            combined_flux[mask] = interpolated_flux
            combined_flux[mask_overlap] = mean_flux_between_order 
            combined_flux[mask_overlap_prev] = np.asarray(mean_flux[prev_order]) 
            combined_orders.extend([order] * np.sum(mask)) #thats not true at the moment
            

    combined_orders = np.asarray(combined_orders)

    return np.flip(common_wavelength), np.flip(combined_flux), combined_orders#, interpolated_flux_1, interpolated_flux_2


#store_spectrum in cross correlation method

def store_spectrum_orders(wavelength_spec, spectrum, save_dir=None, file_name=None):
    num_orders = np.shape(wavelength_spec)[0]  # Number of spectral orders

    # Combine all orders into a single array
    data = np.empty((0, 3))
    for order_idx in range(num_orders):
        # Get the wavelength and flux values for the current order
        wavelength = wavelength_spec[order_idx]
        flux = spectrum[order_idx]
        #print(np.shape(wavelength),np.shape(flux),order_idx)
        # Combine order index, wavelength, and flux into a single array
        order_data = np.column_stack((np.full_like(wavelength, order_idx), wavelength, flux))

        # Append order data to the main data array
        data = np.concatenate((data, order_data), axis=0)

    # Define the file path and name for all orders
    file_path = os.path.join(save_dir, file_name)

    # Create the directory if it does not exist
    os.makedirs(os.path.dirname(file_path), exist_ok=True)

    # Save the data to the text file
    np.savetxt(file_path, data, delimiter="\t")

    print(f"The spectrum for all orders has been stored in the file: {file_path}")




    

# Fits files 

def collect_fits_files(directory, extension = None):
    
    """Collects all fits files in a directory. 
    
    Parameters:
    directory (str): 
        directory name
    
    Returns:
    fits_files (list): 
        list of fits files. 
    
    """
    
    fits_files = []
    
    if extension == None:
        extension = ('.fits')
    
    for filename in os.listdir(directory):
        if filename.endswith(extension):
            fits_files.append(os.path.join(directory, filename))
    return sorted(fits_files) 

def separate_fits_files(file_list):
    
    """Separates all fits files in a directory. 
    
    Parameters:
    file_list (list): 
        list of fits files
        
    
    Returns:
    e2ds_files (list): 
        list of sorted e2ds fits files 
    ccf_files (list): 
        list of sorted ccf fits files 
    
    """
    e2ds_files = []
    ccf_files = []

    for file_path in file_list:
        file_name, file_ext = os.path.splitext(file_path)
        if file_ext == ".fits":
            if "e2ds" in file_name:
                e2ds_files.append(file_path)
            elif "ccf" in file_name:
                ccf_files.append(file_path)

    return sorted(e2ds_files), sorted(ccf_files)

def natural_sort_key(extension):
    
    """
    Custom sorting key for ESO DRS HARPS fits files. Is used to sort 
    the wavelength coefficients for the calculation of the wavelengths.
    Coeffcients usually in header["ESO DRS CAL TH COEFF LL*"]. The sorting
    ensures that the appropriate coefficients are used for a specific 
    wavelength order. The key can be used in the python inbuild sorted() 
    function. 
    
    Parameter: 
    extension 
    
    Returns:
        number or extension 
    
    """
    
    # Define a custom sorting key function
    # Extract the numeric part from the extension
    match = re.search(r'\d+', extension)
    if match:
        number = int(match.group())
        return number
    return extension


# plot spectra and save them to directory

def plot_arrays(array1, array2, save_dir=None, save_filename=None):
    
    #The num_plots variable calculates the number of plots, which is the first dimension of the input arrays. The num_rows variable calculates the number of rows needed for the grid, which is half the number of plots rounded up to the nearest integer.
    #The subplots function now takes num_rows and 2 as arguments to create a grid with the desired number of rows and 2 columns.
    #In the loop, we calculate the row and column indices of each subplot and access them using axs[row, col]. The rest of the code is similar to the previous implementation.
    
    num_plots = array1.shape[0]
    num_cols = 5
    num_rows = int(np.ceil(array1.shape[0]/num_cols))
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(15, 30))
    
    for i in range(num_plots):
        row = i // num_cols
        col = i % num_cols
        axs[row, col].plot(array2[i,:], array1[i,:])
        axs[row, col].set_xlabel('Wavelength')
        axs[row, col].set_ylabel('Flux')
        axs[row, col].set_title('Grating Order {}'.format(i+1))
    
    #plt.tight_layout()
    
    if save_dir and save_filename:
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        plt.savefig(os.path.join(save_dir, save_filename))
        
    plt.close()
    #else:
    #    plt.show()
    
    
def plot_rvs(rv,rv_err,bjd,save_dir=None, save_filename=None,ms = True):
    
    """ Function to plot the RV variation against barycentric julian date. 
    
    Arguments: 
    
    rv (ndarry)
        1-d Radial Velocity in km/s array of length nights of observation
    rv_err (ndarray)
        1-d Radial Velocity Uncertainty in km/s array of length nights of observation
    bjd (ndarray)
        1-d Barycentric Julian Date in days array of length nights of observation
    save_dir (string)
        Directory where the figure is save. None by default
    save_filename=None
        Filename of the figure in save_dir. None by default
    ms (bool)
        meter per second conversion by default true
      
        
    Returns
    figure of rv variation
    
    """
    
    plt.figure()
    
    rv_variation = rv - np.mean(rv)
    if ms == True:
        rv_variation = rv_variation * 1000
        rv_err = rv_err * 1000
    reduced_bjd = bjd - np.min(bjd)
    plt.errorbar(reduced_bjd,rv_variation,yerr = rv_err, marker = "s", linestyle = "None") 
    plt.xlabel("Time [days]")
    plt.ylabel("RV variation [km/s]")
    
    if save_dir and save_filename:
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        plt.savefig(os.path.join(save_dir, save_filename))
        
        
def plot_lomb_scargle(bjd,rv,rv_err, false_alarm_probability = 0.01 , save_dir = None, save_filename = None): 
    
    """ Function to plot a Lomb Scargle Peridogram (based on astropy). Power is plottet against time in days. 
    
    Arguments: 
    
    rv (ndarry)
        1-d Radial Velocity in km/s array of length nights of observation
    rv_err (ndarray)
        1-d Radial Velocity Uncertainty in km/s array of length nights of observation
    bjd (ndarray)
        1-d Barycentric Julian Date in days array of length nights of observation
    false_alarm_level (float)
        Probability of detecting a false alarm. By default 1% 
    save_dir (string)
        Directory where the figure is save. None by default
    save_filename=None
        Filename of the figure in save_dir. None by default
 
        
    Returns
    figure of Lomb Scargle Periodogram
    
    """
    
    rv_variation = rv - np.mean(rv)
    reduced_bjd = bjd - np.min(bjd)
    
    ls = LombScargle(reduced_bjd, rv_variation, rv_err) #is that error
    frequency, power = ls.autopower()
    
    periods = 1/frequency
    false_alarm = ls.false_alarm_level(false_alarm_probability) 
    best_period = periods[np.argmax(power)]
    plt.plot(periods, power, label = "Lomb Scargle Periodogram") 
    x_fa = np.linspace(np.min(periods),np.max(periods),100)
    plt.plot(x_fa,np.ones(len(x_fa)) *false_alarm, "--", label = "{:2f}-False alarm probability".format(false_alarm))
    plt.text(np.max(periods)/3, np.max(power), "Period = {:.2f} d".format(best_period)) 
    plt.xlabel("Time [days]")
    plt.ylabel("Power")
    plt.legend(loc = "lower right", frameon = False)
    
    if save_dir and save_filename:
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        plt.savefig(os.path.join(save_dir, save_filename))
    
    
    
from astroquery.simbad import Simbad
def query_tic_name(target_name):
    custom_simbad = Simbad()
    custom_simbad.add_votable_fields("ids")
    
    result_table = custom_simbad.query_object(target_name)
    
    if result_table is not None and len(result_table) > 0:
        ids = result_table['IDS'][0].split("|")
        tic_name = next((id for id in ids if id.startswith("TIC ")), None)
        if tic_name is not None:
            tic_name = tic_name.split("TIC ")[1]
        return 'TIC'+tic_name
    
    return "Unknown"



def query_gliese_name(target_name):
    custom_simbad = Simbad()
    custom_simbad.add_votable_fields("ids")
    
    result_table = custom_simbad.query_object(target_name)
    
    if result_table is not None and len(result_table) > 0:
        ids = result_table['IDS'][0].split("|")
        gliese_name = next((id for id in ids if id.startswith("GJ ")), None)
        if gliese_name is not None:
            gliese_name = gliese_name.split("GJ ")[1]
            return 'GJ'+gliese_name
    
    return "Unknown"

def query_HD_name(target_name):
    
    custom_simbad = Simbad()
    custom_simbad.add_votable_fields("ids")
    
    result_table = custom_simbad.query_object(target_name)
    
    if result_table is not None and len(result_table) > 0:
        ids = result_table['IDS'][0].split("|")
        draper_name = next((id for id in ids if id.startswith("HD ")), None)
        if draper_name is not None:
            draper_name = draper_name.split("HD ")[1]
            return 'HD'+draper_name
    return "Unknown"



from astropy.io import fits

def extract_target_info_astropy(fits_file):
    # Open the FITS file
    with fits.open(fits_file) as hdul:
        # Get the primary header
        header = hdul[0].header
        
        #print(header)

        # Extract the target name and observed night
        target_name = header.get('OBJECT', 'Unknown')
        observed_night = header.get('DATE-OBS', 'Unknown')
        
        if target_name == "Unknwon" or observed_night == "Unknown":
            target_name, observed_night = extract_target_info(fits_file)
            

    return target_name, observed_night

import os

def extract_target_info(fits_file):
    # Extract target name and observed night from file name
    file_name = os.path.splitext(os.path.basename(fits_file))[0]
    target_name, observed_night = file_name.split('_')

    return target_name, observed_night