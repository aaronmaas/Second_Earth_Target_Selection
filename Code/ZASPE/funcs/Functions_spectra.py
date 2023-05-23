#import relevant libraries

import matplotlib.pyplot as plt
import os
import numpy as np
from astropy.stats import sigma_clip
from astropy.io import fits
from scipy.signal import correlate, correlate2d
from scipy.interpolate import interp1d
from scipy import interpolate 
import scipy.constants as c

c_km = c.c/1000 #speed of light in km/s 




#Combining and aligning Spectra with RV correction 

#

def interpolate_flux(wave, flux, wave_ref, RV):
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
        wave_corr = wave[ind] / (1 + RV[ind] / c_km)
        for order in range(np.shape(wave[ind])[0]):
            interp_func = interpolate.interp1d(wave_corr[order], flux[ind][order], fill_value='extrapolate')
            #the interpolation function takes the reference wavelength
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
    num_nights = wave_corr.shape[0]
    num_orders = wave_corr.shape[1]
    num_pixels = wave_corr.shape[2]

    plt.figure()
    for i in range(num_nights):
        x = wave_corr[i][order]
        y = flux[i][order]
        plt.plot(x, y)
    plt.show()



def plot_median_spectrum(wave_ref_corr, flux_median, order):
    """
    Plot the median spectrum.

    Args:
        wave_ref_corr (ndarray): Corrected reference wavelength values.
        flux_median (ndarray): Median flux values.
    """
    
    plt.plot(wave_ref_corr[order], flux_median[order], '--')
    plt.show()



#Combining and aligning Spectra with cross correlation


def combine_spectra(nights, order, saturation_threshold=None, sigma=10.0, CERES_spectra=True, blaze_corrected=False):
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
            print("Currently, only CERES spectra are supported.")
            return None

        order_data = data[extension][order]
        wavelength_spectrum = data[wavelength_extension][order]

        # Append the order spectrum to the list
        spectra.append(order_data)
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
    
    reference_spectrum = clipped_spectra[0]
    aligned_spectra = [reference_spectrum]

    for spectrum in clipped_spectra[1:]:
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
    correlation_half_range = 2  # Number of indices on each side of the deepest absorption line for correlation
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

def collect_fits_files(directory):
    
    """Collects all fits files in a directory. 
    
    Parameters:
    directory (str): 
        directory name
    
    Returns:
    fits_files (list): 
        list of fits files. 
    
    """
    
    fits_files = []
    for filename in os.listdir(directory):
        if filename.endswith('.fits'):
            fits_files.append(os.path.join(directory, filename))
    return fits_files

# plot spectra and save them to directory

def plot_arrays(array1, array2, save_dir=None, save_filename=None):
    
    #The num_plots variable calculates the number of plots, which is the first dimension of the input arrays. The num_rows variable calculates the number of rows needed for the grid, which is half the number of plots rounded up to the nearest integer.
    #The subplots function now takes num_rows and 2 as arguments to create a grid with the desired number of rows and 2 columns.
    #In the loop, we calculate the row and column indices of each subplot and access them using axs[row, col]. The rest of the code is similar to the previous implementation.
    
    num_plots = array1.shape[0]
    num_rows = 14
    num_cols = 5
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



from astropy.io import fits

def extract_target_info_astropy(fits_file):
    # Open the FITS file
    with fits.open(fits_file) as hdul:
        # Get the primary header
        header = hdul[0].header
        
        print(header)

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