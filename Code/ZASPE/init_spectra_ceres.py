import os
import nbformat as nbf
from funcs import align_spectra
from funcs import combine_spectra
from funcs import collect_fits_files
import sys

#def plot_spectra_function1(star_name):
    # Write your code here to plot spectra using function 1

#def plot_spectra_function2(star_name):
    # Write your code here to plot spectra using function 2

# Other necessary functions and imports

def create_notebook(star_name, directory_name):
    nb = nbf.v4.new_notebook()
    nb.cells = []
    
    #
    header_cell = nbf.v4.new_markdown_cell(source=f"## Analysis of {star_name} Spectra")
    nb.cells.append(header_cell)
    
    # Add a Markdown cell with the description for the first section
    description1_cell = nbf.v4.new_markdown_cell(source="In this step, we will plot all available data for the Star. Please check if the provided targetname and spectra directory are correct.")
    nb.cells.append(description1_cell)
    
    code_cell0 = nbf.v4.new_code_cell(source=f"""star_name = '{star_name}'
directory_name = '{directory_name}'
    """)
    nb.cells.append(code_cell0) 
    # Add code cells and Markdown cells for the first section
    code_cell = nbf.v4.new_code_cell(source="""
from funcs import collect_fits_files
import numpy as np

#collect all fits files in the given directory
fits_files = collect_fits_files(directory_name, extension = 'blaze_correction.fits')
    """)

    nb.cells.append(code_cell)
    code02_cell = nbf.v4.new_code_cell(source=f"""
from funcs import plot_arrays, query_tic_name,query_gliese_name,query_HD_name, extract_target_info, extract_target_info_astropy
import matplotlib.pyplot as plt
from astropy.io import fits 
from tqdm.notebook import tqdm
from time import sleep

for i in tqdm(range(len(fits_files))): 
    spectrum = fits.getdata(fits_files[i])
    flux = spectrum[3,:] #CERES pipeline blaze corrected flux
    wavelength = spectrum[0,:]
    string = fits_files[i] 
    cut_string = string.split("/")[-1][:38] + ".png"
    
    target_name, observed_night = extract_target_info_astropy(fits_files[i])
    t1 = target_name 
    if target_name.startswith("GJ") == False:
        target_name = query_gliese_name(target_name)    #try to get gliese id
    if target_name == "Unknown":
        target_name = query_HD_name(t1)   #try to get HD id
    if target_name == "Unknown":
        target_name = query_tic_name(t1) # try to get TIC id
      
      
    file_name = target_name + observed_night + ".png"
    
    
    
    save_dir = "Plots/" + target_name
    plot_arrays(flux,wavelength,save_dir=save_dir, save_filename=file_name)
    
    sleep(3)
    """)
    nb.cells.append(code02_cell)
    
    
    
    
    code03_cell = nbf.v4.new_code_cell(source=f"""
#example plot
#add the file you intend to show
example_file = save_dir + "/" + file_name
from IPython.display import display, Image
display(Image(filename='example_file'))

    """)
    nb.cells.append(code03_cell)
    

    # Add a Markdown cell with the header for the first section
    header1_cell = nbf.v4.new_markdown_cell(source="## Combine and Align Spectra")
    nb.cells.append(header1_cell)
    
    # Add a Markdown cell with the description for the first section
    description2_cell = nbf.v4.new_markdown_cell(source="In this step, we will combine and align all available spectra in the same time frame while removing artifacts that we dont want to have in the combined spectra, like saturation or cosmic rays. We will use two methods to combine the spectra - RV - correction method and Cross-Correlation Method.")
    nb.cells.append(description2_cell)

    header2_cell = nbf.v4.new_markdown_cell(source="### RV correction method")
    nb.cells.append(header2_cell)
    
    description3_cell = nbf.v4.new_markdown_cell(source="The RV correction is performed using the measured RVs by Harps. The wavelengths are then corrected for RV shift due to different time frames and moved to the first obervation night.")
    nb.cells.append(description3_cell)

    code_cell22 = nbf.v4.new_code_cell(source=f"""  
import scipy.constants as c

c_km = c.c/1000 
    
    
%matplotlib notebook
from funcs import plot_fluxes, plot_median_spectrum, interpolate_flux, get_wavelength_flux_RV, save_spectrum_as_text 

#get fit files
filelist = fits_files

#get observed wavelength, flux and eflux; extract RV and calculate corrected wavelength
wave, flux, eflux, RV, wave_corr, combined_orders = get_wavelength_flux_RV(filelist, ceres = False)

#get reference wavelength from first observation night 
wave_ref_corr = wave[0] / (1 + RV[0] / c_km) 

#interpolate flux at corrected wavelength and on the corrected wavelegnths reference
flux_interp = interpolate_flux(wave, flux, wave_ref_corr, RV)

#median flux for specific order
flux_median = np.median(flux_interp, axis=0)

#controll plots for specific order
norder = 35 
plot_fluxes(wave_corr, flux, norder) 
plot_median_spectrum(wave_ref_corr, flux_median+3000, norder)   

#save the spectrum as txt file for ZASPE
save_directory = "/" + star_name + "/" + star_name + "_RV_method.txt" 
save_spectrum_as_text(combined_orders, wave_ref_corr, flux_median, save_directory)

    """)
    nb.cells.append(code_cell22)
    
    header4_cell = nbf.v4.new_markdown_cell(source="### Cross-correlation method")
    nb.cells.append(header4_cell)
    
    description4_cell = nbf.v4.new_markdown_cell(source="The cross-correlation method takes the median flux from the RV method as star templat to correlate with. It then performes combines and alignes the spectra via cross correlation to the reference template. At the moment there are some issus with the cc method why using only RV method is recommended. In principle should the RV method result in the most accurate shift to the assigned reference spectrum.")
    nb.cells.append(description4_cell)
    
    code_cell21 = nbf.v4.new_code_cell(source=f"""    
from funcs import align_spectra
from funcs import combine_spectra
from astropy.io import fits 
import matplotlib.pyplot as plt

num_orders = np.shape(fits.getdata(fits_files[0]))[1]

#template is the RV corrected version
template = flux_median 

#Specify the spectral order to combine
combined_spectra = []
for i in range(num_orders):
    combined_spectrum = combine_spectra(fits_files, i, template, blaze_corrected=True)
    combined_spectra.append(np.array(combined_spectrum))

combined_spectra = np.array(combined_spectra)

#get wavelength for plotting
spectrum = fits.getdata(fits_files[0])
flux = spectrum[1,:]
wavelength = spectrum[0,:]
blaze_corrected_flux = spectrum[3,:]
    """)
    nb.cells.append(code_cell21)
    
    code_cell3 = nbf.v4.new_code_cell(source=f"""
%matplotlib notebook 

import matplotlib.pyplot as plt 
plt.xlabel(r"Wavelength $[A]$")
plt.ylabel(r"Flux")
for i in range(len(wavelength)):
    plt.plot(wavelength[i], combined_spectra[:,2][i])
""")
    nb.cells.append(code_cell3) 
    
    
    code_cell4 = nbf.v4.new_code_cell(source=f"""
    
from funcs import store_spectrum_orders
save_dir = "Spectra/" + star_name + "/"

#Cross correlation method    
file_name = star_name + "_CCM" + ".txt"

store_spectrum_orders(combined_wavelength,combined_spectra, save_dir = save_dir, file_name = file_name )

""")
    nb.cells.append(code_cell4)     
    
    

    # Save the notebook with a suitable name
    notebook_name = f"{star_name}_spectra.ipynb"
    with open(notebook_name, 'w', encoding='utf-8') as f:
        nbf.write(nb, f)
      


if __name__ == "__main__":
    # Check if the required arguments are provided
    if len(sys.argv) < 3:
        print("Usage: python init_spectra.py [star_name] [directory_name]")
        sys.exit(1)

    # Extract the star_name and directory_name from command-line arguments
    star_name = sys.argv[1]
    directory_name = sys.argv[2]

    create_notebook(star_name, directory_name)
