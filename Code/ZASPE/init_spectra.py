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
fits_files = collect_fits_files(directory_name)
    """)

    nb.cells.append(code_cell)
    code02_cell = nbf.v4.new_code_cell(source=f"""
from funcs import plot_arrays, query_tic_name, extract_target_info
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
    
    target_name, observed_night = extract_target_info(fits_files[i])
   
    if target_name.startswith("TIC") == False:
        target_name = query_tic_name(target_name)
      
    
    file_name = target_name + observed_night + ".png"
    
    
    
    save_dir = "/home/aaron/Desktop/Second_Earth_Target_Selection/Code/ZASPE/Plots/" + target_name
    plot_arrays(flux,wavelength,save_dir=save_dir, save_filename=file_name)
    
    sleep(3)
    """)
    nb.cells.append(code02_cell)
    
    
    
    
    code03_cell = nbf.v4.new_code_cell(source=f"""
#example plot
#add the file you intend to show
from IPython.display import display, Image
display(Image(filename='Plots/TIC170729775/TIC1707297752021-03-14T00:09:36.136sp.png'))

    """)
    nb.cells.append(code03_cell)
    
    
    

    # Add a Markdown cell with the header for the first section
    header1_cell = nbf.v4.new_markdown_cell(source="## Combine and Align Spectra")
    nb.cells.append(header1_cell)
    
    # Add a Markdown cell with the description for the first section
    description2_cell = nbf.v4.new_markdown_cell(source="In this step, we will combine and align all available spectra in the same time frame while removing artifacts that we dont want to have in the combined spectra, like saturation or cosmic rays. ")
    nb.cells.append(description2_cell)
    
    
    code_cell21 = nbf.v4.new_code_cell(source=f"""    
from funcs import align_spectra
from funcs import combine_spectra
from astropy.io import fits 
import matplotlib.pyplot as plt


#Specify the spectral order to combine
combined_spectra = []
for i in range(70):
    combined_spectrum = combine_spectra(fits_files, i, blaze_corrected=True)
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
    plt.plot(wavelength[i], combined_spectra[i])
""")
    nb.cells.append(code_cell3) 
    # Add more Markdown cells and code cells for additional sections as needed

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
