import os
import shutil
import glob

def move_fits_files(source_dir, destination_dir):
    # Create the destination directory if it doesn't exist
    if not os.path.exists(destination_dir):
        os.makedirs(destination_dir)

    # Get a list of all subdirectories within the source directory
    subdirectories = [subdir for subdir in os.listdir(source_dir) if os.path.isdir(os.path.join(source_dir, subdir))]

    # Iterate over each subdirectory
    for subdir in subdirectories:
        subdir_path = os.path.join(source_dir, subdir)

        # Get a list of all FITS files in the subdirectory
        fits_files = glob.glob(os.path.join(subdir_path, '*.fits'))

        # Move each FITS file to the destination directory
        for fits_file in fits_files:
            shutil.move(fits_file, destination_dir)

# Specify the source directory containing the subdirectories with FITS files
source_directory = '/home/aaron/Desktop/ZASPE/Spectra/HD10700/data/reduced'

# Specify the destination directory where the FITS files will be moved
destination_directory = '/home/aaron/Desktop/ZASPE/Spectra/HD10700'

# Move the FITS files from the source directory to the destination directory
move_fits_files(source_directory, destination_directory)

