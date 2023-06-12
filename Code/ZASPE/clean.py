import os
import tarfile
import glob

def extract_tar_files(folder_path):
    tar_files = glob.glob(os.path.join(folder_path, '*.tar'))
    
    for tar_file in tar_files:
        with tarfile.open(tar_file, 'r') as tar:
            tar.extractall(path=folder_path)
        
        os.remove(tar_file)

def delete_non_e2ds_and_CCF_files(folder_path):
    fits_files = glob.glob(os.path.join(folder_path, '*_e2ds_*.fits'))
    fits_files_ccfs = glob.glob(os.path.join(folder_path, '*ccf*.fits')) 
    
    for file in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file)
        if file_path not in fits_files:
        	if file_path not in fits_files_ccfs:
            		os.remove(file_path)

# Specify the folder path where the .tar files are located
folder_path = "/home/aaron/Desktop/ZASPE/Spectra/HD10700"

# Extract .tar files
extract_tar_files(folder_path)

# Delete non-e2ds files
delete_non_e2ds_and_CCF_files(folder_path)

