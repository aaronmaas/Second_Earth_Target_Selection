# __init__.py

from .Functions_spectra import combine_spectra
from .Functions_spectra import align_spectra
from .Functions_spectra import collect_fits_files
from .Functions_spectra import plot_arrays
from .Functions_spectra import extract_target_info_astropy
from .Functions_spectra import extract_target_info
from .Functions_spectra import query_tic_name
from .Functions_spectra import plot_fluxes
from .Functions_spectra import plot_median_spectrum
from .Functions_spectra import interpolate_flux
from .Functions_spectra import store_spectrum_orders
from .Functions_spectra import combine_spectra_interpolation
 

__all__ = ['combine_spectra', 'align_spectra', 'collect_fits_files', 'plot_arrays', 'extract_target_info_astropy', 'extract_target_info', 'query_tic_name','plot_fluxes', 'plot_median_spectrum', 'interpolate_flux','store_spectrum_orders', 'combine_spectra_interpolation']
 
 