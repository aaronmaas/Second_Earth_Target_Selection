# Second_Earth_Target_Selection

This repository contains the neccesarry modules and methods to excess all results from following work in Progress:

1. Maas, A. J. (2023).  *Spectral Frontier: Exploring the Rainbow of Possibilities
Target Selection Pipeline in Search of Earth 2.0*
2. Maas, A. J. (2023) *Spectral Frontier: Exploring the Rainbow of Possibilities
Pipeline revisited - Activity set the limits in Search of Earth2.0* 

Atm the target selection pipeline is focused on combining and alining spectra. All methods for these operations can be found in
\code\ZASPE\. 


#### How to navigate this repository

**Folders:**

- Code/ contains all the notebooks and scripts 
- also contains the modules used in both reports, stored in `/ZASPE/funcs/`
- data/ contains nothing atm
- results/ contain the combined spectra (still need to move them), figures and values for the reports and important calculations 
- images/ conatains presentations and usefull images 

## Project 1: Spectral Frontier: Exploring the Rainbow of Possibilities -  Target Selection Pipeline in Search of Earth 2.0

Notebooks for 10 stars the analysis pipeline as been testet on: \Code\ZASPE\Starname_spectra.ipynb 

To star with the analysis go into directory \Code\ZASPE\ and run 

python init_spectra.py starname directory_name 

starname = objects name (standart is GL or HD) and directory_name = directory where the spectra of the star are stored. It creates a notebook that takes about 15 minutes for 
70 observation nights to run. It will fit blaze functions to the echelle spectra and remove them, perform basic RV analyisis (Lomb Scargle) and combine and align spectra with 
two different methods (Cross-Correlation and RV correction). 

After the combining of the spectra  are fet to ZASPE code (Brahm et al 2017) to obtain stellar intrinsic parameters such as effective temperature, surface gravity, 
radial velocity, rotational velcotiy and metalicity. 

The obtained parameters are then forwarded to PARSEC code \Code\PARSEC\ which performs MCMC to obtain age, mass and interstellar extinction. 

## Project 2: Spectral Frontier: Exploring the Rainbow of Possibilities - Pipeline revisited - Activity set the limits in Search of Earth2.0

- will update the Pipeline for RV signals arising due to activity 


## Required Python packages

- numpy
- scipy
- astropy
- matplotlib
- emcee
- corner 
- pyutil in \code\ZASPE\  dowload project from here https://deepnote.com/workspace/first-deepnote-workspace-1d2a-8723-eaeb4134-8be9-4e57-8654-97dacf9258c3/project/pyutil-e1c1b62c-c84a-4c5c-94ef-537399940de7/notebook/pyutil%2FblazeFitNB-16e76c120c514db699b93011c648e105 

To run ZASPE code: python2.0 (external not in this repo)
