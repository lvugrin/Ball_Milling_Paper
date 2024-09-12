# Ball_Milling_Paper-Experimental part

## Description
This project includes MATLAB scripts, data used for powder X-ray diffraction (PXRD) and Fourier-transformed infrared spectroscopy (FTIR) analysis as part of the research paper titled: **Molecular Dynamics Modelling of Mechanical Activation and Catalysis by Liquid Additives in Ball-Milling Mechanochemistry**. The project is structured into four main directories: **"FTIR data"**, **"PXRD data"**, **"Raman data"** and **"MATLAB scripts"**, each containing files in formats such as `.ASC`, `.m` or `.xy`. The figures produced using the data are created with MATLAB software and saved in `.fig` or `.png` format.

**Project Structure**  
FTIR data - Contains FTIR spectra in .ASC format, which can be analyzed using the provided MATLAB scripts.  

PXRD data - Contains X-ray diffraction data in .xy format for further analysis and plotting.  

Raman data - Contains real-time Raman spectra for further analysis.  

MATLAB scripts - This directory contains MATLAB scripts developed for processing, analyzing, and visualizing FTIR, PXRD, and Raman data.  

**Key MATLAB Scripts**

- **2d Time-Resolved Raman plot.m** was developed for analyzing and processing Raman data collected in a real time. This script allow plotting of two-dimensional Raman plot for detecting changes in the material's vibrational changes during milling.

- **Least Square Analysis method.m** implements the Least Squares Method for detailed analysis of the Raman spectra and enables extracting more accurate information about the Raman spectra and components' behaviour during time.

- **Jar_substracting.m** was used to subtract background caused by the milling vessel. It helps remove contributions from the milling jar in order to better isolate the sample's signal during Raman measurements. 

- **FTIR plot.m** and **PXRD plot.m** was designed for analyzing and plotting raw data, visualizing FTIR spectra and X-ray diffraction patterns.

## Contact Information
Leonarda Vugrin  
lvugrin@irb.hr    
Ruđer Bošković Institute (RBI), Divison of Physical Chemistry, Laboratory for solid-state synthesis and catalysis, Ivan Halasz Group (Milling and Monitoring), Bijenička c. 54, 10000 Zagreb, Croatia  

## Authors and acknowledgment
Leonarda Vugrin  
lvugrin@irb.hr   
Ruđer Bošković Institute (RBI), Divison of Physical Chemistry, Laboratory for solid-state synthesis and catalysis, Ivan Halasz Group (Milling and Monitoring), Bijenička c. 54, 10000 Zagreb, Croatia

## License
This project is licensed under the Creative Commons Attribution 4.0 International License. See the LICENSE file for details.
