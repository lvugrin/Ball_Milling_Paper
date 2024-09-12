# Ball_Milling_Paper
# Experimental part

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


# Simulation part

## Description
This project has two python scripts which are used in the paper titled : Molecular Dynamics Modelling of Mechanical Activation and Catalysis by Liquid Additives in Ball-Milling Mechanochemistry. There are two folders "Cluster_Analysis" and "Counting_Complexes with one Python script each of the same name in py format.


- The Cluster_Analysis.py script was developed to track the number of roaming entities and agglomerates with time in the simulations. An agglomerate is defined as a collection of entities ( \ce{K} , \ce{Cl} and 18c6) which have more than 5 neighbours within 10 \AA. Roaming entities are the ones which are not part of an agglomerate.


- The Counting_Complexes.py script was written to count the number of complexes formed by a \ce{K+} ion and a 18c6 molecule using the MDAnalysis suite. A complex was identified when a \ce{K+} ion was within 3.4 Å of all oxygen atoms in a single 18c6 residue.

## Contact Information
Leonarda Vugrin  
lvugrin@irb.hr    
Ruđer Bošković Institute (RBI), Divison of Physical Chemistry, Laboratory for solid-state synthesis and catalysis, Ivan Halasz Group (Milling and Monitoring), Bijenička c. 54, 10000 Zagreb, Croatia      

Rupam Gayen (rupam.gayen@fau.de)
Friedrich-Alexander-University Erlangen-Nürnberg (FAU), Faculty of Science, Department of Physics, PULS Group, Interdisciplinary Center for Nanostructured Films (IZNF), Cauerstrasse 3, 91058, Germany

## Authors and acknowledgment
Leonarda Vugrin  
lvugrin@irb.hr   
Ruđer Bošković Institute (RBI), Divison of Physical Chemistry, Laboratory for solid-state synthesis and catalysis, Ivan Halasz Group (Milling and Monitoring), Bijenička c. 54, 10000 Zagreb, Croatia  

Rupam Gayen
Friedrich-Alexander-University Erlangen-Nürnberg (FAU), Faculty of Science, Department of Physics, PULS Group, Interdisciplinary   Center for Nanostructured Films (IZNF), Cauerstrasse 3, 91058, Germany  

## License
This project is licensed under the Creative Commons Attribution 4.0 International License. See the LICENSE file for details.
