# Opioids-depress-breathing-through-two-small-brainstem-sites-Bachmutskyetal.2020-
These scripts execute analyses in 'Opioids depress breathing through two small brainstem sites' (Bachmutsky et al. 2020)  

All code is written in MATLAB. In order to use code, download or clone repository and add to MATLAB path. Within-function paths referencing locations of example data will have to be updated in order to run functions appropriately.  

Missing Example Extracellular Data:  
The 'In Vitro Extracellular Analysis Fig. 4' folder contains functions that reference example data files too large (exceeding 100 MB in size) to upload to the Github repository. Example data is therefore available upon request to the corresponding author (Kevin.Yackle@ucsf.edu) or first author (Iris.Bachmutsky@ucsf.edu).  

Citations to Outside Code:  
Functions in the 'From Outside Sources' folder contain various contributions from MathExchange and other places that are used within our custom scripts. License files have been kept when possible. These functions are cited as follows:  
Bernal, Jose Maria Garcia-Valdecasas (2011). colorGradient, MATLAB Central File Exchange. Retrieved Feb, 2019.  
Greene, Chad (2019). hex2rgb, MATLAB Central File Exchange. Retrieved June, 2019.  
Bockstege, John (2012). jbfill, MATLAB Central File EXchange. Retrieved June, 2019.  
Hughley, Jake (2018). nestedSortStruct, MATLAB Central File Exchange. Retrieved June, 2019.  
Lansey, Jonathan C. (2015). nhist, MATLAB Central File Exchange. Retrieved May, 2019.  
Llimona, Quim, P 2018. suptitle. DOI: https://github.com/lemonzi/matlab/tree/master/suptitle.  
Kumpulainen, Pekka (2016). tight_subplot, MATLAB Central File Exchange. Retrieved March, 2019.  

Peak finding algorithm for Extracellular data, called in extractpeaks.m is based on:  
  Du, P., Kibbe, W.A. and Lin, S.M. (2006) Improved peak detection in mass spectrum by incorporating continuous wavelet transform-based pattern matching, Bioinformatics, 22, 2059-2065. DOI: https://www.bioconductor.org/packages/release/bioc/html/MassSpecWavelet.html.
