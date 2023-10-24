# PAPPy
 Positron Annihilation Program in Python for data analysis in positron annihilation settings. 


The main help file for using (and some editing) the program is in PDF form, title "PAPPy_Help_File.pdf." 

The main program is contained in the PAPPY.py file. Other files are imported into that file. 
analysisFunctions.py is where other users can add custom analysis buttons for extending the range of the application. 
pasdatalib.py contains functions for data manipulation. (Where we should put the curve alignment stuff.)
pasphysics.py contains bound guess, background calculations, and numerical integration functions. 

Wheelhouse is an attempt to include required libraries in the folder so the user does not need them to be installed on their computer. 




## CURRENT PROJECTS

- Make initial guess for peak location based off energy rather than channels. The reason for doing this is that different detectors measure a 511 keV peak in different channels.
-  Using a gaussian fit to find peak center, align all the samples to be centered in the same location. (Add a button to toggle this on or off, because it may not be desired in all settings applications.)
-  Fix the bug with the error function center. When using the MIRION detector which measures 511 keV at ~5500 channels, the error function does not center in the curve and is rendered useless.
-  Look into uncertainty calculation. Some mentioned that it is off by a factor of 2.
-  Look into integration bounds. Some mentioned that it includes one more channel than necessary. (May not be a big deal but compared to PASDA that was what they noticed.)

## current collaborators: 
- Addison Ballif
- Kylee Thomas