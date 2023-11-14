# PAPPy
 Positron Annihilation Program in Python for data analysis in positron annihilation settings. 


The main help file for using (and some editing) the program is in PDF form, title "PAPPy_Help_File.pdf." 

The main program is contained in the PAPPY.py file. Other files are imported into that file. 
analysisFunctions.py is where other users can add custom analysis buttons for extending the range of the application. 
pasdatalib.py contains functions for data manipulation. (Where we should put the curve alignment stuff.)
pasphysics.py contains bound guess, background calculations, and numerical integration functions. 

Wheelhouse is an attempt to include required libraries in the folder so the user does not need them to be installed on their computer. 


PAPPy1 is the original code. PAPPy2 is an attempt to update the code. The main structure is being changed to do calculation based on energy rather than channels. Because this is going to be a large remodel, PAPPy1 is provided as the current usable version. 
