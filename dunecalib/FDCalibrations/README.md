# What is in the directory

This code is standalone and is not compiled with the rest of LArSoft. Instructions on how to use it are in howToUse.txt.

It is code that will perform a selection for stopping muons within the DUNE FD, giving purity and efficiency, using linear cuts and a BDT which must be trained. 

For the files with MBM, it runs the modified box model for a chosen Ccal and provides a chi2 value. The other files use the Recombination Scaling model, and provides three fit parameters for the fit to the ratio equation. You must run it twice to get correct graphical plots - once to get the correct parameters, then again once you have keyed them in to use them. 

The default to use is 10kt, though the files for the 1x2x6 simulation are present in a directory.
