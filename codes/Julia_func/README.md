# processMOPS
codes for plotting MOPS climatology
1. *plot_jointpdfMOPS*: plotting function for the joint pdf with histograms on x, y axes. 
2. *plot_MOPS_3var_hist*: preps the data for jointpdfMOPS
3. *plot_MOPlineAverage*: plots the average MOP line: this is in its infancy stages...
4. *get_MOPname*: give it a latitude and it will spit out the name of the nearest MOP line
5. *get_EfluxMOP*: spits out the shoreward (and alongshore) energy flux from the MOP 1D spectra and a1, b1 fourier components
6. *read_MOPline*: give it bounding dates and a MOP name, and it spits out all the info necessary for these processing codes
