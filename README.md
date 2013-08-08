FSI-Foil
========

Soft Wing Sail Simulation Tool (Designed to solve 2D fluid structure simulation problem around a 2D double layered sail profile)

All simulations should be run from the Controller files.  In these it is possible to control most design parameters.  For any more advanced changes, have a look in the respective classes.  The code is extensively commented, if there are any questions, please don't hesitate to contact me (gplus.to/markdekker).

This work formed part of a Masters Thesis at TU Delft in the Netherlands, so if there are any questions on the theoretical basis of the work or any background information, it would be possible to get a copy of the thesis report.

This software also contains code from other projects/functions:

The code was written in Matlab and requires XFOIL to run (included here).  

The XFOIL class, that calls xfoil.exe is a heavily modified version of Rafael Oliveira's "XFOIL - Matlab Interface" (http://www.mathworks.com/matlabcentral/fileexchange/30478-xfoil-matlab-interface).  

The CALFEMNonL class contains the nonlinear beam element formulation and solver routines from the CALFEM 4.3 toolbox (http://sourceforge.net/projects/calfem/).  

The PlotTool class also incorporates the hline and vline functions by Brandon Kuczenski (http://www.mathworks.com/matlabcentral/fileexchange/1039-hline-and-vline).
