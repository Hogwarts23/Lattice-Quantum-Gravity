# Lattice-Quantum-Gravity
Development starting from January, 2020

## User Guide
Run `getcorrelator.py` and `get2particlecor.py` to read the configurations from directory `4b0` and then store the correlators and two-particle correlators for each configuration in directory `correlatordata` with file name `allcorrelators_m0=0.xxxxxx.npy` and `alltwoparticlecorrelators_m0=0.xxxxxx.npy`.

Run `main3.py` to read the correlator information, fit the parameters and store the mass parameter.

## File Description
-mass: The masses are stored in this directory
- qgrmodel.py : Define class QGrModel(), which read configurations and compute correlators
- shell.py : Define class shell(), which is a class that stores the origin, distance from origin number of simplices and indices of simplices of the shelling at some distance. The total shelling information is a lists of shell class objects.
- correlatefitter.py : Define chi square and the correlated fitter (using `scipy.optimize.curve_fit()`
- getcorrelator.py : Use `qgrmodel.py` and `shell.py` to compute the correlator for each configuration in directory `4b0` and then store the 575 correlators in a matrix in file `./correlatordata/allcorrelators.npy`
- main.py : Read the correlator information from `./correlatordata/allcorrelators.npy` and then fit the parameters using functions from `correlatefitter.py`
-smearedsource.py: This file compares the "smeared source" with the normal method and plot the figure.
-main3.py: Run this file to read correlator information, fit the data and store the mass parameter.
-plotcorrelator.py: Plot the correlators.
-binstd.py: Plot the standard error as a function of bin size.
-comparechisquare.py: Calculate the chi squares for different fitting distance ranges.
