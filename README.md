# Lattice-Quantum-Gravity
Development starting from January, 2020, with the help of Jack and Judah

## User Guide
Run `getcorrelator.py` to read the configurations from directory `4b0` and then store the correlator for each configuration in directory `correlatordata` with file name `allcorrelators.npy`

Run `main.py` to read the correlator information, fit the parameters and then plot the graph

## File Description
- qgrmodel.py : Define class QGrModel(), which read configurations and compute correlators
- shell.py : Define class shell(), which is a class that stores the origin, distance from origin number of simplices and indices of simplices of the shelling at some distance. The total shelling information is a lists of shell class objects.
- correlatefitter.py : Define chi square and the correlated fitter (using `scipy.optimize.curve_fit()`
- getcorrelator.py : Use `qgrmodel.py` and `shell.py` to compute the correlator for each configuration in directory `4b0` and then store the 575 correlators in a matrix in file `./correlatordata/allcorrelators.npy`
- main.py : Read the correlator information from `./correlatordata/allcorrelators.npy` and then fit the parameters using functions from `correlatefitter.py`
