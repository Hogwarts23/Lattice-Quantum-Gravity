from scipy import optimize
import astropy.stats as aspys
import numpy as np
import numpy.linalg as lin
from scipy.optimize import curve_fit
import math
from pathlib import Path
import matplotlib.pyplot as plt

maxlen = 50
corperconfig = 20


def bindata(path,maxlen,corperconfig,binsize,bmass):#if ave = 1, return the mean of the fitted parameter, if ave = 0, return original fitted para
	totalcortmp = np.load(Path(path))
	nonzero = len(totalcortmp[:,0])
	numcortmp = len(totalcortmp[0,:])
	if binsize == 1:
		totalcor = totalcortmp
		numcor = numcortmp
	else:
		numcor = math.floor(numcortmp/binsize)
		totalcor = np.zeros((nonzero,numcor))
		for j in range(int(numcor)):
			for i in range(nonzero):
				totalcor[i,j] = np.mean(totalcortmp[i,binsize*j:binsize*j+binsize])
	#print(totalcor[1,:])
	return totalcor

jstderr = np.zeros(40)

bmass = 0.001
for binsize in range(1,41):
	totalcor = bindata(Path('./correlatordata/allcorrelators_m0=0.001000.npy'),maxlen,corperconfig,binsize,bmass)
	jstderr[binsize-1] = np.std(totalcor[2,:])/np.sqrt(len(totalcor[2,:])-1)

plt.figure()
x = np.arange(1,41,1)
plt.plot(x,jstderr)
plt.xlabel('Bin size')
plt.ylabel('Jackknife std error')
plt.title('4b0, m_0 = 0.001')
plt.grid()
plt.savefig('binstderr.png')



