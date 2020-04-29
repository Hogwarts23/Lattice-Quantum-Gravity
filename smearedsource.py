from scipy import optimize
import astropy.stats as aspys
import numpy as np
import numpy.linalg as lin
from scipy.optimize import curve_fit
import math
from pathlib import Path
import matplotlib.pyplot as plt
import correlatefitter as cf

def smearefit(func,path,maxlen,corperconfig,binsize,bmass,start,end,numpara,form,ave):#if ave = 1, return the mean of the fitted parameter, if ave = 0, return original fitted para
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


	resamplemean = np.zeros((nonzero,numcor))
	for i in range(nonzero):
		resamplemean[i,:] = jackkniferesamplemean(totalcor[i,:])

	#kv = kvalue(totalcortmp,partresamplemean,start,end)

	xdata = np.arange(start,end)

	logmean = np.log(resamplemean[start:end,:])
	logmcov = np.cov(logmean, bias=True)*(numcor-1)
	# kv = kvalue(np.log(totalcortmp[0:20,:]),logmean,start,end)
	# if binsize == 1:
	# 	for i in range(end-start):
	# 		for j in range(end-start):
	# 			logmcov[i,j] = logmcov[i,j]*kv[i]*kv[j]
	p0l = np.array([0,bmass,0])
	fittedpara1, allchisq = fit(func,logmean,logmcov,numcor,start,end,numpara,p0l)

	totalcor = totalcor**2
	resamplemean = np.zeros((nonzero,numcor))
	for i in range(nonzero):
		resamplemean[i,:] = jackkniferesamplemean(totalcor[i,:])

	#kv = kvalue(totalcortmp,partresamplemean,start,end)

	xdata = np.arange(start,end)

	logmean = np.log(resamplemean[start:end,:])
	logmcov = np.cov(logmean, bias=True)*(numcor-1)
	# kv = kvalue(np.log(totalcortmp[0:20,:]),logmean,start,end)
	# if binsize == 1:
	# 	for i in range(end-start):
	# 		for j in range(end-start):
	# 			logmcov[i,j] = logmcov[i,j]*kv[i]*kv[j]
	p0l = np.array([0,bmass,0])
	fittedpara2, allchisq2 = fit(func,logmean,logmcov,numcor,start,end,numpara,p0l)

	m = np.mean(fittedpara1[1,:])
	M = np.mean(fittedpara2[1,:])
	E = 2*m - M
	return m,M,E