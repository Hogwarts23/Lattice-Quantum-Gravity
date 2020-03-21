#correlater fitter
from scipy import optimize
import astropy.stats as aspys
import numpy as np
import numpy.linalg as lin
from scipy.optimize import curve_fit
import math
from pathlib import Path
import matplotlib.pyplot as plt
def chisq(xdata,ydata,invm,func,para):
	D = len(ydata)
	nump = len(para)
	Chisq = 0
	for i in range(D):
		for j in range(D):
				Chisq =  Chisq + (func(xdata[i],para[0],para[1],para[2])-ydata[i])*(func(xdata[j],para[0],para[1],para[2])- ydata[j])* invm[i,j]
	Chisq = Chisq/(D-nump)
	return Chisq

def inv_cov(y):
	tmp = np.cov(y)
	#a,b = lin.eig(tmp)
	#print(a,'\n')
	#print(tmp-cov_matrix(y))
	#tmp = cov_matrix(y)
	tmpinv = lin.inv(tmp)
	return tmpinv

def jackkniferesamplemean(x):
	s = np.sum(x)
	lent = len(x)
	y = np.zeros(lent)
	for i in range(lent):
		y[i] = (s-x[i])/(lent-1)
	return y


def plotdata(totalcor,plotlength):
	finalcorrelator = np.zeros(plotlength)
	numcor = np.len(totalcor[0,:])
	err = np.zeros(plotlength)
	for i in range(plotlength):
		estimate, bias, stderr, conf_interval = aspys.jackknife_stats(totalcor[i,:],np.mean,0.95)
		finalcorrelator[i] = estimate
		err[i] = stderr
	x = np.arange(1,plotlength)
	y1 = finalcorrelator[1:plotlength]
	plt.errorbar(x,y1,yerr = err[1:plotlength], ecolor = 'red',linewidth = 1, capsize = 1)

def plotfit(finalpara,form,plotlength):
	x = np.arange(1,plotlength)
	if form == 'log':
		y = np.exp(fittedfunc(x,finalpara[0],finalpara[1],finalpara[2]))
		plt.plot(x,y,'red')
	elif form =='exp':
		y = fittedfunc2(x,finalpara[0],finalpara[1],finalpara[2])
		plt.plot(x,y,'red')
	elif form == 'log/r':
		y = np.zeros(plotlength)
		for i in range(plotlength):
			y[i] = np.exp(fittedfunc3(x[i],finalpara[0],finalpara[1],finalpara[2])*x[i])
		plt.plot(x,y,'red')

def fit(func,y,mcov,numcor,start,end,numpara,p=np.array([1,1,1])):
	invm = lin.inv(mcov)
	#print(invm)
	#find the number of parameters automatically
	#define the array that contains parameters and \chi^2 for all configurations
	fittedpara = np.zeros((numpara,numcor))
	allchisq = np.zeros(numcor)
	#xdata = np.arange(nonzero)
	xdata = np.arange(start,end)
	for j in range(numcor):
		ydata = y[:,j]
		popt,pcov = curve_fit(func,xdata,ydata,sigma = mcov,p0 = p)
		fittedpara[:,j] = popt
		allchisq[j] = chisq(xdata,ydata,invm,func,popt)
		#fig = plt.figure()
		#xxdata = np.arange(nonzero)
		#yydata = fittedfunc(xxdata,popt[0],popt[1],popt[2])
		#plt.plot(xxdata,yydata)
		#print(popt)
		#for i in range(numpara)
	#print(numpara)
	return fittedpara, allchisq

def kvalue(unbinned,binnedresamplemean,start,end):
	length = end-start
	width = len(unbinned[0,:])
	widthb = len(binnedresamplemean[0,:])
	unbinnedresamplemean = np.zeros((length,width))
	for i in range(length):
		unbinnedresamplemean[i,:] = jackkniferesamplemean(unbinned[i+start,:])
	kv = np.zeros(length)
	std1 = np.zeros(length)
	std2 = np.zeros(length)
	for i in range(length):
		std1[i] = np.std(binnedresamplemean[i,:])*np.sqrt(widthb-1)
		std2[i] = np.std(unbinnedresamplemean[i,:])*np.sqrt(width-1)
		kv[i] = std1[i]/std2[i]
	return kv


def fits(func,path,maxlen,corperconfig,binsize,bmass,start,end,numpara,form):
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


	resamplemean = np.zeros((nonzero,numcor))
	for i in range(nonzero):
		resamplemean[i,:] = jackkniferesamplemean(totalcor[i,:])

	#kv = kvalue(totalcortmp,partresamplemean,start,end)

	xdata = np.arange(start,end)
	if form == 'exp':
		partresamplemean = resamplemean[start:end,:]
		mcov = np.cov(partresamplemean, bias=True)*(numcor-1)
		kv = kvalue(totalcortmp,partresamplemean,start,end)
		if binsize == 1:
			for i in range(end-start):
				for j in range(end-start):
					mcov[i,j] = mcov[i,j]*kv[i]*kv[j]
		p0 = np.array([1,bmass,0])
		fittedpara, allchisq = fit(func,partresamplemean,mcov,numcor,start,end,numpara,p0)
	elif form == 'log':
		logmean = np.log(resamplemean[start:end,:])
		logmcov = np.cov(logmean, bias=True)*(numcor-1)
		kv = kvalue(np.log(totalcortmp[0:20,:]),logmean,start,end)
		if binsize == 1:
			for i in range(end-start):
				for j in range(end-start):
					logmcov[i,j] = logmcov[i,j]*kv[i]*kv[j]
		p0l = np.array([0,bmass,0])
		fittedpara, allchisq = fit(func,logmean,logmcov,numcor,start,end,numpara,p0l)
	elif form == 'log/r':
		logrmean = np.zeros((end-start,numcor))
		for i in range(end-start):
			for j in range(numcor):
				logrmean[i,j] = np.log(resamplemean[i+start,j])/xdata[i]
		logrcov = np.cov(logrmean, bias=True)*(numcor-1)
		p0lr = np.array([0,bmass,0])
		fittedpara, allchisq = fit(func,logrmean,logrmcov,numcor,start,end,numpara,p0lr)

	finalpara = np.zeros(numpara)
	parastderr = np.zeros(numpara)
	finalchisq = np.mean(allchisq)
	for i in range(numpara):
		finalpara[i] = np.mean(fittedpara[i,:])
		parastderr[i] = np.std(fittedpara[i,:])*np.sqrt(numcor - 1)

	return totalcor,finalpara,parastderr,finalchisq

