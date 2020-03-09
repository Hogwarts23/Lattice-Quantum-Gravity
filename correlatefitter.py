#correlater fitter
from scipy import optimize
import numpy as np
import numpy.linalg as lin
from scipy.optimize import curve_fit
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