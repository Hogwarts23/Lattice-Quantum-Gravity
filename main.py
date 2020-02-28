import numpy as np
import matplotlib.pyplot as plt
import astropy.stats as aspys
import correlatefitter
import time
from inspect import signature
import numpy.linalg as lin
import random
#import sys
#sys.stdout.flush()
def fittedfunc(x,a,m,b):
	return a - m * x - b*np.log(x)
def fittedfunc2(x,a,m,b):
	return a*np.exp(-m*x)+b

def f1():
	t = time.time()
	time.sleep(.2)
	print(time.time()-t)

totalcor = np.load('./correlatordata/allcorrelators.npy')
nonzero = len(totalcor[:,0])
numfile = len(totalcor[0,:])

#for i in range(maxlen):
#	for j in range(maxwid):
#		if totalcor[i,j] < 0:
#			print('nagative')
#			break

#find the jackknife average of
finalcorrelator = np.zeros(nonzero)
err = np.zeros(nonzero)
for i in range(nonzero):
	estimate, bias, stderr, conf_interval = aspys.jackknife_stats(totalcor[i,:],np.mean,0.95)
	finalcorrelator[i] = estimate
	err[i] = stderr


resamplemean = np.zeros((nonzero,numfile))
for i in range(nonzero):
	resamples = aspys.jackknife_resampling(totalcor[i,:])
	for j in range(numfile):
		resamplemean[i,j] = np.mean(resamples[j,:])

#print(invm)
#mcov = np.cov(resamplemean, bias=True)*(numfile-1)


#for start in range(1,10):
#	for end in range(11,20):
#		if end-start<6:
#			continue

start = 3
end = 16
p0 = np.array([1,1,0]) # First guess of the parameters
#choosing a range of distance to do the fitting
partresamplemean = resamplemean[start:end,:]
mcov = np.cov(partresamplemean, bias=True)*(numfile-1)
sig = signature(fittedfunc2)
numpara = len(sig.parameters) - 1
fittedpara, allchisq = correlatefitter.fit(fittedfunc2,partresamplemean,mcov,numfile,start,end,numpara,p0)
#print(fittedpara)
#print(allchisq)
finalpara = np.zeros(numpara)
parastderr = np.zeros(numpara)
finalchisq = np.mean(allchisq)
for i in range(numpara):
	finalpara[i] = np.mean(fittedpara[i,:])
	parastderr[i] = np.std(fittedpara[i,:])
#		print(finalchisq/(end-start),'start is ',start,'end is ',end)

#exit()
nonzero1 = 33
fig1 = plt.figure()
#Plot the correlator with error bars
x = np.arange(1,nonzero1)
y = fittedfunc2(x,finalpara[0],finalpara[1],finalpara[2])
#y = np.exp(fittedfunc(x,finalpara[0],finalpara[1],finalpara[2]))
y1 = finalcorrelator[1:nonzero1]
#err1 = err[0:nonzero]
plt.errorbar(x,y1,yerr = err[1:nonzero1], ecolor = 'red',linewidth = 1, capsize = 1)
#plt.yscale('log')
plt.plot(x,y,'red')
plt.grid(True)
plt.xlabel('Geodesic Distance')
plt.ylabel('Correlator')
plt.title('start %d, end %d'%(start,end))
plt.show()


exit()
i=60
nonzero1 = 33
fig = plt.figure()
while(i>0):
	i = i-1
	j=random.randint(0,numfile-1)
	xdata = np.arange(1,nonzero1)
	ydata1 = np.log(totalcor[1:nonzero1,j])
	ydata2 = fittedfunc(xdata,fittedpara[0,j],fittedpara[1,j],fittedpara[2,j])
	plt.plot(xdata,ydata1)
	plt.plot(xdata,ydata2)
	plt.savefig('./figtrial/fig %d' % i)
	plt.clf()


exit()
fig1 = plt.figure()
#Plot the correlator with error bars
x = np.arange(1,nonzero)
y = fittedfunc(x,finalpara[0],finalpara[1],finalpara[2])
y1 = np.log(finalcorrelator[1:nonzero])
err1 = err[1:nonzero]
plt.errorbar(x,y1,yerr = err1, ecolor = 'red',linewidth = 1, capsize = 1)
#plt.yscale('log')
#plt.plot(x,y)
plt.grid(True)
plt.xlabel('Geodesic Distance')
plt.ylabel('Correlator in log')
plt.show()

