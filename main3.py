import numpy as np
import matplotlib.pyplot as plt
import correlatefitter as cf
import time
from inspect import signature
import numpy.linalg as lin
import random
from pathlib import Path
maxlen = 50
corperconfig = 20
binsize = 20
start=2
end=9

def fittedfunc(x,a,m,b):
	return a - m * x - b*np.log(x)

def fittedfunc2(x,a,m,b):
	return a*np.exp(-m*x)/(x**b)

def bindingenergy(x,a1,m,b1,a2,M,b2):
	return (2*a1-a2-2*m*x+M*x+(2*b1-b2)*np.log(x))/x

sig = signature(fittedfunc2)
numpara = len(sig.parameters) - 1

bmass = 0.001

totalcor,finalpara1,parastderr,finalchisq = cf.fits(fittedfunc,Path('./correlatordata/allcorrelators_m0=%f.npy'%bmass),maxlen,corperconfig,binsize,bmass,start,end,numpara,'log')
m=finalpara1[1] 
totalcor,finalpara2,parastderr,finalchisq = cf.fits(fittedfunc,Path('./correlatordata/alltwoparticlecorrelators_m0=%f.npy'%bmass),maxlen,corperconfig,binsize,bmass,start,end,numpara,'log')
M=finalpara2[1]
print(m)
print(M)
print(bindingenergy(1000000,finalpara1[0],m,finalpara1[2],finalpara2[0],M,finalpara2[2]))
print(bindingenergy(50,finalpara1[0],m,finalpara1[2],finalpara2[0],M,finalpara2[2]))
print(bindingenergy(10000000,finalpara1[0],m,finalpara1[2],finalpara2[0],M,finalpara2[2]))


exit()
nonzero1 = 33
fig1 = plt.figure()
#Plot the correlator with error bars
x = np.arange(1,nonzero1)
#y = fittedfunc2(x,finalpara[0],finalpara[1],finalpara[2])
y = np.exp(fittedfunc(x,finalpara[0],finalpara[1],finalpara[2]))
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
plotlenth = 33
plotdata(totalcortmp1,plotlenth)
plotfit(finalpara,form,plotlenth)