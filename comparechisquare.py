import numpy as np
import matplotlib.pyplot as plt
import correlatefitter as cf
from pathlib import Path
maxlen = 50
corperconfig = 20
binsize = 20
start=2
end=9
numpara = 3
m=0.001
def fittedfunc(x,a,m,b):
	return a - m * x - b*np.log(x)

totalcor,finalpara1,parastderr1,finalchisq1 = cf.fits(fittedfunc,Path('./correlatordata/allcorrelators_m0=%f.npy'%m),maxlen,corperconfig,binsize,m,start,end,numpara,'log',1)
numfile = len(totalcor[1,:])
print('for %d tp %d, the average chi square is %6f'%(start,end,finalchisq1))
start = 2
end = 8
totalcor,finalpara1,parastderr1,finalchisq1 = cf.fits(fittedfunc,Path('./correlatordata/allcorrelators_m0=%f.npy'%m),maxlen,corperconfig,binsize,m,start,end,numpara,'log',1)
numfile = len(totalcor[1,:])
print('for %d tp %d, the average chi square is %6f'%(start,end,finalchisq1))
start = 3
end=8
totalcor,finalpara1,parastderr1,finalchisq1 = cf.fits(fittedfunc,Path('./correlatordata/allcorrelators_m0=%f.npy'%m),maxlen,corperconfig,binsize,m,start,end,numpara,'log',1)
numfile = len(totalcor[1,:])
print('for %d tp %d, the average chi square is %6f'%(start,end,finalchisq1))
start = 3
end = 7
totalcor,finalpara1,parastderr1,finalchisq1 = cf.fits(fittedfunc,Path('./correlatordata/allcorrelators_m0=%f.npy'%m),maxlen,corperconfig,binsize,m,start,end,numpara,'log',1)
numfile = len(totalcor[1,:])
print('for %d tp %d, the average chi square is %6f'%(start,end,finalchisq1))
start = 2
end=7
totalcor,finalpara1,parastderr1,finalchisq1 = cf.fits(fittedfunc,Path('./correlatordata/allcorrelators_m0=%f.npy'%m),maxlen,corperconfig,binsize,m,start,end,numpara,'log',1)
numfile = len(totalcor[1,:])
print('for %d tp %d, the average chi square is %6f'%(start,end,finalchisq1))




