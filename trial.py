import numpy as np
import os
from pathlib import Path
import astropy.stats as aspys
import time
import math
import correlatefitter as cf
from scipy.optimize import curve_fit
import threading

def fun_timer():
    print('Hello Timer!')
    timer = threading.Timer(2, fun_timer)
    timer.start()

def func(x,a,m,b):
	return a - m * x - b*np.log(x)

timer = threading.Timer(2, fun_timer)
timer.start()
t0 = time.time()
while time.time()-t0<20:
	print('Hi')
	time.sleep(1)

# bmass = 0.001
# binsize = 20
# totalcortmp = np.load(Path('./correlatordata/allcorrelators_m0=%f.npy'%bmass))
# nonzero = len(totalcortmp[:,0])
# numcortmp = len(totalcortmp[0,:])
# if binsize == 1:
# 	totalcor = totalcortmp
# 	numcor = numcortmp
# else:
# 	numcor = math.floor(numcortmp/binsize)
# 	totalcor = np.zeros((nonzero,numcor))
# 	for j in range(int(numcor)):
# 		for i in range(nonzero):
# 			totalcor[i,j] = np.mean(totalcortmp[i,binsize*j:binsize*j+binsize])

# resamplemean = np.zeros((nonzero,numcor))
# for i in range(nonzero):
# 	resamples = aspys.jackknife_resampling(totalcor[i,:])
# 	for j in range(int(numcor)):
# 		resamplemean[i,j] = np.mean(resamples[j,:])
# start=2
# end=9
# kv = cf.kvalue(totalcortmp,resamplemean[start:end,:],start,end)
# print(kv)


# x = np.arange(10)
# y = cf.jackkniferesamplemean(x)
# estimate, bias, stderr, conf_interval = aspys.jackknife_stats(x,np.mean,0.95)

# print(estimate)
# print(stderr)
# print(np.mean(x))
# print(np.std(x))
# print(np.mean(y))
# print(np.std(y))
# print(np.std(y)*np.sqrt(9))
# print(stderr*np.sqrt(9))

# xdata = np.arange(2,9)
# ydata = np.load('testt.npy')
# ydata = np.log(ydata[2:9])
# popt,pcov = curve_fit(func,xdata,ydata,p0 = [1,0.001,0])
# print(popt)

# for i in np.arange(0.03,0.05,0.001):
# 	print('%6f'%i)

# exit()
# m = np.load(Path('renorm_masses_sqcorr_bin20_4kb0_n11500_m00p001-mf0p05_be0-f2-7_rm0-f2-7.npy'))
# print(np.mean(m[23,:]))