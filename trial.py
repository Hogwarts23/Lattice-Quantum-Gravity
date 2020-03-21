import numpy as np
import os
from pathlib import Path
import astropy.stats as aspys
import time
import math
import correlatefitter as cf
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
# estimate, bias, stderr, conf_interval = aspys.jackknife_stats(x,np.mean,0.95)

# print(estimate)
# print(stderr)
# print(np.mean(x))
# print(np.std(x))
# print(stderr*np.sqrt(9))
# for x in range(2,4,1):
# 	print(x)

# print("********************")

# print(10/3)


exit()
m = np.load(Path('renorm_masses_sqcorr_bin20_4kb0_n11500_m00p001-mf0p05_be0-f2-7_rm0-f2-7.npy'))
print(np.mean(m[23,:]))