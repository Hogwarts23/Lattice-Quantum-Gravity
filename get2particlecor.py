import qgrmodel as qm
import os
import astropy.stats as aspys
import numpy as np
from pathlib import Path
maxlen = 50
corperconfig = 20
bmass = 0.001
#get the configuration directory
dirs = os.listdir(Path("./4b0"))
lent = len(dirs)
#define a matrix to store correlators from all the files
#with rows be the distance and columns be the configuration number
totalcor2 = np.zeros((maxlen,lent*corperconfig))
for file in range(lent):
	m = qm.QGrModel(Path("./4b0")/dirs[file],bmass**2)
	m.twoparticlecorrelators(corperconfig)
	#print(dirs[file],l)
	for i in range(maxlen):
		for k in range(corperconfig):
			totalcor2[i,file*20+k] = m.corr2[i,k]
	#print(col-file) cx

# find the maximum distance and truncate the zero part
#nonzero = 0
#for j in range(lent):
#	for i in range(maxlen):
#		if totalcor[i,j] == 0:
#			if i > nonzero:
#				nonzero = i
#			break
#print('The maximun length is',nonzero)
#totalcor = totalcor[0:nonzero,:]
#save the correlator info to a binary file
np.save(Path('./correlatordata/alltwoparticlecorrelators_m0=%f'%bmass),totalcor2)
