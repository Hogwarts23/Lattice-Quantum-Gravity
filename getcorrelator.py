import qgrmodel as qm
import os
import astropy.stats as aspys
import numpy as np
maxlen = 200
triall = 2
#get the configuration directory
dirs = os.listdir('./4b0')
lent = len(dirs)
#define a matrix to store correlators from all the files
#with rows be the distance and columns be the configuration number
totalcor = np.zeros((maxlen,lent))
for file in range(lent):
	m = qm.QGrModel('./4b0/'+dirs[file])
	m.correlator()
	l = len(m.corr)
	#print(dirs[file],l)
	for i in range(l):
		totalcor[i,file] = m.corr[i]
	#print(col-file)

# find the maximum distance and truncate the zero part
nonzero = 0
for j in range(lent):
	for i in range(maxlen):
		if totalcor[i,j] == 0:
			if i > nonzero:
				nonzero = i
			break
print('The maximun length is',nonzero)
totalcor = totalcor[0:nonzero,:]
#save the correlator info to a binary file
np.save('./correlatordata/allcorrelators',totalcor)

