import qgrmodel as qm
import os
import astropy.stats as aspys
import numpy as np
from pathlib import Path
maxlen = 50
corperconfig = 20
for mass in np.arange(0.050,0.051,0.001):
	bmass = mass
	#get the configuration directory
	dirs = os.listdir(Path("./4b0"))
	lent = len(dirs)
	#define a matrix to store correlators from all the files
	#with rows be the distance and columns be the configuration number
	totalcor = np.zeros((maxlen,lent*corperconfig))
	for file in range(lent):
		m = qm.QGrModel(Path("./4b0")/dirs[file],bmass**2)
		m.correlators(corperconfig)
		#print(dirs[file],l)
		for i in range(maxlen):
			for k in range(corperconfig):
				totalcor[i,file*20+k] = m.corr[i,k]
		#print(col-file) cx

	np.save(Path('./correlatordata/allcorrelators_m0=%6f'%bmass),totalcor)

