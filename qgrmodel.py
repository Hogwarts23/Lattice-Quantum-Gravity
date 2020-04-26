import numpy as np
from shell import shelling
import math
from scipy.sparse import csc_matrix, linalg as sla
from scipy.sparse.linalg import inv
import scipy.linalg as lin
import numpy.linalg as lin1
###########################################
import matplotlib.pyplot as plt
import math
import time
class QGrModel():
	def __init__(self, tfname,bmass2):
		self.newton = 6.67430*10**-11
		self.baremass2 = bmass2 #the bare mass square
		self.fname = tfname
		self.nlat = None
		self.label = None
		self.contentN4inf = None
		self.matrix = None
		self.lu = None
		self.prop = None
		self.shellinfo = None
		self.corr = None
		self.corr2 = None
		self.maxlen = 50


	def readpartsfromfile(self, keyword):
		file = open(self.fname,'r')
		content = file.read()
		post = content.find(keyword)
		past = content.find("\n\n",post+len(keyword))
		content = content[post+len(keyword):past]
		return content

	def findcol(self, row, collabel, length):
#		if collabel <= length-1 and collabel == self.label[collabel]:
#			return collabel
		l = row
		r = min(collabel,length-1)
		#Trying to bias the binary search
		ratio = (length-1)/self.label[length-1]
		while l <= r:
			mid = l + math.floor((r - l)/2)
			if self.label[mid] == collabel:
				return mid
			elif self.label[mid] < collabel:
				l = mid + 1
			else:
				r = mid - 1
#		for i in range(min(collabel,length-1),row,-1):
#			if label[i] == collabel:
#				return i

	def cstructpropa(self):
		content = self.readpartsfromfile('N4inf: ').split()
		ind1 = 1
		ind2 = 0
		length = int(content[0])
		self.nlat = length
		self.label = np.zeros(length)
		matrix = np.diag((5+self.baremass2)*np.ones(length))
		while ind2 < length:
			self.label[ind2] = int(content[ind1])
			ind1 = ind1 + 6
			ind2 = ind2 + 1
		row = 0
		ind = 1
		while row < length:
			for i in range(1,6):# range(1,6) gives you 1 2 3 4 5
				collabel = int(content[ind+i])
				if self.label[row] < collabel:
					col = self.findcol(row, collabel,length)
					matrix[row,col] = matrix[row,col] - 1
					matrix[col,row] = matrix[col,row] - 1
			row = row + 1
			ind = ind + 6
		self.contentN4inf = content
		self.matrix = matrix
#		np.savetxt('test1',matrix[0,:])
#		np.savetxt('test2',matrix[len(matrix[0,:])-1,:])

	def cstructshelling(self,origin): #origin is the row number
		if self.label is None:
			self.cstructpropa()
		length = self.nlat
		#set up the first element
		tmpdistance = np.zeros(length)+self.label[length-1]+1
		shell = [] 
		tmp = shelling(origin)
		tmp.addlat(origin)
		shell.append(tmp)
		tmpdistance[origin] = 0
		#construct the matrix
		nleft = length-1 # - 1 because we have already add the origin to the shelling
		shellnum = 0
		while nleft != 0:
			shellnum = shellnum + 1
			tmp = shelling(origin,shellnum)
			tmplat = shell[shellnum-1].shelllattice
			for i in tmplat:
				for j in range(1,6):
					col = self.findcol(0,int(self.contentN4inf[int(i*6+1+j)]),length)
					if tmpdistance[col] > shellnum:
						tmp.addlat(col)
						tmpdistance[col] = shellnum
						nleft = nleft-1
			shell.append(tmp)
		self.shellinfo = shell

	def plotdistribution(self):
		lent = len(self.shellinfo)
		x = np.array(range(lent))
		y = np.zeros(lent)
		for i in x:
			y[i] = self.shellinfo[i].numlat
		plt.figure()
		p1 = plt.plot(x,y)
		plt.xlabel(r"Distance")
		plt.ylabel(r"Number of Simplices")
		plt.show()

	def findmax(self,lst):
		maximum = lst[0]
		index = 0
		ind = -1
		for i in lst:
			ind = ind + 1
			if i > maximum:
				maximum = i
				index = ind
		return index, maximum

	def findgeocenter(self): # return the distance
		self.cstructshelling(0)
		lent = len(self.shellinfo)
		y = np.zeros(lent)
		for i in range(lent):
			y[i] = self.shellinfo[i].numlat
		dis,maximum = self.findmax(y)
		return dis

	def correlator(self):
		pro = self.prop
		#print(pro)
		#np.savetxt('test1',pro)
		num = 0
		cor = np.zeros(self.maxlen*4)
		for shell in self.shellinfo:
			s = 0
			for index in shell.shelllattice:
				s = s + pro[int(index)]
			s = s/shell.numlat
			cor[num] = s
			num = num + 1
		return cor[0:self.maxlen]

	def twoparticlecorrelator(self):
		pro2 = self.prop**2
		#print(pro)
		#np.savetxt('test1',pro)
		num = 0
		cor2 = np.zeros(self.maxlen*4)
		for shell in self.shellinfo:
			s = 0
			for index in shell.shelllattice:
				s = s + pro2[int(index)]
			s = s/shell.numlat
			cor2[num] = s
			num = num + 1
		return cor2[0:self.maxlen]

	def correlators(self,numberofcor):
		if self.label is None:
			self.cstructpropa()
		totalcor = np.zeros((self.maxlen,numberofcor))
		dis = self.findgeocenter()
		lst = self.shellinfo[dis].shelllattice
		for i in range(numberofcor):
			#t2 = time.time()
			self.cstructshelling(int(lst[i]))
			#print('number is',int(lst[i]))
			#print(time.time()-t2)
			ori = self.shellinfo[0].origin
			self.propagator(ori)
			cor = self.correlator()
			#print(cor)
			for j in range(self.maxlen):
				totalcor[j,i] = cor[j]
		self.corr = totalcor

	def twoparticlecorrelators(self,numberofcor):
		if self.label is None:
			self.cstructpropa()
		totalcor2 = np.zeros((self.maxlen,numberofcor))
		dis = self.findgeocenter()
		lst = self.shellinfo[dis].shelllattice
		for i in range(numberofcor):
			#t2 = time.time()
			self.cstructshelling(int(lst[i]))
			#print('number is',int(lst[i]))
			#print(time.time()-t2)
			ori = self.shellinfo[0].origin
			self.propagator(ori)
			cor2 = self.twoparticlecorrelator()
			#print(cor)
			for j in range(self.maxlen):
				totalcor2[j,i] = cor2[j]
		self.corr2 = totalcor2

	def propagator(self,source):
		#t1 = time.time()
		if self.lu is None:
			smat = csc_matrix(self.matrix)
			lu = sla.splu(smat)
			self.lu = lu
		#print(time.time()-t1)
		tmp = np.zeros(self.nlat)
		tmp[source] = 1
		x = self.lu.solve(tmp)
		#print(x)
		self.prop = x


	def sinv(self):
#		#t1 = time.time()
		if self.label is None:
			self.cstructpropa()
#		#print(time.time()-t1)
		smat = csc_matrix(self.matrix)
		self.sinverse = inv(smat)
		#matrix = self.sinverse.toarray()
		#np.savetxt('test11',matrix)
		return self.sinverse

# m = QGrModel('./4b0/l4000k16_h0_all_99404',0.000001)

# m.cstructshelling(0)
# m.propagator(0)
# x = m.correlator()

# np.save('testt',x)


# m.twoparticlecorrelators(1)
# t = time.time()
# m.correlators(20)
# print(time.time()-t)
# #x = m.sinv()
#y = x.toarray()
#print(y[:,3042])
#print(time.time()-t)
#print(m.corr)
#m.cstructshelling(0)
#s=0
#for i in m.shellinfo:
#	s = s + i.numlat
#	print(min(i.shelllattice),i.numlat)
#print(len(m.shellinfo))
#print('after shelling')
#m.reshelling()
#for i in m.shellinfo:
#	print(min(i.shelllattice),i.numlat)
#print(len(m.shellinfo))
#print(s)
#print the shelling info
#for i in range(len(m.shellinfo)):
#	print(sorted(m.shellinfo[i].shelllattice))

#m.plotdistribution()





