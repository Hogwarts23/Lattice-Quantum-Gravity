import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
onemass = np.load('onemass_4b0_nobin.npy')
twomass = np.load('twomass_4b0_nobin.npy')
N=len(onemass[1,:])
m0 = np.arange(0.001,0.051,0.001)
m = np.mean(onemass,axis=1)
mstd = np.std(onemass,axis=1)*np.sqrt(N-1)
M = np.mean(twomass,axis=1)
binding = 2*m - M
print(N)

# fig = plt.figure()
# l1=plt.plot(m0,m,'blue',label='1-particle-mass')
# l2=plt.plot(m0,M,'green',label='2-particle-mass')
# l3=plt.plot(m0,binding,'red',label='binding energy')
# plt.grid(ls='--')
# plt.legend()
# plt.xlabel('bare mass m_0')
# plt.ylabel('mass parameters')
# plt.title('Fitting distance 3 to 8(ensemble 4b-0.6)')
# plt.savefig('mass4b-0.6.png')


alljm = np.load('renorm_masses_sqcorr_bin20_4kb0_n11500_m00p001-mf0p05_be0-f2-7_rm0-f2-7.npy')
N1 = len(alljm[1,:])
print(N1)
jm=np.mean(alljm,axis = 1)
jstd=np.std(alljm,axis = 1)*np.sqrt(N1-1)
fig = plt.figure()
plt.errorbar(m0,m,label='My mass',yerr=mstd)
plt.errorbar(m0,jm,label='Judah\'s mass',yerr=jstd)
plt.grid(ls='--')
plt.legend()
plt.xlabel('bare mass m_0')
plt.ylabel('mass parameters')
plt.title('Comparison with Judah\'s results(4b0)')
plt.savefig('compareresults_nobin.png')
