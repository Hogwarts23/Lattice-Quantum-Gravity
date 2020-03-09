import numpy as np
import os

x = np.load('renorm_masses_sqcorr_bin20_4kb0_n11500_m00p001-mf0p05_be0-f2-7_rm0-f2-7.npy')
print(np.mean(x[0,:]))