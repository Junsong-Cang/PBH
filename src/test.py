
n = 10000

np_file = '/Users/cangtao/Desktop/tmp.npz'
h5_file = '/Users/cangtao/Desktop/tmp.h5'

import numpy as np
import h5py
from PyLab import *

data = np.zeros((n,n))
s = np.shape(data)
print(s)

np.savez(np_file, d = data)

f = h5py.File(h5_file, 'w')
f.create_dataset('d', data = data)
f.close()

# testing data loading speed
t1 = TimeNow()
d1 = np.load(np_file)['d']
Timer(t1)

t2 = TimeNow()
f = h5py.File(h5_file, 'r')
x1 = f['d'][:]
f.close()
Timer(t2)
