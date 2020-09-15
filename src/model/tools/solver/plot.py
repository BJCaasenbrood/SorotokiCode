import os
import numpy as np
import matplotlib.pyplot as plt

cmd = './main config.txt'
os.system(cmd)

q = np.loadtxt("state.log", dtype='f', delimiter=',')
tau = np.loadtxt("tau.log", dtype='f', delimiter=',')

# recover time vector
t = q[:,0]

fig, axs = plt.subplots(2)

for i in range(1,8):
	axs[0].plot(t,q[:,i])

axs[1].plot(t, tau[:,1])	

plt.show()


