__author__ = 'mingpu'

import glob
import os
import numpy as np
import math
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt

N1 = 4
N2 = 16
Ne = N1 * N2 * 8 / 16
t1 = 1.0
t2 = 1.2
# h_field = 0.0
k_t1 = 0.0
k_t2 = 0.0
H = np.zeros((N1, N2))
for i in range(N1):
    for j in range(N2):
        kx = 2 * np.pi * (i + k_t1) / (N1 * 1.0)
        ky = 2 * np.pi * (j + k_t2) / (N2 * 1.0)
        # k_t1 = 2 * np.pi * k_t1 / (N1 * 1.0)
        # k_t2 = 2 * np.pi * k_t2 / (N1 * 1.0)
        tem_1 = -4.0 * t2 * (np.cos(kx) * np.cos(ky))
        tem_2 = -2.0 * t1 * (np.cos(kx) + np.cos(ky))
        H[i, j] = tem_1 + tem_2
H = np.reshape(H, N1 * N2, 1)
H.sort()
print H
E_fermi = H[Ne-1]
print E_fermi
a = [E_fermi, E_fermi]
b = [0, 0.5]
# print "%.15f" % (2.0 * sum(H1))
# print math.cos(math.pi)
# print math.sin(math.pi/2.0)
# print np.cos(math.pi)
# print np.mod(-1, 3)
n, bins, patches = plt.hist(H, 400, normed = True, facecolor='green')
plt.plot(a, b, 'r-')
#plt.show()
plt.savefig("xx.pdf", format="pdf", dpi=1000, bbox_inches='tight')

plt.figure(figsize=(10, 10))

a = [-1.0 * np.pi, 1.0 * np.pi]
b = [np.pi,  np.pi]
a = np.array(a)
b = np.array(b)
plt.plot(a, b, 'k-')
plt.plot(a, -1.0 * b, 'k-')
plt.plot(b, a, 'k-')
plt.plot(-1.0 * b, a, 'k-')

a = [-1.0 * np.pi, 0]
b = [0, np.pi]
a = np.array(a)
b = np.array(b)
plt.plot(a, b, 'b-')
plt.plot(a, -1.0 * b, 'b-')
plt.plot(-1.0 * a, b, 'b-')
plt.plot(-1.0 * a, -1.0 * b, 'b-')


k_x_all = []
k_y_all = []
for i in range(-N1/2, N1/2, 10):
    for j in range(-N2/2, N2/2, 10):
        kx = 2 * np.pi * (i + k_t1) / (N1 * 1.0)
        ky = 2 * np.pi * (j + k_t2) / (N2 * 1.0)
        # k_t1 = 2 * np.pi * k_t1 / (N1 * 1.0)
        # k_t2 = 2 * np.pi * k_t2 / (N1 * 1.0)
        tem_1 = -4.0 * t2 * (np.cos(kx) * np.cos(ky))
        tem_2 = -2.0 * t1 * (np.cos(kx) + np.cos(ky))
        if tem_1 + tem_2 < E_fermi:
            # plt.plot(kx, ky, 'ro')
            k_x_all.append(kx)
            k_y_all.append(ky)
        else:
            pass
            # plt.plot(kx, ky, 'bs')
plt.plot(k_x_all, k_y_all, 'ro', markersize=1)
#plt.show()
plt.savefig("yy.pdf", format="pdf", dpi=1000, bbox_inches='tight')
