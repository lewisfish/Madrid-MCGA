#!/usr/bin/python

import numpy as np
# import matplotlib.pylab as plt
# from scipy.optimize import curve_fit

f = open('image.dat', 'r')
im = np.genfromtxt('image.dat')

f.close()

f = open('test.dat', 'w')
for i in range(len(im[100, :])):
    f.write(str(im[100, i]) + '\n')
f.close()
# def gauss(x, *p):
#     A, mu, sigma = p
#     return A * np.exp(-(x - mu)**2 / (2. * sigma**2))

# hist, bin_edges = np.histogram(im[100, :])
# bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2

# p0 = [1500., 100., 50.]

# coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)

# hist_fit = gauss(bin_centres, *coeff)

# plt.plot(bin_centres, hist, label='Test data')
# plt.plot(bin_centres, hist_fit, label='Fitted data')

# plt.plot(im[100, :])
# plt.show()
