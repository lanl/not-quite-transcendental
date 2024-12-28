#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib as mpl

rc('font',size=18)

def test_func(x):
  return 1.0 + 2.0 * np.power(x, 2.0) + 3.0 * np.power(x, 3.0)

def pressure(rho):
  gamma1 = 5.0 / 3.0
  gamma2 = 2.5
  return 1.0 + np.power(rho, gamma1) + np.power(rho, gamma2) + rho

# NQT
def log2o1(x):
  m, e = np.frexp(x)
  return 2.0 * (m - 1.0) + e

def log10o1(x):
  invlog2o10 = 1.0 / np.log2(10)
  return invlog2o10 * log2o1(x)

def pow2o1(x):
  flr = np.floor(x)
  remainder = x - flr
  mantissa = 0.5 * (remainder + 1.0)
  exponent = flr + 1
  exponent = exponent.astype(int)
  return np.ldexp(mantissa, exponent)

def pow10o1(x):
  log2o10 = np.log2(10.0)
  return pow2o1(log2o10 * x)

def log2o2(x):
    m, e = np.frexp(x)
    return e - (4./3.) * (m - 2) * (m - 1)

def log10o2(x):
    invlog2o10 = 1.0 / np.log2(10)
    return invlog2o10 * log2o2(x) 

def pow2o2(x):
    flr = np.floor(x)
    lm = x - flr - 1
    mantissa = 0.5 * (3 - np.sqrt(1 - 3 * lm))
    exponent = (flr + 1).astype(int)
    return np.ldexp(mantissa, exponent)

def pow10o2(x):
    log2o10 = np.log2(10.0)
    return pow2o2(log2o10 * x)

lxmin = 4.0
lxmax = 12.0
xmin = 10**lxmin
xmax = 10**lxmax
nqt_min_o1 = log10o1(xmin)
nqt_max_o1 = log10o1(xmax)

nqt_min_o2 = log10o2(xmin)
nqt_max_o2 = log10o2(xmax)

nfine = 10000
lfine = np.linspace(lxmin, lxmax, nfine)
fine = 10**lfine
truth = pressure(fine)

#resolutions = [32, 64, 128, 256, 512, 1024]
resolutions = [ 256, 512, 1024]

o2norm = 0.5*1e-2*256*256

fig,axarr = plt.subplots(1, 3,figsize=(16,6),sharey=True)
for i,n in enumerate(resolutions):
    linear_grid = np.linspace(xmin, xmax, n)
    log_grid = np.linspace(lxmin, lxmax, n)
    log_vals = pressure(np.power(10.0, log_grid))
    
    nqt_grid_o1 = np.linspace(nqt_min_o1, nqt_max_o1, n)
    nqt_vals_o1 = pressure(pow10o1(nqt_grid_o1))
    
    nqt_grid_o2 = np.linspace(nqt_min_o2, nqt_max_o2, n)
    nqt_vals_o2 = pressure(pow10o2(nqt_grid_o2))

    log_val = np.interp(np.log10(fine), log_grid, log_vals)
    
    nqt_val_o1 = np.interp(log10o1(fine), nqt_grid_o1, nqt_vals_o1)
    nqt_val_o2 = np.interp(log10o2(fine), nqt_grid_o2, nqt_vals_o2)

    errors_log = abs(truth - log_val) / truth
    errors_nqt_o1 = abs(truth - nqt_val_o1) / truth
    errors_nqt_o2 = abs(truth - nqt_val_o2) / truth

    #plt.loglog(fine, errors_lin)
    axarr[0].loglog(fine, errors_log)
    axarr[1].loglog(fine, errors_nqt_o1)
    axarr[2].loglog(fine, errors_nqt_o2,label='N = {}'.format(n))

for n in resolutions:
    for ax in axarr:
        if n == resolutions[-1]:
            ax.plot(fine, np.ones_like(fine)*o2norm/(n*n), color='r', linestyle='--',label=r'$1/N^2$')
        else:
            ax.plot(fine, np.ones_like(fine)*o2norm/(n*n), linestyle='--',color='r')

axarr[0].set_title('true log')
axarr[1].set_title('NQTo1')
axarr[2].set_title('NQTo2')

axarr[0].set_ylabel('relative error')
axarr[0].set_ylim(1e-5,3)
for ax in axarr:
    ax.set_xlabel('fake density')

axarr[2].legend()

plt.savefig('nqt_convergence.png',bbox_inches='tight')
plt.savefig('nqt_convergence.pdf',bbox_inches='tight')
