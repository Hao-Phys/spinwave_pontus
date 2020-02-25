#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : main_HF.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 02.09.2020
# Last Modified Date: 02.24.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>
import numpy as np
import LSW
import parameters as paras
import cob
import pycuba

spinsize = paras.spinsize
spin_eff = paras.spin_eff
J_ex = paras.J_ex

def print_header(name):
  print('-------------------- %s test -------------------' % name)

def print_results(name, results):
  keys = ['nregions', 'neval', 'fail']
  text = ["%s %d" % (k, results[k]) for k in keys if k in results]
  print("%s RESULT:\t" % name.upper() + "\t".join(text))
  for comp in results['results']:
    print("%s RESULT:\t" % name.upper() + \
	"%(integral).8f +- %(error).8f\tp = %(prob).3f\n" % comp)

def integrand_HF(ndim, xx, ncomp, ff, userdata):
    
    q1, q2 = [xx[i] for i in range(ndim.contents.value)]
    q = np.array([q1, q2])
    k = np.zeros(2)
    k[0], k[1] = cob.k12tokxy(q[0], q[1])
    omegak = LSW.eigenvalue(q)/(3.0*J_ex*spin_eff)
    gammak = LSW.gammak(k) 
    
    ff[0] = 1.0/omegak
    ff[1] = (gammak)/omegak    
    ff[2] = (gammak)**2/omegak    

    return 0
        
        
res = pycuba.Vegas(integrand_HF, 2, epsabs = 1e-4, \
                   verbose=2, ncomp=3, maxeval=1000000)
    
print_results('Vegas', res)

rres = res.get('results')
len_rres = len(rres)
iintegrals = np.zeros(len_rres)
for flag in range(len_rres):
    iintegrals[flag] = rres[flag].get('integral')
    
print(iintegrals)
fname = 'cintegrals.txt'
np.savetxt(fname, iintegrals)
