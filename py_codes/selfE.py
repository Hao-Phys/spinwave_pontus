#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : selfE.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 02.09.2020
# Last Modified Date: 02.09.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 23:00:46 2020

@author: hao
"""

import numpy as np
import LSW
import parameters as paras
import vertexfunctions as vfun
import cob
import pycuba

# the convergence factor for the integrand of cubic vertex
cov_fact = 0.1
spinsize = paras.spinsize
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
        
        
def Sigma_decay(omega, k):
    
    def integrand_cubic_1(ndim, xx, ncomp, ff, userdata):
        
        q1, q2 = [xx[i] for i in range(ndim.contents.value)]
        q = np.array([q1, q2])
        kmq = k - q
        
        omega_q = LSW.eigenvalue(q)
        omega_kmq = LSW.eigenvalue(kmq)
        
        gamma1 = vfun.cubic_gamma1(q, kmq, k)
        integrand = 0.5*(gamma1 * gamma1.conj())/(omega-omega_q-omega_kmq+1j*cov_fact)
        
        ff[0] = np.real(integrand)
        ff[1] = np.imag(integrand)
        
        return 0
    
    print_header('Vegas')
    
    res = pycuba.Vegas(integrand_cubic_1, 2, epsabs = 1e-3, \
                       verbose=2, ncomp=2, maxeval=1000000)
        
    print_results('Vegas', res)
    
    rres = res.get('results')
    len_rres = len(rres)
    iintegrals = np.zeros(len_rres)
    for flag in range(len_rres):
        iintegrals[flag] = rres[flag].get('integral')
        
    return iintegrals
        
def Sigma_source(omega, k):
    
    def integrand_cubic_2(ndim, xx, ncomp, ff, userdata):
        
        q1, q2 = [xx[i] for i in range(ndim.contents.value)]
        q = np.array([q1, q2])
        mkmq = - k - q
        
        omega_q = LSW.eigenvalue(q)
        omega_mkmq = LSW.eigenvalue(mkmq)
        
        gamma2 = vfun.cubic_gamma2(q, mkmq, k)
        integrand = - 0.5*(gamma2 * gamma2.conj())/(omega+omega_q+omega_mkmq)
        
        ff[0] = np.real(integrand)
        ff[1] = np.imag(integrand)
        
        return 0
    
    print_header('Vegas')
    
    res = pycuba.Vegas(integrand_cubic_2, 2, epsabs = 1e-3, \
                       verbose=2, ncomp=2, maxeval=100000)
        
    print_results('Vegas', res)
    
    rres = res.get('results')
    len_rres = len(rres)
    iintegrals = np.zeros(len_rres)
    for flag in range(len_rres):
        iintegrals[flag] = rres[flag].get('integral')
        
    return iintegrals    


# constants passed to hf energy
J_ex = paras.J_ex
Xiplus = paras.Xiplus
Ximinus = paras.Ximinus
Xizz = paras.Xizz
Xiprime = paras.Xiprime
    
def hartree_fock(q):
    
    k = np.zeros((2, 1))
    k[0], k[1] = cob.k12tokxy(q[0], q[1])
    gamma_k = LSW.gammak(k)
    deltaAk = 1.5*J_ex * (Xiplus*gamma_k + Xizz)
    deltaBk = - 3.0*J_ex * (0.5*Ximinus*gamma_k + Xiprime)
    
    uk, vk = LSW.eigenvector(q)
    epsilonk4 = (uk**2 + vk**2) * deltaAk - 2.0*uk*vk*deltaBk
    
    return epsilonk4

        
