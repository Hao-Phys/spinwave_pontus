#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : LSW.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 02.24.2020
# Last Modified Date: 02.24.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 22:29:44 2020

@author: hao
"""

import numpy as np
import parameters as para
import cob

J_ex = para.J_ex
Delta = para.Delta
spinsize = para.spinsize
spin_eff = para.spin_eff
delta_S = para.delta_S

def gammak(k):
    gamma_k = (1.0/3.0) * (np.cos(k[0]) + 2.0*np.cos(0.5*k[0]) \
                        * np.cos(0.5*np.sqrt(3.0)*k[1]))
    return gamma_k

def eigenvalue(q):
    k = np.zeros((2, 1))
    k[0], k[1] = cob.k12tokxy(q[0], q[1])
    gamma_k = gammak(k)
    
# =============================================================================
#     Ak = 3.0*J_ex*spinsize*(1.0+(Delta-0.5)*gamma_k)
#     Bk = 3.0*J_ex*spinsize*(Delta+0.5)*gamma_k
# =============================================================================
    
    # omegak = np.sqrt( (1.0-gamma_k) * (1.0+2.0*Delta*gamma_k) )
    # epsilonk = 3.0*J_ex*spinsize*omegak
    omegak = np.sqrt((1.0-gamma_k+delta_S/spin_eff) \
            * (1.0+2.0*Delta*gamma_k+delta_S/spin_eff))
    
    epsilonk = 3.0*J_ex*spin_eff*omegak

    return epsilonk
    
def eigenvector(q):
    k = np.zeros((2, 1))
    k[0], k[1] = cob.k12tokxy(q[0], q[1])
    gamma_k = gammak(k)
    
    # Ak = 3.0*J_ex*spinsize*(1.0+(Delta-0.5)*gamma_k)
    # Bk = 3.0*J_ex*spinsize*(Delta+0.5)*gamma_k
    
    Ak = 3.0*J_ex*(1.0+delta_S/spin_eff+(Delta-0.5)*gamma_k)
    Bk = 3.0*J_ex*spin_eff*(Delta+0.5)*gamma_k

    omegak = np.sqrt((1.0-gamma_k+delta_S/spin_eff) \
            * (1.0+2.0*Delta*gamma_k+delta_S/spin_eff))
    epsilonk = 3.0*J_ex*spin_eff*omegak
    # omegak = np.sqrt( (1.0-gamma_k) * (1.0+2.0*Delta*gamma_k) )
    # epsilonk = 3.0*J_ex*spinsize*omegak
    
    uk = np.sqrt( 0.5*(abs(Ak/epsilonk)+1.0) )
    vk = Bk/abs(Bk) * np.sqrt( 0.5*(abs(Ak/epsilonk)-1.0) )
    
    return uk, vk

    
