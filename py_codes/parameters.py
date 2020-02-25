#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : parameters.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 02.09.2020
# Last Modified Date: 02.24.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 22:15:52 2020

@author: hao
"""
import numpy as np
import sys

print(sys.path[0])
# model parameters
J_ex = 1.0
Delta = 1.0 #0.89
spinsize = 0.5
spin_eff = 0.5
delta_S = spinsize - spin_eff

# contractions evaluated by C & Z
# c0 = 1.5747334
# c1 = -0.1042539
# c2 = 0.3444458

# characteristic parameter of new expansion
# C4 = 1.0/8.0
cwd = sys.path[0]
expan_type = 'taylor' #'new_expan' 

if (expan_type == 'taylor'):
    C = -0.25/spinsize
    fname = cwd + '/cintegrals_heis.txt'
    cvalues = np.loadtxt(fname)
elif (expan_type == 'new_expan'):
    C = np.sqrt(1.0 - 0.5/spinsize) - 1.0
    fname = cwd + '/cintegrals_new.txt'
    cvalues = np.loadtxt(fname)
else:
    print('wrong expansion type, exiting...')
    exit()

c0 = cvalues[0]
c1 = cvalues[1]
c2 = cvalues[2]
#-0.5*spinsize*(sqrt(1.0-1.0/(2.0*spinsize))-1.0)

# contractions defined in note of Pontus
# nbar = 0.5*(c0+(Delta-0.5)*c1)-0.5
# mbar = 0.5*(c1+(Delta-0.5)*c2)
nbar = 0.5*( (1.0+delta_S/spin_eff)*c0 + (Delta-0.5)*c1 )-0.5
mbar = 0.5*( (1.0+delta_S/spin_eff)*c1 + (Delta-0.5)*c2 )
cdeltabar = 0.5*(Delta+0.5)*c2
ldeltabar = 0.5*(Delta+0.5)*c1
# Ximinus = -2.0*cdeltabar + 16.0*C4*(Delta+0.5)*nbar \
          # -8.0*C4*(Delta-0.5)*ldeltabar
# Xiplus  = -2.0*mbar + 8.0*C4*(Delta+0.5)*ldeltabar \
          # -16.0*C4*(Delta-0.5)*nbar
# Xizz    = -2.0*nbar + 16.0*C4*(Delta+0.5)*cdeltabar \
          # -16.0*C4*(Delta-0.5)*mbar
# Xiprime = 4.0*C4*((Delta+0.5)*mbar-(Delta-0.5)*cdeltabar)         
Ximinus = -cdeltabar - 4.0*spinsize*C*(Delta+0.5)*nbar \
        + 2.0*spinsize*C*(Delta-0.5)*ldeltabar
Xiplus  = -mbar - 2.0*spinsize*C*(Delta+0.5)*ldeltabar \
        + 4.0*spinsize*C*(Delta-0.5)*nbar
Xizz    = -nbar - 4.0*spinsize*C*(Delta+0.5)*cdeltabar \
        + 4.0*spinsize*C*(Delta-0.5)*mbar
Xiprime = 2.0*spinsize*C*(-(Delta+0.5)*mbar + (Delta-0.5)*cdeltabar)
