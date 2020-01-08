#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 22:15:52 2020

@author: hao
"""

# model parameters
J_ex = 1.0
Delta = 1.0 #0.89
spinsize = 0.5

# contractions evaluated by C & Z
c0 = 1.5747334
c1 = -0.1042539
c2 = 0.3444458

# characteristic parameter of new expansion
C4 = 1.0/8.0
#-0.5*spinsize*(sqrt(1.0-1.0/(2.0*spinsize))-1.0)

# contractions defined in note of Pontus
nbar = 0.5*(c0+(Delta-0.5)*c1)-0.5
mbar = 0.5*(c1+(Delta-0.5)*c2)
cdeltabar = 0.5*(Delta+0.5)*c2
ldeltabar = 0.5*(Delta+0.5)*c1
Ximinus = -2.0*cdeltabar + 16.0*C4*(Delta+0.5)*nbar \
          -8.0*C4*(Delta-0.5)*ldeltabar
Xiplus  = -2.0*mbar + 8.0*C4*(Delta+0.5)*ldeltabar \
          -16.0*C4*(Delta-0.5)*nbar
Xizz    = -2.0*nbar + 16.0*C4*(Delta+0.5)*cdeltabar \
          -16.0*C4*(Delta-0.5)*mbar
Xiprime = 4.0*C4*((Delta+0.5)*mbar-(Delta-0.5)*cdeltabar)         