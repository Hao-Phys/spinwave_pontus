#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 23:36:46 2020

@author: hao
"""

"""
Run this script as
python3 main_onshell.py input.txt
"""

import numpy as np
import LSW
import selfE
import sys
import time

st = time.time()
inFile = sys.argv[1]
inM = np.loadtxt(inFile)
num = inM[0]

q = np.array([inM[1], inM[2]])
omegaq_lsw = LSW.eigenvalue(q)

outFile = 'selfE.txt'
f = open(outFile, 'w')
"""
selfE.txt 
ith q1 q2 omega_lsw Re(decay) Im(decay) Re(source) Im(source) hf omega_nlsw decay_rate 
"""
f.write('%4d' %num)
f.write('%8.3f' %q[0])
f.write('%8.3f' %q[1])
f.write('%8.3f' %omegaq_lsw)
res1 = selfE.Sigma_decay(omegaq_lsw, q)
res2 = selfE.Sigma_source(omegaq_lsw, q)
res3 = selfE.hartree_fock(q)

f.write('%20.10f' %res1[0])
f.write('%20.10f' %res1[1])
f.write('%20.10f' %res2[0])
f.write('%20.10f' %res2[1])
f.write('%20.10f' %res3)
tmp = omegaq_lsw+res1[0]+res2[0]+res3
f.write('%20.10f' %tmp)
tmp1 = res1[1]+res2[1]
f.write('%20.10f\n' %tmp1)
f.close()
et = time.time()
print('time elapse is', et-st, 's')
