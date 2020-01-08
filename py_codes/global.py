#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 22:24:14 2020

@author: hao
"""

import numpy as np

def k12tokxy(k1, k2):
    kx = 2.0*np.pi * k1
    ky = 4.0*np.pi / np.sqrt(3.0) * (k2 - 0.5*k1)
    return kx, ky

def kxytok12(kx, ky):
    k1 = kx/(2.0*np.pi)
    k2 = np.sqrt(3.0)/(4.0*np.pi) * ky + kx/(4.0*np.pi)
    return k1, k2