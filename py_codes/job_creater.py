#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 00:24:24 2020

@author: hao
"""

import numpy as np
import os

steps = 50
counter = 0

qxstart = 0
qystart = 0
qxend = 4*np.pi/3.0
qyend = 0
dqx = (qxend-qxstart)/steps
dqy = (qyend-qystart)/steps

for j in range(steps):
    