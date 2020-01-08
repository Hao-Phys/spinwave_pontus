#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 22:43:11 2020

@author: hao
"""

import numpy as np
import LSW
import cob
import parameters as paras

J_ex = paras.J_ex
spinsize = paras.spinsize

def cubic_gamma1(q1, q2, q3):
    
    u1, v1 = LSW.eigenvector(q1)
    u2, v2 = LSW.eigenvector(q2)
    u3, v3 = LSW.eigenvector(q3)
    
    k1 = np.zeros((2, 1))
    k2 = np.zeros((2, 1))
    k3 = np.zeros((2, 1))
    
    k1[0], k1[1] = cob.k12tokxy(q1[0], q1[1])
    k2[0], k2[1] = cob.k12tokxy(q2[0], q2[1])
    k3[0], k3[1] = cob.k12tokxy(q3[0], q3[1])
    
    gammabar1 = (1.0/3.0)*(np.sin(k1[0]) - 2.0*np.sin(0.5*k1[0]) \
                         * np.cos(0.5*np.sqrt(3.0)*k1[1]))
    gammabar2 = (1.0/3.0)*(np.sin(k2[0]) - 2.0*np.sin(0.5*k2[0]) \
                         * np.cos(0.5*np.sqrt(3.0)*k2[1]))
    gammabar3 = (1.0/3.0)*(np.sin(k3[0]) - 2.0*np.sin(0.5*k3[0]) \
                         * np.cos(0.5*np.sqrt(3.0)*k3[1]))
        
    tgamma1 = gammabar1 * (u1+v1) * (u2*u3 + v2*v3) \
            + gammabar2 * (u2+v2) * (u1*u3 + v1*v3) \
            - gammabar3 * (u3+v3) * (u1*v2 + v1*u2)
            
    gamma1 = 3.0*1j*J_ex*np.sqrt(1.5*spinsize)*tgamma1
    
    return gamma1


def cubic_gamma2(q1, q2, q3):
    
    u1, v1 = LSW.eigenvector(q1)
    u2, v2 = LSW.eigenvector(q2)
    u3, v3 = LSW.eigenvector(q3)
    
    k1 = np.zeros((2, 1))
    k2 = np.zeros((2, 1))
    k3 = np.zeros((2, 1))
    
    k1[0], k1[1] = cob.k12tokxy(q1[0], q1[1])
    k2[0], k2[1] = cob.k12tokxy(q2[0], q2[1])
    k3[0], k3[1] = cob.k12tokxy(q3[0], q3[1])
    
    gammabar1 = (1.0/3.0)*(np.sin(k1[0]) - 2.0*np.sin(0.5*k1[0]) \
                         * np.cos(0.5*np.sqrt(3.0)*k1[1]))
    gammabar2 = (1.0/3.0)*(np.sin(k2[0]) - 2.0*np.sin(0.5*k2[0]) \
                         * np.cos(0.5*np.sqrt(3.0)*k2[1]))
    gammabar3 = (1.0/3.0)*(np.sin(k3[0]) - 2.0*np.sin(0.5*k3[0]) \
                         * np.cos(0.5*np.sqrt(3.0)*k3[1]))
        
    tgamma2 = gammabar1 * (u1+v1) * (u2*v3 + v2*u3) \
            + gammabar2 * (u2+v2) * (u1*v3 + v1*u3) \
            + gammabar3 * (u3+v3) * (u1*v2 + v1*u2)
            
    gamma2 = 3.0*1j*J_ex*np.sqrt(1.5*spinsize)*tgamma2
    
    return gamma2        
    
    
    
    
    