#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 21:15:53 2023

@author: julio
"""

import numpy as np
import pandas as pd
from crout_factorization_generalization import Crout_generalization 

def f(x): # boundary condition function
    return np.sin(np.pi*x)

def g(x):
    return 0*x

def finite_difference_method_hyperbolic(l,T,     # l is the x endpoint, T is the maximum time t
                                        alpha,   # constant in the parabolic partial d.e.
                                        N,m,     # integers defining the grid
                                        f,g):
    h=l/m
    k=T/N
    Lambda=k*alpha/h
    x=np.linspace(0,l,m+1)
    t=np.linspace(0,T,N+1)
    A=np.zeros((m-1,m-1))
    w=np.zeros(m-1)
    #Lambda=alpha**2 *k/(h**2)
    w0=f(x[1:-1])
    w1=(1-Lambda**2)*f(x[1:-1])+(Lambda**2)/2*(f(x[2:])+f(x[0:-2]))+k*g(x[1:-1])
    for i in range(m-1):# defines Matrix A
        A[i,i]=2*(1-Lambda**2)
        if i<m-2:
            A[i,i+1]=Lambda**2
            A[i+1,i]=Lambda**2
    #j=0
    for tj in t[2:]:
        
        w_j=A@w1-w0
        print(tj)
        print(w_j)
        w0=w1
        w1=w_j
    return A, w_j, w0, w,x

A,w_j,w0,w,x=finite_difference_method_hyperbolic(1,1,2,20,10,f,g)