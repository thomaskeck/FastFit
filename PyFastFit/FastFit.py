#!/usr/bin/env python3

import numpy as np

import ctypes
import ctypes.util
c_double_p = ctypes.POINTER(ctypes.c_double)
c_float_p = ctypes.POINTER(ctypes.c_float)
c_uint_p = ctypes.POINTER(ctypes.c_uint)

import os

# Apparently find_library does not work as I expected
#FastFit_library_path = ctypes.util.find_library('FastFit_CInterface')
#print('Try to load ', FastFit_library_path)
#FastFit_library = ctypes.cdll.LoadLibrary(FastFit_library_path)

FastFit_library =  ctypes.cdll.LoadLibrary(os.getcwd() + '/libFastFit_CInterface.so')
print('Loaded ', FastFit_library)
    

FastFit_library.Create.restype = ctypes.c_void_p
FastFit_library.Create.argtypes = [ctypes.c_uint, ctypes.c_double]

FastFit_library.Delete.argtypes = [ctypes.c_void_p]

FastFit_library.fit.argtypes = [ctypes.c_void_p, ctypes.c_uint]

FastFit_library.SetIPProfile.argtypes = [ctypes.c_void_p, c_double_p, c_double_p]

FastFit_library.SetDaughter.argtypes = [ctypes.c_void_p, ctypes.c_uint, ctypes.c_int, c_double_p, c_double_p, c_double_p]

FastFit_library.getChi2.argtypes = [ctypes.c_void_p]
FastFit_library.getChi2.restype = ctypes.c_double

FastFit_library.GetVertex.argtypes = [ctypes.c_void_p, ctypes.c_uint]
FastFit_library.GetVertex.restype = ctypes.c_double

FastFit_library.GetDaughterMomentum.argtypes = [ctypes.c_void_p, ctypes.c_uint, ctypes.c_uint]
FastFit_library.GetDaughterMomentum.restype = ctypes.c_double

FastFit_library.GetDaughterVariance.argtypes = [ctypes.c_void_p, ctypes.c_uint, ctypes.c_uint, ctypes.c_uint]
FastFit_library.GetDaughterVariance.restype = ctypes.c_double

FastFit_library.GetVariance.argtypes = [ctypes.c_void_p, ctypes.c_uint, ctypes.c_uint]
FastFit_library.GetVariance.restype = ctypes.c_double

FastFit_library.getNDF.argtypes = [ctypes.c_void_p]
FastFit_library.getNDF.restype = ctypes.c_uint


class FastFit(object):
    def __init__(self, numberOfDaughters, magneticField=1.5):
        self.fitter = FastFit_library.Create(int(numberOfDaughters), float(magneticField));
    
    def setIPProfile(self, ip_vertex, ip_variance):
        ip_vertex = np.require(ip_vertex, dtype=np.float64, requirements=['A', 'W', 'C', 'O'])
        ip_variance = np.require(ip_variance, dtype=np.float64, requirements=['A', 'W', 'C', 'O'])
        FastFit_library.SetIPProfile(self.fitter, ip_vertex.ctypes.data_as(c_double_p), ip_variance.ctypes.data_as(c_double_p))

    def setDaughter(self, daughter, charge, momentum, position, variance):
        momentum = np.require(momentum, dtype=np.float64, requirements=['A', 'W', 'C', 'O'])
        position = np.require(position, dtype=np.float64, requirements=['A', 'W', 'C', 'O'])
        variance = np.require(variance, dtype=np.float64, requirements=['A', 'W', 'C', 'O'])
        FastFit_library.SetDaughter(self.fitter, int(daughter), int(charge),
                                    momentum.ctypes.data_as(c_double_p),
                                    position.ctypes.data_as(c_double_p),
                                    variance.ctypes.data_as(c_double_p))

    def fit(self, numberOfIterations=3):
        FastFit_library.fit(self.fitter, int(numberOfIterations))
    
    def getVertex(self):
        return np.array([FastFit_library.GetVertex(self.fitter, int(i)) for i in range(3)])
    
    def getDaughterMomentum(self, daughter):
        return np.array([FastFit_library.GetDaughterMomentum(self.fitter, int(daughter), int(i)) for i in range(3)])
    
    def getDaughterVariance(self, daughter):
        return np.array([[FastFit_library.GetDaughterVariance(self.fitter, int(daughter), int(i), int(j)) for i in range(6)] for j in range(6)])
    
    def getVariance(self):
        return np.array([[FastFit_library.GetVariance(self.fitter, int(i), int(j)) for i in range(6)] for j in range(6)])

    def getChi2(self):
        return FastFit_library.getChi2(self.fitter)
    
    def getNDF(self):
        return FastFit_library.getNDF(self.fitter)

    def __del__(self):
        FastFit_library.Delete(self.fitter)
