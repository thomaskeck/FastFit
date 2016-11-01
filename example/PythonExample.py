#!/usr/bin/env python3
# Thomas Keck

from PyFastFit import FastFit
import numpy as np

if __name__ == '__main__':

    fitter = FastFit.FastFit(2)
    fitter.setDaughter(0,  1, np.array([ 1.0, 0.0, 0.0]), np.array([-1.0, 0.0, 0.0]), np.diag([0.01] * 7))
    fitter.setDaughter(1, -1, np.array([-1.0, 0.0, 0.0]), np.array([ 1.0, 0.0, 0.0]), np.diag([0.01] * 7))
    fitter.fit(3, 1.5)

    print("Vertex", fitter.getVertex())
    print("Chi2 / NDF", fitter.getChi2(), " / ", fitter.getNDF())
    print("Momentum Daughter 1", fitter.getDaughterMomentum(0))
    print("Momentum Daughter 2", fitter.getDaughterMomentum(1))
    print("Momentum Variance 1", fitter.getDaughterVariance(0))
    print("Momentum Variance 2", fitter.getDaughterVariance(1))