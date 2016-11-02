#!/usr/bin/env python3
# Thomas Keck

import os
import sys

from PyFastFit import FastFit

import numpy as np
import pandas
import matplotlib.pyplot as plt
import scipy.stats

if __name__ == '__main__':

    files = [(2, 'D0:2.pickle'), (3, 'D+:3.pickle'), (4, 'D0:4.pickle'), (5, 'D+:5.pickle')]

    for n, f in files:
        df = pandas.read_pickle(os.path.join('data', f))
        fitter = FastFit.FastFit(n, 1.5)

        df['fastfit_dx'] = 0.0
        df['fastfit_dy'] = 0.0
        df['fastfit_dz'] = 0.0
        
        df['fastfit_chiProb'] = 0.0

        for i in range(n):
            df['fastfit_px_{}'.format(i)] = 0.0
            df['fastfit_py_{}'.format(i)] = 0.0
            df['fastfit_pz_{}'.format(i)] = 0.0

        # TODO Variance is not yet correctly calculated by FastFit
        for j in range(7):
            for k in range(7):
                df['fastfit_momVertCovM{}{}'.format(j, k)] = 0.0

        df = df[:1000]

        for index, decay in df.iterrows():
            sys.stdout.write("\r{}".format(index))
            for i in range(n):
                charge = int(decay['charge_{}'.format(i)])
                momentum = np.array([decay['px_{}'.format(i)], decay['py_{}'.format(i)], decay['pz_{}'.format(i)]])
                position = np.array([decay['dx_{}'.format(i)], decay['dy_{}'.format(i)], decay['dz_{}'.format(i)]])
                variance = np.zeros((6, 6))

                m = [4, 5, 6, 0, 1, 2]
                for j in range(6):
                    for k in range(6):
                        variance[j, k] = decay['momVertCovM{}{}_{}'.format(m[j], m[k], i)]
                fitter.setDaughter(i, charge, momentum, position, variance)

            fitter.fit(3)

            vertex = fitter.getVertex()
            df.loc[index, 'fastfit_dx'] = vertex[0]
            df.loc[index, 'fastfit_dy'] = vertex[1]
            df.loc[index, 'fastfit_dz'] = vertex[2]
            
            df.loc[index, 'fastfit_chiProb'] = 1.0 - scipy.stats.chi2.cdf(fitter.getChi2(), fitter.getNDF())

            for i in range(n):
                momentum = fitter.getDaughterMomentum(i)
                df.loc[index, 'fastfit_px_{}'.format(i)] = momentum[0]
                df.loc[index, 'fastfit_py_{}'.format(i)] = momentum[1]
                df.loc[index, 'fastfit_pz_{}'.format(i)] = momentum[2]

        r = (-0.05, 0.05)
        f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
        f.suptitle("Deviation from MC truth for fit with {} daughters".format(n))
        ax1.set_title("X")
        df['mcDX'].hist(bins=100, range=r, label='Trivial', alpha=0.3, ax=ax1, color='g')
        (df['fastfit_dx'] - df['mcDX']).hist(bins=100, range=r, label='FastFit', alpha=0.3, ax=ax1, color='r')
        (df['dx'] - df['mcDX']).hist(bins=100, range=r, label='KFit', alpha=0.3, ax=ax1, color='b')
        ax2.set_title("Y")
        df['mcDY'].hist(bins=100, range=r, label='Trivial', alpha=0.3, ax=ax2, color='g')
        (df['fastfit_dy'] - df['mcDY']).hist(bins=100, range=r, label='FastFit', alpha=0.3, ax=ax2, color='r')
        (df['dy'] - df['mcDY']).hist(bins=100, range=r, label='KFit', alpha=0.3, ax=ax2, color='b')
        ax3.set_title("Z")
        df['mcDZ'].hist(bins=100, range=r, label='Trivial', alpha=0.3, ax=ax3, color='g')
        (df['fastfit_dz'] - df['mcDZ']).hist(bins=100, range=r, label='FastFit', alpha=0.3, ax=ax3, color='r')
        (df['dz'] - df['mcDZ']).hist(bins=100, range=r, label='KFit', alpha=0.3, ax=ax3, color='b')
        plt.legend()
        plt.savefig('MC_D{}.png'.format(n))
        plt.clf()
        
        r = (-0.01, 0.01)
        f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
        f.suptitle("Deviation from KFit for fit with {} daughters".format(n))
        ax1.set_title("X")
        (df['fastfit_dx'] - df['dx']).hist(bins=100, range=r, label='FastFit', alpha=0.3, ax=ax1, color='r')
        ax2.set_title("Y")
        (df['fastfit_dy'] - df['dy']).hist(bins=100, range=r, label='FastFit', alpha=0.3, ax=ax2, color='r')
        ax3.set_title("Z")
        (df['fastfit_dz'] - df['dz']).hist(bins=100, range=r, label='FastFit', alpha=0.3, ax=ax3, color='r')
        plt.savefig('KF_D{}.png'.format(n))
        plt.clf()
        
        r = (0.0, 1.0)
        f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
        f.suptitle("Chi2: Distribution and Deviation from KFit")
        ax1.set_title("Chi2")
        df['fastfit_chiProb'].hist(bins=100, range=r, label='FastFit', alpha=0.3, ax=ax1, color='r')
        df['chiProb'].hist(bins=100, range=r, label='KFit', alpha=0.3, ax=ax1, color='b')
        ax2.set_title("Chi")
        (df['fastfit_chiProb'] - df['chiProb']).hist(bins=100, range=r, label='Deviation', alpha=0.3, ax=ax2, color='r')
        plt.savefig('ChiProb_D{}.png'.format(n))
        plt.clf()
