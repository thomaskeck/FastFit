#!/usr/bin/env python3
# Thomas Keck

from PyFastFit import FastFit
import numpy as np
import pandas
import os
import matplotlib.pyplot as plt

if __name__ == '__main__':

    files = [(2, 'D0:2.pickle'), (3, 'D+:3.pickle'), (4, 'D0:4.pickle'), (5, 'D+:5.pickle')]

    for n, f in files:
        df = pandas.read_pickle(os.path.join('data', f))
        fitter = FastFit.FastFit(n, 1.5)

        df['fastfit_dx'] = 0.0
        df['fastfit_dy'] = 0.0
        df['fastfit_dz'] = 0.0

        for i in range(n):
            df['fastfit_px_{}'.format(i)] = 0.0
            df['fastfit_py_{}'.format(i)] = 0.0
            df['fastfit_pz_{}'.format(i)] = 0.0

        # TODO Variance is not yet correctly calculated by FastFit
        for j in range(7):
            for k in range(7):
                df['fastfit_momVertCovM_{}{}'.format(j, k)] = 0.0

        for index, decay in df.iterrows():
            for i in range(n):
                charge = int(decay['charge_{}'.format(i)])
                momentum = np.array([decay['px_{}'.format(i)], decay['py_{}'.format(i)], decay['pz_{}'.format(i)]])
                position = np.array([decay['dx_{}'.format(i)], decay['dy_{}'.format(i)], decay['dz_{}'.format(i)]])
                variance = np.zeros((7, 7))
                for j in range(7):
                    for k in range(7):
                        variance[j, k] = decay['momVertCovM_{}{}_{}'.format(j, k, i)]
                fitter.setDaughter(i, charge, momentum, position, variance)

            fitter.fit(3)

            vertex = fitter.getVertex()
            decay['fastfit_x'] = vertex[0]
            decay['fastfit_y'] = vertex[1]
            decay['fastfit_z'] = vertex[2]

            for i in range(n):
                momentum = fitter.getDaughterMomentum(i)
                df['fastfit_px_{}'.format(i)] = momentum[0]
                df['fastfit_py_{}'.format(i)] = momentum[1]
                df['fastfit_pz_{}'.format(i)] = momentum[2]

        f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
        ax1.set_title("X")
        (df['fastfit_dx'] - df['mcDX']).hist(bins=100, range=(-0.2, 0.2), label='FastFit', alpha=0.5, ax=ax1, c='r')
        (df['dx'] - df['mcDX']).hist(bins=100, range=(-0.2, 0.2), label='KFit', alpha=0.5, ax=ax1, c='b')
        ax2.set_title("Y")
        (df['fastfit_dy'] - df['mcDY']).hist(bins=100, range=(-0.2, 0.2), label='FastFit', alpha=0.5, ax=ax2, c='r')
        (df['dy'] - df['mcDY']).hist(bins=100, range=(-0.2, 0.2), label='KFit', alpha=0.5, ax=ax2, c='b')
        ax3.set_title("Z")
        (df['fastfit_dz'] - df['mcDZ']).hist(bins=100, range=(-0.2, 0.2), label='FastFit', alpha=0.5, ax=ax2, c='r')
        (df['dz'] - df['mcDZ']).hist(bins=100, range=(-0.2, 0.2), label='KFit', alpha=0.5, ax=ax2, c='b')
        #plt.savefig(f + '.png')
        plt.show()
        plt.clf()
