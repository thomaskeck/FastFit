#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from basf2 import *
from modularAnalysis import *

# Define some nicer variables names, so the final pandas data frames are easily readable
import variables as v
variables = ['dx', 'dy', 'dz', 'px', 'py', 'pz', 'mcDX', 'mcDY', 'mcDZ', 'mcE', 'mcPX', 'mcPY', 'mcPZ', 'charge']
for i in range(7):
    for j in range(7):
        v.variables.addAlias('momVertCovM{}{}'.format(i,j), 'momVertCovM({},{})'.format(i,j))
        variables.append('momVertCovM{}{}'.format(i,j))
for daughter in range(5):
    for variable in variables:
        v.variables.addAlias('{}_{}'.format(variable, daughter), 'daughter({}, {})'.format(daughter, variable))

# We load some charged Y4S events
inputMdstList('default', [
    '/storage/jbod/tkeck/MC6/charged/sub00/mdst_000001_prod00000189_task00000001.root',
    '/storage/jbod/tkeck/MC6/charged/sub00/mdst_000003_prod00000189_task00000003.root',
    '/storage/jbod/tkeck/MC6/charged/sub00/mdst_000004_prod00000189_task00000004.root',
    '/storage/jbod/tkeck/MC6/charged/sub00/mdst_000005_prod00000189_task00000005.root',
    '/storage/jbod/tkeck/MC6/charged/sub00/mdst_000007_prod00000189_task00000007.root',
    '/storage/jbod/tkeck/MC6/charged/sub00/mdst_000008_prod00000189_task00000008.root',
    '/storage/jbod/tkeck/MC6/charged/sub00/mdst_000009_prod00000189_task00000009.root',
    '/storage/jbod/tkeck/MC6/charged/sub00/mdst_000011_prod00000189_task00000011.root',
    '/storage/jbod/tkeck/MC6/charged/sub00/mdst_000012_prod00000189_task00000012.root',
    '/storage/jbod/tkeck/MC6/charged/sub00/mdst_000014_prod00000189_task00000014.root',
    '/storage/jbod/tkeck/MC6/mixed/sub00/mdst_000001_prod00000188_task00000001.root',
    '/storage/jbod/tkeck/MC6/mixed/sub00/mdst_000003_prod00000188_task00000003.root',
    '/storage/jbod/tkeck/MC6/mixed/sub00/mdst_000004_prod00000188_task00000004.root',
    '/storage/jbod/tkeck/MC6/mixed/sub00/mdst_000006_prod00000188_task00000006.root',
    '/storage/jbod/tkeck/MC6/mixed/sub00/mdst_000008_prod00000188_task00000008.root',
    '/storage/jbod/tkeck/MC6/mixed/sub00/mdst_000009_prod00000188_task00000009.root',
    '/storage/jbod/tkeck/MC6/mixed/sub00/mdst_000010_prod00000188_task00000010.root',
    '/storage/jbod/tkeck/MC6/mixed/sub00/mdst_000011_prod00000188_task00000011.root',
    '/storage/jbod/tkeck/MC6/mixed/sub00/mdst_000012_prod00000188_task00000012.root',
    '/storage/jbod/tkeck/MC6/mixed/sub00/mdst_000013_prod00000188_task00000013.root',
    ])
analysis_main.add_module('ProgressBar')

fillParticleLists([('K+', 'Kid > 0.5'), ('pi+', '')])

# We try to reconstruct 6 different D Mesons, with 2 to 5 daughters
pLists = ['D0:2', 'D+:3', 'D0:4', 'D+:5']

for i, pList in enumerate(pLists):
    reconstructDecay(pList + ' -> K- pi+ ' + ' '.join(['pi+' if j % 2 == 0 else 'pi-' for j in range(i)]), '1.5 < M < 2.0')
    matchMCTruth(pList)
    applyCuts(pList, "isSignal == 1")
    vertexKFit(pList, -2)
    ntuple_variables = ['M', 'chiProb'] + variables
    for j in range(i+2):
        ntuple_variables += [v + '_{}'.format(j) for v in variables]
    variablesToNTuple(pList, ntuple_variables, filename=pList + '.root')

process(analysis_main)
print(statistics)

import root_pandas
for pList in pLists:
    df = root_pandas.read_root(pList + '.root')
    df.to_pickle(pList + '.pickle')
