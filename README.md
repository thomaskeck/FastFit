# FastFit
A fast vertex fitter based on the eigen library and the kalman filter aglorithm

This is mostly a test-project for a faster vertex fitter for the Belle II experiment

# Installation
   * cmake .
   * make
   * make install


# Installation for the Belle 2 Analysis Software Framework
   * Setup basf2
   * git clone https://github.com/thomaskeck/FastFit.git
   * cd FastFit
   * cmake -DCMAKE_INSTALL_PREFIX:PATH=$BELLE2_EXTERNALS_DIR/$BELLE2_EXTERNALS_SUBDIR .
   * make
   * make install
   * cd $BELLE2_LOCAL_DIR
   * scons


# Other implementations
The Belle and Belle II experiment employ the following vertex fitter implementations
   * RAVE https://github.com/newtrino/rave
   * KFitter (no online resource)

# Further Reading
This work is based on
   * http://web-docs.gsi.de/~ikisel/reco/Methods/Billoir_VertexFitting_NIM_A311_1992.pdf (See also Erratum from 1994. There was a sign error.)
   * http://web-docs.gsi.de/~ikisel/reco/Methods/Fruehwirth-Kalman-NIMA262-1987.pdf
   * http://www-hera-b.desy.de/general/publications/pub_notes/95-013.ps.gz
