/**
 * Thomas Keck 2016
 */

extern "C" {

    void* Create(unsigned int numberOfDaughters);

    void Delete(void*);

    bool fit(void*, unsigned int maximumNumberOfFitIterations, double magnetic_field);
    
    void SetDaughter(void*, unsigned int i, int charge, double* momentum, double* vertex, double* error);

    double getChi2(void*);

    unsigned int getNDF(void*);

    double GetDaughterMomentum(void*, unsigned int i, unsigned int component);

    double GetDaughterVariance(void*, unsigned int i, unsigned int component_i, unsigned int component_j);

    double GetVertex(void*, unsigned int component);

}
