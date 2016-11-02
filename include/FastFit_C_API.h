/**
 * Thomas Keck 2016
 */

extern "C" {

    void* Create(unsigned int numberOfDaughters, double magnetic_field);

    void Delete(void*);

    bool fit(void*, unsigned int maximumNumberOfFitIterations);
    
    void SetIPProfile(void*, double* vertex, double* variance);
    
    void SetDaughter(void*, unsigned int i, int charge, double* momentum, double* vertex, double* error);

    double getChi2(void*);

    unsigned int getNDF(void*);

    double GetDaughterMomentum(void*, unsigned int i, unsigned int component);

    double GetDaughterVariance(void*, unsigned int i, unsigned int component_i, unsigned int component_j);
    
    double GetVariance(void*, unsigned int component_i, unsigned int component_j);

    double GetVertex(void*, unsigned int component);

}
