/**
 * Thomas Keck 2016
*/

#include "FastFit.h"

extern "C" {

    void* Create(unsigned int numberOfDaughters, double magnetic_field) {
        FastFit *ff = new(std::nothrow) FastFit(numberOfDaughters, magnetic_field);
        return ff;
    }

    void Delete(void* ptr) {
        delete reinterpret_cast<FastFit*>(ptr);
    }

    bool fit(void* ptr, unsigned int maximumNumberOfFitIterations) {
        FastFit *ff = reinterpret_cast<FastFit*>(ptr);
        return ff->fit(maximumNumberOfFitIterations);
    }
    
    void SetDaughter(void* ptr, unsigned int i, int charge, double* momentum, double* vertex, double* error) {
        FastFit *ff = reinterpret_cast<FastFit*>(ptr);
        std::vector<std::vector<double>> error_2d(7);
        for(unsigned int j = 0; j < 6; ++j) {
            std::vector<double> temp(6);
            for(unsigned int k = 0; k < 6; ++k) {
                temp[k] = error[j*6 + k];
            }
            error_2d[j] = temp;
        }
        ff->SetDaughter(i, charge, momentum, vertex, error_2d);
    }

    double getChi2(void* ptr) {
        FastFit *ff = reinterpret_cast<FastFit*>(ptr);
        return ff->getChi2();
    }

    unsigned int getNDF(void* ptr) {
        FastFit *ff = reinterpret_cast<FastFit*>(ptr);
        return ff->getNDF();
    }

    double GetDaughterMomentum(void* ptr, unsigned int i, unsigned int component) {
        FastFit *ff = reinterpret_cast<FastFit*>(ptr);
        return ff->GetDaughterMomentum(i, component);
    }

    double GetDaughterVariance(void* ptr, unsigned int i, unsigned int component_i, unsigned int component_j) {
        FastFit *ff = reinterpret_cast<FastFit*>(ptr);
        return ff->GetDaughterVariance(i, component_i, component_j);
    }

    double GetVariance(void* ptr, unsigned int component_i, unsigned int component_j) {
        FastFit *ff = reinterpret_cast<FastFit*>(ptr);
        return ff->GetVariance(component_i, component_j);
    }

    double GetVertex(void* ptr, unsigned int component) {
        FastFit *ff = reinterpret_cast<FastFit*>(ptr);
        return ff->GetVertex(component);
    }

}
