#ifndef DISPERSION_HH
#define DISPERSION_HH

#include "plasma_params.hh"
#include <complex>
#include <vector>
#include <tuple> // Añadir esta inclusión

class DispersionRelation {
public:
    DispersionRelation(const PlasmaParams& params);
    
    // Relaciones de dispersión para diferentes modos
    std::complex<double> right_hand_mode(double omega) const;
    std::complex<double> left_hand_mode(double omega) const;
    std::complex<double> ordinary_mode(double omega) const;
    std::complex<double> extraordinary_mode(double omega) const;
    
    // Calcular curvas de dispersión (nueva versión que devuelve resultados)
    std::tuple<std::vector<double>, 
               std::vector<std::complex<double>>,
               std::vector<std::complex<double>>,
               std::vector<std::complex<double>>,
               std::vector<std::complex<double>>> 
    calculate_dispersion_curves(double omega_min, double omega_max, int points) const;
    
private:
    PlasmaParams params_;
};

#endif
