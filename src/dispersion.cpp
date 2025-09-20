#include "dispersion.hh"
#include <cmath>
#include <iostream>
#include <tuple>

DispersionRelation::DispersionRelation(const PlasmaParams& params)
    : params_(params) {}

double PlasmaParams::electron_plasma_frequency() const {
    return std::sqrt(electron_density * ELEMENTARY_CHARGE * ELEMENTARY_CHARGE / 
                    (ELECTRON_MASS * VACUUM_PERMITTIVITY));
}

double PlasmaParams::electron_cyclotron_frequency() const {
    return ELEMENTARY_CHARGE * magnetic_field / ELECTRON_MASS;
}

std::complex<double> DispersionRelation::right_hand_mode(double omega) const {
    double omega_pe = params_.electron_plasma_frequency();
    double omega_ce = params_.electron_cyclotron_frequency();
    double nu = params_.collision_frequency;
    
    std::complex<double> denominator(omega * (omega + omega_ce), nu * omega);
    std::complex<double> k = omega / params_.LIGHT_SPEED * 
                            std::sqrt(1.0 - (omega_pe * omega_pe) / denominator);
    
    return k;
}

std::complex<double> DispersionRelation::left_hand_mode(double omega) const {
    double omega_pe = params_.electron_plasma_frequency();
    double omega_ce = params_.electron_cyclotron_frequency();
    double nu = params_.collision_frequency;
    
    std::complex<double> denominator(omega * (omega - omega_ce), nu * omega);
    std::complex<double> k = omega / params_.LIGHT_SPEED * 
                            std::sqrt(1.0 - (omega_pe * omega_pe) / denominator);
    
    return k;
}

std::complex<double> DispersionRelation::ordinary_mode(double omega) const {
    double omega_pe = params_.electron_plasma_frequency();
    double nu = params_.collision_frequency;
    
    std::complex<double> denominator(omega * omega, nu * omega);
    std::complex<double> k = omega / params_.LIGHT_SPEED * 
                            std::sqrt(1.0 - (omega_pe * omega_pe) / denominator);
    
    return k;
}

std::complex<double> DispersionRelation::extraordinary_mode(double omega) const {
    double omega_pe = params_.electron_plasma_frequency();
    double omega_ce = params_.electron_cyclotron_frequency();
    double nu = params_.collision_frequency;
    
    // Relación más precisa para el modo X
    std::complex<double> omega_complex(omega, -nu); // ω → ω - iν
    
    std::complex<double> S = 1.0 - (omega_pe * omega_pe) / 
                            (omega_complex * omega_complex - omega_ce * omega_ce);
    std::complex<double> D = (omega_ce * omega_pe * omega_pe) / 
                            (omega_complex * (omega_complex * omega_complex - omega_ce * omega_ce));
    std::complex<double> P = 1.0 - (omega_pe * omega_pe) / (omega_complex * omega_complex);
    
    // Relación de dispersión completa para el modo X
    std::complex<double> numerator = S * S - D * D;
    std::complex<double> denominator = S;
    
    std::complex<double> n_squared = (numerator / denominator); // Índice de refracción al cuadrado
    std::complex<double> k = (omega_complex / params_.LIGHT_SPEED) * std::sqrt(n_squared);
    
    return k;
}

// Nueva implementación de calculate_dispersion_curves
std::tuple<std::vector<double>, 
           std::vector<std::complex<double>>,
           std::vector<std::complex<double>>,
           std::vector<std::complex<double>>,
           std::vector<std::complex<double>>> 
DispersionRelation::calculate_dispersion_curves(double omega_min, double omega_max, int points) const {
    // Crear vectores locales para los resultados
    std::vector<double> frequencies(points);
    std::vector<std::complex<double>> k_R(points);
    std::vector<std::complex<double>> k_L(points);
    std::vector<std::complex<double>> k_O(points);
    std::vector<std::complex<double>> k_X(points);
    
    double delta_omega = (omega_max - omega_min) / (points - 1);
    
    for (int i = 0; i < points; ++i) {
        double omega = omega_min + i * delta_omega;
        frequencies[i] = omega;
        k_R[i] = right_hand_mode(omega);
        k_L[i] = left_hand_mode(omega);
        k_O[i] = ordinary_mode(omega);
        k_X[i] = extraordinary_mode(omega);
    }
    
    // Devolver los resultados como una tupla
    return std::make_tuple(frequencies, k_R, k_L, k_O, k_X);
}
