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
double PlasmaParams::ion_plasma_frequency() const {
    return std::sqrt(ion_density * ELEMENTARY_CHARGE * ELEMENTARY_CHARGE / 
                    (ion_mass * VACUUM_PERMITTIVITY));
}

double PlasmaParams::ion_cyclotron_frequency() const {
    return ELEMENTARY_CHARGE * magnetic_field / ion_mass;
}

std::complex<double> DispersionRelation::right_hand_mode(double omega) const {
    double omega_pe = params_.electron_plasma_frequency();
    double omega_ce = params_.electron_cyclotron_frequency();
    double omega_pi = params_.ion_plasma_frequency();
    double omega_ci = params_.ion_cyclotron_frequency();
    double nu_e = params_.collision_frequency;
    double nu_i = params_.ion_collision_frequency;
    
    // Usar frecuencias complejas para incluir colisiones
    std::complex<double> omega_e(omega, -nu_e);
    std::complex<double> omega_i(omega, -nu_i);
    
    // Contribución de electrones e iones
    std::complex<double> electron_term = (omega_pe * omega_pe) / 
                                       (omega_e * (omega_e + omega_ce));
    std::complex<double> ion_term = (omega_pi * omega_pi) / 
                                  (omega_i * (omega_i - omega_ci));
    
    std::complex<double> k = omega / params_.LIGHT_SPEED * 
                            std::sqrt(1.0 - electron_term - ion_term);
    
    return k;
}

std::complex<double> DispersionRelation::left_hand_mode(double omega) const {
    double omega_pe = params_.electron_plasma_frequency();
    double omega_ce = params_.electron_cyclotron_frequency();
    double omega_pi = params_.ion_plasma_frequency();
    double omega_ci = params_.ion_cyclotron_frequency();
    double nu_e = params_.collision_frequency;
    double nu_i = params_.ion_collision_frequency;
    
    std::complex<double> omega_e(omega, -nu_e);
    std::complex<double> omega_i(omega, -nu_i);
    
    std::complex<double> electron_term = (omega_pe * omega_pe) / 
                                       (omega_e * (omega_e - omega_ce));
    std::complex<double> ion_term = (omega_pi * omega_pi) / 
                                  (omega_i * (omega_i + omega_ci));
    
    std::complex<double> k = omega / params_.LIGHT_SPEED * 
                            std::sqrt(1.0 - electron_term - ion_term);
    
    return k;
}

std::complex<double> DispersionRelation::ordinary_mode(double omega) const {
    double omega_pe = params_.electron_plasma_frequency();
    double omega_pi = params_.ion_plasma_frequency();
    double nu_e = params_.collision_frequency;
    double nu_i = params_.ion_collision_frequency;
    
    std::complex<double> omega_e(omega, -nu_e);
    std::complex<double> omega_i(omega, -nu_i);
    
    // Contribución de ambas especies
    std::complex<double> electron_term = (omega_pe * omega_pe) / (omega_e * omega_e);
    std::complex<double> ion_term = (omega_pi * omega_pi) / (omega_i * omega_i);
    
    std::complex<double> k = omega / params_.LIGHT_SPEED * 
                            std::sqrt(1.0 - electron_term - ion_term);
    
    return k;
}

std::complex<double> DispersionRelation::extraordinary_mode(double omega) const {
    double omega_pe = params_.electron_plasma_frequency();
    double omega_ce = params_.electron_cyclotron_frequency();
    double omega_pi = params_.ion_plasma_frequency();
    double omega_ci = params_.ion_cyclotron_frequency();
    double nu_e = params_.collision_frequency;
    double nu_i = params_.ion_collision_frequency;
    
    std::complex<double> omega_e(omega, -nu_e);
    std::complex<double> omega_i(omega, -nu_i);
    
    // Componentes del tensor dieléctrico con ambas especies
    std::complex<double> S = 1.0 - (omega_pe * omega_pe) / 
                            (omega_e * omega_e - omega_ce * omega_ce) -
                            (omega_pi * omega_pi) / 
                            (omega_i * omega_i - omega_ci * omega_ci);
    
    std::complex<double> D = (omega_ce * omega_pe * omega_pe) / 
                            (omega_e * (omega_e * omega_e - omega_ce * omega_ce)) +
                            (omega_ci * omega_pi * omega_pi) / 
                            (omega_i * (omega_i * omega_i - omega_ci * omega_ci));
    
    std::complex<double> P = 1.0 - (omega_pe * omega_pe) / (omega_e * omega_e) -
                            (omega_pi * omega_pi) / (omega_i * omega_i);
    
    // Fórmula completa para el modo X
    std::complex<double> n_squared = (S * S - D * D) / S;
    std::complex<double> k = (omega / params_.LIGHT_SPEED) * std::sqrt(n_squared);
    
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
