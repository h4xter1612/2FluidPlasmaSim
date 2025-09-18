#include "simulator.hh"
#include "dispersion.hh"
#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

TwoFluidSimulator::TwoFluidSimulator(const PlasmaParams& params)
    : params_(params), dispersion_(params) {}

void TwoFluidSimulator::initialize() {
    // Inicializar la malla espacial
    z_grid_.resize(params_.grid_points);
    double dz = params_.length / (params_.grid_points - 1);
    for (int i = 0; i < params_.grid_points; ++i) {
        z_grid_[i] = i * dz;
    }
    
    // Inicializar campos y corrientes
    fields_.resize(6, std::vector<double>(params_.grid_points, 0.0));
    currents_.resize(3, std::vector<double>(params_.grid_points, 0.0));
}

void TwoFluidSimulator::run_timesteps(int num_steps, double dt) {
    for (int step = 0; step < num_steps; ++step) {
        // Aplicar excitación
        if (!current_mode_.empty()) {
            excite_mode(current_mode_, current_frequency_, current_amplitude_);
        }
        
        // Actualizar campos y corrientes
        update_currents(dt);
        update_fields(dt);
        apply_collisions(dt);
        apply_boundary_conditions();
    }
}

void TwoFluidSimulator::excite_mode(const std::string& mode_type, double frequency, double amplitude) {
    double omega = 2.0 * M_PI * frequency;
    double t = 0.0; // En una implementación real, esto debería ser el tiempo actual
    
    // Excitar el modo apropiado
    if (mode_type == "R") {
        // Modo derecho
        fields_[0][0] = amplitude * std::cos(omega * t); // Ex
        fields_[1][0] = amplitude * std::sin(omega * t); // Ey
    } else if (mode_type == "L") {
        // Modo izquierdo
        fields_[0][0] = amplitude * std::cos(omega * t); // Ex
        fields_[1][0] = -amplitude * std::sin(omega * t); // Ey
    } else if (mode_type == "O") {
        // Modo ordinario
        fields_[2][0] = amplitude * std::cos(omega * t); // Ez
    } else if (mode_type == "X") {
        // Modo extraordinario
        fields_[0][0] = amplitude * std::cos(omega * t); // Ex
        fields_[3][0] = amplitude * std::sin(omega * t); // Bx
    }
}

void TwoFluidSimulator::update_fields(double dt) {
    // Implementación simplificada de la actualización de campos
    // En una implementación real, se usarían las ecuaciones de Maxwell
    double c = params_.LIGHT_SPEED;
    double epsilon0 = params_.VACUUM_PERMITTIVITY;
    
    for (int i = 1; i < params_.grid_points - 1; ++i) {
        // Actualizar campos eléctricos
        fields_[0][i] += dt * (c * c * (fields_[5][i+1] - fields_[5][i-1]) / (2.0 * (z_grid_[1] - z_grid_[0])) - 
                              currents_[0][i] / epsilon0);
        fields_[1][i] += dt * (c * c * (fields_[3][i+1] - fields_[3][i-1]) / (2.0 * (z_grid_[1] - z_grid_[0])) - 
                              currents_[1][i] / epsilon0);
        fields_[2][i] += dt * (-currents_[2][i] / epsilon0);
        
        // Actualizar campos magnéticos
        fields_[3][i] += dt * ((fields_[1][i+1] - fields_[1][i-1]) / (2.0 * (z_grid_[1] - z_grid_[0])));
        fields_[4][i] += dt * (-(fields_[0][i+1] - fields_[0][i-1]) / (2.0 * (z_grid_[1] - z_grid_[0])));
        fields_[5][i] += dt * ((fields_[0][i+1] - fields_[0][i-1]) / (2.0 * (z_grid_[1] - z_grid_[0])));
    }
}

void TwoFluidSimulator::update_currents(double dt) {
    double omega_pe = params_.electron_plasma_frequency();
    double omega_ce = params_.electron_cyclotron_frequency();
    double epsilon0 = params_.VACUUM_PERMITTIVITY;
    
    for (int i = 0; i < params_.grid_points; ++i) {
        // Ecuaciones de momento para electrones
        currents_[0][i] += dt * (epsilon0 * omega_pe * omega_pe * fields_[0][i] + 
                                omega_ce * currents_[1][i]);
        currents_[1][i] += dt * (epsilon0 * omega_pe * omega_pe * fields_[1][i] - 
                                omega_ce * currents_[0][i]);
        currents_[2][i] += dt * (epsilon0 * omega_pe * omega_pe * fields_[2][i]);
    }
}

void TwoFluidSimulator::apply_collisions(double dt) {
    double nu = params_.collision_frequency;
    
    for (int i = 0; i < params_.grid_points; ++i) {
        currents_[0][i] *= std::exp(-nu * dt);
        currents_[1][i] *= std::exp(-nu * dt);
        currents_[2][i] *= std::exp(-nu * dt);
    }
}

void TwoFluidSimulator::apply_boundary_conditions() {
    // Condiciones de contorno periódicas
    for (int i = 0; i < 6; ++i) {
        fields_[i][0] = fields_[i][params_.grid_points - 2];
        fields_[i][params_.grid_points - 1] = fields_[i][1];
    }
    
    for (int i = 0; i < 3; ++i) {
        currents_[i][0] = currents_[i][params_.grid_points - 2];
        currents_[i][params_.grid_points - 1] = currents_[i][1];
    }
}

void TwoFluidSimulator::export_field_data(const std::string& filename) const {
    std::ofstream file(filename);
    file << "z,Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz\n";
    
    for (int i = 0; i < params_.grid_points; ++i) {
        file << z_grid_[i];
        for (int j = 0; j < 6; ++j) file << "," << fields_[j][i];
        for (int j = 0; j < 3; ++j) file << "," << currents_[j][i];
        file << "\n";
    }
    
    file.close();
}
void TwoFluidSimulator::export_dispersion_data(const std::string& filename) const {
    // Calcular curvas de dispersión y obtener los resultados
    auto [frequencies, k_R, k_L, k_O, k_X] = dispersion_.calculate_dispersion_curves(1e6, 1e11, 1000);
    
    std::ofstream file(filename);
    file << "frequency,Re(k_R),Im(k_R),Re(k_L),Im(k_L),Re(k_O),Im(k_O),Re(k_X),Im(k_X)\n";
    
    for (size_t i = 0; i < frequencies.size(); ++i) {
        file << frequencies[i] << "," 
             << k_R[i].real() << "," << k_R[i].imag() << ","
             << k_L[i].real() << "," << k_L[i].imag() << ","
             << k_O[i].real() << "," << k_O[i].imag() << ","
             << k_X[i].real() << "," << k_X[i].imag() << "\n";
    }
    
    file.close();
}
