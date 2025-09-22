#include "simulator_1d.hh"
#include "dispersion.hh"
#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

Simulator1D::Simulator1D(const PlasmaParams& params)
    : SimulatorBase(params)  // Llama al constructor de la clase base
{
    // Inicializaciones adicionales si son necesarias
}

void Simulator1D::initialize() {
    // Inicializar la malla espacial
    z_grid_.resize(params_.grid_points);
    double dz = params_.length / (params_.grid_points - 1);
    for (int i = 0; i < params_.grid_points; ++i) {
        z_grid_[i] = i * dz;
    }

    // Inicializar campos y corrientes a cero
    fields_.resize(6, std::vector<double>(params_.grid_points, 0.0));
    currents_.resize(3, std::vector<double>(params_.grid_points, 0.0));
    currents_ion_.resize(3, std::vector<double>(params_.grid_points, 0.0));
    fields_prev_.resize(6, std::vector<double>(params_.grid_points, 0.0));
    currents_prev_.resize(3, std::vector<double>(params_.grid_points, 0.0));
    currents_ion_prev_.resize(3, std::vector<double>(params_.grid_points, 0.0));
}

void Simulator1D::excite_mode(const std::string& mode_type, double frequency, double amplitude, double time) {
    double omega = 2.0 * M_PI * frequency;
    double c = params_.LIGHT_SPEED;

    // Excitar el modo apropiado
    if (mode_type == "R") {
        fields_[0][0] = amplitude * std::cos(omega * time);
        fields_[1][0] = amplitude * std::sin(omega * time);
        fields_[3][0] = -amplitude * std::sin(omega * time) / c;
        fields_[4][0] = amplitude * std::cos(omega * time) / c;
    } else if (mode_type == "L") {
        fields_[0][0] = amplitude * std::cos(omega * time);
        fields_[1][0] = -amplitude * std::sin(omega * time);
        fields_[3][0] = amplitude * std::sin(omega * time) / c;
        fields_[4][0] = amplitude * std::cos(omega * time) / c;
    } else if (mode_type == "O") {
        fields_[2][0] = amplitude * std::cos(omega * time);
    } else if (mode_type == "X") {
        fields_[0][0] = amplitude * std::cos(omega * time);
        fields_[3][0] = amplitude * std::sin(omega * time);
    }
}

void Simulator1D::run_timesteps(int num_steps, double dt) {
    double current_time = 0.0;
    double c = params_.LIGHT_SPEED;
    double dz = z_grid_[1] - z_grid_[0];
    
    // Verificar condición CFL
    double cfl_dt = dz / c;
    if (dt > cfl_dt * 0.1) {
        std::cout << "Advertencia: dt = " << dt << " puede ser demasiado grande. CFL recomienda dt < " << cfl_dt * 0.1 << std::endl;
    }

    // Archivo para guardar evolución de energía
    std::ofstream energy_file("data/energy_evolution.csv");
    energy_file << "time,em_energy,electron_energy,ion_energy,total_energy\n";
    
    // Energía inicial
    auto initial_energy = calculate_energy_densities();
    energy_file << current_time << "," << initial_energy[0] << "," 
                << initial_energy[1] << "," << initial_energy[2] << "," 
                << initial_energy[3] << "\n";

    for (int step = 0; step < num_steps; ++step) {
        if (!current_mode_.empty()) {
            excite_mode(current_mode_, current_frequency_, current_amplitude_, current_time);
        }

        // Guardar snapshot en intervalos regulares
        if (step_count_ % save_interval_ == 0) {
            std::string filename = "data/snap/field_data_" + current_mode_ + "_" + 
                                  std::to_string(step_count_ / save_interval_) + ".csv";
            export_field_data(filename);
        }

        if (step_count_ % save_interval_ == 0) {
            // Guardar energía cada cierto intervalo
            auto energies = calculate_energy_densities();
            energy_file << current_time << "," << energies[0] << "," 
                       << energies[1] << "," << energies[2] << "," 
                       << energies[3] << "\n";
        }
        
        update_system_rk4(dt);
        apply_collisions(dt);
        apply_boundary_conditions_pml();
        
        current_time += dt;
        step_count_++;
    }
    
    // Guardar el estado final
    std::string filename = "data/field_data_" + current_mode_ + "_final.csv";
    export_field_data(filename);
}

void Simulator1D::export_field_data_binary(const std::string& filename) const {
    std::ofstream file(filename, std::ios::binary);
    
    // Escribir dimensiones
    int rows = params_.grid_points;
    int cols = 10;
    file.write(reinterpret_cast<const char*>(&rows), sizeof(int));
    file.write(reinterpret_cast<const char*>(&cols), sizeof(int));
    
    // Escribir datos
    for (int i = 0; i < params_.grid_points; ++i) {
        double z_val = z_grid_[i];
        file.write(reinterpret_cast<const char*>(&z_val), sizeof(double));
        
        for (int j = 0; j < 6; ++j) {
            file.write(reinterpret_cast<const char*>(&fields_[j][i]), sizeof(double));
        }
        
        for (int j = 0; j < 3; ++j) {
            file.write(reinterpret_cast<const char*>(&currents_[j][i]), sizeof(double));
        }
    }
    
    file.close();
}

void Simulator1D::compute_derivatives(
    const std::vector<std::vector<double>>& fields,
    const std::vector<std::vector<double>>& currents,
    const std::vector<std::vector<double>>& currents_ion,
    std::vector<std::vector<double>>& dfields_dt,
    std::vector<std::vector<double>>& dcurrents_dt,
    std::vector<std::vector<double>>& dcurrents_dt_ion) {

    double c = params_.LIGHT_SPEED;
    double epsilon0 = params_.VACUUM_PERMITTIVITY;
    double dz = z_grid_[1] - z_grid_[0];
    
    double omega_pe = params_.electron_plasma_frequency();
    double omega_ce = params_.electron_cyclotron_frequency();
    double nu = params_.collision_frequency;
    double omega_pi = params_.ion_plasma_frequency();
    double omega_ci = params_.ion_cyclotron_frequency();
    double nu_i = params_.ion_collision_frequency;

    // Inicializar derivadas a cero
    for (auto& vec : dfields_dt) std::fill(vec.begin(), vec.end(), 0.0);
    for (auto& vec : dcurrents_dt) std::fill(vec.begin(), vec.end(), 0.0);
    for (auto& vec : dcurrents_dt_ion) std::fill(vec.begin(), vec.end(), 0.0);

    // Diferencias finitas de 4to orden para derivadas espaciales
    for (int i = 2; i < params_.grid_points - 2; ++i) {
        // Calcular derivadas espaciales
        double dEx_dz = (-fields[0][i+2] + 8*fields[0][i+1] - 8*fields[0][i-1] + fields[0][i-2]) / (12.0 * dz);
        double dEy_dz = (-fields[1][i+2] + 8*fields[1][i+1] - 8*fields[1][i-1] + fields[1][i-2]) / (12.0 * dz);
        double dBx_dz = (-fields[3][i+2] + 8*fields[3][i+1] - 8*fields[3][i-1] + fields[3][i-2]) / (12.0 * dz);
        double dBy_dz = (-fields[4][i+2] + 8*fields[4][i+1] - 8*fields[4][i-1] + fields[4][i-2]) / (12.0 * dz);

        // Ecuaciones de Maxwell para 1D
        dfields_dt[0][i] = c * c * dBy_dz -  (currents[0][i] + currents_ion[0][i]) / epsilon0;
        dfields_dt[1][i] = -c * c * dBx_dz - (currents[1][i] + currents_ion[1][i]) / epsilon0;
        dfields_dt[2][i] = -currents[2][i] / epsilon0;

        dfields_dt[3][i] = -dEy_dz;
        dfields_dt[4][i] = dEx_dz;
        dfields_dt[5][i] = 0.0;

        // Ecuaciones de momento para corrientes electronicas
        dcurrents_dt[0][i] = epsilon0 * omega_pe * omega_pe * fields[0][i] + 
            omega_ce * currents[1][i] - nu * currents[0][i];
        dcurrents_dt[1][i] = epsilon0 * omega_pe * omega_pe * fields[1][i] - 
            omega_ce * currents[0][i] - nu * currents[1][i];
        dcurrents_dt[2][i] = epsilon0 * omega_pe * omega_pe * fields[2][i] - 
            nu * currents[2][i];
        
        // Ecuaciones de momento para iones
        dcurrents_dt_ion[0][i] = epsilon0 * omega_pi * omega_pi * fields[0][i] - 
                                omega_ci * currents_ion[1][i] - nu_i * currents_ion[0][i];
        dcurrents_dt_ion[1][i] = epsilon0 * omega_pi * omega_pi * fields[1][i] + 
                                omega_ci * currents_ion[0][i] - nu_i * currents_ion[1][i];
        dcurrents_dt_ion[2][i] = epsilon0 * omega_pi * omega_pi * fields[2][i] - 
                                nu_i * currents_ion[2][i];
    }
}

void Simulator1D::update_system_rk4(double dt) {
    // Vectores para almacenar las k (pendientes) de RK4
    std::vector<std::vector<double>> k1_fields(6, std::vector<double>(params_.grid_points));
    std::vector<std::vector<double>> k2_fields(6, std::vector<double>(params_.grid_points));
    std::vector<std::vector<double>> k3_fields(6, std::vector<double>(params_.grid_points));
    std::vector<std::vector<double>> k4_fields(6, std::vector<double>(params_.grid_points));
    
    std::vector<std::vector<double>> k1_currents(3, std::vector<double>(params_.grid_points));
    std::vector<std::vector<double>> k2_currents(3, std::vector<double>(params_.grid_points));
    std::vector<std::vector<double>> k3_currents(3, std::vector<double>(params_.grid_points));
    std::vector<std::vector<double>> k4_currents(3, std::vector<double>(params_.grid_points));
    
    std::vector<std::vector<double>> k1_currents_ion(3, std::vector<double>(params_.grid_points));
    std::vector<std::vector<double>> k2_currents_ion(3, std::vector<double>(params_.grid_points));
    std::vector<std::vector<double>> k3_currents_ion(3, std::vector<double>(params_.grid_points));
    std::vector<std::vector<double>> k4_currents_ion(3, std::vector<double>(params_.grid_points));
    
    // Vectores temporales para almacenar estados intermedios
    std::vector<std::vector<double>> fields_temp = fields_;
    std::vector<std::vector<double>> currents_temp = currents_;
    std::vector<std::vector<double>> currents_ion_temp = currents_ion_;

    // --- Primer paso de RK4 (k1) ---
    compute_derivatives(fields_, currents_, currents_ion_, k1_fields, k1_currents, k1_currents_ion);

    // --- Segundo paso de RK4 (k2) ---
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < params_.grid_points; ++j) {
            fields_temp[i][j] = fields_[i][j] + 0.5 * dt * k1_fields[i][j];
        }
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < params_.grid_points; ++j) {
            currents_temp[i][j] = currents_[i][j] + 0.5 * dt * k1_currents[i][j];
        }
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < params_.grid_points; ++j) {
            currents_ion_temp[i][j] = currents_ion_[i][j] + 0.5 * dt * k1_currents_ion[i][j];
        }
    }
    compute_derivatives(fields_temp, currents_temp, currents_ion_temp, k2_fields, k2_currents, k2_currents_ion);

    // --- Tercer paso de RK4 (k3) ---
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < params_.grid_points; ++j) {
            fields_temp[i][j] = fields_[i][j] + 0.5 * dt * k2_fields[i][j];
        }
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < params_.grid_points; ++j) {
            currents_temp[i][j] = currents_[i][j] + 0.5 * dt * k2_currents[i][j];
        }
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < params_.grid_points; ++j) {
            currents_ion_temp[i][j] = currents_ion_[i][j] + 0.5 * dt * k2_currents_ion[i][j];
        }
    }
    compute_derivatives(fields_temp, currents_temp, currents_ion_temp, k3_fields, k3_currents, k3_currents_ion);

    // --- Cuarto paso de RK4 (k4) ---
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < params_.grid_points; ++j) {
            fields_temp[i][j] = fields_[i][j] + dt * k3_fields[i][j];
        }
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < params_.grid_points; ++j) {
            currents_temp[i][j] = currents_[i][j] + dt * k3_currents[i][j];
        }
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < params_.grid_points; ++j) {
            currents_ion_temp[i][j] = currents_ion_[i][j] + dt * k3_currents_ion[i][j];
        }
    }
    compute_derivatives(fields_temp, currents_temp, currents_ion_temp, k4_fields, k4_currents, k4_currents_ion);

    // --- Combinar resultados para campos ---
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < params_.grid_points; ++j) {
            fields_[i][j] += dt * (k1_fields[i][j] + 2*k2_fields[i][j] + 
                2*k3_fields[i][j] + k4_fields[i][j]) / 6.0;
        }
    }

    // --- Combinar resultados para corrientes ---
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < params_.grid_points; ++j) {
            currents_[i][j] += dt * (k1_currents[i][j] + 2*k2_currents[i][j] + 
                2*k3_currents[i][j] + k4_currents[i][j]) / 6.0;
        }
    }
    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < params_.grid_points; ++j) {
            currents_ion_[i][j] += dt * (k1_currents_ion[i][j] + 2*k2_currents_ion[i][j] + 
                2*k3_currents_ion[i][j] + k4_currents_ion[i][j]) / 6.0;
        }
    }
    
    // Verificar valores problemáticos
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < params_.grid_points; ++j) {
            if (std::isnan(fields_[i][j]) || std::isinf(fields_[i][j])) {
                fields_[i][j] = 0.0;
            }
        }
    }
    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < params_.grid_points; ++j) {
            if (std::isnan(currents_[i][j]) || std::isinf(currents_[i][j])) {
                currents_[i][j] = 0.0;
            }
        }
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < params_.grid_points; ++j) {
            if (std::isnan(currents_ion_[i][j]) || std::isinf(currents_ion_[i][j])) {
                currents_ion_[i][j] = 0.0;
            }
        }
    }
}

void Simulator1D::apply_collisions(double dt) {
    double nu = params_.collision_frequency;
    double nu_i = params_.ion_collision_frequency;

    // Aplicar amortiguamiento exponencial debido a colisiones
    for (int i = 0; i < params_.grid_points; ++i) {
        currents_[0][i] *= std::exp(-nu * dt);
        currents_[1][i] *= std::exp(-nu * dt);
        currents_[2][i] *= std::exp(-nu * dt);
        
        currents_ion_[0][i] *= std::exp(-nu_i * dt);
        currents_ion_[1][i] *= std::exp(-nu_i * dt);
        currents_ion_[2][i] *= std::exp(-nu_i * dt);
    }
}

void Simulator1D::apply_boundary_conditions_pml() {
    int pml_width = 10;
    double sigma_max = 0.001;

    for (int i = 0; i < pml_width; ++i) {
        double sigma = sigma_max * std::pow((double)(pml_width - i) / pml_width, 4);

        for (int j = 0; j < 6; ++j) {
            fields_[j][i] *= std::exp(-sigma);
            fields_[j][params_.grid_points - 1 - i] *= std::exp(-sigma);
        }
    }
}

std::vector<double> Simulator1D::calculate_energy_densities() const {
    std::vector<double> energies(4, 0.0); // [EM, electron, ion, total]
    double dz = z_grid_[1] - z_grid_[0];
    
    // Constantes físicas
    double epsilon0 = params_.VACUUM_PERMITTIVITY;
    double mu0 = 1.0 / (epsilon0 * params_.LIGHT_SPEED * params_.LIGHT_SPEED);
    double omega_pe = params_.electron_plasma_frequency();
    double omega_pi = params_.ion_plasma_frequency();
    
    // Verificar que las frecuencias de plasma no sean cero
    if (omega_pe == 0.0) omega_pe = 1e-10; // Evitar división por cero
    if (omega_pi == 0.0) omega_pi = 1e-10;
    
    for (int i = 0; i < params_.grid_points; ++i) {
        // 1. Energía electromagnética por unidad de volumen: (ε₀E² + B²/μ₀)/2
        double E_sq = 0.0, B_sq = 0.0;
        for (int comp = 0; comp < 3; ++comp) {
            E_sq += fields_[comp][i] * fields_[comp][i];
            B_sq += fields_[comp+3][i] * fields_[comp+3][i];
        }
        double em_energy_density = 0.5 * (epsilon0 * E_sq + B_sq / mu0);
        
        // 2. Energía cinética electrónica por unidad de volumen: Jₑ²/(2ε₀ωₚₑ²)
        double Je_sq = 0.0;
        for (int comp = 0; comp < 3; ++comp) {
            Je_sq += currents_[comp][i] * currents_[comp][i];
        }
        double electron_energy_density = Je_sq / (2.0 * epsilon0 * omega_pe * omega_pe);
        
        // 3. Energía cinética iónica por unidad de volumen: Jᵢ²/(2ε₀ωₚᵢ²)
        double Ji_sq = 0.0;
        for (int comp = 0; comp < 3; ++comp) {
            Ji_sq += currents_ion_[comp][i] * currents_ion_[comp][i];
        }
        double ion_energy_density = Ji_sq / (2.0 * epsilon0 * omega_pi * omega_pi);
        
        // Integrar sobre el volumen (1D: multiplicar por dz)
        energies[0] += em_energy_density * dz;
        energies[1] += electron_energy_density * dz;
        energies[2] += ion_energy_density * dz;
    }
    
    energies[3] = energies[0] + energies[1] + energies[2];
    
    return energies;
}

void Simulator1D::export_field_data(const std::string& filename) const {
    std::ofstream file(filename);
    file << "z,Ex,Ey,Ez,Bx,By,Bz,Jx_e,Jy_e,Jz_e,Jx_i,Jy_i,Jz_i\n";

    for (int i = 0; i < params_.grid_points; ++i) {
        file << z_grid_[i];
        for (int j = 0; j < 6; ++j) file << "," << fields_[j][i];
        for (int j = 0; j < 3; ++j) file << "," << currents_[j][i];
        for (int j = 0; j < 3; ++j) file << "," << currents_ion_[j][i];
        file << "\n";
    }

    file.close();
}

void Simulator1D::export_dispersion_data(const std::string& filename) const {
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
