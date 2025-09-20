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
    : params_(params), dispersion_(params), step_count_(0) {}

void TwoFluidSimulator::initialize() {
    // Inicializar la malla espacial
    z_grid_.resize(params_.grid_points);
    double dz = params_.length / (params_.grid_points - 1);
    for (int i = 0; i < params_.grid_points; ++i) {
        z_grid_[i] = i * dz;
    }

    // Inicializar campos y corrientes a cero
    fields_.resize(6, std::vector<double>(params_.grid_points, 0.0));
    currents_.resize(3, std::vector<double>(params_.grid_points, 0.0));
    currents_ion_.resize(3, std::vector<double>(params_.grid_points, 0.0)); // Nuevo
    fields_prev_.resize(6, std::vector<double>(params_.grid_points, 0.0));
    currents_prev_.resize(3, std::vector<double>(params_.grid_points, 0.0));
    currents_ion_prev_.resize(3, std::vector<double>(params_.grid_points, 0.0)); // Nuevo
}

void TwoFluidSimulator::excite_mode(const std::string& mode_type, double frequency, double amplitude, double time) {
    double omega = 2.0 * M_PI * frequency;
    double c = params_.LIGHT_SPEED;

    // Excitar el modo apropiado
    if (mode_type == "R") {
        // Modo derecho - campos electromagnéticos consistentes
        fields_[0][0] = amplitude * std::cos(omega * time); // Ex
        fields_[1][0] = amplitude * std::sin(omega * time); // Ey
        
        // Campos magnéticos correspondientes (relación E/B = c para onda plana)
        fields_[3][0] = -amplitude * std::sin(omega * time) / c; // Bx
        fields_[4][0] = amplitude * std::cos(omega * time) / c;  // By
        //
        // std::cout << "DEBUG: Inicializando Bx=" << fields_[3][0] 
        //           << ", By=" << fields_[4][0] << std::endl;
    } else if (mode_type == "L") {
        // Modo izquierdo
        fields_[0][0] = amplitude * std::cos(omega * time); // Ex
        fields_[1][0] = -amplitude * std::sin(omega * time); // Ey
        fields_[3][0] = amplitude * std::sin(omega * time) / c; // Bx
        fields_[4][0] = amplitude * std::cos(omega * time) / c; // By
        //
        // std::cout << "DEBUG: Inicializando Bx=" << fields_[3][0]
        //     << ", By=" << fields_[4][0] << std::endl;
    } else if (mode_type == "O") {
        // Modo ordinario (longitudinal)
        fields_[2][0] = amplitude * std::cos(omega * time); // Ez
    } else if (mode_type == "X") {
        // Modo extraordinario
        fields_[0][0] = amplitude * std::cos(omega * time); // Ex
        fields_[3][0] = amplitude * std::sin(omega * time); // Bx
    }
}

void TwoFluidSimulator::run_timesteps(int num_steps, double dt) {
    double current_time = 0.0;
    double c = params_.LIGHT_SPEED;
    double dz = z_grid_[1] - z_grid_[0];
    
    // Verificar condición CFL
    double cfl_dt = dz / c;
    if (dt > cfl_dt * 0.1) {
        std::cout << "Advertencia: dt = " << dt << " puede ser demasiado grande. CFL recomienda dt < " << cfl_dt * 0.1 << std::endl;
    }
    
    for (int step = 0; step < num_steps; ++step) {
        if (!current_mode_.empty()) {
            excite_mode(current_mode_, current_frequency_, current_amplitude_, current_time);
        }
        // Diagnóstico más detallado
        // if (step % 100 == 0) {
        //     std::cout << "Paso " << step << ", T = " << current_time << ": "
        //               << "Ex[0] = " << fields_[0][0] << ", "
        //               << "Ey[0] = " << fields_[1][0] << ", "
        //               << "Bx[0] = " << fields_[3][0] << ", "
        //               << "By[0] = " << fields_[4][0] << ", "
        //               << "Ex[100] = " << fields_[0][100] << ", "
        //               << "Ey[100] = " << fields_[1][100] << std::endl;
        // }
        // // DIAGNÓSTICO: Verificar valores de campos magnéticos
        // if (step % 100 == 0) {
        //     std::cout << "Paso " << step << ": "
        //               << "Bx[0] = " << fields_[3][0] << ", "
        //               << "By[0] = " << fields_[4][0] << ", "
        //               << "Bx[100] = " << fields_[3][100] << ", "
        //               << "By[100] = " << fields_[4][100] << std::endl;
        // }

        // Guardar snapshot en intervalos regulares
        if (step_count_ % save_interval_ == 0) {
            std::string filename = "data/snap/field_data_" + current_mode_ + "_" + 
                                  std::to_string(step_count_ / save_interval_) + ".csv";
            export_field_data(filename);
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

void TwoFluidSimulator::export_field_data_binary(const std::string& filename) const {
    std::ofstream file(filename, std::ios::binary);
    
    // Escribir dimensiones
    int rows = params_.grid_points;
    int cols = 10; // z + 6 campos + 3 corrientes
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

void TwoFluidSimulator::compute_derivatives(
    const std::vector<std::vector<double>>& fields,
    const std::vector<std::vector<double>>& currents,
    const std::vector<std::vector<double>>& currents_ion,
    std::vector<std::vector<double>>& dfields_dt,
    std::vector<std::vector<double>>& dcurrents_dt,
    std::vector<std::vector<double>>& dcurrents_dt_ion) {

    double c = params_.LIGHT_SPEED;
    double epsilon0 = params_.VACUUM_PERMITTIVITY;
    double dz = z_grid_[1] - z_grid_[0];
    // Electronic parameters
    double omega_pe = params_.electron_plasma_frequency();
    double omega_ce = params_.electron_cyclotron_frequency();
    double nu = params_.collision_frequency;
    // Ionic parameters
    double omega_pi = params_.ion_plasma_frequency();
    double omega_ci = params_.ion_cyclotron_frequency();
    double nu_i = params_.ion_collision_frequency;

    // Inicializar derivadas a cero
    for (auto& vec : dfields_dt) std::fill(vec.begin(), vec.end(), 0.0);
    for (auto& vec : dcurrents_dt) std::fill(vec.begin(), vec.end(), 0.0);
    for (auto& vec : dcurrents_dt_ion) std::fill(vec.begin(), vec.end(), 0.0);  // Ionic variable

    // Diferencias finitas de 4to orden para derivadas espaciales
    for (int i = 2; i < params_.grid_points - 2; ++i) {
        // Calcular derivadas espaciales
        double dEx_dz = (-fields[0][i+2] + 8*fields[0][i+1] - 8*fields[0][i-1] + fields[0][i-2]) / (12.0 * dz);
        double dEy_dz = (-fields[1][i+2] + 8*fields[1][i+1] - 8*fields[1][i-1] + fields[1][i-2]) / (12.0 * dz);
        double dBx_dz = (-fields[3][i+2] + 8*fields[3][i+1] - 8*fields[3][i-1] + fields[3][i-2]) / (12.0 * dz);
        double dBy_dz = (-fields[4][i+2] + 8*fields[4][i+1] - 8*fields[4][i-1] + fields[4][i-2]) / (12.0 * dz);

        // Ecuaciones de Maxwell para 1D
        dfields_dt[0][i] = c * c * dBy_dz -  (currents[0][i] + currents_ion[0][i]) / epsilon0;  // ∂Ex/∂t = c²∂By/∂z - Jx/ε₀
        dfields_dt[1][i] = -c * c * dBx_dz - (currents[1][i] + currents_ion[1][i]) / epsilon0; // ∂Ey/∂t = -c²∂Bx/∂z - Jy/ε₀
        dfields_dt[2][i] = -currents[2][i] / epsilon0;                  // ∂Ez/∂t = -Jz/ε₀

        dfields_dt[3][i] = -dEy_dz;  // ∂Bx/∂t = -∂Ey/∂z
        dfields_dt[4][i] = dEx_dz;   // ∂By/∂t = ∂Ex/∂z
        dfields_dt[5][i] = 0.0;      // ∂Bz/∂t = 0

        // Ecuaciones de momento para corrientes electronicas
        dcurrents_dt[0][i] = epsilon0 * omega_pe * omega_pe * fields[0][i] + 
            omega_ce * currents[1][i] - nu * currents[0][i];
        dcurrents_dt[1][i] = epsilon0 * omega_pe * omega_pe * fields[1][i] - 
            omega_ce * currents[0][i] - nu * currents[1][i];
        dcurrents_dt[2][i] = epsilon0 * omega_pe * omega_pe * fields[2][i] - 
            nu * currents[2][i];
        // ECUACIONES DE MOMENTO PARA IONES (nuevas)
        // Nota: signo opuesto para el término de ciclotrón
        dcurrents_dt_ion[0][i] = epsilon0 * omega_pi * omega_pi * fields[0][i] - 
                                omega_ci * currents_ion[1][i] - nu_i * currents_ion[0][i];
        dcurrents_dt_ion[1][i] = epsilon0 * omega_pi * omega_pi * fields[1][i] + 
                                omega_ci * currents_ion[0][i] - nu_i * currents_ion[1][i];
        dcurrents_dt_ion[2][i] = epsilon0 * omega_pi * omega_pi * fields[2][i] - 
                                nu_i * currents_ion[2][i];
            
        // Diagnóstico adicional para detectar problemas
        if (i == 100 && std::isnan(dfields_dt[0][i])) {
            std::cout << "NaN detectado en derivadas: " 
                      << "dBy_dz = " << dBy_dz << ", "
                      << "currents[0][i] = " << currents[0][i] << std::endl;
        }
    }
}

void TwoFluidSimulator::update_system_rk4(double dt) {
    // Vectores para almacenar las k (pendientes) de RK4
    std::vector<std::vector<double>> k1_fields(6, std::vector<double>(params_.grid_points));
    std::vector<std::vector<double>> k2_fields(6, std::vector<double>(params_.grid_points));
    std::vector<std::vector<double>> k3_fields(6, std::vector<double>(params_.grid_points));
    std::vector<std::vector<double>> k4_fields(6, std::vector<double>(params_.grid_points));
    // Electrons
    std::vector<std::vector<double>> k1_currents(3, std::vector<double>(params_.grid_points));
    std::vector<std::vector<double>> k2_currents(3, std::vector<double>(params_.grid_points));
    std::vector<std::vector<double>> k3_currents(3, std::vector<double>(params_.grid_points));
    std::vector<std::vector<double>> k4_currents(3, std::vector<double>(params_.grid_points));
    // Ions
    // Añadir implementación similar para iones:
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
    // Después de actualizar campos y corrientes, verificar valores
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < params_.grid_points; ++j) {
            if (std::isnan(fields_[i][j]) || std::isinf(fields_[i][j])) {
                fields_[i][j] = 0.0;  // Resetear valores problemáticos
            }
        }
    }
    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < params_.grid_points; ++j) {
            if (std::isnan(currents_[i][j]) || std::isinf(currents_[i][j])) {
                currents_[i][j] = 0.0;  // Resetear valores problemáticos
            }
        }
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < params_.grid_points; ++j) {
            if (std::isnan(currents_ion_[i][j]) || std::isinf(currents_ion_[i][j])) {
                currents_ion_[i][j] = 0.0;  // Resetear valores problemáticos
            }
        }
    }

}

void TwoFluidSimulator::apply_collisions(double dt) {
    double nu = params_.collision_frequency;
    double nu_i = params_.ion_collision_frequency;

    // Aplicar amortiguamiento exponencial debido a colisiones
    // Electrons
    for (int i = 0; i < params_.grid_points; ++i) {
        currents_[0][i] *= std::exp(-nu * dt);
        currents_[1][i] *= std::exp(-nu * dt);
        currents_[2][i] *= std::exp(-nu * dt);
    }
    // Ions
    for (int i = 0; i < params_.grid_points; ++i) {
        currents_ion_[0][i] *= std::exp(-nu_i * dt);
        currents_ion_[1][i] *= std::exp(-nu_i * dt);
        currents_ion_[2][i] *= std::exp(-nu_i * dt);
    }
}

void TwoFluidSimulator::apply_boundary_conditions_pml() {
    int pml_width = 10;  // Ancho de la capa PML
    double sigma_max = 0.001;  // Reducir la atenuación máxima

    for (int i = 0; i < pml_width; ++i) {
        // Calcular perfil de conductividad (suave)
        double sigma = sigma_max * std::pow((double)(pml_width - i) / pml_width, 4);

        // Aplicar atenuación PML a los campos (solo en los bordes)
        for (int j = 0; j < 6; ++j) {
            fields_[j][i] *= std::exp(-sigma);  // Borde izquierdo
            fields_[j][params_.grid_points - 1 - i] *= std::exp(-sigma);  // Borde derecho
        }
    }
}

void TwoFluidSimulator::export_field_data(const std::string& filename) const {
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
