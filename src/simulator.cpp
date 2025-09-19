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

    // Inicializar campos y corrientes a cero
    fields_.resize(6, std::vector<double>(params_.grid_points, 0.0));
    currents_.resize(3, std::vector<double>(params_.grid_points, 0.0));
    fields_prev_.resize(6, std::vector<double>(params_.grid_points, 0.0));
    currents_prev_.resize(3, std::vector<double>(params_.grid_points, 0.0));
}

void TwoFluidSimulator::excite_mode(const std::string& mode_type, double frequency, double amplitude, double time) {
    double omega = 2.0 * M_PI * frequency;

    // Excitar el modo apropiado
    if (mode_type == "R") {
        // Modo derecho
        fields_[0][0] = amplitude * std::cos(omega * time); // Ex
        fields_[1][0] = amplitude * std::sin(omega * time); // Ey
    } else if (mode_type == "L") {
        // Modo izquierdo
        fields_[0][0] = amplitude * std::cos(omega * time); // Ex
        fields_[1][0] = -amplitude * std::sin(omega * time); // Ey
    } else if (mode_type == "O") {
        // Modo ordinario
        fields_[2][0] = amplitude * std::cos(omega * time); // Ez
    } else if (mode_type == "X") {
        // Modo extraordinario
        fields_[0][0] = amplitude * std::cos(omega * time); // Ex
        fields_[3][0] = amplitude * std::sin(omega * time); // Bx
    }
}


void TwoFluidSimulator::run_timesteps(int num_steps, double dt) {
    double current_time = 0.0;

    for (int step = 0; step < num_steps; ++step) {
        if (!current_mode_.empty()) {
            excite_mode(current_mode_, current_frequency_, current_amplitude_, current_time);
        }

        // Guardar estado actual para posibles reinicios
        fields_prev_ = fields_;
        currents_prev_ = currents_;

        try {
            // Usar RK4 para integración temporal
            update_system_rk4(dt);

            // Aplicar colisiones
            apply_collisions(dt);

            // Aplicar condiciones de contorno PML
            apply_boundary_conditions_pml();

        } catch (const std::exception& e) {
            // En caso de inestabilidad, restaurar estado anterior
            std::cerr << "Advertencia: Inestabilidad detectada. Restaurando estado anterior." << std::endl;
            fields_ = fields_prev_;
            currents_ = currents_prev_;
            break;
        }

        current_time += dt;
    }
}


void TwoFluidSimulator::compute_derivatives(
    const std::vector<std::vector<double>>& fields,
    const std::vector<std::vector<double>>& currents,
    std::vector<std::vector<double>>& dfields_dt,
    std::vector<std::vector<double>>& dcurrents_dt) {

    double c = params_.LIGHT_SPEED;
    double epsilon0 = params_.VACUUM_PERMITTIVITY;
    double dz = z_grid_[1] - z_grid_[0];
    double omega_pe = params_.electron_plasma_frequency();
    double omega_ce = params_.electron_cyclotron_frequency();
    double nu = params_.collision_frequency;

    // Inicializar derivadas a cero
    for (auto& vec : dfields_dt) std::fill(vec.begin(), vec.end(), 0.0);
    for (auto& vec : dcurrents_dt) std::fill(vec.begin(), vec.end(), 0.0);

    // Diferencias finitas de 4to orden para derivadas espaciales
    for (int i = 2; i < params_.grid_points - 2; ++i) {
        // Calcular todas las derivadas espaciales necesarias
        double dEx_dz = (-fields[0][i+2] + 8*fields[0][i+1] - 8*fields[0][i-1] + fields[0][i-2]) / (12.0 * dz);
        double dEy_dz = (-fields[1][i+2] + 8*fields[1][i+1] - 8*fields[1][i-1] + fields[1][i-2]) / (12.0 * dz);
        double dBx_dz = (-fields[3][i+2] + 8*fields[3][i+1] - 8*fields[3][i-1] + fields[3][i-2]) / (12.0 * dz);
        double dBy_dz = (-fields[4][i+2] + 8*fields[4][i+1] - 8*fields[4][i-1] + fields[4][i-2]) / (12.0 * dz);

        // Ecuaciones de Maxwell CORREGIDAS para 1D
        dfields_dt[0][i] = c * c * dBy_dz - currents[0][i] / epsilon0;  // dEx/dt = c²(∂By/∂z) - Jx/ε₀
        dfields_dt[1][i] = -c * c * dBx_dz - currents[1][i] / epsilon0; // dEy/dt = -c²(∂Bx/∂z) - Jy/ε₀
        dfields_dt[2][i] = -currents[2][i] / epsilon0;                  // dEz/dt = -Jz/ε₀

        dfields_dt[3][i] = -dEy_dz;  // dBx/dt = -∂Ey/∂z
        dfields_dt[4][i] = dEx_dz;   // dBy/dt = ∂Ex/∂z
        dfields_dt[5][i] = 0.0;      // dBz/dt = 0

        // Ecuaciones de momento para corrientes (respuesta del plasma)
        dcurrents_dt[0][i] = epsilon0 * omega_pe * omega_pe * fields[0][i] + 
            omega_ce * currents[1][i] - nu * currents[0][i];
        dcurrents_dt[1][i] = epsilon0 * omega_pe * omega_pe * fields[1][i] - 
            omega_ce * currents[0][i] - nu * currents[1][i];
        dcurrents_dt[2][i] = epsilon0 * omega_pe * omega_pe * fields[2][i] - 
            nu * currents[2][i];
    }
}

void TwoFluidSimulator::update_system_rk4(double dt) {
    // Vectores para almacenar las k (pendientes) de RK4
    std::vector<std::vector<double>> k1_fields(6, std::vector<double>(params_.grid_points));
    std::vector<std::vector<double>> k1_currents(3, std::vector<double>(params_.grid_points));
    std::vector<std::vector<double>> k2_fields(6, std::vector<double>(params_.grid_points));
    std::vector<std::vector<double>> k2_currents(3, std::vector<double>(params_.grid_points));
    std::vector<std::vector<double>> k3_fields(6, std::vector<double>(params_.grid_points));
    std::vector<std::vector<double>> k3_currents(3, std::vector<double>(params_.grid_points));
    std::vector<std::vector<double>> k4_fields(6, std::vector<double>(params_.grid_points));
    std::vector<std::vector<double>> k4_currents(3, std::vector<double>(params_.grid_points));

    // Vectores temporales para almacenar estados intermedios
    std::vector<std::vector<double>> fields_temp = fields_;
    std::vector<std::vector<double>> currents_temp = currents_;

    // --- Primer paso de RK4 (k1) ---
    compute_derivatives(fields_, currents_, k1_fields, k1_currents);

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

    compute_derivatives(fields_temp, currents_temp, k2_fields, k2_currents);

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

    compute_derivatives(fields_temp, currents_temp, k3_fields, k3_currents);

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

    compute_derivatives(fields_temp, currents_temp, k4_fields, k4_currents);

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
}

void TwoFluidSimulator::apply_collisions(double dt) {
    double nu = params_.collision_frequency;

    // Aplicar amortiguamiento exponencial debido a colisiones
    for (int i = 0; i < params_.grid_points; ++i) {
        currents_[0][i] *= std::exp(-nu * dt);
        currents_[1][i] *= std::exp(-nu * dt);
        currents_[2][i] *= std::exp(-nu * dt);
    }
}

void TwoFluidSimulator::apply_boundary_conditions_pml() {
    int pml_width = 20;  // Ancho de la capa PML
    double sigma_max = 0.01;  // Reducir la atenuación máxima

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
