#include "simulator_2d.hh"
#include "dispersion.hh"
#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

Simulator2D::Simulator2D(const PlasmaParams& params)
    : SimulatorBase(params) {}

void Simulator2D::initialize() {
    // Inicializar malla 2D
    x_grid_.resize(params_.nx);
    y_grid_.resize(params_.ny);
    
    double dx = params_.length_x / (params_.nx - 1);
    double dy = params_.length_y / (params_.ny - 1);
    
    for (int i = 0; i < params_.nx; ++i) {
        x_grid_[i] = i * dx;
    }
    for (int j = 0; j < params_.ny; ++j) {
        y_grid_[j] = j * dy;
    }

    // Inicializar campos y corrientes (6 campos + 3 corrientes electrones + 3 corrientes iones)
    fields_.resize(6, std::vector<std::vector<double>>(params_.nx, std::vector<double>(params_.ny, 0.0)));
    currents_.resize(3, std::vector<std::vector<double>>(params_.nx, std::vector<double>(params_.ny, 0.0)));
    currents_ion_.resize(3, std::vector<std::vector<double>>(params_.nx, std::vector<double>(params_.ny, 0.0)));
}

void Simulator2D::excite_mode(const std::string& mode_type, double frequency, double amplitude, double time) {
    double omega = 2.0 * M_PI * frequency;
    double c = params_.LIGHT_SPEED;

    // Excitar en el borde izquierdo (x=0) para todo y
    for (int j = 0; j < params_.ny; ++j) {
        if (mode_type == "R") {
            fields_[0][0][j] = amplitude * std::cos(omega * time);
            fields_[1][0][j] = amplitude * std::sin(omega * time);
            fields_[3][0][j] = -amplitude * std::sin(omega * time) / c;
            fields_[4][0][j] = amplitude * std::cos(omega * time) / c;
        } else if (mode_type == "L") {
            fields_[0][0][j] = amplitude * std::cos(omega * time);
            fields_[1][0][j] = -amplitude * std::sin(omega * time);
            fields_[3][0][j] = amplitude * std::sin(omega * time) / c;
            fields_[4][0][j] = amplitude * std::cos(omega * time) / c;
        } else if (mode_type == "O") {
            fields_[2][0][j] = amplitude * std::cos(omega * time);
        } else if (mode_type == "X") {
            fields_[0][0][j] = amplitude * std::cos(omega * time);
            fields_[3][0][j] = amplitude * std::sin(omega * time);
        }
    }
}

void Simulator2D::run_timesteps(int num_steps, double dt) {
    double current_time = 0.0;
    double c = params_.LIGHT_SPEED;
    double dx = x_grid_[1] - x_grid_[0];
    double dy = y_grid_[1] - y_grid_[0];
    
    // Verificar condición CFL 2D
    double cfl_dt = 1.0 / (c * std::sqrt(1.0/(dx*dx) + 1.0/(dy*dy)));
    if (dt > cfl_dt * 0.1) {
        std::cout << "Advertencia: dt = " << dt << " puede ser demasiado grande. CFL recomienda dt < " << cfl_dt * 0.1 << std::endl;
    }
    
    for (int step = 0; step < num_steps; ++step) {
        if (!current_mode_.empty()) {
            excite_mode(current_mode_, current_frequency_, current_amplitude_, current_time);
        }

        if (step_count_ % save_interval_ == 0) {
            std::string filename = "data/snap/field_data_2d_" + current_mode_ + "_" + 
                                  std::to_string(step_count_ / save_interval_) + ".csv";
            export_field_data(filename);
        }
        
        update_system_rk4(dt);
        apply_collisions(dt);
        apply_boundary_conditions_pml();
        
        current_time += dt;
        step_count_++;
    }
    
    std::string filename = "data/field_data_2d_" + current_mode_ + "_final.csv";
    export_field_data(filename);
}

void Simulator2D::compute_derivatives(
    const std::vector<std::vector<std::vector<double>>>& fields,
    const std::vector<std::vector<std::vector<double>>>& currents,
    const std::vector<std::vector<std::vector<double>>>& currents_ion,
    std::vector<std::vector<std::vector<double>>>& dfields_dt,
    std::vector<std::vector<std::vector<double>>>& dcurrents_dt,
    std::vector<std::vector<std::vector<double>>>& dcurrents_dt_ion) {

    double c = params_.LIGHT_SPEED;
    double epsilon0 = params_.VACUUM_PERMITTIVITY;
    double dx = x_grid_[1] - x_grid_[0];
    double dy = y_grid_[1] - y_grid_[0];
    
    double omega_pe = params_.electron_plasma_frequency();
    double omega_ce = params_.electron_cyclotron_frequency();
    double nu = params_.collision_frequency;
    double omega_pi = params_.ion_plasma_frequency();
    double omega_ci = params_.ion_cyclotron_frequency();
    double nu_i = params_.ion_collision_frequency;

    // Inicializar derivadas a cero
    for (auto& comp : dfields_dt)
        for (auto& row : comp)
            std::fill(row.begin(), row.end(), 0.0);
            
    for (auto& comp : dcurrents_dt)
        for (auto& row : comp)
            std::fill(row.begin(), row.end(), 0.0);
            
    for (auto& comp : dcurrents_dt_ion)
        for (auto& row : comp)
            std::fill(row.begin(), row.end(), 0.0);

    // Diferencias finitas de 4to orden
    for (int i = 2; i < params_.nx - 2; ++i) {
        for (int j = 2; j < params_.ny - 2; ++j) {
            // Derivadas de E
            double dEx_dy = (-fields[0][i][j+2] + 8*fields[0][i][j+1] - 8*fields[0][i][j-1] + fields[0][i][j-2]) / (12.0 * dy);
            double dEy_dx = (-fields[1][i+2][j] + 8*fields[1][i+1][j] - 8*fields[1][i-1][j] + fields[1][i-2][j]) / (12.0 * dx);
            double dEz_dx = (-fields[2][i+2][j] + 8*fields[2][i+1][j] - 8*fields[2][i-1][j] + fields[2][i-2][j]) / (12.0 * dx);
            double dEz_dy = (-fields[2][i][j+2] + 8*fields[2][i][j+1] - 8*fields[2][i][j-1] + fields[2][i][j-2]) / (12.0 * dy);

            // Derivadas de B
            double dBx_dy = (-fields[3][i][j+2] + 8*fields[3][i][j+1] - 8*fields[3][i][j-1] + fields[3][i][j-2]) / (12.0 * dy);
            double dBy_dx = (-fields[4][i+2][j] + 8*fields[4][i+1][j] - 8*fields[4][i-1][j] + fields[4][i-2][j]) / (12.0 * dx);
            double dBz_dx = (-fields[5][i+2][j] + 8*fields[5][i+1][j] - 8*fields[5][i-1][j] + fields[5][i-2][j]) / (12.0 * dx);
            double dBz_dy = (-fields[5][i][j+2] + 8*fields[5][i][j+1] - 8*fields[5][i][j-1] + fields[5][i][j-2]) / (12.0 * dy);

            // Ley de Faraday para B
            dfields_dt[3][i][j] = -dEz_dy;  // dBx/dt
            dfields_dt[4][i][j] = dEz_dx;   // dBy/dt
            dfields_dt[5][i][j] = dEx_dy - dEy_dx;  // dBz/dt

            // Ley de Ampère-Maxwell para E
            dfields_dt[0][i][j] = c*c * dBz_dy - (currents[0][i][j] + currents_ion[0][i][j])/epsilon0;  // dEx/dt
            dfields_dt[1][i][j] = -c*c * dBz_dx - (currents[1][i][j] + currents_ion[1][i][j])/epsilon0; // dEy/dt
            dfields_dt[2][i][j] = c*c * (dBy_dx - dBx_dy) - (currents[2][i][j] + currents_ion[2][i][j])/epsilon0; // dEz/dt

            // Ecuaciones de momento para electrones
            dcurrents_dt[0][i][j] = epsilon0*omega_pe*omega_pe*fields[0][i][j] + 
                                   omega_ce*currents[1][i][j] - nu*currents[0][i][j];
            dcurrents_dt[1][i][j] = epsilon0*omega_pe*omega_pe*fields[1][i][j] - 
                                   omega_ce*currents[0][i][j] - nu*currents[1][i][j];
            dcurrents_dt[2][i][j] = epsilon0*omega_pe*omega_pe*fields[2][i][j] - 
                                   nu*currents[2][i][j];

            // Ecuaciones de momento para iones
            dcurrents_dt_ion[0][i][j] = epsilon0*omega_pi*omega_pi*fields[0][i][j] - 
                                      omega_ci*currents_ion[1][i][j] - nu_i*currents_ion[0][i][j];
            dcurrents_dt_ion[1][i][j] = epsilon0*omega_pi*omega_pi*fields[1][i][j] + 
                                      omega_ci*currents_ion[0][i][j] - nu_i*currents_ion[1][i][j];
            dcurrents_dt_ion[2][i][j] = epsilon0*omega_pi*omega_pi*fields[2][i][j] - 
                                      nu_i*currents_ion[2][i][j];
        }
    }
}

void Simulator2D::update_system_rk4(double dt) {
    // Vectores para almacenar las k (pendientes) de RK4
    std::vector<std::vector<std::vector<double>>> k1_fields(6, 
        std::vector<std::vector<double>>(params_.nx, std::vector<double>(params_.ny)));
    std::vector<std::vector<std::vector<double>>> k2_fields(6, 
        std::vector<std::vector<double>>(params_.nx, std::vector<double>(params_.ny)));
    std::vector<std::vector<std::vector<double>>> k3_fields(6, 
        std::vector<std::vector<double>>(params_.nx, std::vector<double>(params_.ny)));
    std::vector<std::vector<std::vector<double>>> k4_fields(6, 
        std::vector<std::vector<double>>(params_.nx, std::vector<double>(params_.ny)));
    
    std::vector<std::vector<std::vector<double>>> k1_currents(3, 
        std::vector<std::vector<double>>(params_.nx, std::vector<double>(params_.ny)));
    std::vector<std::vector<std::vector<double>>> k2_currents(3, 
        std::vector<std::vector<double>>(params_.nx, std::vector<double>(params_.ny)));
    std::vector<std::vector<std::vector<double>>> k3_currents(3, 
        std::vector<std::vector<double>>(params_.nx, std::vector<double>(params_.ny)));
    std::vector<std::vector<std::vector<double>>> k4_currents(3, 
        std::vector<std::vector<double>>(params_.nx, std::vector<double>(params_.ny)));
    
    std::vector<std::vector<std::vector<double>>> k1_currents_ion(3, 
        std::vector<std::vector<double>>(params_.nx, std::vector<double>(params_.ny)));
    std::vector<std::vector<std::vector<double>>> k2_currents_ion(3, 
        std::vector<std::vector<double>>(params_.nx, std::vector<double>(params_.ny)));
    std::vector<std::vector<std::vector<double>>> k3_currents_ion(3, 
        std::vector<std::vector<double>>(params_.nx, std::vector<double>(params_.ny)));
    std::vector<std::vector<std::vector<double>>> k4_currents_ion(3, 
        std::vector<std::vector<double>>(params_.nx, std::vector<double>(params_.ny)));
    
    // Vectores temporales para almacenar estados intermedios
    auto fields_temp = fields_;
    auto currents_temp = currents_;
    auto currents_ion_temp = currents_ion_;

    // --- Primer paso de RK4 (k1) ---
    compute_derivatives(fields_, currents_, currents_ion_, k1_fields, k1_currents, k1_currents_ion);

    // --- Segundo paso de RK4 (k2) ---
    for (int comp = 0; comp < 6; ++comp) {
        for (int i = 0; i < params_.nx; ++i) {
            for (int j = 0; j < params_.ny; ++j) {
                fields_temp[comp][i][j] = fields_[comp][i][j] + 0.5 * dt * k1_fields[comp][i][j];
            }
        }
    }
    for (int comp = 0; comp < 3; ++comp) {
        for (int i = 0; i < params_.nx; ++i) {
            for (int j = 0; j < params_.ny; ++j) {
                currents_temp[comp][i][j] = currents_[comp][i][j] + 0.5 * dt * k1_currents[comp][i][j];
            }
        }
    }
    for (int comp = 0; comp < 3; ++comp) {
        for (int i = 0; i < params_.nx; ++i) {
            for (int j = 0; j < params_.ny; ++j) {
                currents_ion_temp[comp][i][j] = currents_ion_[comp][i][j] + 0.5 * dt * k1_currents_ion[comp][i][j];
            }
        }
    }
    compute_derivatives(fields_temp, currents_temp, currents_ion_temp, k2_fields, k2_currents, k2_currents_ion);

    // --- Tercer paso de RK4 (k3) ---
    for (int comp = 0; comp < 6; ++comp) {
        for (int i = 0; i < params_.nx; ++i) {
            for (int j = 0; j < params_.ny; ++j) {
                fields_temp[comp][i][j] = fields_[comp][i][j] + 0.5 * dt * k2_fields[comp][i][j];
            }
        }
    }
    for (int comp = 0; comp < 3; ++comp) {
        for (int i = 0; i < params_.nx; ++i) {
            for (int j = 0; j < params_.ny; ++j) {
                currents_temp[comp][i][j] = currents_[comp][i][j] + 0.5 * dt * k2_currents[comp][i][j];
            }
        }
    }
    for (int comp = 0; comp < 3; ++comp) {
        for (int i = 0; i < params_.nx; ++i) {
            for (int j = 0; j < params_.ny; ++j) {
                currents_ion_temp[comp][i][j] = currents_ion_[comp][i][j] + 0.5 * dt * k2_currents_ion[comp][i][j];
            }
        }
    }
    compute_derivatives(fields_temp, currents_temp, currents_ion_temp, k3_fields, k3_currents, k3_currents_ion);

    // --- Cuarto paso de RK4 (k4) ---
    for (int comp = 0; comp < 6; ++comp) {
        for (int i = 0; i < params_.nx; ++i) {
            for (int j = 0; j < params_.ny; ++j) {
                fields_temp[comp][i][j] = fields_[comp][i][j] + dt * k3_fields[comp][i][j];
            }
        }
    }
    for (int comp = 0; comp < 3; ++comp) {
        for (int i = 0; i < params_.nx; ++i) {
            for (int j = 0; j < params_.ny; ++j) {
                currents_temp[comp][i][j] = currents_[comp][i][j] + dt * k3_currents[comp][i][j];
            }
        }
    }
    for (int comp = 0; comp < 3; ++comp) {
        for (int i = 0; i < params_.nx; ++i) {
            for (int j = 0; j < params_.ny; ++j) {
                currents_ion_temp[comp][i][j] = currents_ion_[comp][i][j] + dt * k3_currents_ion[comp][i][j];
            }
        }
    }
    compute_derivatives(fields_temp, currents_temp, currents_ion_temp, k4_fields, k4_currents, k4_currents_ion);

    // --- Combinar resultados para campos ---
    for (int comp = 0; comp < 6; ++comp) {
        for (int i = 0; i < params_.nx; ++i) {
            for (int j = 0; j < params_.ny; ++j) {
                fields_[comp][i][j] += dt * (k1_fields[comp][i][j] + 2*k2_fields[comp][i][j] + 
                    2*k3_fields[comp][i][j] + k4_fields[comp][i][j]) / 6.0;
            }
        }
    }

    // --- Combinar resultados para corrientes ---
    for (int comp = 0; comp < 3; ++comp) {
        for (int i = 0; i < params_.nx; ++i) {
            for (int j = 0; j < params_.ny; ++j) {
                currents_[comp][i][j] += dt * (k1_currents[comp][i][j] + 2*k2_currents[comp][i][j] + 
                    2*k3_currents[comp][i][j] + k4_currents[comp][i][j]) / 6.0;
            }
        }
    }
    
    for (int comp = 0; comp < 3; ++comp) {
        for (int i = 0; i < params_.nx; ++i) {
            for (int j = 0; j < params_.ny; ++j) {
                currents_ion_[comp][i][j] += dt * (k1_currents_ion[comp][i][j] + 2*k2_currents_ion[comp][i][j] + 
                    2*k3_currents_ion[comp][i][j] + k4_currents_ion[comp][i][j]) / 6.0;
            }
        }
    }
    
    // Verificar valores problemáticos
    for (int comp = 0; comp < 6; ++comp) {
        for (int i = 0; i < params_.nx; ++i) {
            for (int j = 0; j < params_.ny; ++j) {
                if (std::isnan(fields_[comp][i][j]) || std::isinf(fields_[comp][i][j])) {
                    fields_[comp][i][j] = 0.0;
                }
            }
        }
    }
    
    for (int comp = 0; comp < 3; ++comp) {
        for (int i = 0; i < params_.nx; ++i) {
            for (int j = 0; j < params_.ny; ++j) {
                if (std::isnan(currents_[comp][i][j]) || std::isinf(currents_[comp][i][j])) {
                    currents_[comp][i][j] = 0.0;
                }
            }
        }
    }
    for (int comp = 0; comp < 3; ++comp) {
        for (int i = 0; i < params_.nx; ++i) {
            for (int j = 0; j < params_.ny; ++j) {
                if (std::isnan(currents_ion_[comp][i][j]) || std::isinf(currents_ion_[comp][i][j])) {
                    currents_ion_[comp][i][j] = 0.0;
                }
            }
        }
    }
}

void Simulator2D::apply_collisions(double dt) {
    double nu = params_.collision_frequency;
    double nu_i = params_.ion_collision_frequency;

    for (int i = 0; i < params_.nx; ++i) {
        for (int j = 0; j < params_.ny; ++j) {
            for (int comp = 0; comp < 3; ++comp) {
                currents_[comp][i][j] *= std::exp(-nu * dt);
                currents_ion_[comp][i][j] *= std::exp(-nu_i * dt);
            }
        }
    }
}

void Simulator2D::apply_boundary_conditions_pml() {
    int pml_width = 10;
    double sigma_max = 0.001;

    // Capas PML en los bordes x
    for (int i = 0; i < pml_width; ++i) {
        double sigma = sigma_max * std::pow((double)(pml_width - i) / pml_width, 4);
        for (int j = 0; j < params_.ny; ++j) {
            for (int comp = 0; comp < 6; ++comp) {
                fields_[comp][i][j] *= std::exp(-sigma);
                fields_[comp][params_.nx-1-i][j] *= std::exp(-sigma);
            }
        }
    }

    // Capas PML en los bordes y
    for (int j = 0; j < pml_width; ++j) {
        double sigma = sigma_max * std::pow((double)(pml_width - j) / pml_width, 4);
        for (int i = 0; i < params_.nx; ++i) {
            for (int comp = 0; comp < 6; ++comp) {
                fields_[comp][i][j] *= std::exp(-sigma);
                fields_[comp][i][params_.ny-1-j] *= std::exp(-sigma);
            }
        }
    }
}

void Simulator2D::export_field_data(const std::string& filename) const {
    std::ofstream file(filename);
    file << "x,y,Ex,Ey,Ez,Bx,By,Bz,Jx_e,Jy_e,Jz_e,Jx_i,Jy_i,Jz_i\n";

    for (int i = 0; i < params_.nx; ++i) {
        for (int j = 0; j < params_.ny; ++j) {
            file << x_grid_[i] << "," << y_grid_[j];
            for (int comp = 0; comp < 6; ++comp) file << "," << fields_[comp][i][j];
            for (int comp = 0; comp < 3; ++comp) file << "," << currents_[comp][i][j];
            for (int comp = 0; comp < 3; ++comp) file << "," << currents_ion_[comp][i][j];
            file << "\n";
        }
    }
    file.close();
}

void Simulator2D::export_dispersion_data(const std::string& filename) const {
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
