#ifndef SIMULATOR_HH
#define SIMULATOR_HH

#include "plasma_params.hh"
#include "dispersion.hh"
#include <vector>
#include <string>

class TwoFluidSimulator {
public:
    TwoFluidSimulator(const PlasmaParams& params);
    
    void initialize();
    void run_timesteps(int num_steps, double dt);
    void excite_mode(const std::string& mode_type, double frequency, double amplitude, double time);
    void export_field_data(const std::string& filename) const;
    void export_dispersion_data(const std::string& filename) const;
    
private:
    PlasmaParams params_;
    DispersionRelation dispersion_;
    
    std::vector<double> z_grid_;
    std::vector<std::vector<double>> fields_; // Ex, Ey, Ez, Bx, By, Bz
    std::vector<std::vector<double>> currents_; // Jx, Jy, Jz
    std::vector<std::vector<double>> fields_prev_; // Campos en paso anterior (para RK4)
    std::vector<std::vector<double>> currents_prev_; // Corrientes en paso anterior (para RK4)
    
    std::string current_mode_;
    double current_frequency_;
    double current_amplitude_;
    
    // Métodos internos para RK4
    void update_system_rk4(double dt);  // Update fields and currents
    void apply_boundary_conditions_pml();
    void apply_collisions(double dt);
    
    // Helper functions para RK4
    void compute_derivatives(
        const std::vector<std::vector<double>>& fields,
        const std::vector<std::vector<double>>& currents,
        std::vector<std::vector<double>>& dfields_dt,
        std::vector<std::vector<double>>& dcurrents_dt);
};

#endif
