#ifndef SIMULATOR_1D_HH
#define SIMULATOR_1D_HH

#include "simulator_base.hh"
#include <vector>

class Simulator1D : public SimulatorBase {
public:
    Simulator1D(const PlasmaParams& params);
    
    void initialize() override;
    void run_timesteps(int num_steps, double dt) override;
    void excite_mode(const std::string& mode_type, double frequency, double amplitude, double time) override;
    void export_field_data(const std::string& filename) const override;
    void export_dispersion_data(const std::string& filename) const override;
    void export_field_data_binary(const std::string& filename) const;

    std::vector<double> calculate_energy_densities() const;

private:
    std::vector<double> z_grid_;
    std::vector<std::vector<double>> fields_;
    std::vector<std::vector<double>> currents_;
    std::vector<std::vector<double>> currents_ion_;
    std::vector<std::vector<double>> fields_prev_;
    std::vector<std::vector<double>> currents_prev_;
    std::vector<std::vector<double>> currents_ion_prev_;
    
    void update_system_rk4(double dt);
    void apply_collisions(double dt);
    void apply_boundary_conditions_pml();
    
    void compute_derivatives(
        const std::vector<std::vector<double>>& fields,
        const std::vector<std::vector<double>>& currents,
        const std::vector<std::vector<double>>& currents_ion,
        std::vector<std::vector<double>>& dfields_dt,
        std::vector<std::vector<double>>& dcurrents_dt,
        std::vector<std::vector<double>>& dcurrents_dt_ion);
};

#endif
