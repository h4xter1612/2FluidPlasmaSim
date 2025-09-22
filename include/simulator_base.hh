#ifndef SIMULATOR_BASE_HH
#define SIMULATOR_BASE_HH

#include "plasma_params.hh"
#include "dispersion.hh"
#include <vector>
#include <string>
#include <memory>

class SimulatorBase {
public:
    SimulatorBase(const PlasmaParams& params);
    virtual ~SimulatorBase() = default;
    
    // Métodos virtuales puros
    virtual void initialize() = 0;
    virtual void run_timesteps(int num_steps, double dt) = 0;
    virtual void excite_mode(const std::string& mode_type, double frequency, double amplitude, double time) = 0;
    virtual void export_field_data(const std::string& filename) const = 0;
    virtual void export_dispersion_data(const std::string& filename) const = 0;
    
    // Métodos con implementación por defecto
    virtual void set_mode(const std::string& mode) { current_mode_ = mode; }
    virtual void set_frequency(double freq) { current_frequency_ = freq; }
    virtual void set_amplitude(double amp) { current_amplitude_ = amp; }
    virtual void set_save_interval(int interval) { save_interval_ = interval; }

protected:
    PlasmaParams params_;
    DispersionRelation dispersion_;
    
    std::string current_mode_;
    double current_frequency_;
    double current_amplitude_;
    int save_interval_ = 100;
    int step_count_ = 0;
};

#endif
