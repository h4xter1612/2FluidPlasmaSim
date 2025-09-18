#ifndef SIMULATOR_HH
#define SIMULATOR_HH

#include "plasma_params.hh"
#include "dispersion.hh"
#include <vector>
#include <string>

class TwoFluidSimulator {
public:
    TwoFluidSimulator(const PlasmaParams& params);
    
    // Inicializar la simulación
    void initialize();
    
    // Ejecutar pasos de tiempo
    void run_timesteps(int num_steps, double dt);
    
    // Excitar un modo específico
    void excite_mode(const std::string& mode_type, double frequency, double amplitude);
    
    // Exportar datos para visualización
    void export_field_data(const std::string& filename) const;
    void export_dispersion_data(const std::string& filename) const;
    
private:
    PlasmaParams params_;
    DispersionRelation dispersion_;
    
    // Campos y variables de estado
    std::vector<double> z_grid_;
    std::vector<std::vector<double>> fields_; // Ex, Ey, Ez, Bx, By, Bz
    std::vector<std::vector<double>> currents_; // Jx, Jy, Jz
    
    // Parámetros de la excitación
    std::string current_mode_;
    double current_frequency_;
    double current_amplitude_;
    
    // Métodos internos
    void update_fields(double dt);
    void update_currents(double dt);
    void apply_boundary_conditions();
    void apply_collisions(double dt);
};

#endif
