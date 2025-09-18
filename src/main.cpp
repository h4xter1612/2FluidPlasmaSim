#include "simulator.hh"
#include <iostream>

int main() {
    // Configurar parámetros del plasma
    PlasmaParams params;
    params.electron_density = 1e18;      // m^{-3}
    params.magnetic_field = 0.1;         // Tesla
    params.collision_frequency = 1e7;    // s^{-1}
    params.length = 1.0;                 // m
    params.grid_points = 1000;
    
    // Crear y configurar simulador
    TwoFluidSimulator simulator(params);
    simulator.initialize();
    
    // Data folder
    std::string data = "data/";
    
    // Calcular y exportar relaciones de dispersión
    simulator.export_dispersion_data(data+"dispersion_data.csv");
   
    // Simular diferentes modos
    std::cout << "Simulating mode R...\n";
    simulator.excite_mode("R", 5e9, 1.0);
    simulator.run_timesteps(1000, 1e-12);
    simulator.export_field_data(data+"field_data_R.csv");
    
    std::cout << "Simulating mode L...\n";
    simulator.excite_mode("L", 5e9, 1.0);
    simulator.run_timesteps(1000, 1e-12);
    simulator.export_field_data(data+"field_data_L.csv");
    
    std::cout << "Simulating mode O...\n";
    simulator.excite_mode("O", 5e9, 1.0);
    simulator.run_timesteps(1000, 1e-12);
    simulator.export_field_data(data+"field_data_O.csv");
    
    std::cout << "Simulating mode X...\n";
    simulator.excite_mode("X", 5e9, 1.0);
    simulator.run_timesteps(1000, 1e-12);
    simulator.export_field_data(data+"field_data_X.csv");
    
    std::cout << "Simulaction completed. Use Python for results visualization.\n";
    return 0;
}
