#include "simulator.hh"
#include <iostream>
#include <cmath>
#include <string>

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

    std::string dataf = "data/";

    // Calcular y exportar relaciones de dispersión
    simulator.export_dispersion_data(dataf+"dispersion_data.csv");
    std::cout << "Dispersion data exported to dispersion_data.csv\n";

    // Calcular el paso de tiempo máximo permitido por CFL
    double dz = params.length / (params.grid_points - 1);
    double dt = 0.1 * dz / params.LIGHT_SPEED;  // 10% del límite CFL
    
    // Courant-Friedrichs-Levy Condition
    std::cout << "Using timestep dt = " << dt << " s (CFL recommends < " << dz/params.LIGHT_SPEED << " s)\n";

    simulator.set_save_interval(100);  // Guardar cada 100 pasos

    // Simular modo R
    std::cout << "Simulating Mode R...\n";
    simulator.set_mode("R");
    simulator.set_frequency(5e9);  // 5 GHz
    simulator.set_amplitude(1.0);
    simulator.run_timesteps(10000, dt);
    simulator.export_field_data(dataf+"field_data_R.csv");
    std::cout << "Data of mode R exported.\n";

    // Reinicializar para la próxima simulación
    simulator.initialize();
    
    // Simular modo L
    std::cout << "Simulating Mode L...\n";
    simulator.set_mode("L");
    simulator.set_frequency(5e9);
    simulator.set_amplitude(1.0);
    simulator.run_timesteps(10000, dt);
    simulator.export_field_data(dataf+"field_data_L.csv");
    std::cout << "Data of mode L exported.\n";

    // Reinicializar para la próxima simulación
    simulator.initialize();
    
    // Simular modo O
    std::cout << "Simulating Mode O...\n";
    simulator.set_mode("O");
    simulator.set_frequency(5e9);
    simulator.set_amplitude(1.0);
    simulator.run_timesteps(10000, dt);
    simulator.export_field_data(dataf+"field_data_O.csv");
    std::cout << "Data of mode O exported.\n";

    // Reinicializar para la próxima simulación
    simulator.initialize();
    
    // Simular modo X
    std::cout << "Simulating Mode X...\n";
    simulator.set_mode("X");
    simulator.set_frequency(5e9);
    simulator.set_amplitude(1.0);
    simulator.run_timesteps(10000, dt);
    simulator.export_field_data(dataf+"field_data_X.csv");
    std::cout << "Data of mode X exported.\n";

    std::cout << "All simulations completed. Use Python for results visualization.\n";
    
    return 0;
}
