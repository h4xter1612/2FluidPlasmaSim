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
    std::cout << "Datos de dispersion exportados a dispersion_data.csv\n";

    // Calcular el paso de tiempo máximo permitido por CFL
    double dz = params.length / (params.grid_points - 1);
    double dt = 0.1 * dz / params.LIGHT_SPEED;  // 10% del límite CFL
    
    std::cout << "Usando paso de tiempo dt = " << dt << " s (CFL recomienda < " << dz/params.LIGHT_SPEED << " s)\n";

    simulator.set_save_interval(100);  // Guardar cada 100 pasos

    // Simular modo R
    std::cout << "Simulando modo R...\n";
    simulator.set_mode("R");
    simulator.set_frequency(5e9);  // 5 GHz
    simulator.set_amplitude(1.0);
    simulator.run_timesteps(10000, dt);
    simulator.export_field_data(dataf+"field_data_R.csv");
    std::cout << "Datos del modo R exportados a field_data_R.csv\n";

    // Reinicializar para la próxima simulación
    simulator.initialize();
    
    // Simular modo L
    std::cout << "Simulando modo L...\n";
    simulator.set_mode("L");
    simulator.set_frequency(5e9);
    simulator.set_amplitude(1.0);
    simulator.run_timesteps(10000, dt);
    simulator.export_field_data(dataf+"field_data_L.csv");
    std::cout << "Datos del modo L exportados a field_data_L.csv\n";

    // Reinicializar para la próxima simulación
    simulator.initialize();
    
    // Simular modo O
    std::cout << "Simulando modo O...\n";
    simulator.set_mode("O");
    simulator.set_frequency(5e9);
    simulator.set_amplitude(1.0);
    simulator.run_timesteps(10000, dt);
    simulator.export_field_data(dataf+"field_data_O.csv");
    std::cout << "Datos del modo O exportados a field_data_O.csv\n";

    // Reinicializar para la próxima simulación
    simulator.initialize();
    
    // Simular modo X
    std::cout << "Simulando modo X...\n";
    simulator.set_mode("X");
    simulator.set_frequency(5e9);
    simulator.set_amplitude(1.0);
    simulator.run_timesteps(10000, dt);
    simulator.export_field_data(dataf+"field_data_X.csv");
    std::cout << "Datos del modo X exportados to field_data_X.csv\n";

    std::cout << "Todas las simulaciones completadas. Use Python para visualizar los resultados.\n";
    
    return 0;
}
