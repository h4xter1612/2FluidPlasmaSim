#include "simulator_1d.hh"
#include "simulator_2d.hh"
#include <iostream>
#include <cmath>
#include <string>
#include <memory>
#include <algorithm>

int main(int argc, char* argv[]) {
    // Determinar dimensionalidad desde argumentos
    int dimension = 1;
    if (argc > 1) {
        dimension = std::stoi(argv[1]);
    }
    
    // Configurar parámetros del plasma
    PlasmaParams params;
    params.dimension = dimension;
    std::string dataf;
    
    if (dimension == 1) {
        dataf = "data/field_data";
        params.nx = 1000;
        params.ny = 1;
        params.length_x = 1.0;
        params.length_y = 0.0;
        params.dx = params.length_x / (params.nx - 1);
        params.dy = 0.0;
        // Mantener compatibilidad con código 1D existente
        params.grid_points = params.nx;
        params.length = params.length_x;
    } else {
        dataf = "data/field_data_2d";
        // Para 2D, usar una malla más pequeña por razones de rendimiento
        params.nx = 200;
        params.ny = 200;
        params.length_x = 1.0;
        params.length_y = 1.0;
        params.dx = params.length_x / (params.nx - 1);
        params.dy = params.length_y / (params.ny - 1);
        // Para compatibilidad con código 1D existente
        params.grid_points = params.nx;
        params.length = params.length_x;
    }

    params.electron_density = 1e18;
    params.magnetic_field = 0.1;
    params.collision_frequency = 1e7;
    params.ion_density = 1e18;
    params.ion_mass = params.PROTON_MASS;
    params.ion_collision_frequency = 1e6;

    // Crear el simulador apropiado
    std::unique_ptr<SimulatorBase> simulator;
    
    if (dimension == 1) {
        simulator = std::make_unique<Simulator1D>(params);
        std::cout << "Running 1D simulation\n";
    } else {
        simulator = std::make_unique<Simulator2D>(params);
        std::cout << "Running 2D simulation\n";
    }
    
    simulator->initialize();


    // Calcular y exportar relaciones de dispersión
    simulator->export_dispersion_data("data/dispersion_data.csv");
    std::cout << "Dispersion data exported to dispersion_data.csv\n";

    // Calcular el paso de tiempo máximo permitido por CFL
    double dt;
    if (dimension == 1) {
        double dz = params.length_x / (params.nx - 1);
        dt = 0.01 * dz / params.LIGHT_SPEED;
        std::cout << "Using timestep dt = " << dt << " s (CFL recommends < " << dz/params.LIGHT_SPEED << " s)\n";
    } else {
        // Condición CFL para 2D
        double dx = params.length_x / (params.nx - 1);
        double dy = params.length_y / (params.ny - 1);
        double cfl_dt = 1.0 / (params.LIGHT_SPEED * std::sqrt(1.0/(dx*dx) + 1.0/(dy*dy)));
        dt = 0.01 * cfl_dt;
        std::cout << "Using timestep dt = " << dt << " s (CFL recommends < " << cfl_dt << " s)\n";
    }

    simulator->set_save_interval(100);

    // Para 2D, reducir el número de pasos de tiempo
    int num_steps = (dimension == 1) ? 10000 : 1000;

    // Simular modo R
    std::cout << "Simulating Mode R...\n";
    simulator->set_mode("R");
    simulator->set_frequency(5e9);
    simulator->set_amplitude(1.0);
    simulator->run_timesteps(num_steps, dt);
    simulator->export_field_data(dataf+"_R.csv");
    std::cout << "Data of mode R exported.\n";

    // Reinicializar para la próxima simulación
    simulator->initialize();
    
    // Simular modo L
    std::cout << "Simulating Mode L...\n";
    simulator->set_mode("L");
    simulator->set_frequency(5e9);
    simulator->set_amplitude(1.0);
    simulator->run_timesteps(num_steps, dt);
    simulator->export_field_data(dataf+"_L.csv");
    std::cout << "Data of mode L exported.\n";

    // Simular modo O
    std::cout << "Simulating Mode O...\n";
    simulator->set_mode("O");
    simulator->set_frequency(5e9);
    simulator->set_amplitude(1.0);
    simulator->run_timesteps(num_steps, dt);
    simulator->export_field_data(dataf+"_O.csv");
    std::cout << "Data of mode O exported.\n";

    // Simular modo X
    std::cout << "Simulating Mode X...\n";
    simulator->set_mode("X");
    simulator->set_frequency(5e9);
    simulator->set_amplitude(1.0);
    simulator->run_timesteps(num_steps, dt);
    simulator->export_field_data(dataf+"_X.csv");
    std::cout << "Data of mode X exported.\n";

    std::cout << "All simulations completed. Use Python for results visualization.\n";
    
    return 0;
}
