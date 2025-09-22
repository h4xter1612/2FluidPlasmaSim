#ifndef PLASMA_PARAMS_HH
#define PLASMA_PARAMS_HH

struct PlasmaParams {
    double electron_density;    // n_e (m^{-3})
    double magnetic_field;      // B_0 (Tesla)
    double collision_frequency; // ν (s^{-1})
    double length;              // Longitud del dominio (m)
    int grid_points;            // Puntos de la malla
    double ion_density;                   // n_i (m^{-3})
    double ion_temperature;               // T_i (eV) - para extensiones futuras
    double ion_collision_frequency;       // ν_i (s^{-1})
    double ion_mass;                      // m_i (kg) - typically proton mass

    // 2D PlasmaParams
    int dimension;           // 1 o 2
    int nx, ny;              // Puntos de la malla en x e y
    double dx, dy;           // Espaciado en x e y
    double length_x, length_y; // Dimensiones del dominio
    
    // Constantes físicas
    static constexpr double PROTON_MASS = 1.6726e-27;
    static constexpr double ELECTRON_MASS = 9.109e-31;
    static constexpr double ELEMENTARY_CHARGE = 1.602e-19;
    static constexpr double LIGHT_SPEED = 2.998e8;
    static constexpr double VACUUM_PERMITTIVITY = 8.854e-12;
    
    // Métodos para calcular frecuencias características
    double electron_plasma_frequency() const;
    double electron_cyclotron_frequency() const;
    // Métodos para frecuencias iónicas
    double ion_plasma_frequency() const;
    double ion_cyclotron_frequency() const;
    // Métodos auxiliares
    bool is_1d() const { return dimension == 1; }
    bool is_2d() const { return dimension == 2; }
};

#endif
