import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd
import glob
import os
import re

def numeric_sort_key(path):
    # Busca el último número en el nombre del archivo
    match = re.search(r'_(\d+)\.csv$', path)
    return int(match.group(1)) if match else -1

def animate_field_propagation(mode):
    # Buscar todos los archivos de datos para este modo
    data_files = sorted(glob.glob(f'data/snap/field_data_{mode}_*.csv'))

    # Ordenar por el número al final
    data_files = sorted(data_files, key=numeric_sort_key)
    
    if not data_files:
        print(f"No data files found for mode {mode}")
        return
    
    # Leer todos los datos
    all_data = []
    for file in data_files:
        data = pd.read_csv(file)
        all_data.append(data)
    
    # Crear figura con subplots
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))
    
    # Configurar ejes para campos eléctricos
    axes[0].set_xlim(0, all_data[0]['z'].max())
    axes[0].set_ylim(-1.5, 1.5)
    axes[0].set_ylabel('Electric Field (V/m)')
    axes[0].grid(True)
    axes[0].set_title(f'Propagation of Mode {mode}')
    
    # Configurar ejes para campos magnéticos
    axes[1].set_xlim(0, all_data[0]['z'].max())
    
    # Establecer límites adecuados para campos magnéticos
    if mode in ['R', 'L', 'O']:
        axes[1].set_ylim(-2e-9, 2e-9)  # Para modos R, L, O
    else:  # Modo X
        axes[1].set_ylim(-1.0, 1.0)    # Para modo X (Bx más grande)
    
    axes[1].set_ylabel('Magnetic Field (T)')
    axes[1].grid(True)
    
    # Configurar ejes para corrientes
    axes[2].set_xlim(0, all_data[0]['z'].max())
    
    # Establecer límites adecuados para corrientes
    if mode in ['R', 'L', 'O']:
        axes[2].set_ylim(-0.2, 0.2)  # Para modos R, L, O
    else:  # Modo X
        axes[2].set_ylim(-2e6, 2e6)  # Para modo X (corrientes más grandes)
    
    axes[2].set_xlabel('Position (m)')
    axes[2].set_ylabel('Current Density (A/m²)')
    axes[2].grid(True)
    
    # Inicializar líneas
    # Campos eléctricos
    line_ex, = axes[0].plot([], [], 'b-', label='Ex')
    line_ey, = axes[0].plot([], [], 'r-', label='Ey')
    axes[0].legend()
    
    # Campos magnéticos
    line_bx, = axes[1].plot([], [], 'c-', label='Bx')
    line_by, = axes[1].plot([], [], 'm-', label='By')
    axes[1].legend()
    
    # Corrientes
    line_jx_e, = axes[2].plot([], [], 'b-', label='Jx_e')
    line_jx_i, = axes[2].plot([], [], 'r-', label='Jx_i')
    axes[2].legend()
    
    # Texto para mostrar el tiempo
    time_text = axes[0].text(0.02, 0.95, '', transform=axes[0].transAxes)
    
    # Función de inicialización
    def init():
        line_ex.set_data([], [])
        line_ey.set_data([], [])
        line_bx.set_data([], [])
        line_by.set_data([], [])
        line_jx_e.set_data([], [])
        line_jx_i.set_data([], [])
        return line_ex, line_ey, line_bx, line_by, line_jx_e, line_jx_i, time_text
    
    # Función de animación
    def animate(i):
        data = all_data[i]
        time = i * 3.3389e-11  # Aproximación del tiempo
        
        line_ex.set_data(data['z'], data['Ex'])
        line_ey.set_data(data['z'], data['Ey'])
        line_bx.set_data(data['z'], data['Bx'])
        line_by.set_data(data['z'], data['By'])
        line_jx_e.set_data(data['z'], data['Jx_e'])
        line_jx_i.set_data(data['z'], data['Jx_i'])
        time_text.set_text(f'Time = {time:.2e} s')
        
        return line_ex, line_ey, line_bx, line_by, line_jx_e, line_jx_i, time_text
    
    # Crear animación
    ani = FuncAnimation(fig, animate, frames=len(all_data),
                        init_func=init, interval=200, blit=True)
    
    plt.tight_layout()
    
    # Guardar animación
    try:
        os.makedirs('animations', exist_ok=True)
        ani.save(f'animations/propagation_{mode}.gif', writer='pillow', fps=5)
        print(f"Animation saved as animations/propagation_{mode}.gif")
    except Exception as e:
        print(f"Error saving animation: {e}")
    
    plt.show()

if __name__ == "__main__":
    for mode in ['R', 'L', 'O', 'X']:
        print(f"Generating animation for mode {mode}...")
        animate_field_propagation(mode)
