# animate_propagation_fixed.py
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
        print(f"No se encontraron archivos de datos para el modo {mode}")
        return
    
    # Leer todos los datos
    all_data = []
    for file in data_files:
        data = pd.read_csv(file)
        all_data.append(data)
    
    # Crear figura con subplots independientes
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
    
    # Configurar ejes para campos eléctricos
    ax1.set_xlim(0, all_data[0]['z'].max())
    ax1.set_ylim(-1.5, 1.5)
    ax1.set_ylabel('Campo Eléctrico (V/m)')
    ax1.grid(True)
    ax1.set_title(f'Propagación del Modo {mode}')
    
    # Configurar ejes para campos magnéticos (con escala apropiada)
    ax2.set_xlim(0, all_data[0]['z'].max())
    
    # Establecer límites adecuados para campos magnéticos
    if mode in ['R', 'L', 'O']:
        ax2.set_ylim(-2e-9, 2e-9)  # Para modos R, L, O
    else:  # Modo X
        ax2.set_ylim(-1.0, 1.0)    # Para modo X (Bx más grande)
    
    ax2.set_xlabel('Posición (m)')
    ax2.set_ylabel('Campo Magnético (T)')
    ax2.grid(True)
    
    # Inicializar líneas
    if mode in ['R', 'L']:
        # Campos eléctricos
        line_ex, = ax1.plot([], [], 'b-', label='Ex')
        line_ey, = ax1.plot([], [], 'r-', label='Ey')
        ax1.legend()
        
        # Campos magnéticos
        line_bx, = ax2.plot([], [], 'c-', label='Bx')
        line_by, = ax2.plot([], [], 'm-', label='By')
        ax2.legend()
        
    elif mode == 'O':
        # Modo ordinario (solo Ez)
        line_ez, = ax1.plot([], [], 'g-', label='Ez')
        ax1.legend()
        
        # Para modo O, los campos magnéticos son pequeños
        line_bx, = ax2.plot([], [], 'c-', label='Bx')
        line_by, = ax2.plot([], [], 'm-', label='By')
        ax2.legend()
        
    else:  # Modo X
        # Campos eléctricos
        line_ex, = ax1.plot([], [], 'b-', label='Ex')
        ax1.legend()
        
        # Campos magnéticos (Bx es dominante en modo X)
        line_bx, = ax2.plot([], [], 'r-', label='Bx')
        line_by, = ax2.plot([], [], 'm-', label='By')
        ax2.legend()
    
    # Texto para mostrar el tiempo
    time_text = ax1.text(0.02, 0.95, '', transform=ax1.transAxes)
    
    # Función de inicialización
    def init():
        if mode in ['R', 'L']:
            line_ex.set_data([], [])
            line_ey.set_data([], [])
            line_bx.set_data([], [])
            line_by.set_data([], [])
            return line_ex, line_ey, line_bx, line_by, time_text
        elif mode == 'O':
            line_ez.set_data([], [])
            line_bx.set_data([], [])
            line_by.set_data([], [])
            return line_ez, line_bx, line_by, time_text
        else:
            line_ex.set_data([], [])
            line_bx.set_data([], [])
            line_by.set_data([], [])
            return line_ex, line_bx, line_by, time_text
    
    # Función de animación
    def animate(i):
        data = all_data[i]
        time = i * 3.3389e-11  # Aproximación del tiempo
        
        if mode in ['R', 'L']:
            line_ex.set_data(data['z'], data['Ex'])
            line_ey.set_data(data['z'], data['Ey'])
            line_bx.set_data(data['z'], data['Bx'])
            line_by.set_data(data['z'], data['By'])
            time_text.set_text(f'Tiempo = {time:.2e} s')
            return line_ex, line_ey, line_bx, line_by, time_text
            
        elif mode == 'O':
            line_ez.set_data(data['z'], data['Ez'])
            line_bx.set_data(data['z'], data['Bx'])
            line_by.set_data(data['z'], data['By'])
            time_text.set_text(f'Tiempo = {time:.2e} s')
            return line_ez, line_bx, line_by, time_text
            
        else:  # Modo X
            line_ex.set_data(data['z'], data['Ex'])
            line_bx.set_data(data['z'], data['Bx'])
            line_by.set_data(data['z'], data['By'])
            time_text.set_text(f'Tiempo = {time:.2e} s')
            return line_ex, line_bx, line_by, time_text
    
    # Crear animación
    ani = FuncAnimation(fig, animate, frames=len(all_data),
                        init_func=init, interval=200, blit=True)
    
    plt.tight_layout()
    
    # Guardar animación
    try:
        os.makedirs('animations', exist_ok=True)
        ani.save(f'animations/propagation_{mode}.gif', writer='pillow', fps=5)
        print(f"Animación guardada como animations/propagation_{mode}.gif")
    except Exception as e:
        print(f"Error al guardar la animación: {e}")
    
    plt.show()

if __name__ == "__main__":
    for mode in ['R', 'L', 'O', 'X']:
        print(f"Generando animación para el modo {mode}...")
        animate_field_propagation(mode)
