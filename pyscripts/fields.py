import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

def plot_field_data(mode):
    # Leer datos de campo
    data = pd.read_csv('data/'+f'field_data_{mode}.csv')
    
    # Crear figura
    fig, axes = plt.subplots(3, 2, figsize=(15, 12))
    
    # Campos eléctricos
    axes[0,0].plot(data['z'], data['Ex'], 'b-', label='Ex')
    axes[0,0].plot(data['z'], data['Ey'], 'r-', label='Ey')
    axes[0,0].plot(data['z'], data['Ez'], 'g-', label='Ez')
    axes[0,0].set_title(f'Electric Fields - Mode {mode}')
    axes[0,0].set_xlabel('Position (m)')
    axes[0,0].set_ylabel('Field (V/m)')
    axes[0,0].legend()
    axes[0,0].grid(True)
    
    # Campos magnéticos
    axes[0,1].plot(data['z'], data['Bx'], 'b-', label='Bx')
    axes[0,1].plot(data['z'], data['By'], 'r-', label='By')
    axes[0,1].plot(data['z'], data['Bz'], 'g-', label='Bz')
    axes[0,1].set_title(f'Magnetic Fields - Mode {mode}')
    axes[0,1].set_xlabel('Position (m)')
    axes[0,1].set_ylabel('Field (T)')
    axes[0,1].legend()
    axes[0,1].grid(True)
    
    # Corrientes electrónicas
    axes[1,0].plot(data['z'], data['Jx_e'], 'b-', label='Jx_e')
    axes[1,0].plot(data['z'], data['Jy_e'], 'r-', label='Jy_e')
    axes[1,0].plot(data['z'], data['Jz_e'], 'g-', label='Jz_e')
    axes[1,0].set_title(f'Electron Currents - Mode {mode}')
    axes[1,0].set_xlabel('Position (m)')
    axes[1,0].set_ylabel('Current (A/m²)')
    axes[1,0].legend()
    axes[1,0].grid(True)
    
    # Corrientes iónicas
    axes[1,1].plot(data['z'], data['Jx_i'], 'b-', label='Jx_i')
    axes[1,1].plot(data['z'], data['Jy_i'], 'r-', label='Jy_i')
    axes[1,1].plot(data['z'], data['Jz_i'], 'g-', label='Jz_i')
    axes[1,1].set_title(f'Ion Currents - Mode {mode}')
    axes[1,1].set_xlabel('Position (m)')
    axes[1,1].set_ylabel('Current (A/m²)')
    axes[1,1].legend()
    axes[1,1].grid(True)
    
    # Diagrama de fases para modos circulares
    if mode in ['R', 'L']:
        axes[2,0].plot(data['Ex'], data['Ey'], 'b-')
        axes[2,0].set_title(f'Phase Diagram Ex-Ey - Mode {mode}')
        axes[2,0].set_xlabel('Ex (V/m)')
        axes[2,0].set_ylabel('Ey (V/m)')
        axes[2,0].grid(True)
        axes[2,0].axis('equal')
        
        # Comparación de corrientes electrónicas e iónicas
        axes[2,1].plot(data['z'], np.abs(data['Jx_e']), 'b-', label='|Jx_e|')
        axes[2,1].plot(data['z'], np.abs(data['Jx_i']), 'r-', label='|Jx_i|')
        axes[2,1].set_title(f'Current Comparison - Mode {mode}')
        axes[2,1].set_xlabel('Position (m)')
        axes[2,1].set_ylabel('Current Magnitude (A/m²)')
        axes[2,1].legend()
        axes[2,1].grid(True)
        axes[2,1].set_yscale('log')
    else:
        # Para modos O y X, mostrar solo el diagrama de fases
        axes[2,0].remove()
        axes[2,1].remove()
    
    plt.tight_layout()
    os.makedirs('plots', exist_ok=True)
    plt.savefig('plots/'+f'field_plot_{mode}.png', dpi=300)
    plt.show()

if __name__ == "__main__":
    for mode in ['R', 'L', 'O', 'X']:
        plot_field_data(mode)
