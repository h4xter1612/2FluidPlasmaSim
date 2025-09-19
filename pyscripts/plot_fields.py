import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_field_data(mode):
    # Leer datos de campo
    data = pd.read_csv('data/'+f'field_data_{mode}.csv')
    
    # Crear figura
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
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
    
    # Corrientes
    axes[1,0].plot(data['z'], data['Jx'], 'b-', label='Jx')
    axes[1,0].plot(data['z'], data['Jy'], 'r-', label='Jy')
    axes[1,0].plot(data['z'], data['Jz'], 'g-', label='Jz')
    axes[1,0].set_title(f'Corrientes - Modo {mode}')
    axes[1,0].set_xlabel('Position (m)')
    axes[1,0].set_ylabel('Current (A/m²)')
    axes[1,0].legend()
    axes[1,0].grid(True)
    
    # Diagrama de fases para modos circulares
    if mode in ['R', 'L']:
        axes[1,1].plot(data['Ex'], data['Ey'], 'b-')
        axes[1,1].set_title(f'Phase Diagram Ex-Ey - Mode {mode}')
        axes[1,1].set_xlabel('Ex (V/m)')
        axes[1,1].set_ylabel('Ey (V/m)')
        axes[1,1].grid(True)
        axes[1,1].axis('equal')
    else:
        axes[1,1].remove()
    
    plt.tight_layout()
    plt.savefig('data/'+f'field_plot_{mode}.png', dpi=300)
    plt.show()

if __name__ == "__main__":
    for mode in ['R', 'L', 'O', 'X']:
        plot_field_data(mode)
