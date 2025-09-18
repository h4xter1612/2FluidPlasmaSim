import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd

def animate_field_propagation(mode):
    # Leer datos de campo
    data = pd.read_csv('data/'+f'field_data_{mode}.csv')
    
    # Crear animación
    fig, ax = plt.subplots(figsize=(10, 6))
    
    if mode in ['R', 'L']:
        line, = ax.plot(data['z'], data['Ex'], 'b-', label='Ex')
        line2, = ax.plot(data['z'], data['Ey'], 'r-', label='Ey')
        ax.set_ylabel('Campo (V/m)')
    elif mode == 'O':
        line, = ax.plot(data['z'], data['Ez'], 'g-', label='Ez')
        ax.set_ylabel('Campo (V/m)')
    else:  # Modo X
        line, = ax.plot(data['z'], data['Ex'], 'b-', label='Ex')
        line2, = ax.plot(data['z'], data['Bx'], 'r-', label='Bx')
        ax.set_ylabel('Campo')
    
    ax.set_xlabel('Posición (m)')
    ax.set_title(f'Propagación del Modo {mode}')
    ax.grid(True)
    ax.legend()
    
    def update(frame):
        # Esta función se actualizaría con datos de múltiples frames en una simulación real
        # Aquí solo mostramos un frame estático como ejemplo
        return line,
    
    ani = FuncAnimation(fig, update, frames=1, interval=100, blit=True)
    plt.show()
    
    # Guardar animación (requiere ffmpeg)
    # ani.save(f'propagation_{mode}.mp4', writer='ffmpeg', fps=10)

if __name__ == "__main__":
    for mode in ['R', 'L', 'O', 'X']:
        animate_field_propagation(mode)
