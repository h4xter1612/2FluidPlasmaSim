import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd
import glob
import os
import re
from mpl_toolkits.axes_grid1 import make_axes_locatable

def numeric_sort_key(path):
    # Busca el último número en el nombre del archivo
    match = re.search(r'_(\d+)\.csv$', path)
    return int(match.group(1)) if match else -1

def animate_field_propagation_2d(mode):
    # Buscar todos los archivos de datos para este modo
    data_files = sorted(glob.glob(f'data/snap/field_data_2d_{mode}_*.csv'))
    
    # Ordenar por el número al final
    data_files = sorted(data_files, key=numeric_sort_key)
    
    if not data_files:
        print(f"No data files found for mode {mode}")
        return
    
    # Leer el primer archivo para obtener la estructura de la malla
    first_data = pd.read_csv(data_files[0])
    x_unique = np.sort(first_data['x'].unique())
    y_unique = np.sort(first_data['y'].unique())
    nx = len(x_unique)
    ny = len(y_unique)
    
    # Crear figura con subplots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    # Configurar títulos
    titles = [
        f'Ex - Mode {mode}',
        f'Ey - Mode {mode}',
        f'Bx - Mode {mode}',
        f'By - Mode {mode}'
    ]
    
    # Inicializar imágenes
    images = []
    for i, ax in enumerate(axes):
        im = ax.imshow(np.zeros((nx, ny)).T, origin='lower', 
                      extent=[x_unique[0], x_unique[-1], y_unique[0], y_unique[-1]],
                      cmap='RdBu_r', animated=True)
        ax.set_title(titles[i])
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        
        # Añadir barra de color
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax)
        
        images.append(im)
    
    # Texto para mostrar el tiempo
    time_text = axes[0].text(0.02, 0.95, '', transform=axes[0].transAxes, color='white',
                           bbox=dict(facecolor='black', alpha=0.7))
    
    # Función de inicialización
    def init():
        for im in images:
            im.set_array(np.zeros((nx, ny)).T)
        time_text.set_text('')
        return images + [time_text]
    
    # Función de animación
    def animate(i):
        data = pd.read_csv(data_files[i])
        
        # Reorganizar datos en matrices 2D
        Ex = data['Ex'].values.reshape(nx, ny)
        Ey = data['Ey'].values.reshape(nx, ny)
        Bx = data['Bx'].values.reshape(nx, ny)
        By = data['By'].values.reshape(nx, ny)
        
        # Actualizar imágenes
        images[0].set_array(Ex.T)
        images[1].set_array(Ey.T)
        images[2].set_array(Bx.T)
        images[3].set_array(By.T)
        
        # Obtener el tiempo del nombre del archivo
        time_match = re.search(r'_(\d+)\.csv$', data_files[i])
        time_step = int(time_match.group(1)) if time_match else i
        time = time_step * 3.3389e-11  # Aproximación del tiempo
        
        time_text.set_text(f'Time = {time:.2e} s')
        
        return images + [time_text]
    
    # Crear animación
    ani = FuncAnimation(fig, animate, frames=len(data_files),
                        init_func=init, interval=200, blit=True)
    
    plt.tight_layout()
    
    # Guardar animación
    try:
        os.makedirs('animations_2d', exist_ok=True)
        ani.save(f'animations_2d/propagation_2d_{mode}.gif', writer='pillow', fps=5)
        print(f"Animation saved as animations_2d/propagation_2d_{mode}.gif")
    except Exception as e:
        print(f"Error saving animation: {e}")
    
    plt.show()

def create_propagation_montage_2d(mode):
    # Buscar todos los archivos de datos para este modo
    data_files = sorted(glob.glob(f'data/snap/field_data_2d_{mode}_*.csv'))
    
    # Ordenar por el número al final
    data_files = sorted(data_files, key=numeric_sort_key)
    
    if not data_files:
        print(f"No data files found for mode {mode}")
        return
    
    # Seleccionar un subconjunto de archivos para el montaje
    step = max(1, len(data_files) // 9)  # 9 imágenes en el montaje
    selected_files = data_files[::step][:9]
    
    # Leer el primer archivo para obtener la estructura de la malla
    first_data = pd.read_csv(selected_files[0])
    x_unique = np.sort(first_data['x'].unique())
    y_unique = np.sort(first_data['y'].unique())
    nx = len(x_unique)
    ny = len(y_unique)
    
    # Crear figura con subplots
    fig, axes = plt.subplots(3, 3, figsize=(15, 15))
    axes = axes.flatten()
    
    # Plot cada paso de tiempo
    for i, (file, ax) in enumerate(zip(selected_files, axes)):
        data = pd.read_csv(file)
        Ex = data['Ex'].values.reshape(nx, ny)
        
        im = ax.imshow(Ex.T, origin='lower', 
                      extent=[x_unique[0], x_unique[-1], y_unique[0], y_unique[-1]],
                      cmap='RdBu_r')
        ax.set_title(f'Time step {i*step}')
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        
        # Añadir barra de color
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax)
    
    plt.suptitle(f'Propagation of Ex - Mode {mode}', fontsize=16)
    plt.tight_layout()
    
    # Guardar montaje
    os.makedirs('montages_2d', exist_ok=True)
    plt.savefig(f'montages_2d/propagation_montage_2d_{mode}.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    for mode in ['R', 'L', 'O', 'X']:
        print(f"Generating 2D animation for mode {mode}...")
        animate_field_propagation_2d(mode)
        create_propagation_montage_2d(mode)
