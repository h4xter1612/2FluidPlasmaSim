import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot_field_data_2d(mode):
    # Leer datos de campo 2D
    data = pd.read_csv(f'data/field_data_2d_{mode}.csv')
    
    # Obtener dimensiones de la malla
    x_unique = np.sort(data['x'].unique())
    y_unique = np.sort(data['y'].unique())
    nx = len(x_unique)
    ny = len(y_unique)
    
    # Reorganizar datos en matrices 2D
    Ex = data['Ex'].values.reshape(nx, ny)
    Ey = data['Ey'].values.reshape(nx, ny)
    Ez = data['Ez'].values.reshape(nx, ny)
    Bx = data['Bx'].values.reshape(nx, ny)
    By = data['By'].values.reshape(nx, ny)
    Bz = data['Bz'].values.reshape(nx, ny)
    Jx_e = data['Jx_e'].values.reshape(nx, ny)
    Jy_e = data['Jy_e'].values.reshape(nx, ny)
    Jz_e = data['Jz_e'].values.reshape(nx, ny)
    Jx_i = data['Jx_i'].values.reshape(nx, ny)
    Jy_i = data['Jy_i'].values.reshape(nx, ny)
    Jz_i = data['Jz_i'].values.reshape(nx, ny)
    
    # Crear figura con subplots
    fig, axes = plt.subplots(4, 3, figsize=(18, 20))
    axes = axes.flatten()
    
    # Lista de campos y títulos
    fields = [Ex, Ey, Ez, Bx, By, Bz, Jx_e, Jy_e, Jz_e, Jx_i, Jy_i, Jz_i]
    titles = [
        f'Ex - Mode {mode}', f'Ey - Mode {mode}', f'Ez - Mode {mode}',
        f'Bx - Mode {mode}', f'By - Mode {mode}', f'Bz - Mode {mode}',
        f'Jx_e - Mode {mode}', f'Jy_e - Mode {mode}', f'Jz_e - Mode {mode}',
        f'Jx_i - Mode {mode}', f'Jy_i - Mode {mode}', f'Jz_i - Mode {mode}'
    ]
    
    # Plot cada campo
    for i, (field, title) in enumerate(zip(fields, titles)):
        im = axes[i].imshow(field.T, origin='lower', 
                          extent=[x_unique[0], x_unique[-1], y_unique[0], y_unique[-1]],
                          cmap='RdBu_r')
        axes[i].set_title(title)
        axes[i].set_xlabel('x (m)')
        axes[i].set_ylabel('y (m)')
        
        # Añadir barra de color
        divider = make_axes_locatable(axes[i])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax)
    
    plt.tight_layout()
    os.makedirs('plots_2d', exist_ok=True)
    plt.savefig(f'plots_2d/field_plot_2d_{mode}.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_field_comparison_2d(mode):
    # Leer datos de campo 2D
    data = pd.read_csv(f'data/field_data_2d_{mode}.csv')
    
    # Obtener dimensiones de la malla
    x_unique = np.sort(data['x'].unique())
    y_unique = np.sort(data['y'].unique())
    nx = len(x_unique)
    ny = len(y_unique)
    
    # Reorganizar datos en matrices 2D
    Ex = data['Ex'].values.reshape(nx, ny)
    Ey = data['Ey'].values.reshape(nx, ny)
    Bx = data['Bx'].values.reshape(nx, ny)
    By = data['By'].values.reshape(nx, ny)
    Jx_e = data['Jx_e'].values.reshape(nx, ny)
    Jy_e = data['Jy_e'].values.reshape(nx, ny)
    Jx_i = data['Jx_i'].values.reshape(nx, ny)
    Jy_i = data['Jy_i'].values.reshape(nx, ny)
    
    # Crear figura con subplots
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    
    # Campos eléctricos
    im1 = axes[0, 0].imshow(Ex.T, origin='lower', 
                          extent=[x_unique[0], x_unique[-1], y_unique[0], y_unique[-1]],
                          cmap='RdBu_r')
    axes[0, 0].set_title(f'Ex - Mode {mode}')
    axes[0, 0].set_ylabel('y (m)')
    plt.colorbar(im1, ax=axes[0, 0])
    
    im2 = axes[0, 1].imshow(Ey.T, origin='lower', 
                          extent=[x_unique[0], x_unique[-1], y_unique[0], y_unique[-1]],
                          cmap='RdBu_r')
    axes[0, 1].set_title(f'Ey - Mode {mode}')
    plt.colorbar(im2, ax=axes[0, 1])
    
    # Campos magnéticos
    im3 = axes[0, 2].imshow(Bx.T, origin='lower', 
                          extent=[x_unique[0], x_unique[-1], y_unique[0], y_unique[-1]],
                          cmap='RdBu_r')
    axes[0, 2].set_title(f'Bx - Mode {mode}')
    plt.colorbar(im3, ax=axes[0, 2])
    
    im4 = axes[0, 3].imshow(By.T, origin='lower', 
                          extent=[x_unique[0], x_unique[-1], y_unique[0], y_unique[-1]],
                          cmap='RdBu_r')
    axes[0, 3].set_title(f'By - Mode {mode}')
    plt.colorbar(im4, ax=axes[0, 3])
    
    # Corrientes electrónicas
    im5 = axes[1, 0].imshow(Jx_e.T, origin='lower', 
                          extent=[x_unique[0], x_unique[-1], y_unique[0], y_unique[-1]],
                          cmap='RdBu_r')
    axes[1, 0].set_title(f'Jx_e - Mode {mode}')
    axes[1, 0].set_xlabel('x (m)')
    axes[1, 0].set_ylabel('y (m)')
    plt.colorbar(im5, ax=axes[1, 0])
    
    im6 = axes[1, 1].imshow(Jy_e.T, origin='lower', 
                          extent=[x_unique[0], x_unique[-1], y_unique[0], y_unique[-1]],
                          cmap='RdBu_r')
    axes[1, 1].set_title(f'Jy_e - Mode {mode}')
    axes[1, 1].set_xlabel('x (m)')
    plt.colorbar(im6, ax=axes[1, 1])
    
    # Corrientes iónicas
    im7 = axes[1, 2].imshow(Jx_i.T, origin='lower', 
                          extent=[x_unique[0], x_unique[-1], y_unique[0], y_unique[-1]],
                          cmap='RdBu_r')
    axes[1, 2].set_title(f'Jx_i - Mode {mode}')
    axes[1, 2].set_xlabel('x (m)')
    plt.colorbar(im7, ax=axes[1, 2])
    
    im8 = axes[1, 3].imshow(Jy_i.T, origin='lower', 
                          extent=[x_unique[0], x_unique[-1], y_unique[0], y_unique[-1]],
                          cmap='RdBu_r')
    axes[1, 3].set_title(f'Jy_i - Mode {mode}')
    axes[1, 3].set_xlabel('x (m)')
    plt.colorbar(im8, ax=axes[1, 3])
    
    plt.tight_layout()
    os.makedirs('plots_2d', exist_ok=True)
    plt.savefig(f'plots_2d/field_comparison_2d_{mode}.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    for mode in ['R', 'L', 'O', 'X']:
        try:
            print(f"Plotting 2D field data for mode {mode}...")
            plot_field_data_2d(mode)
            plot_field_comparison_2d(mode)
        except Exception as e:
            print(f"Error plotting mode {mode}: {e}")
