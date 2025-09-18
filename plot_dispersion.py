import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_dispersion_relations():
    # Leer datos de dispersión
    data = pd.read_csv('data/dispersion_data.csv')
    
    # Crear figura
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Modo R
    axes[0,0].plot(data['frequency'], data['Re(k_R)'], 'b-', label='Parte real')
    axes[0,0].plot(data['frequency'], data['Im(k_R)'], 'r--', label='Parte imaginaria')
    axes[0,0].set_title('Modo R (Right-hand)')
    axes[0,0].set_xlabel('Frecuencia (Hz)')
    axes[0,0].set_ylabel('k (m⁻¹)')
    axes[0,0].legend()
    axes[0,0].grid(True)
    axes[0,0].set_xscale('log')
    
    # Modo L
    axes[0,1].plot(data['frequency'], data['Re(k_L)'], 'b-', label='Parte real')
    axes[0,1].plot(data['frequency'], data['Im(k_L)'], 'r--', label='Parte imaginaria')
    axes[0,1].set_title('Modo L (Left-hand)')
    axes[0,1].set_xlabel('Frecuencia (Hz)')
    axes[0,1].set_ylabel('k (m⁻¹)')
    axes[0,1].legend()
    axes[0,1].grid(True)
    axes[0,1].set_xscale('log')
    
    # Modo O
    axes[1,0].plot(data['frequency'], data['Re(k_O)'], 'b-', label='Parte real')
    axes[1,0].plot(data['frequency'], data['Im(k_O)'], 'r--', label='Parte imaginaria')
    axes[1,0].set_title('Modo O (Ordinary)')
    axes[1,0].set_xlabel('Frecuencia (Hz)')
    axes[1,0].set_ylabel('k (m⁻¹)')
    axes[1,0].legend()
    axes[1,0].grid(True)
    axes[1,0].set_xscale('log')
    
    # Modo X
    axes[1,1].plot(data['frequency'], data['Re(k_X)'], 'b-', label='Parte real')
    axes[1,1].plot(data['frequency'], data['Im(k_X)'], 'r--', label='Parte imaginaria')
    axes[1,1].set_title('Modo X (Extraordinary)')
    axes[1,1].set_xlabel('Frecuencia (Hz)')
    axes[1,1].set_ylabel('k (m⁻¹)')
    axes[1,1].legend()
    axes[1,1].grid(True)
    axes[1,1].set_xscale('log')
    
    plt.tight_layout()
    plt.savefig('data/dispersion_relations.png', dpi=300)
    plt.show()

if __name__ == "__main__":
    plot_dispersion_relations()
