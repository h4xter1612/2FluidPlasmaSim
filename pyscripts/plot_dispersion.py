import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_dispersion_relations():
    # Leer datos de dispersión
    data = pd.read_csv('data/dispersion_data.csv')
    
    # Crear figura
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Mode R
    axes[0,0].plot(data['frequency'], data['Re(k_R)'], 'b-', label='Real part')
    axes[0,0].plot(data['frequency'], data['Im(k_R)'], 'r--', label='Imaginary part')
    axes[0,0].set_title('Mode R (Right-hand)')
    axes[0,0].set_xlabel('Frequency (Hz)')
    axes[0,0].set_ylabel('k (m⁻¹)')
    axes[0,0].legend()
    axes[0,0].grid(True)
    axes[0,0].set_xscale('log')
    
    # Mode L
    axes[0,1].plot(data['frequency'], data['Re(k_L)'], 'b-', label='Real part')
    axes[0,1].plot(data['frequency'], data['Im(k_L)'], 'r--', label='Imaginary part')
    axes[0,1].set_title('Mode L (Left-hand)')
    axes[0,1].set_xlabel('Frequency (Hz)')
    axes[0,1].set_ylabel('k (m⁻¹)')
    axes[0,1].legend()
    axes[0,1].grid(True)
    axes[0,1].set_xscale('log')
    
    # Mode O
    axes[1,0].plot(data['frequency'], data['Re(k_O)'], 'b-', label='Real part')
    axes[1,0].plot(data['frequency'], data['Im(k_O)'], 'r--', label='Imaginary part')
    axes[1,0].set_title('Mode O (Ordinary)')
    axes[1,0].set_xlabel('Frequency (Hz)')
    axes[1,0].set_ylabel('k (m⁻¹)')
    axes[1,0].legend()
    axes[1,0].grid(True)
    axes[1,0].set_xscale('log')
    
    # Mode X
    axes[1,1].plot(data['frequency'], data['Re(k_X)'], 'b-', label='Real part')
    axes[1,1].plot(data['frequency'], data['Im(k_X)'], 'r--', label='Imaginary part')
    axes[1,1].set_title('Mode X (Extraordinary)')
    axes[1,1].set_xlabel('Frequency (Hz)')
    axes[1,1].set_ylabel('k (m⁻¹)')
    axes[1,1].legend()
    axes[1,1].grid(True)
    axes[1,1].set_xscale('log')
    
    plt.tight_layout()
    plt.savefig('data/dispersion_relations.png', dpi=300)
    plt.show()

if __name__ == "__main__":
    plot_dispersion_relations()
