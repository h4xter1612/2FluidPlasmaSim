import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

def detailed_energy_analysis():
    data = pd.read_csv('data/energy_evolution.csv')
    
    # Análisis detallado de la evolución temporal
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    
    # 1. Energía total con derivada temporal
    axes[0,0].plot(data['time'], data['total_energy'], 'k-', linewidth=2, label='Total Energy')
    axes[0,0].set_xlabel('Time (s)')
    axes[0,0].set_ylabel('Total Energy (J/m)')
    axes[0,0].set_title('Total Energy Evolution')
    axes[0,0].grid(True, alpha=0.3)
    
    # 2. Derivada temporal de la energía (¡esto es clave!)
    time_diff = np.diff(data['time'])
    energy_diff = np.diff(data['total_energy'])
    energy_derivative = energy_diff / time_diff
    
    axes[0,1].plot(data['time'][1:], energy_derivative, 'r-', linewidth=2)
    axes[0,1].set_xlabel('Time (s)')
    axes[0,1].set_ylabel('dE/dt (W/m)')
    axes[0,1].set_title('Energy Rate of Change')
    axes[0,1].grid(True, alpha=0.3)
    
    # 3. Energías individuales (escala log)
    axes[0,2].semilogy(data['time'], data['em_energy'] + 1e-20, 'b-', label='EM')
    axes[0,2].semilogy(data['time'], data['electron_energy'] + 1e-20, 'r-', label='Electron')
    axes[0,2].semilogy(data['time'], data['ion_energy'] + 1e-20, 'g-', label='Ion')
    axes[0,2].set_xlabel('Time (s)')
    axes[0,2].set_ylabel('Energy (J/m)')
    axes[0,2].set_title('Component Energies (Log Scale)')
    axes[0,2].legend()
    axes[0,2].grid(True, alpha=0.3)
    
    # 4. Conservación de energía por componentes
    total_initial = data['total_energy'].iloc[0]
    if total_initial > 0:
        relative_em = (data['em_energy'] - data['em_energy'].iloc[0]) / total_initial
        relative_elec = (data['electron_energy'] - data['electron_energy'].iloc[0]) / total_initial
        relative_ion = (data['ion_energy'] - data['ion_energy'].iloc[0]) / total_initial
        relative_total = (data['total_energy'] - data['total_energy'].iloc[0]) / total_initial
        
        axes[1,0].plot(data['time'], relative_em * 100, 'b-', label='ΔEM')
        axes[1,0].plot(data['time'], relative_elec * 100, 'r-', label='ΔElectron')
        axes[1,0].plot(data['time'], relative_ion * 100, 'g-', label='ΔIon')
        axes[1,0].plot(data['time'], relative_total * 100, 'k-', label='ΔTotal', linewidth=2)
        axes[1,0].set_xlabel('Time (s)')
        axes[1,0].set_ylabel('Energy Change (% of Initial)')
        axes[1,0].set_title('Energy Conservation by Component')
        axes[1,0].legend()
        axes[1,0].grid(True, alpha=0.3)
    
    # 5. Análisis espectral de las fluctuaciones
    from scipy import signal
    if len(data) > 10:
        total_energy = data['total_energy'].values
        # Remover tendencia lineal
        detrended = signal.detrend(total_energy)
        # Calcular FFT
        sample_rate = 1.0 / np.mean(np.diff(data['time']))
        freqs = np.fft.fftfreq(len(detrended), 1/sample_rate)
        fft_vals = np.abs(np.fft.fft(detrended))
        
        positive_freqs = freqs[:len(freqs)//2]
        positive_fft = fft_vals[:len(fft_vals)//2]
        
        axes[1,1].semilogy(positive_freqs, positive_fft)
        axes[1,1].set_xlabel('Frequency (Hz)')
        axes[1,1].set_ylabel('FFT Magnitude')
        axes[1,1].set_title('Spectral Analysis of Energy Fluctuations')
        axes[1,1].grid(True, alpha=0.3)
    
    # 6. Estadísticas de conservación
    axes[1,2].text(0.1, 0.9, f"Initial Energy: {data['total_energy'].iloc[0]:.6e} J/m", 
                   transform=axes[1,2].transAxes, fontsize=12)
    axes[1,2].text(0.1, 0.7, f"Final Energy: {data['total_energy'].iloc[-1]:.6e} J/m", 
                   transform=axes[1,2].transAxes, fontsize=12)
    axes[1,2].text(0.1, 0.5, f"Max Variation: {((data['total_energy'].max() - data['total_energy'].min()) / data['total_energy'].iloc[0]) * 100:.2f}%", 
                   transform=axes[1,2].transAxes, fontsize=12)
    axes[1,2].text(0.1, 0.3, f"Std Dev: {data['total_energy'].std():.6e} J/m", 
                   transform=axes[1,2].transAxes, fontsize=12)
    axes[1,2].set_xlim(0, 1)
    axes[1,2].set_ylim(0, 1)
    axes[1,2].set_title('Energy Statistics')
    axes[1,2].axis('off')
    
    plt.tight_layout()
    os.makedirs('plots', exist_ok=True)
    plt.savefig('plots/detailed_energy_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Diagnóstico automático
    print("=== ENERGY CONSERVATION DIAGNOSIS ===")
    
    # Verificar si hay crecimiento exponencial
    energy_growth = data['total_energy'].iloc[-1] / data['total_energy'].iloc[10]  # Ignorar primeros puntos
    if energy_growth > 10:
        print("❌ CRITICAL: Exponential energy growth detected!")
        print("   This indicates numerical instability.")
        print("   Solutions: Reduce time step, check CFL condition")
    
    # Verificar fluctuaciones
    energy_std = data['total_energy'].std()
    energy_mean = data['total_energy'].mean()
    if energy_std / energy_mean > 0.1:  # Más del 10% de fluctuación
        print("❌ HIGH fluctuations in total energy")
        print("   Possible causes:")
        print("   - Time step too large")
        print("   - PML boundary conditions too strong")
        print("   - Numerical scheme instability")
    
    # Verificar si las componentes suman correctamente
    computed_total = data['em_energy'] + data['electron_energy'] + data['ion_energy']
    discrepancy = np.max(np.abs(computed_total - data['total_energy']))
    if discrepancy > 1e-10:
        print(f"⚠️  Energy components don't sum correctly (max error: {discrepancy:.2e})")
    
    print("=== RECOMMENDED ACTIONS ===")
    print("1. Reduce time step by factor of 10")
    print("2. Check CFL condition: dt < dx/c")
    print("3. Verify PML parameters are not too aggressive")
    print("4. Consider using a more stable time integration scheme")

if __name__ == "__main__":
    detailed_energy_analysis()
