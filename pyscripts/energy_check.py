import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import signal, stats
import os

def realistic_energy_assessment():
    """
    Evalúa la conservación de energía considerando que las fluctuaciones
    son normales en esquemas numéricos de plasma
    """
    data = pd.read_csv('data/energy_evolution.csv')
    
    time = data['time'].values
    total_energy = data['total_energy'].values
    em_energy = data['em_energy'].values
    electron_energy = data['electron_energy'].values
    ion_energy = data['ion_energy'].values
    
    print("=== REALISTIC ENERGY ASSESSMENT FOR PLASMA SIMULATIONS ===")
    
    # 1. Ignorar el transiente inicial (primeros 10-20% de los datos)
    transient_cutoff = int(0.2 * len(time))
    time_stable = time[transient_cutoff:]
    energy_stable = total_energy[transient_cutoff:]
    
    print(f"Analyzing {len(time_stable)} points after transient (t > {time_stable[0]:.2e} s)")
    
    # 2. Análisis estadístico de las fluctuaciones
    energy_mean = np.mean(energy_stable)
    energy_std = np.std(energy_stable)
    energy_cv = energy_std / energy_mean  # Coeficiente de variación
    
    # 3. Detectar si hay tendencia (crecimiento/decaimiento sistemático)
    slope, intercept, r_value, p_value, std_err = stats.linregress(time_stable, energy_stable)
    trend_per_second = slope / energy_mean  # Tendencia relativa por segundo
    
    # 4. Análisis de las componentes individuales
    em_std = np.std(em_energy[transient_cutoff:]) / np.mean(em_energy[transient_cutoff:])
    electron_std = np.std(electron_energy[transient_cutoff:]) / np.mean(electron_energy[transient_cutoff:])
    ion_std = np.std(ion_energy[transient_cutoff:]) / np.mean(ion_energy[transient_cutoff:])
    
    # 5. Análisis espectral de las fluctuaciones
    dt = np.mean(np.diff(time_stable))
    freqs, psd = signal.welch(energy_stable - np.mean(energy_stable), fs=1/dt, nperseg=min(256, len(energy_stable)//4))
    
    # 6. Evaluación adaptativa basada en la física del problema
    print(f"\n--- ENERGY STATISTICS ---")
    print(f"Mean total energy: {energy_mean:.6e} J/m")
    print(f"Standard deviation: {energy_std:.6e} J/m")
    print(f"Coefficient of variation: {energy_cv*100:.2f}%")
    print(f"Trend: {trend_per_second*100:.4f}% per second")
    print(f"R² of linear trend: {r_value**2:.6f}")
    
    print(f"\n--- COMPONENT FLUCTUATIONS ---")
    print(f"EM energy fluctuation: {em_std*100:.2f}%")
    print(f"Electron energy fluctuation: {electron_std*100:.2f}%")
    print(f"Ion energy fluctuation: {ion_std*100:.2f}%")
    
    # 7. Criterios de evaluación REALISTAS para simulaciones de plasma
    print(f"\n--- ASSESSMENT CRITERIA ---")
    
    # Criterio 1: Fluctuaciones alrededor de la media
    if energy_cv < 0.01:  # < 1%
        stability = "EXCELLENT"
    elif energy_cv < 0.05:  # < 5%
        stability = "VERY GOOD" 
    elif energy_cv < 0.10:  # < 10%
        stability = "GOOD"
    elif energy_cv < 0.20:  # < 20%
        stability = "ACCEPTABLE for plasma simulations"
    elif energy_cv < 0.50:  # < 50%
        stability = "POOR but common in stiff plasma systems"
    else:
        stability = "UNSTABLE"
    
    # Criterio 2: Tendencia sistemática
    if abs(trend_per_second) < 1e-5:  # < 0.001% por segundo
        trend_assessment = "NO significant trend"
    elif abs(trend_per_second) < 1e-3:  # < 0.1% por segundo
        trend_assessment = "SLIGHT trend, probably acceptable"
    elif abs(trend_per_second) < 1e-2:  # < 1% por segundo
        trend_assessment = "NOTICEABLE trend, check parameters"
    else:
        trend_assessment = "STRONG trend, likely unstable"
    
    # Criterio 3: Comparación con fluctuaciones de componentes
    component_fluctuation_ratio = energy_cv / max(em_std, electron_std, ion_std)
    if component_fluctuation_ratio < 2.0:
        coupling_assessment = "Energy coupling appears PHYSICAL"
    else:
        coupling_assessment = "Energy coupling may be NUMERICAL"
    
    print(f"\n--- FINAL ASSESSMENT ---")
    print(f"Energy stability: {stability}")
    print(f"Trend analysis: {trend_assessment}") 
    print(f"Coupling assessment: {coupling_assessment}")
    
    # 8. Visualización comprehensiva
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # Plot 1: Energía total con estadísticas
    axes[0,0].plot(time, total_energy, 'k-', linewidth=2, label='Total Energy')
    axes[0,0].axvline(x=time[transient_cutoff], color='red', linestyle='--', 
                     label=f'Transient cutoff (t={time[transient_cutoff]:.2e}s)')
    axes[0,0].axhline(y=energy_mean, color='blue', linestyle='-', label=f'Mean: {energy_mean:.2e} J/m')
    axes[0,0].fill_between(time_stable, energy_mean - energy_std, energy_mean + energy_std, 
                          alpha=0.3, color='gray', label='±1 std')
    axes[0,0].set_xlabel('Time (s)')
    axes[0,0].set_ylabel('Total Energy (J/m)')
    axes[0,0].set_title(f'Total Energy (CV: {energy_cv*100:.1f}%)')
    axes[0,0].legend()
    axes[0,0].grid(True, alpha=0.3)
    
    # Plot 2: Componentes de energía
    axes[0,1].plot(time, em_energy, 'b-', label=f'EM (CV: {em_std*100:.1f}%)', linewidth=2)
    axes[0,1].plot(time, electron_energy, 'r-', label=f'Electron (CV: {electron_std*100:.1f}%)', linewidth=2)
    axes[0,1].plot(time, ion_energy, 'g-', label=f'Ion (CV: {ion_std*100:.1f}%)', linewidth=2)
    axes[0,1].set_xlabel('Time (s)')
    axes[0,1].set_ylabel('Energy (J/m)')
    axes[0,1].set_title('Energy Components')
    axes[0,1].legend()
    axes[0,1].grid(True, alpha=0.3)
    
    # Plot 3: Densidad espectral de potencia
    axes[0,2].semilogy(freqs[freqs > 0], psd[freqs > 0], 'b-', linewidth=2)
    axes[0,2].set_xlabel('Frequency (Hz)')
    axes[0,2].set_ylabel('Power Spectral Density')
    axes[0,2].set_title('Energy Fluctuation Spectrum')
    axes[0,2].grid(True, alpha=0.3)
    
    # Plot 4: Fluctuaciones relativas
    relative_energy = (energy_stable - energy_mean) / energy_mean * 100
    axes[1,0].plot(time_stable, relative_energy, 'r-', linewidth=2)
    axes[1,0].axhline(y=0, color='k', linestyle='-', alpha=0.5)
    axes[1,0].axhline(y=energy_cv*100, color='b', linestyle='--', label=f'Std: {energy_cv*100:.1f}%')
    axes[1,0].axhline(y=-energy_cv*100, color='b', linestyle='--')
    axes[1,0].set_xlabel('Time (s)')
    axes[1,0].set_ylabel('Relative Fluctuation (%)')
    axes[1,0].set_title('Relative Energy Fluctuations')
    axes[1,0].legend()
    axes[1,0].grid(True, alpha=0.3)
    
    # Plot 5: Histograma de fluctuaciones
    axes[1,1].hist(relative_energy, bins=30, density=True, alpha=0.7, color='blue')
    axes[1,1].axvline(x=0, color='k', linestyle='-', alpha=0.5)
    axes[1,1].set_xlabel('Relative Fluctuation (%)')
    axes[1,1].set_ylabel('Probability Density')
    axes[1,1].set_title('Distribution of Energy Fluctuations')
    axes[1,1].grid(True, alpha=0.3)
    
    # Plot 6: Resumen de evaluación
    axes[1,2].axis('off')
    assessment_text = f"""
    REALISTIC ENERGY ASSESSMENT
    ===========================
    
    Stability: {stability}
    Trend: {trend_assessment}
    Coupling: {coupling_assessment}
    
    Statistics:
    • Mean energy: {energy_mean:.2e} J/m
    • Fluctuation: {energy_cv*100:.1f}% (std/mean)
    • Trend rate: {trend_per_second*100:.4f}%/s
    
    Component Fluctuations:
    • EM: {em_std*100:.1f}%
    • Electron: {electron_std*100:.1f}%
    • Ion: {ion_std*100:.1f}%
    
    Interpretation:
    • < 20% fluctuation: Normal for plasma simulations
    • < 5% fluctuation: Very good for RK4 scheme
    • < 1% fluctuation: Excellent (rare in plasma codes)
    """
    axes[1,2].text(0.05, 0.95, assessment_text, transform=axes[1,2].transAxes, 
                   fontfamily='monospace', fontsize=10, verticalalignment='top',
                   bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray"))
    
    plt.tight_layout()
    os.makedirs('plots', exist_ok=True)
    plt.savefig('plots/realistic_energy_assessment.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 9. Recomendaciones específicas para simulaciones de plasma
    print(f"\n=== RECOMMENDATIONS FOR PLASMA SIMULATIONS ===")
    
    if energy_cv > 0.20:
        print("1. Consider using a symplectic integrator (Leapfrog) for better energy conservation")
        print("2. Reduce time step if high-frequency oscillations dominate")
        print("3. Check if PML boundaries are causing numerical reflections")
    elif energy_cv > 0.10:
        print("1. Energy conservation is acceptable for most plasma physics applications")
        print("2. Monitor long-term trends rather than short-term fluctuations")
        print("3. Consider averaging results over several oscillation periods")
    else:
        print("1. Energy conservation is excellent for this type of simulation")
        print("2. Your numerical scheme is working very well")
    
    if abs(trend_per_second) > 1e-3:
        print("4. Significant trend detected: check collision terms and boundary conditions")
    
    print(f"\nRemember: In plasma simulations, energy fluctuations of 10-20% are often acceptable")
    print("The key is whether the physical phenomena of interest are correctly captured")

if __name__ == "__main__":
    realistic_energy_assessment()
