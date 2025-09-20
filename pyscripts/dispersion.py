import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

def plot_dispersion_relations(data):
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
    plt.savefig('plots/dispersion_relations.png', dpi=300)
    plt.show()

def analyze_cutoffs_resonances(data):
    # Encontrar cutoffs (donde k ≈ 0)
    cutoff_R = data['frequency'][np.argmin(np.abs(data['Re(k_R)']))]
    cutoff_L = data['frequency'][np.argmin(np.abs(data['Re(k_L)']))]
    
    print(f"Cutoff modo R: {cutoff_R:.2e} Hz")
    print(f"Cutoff modo L: {cutoff_L:.2e} Hz")

def infer_L_from_fieldfile(mode, folder='data'):
    """
    Read data/field_data_{mode}.csv and return L = z_max - z_min (m).
    """
    path = os.path.join(folder, f'field_data_{mode}.csv')
    if not os.path.isfile(path):
        raise FileNotFoundError(f"{path} not found.")
    df = pd.read_csv(path)
    if 'z' not in df.columns:
        raise ValueError(f"Column 'z' not found in {path}.")
    zmin, zmax = float(df['z'].min()), float(df['z'].max())
    return zmax - zmin, path

def faraday_from_dispersion_with_inferred_L(data,
                                            mode_for_L='R', folder='data'):
    """
    Compute Faraday rotation from dispersion CSV, inferring L from field_data_{mode_for_L}.csv.
    Returns DataFrame with frequency, theta_deg, ellipticity and also rot_per_m (deg/m).
    """
    disp = data
    L, src = infer_L_from_fieldfile(mode_for_L, folder=folder)
    # print(f"Inferred L = {L:.6g} m from {src}")
    kR = disp['Re(k_R)'].to_numpy(dtype=float) + 1j * disp['Im(k_R)'].to_numpy(dtype=float)
    kL = disp['Re(k_L)'].to_numpy(dtype=float) + 1j * disp['Im(k_L)'].to_numpy(dtype=float)
    dk = kR - kL
    theta_rad = 0.5 * np.real(dk) * L          # rad
    theta_deg = theta_rad * 180.0 / np.pi     # deg
    ellipticity = 0.5 * np.imag(dk) * L
    # rotation per meter (deg/m) from dispersion: 0.5 * Re(dk) in rad/m converted to deg/m
    rot_per_m_deg = 0.5 * np.real(dk) * (180.0 / np.pi)
    results = pd.DataFrame({
        'frequency': disp['frequency'],
        'theta_deg': theta_deg,
        'ellipticity': ellipticity,
        'rot_per_m_deg': rot_per_m_deg
    })
    results.to_csv('data/faraday_from_dispersion_inferred_L.csv', index=False)
    return results, L

def compute_stokes_angle_from_fieldfile(mode, folder='data'):
    """
    Compute averaged Stokes-derived polarization angle psi (rad) from a single field file.
    Uses columns Ex, Ey (real). Returns psi (rad) and psi_deg.
    """
    path = os.path.join(folder, f'field_data_{mode}.csv')
    df = pd.read_csv(path)
    if not {'Ex','Ey'}.issubset(df.columns):
        raise ValueError(f"Columns Ex and Ey required in {path}.")
    # drop nan rows
    df = df.dropna(subset=['Ex','Ey'])
    Ex = df['Ex'].to_numpy(dtype=float)
    Ey = df['Ey'].to_numpy(dtype=float)
    # averaged Stokes-like quantities
    S1 = np.mean(Ex**2 - Ey**2)
    S2 = 2.0 * np.mean(Ex * Ey)
    psi = 0.5 * np.arctan2(S2, S1)  # radians
    return {'psi_rad': psi, 'psi_deg': psi * 180.0 / np.pi, 'npoints': len(df)}

def rotation_from_fields_vs_dispersion(mode='R', dispersion_csv='data/dispersion_data.csv', folder='data'):
    """
    Compare field-measured rotation per meter with dispersion-predicted rotation per meter.
    - Field-measured rot_per_m: compute psi at two reference files or use single file (psi_total / L).
      Here we assume a single file containing the full z-range; rotation per meter = psi_deg / L.
    - Dispersion-predicted rot_per_m: from 0.5*Re(kR-kL) converted to deg/m (independent of L).
    """
    # 1) infer L from the field file
    L, src = infer_L_from_fieldfile(mode, folder=folder)
    # print(f"Inferred L = {L:.6g} m from {src}")
    # 2) compute psi from fields (averaged over z)
    st = compute_stokes_angle_from_fieldfile(mode, folder=folder)
    psi_deg = st['psi_deg']
    psi_rad = st['psi_rad']
    # Field-measured rotation per meter (deg/m)
    field_rot_per_m = psi_deg / L
    # 3) dispersion prediction (rot per meter in deg/m) - independent of L
    disp = pd.read_csv(dispersion_csv)
    kR = disp['Re(k_R)'].to_numpy(dtype=float) + 1j * disp['Im(k_R)'].to_numpy(dtype=float)
    kL = disp['Re(k_L)'].to_numpy(dtype=float) + 1j * disp['Im(k_L)'].to_numpy(dtype=float)
    rot_per_m_deg_disp = 0.5 * np.real(kR - kL) * (180.0 / np.pi)  # deg per meter, array over frequency
    # Save small summary
    summary = {
        'mode': mode,
        'L_m': L,
        'psi_deg_avg': psi_deg,
        'field_rot_per_m_deg': field_rot_per_m
    }
    # print("Field-averaged psi = {:.4f} deg (n={} points)".format(psi_deg, st['npoints']))
    # print("Field-measured rotation per meter (deg/m) = {:.4e}".format(field_rot_per_m))
    # Optionally save rot_per_m dispersion vs frequency
    out = pd.DataFrame({
        'frequency': disp['frequency'],
        'rot_per_m_deg_disp': rot_per_m_deg_disp
    })
    out.to_csv(f'data/rot_per_m_dispersion_{mode}.csv', index=False)
    # Plot comparison: single horizontal line (field) vs dispersion curve (freq)
    plt.figure(figsize=(8,4))
    plt.plot(out['frequency'], out['rot_per_m_deg_disp'], label='Dispersion predicted (deg/m)')
    plt.axhline(field_rot_per_m, color='k', linestyle='--', label='Field-measured (deg/m)')
    plt.xscale('log')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Rotation (deg/m)')
    plt.legend()
    plt.grid(True)
    plt.title(f'Rotation per meter: mode {mode} (L={L:.3g} m)')
    plt.tight_layout()
    plt.savefig(f'plots/rot_per_m_comparison_{mode}.png', dpi=300)
    plt.show()
    return summary, out

def FR_print_results(summary, L, src):
    print("\n" + "="*60)
    print("            FARADAY ROTATION ANALYSIS REPORT")
    print("="*60)
    print(f"Field data source     : {src}")
    print(f"Inferred length (L)   : {L:.4f} m")
    print(f"Mode                  : {summary['mode']}")
    print("-"*60)
    print(f"Field-averaged ψ      : {summary['psi_deg_avg']:.3f} deg")
    print(f"Field rotation / m    : {summary['field_rot_per_m_deg']:.3f} deg/m")
    print("-"*60)
    print("Dictionary summary object:")
    print(summary)
    print("="*60 + "\n")

def plot_current_ratio_analysis():
    """
    Analizar la relación entre corrientes electrónicas e iónicas para todos los modos
    """
    modes = ['R', 'L', 'O', 'X']
    ratios = []
    
    plt.figure(figsize=(10, 6))
    
    for i, mode in enumerate(modes):
        try:
            data = pd.read_csv(f'data/field_data_{mode}.csv')
            
            # Calcular la relación máxima entre corrientes
            j_e_max = np.max(np.abs(data[['Jx_e', 'Jy_e', 'Jz_e']].values))
            j_i_max = np.max(np.abs(data[['Jx_i', 'Jy_i', 'Jz_i']].values))
            ratio = j_e_max / j_i_max if j_i_max > 0 else float('inf')
            ratios.append(ratio)
            
            # Graficar la relación a lo largo de z
            jx_ratio = np.abs(data['Jx_e']) / np.abs(data['Jx_i'] + 1e-30)  # Evitar división por cero
            plt.plot(data['z'], jx_ratio, label=f'Mode {mode} (avg: {np.mean(jx_ratio):.1f})')
            
        except FileNotFoundError:
            print(f"File for mode {mode} not found")
            ratios.append(0)
    
    plt.xlabel('Position (m)')
    plt.ylabel('|J_e| / |J_i| Ratio')
    plt.title('Electron to Ion Current Ratio')
    plt.legend()
    plt.grid(True)
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig('plots/current_ratio_analysis.png', dpi=300)
    plt.show()
    
    # Imprimir resumen
    print("\nCurrent Ratio Summary:")
    for mode, ratio in zip(modes, ratios):
        print(f"Mode {mode}: {ratio:.1f}")

if __name__ == "__main__":
    os.makedirs('animations', exist_ok=True)
    os.makedirs('plots', exist_ok=True)
    
    # Leer datos de dispersión
    data = pd.read_csv('data/dispersion_data.csv')
    
    # Predict Faraday using inferred L
    res_disp, L = faraday_from_dispersion_with_inferred_L(
        data, mode_for_L='R', folder='data'
    )
    
    # Compare measured vs predicted
    summary, disp_vs_freq = rotation_from_fields_vs_dispersion(
        mode='R', dispersion_csv='data/dispersion_data.csv', folder='data'
    )
    
    # Output
    _, src = infer_L_from_fieldfile('R', folder='data')
    FR_print_results(summary, L, src)
    
    # Análisis de relaciones de dispersión
    plot_dispersion_relations(data)
    analyze_cutoffs_resonances(data)
    
    # Análisis de relación de corrientes
    plot_current_ratio_analysis()
