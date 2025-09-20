# check_fields.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def check_magnetic_fields(mode):
    try:
        data = pd.read_csv(f'data/field_data_{mode}.csv')
        
        print(f"\n=== ANALYSIS OF MODE {mode} ===")
        
        # Campos eléctricos
        print(f"Range of Ex: [{data['Ex'].min():.3e}, {data['Ex'].max():.3e}]")
        print(f"Range of Ey: [{data['Ey'].min():.3e}, {data['Ey'].max():.3e}]")
        print(f"Range of Ez: [{data['Ez'].min():.3e}, {data['Ez'].max():.3e}]")
        
        # Campos magnéticos
        print(f"Range of Bx: [{data['Bx'].min():.3e}, {data['Bx'].max():.3e}]")
        print(f"Range of By: [{data['By'].min():.3e}, {data['By'].max():.3e}]")
        print(f"Range of Bz: [{data['Bz'].min():.3e}, {data['Bz'].max():.3e}]")
        
        # Corrientes electrónicas
        print(f"Range of Jx_e: [{data['Jx_e'].min():.3e}, {data['Jx_e'].max():.3e}]")
        print(f"Range of Jy_e: [{data['Jy_e'].min():.3e}, {data['Jy_e'].max():.3e}]")
        print(f"Range of Jz_e: [{data['Jz_e'].min():.3e}, {data['Jz_e'].max():.3e}]")
        
        # Corrientes iónicas
        print(f"Range of Jx_i: [{data['Jx_i'].min():.3e}, {data['Jx_i'].max():.3e}]")
        print(f"Range of Jy_i: [{data['Jy_i'].min():.3e}, {data['Jy_i'].max():.3e}]")
        print(f"Range of Jz_i: [{data['Jz_i'].min():.3e}, {data['Jz_i'].max():.3e}]")
        
        # Verificar si los campos magnéticos son esencialmente cero
        bx_max = np.max(np.abs(data['Bx']))
        by_max = np.max(np.abs(data['By']))
        bz_max = np.max(np.abs(data['Bz']))
        
        print(f"Absolute maximum value - Bx: {bx_max:.3e}, By: {by_max:.3e}, Bz: {bz_max:.3e}")
        
        if bx_max < 1e-15 and by_max < 1e-15:
            print("❌ PROBLEM: Magnetic fields essentially zero")
        else:
            print("✅ Magnetic fields seem correct")
            
        # Comparar magnitudes de corrientes electrónicas e iónicas
        j_e_max = np.max(np.abs(data[['Jx_e', 'Jy_e', 'Jz_e']].values))
        j_i_max = np.max(np.abs(data[['Jx_i', 'Jy_i', 'Jz_i']].values))
        
        print(f"Max electron current: {j_e_max:.3e}")
        print(f"Max ion current: {j_i_max:.3e}")
        print(f"Ratio (J_e/J_i): {j_e_max/j_i_max:.3f}")
        
        # Graficar para inspección visual
        plt.figure(figsize=(15, 10))
        
        # Campos eléctricos
        plt.subplot(2, 2, 1)
        plt.plot(data['z'], data['Ex'], 'b-', label='Ex')
        plt.plot(data['z'], data['Ey'], 'r-', label='Ey')
        plt.plot(data['z'], data['Ez'], 'g-', label='Ez')
        plt.xlabel('Position (m)')
        plt.ylabel('Electric Field (V/m)')
        plt.legend()
        plt.title(f'Electric Fields - Mode {mode}')
        
        # Campos magnéticos
        plt.subplot(2, 2, 2)
        plt.plot(data['z'], data['Bx'], 'c-', label='Bx')
        plt.plot(data['z'], data['By'], 'm-', label='By')
        plt.plot(data['z'], data['Bz'], 'y-', label='Bz')
        plt.xlabel('Position (m)')
        plt.ylabel('Magnetic Field (T)')
        plt.legend()
        plt.title(f'Magnetic Fields - Mode {mode}')
        
        # Corrientes electrónicas
        plt.subplot(2, 2, 3)
        plt.plot(data['z'], data['Jx_e'], 'b-', label='Jx_e')
        plt.plot(data['z'], data['Jy_e'], 'r-', label='Jy_e')
        plt.plot(data['z'], data['Jz_e'], 'g-', label='Jz_e')
        plt.xlabel('Position (m)')
        plt.ylabel('Electron Current Density (A/m²)')
        plt.legend()
        plt.title(f'Electron Currents - Mode {mode}')
        
        # Corrientes iónicas
        plt.subplot(2, 2, 4)
        plt.plot(data['z'], data['Jx_i'], 'b-', label='Jx_i')
        plt.plot(data['z'], data['Jy_i'], 'r-', label='Jy_i')
        plt.plot(data['z'], data['Jz_i'], 'g-', label='Jz_i')
        plt.xlabel('Position (m)')
        plt.ylabel('Ion Current Density (A/m²)')
        plt.legend()
        plt.title(f'Ion Currents - Mode {mode}')
        
        plt.tight_layout()
        os.makedirs('plots', exist_ok=True)
        plt.savefig(f'plots/field_check_{mode}.png', dpi=300)
        plt.show()
        
    except FileNotFoundError:
        print(f"File for mode {mode} not found")
    except KeyError as e:
        print(f"Column {e} not found in the data file. Check your CSV format.")

if __name__ == "__main__":
    for mode in ['R', 'L', 'O', 'X']:
        check_magnetic_fields(mode)
