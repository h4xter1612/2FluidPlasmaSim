# check_fields.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def check_magnetic_fields(mode):
    try:
        data = pd.read_csv(f'data/field_data_{mode}.csv')
        
        print(f"\n=== ANALYSIS OF MODE {mode} ===")
        print(f"Range of Bx: [{data['Bx'].min():.3e}, {data['Bx'].max():.3e}]")
        print(f"Range of By: [{data['By'].min():.3e}, {data['By'].max():.3e}]")
        print(f"Range of Bz: [{data['Bz'].min():.3e}, {data['Bz'].max():.3e}]")
        
        # Verificar si los campos magnéticos son esencialmente cero
        bx_max = np.max(np.abs(data['Bx']))
        by_max = np.max(np.abs(data['By']))
        bz_max = np.max(np.abs(data['Bz']))
        
        print(f"Absolute maximum value - Bx: {bx_max:.3e}, By: {by_max:.3e}, Bz: {bz_max:.3e}")
        
        if bx_max < 1e-15 and by_max < 1e-15:
            print("❌ PROBLEM: Magnetic fields essentially zero")
        else:
            print("✅ Magnetic fields seem correct")
            
        # Graficar para inspección visual
        plt.figure(figsize=(12, 6))
        plt.subplot(1, 2, 1)
        plt.plot(data['z'], data['Ex'], 'b-', label='Ex')
        plt.plot(data['z'], data['Ey'], 'r-', label='Ey')
        plt.xlabel('Position (m)')
        plt.ylabel('Electric Field (V/m)')
        plt.legend()
        plt.title(f'Electric Fields - Mode {mode}')
        
        plt.subplot(1, 2, 2)
        plt.plot(data['z'], data['Bx'], 'c-', label='Bx')
        plt.plot(data['z'], data['By'], 'm-', label='By')
        plt.xlabel('Position (m)')
        plt.ylabel('Magnetic Field (T)')
        plt.legend()
        plt.title(f'Magnetic Fields - Mode {mode}')
        
        plt.tight_layout()
        os.makedirs('plots', exist_ok=True)
        plt.savefig(f'plots/field_check_{mode}.png')
        plt.show()
        
    except FileNotFoundError:
        print(f"File for mode {mode} not found")

if __name__ == "__main__":
    for mode in ['R', 'L', 'O', 'X']:
        check_magnetic_fields(mode)
