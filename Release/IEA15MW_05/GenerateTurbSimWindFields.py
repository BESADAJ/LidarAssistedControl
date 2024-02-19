import os
import numpy as np
import shutil
import sys

sys.path.append('..\PythonFunctions')

from ManipulateFastInputFile import ManipulateFastInputFile

# Parameters (can be adjusted, but will provide different results)
URef_vector = np.arange(4, 25, 2)  # [m/s] range of wind speeds (operation points)
n_Seed = 6  # [-] number of stochastic turbulence field seeds

# Seed Matrix Definition
Seed_vector = np.arange(1, n_Seed + 1)
n_URef = len(URef_vector)
Seed_matrix = np.outer(URef_vector, np.ones(n_Seed)) * 100 + np.outer(np.ones(n_URef), Seed_vector)

# Files (should not be changed)
TurbSimExeFile = 'TurbSim_x64.exe'
TurbSimTemplateFile = 'TurbSimInputFileTemplateIEA15MW.inp'

# Generate folder
if not os.path.exists('TurbulentWind'):
    os.mkdir('TurbulentWind')

# Preprocessing: generate turbulent wind field

# Copy the adequate TurbSim version to the example folder
shutil.copyfile(os.path.join('..', 'TurbSim', TurbSimExeFile), os.path.join('TurbulentWind', TurbSimExeFile))

# Generate all wind fields for different URef and RandSeed1
for i_URef in range(n_URef):
    URef = URef_vector[i_URef]
    for i_Seed in range(n_Seed):
        Seed = Seed_matrix[i_URef, i_Seed]
        WindFileName = f'URef_{int(URef):02d}_Seed_{int(Seed):04d}'
        TurbSimInputFile = os.path.join('TurbulentWind', f'{WindFileName}.ipt')
        TurbSimResultFile = os.path.join('TurbulentWind', f'{WindFileName}.wnd')
        if not os.path.exists(TurbSimResultFile):
            shutil.copyfile(TurbSimTemplateFile, TurbSimInputFile)
            # Adjust the TurbSim input file
            ManipulateFastInputFile(TurbSimInputFile, 'URef ', f'{URef:4.1f}')  # adjust URef
            ManipulateFastInputFile(TurbSimInputFile, 'RandSeed1 ', str(Seed))  # adjust seed
            # Generate wind field
            os.system(os.path.join('TurbulentWind', TurbSimExeFile) + ' ' + TurbSimInputFile)

# Clean up
os.remove(os.path.join('TurbulentWind', TurbSimExeFile))
