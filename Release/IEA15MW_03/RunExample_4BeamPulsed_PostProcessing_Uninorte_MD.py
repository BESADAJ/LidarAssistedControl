# IEA15MW_03: IEA 15 MW monopile + realistic wind preview from a
# 4-beam pulsed lidar system measuring at 160 m.
# This script needs to be run after RunExample_4BeamPulsed.py.
# Purpose:
# A postprocessing version without the need to compile DLLs for lidar data
# processing to be used in the LAC Summer Games 2024.
# To implement your own solution, replace line 60
# R_FBFF = CalculateREWSfromLidarData_LDP_v1(FBFF, DT, TMax, LDP)
# with your own function with the same inputs and outputs.
# Result:
# Cost for Summer Games 2024 ("18 m/s hurdles"):  0.444617 m/s

import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.interpolate import interp1d
from CalculateREWSfromLidarData_LDP_MD import CalculateREWSfromLidarData_LDP_MD

sys.path.append('../PythonFunctions')
from PreProcessing.CalculateREWSfromWindField import CalulateREWSfromWindField
from FileOperations.ReadFASTbinaryIntoStruct import ReadFASTbinaryIntoStruct

# Seeds (can be adjusted, but will provide different results)
nSeed = 1
Seed_vec = np.arange(1, nSeed + 1) + 18 * 100

# Parameters postprocessing (can be adjusted, but will provide different results)
t_start = 60                                        # [s] 	ignore data before for STD and spectra
TMax = 660                                          # [s]   total run time, same as in *.fst
DT = 0.0125                                         # [s]   time step, same as in *.fst
R = 120                                             # [m]  	rotor radius to calculate REWS

# Parameter for Cost (Summer Games 2024)
tau = 2                                             # [s]   time to overcome pitch actuator, from Example 1: tau = T_Taylor - T_buffer, since there T_filter = T_scan = 0

#Parameters modification for brute force, T_Buffer, f_cutoff

T_buffer_start = 1.8
T_buffer_reset = T_buffer_start
T_buffer_step = 0.1
T_buffer_count = 1
f_cutoff_start = 0.16
f_cutoff_step = 0.01
f_cutoff_count = 1
CostPlot = np.zeros(T_buffer_count)
BufferPlot = np.zeros(T_buffer_count)

# Configuration from LDP_v1_4BeamPulsed.IN and LDP_v1_4BeamPulsed.IN
for i_f_cutoff in range(f_cutoff_count):
    T_buffer_start=T_buffer_reset
    for iBuffer in range(T_buffer_count):
        LDP = {
            "NumberOfBeams": 4,
            "AngleToCenterline": 19.176,
            "IndexGate": 6,
            "FlagLPF": 1,
            "f_cutoff":  float(f_cutoff_start),
            "T_buffer": float(T_buffer_start),
        }


        # Files (should not be changed)
        SimulationFolderLAC = "SimulationResults_4BeamPulsed"

        # Allocate
        MAE = np.empty(nSeed)

        # Loop over all seeds

        for iSeed in range(nSeed):
            # Load data
            Seed = Seed_vec[iSeed]
            WindFileName = f'URef_18_Seed_{Seed:02d}'
            FASTresultFile = f'{SimulationFolderLAC}/{WindFileName}_FlagLAC_1.outb'
            FBFF = ReadFASTbinaryIntoStruct(FASTresultFile)

            # Calculate REWS
            R_FBFF = CalculateREWSfromLidarData_LDP_MD(FBFF, DT, TMax, LDP)

            # Get REWS from the wind field and interpolate it on the same time vector
            TurbSimResultFile = f'TurbulentWind/URef_18_Seed_{Seed:02d}.wnd'
            REWS_WindField, Time_WindField = CalulateREWSfromWindField(TurbSimResultFile, R, 2)
            REWS_WindField_Fs = interp1d(Time_WindField.ravel(), REWS_WindField.ravel())(R_FBFF['Time'])

            # Calculate mean absolute error
            REWS_WindField_Fs_shifted = interp1d(Time_WindField.ravel() - tau, REWS_WindField.ravel())(R_FBFF['Time'])
            Error = REWS_WindField_Fs_shifted - R_FBFF["REWS_b"]
            MAE[iSeed] = np.mean(np.abs(Error[R_FBFF["Time"] >= t_start]))

            # Plot REWS for absolute error
            plt.figure('REWS seed {}'.format(Seed))
            plt.subplot(311)
            plt.plot(R_FBFF['Time'], REWS_WindField_Fs, label='wind field')
            plt.plot(R_FBFF['Time'], R_FBFF['REWS'], label='lidar estimate')
            plt.ylabel('REWS [m/s]')
            plt.legend()
            plt.grid(True)
            plt.subplot(312)
            plt.plot(R_FBFF['Time'], REWS_WindField_Fs_shifted, label='wind field shifted')
            plt.plot(R_FBFF['Time'], R_FBFF['REWS_b'], label='lidar estimate filtered and buffered')
            plt.ylabel('REWS [m/s]')
            plt.legend()
            plt.grid(True)
            plt.subplot(313)
            plt.plot(R_FBFF['Time'], Error)
            plt.ylabel('error [m/s]')
            plt.xlabel('time [s]')
            plt.grid(True)
        BufferPlot[iBuffer] = T_buffer_start
        T_buffer_start = T_buffer_start + T_buffer_step
        Cost = np.mean(MAE)




        print(f'Cost for Summer Games 2024 ("18 m/s hurdles"): {Cost:.6f}')

        CostPlot[iBuffer]= Cost



    plt.plot(BufferPlot,CostPlot,label=f_cutoff_start )
    plt.legend(title= "f_cutoff [Hz]")
    f_cutoff_start = f_cutoff_start + f_cutoff_step


plt.ylabel('Cost')
plt.xlabel('T_Buffer')
plt.show()


# Calculation of Cost for Summer Games 2024
