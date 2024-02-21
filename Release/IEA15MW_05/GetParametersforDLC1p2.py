import numpy as np
import sys

sys.path.append('..\PythonFunctions')
from GetStatistics import GetStatistics
from CalculateDEL import CalculateDEL
def GetParametersforDLC1p2(simulation_mode):
    assert simulation_mode in ['FeedbackOnly', 'LAC_CircularCW', 'LAC_4BeamPulsed'], "Invalid SimulationMode"

    # Process parameters
    start_time = 60  # [s] time to start evaluation (all signals should be settled)

    # files
    statistics_file = 'Statistics_SteadyStates.csv'

    # Variation
    pre_processing_variation = [['URef', np.arange(4, 21, 4), '%02d'],
                                ['Seed', np.arange(1, 3), '%02d']]

    # template files
    input_files = [['IEA-15-240-RWT-Monopile.fst'],  # main file
                   ['IEA-15-240-RWT-Monopile_ElastoDyn.dat'],  # to adjust initial conditions
                   ['IEA-15-240-RWT_InflowFile.dat']]  # to adjust wind speed

    if simulation_mode == 'LAC_CircularCW':
        input_files.extend([
            ['LidarFile_CircularCW.dat'],  # to adjust wind speed
            ['IEA-15-240-RWT-Monopile_ServoDyn.dat'],  # to link to new WRAPPER_CircularCW.IN
            ['WRAPPER_CircularCW.IN'],  # to link to new FFP_v1_CircularCW.IN and ROSCO_v2d6.IN
            ['ROSCO_v2d6.IN'],  # to enable FF
            ['FFP_v1_CircularCW.IN']  # to adjust f_cutoff and T_buffer
        ])
    elif simulation_mode == 'LAC_4BeamPulsed':
        input_files.extend([
            ['LidarFile_4BeamPulsed.dat'],  # to adjust wind speed
            ['IEA-15-240-RWT-Monopile_ServoDyn.dat'],  # to link to new WRAPPER_4BeamPulsed.IN
            ['WRAPPER_4BeamPulsed.IN'],  # to link to new FFP_v1_4BeamPulsed.IN and ROSCO_v2d6.IN
            ['ROSCO_v2d6.IN'],  # to enable FF
            ['FFP_v1_4BeamPulsed.IN']  # to adjust f_cutoff and T_buffer
        ])

        # new files to be modified
    input_files.extend([
        ['<SimulationName>.fst'],
        ['<SimulationName>_ElastoDyn.dat'],
        ['<SimulationName>_InflowWind.dat']
    ])

    if simulation_mode in ['LAC_CircularCW', 'LAC_4BeamPulsed']:
        input_files.extend([
            ['<SimulationName>_LidarFile.dat'],
            ['<SimulationName>_ServoDyn.dat'],
            ['<SimulationName>_WRAPPER.IN'],
            ['ROSCO_v2d6_LAC.IN'],
            ['<SimulationName>_FFP_v1.IN']
        ])

    # modifications in new files (for all Modes)
    modifications = [
        # OpenFAST: change time and link to new files
        ['1', 'I', 'TMax', str(600 + start_time)],
        ['1', 'I', 'EDFile', input_files[1]],
        ['1', 'I', 'InflowFile', input_files[2]],
        # ElastoDyn: adjust initial conditions
        ['2', 'I', 'BlPitch(1|2|3)',
         lambda variation_values: '%5.2f' % GetStatistics(statistics_file, 'mean_BldPitch1', variation_values[0])],
        ['2', 'I', 'RotSpeed',
         lambda variation_values: '%5.2f' % GetStatistics(statistics_file, 'mean_RotSpeed', variation_values[0])],
        # InflowWind: change wind type and link to turbulent wind file
        ['3', 'I', 'WindType', '4'],
        ['3', 'I', 'FilenameRoot', lambda variation_values: '../TurbulentWind/URef_%02d_Seed_%02d%02d' % (
        variation_values[0], variation_values[0], variation_values[1])]
    ]
    # FFP_v1 Parameters
    u_ref_v = np.arange(10, 25, 2)  # [m/s]
    if simulation_mode == 'LAC_CircularCW':
        f_cutoff_v = [0.1816, 0.2179, 0.2542, 0.2905, 0.3268, 0.3632, 0.3995, 0.4358]  # [rad/s]
        t_buffer_v = [13.0000, 10.4167, 8.5714, 7.1875, 6.1111, 5.2500, 4.5455, 3.9583]  # [s]
    elif simulation_mode == 'LAC_4BeamPulsed':
        f_cutoff_v = [0.0684, 0.0821, 0.0958, 0.1095, 0.1232, 0.1369, 0.1505, 0.1642]  # [rad/s]
        t_buffer_v = [4.5000, 3.3333, 2.5000, 1.8750, 1.3889, 1.0000, 0.6818, 0.4167]  # [s]

    # additional modifications for LAC modes
    if simulation_mode in ['LAC_CircularCW', 'LAC_4BeamPulsed']:
        modifications.extend([
            # OpenFAST: link to new files
            ['1', 'I', 'SWELidarFile', input_files[3][1]],
            ['1', 'I', 'ServoFile', input_files[4][1]],
            # LidarFile: adjust wind speed
            ['4', 'I', 'URef', lambda variation_values: str(variation_values[0])],
            # ServoDyn: link to Wrapper
            ['5', 'I', 'DLL_FileName', 'WRAPPER.dll'],
            ['5', 'I', 'DLL_InFile', input_files[5][1]],
            # Wrapper: link to new ROSCO_v2d6.IN
            ['6', 'R', input_files[6][0], input_files[6][1]],
            ['6', 'R', input_files[7][0], input_files[7][1]],
            # Rosco: enable FF
            ['7', 'I', '! FlagLAC', '1'],
            # FFP_v1: adjust f_cutoff and T_buffer
            ['8', 'I', '! f_cutoff',
             lambda variation_values: '%8.4f' % np.interp(variation_values[0], u_ref_v, f_cutoff_v)],
            ['8', 'I', '! T_buffer',
             lambda variation_values: '%8.4f' % np.interp(variation_values[0], u_ref_v, t_buffer_v)]
        ])

    # PlotTimeResults
    n_u_ref = len(pre_processing_variation[0][1])
    n_seed = len(pre_processing_variation[1][1])
    post_processing_config = {'Plots': {'BasicTimePlot': {}}}
    for i_u_ref in range(n_u_ref):
        id = i_u_ref
        post_processing_config['Plots']['BasicTimePlot'][id] = {
            'Enable': 1,
            'Channels': ['Wind1VelX', 'BldPitch1', 'RotSpeed'],
            'gca': {'xlim': [0, 600, start_time]},
            'IndicesConsideredDataFiles': list(range(1, n_seed + 1)) + [(i_u_ref - 1) * n_seed]
        }

    # CalculateStatistics
    woehler_exponent = 4
    post_processing_config['CalculateStatistics'] = {
        'mean': [lambda data, time: np.mean(data[time >= start_time]), ['Wind1VelX', 'GenPwr']],
        'DEL': [lambda data, time: CalculateDEL(data[time >= start_time], time[time >= start_time], woehler_exponent),
                ['TwrBsMyt']]
    }

    return  pre_processing_variation, input_files, modifications, post_processing_config
