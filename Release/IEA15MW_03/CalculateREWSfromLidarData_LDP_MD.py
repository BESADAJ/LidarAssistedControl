import numpy as np

    #multi-distance data processing function
def CalculateREWSfromLidarData_LDP_MD(FBFF, DT, TMax, LDP):
    # Function to postprocess lidar data and calculate rotor-effective wind speed (REWS)
    # similar to LDP_v1/FFP_v1 without the need for compiling a DLL.

    # Initialize time
    R_FBFF = {}
    R_FBFF["Time"] = np.arange(0, TMax + DT, DT)
    n_t = len(R_FBFF["Time"])

    # Allocate arrays
    R_FBFF["REWS"] = np.full(n_t, np.nan)
    R_FBFF["REWS_f"] = np.full(n_t, np.nan)
    R_FBFF["REWS_b"] = np.full(n_t, np.nan)
    
    # Loop over time
    for i_t in range(n_t):
        Idx = max(i_t - 1, 0)  # Due to co-simulation of OpenFAST and the controller DLLs

        # If there is a new measurement, perform wind field reconstruction
        if FBFF["NEWDATALI"][Idx]:
            v_los_MD = np.array([FBFF["VLOS01LI"][Idx],FBFF["VLOS02LI"][Idx],FBFF["VLOS03LI"][Idx],FBFF["VLOS04LI"][Idx],FBFF["VLOS05LI"][Idx],FBFF["VLOS06LI"][Idx],FBFF["VLOS07LI"][Idx],FBFF["VLOS08LI"][Idx],FBFF["VLOS09LI"][Idx],FBFF["VLOS10LI"][Idx]])
            REWS_MD = WindFieldReconstruction(v_los_MD, LDP["NumberOfBeams"], LDP["AngleToCenterline"])

            #combine to one distance
            [REWS, REWS_MD_shifted] = CombineAndShift(REWS_MD, LDP['MeasurementDistances'], LDP['IdxFirstDistance'], LDP["MeanWindSpeed"], DT);


            # Low-pass filter the REWS
        if LDP["FlagLPF"]:
            REWS_f = LPFilter(REWS, DT, LDP["f_cutoff"])
        else:
            REWS_f = REWS

        # Get buffered and filtered REWS from buffer
        REWS_b = Buffer(REWS_f, DT, LDP["T_buffer"])

        # Store in structure
        R_FBFF["REWS"][i_t] = REWS
        R_FBFF["REWS_f"][i_t] = REWS_f
        R_FBFF["REWS_b"][i_t] = REWS_b

    return R_FBFF


#Line of sight 1*10
def WindFieldReconstruction(v_los_MD, NumberOfBeams, AngleToCenterline):
    NumberOfDistances = 10
    u_est = v_los_MD / np.cos(np.deg2rad(AngleToCenterline))

    # Persistent variable (use list for simplicity)
    if not hasattr(WindFieldReconstruction, 'u_est_Buffer'):
        WindFieldReconstruction.u_est_Buffer =  np.full((NumberOfBeams,NumberOfDistances), np.nan)

    # Estimate u component assuming perfect alignment
    u_est = v_los_MD / np.cos(np.deg2rad(AngleToCenterline))

    # Update buffer for estimated u component
    WindFieldReconstruction.u_est_Buffer = np.vstack((u_est,WindFieldReconstruction.u_est_Buffer[0:NumberOfBeams-1] ))

    # Calculate REWS from mean over all estimated u components (handling NaNs)
    REWS_MD = np.nanmean(WindFieldReconstruction.u_est_Buffer, axis=0)

    return REWS_MD



    return reshaped_arr
def CombineAndShift(REWS_MD, MeasurementDistances, IdxFirstDistance, MeanWindSpeed, DT):
    # internal variables
    nBuffer = 800  # Worst case: T_Taylor/dt with T_Taylor=180/18 s, dt=0.0125 s
    REWS_MD_Buffer = np.full((nBuffer, len(REWS_MD)), np.nan)
    # get numberOfDistances from signal
    NumberOfDistances = len(REWS_MD)

    if REWS_MD_Buffer is None:
        REWS_MD_Buffer = [[REWS_MD[0]] * NumberOfDistances for _ in range(nBuffer)]

    # update FirstInLastOut buffer
    REWS_MD_Buffer[1:nBuffer, :] = np.reshape(REWS_MD, (-1, 1)).T

    REWS_MD_Buffer= REWS_MD_Buffer[:-1]
    # get shifted values
    REWS_MD_shifted = np.empty(NumberOfDistances)
    for iDistance in range(NumberOfDistances):
        T_Taylor = (MeasurementDistances[iDistance] - MeasurementDistances[IdxFirstDistance-1]) / MeanWindSpeed
        Idx = max(0, min(nBuffer - 1, int(np.ceil(T_Taylor / DT))))
        REWS_MD_shifted[iDistance] = REWS_MD_Buffer[Idx, iDistance]

    # combine distances: overall REWS is mean over REWS from considered distances
    REWS = np.mean(REWS_MD_shifted[IdxFirstDistance-1:])

    return REWS, REWS_MD_shifted
def LPFilter(InputSignal, DT, CornerFreq):

    # Initialization of persistent variables
    if not hasattr(LPFilter, 'OutputSignalLast'):
        LPFilter.OutputSignalLast = InputSignal
        LPFilter.InputSignalLast = InputSignal

    # Define coefficients
    a1 = 2 + CornerFreq * DT
    a0 = CornerFreq * DT - 2
    b1 = CornerFreq * DT
    b0 = CornerFreq * DT

    # Filter
    OutputSignal = 1.0 / a1 * (-a0 * LPFilter.OutputSignalLast + b1 * InputSignal + b0 * LPFilter.InputSignalLast)

    # Save signals for next time step
    LPFilter.InputSignalLast = InputSignal
    LPFilter.OutputSignalLast = OutputSignal

    return OutputSignal

def Buffer(REWS, DT, T_buffer):
    if not hasattr(Buffer, 'REWS_f_Buffer'):
        nBuffer = 2000

        Buffer.REWS_f_Buffer = np.full(nBuffer, np.nan)

    # Initialize REWS_f_Buffer
    nBuffer = 2000  # Size of REWS_f_buffer, 25 seconds at 80 Hz
    Buffer.REWS_f_Buffer = np.roll(Buffer.REWS_f_Buffer, 1)
    Buffer.REWS_f_Buffer[0] = REWS

    # Index for entry at T_buffer, minimum 1, maximum nBuffer
    Idx = min(max(int(np.floor(T_buffer / DT-1)), 0), nBuffer-1)

    # Get buffered and filtered REWS from buffer
    REWS_b = Buffer.REWS_f_Buffer[Idx]

    return REWS_b
