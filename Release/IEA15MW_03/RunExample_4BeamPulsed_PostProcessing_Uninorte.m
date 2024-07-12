% IEA15MW_03: IEA 15 MW monopile + realistic wind preview from a 
% 4-beam pulsed lidar system measuring at 160 m.
% This script needs to be run after RunExample_4BeamPulsed.m. 
% Purpose:
% A postprocessing version without the need to compile DLLs for lidar data
% processing to be used in the LAC Summer Games 2024. 
% To implement your own solution, replace line 57
% R_FBFF = CalculateREWSfromLidarData_LDP_v1(FBFF,DT,TMax,LDP);
% with your own function with the same inputs and outputs.
% Result:
% Cost for Summer Games 2024 ("18 m/s hurdles"):  0.444616 m/s 

%% Setup
clearvars;close all;clc;
addpath(genpath('..\WetiMatlabFunctions'))

% Seeds (can be adjusted, but will provide different results)
nSeed               = 6;                        % [-]	number of stochastic turbulence field samples
Seed_vec            = [1:nSeed]+18*100;         % [-]  	vector of seeds

% Parameters postprocessing (can be adjusted, but will provide different results)
t_start             = 60;                       % [s] 	ignore data before for STD and spectra
TMax                = 660;                      % [s]   total run time, same as in *.fst
DT                  = 0.0125;                   % [s]   time step, same as in *.fst
R                   = 120;                      % [m]  	rotor radius to calculate REWS

% Parameter for Cost (Summer Games 2024)
tau                 = 2;                        % [s]   time to overcome pitch actuator, from Example 1: tau = T_Taylor - T_buffer, since there T_filter = T_scan = 0

% configuration from LDP_v1_4BeamPulsed.IN and FFP_v1_4BeamPulsed.IN
LDP.NumberOfBeams       = 4;        % Number of beams measuring at different directions [-]               
LDP.AngleToCenterline   = 19.176;   % Angle around centerline [deg]
LDP.IndexGate           = 6;        % IndexGate
LDP.FlagLPF             = 1;        % Enable low-pass filter (flag) [0/1]
LDP.f_cutoff            = 0.1232;   % Corner frequency (-3dB) of the low-pass filter [rad/s]
LDP.T_buffer            = 1.3889;   % Buffer time for filtered REWS signal [s]

% Files (should not be be changed)
SimulationFolderLAC = 'SimulationResults_4BeamPulsed';

%% Postprocessing: evaluate data

T_buffer_v              = 1.0:0.1:2.5;
f_cutoff_v              = 0.15:0.01:0.18;
n_T_buffer              = length(T_buffer_v);
n_f_cutoff              = length(f_cutoff_v);
Cost_M                  = NaN(n_f_cutoff,n_T_buffer);

for i_f_cutoff = 1:n_f_cutoff
    for i_T_buffer = 1:n_T_buffer

        % Update T_buffer and f_cutoff
        LDP.T_buffer            = T_buffer_v(i_T_buffer);
        LDP.f_cutoff            = f_cutoff_v(i_f_cutoff);

        % Allocation
        MAE                     = NaN(1,nSeed); % mean absolute error [m/s]

        % Loop over all seeds
        for iSeed = 1:nSeed

            % Load data
            Seed                = Seed_vec(iSeed);
    	    WindFileName        = ['URef_18_Seed_',num2str(Seed,'%02d')];
            FASTresultFile      = fullfile(SimulationFolderLAC,[WindFileName,'_FlagLAC_1.outb']);
            FBFF                = ReadFASTbinaryIntoStruct(FASTresultFile);

            % Calculate REWS
            clear CalculateREWSfromLidarData_LDP_v1 % clearing all persistent variables from previous call
            R_FBFF              = CalculateREWSfromLidarData_LDP_v1(FBFF,DT,TMax,LDP);

            % Get REWS from the wind field and interpolate it on the same time vector
            TurbSimResultFile                 	= ['TurbulentWind\URef_18_Seed_',num2str(Seed,'%02d'),'.wnd'];
            [REWS_WindField,Time_WindField]  	= CalculateREWSfromWindField(TurbSimResultFile,R,2);
            REWS_WindField_Fs                   = interp1(Time_WindField,REWS_WindField,R_FBFF.Time);

            % Calculate mean absolute error
            REWS_WindField_Fs_shifted   = interp1(Time_WindField-tau,REWS_WindField,R_FBFF.Time); % shift the REWS from wind field by tau (intented prediction time) into the future (lower times)
            Error                       = REWS_WindField_Fs_shifted-R_FBFF.REWS_b; % error is  REWS from wind field shifted minus REWS from lidar filtered and buffered.
            MAE(iSeed)                  = mean(abs(Error(R_FBFF.Time>=t_start)));  % only consider error after t_start

            % Plot REWS for absolute error
            %     figure('Name',['REWS seed ',num2str(Seed)])
            %     subplot(311)
            %     hold on; grid on; box on
            %     plot(R_FBFF.Time,   REWS_WindField_Fs);
            %     plot(R_FBFF.Time,   R_FBFF.REWS);
            %     ylabel('REWS [m/s]');
            %     legend('wind field','lidar estimate')
            %     subplot(312)
            %     hold on; grid on; box on
            %     plot(R_FBFF.Time,   REWS_WindField_Fs_shifted);
            %     plot(R_FBFF.Time,   R_FBFF.REWS_b);
            %     ylabel('REWS [m/s]');
            %     legend('wind field shifted','lidar estimate filtered and buffered')
            %     subplot(313)
            %     hold on; grid on; box on
            %     plot(R_FBFF.Time,   Error);
            %     ylabel('error [m/s]');
            %     xlabel('time [s]')

        end


        Cost_M(i_f_cutoff,i_T_buffer)      = mean(MAE);
    end

end

%% Calculation of Cost for Summer Games 2024
figure
hold on;grid on;box on
plot(T_buffer_v,Cost_M,'.-')
xlabel('Tbuffer [s]')
ylabel('Cost [m/s]')
legend(f_cutoff_v'+" Hz")
Cost = min(Cost_M(:));

fprintf('Cost for Summer Games 2024 ("18 m/s hurdles"):  %f \n',Cost);