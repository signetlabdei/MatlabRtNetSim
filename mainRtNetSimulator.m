clear
close all
clc

rtFolder = "../qd-realization/src";

addpath("classes",...
    "functions",...
    fullfile("functions", "rayTracerUtils"),...
    rtFolder,...
    fullfile(rtFolder, "utils"))

%% setup
qdCampaignFolder = fullfile(rtFolder, 'custom_scenarios/Journal1Lroom');
scenario = 'refl2_qd0_relTh-40_floorMetal';

params.processRatios = true; % ratios are SNR, SINR, SIR, INR
params.saveHref = false; % warning: the output can be extremely large!

params.bsAnt = Antenna(8, 8, 0.5, 0.5, @(t,f) 1); % Antenna configuration for base stations
params.utAnt = Antenna(4, 4, 0.5, 0.5, @(t,f) 1); % Antenna configuration for users

params.bsIdxs = [2, 4]; % Index of base station nodes, corresponding to NodePosition(idx).dat
params.utIdxs = [1, 3]; % Index of user nodes, corresponding to NodePosition(idx).dat

params.bsSectorDir = 0; % Antenna direction for base stations [rad] (untested)
params.bsDowntilt = pi/2; % Antenna downtilt for base stations: pi/2 points at the horizon [rad]
params.utSectorDir = 0; % Antenna direction for users [rad] (untested)
params.utDowntilt = pi/2; % Antenna downtilt for users: pi/2 points at the horizon [rad]

params.Ptx = 30; % Transmitted power [dBm]
params.fc = 60e9; % Carrier frequency [Hz]
params.BW = 29 * 13.889e6; % Bandiwdth [Hz]
params.F = 5; % Noise figure [dB]

params.dataDirection = "DL"; % ["DL", "UL"]
params.bfMode = "SVD"; % ["SVD", "geometric"]

% Reference nodes are the ones whose metrics will be collected
params.txRefIdx = 2; % Node idx of the reference transmitter
params.rxRefIdx = 1; % Node idx of the reference receiver

% Nodes can interfere with the reference ones
% They are defined by vector, coupling corresponding rx/tx by index
% E.g.,
% params.bsInterfIdxs = [1, 1, 2]
% params.utInterfIdxs = [3, 4, 5]
% Will create interference with BS 1 communicating with users 3 and 4, and
% BS 2 communicating with user 5.
% The communication direction is set in params.dataDirection
params.bsInterfIdxs = [4]; % Node idxs of the interfering base station
params.utInterfIdxs = [3]; % Node idxs of the interfering user

%% simulation
params.rtFolder = rtFolder;
out = launchRtNetSimulation(fullfile(qdCampaignFolder, scenario), params);

%% Plot (sanity check)
figure
plot(out.SINR_db, 'DisplayName', 'SINR'); hold on
plot(out.SNR_db, 'DisplayName', 'SNR')
xlabel('t [step]')
ylabel('[dB]')
legend('show')