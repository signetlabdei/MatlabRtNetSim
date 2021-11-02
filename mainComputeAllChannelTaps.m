clear
close all
clc

rtSrcFolder = "../qd-realization/src";

addpath("classes",...
    "functions",...
    fullfile("functions", "rayTracerUtils"),...
    rtSrcFolder,...
    fullfile(rtSrcFolder, "utils"))

%% setup
scenario = fullfile(rtSrcFolder, 'custom_scenarios/BoxConferenceRoom_split_joined');

params.bsAnt = Antenna(8, 8, 0.5, 0.5, @(t,f) 1); % Antenna configuration for base stations
params.utAnt = Antenna(1, 1, 0.5, 0.5, @(t,f) 1); % Antenna configuration for users

params.bsSectorDir = 0; % Antenna direction for base stations [rad] (untested)
params.bsDowntilt = pi/2; % Antenna downtilt for base stations: pi/2 points at the horizon [rad]
params.utSectorDir = 0; % Antenna direction for users [rad] (untested)
params.utDowntilt = pi/2; % Antenna downtilt for users: pi/2 points at the horizon [rad]

params.bfMode = "SVD"; % ["SVD", "geometric", "codebook"]
params.bsCodebookFile = sprintf("codebooks/%dx%d_noTaper.mat", params.bsAnt.M, params.bsAnt.N);
params.utCodebookFile = sprintf("codebooks/%dx%d_noTaper.mat", params.utAnt.M, params.utAnt.N);

%% Simulation
tapMatrix = computeAllChannelTaps(scenario, params, rtSrcFolder);
save(fullfile(scenario, 'allChannelTaps.mat'), 'tapMatrix', 'params', 'scenario')

%% Plot (sanity check)
% extract the first timestep for a pair of nodes
ch = tapMatrix{1, 2}(1);

figure
stem(ch.delay * 1e9, mag2db(abs(ch.iq)))
xlabel('$\tau$ [ns]')
ylabel('Tap magnitude [dB]')
