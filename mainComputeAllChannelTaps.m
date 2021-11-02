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
scenario = fullfile(rtSrcFolder, 'custom_scenarios/Lroom_split_joined');

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
params.rtFolder = rtSrcFolder;
tapMatrix = runAllNodePairs(scenario, params);
save(fullfile(scenario, 'allChannelTaps.mat'), 'tapMatrix', 'params', 'scenario')

%% Plot (sanity check)
% extract the first timestep for a pair of nodes
ch = tapMatrix{1, 2}(1);

figure
stem(ch.delay * 1e9, mag2db(abs(ch.iq)))
xlabel('$\tau$ [ns]')
ylabel('Tap magnitude [dB]')


%% UTILS
function tapMatrix = runAllNodePairs(scenario, params)

params.paraCfg = parameterCfg(scenario);

params.processRatios = false; % ratios are SNR, SINR, SIR, INR
params.saveHref = false; % warning: the output can be extremely large!
params.saveScalarTaps = true; % scalar taps contain delay and complex (IQ) gain

params.Ptx = 0; % Transmitted power [dBm]
params.fc = params.paraCfg.carrierFrequency; % Carrier frequency [Hz] extracted from scenario
params.BW = 1e6; % Bandiwdth [Hz] - dummy value
params.F = 5; % Noise figure [dB] - dummy value

params.dataDirection = "DL"; % ["DL", "UL"]

% Interference is not needed
params.bsInterfIdxs = []; % Node idxs of the interfering base station
params.utInterfIdxs = []; % Node idxs of the interfering user

params.numNodes = params.paraCfg.numberOfNodes;
tapMatrix = cell(params.numNodes, params.numNodes);
for i = 1:params.numNodes
    for j = i+1:params.numNodes
        % TODO: what to do with the main diagonal
        taps = runSingleNodePair(scenario, params, i, j);
        % The channel is assumed to be symmetric for the pair of nodes
        tapMatrix{i, j} = taps;
        tapMatrix{j, i} = taps;
    end
end

end

function taps = runSingleNodePair(scenario, params, i, j)
% assuming: params.dataDirection = "DL"
params.bsIdxs = i; % Index of base station nodes, corresponding to NodePosition(idx).dat
params.utIdxs = j; % Index of user nodes, corresponding to NodePosition(idx).dat

% Reference nodes are the ones whose metrics will be collected
params.txRefIdx = params.bsIdxs; % Node idx of the reference transmitter
params.rxRefIdx = params.utIdxs; % Node idx of the reference receiver

out = launchRtNetSimulation(scenario, params);
taps = out.scalarTaps;

end