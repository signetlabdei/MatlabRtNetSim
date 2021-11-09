function tapMatrix = computeAllChannelTaps(scenario, params)

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

%% UTILS

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