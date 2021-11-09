function tapMatrix = computeAllChannelTaps(scenario, params)
% Compute channel taps for the given scenario based on the parameters.
% Taps are reduced to complex scalars by applying the TX/RX beamforming
% vectors to the channel matrix, for each node pair and for each time
% stamp.
% IN:
%   - scenario: the path to the qd-realization scenario to process
%   - params: the parameters needed by the `launchRtNetSimulation` function
% OUT:
%   tapMatrix: NxN cell array (where N is the number of nodes) containing
%     arrays of structs of length T (where T is the number of time steps). The
%     structs contains fields `delay` and `iq`, numeric arrays representing
%     complex taps.
%     Cells (i,j) and (j,i) contain the same taps (symmetric channel) and cell
%     (i,i) (i.e., the self-channel) is empty.

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