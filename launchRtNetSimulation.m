function out = launchRtNetSimulation(scenarioPath, params)
% Simulate the link-layer based on a given qd-realization scenario.
% The details of the simulation are passed by the params struct.
%
% channel taps for the given scenario based on the parameters.
% Taps are reduced to complex scalars by applying the TX/RX beamforming
% vectors to the channel matrix, for each node pair and for each time
% stamp.
% IN:
%   - scenarioPath: the path to the qd-realization scenario to process
%   - params: the parameters for the link-layer simulation. They comprise:
%     * processRatios (bool): ratios are SNR, SINR, SIR, INR
%     * saveHref (bool): saves all channel matrices of the reference
%     transmission. warning: the output can be extremely large!
%     * saveScalarTaps (bool): scalar taps contain delay and complex (IQ) gain
%     * bsAnt (classes.Antenna): Antenna configuration for base stations
%     * utAnt (classes.Antenna): Antenna configuration for users
%     * bsIdxs: Indexes of base station nodes, corresponding to NodePosition(idx).dat
%     * utIdxs: Indexes of user nodes, corresponding to NodePosition(idx).dat
%     * txRefIdx: Node idx of the reference transmitter
%     * rxRefIdx: Node idx of the reference receiver
%     * dataDirection: "DL" or "UL"
%     * bsInterfIdxs: Node idxs of the interfering base station(s)
%     * utInterfIdxs: Node idxs of the interfering user(s)
%     * bsSectorDir: Antenna azimuth direction for base stations [rad] (untested)
%     * bsDowntilt: Antenna downtilt for base stations: pi/2 points at the horizon [rad]
%     * utSectorDir: Antenna azimuth direction for users [rad] (untested)
%     * utDowntilt: Antenna downtilt for users: pi/2 points at the horizon [rad]
%     * Ptx: Transmitted power [dBm]
%     * fc: Carrier frequency [Hz]
%     * BW: Bandiwdth [Hz]
%     * F: Noise figure [dB]
%     * bfMode: "SVD", "geometric", or "codebook"
%     * bsCodebookFile: path to codebook used by all base stations (if
%     bfMode == "codebook")
%     * utCodebookFile: path to codebook used by all users (if
%     bfMode == "codebook")
% OUT:
%   out (struct): contains relevant link-layer metrics to the reference link:
%     * Prx_dbm (tTot, 1): received power
%     * I_dbm (tTot, 1): interference power
%     * sv1 (tTot, 1): largest singular value of Href
%     * sv2 (tTot, 1): second largest singular value of Href
%     * channelGenTime (tTot, 1): time to generate Href
%     * N_dbm: noise power
%     * fullSimTime: simulation time to run the whole function
%     * INR_db (tTot, 1): Interference-to-Noise ratio (if
%     params.processRatios)
%     * SNR_db (tTot, 1): Signal-to-Noise ratio (if params.processRatios)
%     * SIR_db (tTot, 1): Signal-to-Interference ratio (if params.processRatios)
%     * SINR_db (tTot, 1): Signal-to-(Interference and Noise) ratio (if
%     params.processRatios)
%     * Href (nAntRx, nAntTx, tTot): channel matrices (if params.saveHref)
%     * scalarTaps (tTot, 1): list of struct containing fields 'delay' [s] and
%     'iq' [linear] of the same length (if params.saveScalarTaps)


assert(isValidScenarioPath(scenarioPath),...
    '''%s'' is not a valid scenario path', scenarioPath)

t0 = tic;

% import RT traces
qdFiles = readAllQdFiles(scenarioPath);
nodesPositions = readAllNodesPositions(scenarioPath);

tTot = length(qdFiles{1,2});

% Setup antennas
if params.bfMode == "codebook"
    bsCodebook = load(params.bsCodebookFile);
    utCodebook = load(params.utCodebookFile);
    
    params.bsAnt.codebook = bsCodebook.codebook;
    params.utAnt.codebook = utCodebook.codebook;
end

% BS setup
bss = SimulationNode.empty();
for i = 1:length(params.bsIdxs)
    nodeIdx = params.bsIdxs(i);
    bss(i) = SimulationNode(params.bsAnt, params.bsSectorDir, params.bsDowntilt,...
        nodesPositions{nodeIdx}(1,:), [0,0,0], "",...
        NaN, params.Ptx, params.fc, params.BW);
    bss(i).nodeIdx = nodeIdx;
end

% UT setup
uts = SimulationNode.empty();
for i = 1:length(params.utIdxs)
    nodeIdx = params.utIdxs(i);
    uts(i) = SimulationNode(params.utAnt, params.utSectorDir, params.utDowntilt,...
        nodesPositions{nodeIdx}(1,:), [0, 0, 0], "",...
        params.F, NaN, params.fc, params.BW);
    uts(i).nodeIdx = nodeIdx;
end

[txRef, rxRef] = getTxRx(bss, uts,...
    params.txRefIdx, params.rxRefIdx, params.dataDirection);

N_dbm = rxRef.getNoise();

% Prepare loop
out.Prx_dbm = nan(tTot,1);
out.I_dbm = nan(tTot,1);
out.sv1 = nan(tTot,1);
out.sv2 = nan(tTot,1);
out.channelGenTime= nan(tTot,1);

if params.processRatios
    out.INR_db = nan(tTot,1);
    out.SNR_db = nan(tTot,1);
    out.SIR_db = nan(tTot,1);
    out.SINR_db = nan(tTot,1);
end
if params.saveHref
    out.Href = nan(rxRef.ant.getNumAnt(), txRef.ant.getNumAnt(), tTot);
end

for t = 1:tTot
    uts = updateNodesPositions(nodesPositions, t, uts);
    bss = updateNodesPositions(nodesPositions, t, bss);
    
    [Href, Htime] = getQdChannel(qdFiles, t, txRef, rxRef, params);
    [txBf, rxRefBf] = SimulationNode.getBfVectors(params.bfMode, Href, txRef, rxRef);
    % Store UT ref BF vector towards BS for interference
    % computation
    rxRef.storeBfVector(rxRefBf);
    bfGain = getBfGain(Href, NaN, rxRef.storedBfVec, txBf);
    Prx_dbm = txRef.Ptx + bfGain;
    
    I_dbm = getQdInterference(qdFiles, t, rxRef, uts, bss, params);
    
    out.Prx_dbm(t) = Prx_dbm;
    out.I_dbm(t) = I_dbm;
    [out.sv1(t), out.sv2(t)] = getSvH(Href, 2);
    out.channelGenTime(t) = Htime;
    
    if params.processRatios
        out.INR_db(t) = I_dbm - N_dbm;
        out.SNR_db(t) = Prx_dbm - N_dbm;
        out.SIR_db(t) = Prx_dbm - I_dbm;
        out.SINR_db(t) = Prx_dbm - 10*log10(10^(I_dbm/10) + 10^(N_dbm/10));
    end
    if params.saveHref
        out.Href(:,:,t) = Href;
    end
end

out.N_dbm = N_dbm;
out.fullSimTime = toc(t0);

end