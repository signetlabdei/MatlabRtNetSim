function out = launchRtNetSimulation(scenarioPath, params)


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