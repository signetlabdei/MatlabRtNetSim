function out = launchRtNetSimulation(scenarioPath, params)

addpath("classes",...
    "functions",...
    params.rtFolder,...
    fullfile(params.rtFolder, "utils"))

assert(isValidScenarioPath(scenarioPath),...
    '''%s'' is not a valid scenario path', scenarioPath)

t0 = tic;

% import RT traces
qdFiles = readAllQdFiles(scenarioPath);
nodesPositions = readAllNodesPositions(scenarioPath);

tTot = length(qdFiles{1,2});

% BS setup
bs = SimulationNode.empty();
for i = 1:length(params.bsIdxs)
    nodeIdx = params.bsIdxs(i);
    bs(i) = SimulationNode(params.bsAnt, params.bsSectorDir, params.bsDowntilt,...
        nodesPositions{nodeIdx}(1,:), [0,0,0], "",...
        NaN, params.Ptx, params.fc, params.BW);
    bs(i).nodeIdx = nodeIdx;
end

% UT setup
ut = SimulationNode.empty();
for i = 1:length(params.utIdxs)
    nodeIdx = params.utIdxs(i);
    ut(i) = SimulationNode(params.utAnt, params.utSectorDir, params.utDowntilt,...
        nodesPositions{nodeIdx}(1,:), [0, 0, 0], "",...
        params.F, NaN, params.fc, params.BW);
    ut(i).nodeIdx = nodeIdx;
end

switch(params.dataDirection)
    case "DL"
        txRef = bs([bs.nodeIdx] == params.txRefIdx);
        rxRef = ut([ut.nodeIdx] == params.rxRefIdx);
    case "UL"
        txRef = ut([ut.nodeIdx] == params.txRefIdx);
        rxRef = bs([bs.nodeIdx] == params.rxRefIdx);
    otherwise
        error("Data direction '%s' not recognized",params.dataDirection)
end

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
    ut = updateNodesPositions(nodesPositions, t, ut);
    bs = updateNodesPositions(nodesPositions, t, bs);
    
    [Href, Htime] = getQdChannel(qdFiles, t, txRef, rxRef, params);
    [txBf, rxRefBf] = SimulationNode.getBfVectors(params.bfMode, Href, txRef, rxRef);
    % Store UT ref BF vector towards BS for interference
    % computation
    rxRef.storeBfVector(rxRefBf);
    bfGain = getBfGain(Href, NaN, rxRef.storedBfVec, txBf);
    Prx_dbm = txRef.Ptx + bfGain;
    
    I_dbm = getQdInterference(qdFiles, t, rxRef, ut, bs, params);
    
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