function [H, Htime] = getQdChannel(qdFiles, t, txNode, rxNode, params)

qdFile = qdFiles{txNode.nodeIdx, rxNode.nodeIdx}(t);

txDowntilt = txNode.downtilt;
txSectorDir = txNode.sectorDir; % assuming scalar sectorDir
rxDowntilt = rxNode.downtilt;
rxSectorDir = rxNode.sectorDir; % assuming scalar sectorDir

t0 = tic;

if qdFile.numRays == 0
    H = zeros(rxNode.ant.getNumAnt(), txNode.ant.getNumAnt());
    
else
    txAng.el = deg2rad(qdFile.aodEl');
    txAng.az = deg2rad(qdFile.aodAz');
    rxAng.el = deg2rad(qdFile.aoaEl');
    rxAng.az = deg2rad(qdFile.aoaAz');
    
    txSteerVec = txNode.ant.getSteeringVector(txAng.el, txAng.az, txDowntilt, txSectorDir);
    txElPatt = txNode.ant.getElementPattern(txAng.el, txAng.az, txDowntilt, txSectorDir);
    txVec = txElPatt .* txSteerVec;
    
    rxSteerVec = rxNode.ant.getSteeringVector(rxAng.el, rxAng.az, rxDowntilt, rxSectorDir);
    rxElPatt = rxNode.ant.getElementPattern(rxAng.el, rxAng.az, rxDowntilt, rxSectorDir);
    rxVec = rxElPatt .* rxSteerVec;
    
    signalPhase = exp(1i * (-2*pi * qdFile.delay * params.fc + qdFile.phaseOffset)).';
    pathGain_dB = qdFile.pathGain';
    pathGainMag_lin = 10.^(pathGain_dB / 20);
    H = pathGainMag_lin .* signalPhase .* conj(rxVec) * txVec';
    
end

Htime = toc(t0);

end