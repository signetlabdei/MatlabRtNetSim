function I_dbm = getQdInterference(qdFiles, t, rxRef, ut, bs, params)
I = 0;
switch(params.dataDirection)
    case "DL"
        txInterf = bs([bs.nodeIdx] == params.bsInterfIdxs);
        rxInterf = ut([ut.nodeIdx] == params.utInterfIdxs);
    case "UL"
        txInterf = ut([ut.nodeIdx] == params.utInterfIdxs);
        rxInterf = bs([bs.nodeIdx] == params.bsInterfIdxs);
end

for i = 1:length(txInterf)
    % Comput interferent BF
    H = getQdChannel(qdFiles, t, txInterf(i), rxInterf(i), params);
    txInterfBf = SimulationNode.getBfVectors(params.bfMode, H, txInterf(i), rxInterf(i));
    
    HInterf = getQdChannel(qdFiles, t, txInterf(i), rxRef, params);
    bfGain = getBfGain(HInterf, NaN, rxRef.storedBfVec, txInterfBf);
    I = I + db2pow(txInterf(i).Ptx + bfGain);
end

I_dbm = pow2db(I);    

end