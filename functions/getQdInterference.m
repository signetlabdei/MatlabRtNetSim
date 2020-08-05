function I_dbm = getQdInterference(qdFiles, t, rxRef, uts, bss, params)
assert(length(params.bsInterfIdxs) == length(params.utInterfIdxs),...
    'bsInterfIdxs and utInterfIdxs must have the same length')

I = 0;
for i = 1:length(params.bsInterfIdxs)
    [txInterf, rxInterf] = getTxRx(bss, uts,...
        params.bsInterfIdxs(i), params.utInterfIdxs(i), params.dataDirection);
    
    % Comput interferent BF
    H = getQdChannel(qdFiles, t, txInterf, rxInterf, params);
    txInterfBf = SimulationNode.getBfVectors(params.bfMode, H, txInterf, rxInterf);
    
    HInterf = getQdChannel(qdFiles, t, txInterf, rxRef, params);
    bfGain = getBfGain(HInterf, NaN, rxRef.storedBfVec, txInterfBf);
    I = I + db2pow(txInterf.Ptx + bfGain);
end

I_dbm = pow2db(I);

end