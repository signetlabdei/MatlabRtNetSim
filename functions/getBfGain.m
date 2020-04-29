function bfGain = getBfGain(H,tau,bfRx,bfTx)

% Narrowband
H = sum(H,3);
bfGain = 20*log10(abs(bfRx.' * H * bfTx));

end