function [bfGain, tapsIq] = getBfGain(H,tau,bfRx,bfTx)
% H is 3d with last dimension representing taps
% assume bf vectors are columns
tapsIq = squeeze(sum(bfRx .* H .* bfTx.', [1,2]));
bfGain = 20*log10(abs(sum(tapsIq)));

end