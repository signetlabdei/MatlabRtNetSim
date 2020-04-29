function Gmax = getAntennaElementMaxGain(antConfig)
%GETANTENNAELEMENTMAXGAIN Returns the max antenna element gain due to
%element spacing over the array. Considers single antennas, ULAs and UPAs.

% init
if ~isfield(antConfig,"M")
    antConfig.M = 1;
end
if ~isfield(antConfig,"N")
    antConfig.N = 1;
end

%% Gain Gmax
% Single element
if antConfig.M == 1 && antConfig.N  == 1
    Gmax = 0;
    return
end
% UPA
if antConfig.M > 1 && antConfig.N > 1
    Gmax = 10*log10(4*pi* antConfig.dv * antConfig.dh);
    return;
end
% ULA
if antConfig.M > 1
    spacing = antConfig.dv;
else % antConfig.N >1
    spacing = antConfig.dh;
end 

% As suggested by an RF engineer, for ULAs using truncated waveguide
% antennas, a good approximation to the element gain is obtained by using
% 0.58*lambda in the basic formula. In general, values in [0.5,0.6] make
% sense, depending on the radiation element considered.
Gmax = 10*log10(4*pi * spacing * 0.58);

end