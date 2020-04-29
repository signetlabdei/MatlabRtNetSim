function antConfig = completeAntConfigStruct(antConfig)
%COMPLETEANTCONFIGSTRUCT Complete antConfig structure with default values

% As suggested by an RF engineer, for ULAs using truncated waveguide
% antennas, a good approximation to the element gain is obtained by using
% 0.58*lambda in the basic formula. In general, values in [0.5,0.6] make
% sense, depending on the radiation element considered.
if ~isfield(antConfig,"M")
    antConfig.M = 1;
    antConfig.dv = 0.58;
end
if ~isfield(antConfig,"N")
    antConfig.N = 1;
    antConfig.dh = 0.58;
end
if ~isfield(antConfig,"P")
    antConfig.P = 1;
end

end