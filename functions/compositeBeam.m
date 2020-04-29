function E = compositeBeam(beam, el, az, el0, az0)
E = 0;
for i = 1:length(el0)
    E = E + beam(el, az, el0(i), az0(i));
    E = E + beam(el, az, el0(i), az0(i) + 360); % wrap angles
    E = E + beam(el, az, el0(i), az0(i) - 360); % wrap angles
end

end