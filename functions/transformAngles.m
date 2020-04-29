function [theta1, phi1] = transformAngles(theta, phi, downtilt, sectorDir)
% From eqs. 7.1-6, 7.1-7, 7.1-8 of 3GPP TR 38.901
assert(0 <= downtilt && downtilt <= pi)
assert(all(size(theta) == size(phi)))

% alpha, beta, gamma
a = sectorDir;
b = downtilt - pi/2;
g = 0;

% transformations
rho = getWaveVector(theta, phi, 2*pi);
Rinv = [cos(a)*cos(b), sin(a)*cos(b), -sin(b);...
    cos(a)*sin(b)*sin(g) - sin(a)*cos(g), sin(a)*sin(b)*sin(g) + cos(a)*cos(g), cos(b)*sin(g);...
    cos(a)*sin(b)*cos(g)+sin(a)*sin(g), sin(a)*sin(b)*cos(g) - cos(a)*sin(g), cos(b)*cos(g)];

theta1 = acos([0, 0, 1] * Rinv * rho); % 7.1-7
phi1 = angle([1, 1j, 0] * Rinv * rho); % 7.1-8

% reshape to original form
theta1 = reshape(theta1, size(theta));
phi1 = reshape(phi1, size(phi));

end