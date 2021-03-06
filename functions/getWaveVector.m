function rho = getWaveVector(theta, phi, lambda)
rho = (2 * pi/lambda) * [sin(theta(:)) .* cos(phi(:)),...
    sin(theta(:)) .* sin(phi(:)),...
    cos(theta(:))].'; % 7.1-6
end