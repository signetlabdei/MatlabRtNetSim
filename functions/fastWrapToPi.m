function a = fastWrapToPi(a)
%FASTWRAPTOPI Equivalent to MATLAB's wrapToPi, but much faster.

a = mod(a+pi, 2*pi) - pi;

end