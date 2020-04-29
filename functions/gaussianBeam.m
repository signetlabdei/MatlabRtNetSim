function E = gaussianBeam(el, az, el0, az0, beamwidth, A)
E = A * exp(-((el-el0).^2 + (az-az0).^2) / (beamwidth^2 / log(4)));

end