function A = threeGppElementPattern(tRad, fRad)

tDeg = rad2deg(tRad);
fDeg = rad2deg(fRad);

Av = -min(12 * ((tDeg - 90) / 65) .^ 2, 30);
Ah = -min(12 * (fDeg / 65) .^ 2, 30);
Adb = 8 - min(-(Av + Ah), 30);
A = 10 .^ (Adb / 20);

end