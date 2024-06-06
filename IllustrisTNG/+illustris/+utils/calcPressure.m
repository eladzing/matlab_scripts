function res = calcPressure(nDensity,temperature)
%CALCPRESSURE function to calculate PRESSURE based on temperature (in K) and
% number denstiy (in cm^-3). The result is in cgs units - erg/cm^3 .

global illUnits

res=nDensity.*temperature.*illUnits.physUnits.kb;


end

