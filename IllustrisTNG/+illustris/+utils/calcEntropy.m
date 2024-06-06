function res = calcEntropy(temperature,density)
%CALCENTROPY function to calculate the entropy from the temperature  (in K) 
% and the density in simulation units. The density is converted to number density 
% and the result is given in KeV cm^2 

units;
%global cosmoStruct
global illUnits


ndensity=density.*illUnits.numberDensityFactor;
res=illUnits.physUnits.kb.*temperature./(ndensity).^(2/3)./(illUnits.physUnits.ev.*1000);

end

