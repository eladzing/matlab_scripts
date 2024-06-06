function res = calcTemperature(internalEnergy,electronAbundance)
%CALCTEMPERATURE function to calculate temperature from internal energy and
% electron abundance in the Illustris simulation. Result is in K.
% Internal Energy is in (km/sec)^2 ( per unit mass)

global illUnits
global cosmoStruct

% define default values 
gamma = 5/3; % adiabatic index 
%boltzmann constant erg/K
Xh=cosmoStruct.Hfraction; % hydrogen mass fraction 

if ~exist('electronAbundance','var')
    ep=(Xh^-1-1)/4;
    electronAbundance=1+2*ep;
end

mu= 4./(1+3*Xh+4.*Xh.*electronAbundance).*illUnits.physUnits.mp;

res=(gamma-1).*(internalEnergy.*(illUnits.physUnits.km).^2)./illUnits.physUnits.kb.*mu;

end

