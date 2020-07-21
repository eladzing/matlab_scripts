function  res=calc_illUnits(snap)
%CALC_ILLUNITS calculate the  conversion factors from TNG simulation units to more
%reasonable things, also changes from comoving to physical units. This
%function does the calculation of the units. will also return a vector of
%values. 

% % if no snap given, use z=0;
% if ~exist('snap','var')
%     snap=99;
% end

zFac=illustris.utils.set_redshiftFactor(snap);

global cosmoStruct
%global illUnits 
units;
res.snap=snap;
res.aexp=zFac.length;
res.zred=zFac.zred;
res.massUnit=1e10/cosmoStruct.hub;  % Change simulation units to solar mass 
res.lengthUnit=1/cosmoStruct.hub.*zFac.length; % change units to kpc
res.volumeUnit=(res.lengthUnit).^3; % change volume to kpc^3
res.areaUnit=(res.lengthUnit).^2; % change area to kpc^2

res.densityUnit=res.massUnit./res.volumeUnit; % Change simulation units to solar mass / kpc^3
res.surfaceDensityUnit=res.massUnit./res.areaUnit; % Change simulation units to solar mass / kpc^3

% some simulation use a strange time unit of 0.978 Gyr / h
% we convert this to years 
res.timeUnit= 0.978 * 1.e9 / cosmoStruct.hub;  % in years 

% convert simulation density to number density in cm^-3
res.numberDensityFactor = (res.densityUnit*Units.Ms/Units.kpc^3)./(cosmoStruct.muMass.*Units.mp);

% convert density to Hydrogen number density in cm^-3
res.HydrogenNumberDensityFactor= (res.densityUnit.*Units.Ms./Units.kpc^3).*(cosmoStruct.Hfraction/Units.mp);

% BH mdot
res.BHMdotFactor=res.massUnit./ res.timeUnit;  % This is in solar-mass per year.

% BH Energy - This is in units of 1e53 ergs.
res.BHEnergyFactor=(res.massUnit.*Units.Ms).*(res.lengthUnit.*Units.kpc).^2./(res.timeUnit.*Units.yr)^2./1e53;  

% energy dissipation in gas 
res.EnergyDissipationUnit=(res.massUnit./res.lengthUnit*Units.Ms/Units.kpc)...
    .*(Units.km)^3.*(res.aexp).^(-1)*Units.yr./1e45; % This is in units of 10^45 erg/year 


end