function res = get_phaseSpace(gasStruct,varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


mask=true(size(gasStruct.Masses));
plotFlag;

%% parse varargin
i=1;

while i<=length(varargin)
    switch lower(varargin)
        case 'mask'
            i=i+1;
            mask=varargin{i};
        case 'noplot'
            plotFlag=false;
    


global massUnit
global densityUnit


mass=gasStruct.Masses.*massUnit;
dens=gasStruct.Density.*densityUnit;
temp=illustris.utils.calcTemperature(gasStruct.InternalEnergy,gasStruct.ElectronAbundance);




end

