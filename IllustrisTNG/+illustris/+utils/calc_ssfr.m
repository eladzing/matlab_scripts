function ssfr = calc_ssfr(subStruct,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global illUnits 

base=10^0.5*1e-17;
scat=0.5;

i=1;
while i<=length(varargin)
    switch(lower(varargin{i}))
        case {'base','floor'}
            i=i+1;
            base=varargin{i};
        case {'scatter','scat','fuzz'}
            i=i+1;
            scat=varargin{i};
        otherwise
            error('CALC_SSFR - illegal argument: %s',varargin{i});
    end
    i=i+1;
end

  
stellarMass= subStruct.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit; % stellar mass within 2*rhalf
sfr= subStruct.SubhaloSFRinRad;

ssfr=sfr./stellarMass; % + base.*10.^(scat.*rand(size(sfr)));

 mask=ssfr==0;
 
 ssfr(mask)=base.*10.^(scat.*rand(1,sum(mask)));

end

