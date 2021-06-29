function result = angular_distance(zend,varargin )
%ANGULAR_DISTANCE calculates the angular distance to an object for a given
%redshift in a cosmological model defined by the parameters Omega_m ,
%Omega_l and h at z=0. The distance is found by integration over the
%friedmann equation. Result is in Mpc;


% if length(zend)>1
%     error('enter only one redshift')
% end

%% set defualts: 
global hub 
h0 = hub ;
Omm = 0.3; 
Oml= 0.7; 
res=1e-4;
units;

i=1;
while i<=length(varargin)
    switch lower(varargin{i})
        case {'h0','hubble','hub'}
            i=i+1;
            h0=varargin{i};
        case {'omegam','omm','omega_m'}
            i=i+1;
            Omm=varargin{i};
        case {'omegal','oml','omega_lambda','lambda'}
            i=i+1;
            Oml=varargin{i};
        case{'cosmostruct','cosmo','cosmology'}
            i=i+1;
            cosmoStruct=varargin{i};
            if ~isstruct(cosmoStruct)
                error('comoving distance: cosmology structure must be a structure')
            end
            Omm=cosmoStruct.Omm;
            Oml=cosmoStruct.Oml;
            h0=cosmoStruct.hub;

        
        otherwise
            error('friedeq - illegal argument: %s',varargin{i})
    end
    i=i+1;
end

if isempty(h0)
    h0=0.7;
end

if h0>1
    h0=h0./100;
end

for j=1:length(zend)

zz=0:zend(j)*res:zend(j);

eofz=friedeq(zz,'hub',h0,'Omm',Omm,'Oml',Oml)./(100.*h0);

p1=trapz(zz,1./eofz);

result(j)=(Units.cspeed./Units.km)/(100.*h0)/(1+zend(j)).*p1;




end

