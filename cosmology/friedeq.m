function result = friedeq(zz,varargin)
% Freidmans equation for the hubble parameter
% in units of km/sec/Mpc.
% for a Lambda-CDM concordance model.



%% set defualts:
%global hub
h0 = 0.7 ;
Omm = 0.3;
Oml= 0.7;

if isempty(h0)
    h0=0.7;
end

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
                error('OMEGA_M: cosmology structure must be a structure')
            end
            Omm=cosmoStruct.Omm;
            Oml=cosmoStruct.Oml;
            h0=cosmoStruct.hub;
        otherwise
            error('friedeq - illegal argument: %s',varargin{i})
    end
    i=i+1;
end

if h0>1
    h0=h0./100;
end

zp1=zz+1;

result=h0*100.*sqrt(Omm.*zp1.^3+Oml+(1-Oml-Omm).*zp1.^2);

end