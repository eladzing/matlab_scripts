function  res= redshift2time(zred,varargin)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


%% set defualts:
%global hub
h0 = 0.7 ;
Omm = 0.3;
Oml= 0.7;

if isempty(h0)
    h0=0.7;
end

zend=500;


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
        case {'zend'}
            i=i+1;
            zend=varargin{i};
            
        otherwise
            error('friedeq - illegal argument: %s',varargin{i})
    end
    i=i+1;
end

%% set hubble time 
units; 

H=h0.*100.*Units.km./Units.Mpc; % Hubble 
th=1/H./(1e9.*Units.yr);



%% all time 
zr0=0:0.001:zend;
zpl=zr0+1;

EofZ=sqrt(Omm.*zpl.^3 + Oml + (1-Oml-Omm).*zpl.^2);
integrand=1./(zpl.*EofZ);

time=th.*cumtrapz(zr0,integrand);


%% lookback 

% zr=0:0.0001:zred;
% zpl=zr+1;
% 
% EofZ=sqrt(Omm.*zpl.^3 + Oml + (1-Oml-Omm).*zpl.^2);
% integrand=1./(zpl.*EofZ);
% 
% lookback(1)=th.*trapz(zr,integrand);

res.lookback=interp1(zr0,time,zred,'pchip');

%% age 

% zr=zred:0.0001:200;
% zpl=zr+1;
% 
% EofZ=sqrt(Omm.*zpl.^3 + Oml + (1-Oml-Omm).*zpl.^2);
% integrand=1./(zpl.*EofZ);
% 
% age(1)=th.*trapz(zr,integrand);

res.age=time(end)-res.lookback;



end

