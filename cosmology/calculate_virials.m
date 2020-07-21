function [rvir, mvir, tvir, vvir]=calculate_virials(param,arr,varargin)
% given one of the virial parameters and a choice for delta_vir (defualt=200)
% calculate all other virial quantities

%% units in cgs
% kb= 1.38066e-16; %erg/K
% mu= 0.5926; %mass per particle for primordial composition
% mp= 1.672649e-24; % gram
% pc= 3.0856e18; % parsec
% km= 1e5; %cm
% Ms= 1.9891e33; %gr
% G= 6.67e-8; % Gravitational Constant (cm^3 g^-1 s-2)
% Mpc=1e6*pc;

%delta_vir=337;
zred=0;
hub=0.7;
Omm=0.3;
Oml=1-Omm;
cosmoStructFlag=false;
rhoType='mean';

units;
mumass=Units.muMass;


if nargin<2
    error('calculate virials: not enough arguments, must be at least (mvir,val) or (rvir,val)');
end

delvFlag=false;

i=1;
while i<=length(varargin)
    switch lower(varargin{i})
        case {'h','hub','hubble'}
            i=i+1;
            hub=varargin{i};
        case {'omm', 'omega_matter','om_m'}
            i=i+1;
            Omm=varargin{i};
        case {'oml', 'omega_lambda','om_l'}
            i=i+1;
            Oml=varargin{i};
        case {'zred','z','redshift'}
            i=i+1;
            zred=varargin{i};
        case {'delv','delta','delta_vir','deltavir'}
            i=i+1;
            delta_vir=varargin{i};
            delvFlag=true;
        case {'mean','rho_mean','rhomean'}
            rhoType='mean';
        case {'crit','rho_crit','rhocrit'}
            rhoType='crit';
        case {'cosmostruct','cosmo'}
            i=i+1;
            cosmoSt=varargin{i};
            
            if ~isstruct(cosmoSt)
                error('CALCULATE_VIRIALS: cosmology structure must be a structure')
            end
            mumass=cosmoSt.muMass;
            cosmoStructFlag=true;
        otherwise
            error('calculate virials: illegal argument: %s',varargin{i})
    end
    i=i+1;
end

if ~cosmoStructFlag
    cosmoSt.Omm=Omm;
    cosmoSt.Oml=Oml;
    cosmoSt.hub=hub;
end
    
   

% find mean density of the universe
switch rhoType
    case 'mean'
        rhoRef=rho_mean(zred,'cosmoStruct',cosmoSt);
    case 'crit'
        rhoRef=rho_crit(zred,'cosmoStruct',cosmoSt);
end

if ~delvFlag
    delta_vir=deltavir(zred,cosmoSt.Omm,cosmoSt.Oml);
end


switch param
    case{'mvir','mv'}
        mvir=double(arr);
        rvir=(3.*mvir./(4*pi*delta_vir*rhoRef)).^(1/3); % in Mpc
        vvir=sqrt(Units.G.*(mvir.*Units.Ms)./(rvir.*Units.Mpc))./Units.km; %in km/sec
        tvir=(vvir.*Units.km).^2*mumass*Units.mp/(2*Units.kb); % inKelvin
    case{'rvir','rv'}
        rvir=double(arr);
        mvir=(4*pi/3).*rvir.^3.*(delta_vir*rhoRef);
        vvir=sqrt(Units.G*(mvir.*Units.Ms)./(rvir.*Units.Mpc))./Units.km; %in km/sec
        tvir=(vvir.*Units.km).^2*mumass*Units.mp/(2*Units.kb); % inKelvin
        %     case{'vvir','vv'}
        %         vvir=value;
        %         vvir=sqrt(G*(mvir*Ms)/(rvir.*1e6*pc))/km; %in km/sec
        %         tvir=(vvir*km)^2*mu*mp/(2*kb); % inKelvin
        %     case{'tvir','tv'}
        %         tvir=values;
        %         vvir=sqrt(G*(mvir*Ms)/(rvir.*1e6*pc))/km; %in km/sec
        %         tvir=(vvir*km)^2*mu*mp/(2*kb); % inKelvin
    otherwise
        error('calculate virials: illegal virial parameter')
end
