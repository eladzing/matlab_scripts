function res = mass_entropy_relation(mass,varargin)
%MASS_ENTROPY_RELATION Calculate the virial mass entropy relation Summary of this function goes her
%   Based on the relation K= kb * T_vir  / n^(2/3) where n ~ \delta
%   \rho_ref = const. result is given in Kev cm^2

units;
delta=200;
rhoType='crit';
zred=0;
cosmoFlag=false;
fgas=[];

i=1;
while (i<=length(varargin))
   % fprintf('%s \n',lower(varargin{i}))
    switch lower(varargin{i})
        case{'zred'}
            i=i+1;
            zred=varargin{i};
            
        case{'delta','deltavir','delta_vir'}
            i=i+1;
            delta=varargin{i};
            
            
        case{'rho','rhotype'}
            i=i+1;
            rhoType=varargin{i};
            
        case{'cosmo','cosmology','cosmostruct'}
            i=i+1;
            cosmoStruct=varargin{i};
            cosmoFlag=true;
            
        case{'fg','fgas','gasfraction'}
            i=i+1;
            fgas=varargin{i};
        otherwise
            error('MASS_ENTROPY_RELATION - Illegal argument: %s',varargin{i});
    end
    i=i+1;
end

if ~cosmoFlag
    cosmoStruct=set_LCDM_cosmology('quiet');
end

switch lower(rhoType)
    case 'mean'
        rhoRef=rho_mean(zred,'cosmo',cosmoStruct);
    case 'crit'
        rhoRef=rho_crit(zred,'cosmo',cosmoStruct);
    otherwise
        error('MASS_ENTROPY_RELATION - Illegal density type: %s',rhoType);
end

if ischar(delta)
    if strcmpi(delta,'vir') || strcmpi(delta,'virial')
        delt=deltavir(zred,cosmoStruct.Omm,cosmoStruct.Oml);
    else
        error('MASS_ENTROPY_RELATION - unknown delta type: %s',delta);
    end
elseif isnumeric(delta)
    delt=delta;
else
    error('MASS_ENTROPY_RELATION - unknown delta type: %s',delta);
end

if isempty(fgas)
    fgas=cosmoStruct.Omb./cosmoStruct.Omm;
end


%% actual calculation

preFac=0.5*Units.G*(4*pi/3)^(1/3)*...
    (cosmoStruct.muMass.*Units.mp)^(5/3)*...
    (delt*rhoRef*Units.Ms/Units.Mpc^3)^(-1/3)/(1000*Units.ev);
res=preFac*fgas.^(-2/3).*Units.Ms^(2/3).*(mass).^(2/3);




end

