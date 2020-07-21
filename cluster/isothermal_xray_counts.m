function result=isothermal_xray_counts(x,varargin)
%% calculates the prefactor of the xray ouptut in photon counts 
% for an isothermal sphere model
% x- is the projected distance from the center in units of Rvir
% output is in units of cm^{-2}

%% defaults
units;

zred=0;
mv=1e15;
fc=0.15;
delv=deltavir(zred);
rhoType='mean';

i=1;
while(i< length(varargin))
    switch lower(varargin{i})
        case {'mvir','mv'}
            i=i+1;
            mv=varargin{i};
        case {'deltavir','delv'}
            i=i+1;
            delv=varargin{i};
        case {'rhomean','mean'}
            rhoType='mean';
        case {'rhocrit','crit'}
            rhoType='crit';
        case {'zred','z'}
            i=i+1;
            zred=varargin{i};
        case {'fc','fg'}
            i=i+1;
            fc=varargin{i};
        otherwise
            error('isothermal_xray_counts - Illegal argument: %s',varargin{i})
    end
    i=i+1;
end


if zred>0
  delv=deltavir(zred);
end

switch rhoType
    case 'mean'
        rhoPar=rho_mean(zred).*(Ms/Mpc^3);
    case 'crit'
        rhoPar=rho_crit(zred).*(Ms/Mpc^3);
end

mv=mv.*Ms;
    

%% calculate prefactor 

result =(mm)^-2*(pi/2/(4*pi)^(1/3)/3^(5/3)).*(delv.*rhoPar).^(5/3).*fc.^2 .* mv.^(1/3).*x.^(-3);


