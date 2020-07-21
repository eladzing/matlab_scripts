function res=sigma_factor(arg)

%% pre-factor for the stellar surface density parameter
% to convert to effective ovbservational parameter

if exist('arg','var')
    eff=arg;
else
    eff='half';
end

switch(eff)
    case {'eff','half'}
        res=1/effective_rad(0.5)^2;
    case {'5090','50-90'}
        res=0.8/(effective_rad(0.9)^2-effective_rad(0.5)^2);
    otherwise
        error('sigma_factor: Illegal argument')
end