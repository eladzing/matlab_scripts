function pval=rps_factor_expdisk(varargin)

%% *** obsolete *** 

display('****   OBSOLETE  ****')

%% calculate the prefactor for the rps toy model for exponential disk
% and isothermal ICM 
% 
%default values: 
alfa=0.5; % fudge factor 
fc=0.1;   % cluster gas fraction 
fd=0.1;  % disk gas fraction 
%Mc=1e15;  % cluster mass 
delv = 337; % delta_vir
zred=0.0;
%type='reg';
%Ms=1e10; % disk mass  - stellar
%rd=3;     % scale radius (kpc)
%re_fac=1.678347;
%sig=sigma_effective(Md,rd);
%beta=1;

sigma_s=0;
Mc=0;

i=1;
while i<=length(varargin)
    switch varargin{i}
        case{'sigma_eff'}
            i=i+1;
            sig=varargin{i};
            sigma_s=sig./sigma_factor('half');
        case{'sigma5090'}
            i=i+1;
            sig=varargin{i};
            sigma_s=sig./sigma_factor('5090');
        case{'sigma','sigma_s'}
            i=i+1;
            sigma_s=varargin{i};
        case{'alpha'}
             i=i+1;
            alfa=varargin{i};
        case{'fc'}
             i=i+1;
            fc=varargin{i};
        case{'fd'}
             i=i+1;
            fd=varargin{i};
        case{'Mc','mc','mv'}
             i=i+1;
            Mc=varargin{i};
        %case{'Md','md','Msat'}
        %    Md=varargin{i+1};
        %case{'rd','Rd'}
        %    rd=varargin{i+1};
        %case{'beta'}
        %    beta=varargin{i+1};
        case{'deltavir','delv'}
             i=i+1;
            delv=varargin{i};
        case{'zred'}
             i=i+1;
            zred=varargin{i};
        otherwise
            error('rps_factor - illegal argument: %s',varargin{i});
    end
end

%if(Mc==0)
%    error('rps_factor

sigma_cl=sigma_cluster(Mc,'delv',delv,'zred',zred);

pval=alfa.*(fc/fd).*(sigma_cl./sigma_s)

switch type
    case 'sigma_s'
        
        pval=alfa.*(fc/fd).*
        %pval1=alfa*fc/fd*(Mc/Md)^2*(rd/rc*1e-3)^4;
    case 'sigma_eff'
        pval=beta^2.*(alfa*fc/fd*(Mc/(rc*1e3)^2)^2).*(sig.*2.*pi.*re_fac^2).^-2;
       
    otherwise
        error('illegal type')
end

