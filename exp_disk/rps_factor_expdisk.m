function pval=rps_factor_expdisk(varargin)
%% calculate the prefactor for the rps toy model for exponential disk
% and isothermal ICM 
% 
%default values: 
alfa=0.5; % fudge factor 
fc=0.1;   % cluster gas fraction 
fd=0.1;  % disk gas fraction 
zred=0.0;
rp=1;  % position in cluster in units of Rvir
sigma_s=1;
Mc=1;

delv = deltavir(zred); % delta_vir

sigmaclFlag=false;
McFlag=false;
sigmastFlag=false;


i=1;
while i<=length(varargin)
    switch varargin{i}
        case{'rp','etap','pos'}
            i=i+1;
            rp=varargin{i};
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
            sigmastFlag=true;
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
            McFlag=true;
        case{'sigma_cl','sigma_cluster'}
             i=i+1;
            sigma_cl=varargin{i};
            sigmaclFlag=true;
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
    i=i+1;
end

if sigmaclFlag && McFlag 
    error('rps_factor_expdisk: Mc & Sigma_cl cannot be enterd together');
elseif ~sigmaclFlag && ~McFlag
    error('rps_factor_expdisk: Must enter cluster mass Mc or sigma_cl');
end


if ~sigmastFlag
    error('rps_factor_expdisk: Must enter sigma disk (star,effective,5090)');
end

if McFlag
    sigma_cl=sigma_cluster(Mc,'delv',delv,'zred',zred);
end

pval=alfa.*(fc./fd).*(sigma_cl./sigma_s).^2.*rp.^(-2);

