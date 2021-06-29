function [zeta_vir] = zeta_strip(varargin)
%% function for simple toy model of disk strippng a-la Gunn & Gott 77 

% defualt values 
alpha = 0.5 ;  % stripping efficiency
fg = 0.1 ;     % gas fraction of ICM  
fd = 0.25 ;    % gas fraction of disk 
Mc = 1e14 ;    % Mass of cluster [Msun]
Md = 1e11 ;    % Mass of disk [Msun]
rd = 20 ;      % Outer radius of disk [kpc]
zred = 0.0;     % redshift
delvflag=false;

for i=1:2:nargin
    switch varargin{i}
        case{'alpha'}
            alpha=varargin{i+1};
        case{'fg'}
            fg=varargin{i+1};
        case{'fd'}
            fd=varargin{i+1};
        case{'Mc','mc','mv'}
            Mc=varargin{i+1};
        case{'Md','md'}
            Md=varargin{i+1};
        case{'rd','Rd'}
            rd=varargin{i+1};
        case{'zred'}
            zred=varargin{i+1};
        case{'deltavir'}
            delv=varargin{i+1};
            delvflag=true;
        otherwise
            error('zeta_strip - illegal argument: %s',varargin{i});
    end
end

%% calculate 
units;
Pdisk=2.*fd.*(Md.*Ms).^2/(rd*kpc).^4;

if ~delvflag
    delv=deltavir(zred);
end

Pram=alpha./4.0.*fg.*(Mc.*Ms).^(2/3)*(4*pi/3*delv*rho_mean(zred)*Ms/Mpc^3)^(4/3);

zeta_vir=Pram./Pdisk;
