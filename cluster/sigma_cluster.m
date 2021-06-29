function res=sigma_cluster(Mc,varargin)

%% function for calculating the contrived parameter of cluster surface density
% sigma_cluster which is defined:
% sigma=Mv/(2*pi*rv^2)
% in units of Msun/kpc^2

delv = 337; % delta_vir
zred=0.0;

i=1;
while i<=length(varargin)
    switch varargin{i}
        case{'deltavir','delv'}
            i=i+1;
            delv=varargin{i};
        case{'zred'}
            i=i+1;
            zred=varargin{i};
        otherwise
            error('sigma_cluster: Illegal argument')
    end
    i=i+1;
end

[rc,~,~,~]=calculate_virials('mvir',Mc,'delv',delv,'zred',zred);
rc=rc.*1e3; %convert to kpc
res=Mc./(2.*pi.*rc.^2);

