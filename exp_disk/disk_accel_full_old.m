function [g_tot g_st g_gs]=disk_accel_full(r,varargin)
%% function to calculate the total accelerattion in an exponential disk
% due the stellar and gaseous component. The reuslt is unitless and does
% not include factors of Sigma_0 or G
%% r is in units of the scale radius rd
% beta is the ratio between rd and rg (star & gas scale radii
% fg = gas fraction
% fs = stellar fraction

disp('*** WARNING: the function disk_accel_full is obsolete *** ');


beta=1;
fg=0;

i=1;
while i<=length(varargin)
    switch(varargin{i})
        
        case 'fg'
            i=i+1;
            fg=varargin{i};
        case 'beta'
            i=i+1;
            beta=varargin{i};
        case {'invbeta','beta_inv','Ibeta'}
            i=i+1;
            beta=1/varargin{i};
        otherwise
            error('disk_acccel: Illegal option')
    end
    i=i+1;
end

fs=1-fg;

g_st=fs.*disk_accel(r,'beta',1);
g_gs=fg.*disk_accel(r,'beta',beta);

g_tot=g_st+g_gs;




