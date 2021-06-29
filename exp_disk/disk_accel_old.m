function res=disk_accel(r,varargin)

%% Function calculates the gravitational acceleration exerted by an 
% exponential disk.   
% r - is in units of the scale radius unless stated otherwise 
% r is in units of the scale radius rd 
% beta is the ratio between rd and rg (star & gas scale radii

disp('WARNING: the function disk_accel is obsolete')

beta=1;
i=1;
while i<=length(varargin)
    switch(varargin{i})
        case 'beta'
            i=i+1;
            beta=varargin{i};
        case {'invbeta','beta_inv','Ibeta'}
            i=i+1;
            beta=1/varargin{i};
        otherwise
            error('disk_force_reduced: Illegal option')
    end
    i=i+1;
end

res=beta.^3.*r.*(besseli(0,0.5.*r*beta).*besselk(0,0.5.*r.*beta)-besseli(1,0.5.*r.*beta).*besselk(1,0.5.*r.*beta));

