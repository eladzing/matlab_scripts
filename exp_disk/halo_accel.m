function g_hal=halo_accel(r,mvms,rs,varargin)
%% function to calculate the accelerattion of the halo component by a
% NFW, Hernquist or  isothermal model.
% The reuslt is unitless and does not include factors of M or G
% r is in units of the scale radius of the disk 
% mv ms is the ratio of Mvir to the disk stellar mass
% rs is 
%  1. the scale radius of the hernquist model in units of the disk scale
%  2. the ratio Rd/Rv for the isothermal and NFW case 
% radius
% cv is the concentration parameter (NFW)

model='nfw';

i=1;
cv=10;
%fb=0.25;

while i<=length(varargin)
    switch(varargin{i})
        case {'cv'}
            i=i+1;
            cv=varargin{i}; 
        
        case {'hern','hernquist'}
            model='hern';
 %       case 'jaffe'
 %           model='jaffe';
        case 'isothermal'
            model='isothermal';
        case {'nfw','NFW'}
            model='nfw';
            
        otherwise
            error('bulge_acccel: Illegal option')
    end
    i=i+1;
end

switch model
    case 'hern'
        g_hal=1./(r+rs).^2;  
    case 'isothermal'
         g_hal=1./(r.*rs);
    case 'nfw'
        xx=rs.*r;
        g_hal=((log(1+cv.*xx)-cv.*xx./(1+cv.*xx))./(log(1+cv)-cv./(1+cv)))./r.^2;
    otherwise
        error('bulge_acccel: Illegal model')
end
g_hal=g_hal.*2.*mvms;



