function g_bul=bulge_accel(r,varargin)
%% function to calculate the accelerattion of the bulge component by a
% Hernquist or Jaffe model.
% The reuslt is unitless and does not include factors of M or G
% r is in units of the scale radius of the model 
% xi is the ratio of the disk scale radius to the bulge scale radius 

model='hern';

i=1;
xi=3;
fb=0.25;
while i<=length(varargin)
    switch(varargin{i})
        case {'xi'}
            i=i+1;
            xi=varargin{i}; 
        case {'fb'}
            i=i+1;
           fb=varargin{i};
        case {'hern','hernquist'}
            model='hern';
 %       case 'jaffe'
 %           model='jaffe';
        case 'point'
            model='point';
            
        otherwise
            error('bulge_acccel: Illegal option')
    end
    i=i+1;
end

switch model
    case 'hern'
        g_bul=xi.^2./(1+xi.*r).^2;
%     case 'jaffe'
%         g_bul=1./(r.*(rs+r));
    case 'point'
        g_bul=1./r.^2;
    otherwise
        error('bulge_acccel: Illegal model')
end
g_bul=g_bul.*fb.*2;



