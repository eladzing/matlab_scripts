function bul_pot=bulge_potential(r,varargin)
%% function to calculate the accelerattion of the bulge component by a
% Hernquist or Jaffe model.
% The reuslt is unitless and does not include factors of M or G
% r is in units of the scale radius of the profile 
% rs is the scale radius of the model in units of the scale radius of the profile 
model='hern';

i=1;
rs=0.2;
fb=0.25;
while i<=length(varargin)
    switch(varargin{i})
        case 'rs'
            i=i+1;
            rs=varargin{i};
            case 'fb'
            i=i+1;
            fb=varargin{i};
        case {'hern','hernquist'}
            model='hern';
%        case 'jaffe'
%            model='jaffe';
        case 'point'
            model='point';
            
        otherwise
            error('bulge_acccel: Illegal option')
    end
    i=i+1;
end

switch model
    case 'hern'
        bul_pot=1./(rs+r);
 %   case 'jaffe'
 %       bul_pot=log(r./(rs+r));
    case 'point'
        bul_pot=1./r;
    otherwise
        error('bulge_acccel: Illegal model')
end
bul_pot=bul_pot.*2.*fb;


