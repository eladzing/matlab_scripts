function m_bul=bulge_density(r,varargin)
%% function to calculate the density of the bulge component by a
% Hernquist or Jaffe model.
% The reuslt is unitless and does not include factors of M or G
% r is in units of the scale radius of the disk  
% rs is the scale radius of the model in units of the scale radius of the model 

model='hern';

i=1;
rs=0.2 ;
fb=0.25;
while i<=length(varargin)
    switch(varargin{i})
        case {'rs'}
            i=i+1;
            rs=varargin{i}; 
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
        m_bul=(rs/2*pi)./(r.*(rs+r).^3);
%     case 'jaffe'
%         g_bul=1./(r.*(rs+r));
    case 'point'
        error('bulge_density: noe density for poin-mass model')
    otherwise
        error('bulge_acccel: Illegal model')
end
m_bul=m_bul.*fb;



