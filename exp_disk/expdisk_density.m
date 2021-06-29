function res=expdisk_density(r,type,varargin)
%% function to calculate the total accelerattion in an exponential disk
% due the stellar and gaseous component. The reuslt is unitless and does
% not include factors of Sigma_0 or G
%% r is in units of the scale radius rd
% beta is the ratio between rd and rg (star & gas scale radii

beta=1;
fg=0.1;

i=1;
while i<=length(varargin)
    switch(varargin{i})
      
        case 'fg'
            i=i+1;
            fg=varargin{i};
        case 'beta'
            i=i+1;
            beta=varargin{i};
     
        otherwise
            error('disk_acccel: Illegal option')
    end
    i=i+1;
end



switch type 
    case {'stellar','stars'}
        res=exp(-r);
    case 'gas'
        res=fg.*beta^2.*exp(-beta.*r);
    case{'total','both'}
        res=exp(-r)+fg.*beta^2.*exp(-beta.*r);
    otherwise 
        error('expdisk_density: Illegal type')
end


