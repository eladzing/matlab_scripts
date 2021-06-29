function result = deltavir(zred,varargin)
% find delta vir given redshift and the comsological parameters at z=0 
% based on bryan & Norman 1998   

eps=1e-6;

%Don't want more than 2 arugments altogether
numvarargs=length(varargin);
if numvarargs>2
    error('deltavir: too many inputs, Only 3 are allowed');
end

%set defualt omega_m and omega_lambda at z=0
defvals={0.3 ; 0.7};

%assign the optional values
defvals(1:numvarargs)=varargin;

% transfer to easy to use varaibles
[omm0, oml0]=defvals{:};

omm=omega_m(zred,'omm',omm0,'oml',oml0);

if(omm0+oml0-1.0)<eps %flat universe
    result=(18*pi^2+82.*(omm-1)-39.*(omm-1).^2)./omm;
elseif(omm0+oml0<1.0) %open universe 
    result=(18*pi^2+60.*(omm-1)-32.*(omm-1).^2)./omm;
else
    error('deltavir: illegal cosmology - Omega>1');
end
	     
