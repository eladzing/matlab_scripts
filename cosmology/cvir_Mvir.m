function cvir=cvir_Mvir(Mv,zred,varargin)
% function to realize the cvir-Mvir relation as presented in
% Munoz-Cuartas et al. (2010)
% Mv should be in solarmass NOT solarmass h^-1 !! 
% the fitting relation is for Mv h^-1 so the value of h is an
% optional argument

if nargin<2
    error('cvir_Mvir: two few arguments')
end

%% formula parameters
w=0.029;
m=0.097;
alpha=-110.001;
beta=2469.72;
gamma=16.885;

a=w*zred-m;
b=alpha/(zred+gamma)+beta/(zred+gamma)^2;



%% random parameters
sigfac=0.05;
distType='normal';
randomFlag=false;

hub=0.7;

i=1;
while i<=length(varargin)
    switch(varargin{i})
        case {'h','hub'}
            i=i+1;
            hub=varargin{i};
        case {'random'}
            randomFlag=true;
        case {'type','distribution'}
            i=i+1;
            distType=varargin{i};
            randomFlag=true;
        case {'sigma','sig'}
            i=i+1;
            sigfac=varargin{i};
        otherwise
            error('cvir_Mvir: Illegal argument: %s',varargin{i})
    end
    i=i+1;
end

mv=Mv.*hub;

cvir=10.^(a.*log10(mv)+b);

if randomFlag
    mu=cvir;
    sig=sigfac.*mu;
    
    switch distType
        case {'normal','gaussian'}
            cvir=normrnd(mu,sig);
        case {'uniform','constant'}
            cvir=(mu-sig./2)+sig.*rand(size(mu));
        case {'lognormal'}
            cvir=lognrnd(mu,sig);
        otherwise
            error('cvir_Mvir: Illegal probablity distribution: %s',distType)
    end
    
    
end




