function cvir=cvir_Mvir_200(Mv,zred,type,varargin)
% function to realize the cvir-Mvir relation as presented in
% Dutton and  Maccio  2014
% Mv should be in solarmass NOT solarmass h^-1!! 
% the fitting relation is for Mv h^-1 so the value of h is an
% optional argument

if nargin<2
    error('cvir_Mvir: two few arguments')
end

%% formula parameters
switch type
    case{'200'}
a=0.52+(0.905 - 0.520)*exp(-0.617.*zred.^1.21);
b=-0.101+0.026.*zred;
    case{'vir'}
a=0.537+(1.025 - 0.537)*exp(-0.718.*zred.^1.08);
b=-0.097+0.024.*zred;
    otherwise
        error('CVIR_MVIR_200 - Illegal type (200,vir): %s',type)
end

% w=0.029;
% m=0.097;
% alpha=-110.001;
% beta=2469.72;
% gamma=16.885;
% 
% a=w*zred-m;
% b=alpha/(zred+gamma)+beta/(zred+gamma)^2;



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

%cvir=10.^(a.*log10(mv)+b);

cvir=10.^(a+b.*log10(mv./1e12));


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




