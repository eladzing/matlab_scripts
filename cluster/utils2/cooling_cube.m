function coolin=cooling_cube(boxx,unitArg)
%% calculates cooling rate in each cell 
% in units of   

unitType='custom';

if exist('unitArg','var')
    unitType=unitArg;
end
    


units;
[tlam, zlam, lambda]=read_lambda; %%read interpolation data for cooling function (in log10) 
%vol=(boxx./hub./NCELL).^3;

% read datacubes
n=RHOG(boxx);
tt=log10(T(boxx));
zmet=ZIa(boxx)+ZII(boxx);

%ignore cells which are too cold
mask=tt>4.0;

% calculate cooling function
lamb=zeros(size(tt));
len=size(zmet,1);
for ii=1:len    
    zm=squeeze(zmet(ii,:,:));
    tm=squeeze(tt(ii,:,:));  
    lamb(ii,:,:)=interp2(zlam,tlam,lambda,zm,tm,'spline');
    clear tm zm   
end

lamb22=10.^((lamb+22).*mask).*mask; 
clear mask lamb

coolin=n.^2.*lamb22;

switch unitType
    case 'custom' % units of M_sun*(km/sec)^2 per Mpc^3
        factor=Units.factors.f_cl;
    case 'cgs' 
        factor=(Units.xi/(Units.muMass*Units.mp))^2*(Units.Ms/Units.Mpc^3)^2*1e-22/1e42; %ergs/sec per cm^3 in units of 10^42 erg/sec
    case 'cgs/Mpc' 
        factor=(Units.xi/(Units.muMass*Units.mp))^2*(Units.Ms^2/Units.Mpc^3)*1e-22/1e42; %ergs/sec per Mpc^3  in units of 10^42 erg/sec
    otherwise
        error('COOLING_CUBE - Illegal unit type %s',unitArg);
end
           
coolin=coolin.*factor; %cooling rate per unit volume