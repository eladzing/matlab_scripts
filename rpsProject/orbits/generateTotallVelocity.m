function res = generateTotallVelocity(nOrb,hostMass,massRatio,varargin)
% GENERATETOTALVELOCITY - Generate the total velocity of an Orbit
%   Using the distributions of Jiang et al. 2015, generate random values of
%   the total velocity component (in units of Vvir).
%   nOrb - the number of / array size of values to generate
%   massRatio - mass ratio of satellite to host  - must be constant;
%   hostMass - mass of host halo.
%   values of the radial component are between 0 and 1;

% final halo mass:
mf=hostMass*(1+massRatio);

%% find appropriate B Value
if mf>=1e12 && mf<1e13
    % if massRatio>=0.0001 && massRatio<0.005
    if massRatio<0.005
        gamma=0.109;
        sigma=0.077;
        mu=1.22;
    elseif massRatio>=0.005 && massRatio<0.05
        gamma=0.098;
        sigma=0.073;
        mu=1.181;
    elseif massRatio>=0.05 && massRatio<0.5
        gamma=0.071;
        sigma=0.091;
        mu=1.1;
    else
        error('GENERATERADIALVELOCITY - mass ratio out of range: %s',...
            num2str(massRatio));
    end
elseif mf>=1e13 && mf<1e14
    if massRatio>=0.0001 && massRatio<0.005
        gamma=0.114;
        sigma=0.094;
        mu=1.231;
    elseif massRatio>=0.005 && massRatio<0.05
        gamma=0.087;
        sigma=0.083;
        mu=1.201;
    elseif massRatio>=0.05 && massRatio<0.5
        gamma=0.03;
        sigma=0.139;
        mu=1.1;
    else
        error('GENERATERADIALVELOCITY - mass ratio out of range: %s',...
            num2str(massRatio));
    end
elseif mf>=1e14
    if massRatio>=0.0001 && massRatio<0.005
        gamma=0.11;
        sigma=0.072;
        mu=1.254;
    elseif massRatio>=0.005 && massRatio<0.05
        gamma=0.05;
        sigma=0.118;
        mu=1.236;
    elseif massRatio>=0.05 && massRatio<0.5
        gamma=-0.012;
        sigma=0.187;
        mu=1.084;
    else
        error('GENERATERADIALVELOCITY - mass ratio out of range: %s',...
            num2str(massRatio));
    end
else
    error('GENERATERADIALVELOCITY - host mass out of range: %s',...
        num2str(hostMass));
end


i=1;
while i<=length(varargin)
    switch(lower(varargin{i}))
        case 'gamma'
            i=i+1;
            gamma=varargin{i};
            case 'sigma'
            i=i+1;
            sigma=varargin{i};
            case 'mu'
            i=i+1;
            mu=varargin{i};
        otherwise
            error('GENERATERADIALVELOCITY - unknown argument: %s',varargin{i});
    end
    i=i+1;
end


%% create cdf

vv=0:0.001:3;  % Vr/Vvir

xx=-1.5:0.001:1.5;
xx=xx+mu;
pg=1/(sqrt(2*pi)*sigma).*exp(-(xx-mu).^2./(2*sigma^2));

pdf=zeros(size(vv));
cdf=pdf;
for i=1:length(vv)
    pl=gamma./(pi.*((vv(i)-xx).^2+gamma^2));
    
    pdf(i)=trapz(xx,pg.*pl);
    if i>11
        cdf(i)=trapz(vv(1:i),pdf(1:i));
    else
        cdf(i)=pdf(i);
    end
    
end

cdf=cdf./trapz(vv,pdf);



%% generate random values
[cdf, mask] = unique(cdf);
vv=vv(mask);

rv=rand(1,nOrb);
resTemp=interp1(cdf,vv,rv);
resTemp=resTemp(resTemp>0);

df=nOrb-length(resTemp);
cnt=0;
while df>0 || cnt>20
    cnt=cnt+1;
    
    rv=rand(1,df);
    resT=interp1(cdf,vv,rv);
    resT=resT(resT>0);
        
    resTemp=cat(2,resTemp,resT);
        
    df=nOrb-length(resTemp);
end

res=resTemp;





