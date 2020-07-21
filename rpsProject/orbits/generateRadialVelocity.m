function res = radialVelocity_pdf(nOrb,hostMass,massRatio)
% GENERATERADIALVELOCITY - Generate the radial velocity component of an
% Orbit
%   Using the distributions of Jiang et al. 2015, generate random values of
%   the radial velocity component (in units of Vvir).
%   nOrb - the number of / array size of values to generate
%   massRatio - mass ratio of satellite to host  - must be constant;
%   hostMass - mass of host halo.
%   values of the radial component are between 0 and 1;

% final halo mass:
mf=hostMass*(1+massRatio);

%% find appropriate B Value
if mf>=1e12 && mf<1e13
    %if massRatio>=0.0001 && massRatio<0.005
    if  massRatio<0.005
        bVal=0.049;
    elseif massRatio>=0.005 && massRatio<0.05
        bVal=1.044;
    elseif massRatio>=0.05 && massRatio<0.5
        bVal=2.878;
    else
        error('GENERATERADIALVELOCITY - mass ratio out of range: %s',...
            num2str(massRatio));
    end
elseif mf>=1e13 && mf<1e14
    if massRatio>=0.0001 && massRatio<0.005
        bVal=0.548;
    elseif massRatio>=0.005 && massRatio<0.05
        bVal=1.535;
    elseif massRatio>=0.05 && massRatio<0.5
        bVal=3.946;
    else
        error('GENERATERADIALVELOCITY - mass ratio out of range: %s',...
            num2str(massRatio));
    end
elseif mf>=1e14
    if massRatio>=0.0001 && massRatio<0.005
        bVal=1.229;
    elseif massRatio>=0.005 && massRatio<0.05
        bVal=3.396;
    elseif massRatio>=0.05 && massRatio<0.5
        bVal=2.982;
    else
        error('GENERATERADIALVELOCITY - mass ratio out of range: %s',...
            num2str(massRatio));
    end
else
    error('GENERATERADIALVELOCITY - host mass out of range: %s',...
        num2str(hostMass));
end

%% create cdf

vv=0:0.001:1;  % Vr/Vvir
cdf=(exp(bVal.*vv)-bVal.*vv-1)./(exp(bVal)-bVal-1);

%% generate random values 
rv=rand(1,nOrb);
res=interp1(cdf,vv,rv);





