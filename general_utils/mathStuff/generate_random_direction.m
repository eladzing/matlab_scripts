function res = generate_random_direction(nn)
%GENERATE_RANDOM_DIRECTION generate a random direction characterizeds by
% angles 
%   Detailed explanation goes here

phiRange=[0 2*pi];
psiRange=phiRange;
zRange=[-1 1];
%thetRange=[-pi/2 pi/2];

res.phi=phiRange(1)+diff(phiRange).*rand(nn,1);
res.psi=psiRange(1)+diff(psiRange).*rand(nn,1);
zz=zRange(1)+diff(zRange).*rand(nn,1);

res.theta=acos(zz);



end

