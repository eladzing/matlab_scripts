function res = generate_projectionFac(nn)
%GENERATE_PFAC generate nn randomly selected projection factors


phiRange=[0 2*pi];
psiRange=phiRange;
zRange=[-1 1];
%thetRange=[-pi/2 pi/2];

phi=phiRange(1)+diff(phiRange).*rand(nn,1);
psi=psiRange(1)+diff(psiRange).*rand(nn,1);
zz=zRange(1)+diff(zRange).*rand(nn,1);

theta=acos(zz);

%theta=thetRange(1)+diff(thetRange).*rand(nn,1);
res=zeros(1,length(phi));
for i=1:nn
    res(i)=projectionFac(phi(i),theta(i),psi(i));
end

end

