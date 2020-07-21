function [r50 sig50 sig5090]=sigmaObserved(ms,rd,fb,xi)

eta=logspace(-3,1.5,1e4);

%
f50=0.5.*(1+fb);% half mass
f90=0.9.*(1+fb);% half mass

eta50=zeros(size(ms));
eta90=zeros(size(ms));
for i=1:length(ms)

    ff=1-exp(-eta).*(1+eta)+fb(i).*xi(i).^2.*eta.^2./(1+xi(i).*eta).^2;

    ff50=ff-f50(i);
    ff90=ff-f90(i);
    eta50(i)=interp1(ff50,eta,0);
    eta90(i)=interp1(ff90,eta,0);
end



r50=eta50.*rd;
r90=eta90.*rd;


sig50=0.5.*ms.*(1+fb)./(pi.*r50.^2);
sig5090=0.4.*ms.*(1+fb)./(pi.*r90.^2-pi.*r50.^2);
