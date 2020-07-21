function res=vcSquare(eta,fg,beta,fb,xi)

J=bfunc(eta,1)+fg.*beta.^3.*bfunc(eta,beta)+2.*fb.*xi.^2./(eta.*(1+xi.*eta).^2);


res=eta.^2.*J;

