function for_amri

xx=[0 pi];  %linspace(0,pi,100);
yy=[400 100]; % (200 -350)./pi.*xx+350;

soloint=bvpinit(xx,yy);

sol=bvp4c(@amriODE,@bcs,soloint);

end


function dydx=amriODE(x,y)

a=1;
b=10;
c=10^-9;

dydx=a.*y.^4-b.*cos(x)-c;

end

function res=bcs(ya,yb)

res=[ya(1) yb(1)];

end