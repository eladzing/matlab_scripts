function cv=find_cee(ros,dvir,rou)

rofac=1.989e33./(1e6.*3.0856e18).^3;

fac=3.*ros.*rofac./(rou.*dvir);
c=1:1e-5:100;
a=fac.*(log(1+c)-c./(c+1))./c.^3;
i1=find(a<=1,1,'first');
i2=find(a>1,1,'last');

d1=abs(a(i1)-1);
d2=abs(a(i2)-1);

if d1<d2
    cv=c(i1);
else
    cv=c(i2);
end

