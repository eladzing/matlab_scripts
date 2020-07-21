function result = mbird()

env;

t=T(1);
ro=RHOG(1);
xx=
weit=


tro1=hist2d(t,ro,xx,weit);

t=T(2);
ro=RHOG(2);
xx=
weit=

l1=round(size(t)*0.25)+1;
l2=round(size(t)*0.75);
weit(l1(1):l2(1),l1(2):l2(2),l1(3):l2(3))=0;
xx(l1(1):l2(1),l1(2):l2(2),l1(3):l2(3))=0;



tro2=hist2d(t,ro,xx,weit);
