
env

s=S(1);
ls=log10(s.*7415.967);
locs=find(ls<lim2 & ls>lim1);
[i,j,k]=ind2sub(size(ls),locs);
x1=1*(i/256-0.5);
y1=1*(j/256-0.5);
z1=1*(k/256-0.5);
clear i j k 
%

s=S(2);
ls=log10(s.*7415.967);
locs=find(ls<lim2 & ls>lim1);
[i,j,k]=ind2sub(size(ls),locs);
lo=~((i>64 & i<193) & (j>64 & j<193) & (k>64 & k<193));

i=i.*lo;
j=j.*lo;
k=k.*lo;
x2=2*(i/256-0.5);
y2=2*(j/256-0.5);
z2=2*(k/256-0.5);
clear i j k 
%
%

s=S(4);
ls=log10(s.*7415.967);
locs=find(ls<lim2 & ls>lim1);
[i,j,k]=ind2sub(size(ls),locs);
lo=~((i>64 & i<193) & (j>64 & j<193) & (k>64 & k<193));
i=i.*lo;
j=j.*lo;
k=k.*lo;
x4=4*(i/256-0.5);
y4=4*(j/256-0.5);
z4=4*(k/256-0.5);
clear i j k
%
s=S(8);
ls=log10(s.*7415.967);
locs=find(ls<lim2 & ls>lim1);
[i,j,k]=ind2sub(size(ls),locs);
lo=~((i>64 & i<193) & (j>64 & j<193) & (k>64 & k<193));
i=i.*lo;
%lj=(j<128 | j>192);
j=j.*lo;
%lk=(k<128 | k>192);
k=k.*lo;
x8=8*(i/256-0.5);
y8=8*(j/256-0.5);
z8=8*(k/256-0.5);
clear i j k
%

x=cat(1,x1,x2,x4,x8);
y=cat(1,y1,y2,y4,y8);
z=cat(1,z1,z2,z4,z8);

plot3(x,y,z,'.')