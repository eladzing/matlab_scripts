%% solve rps for NFW model 
zred=0;
mc=1e15;
ms=1e12;

%% find relevent NFW parameters: 
cc=cvir_Mvir(mc,zred);
cs=cvir_Mvir(ms,zred);

ac=nfwA(1,cc);
as=nfwA(1,cs);

%% set rc and ls 
fg=0.1;
fc=0.1;
alfa=1;


rcc=0.5:0.001:3;
ls=0.01:0.0001:2;
lhs=(ls.^2.*(cs.^-1+ls).^2)./nfwA(ls,cs);
lls=zeros(size(rcc));
for i=length(rcc):-1:1

rc=rcc(i);    
%rhs=(alfa*fg*(1+fg)/fc)*(ms./mc)^(2/3)*(ac/as).^2.*rc.*(cc^-1+rc).^2./nfwA(rc,cc);
rhs=(fg*(1+0)/fc/alfa)*(ms./mc)^(2/3)*(ac/as.^2).*rc.*(cc^-1+rc).^2;%./nfwA(rc,cc);

ff=lhs-rhs;

ii=find(ff>0,1,'first');
lls(i)=ls(ii);
end

mms=nfwA(lls,cs)/nfwA(1,cs);



