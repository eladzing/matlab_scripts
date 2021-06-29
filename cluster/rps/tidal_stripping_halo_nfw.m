%% solve rps for NFW model
zred=0;
mc=1e15;
ms=1e12;

mu=(ms/mc);
muu=mu^(1/3);

%% find relevent NFW parameters:
cc=cvir_Mvir(mc,zred);
cs=cvir_Mvir(ms,zred);

ac=nfwA(1,cc);
as=nfwA(1,cs);

%% set rc and ls
fg=0.1;
fc=0.1;
alfa=1;


rcc=0.01:0.001:3;
ls=0.01:0.0001:2;

lls=zeros(size(rcc));
for i=length(rcc):-1:1
    
    rc=rcc(i);
    %rhs=(alfa*fg*(1+fg)/fc)*(ms./mc)^(2/3)*(ac/as).^2.*rc.*(cc^-1+rc).^2./nfwA(rc,cc);
    
    ft1=-nfwA(rc+ls.*muu,cc)./(rc+ls.*muu).^2+nfwA(rc,cc)./rc.^2;
    fs1=-muu.*(ac./as).*nfwA(ls,cs)./ls.^2;
    
    ff1=ft1+fs1;
    ii1=find(ff1>0,1,'first');
    
    ft2=-nfwA(rc-ls.*muu,cc)./(rc-ls.*muu).^2+nfwA(rc,cc)./rc.^2;
    fs2=muu.*(ac./as).*nfwA(ls,cs)./ls.^2;
    
    ff2=ft2+fs2;
    ii2=find(ff2<0,1,'first');
    %%rhs=(fg*(1+0)/fc/alfa)*(ms./mc)^(2/3)*(ac/as.^2).*rc.*(cc^-1+rc).^2;%./nfwA(rc,cc);
    
    %%ff=lhs-rhs;
    
    %%ii=find(ff>1,1,'last');
    if isempty(ii1)
        lt(i,1)=inf;
    else
        lt(i,1)=ls(ii1);
    end
    if isempty(ii2)
        lt(i,2)=inf;
    else
        lt(i,2)=ls(ii2);
    end
end
%mms=nfwA(lls,cs)/nfwA(1,cs);



