cvir1=zeros(16,3);
cvir2=zeros(16,3);

for i=1:16
[rp lrp rot lrot rog lrog rodm lrodm mto lmto rv mv clname]=extract_rhoprof(i);
rv=r200(i);
ft=rho_fit(lrp,lrot,lrodm,lrog,rv);
cofs1=coeffvalues(ft(1).cft);
%cofs2=coeffvalues(ft(2).cft);
conf1=confint(ft(1).cft);
%conf2=confint(ft(2).cft);
cvir1(i,1)=cofs1(2);
cvir1(i,2)=conf1(2,2);
cvir1(i,3)=conf1(1,2);
%cvir2(i,1)=rv./cofs2(2);
%cvir2(i,2)=rv./conf2(2,2);
%cvir2(i,3)=rv./conf2(1,2);


mvir(i)=m200(i);
%disp(cv)
clname

close all
end
%cvir
%close all;figure;scatter(mvir,cvir);
%close all;
%figure;scatter(log10(mvir),cvir);