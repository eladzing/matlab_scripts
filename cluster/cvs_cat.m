cvir1=zeros(16,3);
%cvir2=zeros(16,3);
delv=336.9;
rou=2.776e11*0.7^2*0.3*2e33/(1e6*3.0856e18)^3; %1e-30;%2.76e-30;
disp(rou)

for i=1:16
[rp lrp rot lrot rog lrog rodm lrodm mto lmto rv mv clname]=extract_rhoprof(i);
ft=rho_fit(lrp,lrot,lrodm,lrog,rv);
cofs1=coeffvalues(ft(1).cft);
conf1=confint(ft(1).cft);
ros=10.^cofs1(1);rs=cofs1(2);
rosh=10.^conf1(2,1);rosl=10^conf1(1,1);
rsh=conf1(2,2);rsl=conf1(1,2);

%cofs2=coeffvalues(ft(2).cft);
%conf1=confint(ft(1).cft);
%conf2=confint(ft(2).cft);
cvir1(i,1)=find_cee(ros,delv,rou);%rv./cofs1(2);
cvir1(i,2)=find_cee(rosl,delv,rou);%rv./conf1(2,2);
cvir1(i,3)=find_cee(rosh,delv,rou);%rv./conf1(1,2);
%cvir2(i,1)=rv./cofs2(2);
%cvir2(i,2)=rv./conf2(2,2);
%cvir2(i,3)=rv./conf2(1,2);

cv=cvir1(i,1);
mvir(i,1)=4.*pi.*ros.*rs.^3.*(log(1+cv)-cv./(cv+1));

mvir(i,2)=mv;

disp(cv)
clname

%pause
close all
end
%cvir
%close all;figure;scatter(mvir,cvir);
%close all;
%figure;scatter(log10(mvir),cvir);