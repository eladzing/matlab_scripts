%% explore_rps_disk 

etap=1;
mc=5e14;
fc=0.15;


%disk 
%sig=10^8.*ones(size(mcMat));
Ms=1e10;
fg=0.1:0.05:0.25;
beta=0.1:0.1:1;
BT=0;
mv=1e11;
cv=cvir_Mvir(mv,0);
lambda=0.04;


%pv=zeros(length(etap),length(mc));

[fgMat,betapMat]=meshgrid(fg,beta);



pval=rps_factor_expdisk('sigma_s',sig,'fd',fg,...
    'fc',fc,'mc',mcMat,'etap',etapMat);


eta=0.01:0.01:10;
f1=disk_force_reduced(eta,'beta',beta(1),'fg',fg(1),'BT',BT(1));
[f1Max, id]=max(f1);
%f1(1:id)=f1Max;
mf1=exp_disk_mass(eta,beta(1));

ms=zeros(size(etapMat));
tic
for i=1:length(etap)
    for j=1:length(mc)
        if pval(i,j)<=f1Max
              ms(i,j)=interp1(f1,mf1,pval(i,j));
        else
            ms(i,j)=0;
        end
              
    end
end
toc
%ms(pval>max(f1))=0;

ms=(1-ms);

figure 
imagesc(mc,etap,log10(ms))
set(gca,'Ydir','normal')
colorbar
