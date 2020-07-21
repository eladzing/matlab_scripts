%% explore_rps_disk

sig=8:0.02:10;
mc=13:0.02:16;
fc=0.15;
etap=0.5;
%pv=zeros(length(etap),length(mc));

[mcMat,sigMat]=meshgrid(10.^mc,10.^sig);

%disk
%sig=5.*10^8.*ones(size(mcMat));
fg=0.2;
beta=0.33;
BT=0;

pval=rps_factor_expdisk('sigma_s',sigMat,'fd',fg,...
    'fc',fc,'mc',mcMat,'etap',etap);


eta=0.01:0.01:10/beta;
f1=disk_force_reduced(eta,'beta',beta(1),'fg',fg(1),'BT',BT(1));
[f1Max, id]=max(f1);
%f1(1:id)=f1Max;
mf1=exp_disk_mass(eta,beta(1));

ms=zeros(size(sigMat));

pad=5;

tic
for i=1:length(sig)
    for j=1:length(mc)
        if pval(i,j)<=f1Max
            ind=find(f1>pval(i,j),1,'last');
         
        l2=min(ind+pad,length(f1));
        l1=max(1,ind-pad);
        
        ll=l1:l2;
        
        ms(i,j)=interp1(f1(ll),mf1(ll),pval(i,j),'cubic');
        
        %
         %   ms(i,j)=interp1(f1,mf1,pval(i,j),'cubic');
        else
            ms(i,j)=0;
        end
        
    end
end
toc
%ms(pval>max(f1))=0;

ms=100.*(1-ms);

bartag='$\%M_{strip}$';
load('MyColormaps','avijet');

figure
imagesc(mc,sig,(ms),[0 100])
xlabelmine('$\log(M_c)\,[\mathrm{M_\odot}]$',14)
ylabelmine('$\log(\Sigma_s)\,[\mathrm{M_\odot/kpc^2}]$',14)
set(gca,'Ydir','normal','Fontsize',14)
colormap(avijet);
bar=colorbar;
set(get(bar,'Title'),'String',bartag,'Fontsize',12,'Interpreter','latex');
set(gcf,'Colormap',avijet);
titlemine(sprintf('$\\eta_p=%s,\\,\\beta=%s$',num2str(etap),num2str(beta)));

