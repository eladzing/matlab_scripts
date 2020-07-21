%% plot stripping dependence on cluster position 

ss=7:0.01:10;
sigma=10.^(ss);
etap=0.01:0.01:3;


len=length(sigma)*length(etap);

sig=reshape(repmat(sigma,[length(etap) 1]),[len 1]);
etp=repmat(etap,[1 length(sigma)])';

Mc=1e15;
fg=0.1;
fc=0.1;
beta=0.5;
BT=0.05;

pval=rps_factor_expdisk('etap',etp,'sigma',sig,...
'fd',fg,'fc',fc,'Mc',Mc);

r=10.^(-3:0.001:2);%0.001:0.001:100; 
%r1=0.0001:0.0001:100; 

f=disk_force_reduced(r,'fg',fg,'beta',beta,'BT',BT);
%f1=disk_force_reduced(r1,'fg',fg,'beta',beta,'BT',BT);
[fmax, imax]=max(f);
r=r(imax:end);
f=f(imax:end);
%f(1:imax)=fmax;
mStr=zeros(size(pval));
ind=find(pval<=fmax);

for i=ind %1:length(pval)
       
    %if fmax>pval(i)
        rStr=interp1(f,r,pval(i));
        mStr(i)=exp_disk_mass(rStr,beta);
    %else
    %    mStr(i)=0;
    %end
    
end

mStr=reshape(mStr,[length(etap) length(sigma)]);
mStr=(1-mStr).*100;

figure 
imagesc(ss,etap,mStr);
set(gca,'Ydir','normal','Fontsize',14)
bar=colorbar;
set(bar','Fontsize',14)
set(get(bar,'Title'),'String','$M_{stripped}\,[\%]$','Fontsize',14,'Interpreter','latex');
xlabelmine('$\log(\Sigma_s)\,[\mathrm{M_\odot/kpc^2}]$')
ylabelmine('$r/R_c$')