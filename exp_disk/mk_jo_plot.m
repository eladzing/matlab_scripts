%% recreate Jo's figure

%create SFR before stripping

%rm=0.01:0.01:100;


%ssfr0=zeros(size(cat.Ms));

rm=5.*cat.rd(ind);
%f%or i=1:length(cat.Ms)
ssfr0=ssfr_proxy(cat.Ms(ind),cat.Ms(ind).*cat.fg(ind),rm);%,cat.beta(i),rm));
%end

%create SFR after stripping
mh=10.^(14:0.05:14.6);
etap=1;
fc=0.15;
eta=0.01:0.01:50;
mmh=repmat(mh,length(cat.sigma(ind)),1);
mmh=reshape(mmh,numel(mmh),1);
sig=repmat(cat.sigma(ind),length(mh),1);
sig59=repmat(cat.sigma5090(ind),length(mh),1);
fg=repmat(cat.fg(ind),length(mh),1);
Ms=repmat(cat.Ms(ind),length(mh),1);

bet=repmat(cat.beta(ind),length(mh),1);
bt=repmat(cat.BT(ind),length(mh),1);
fc=fc.*ones(size(mmh));
etap=etap.*ones(size(mmh));
rm=repmat(rm,length(mh),1);
mstrip=zeros(size(sig));
ssfr=zeros(size(mstrip));
pad=7;

pval=rps_factor_expdisk('sigma_s',sig,'fd',fg,...
    'fc',fc,'mc',mmh,'etap',etap);

for i=1:length(mstrip)
    
    f1=disk_force_reduced(eta,'beta',bet(i),'fg',fg(i),...
        'BT',bt(i));
    [f1Max, id]=max(f1);
    
    mf1=exp_disk_mass(eta,bet(i));
    
    
    if pval(i)<=f1Max
        ind=find(f1>pval(i),1,'last');
        
        l2=min(ind+pad,length(f1));
        l1=max(1,ind-pad);
        
        ll=l1:l2;
        
        mstrip(i)=interp1(f1(ll),mf1(ll),pval(i),'cubic');
        
    else
        mstrip(i)=0;
    end
    
    
end


ssfr=ssfr_proxy(Ms,mstrip.*Ms.*fg,rm);%


[bir, binsize, xxlim,yylim]= histogram2d(log10(mmh),log10(sig59),ssfr,'xxlim',[14 14.6],'yylim',[5  8],'bins',[12 12]);
bird=bir(:,:,1)./bir(:,:,2);
xx=xxlim(1)+0.5*binsize(1):binsize(1):xxlim(2)-0.5*binsize(1);
yy=yylim(1)+0.5*binsize(2):binsize(2):yylim(2)-0.5*binsize(2);

figure
imagesc(xx,yy,log10(bird))
set(gca,'Ydir','normal','Fontsize',14)

bar=colorbar;
set(bar,'Fontsize',14)
set(get(bar,'Title'),'String','$\log(\mathrm{SSFR}$','Fontsize',14,'Interpreter','latex');
xlabelmine('$ \log(M_{halo}) \, [\mathrm{M_\odot}]$',14)
ylabelmine('$\Sigma_{50-90}\,[\mathrm{M_\odot/kpc^2}]$',14)


