%% plot stripping catalog

%load('mockCat_1000.mat')

etap=0.2;
mc=1e15;
fc=0.15;

pval=rps_factor_nfw_expdisk('sigma_s',cata.sigma,'fd',cata.fg,...
    'fc',fc,'mc',mc,'etap',etap,'alpha',1);

eta=0.01:0.01:50;
mstrip=zeros(size(cata.sigma));
rstrip=mstrip;
pad=7;

% [rvir,~,~,~]=calculate_virials('mv',cata.Mv);

tic
for i=1:length(mstrip)
    %^if ~mask(i)
    %    continue
    %end
    fd=disk_force_reduced(eta,'beta',cata.beta(i),'fg',cata.fg(i),...
        'BT',cata.BT(i));
    % create mhalo vs. eta from raw data
    mh0=interp1(cata.raw.rr(i,:)./cata.rd(i),cata.raw.mDM(i,:)./cata.Ms(i),eta,'PCHIP');
    
    
    fh=2.*mh0./eta.^2.*cata.beta(i).^2.*exp(-cata.beta(i).*eta);
    % test mdm contribution
    
    % fht=halo_accel(eta,cata.Mv(i)/cata.Ms(i),cata.rd(i)./(1e3*rvir(i)),'nfw','cv',cata.cv(i)).*...
    %     expdisk_density(eta,'gas','beta',cata.beta(i),'fg',1.0);
    
    f1=fd+fh;
    
    [f1Max, id]=max(f1);
    
    mf1=exp_disk_mass(eta,cata.beta(i));
    
    if(any(isnan(f1)) || any(isnan(mf1)) ||any(isnan(pval)))
        i;
    end
    
    if pval(i)<=f1Max
        ind=find(f1>pval(i),1,'last');
        
        l2=min(ind+pad,length(f1));
        l1=max(1,ind-pad);
        
        ll=l1:l2;
        
        
        mstrip(i)=interp1(f1(ll),mf1(ll),pval(i),'PCHIP');
        rstrip(i)=interp1(f1(ll),eta(ll),pval(i),'PCHIP');
    else
        mstrip(i)=0;
        rstrip(i)=0;
    end
    
    
    
    
end
toc

mstrip=(1-mstrip).*100;

%% calculate SSFR

knob=0;
A=(2.5+knob.*0.7).*1e-4;
alfa=1.4+knob.*0.15;

%sfrfac=A*10^(-6*alfa)/alfa^2;
%ssfr0=sfrfac*(cata.sigma.*cata.fg.*cata.beta.^2).^alfa...
%    .*(cata.rd./cata.beta).^2./(cata.Ms.*(1+cata.fb));
sfrFrac=(1-exp(-alfa.*cata.beta.*rstrip).*(1+alfa.*cata.beta.*rstrip));

sfrFrac=100.*(1-sfrFrac);  % reduction in SFR

% figure
% subplot(3,1,1:2)
% plot(log10(cata.sigma),mstrip,'.')
% xlim([6 11])
% grid
% ylabelmine('$M_{strip}\,[\%]$',14)
% set(gca,'Fontsize',14,'box','on');
%
% subplot(3,1,3)
% [nh xo]=hist(log10(cata.sigma),20);
% bar(xo,nh./length(mstrip).*100);
% set(gca,'Fontsize',14,'box','on');
% xlim([6 11])
% xlabelmine('$\Sigma_s\,[\mathrm{M_\odot/kpc^{2}}]$',14)
% ylabelmine('$\%$',14)
hh=fspecial('gaussian',10,3);
totMask=cata.lambdaMask & cata.qmin>=1;

%% plot for stripped mass
re=find_re(cata.fb,cata.xi);
sig=0.5.*cata.Ms.*(1+cata.fb)./(pi.*re'.^2);

[bir, binsize, xxlim,yylim]= histogram2d(log10(cata.sigma(totMask)),mstrip(totMask),ones(size(mstrip(totMask))),'xxlim',[6 11],'yylim',[0 100]);
bird=bir(:,:,1)./length(mstrip(totMask)).*100;
bss=imfilter(bird,hh);
%bs=imfilter(bird,hh,'full');
bs=bss./sum(sum(bss)).*100;

%figure
bb=zeros(size(bird));

for i=1:size(bb,1)
    for j=1:size(bb,2)
        maskTmp=bs>=bs(i,j);
        bb(i,j)=sum(sum(bs(maskTmp)));
    end
end
% contourf(bb)
% contourf(bb,10:10:100)
% contourf(bb,1:10:100)
xx=xxlim(1)+0.5*binsize(1):binsize(1):xxlim(2)-0.5*binsize(1);
yy=yylim(1)+0.5*binsize(2):binsize(2):yylim(2)-0.5*binsize(2);
lc=[0 10 20 30 40 50 60 70 80 90 100] ;%    0:10:100;%lc(end+1)=99.9;
bb(bb==max(max(bb)))=100;

ylab='Stripped Mass $[\%]$';
createfigure_stripped_catalog_with_cumsum(xx, yy, bb, cumsum(sum(bird,2)),ylab)

%% plot for SFR reduction

% [bir, binsize, xxlim,yylim]= histogram2d(log10(sig(totMask)),sfrFrac(totMask),ones(size(mstrip(totMask))),'xxlim',[6 11],'yylim',[0 100]);
% bird=bir(:,:,1)./length(mstrip(totMask)).*100;
% bss=imfilter(bird,hh);
% %bs=imfilter(bird,hh,'full');
% bs=bss./sum(sum(bss)).*100;
% 
% %figure
% bb=zeros(size(bird));
% 
% for i=1:size(bb,1)
%     for j=1:size(bb,2)
%         maskTmp=bs>=bs(i,j);
%         bb(i,j)=sum(sum(bs(maskTmp)));
%     end
% end
% % contourf(bb)
% % contourf(bb,10:10:100)
% % contourf(bb,1:10:100)
% xx=xxlim(1)+0.5*binsize(1):binsize(1):xxlim(2)-0.5*binsize(1);
% yy=yylim(1)+0.5*binsize(2):binsize(2):yylim(2)-0.5*binsize(2);
% lc=[0 10 20 30 40 50 60 70 80 90 100] ;%    0:10:100;%lc(end+1)=99.9;
% bb(bb==max(max(bb)))=100;
% 
% ylab='SFR Reduction $[\%]$';
% createfigure_stripped_catalog_with_cumsum(xx, yy, bb, cumsum(sum(bird,2)),ylab)
% 
% % [C,h]= contourf(xx,yy,bb,lc);
% % clabel(C,h)
% % bar=colorbar;
% % set(bar','Fontsize',14)
% % set(get(bar,'Title'),'String','Population $[\%]$','Fontsize',14,'Interpreter','latex');
% % set(gca,'Ydir','normal','Fontsize',14)
% % xlabelmine('$ \log(\Sigma_s) \, [\mathrm{M_\odot/kpc^2}]$',14)
% % ylabelmine('Stripped Mass $[\%]$',14)
