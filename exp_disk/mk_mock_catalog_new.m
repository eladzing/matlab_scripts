%% make mock catalog
units
%nSample=1;
ngal=100;

% set parameter values
lms=9.*ones(ngal,1); %; +(11.5-9).*rand(ngal,1);
%lms=8.75+(9.5-8.75).*rand(ngal,1);
ms=10.^lms;

fg=0.05+(1-0.05).*rand(ngal,1);
beta=0.4+(1-0.4).*rand(ngal,1);
fb=0+(1-0).*rand(ngal,1);
xi=0.5+(4-0.5).*rand(ngal,1);
zeta=normrnd(0.23,0.071,size(ms));
zeta(zeta<0)=min(zeta(zeta>0));

Mv=halomass_from_moster(ms.*(1+fb),'noshow');
cv=cvir_Mvir(Mv,0,'random');
lambda=lambda_prime(ms);


%% calculate rd 
tic
res=rscale_mmw_array(ms,'fg',fg,'beta',beta,'fb',fb,'xi',xi,...
    'Mv',Mv,'cv',cv,'lambda',lambda,'noshow');
toc


% check stability  according to MMW98 and Toomre Q
%eta=0.01:0.01:50;
eps=zeros(size(ms));
qmin=eps;
qsub=eps;

for i=1:length(ms)
    ind=find(res.rr(i,:)./res.rd(i)<10);
    % eps a al Efstathiou et al 1982
    vmax=max(res.vc(i,ind)); % in km/sec
    eps(i)=vmax/sqrt(GG*ms(i)/res.rd(i));
    
    %Q Toomre
    Q=QToomre(res.rr(i,ind),res.vc(i,ind),ms(i),res.rd(i),1,zeta(i),GG);
    qmin(i)=min(Q);
    subInd=find(Q<1);
    if ~isempty(subInd)
        qsub(i)=diff(res.rr(i,[subInd(1)+1 subInd(end)+1]))./res.rd(i);
    end
    % loglog(res.rr(i,ind(1)+1:ind(end)-1)./res.rd(i),Q)
    %pause
end

%maskStab=eps>1.1;
%mask=maskStab & res.lambdaMask;
mask=true(size(ms));
%% build catalog;
cata.Ms=ms(mask);
cata.rd=res.rd(mask);
cata.Mv=res.Mv(mask);
cata.md=res.md(mask);
cata.BT=res.BT(mask);
cata.fb=fb(mask);
cata.fg=fg(mask);
cata.xi=xi(mask);
cata.cv=res.cv(mask);
cata.lambda=res.lambda(mask);
cata.sigma=ms(mask)./(2*pi*res.rd(mask).^2);
cata.sigmaeff=ms(mask)./(2*pi*res.rd(mask).^2).*sigma_factor('half');
cata.sigma5090=ms(mask)./(2*pi*res.rd(mask).^2).*sigma_factor('5090');
cata.beta=beta(mask);
cata.qmin=qmin(mask);
cata.qsub=qsub(mask);
cata.eps=eps(mask);
cata.lambdaMask=res.lambdaMask;
cata.raw=res;


%% older
% mms=repmat(ms,nSample,1);
% 
% 
% %beta=0.5.*ones(size(mms));
% beta=0.4+(1-0.4).*rand(size(mms)); % up to x2.5
% fg=0.05+(0.25-0.05).*rand(size(mms));
% fb=0.0+(0.25-0.0).*rand(size(mms));
% xi=0.8+(4-0.8).*rand(size(mms));
% zeta=normrnd(0.23,0.071,size(mms));
% zeta(zeta<0)=min(zeta(zeta>0));
% 
% 
% % find halo mass according to stellar to Halo mass
% % relation of Moster et al. 09
% Mv=halomass_from_moster(mms.*(1+fb),'noshow');
% 
% tic
% res=rscale_mmw_array(mms,'beta',beta,'fg',fg,'fb',fb,'xi',xi,...
%     'Mv',Mv,'noshow');
% toc


%% plot distributions

% hfh1=figure;
% [nr xr]=hist(res.rd,50);
% bar(xr,nr./ngal.*100)
% xlabelmine('$R_d\,[\mathrm{kpc}]$',14)
% ylabelmine('$\%$')
%
% hfh2=figure;
% [ns xs]=hist(log10(mms),50);
% bar(xs,ns./ngal.*100)
% xlabelmine('$\log(\Sigma_s) \,[\mathrm{M_\odot/kpc^2}]$',14)
% ylabelmine('$\%$')
%
% %% 2d histogram
% yylim=[0 10];
% [bir, binsize, xxlim,yylim]= histogram2d(log10(mms),res.rd,ones(size(mms)),'yylim',yylim);
% bird=bir(:,:,1)./ngal.*100;
%
%
% %hh=ones(10);
% hh=fspecial('gaussian',10,3);
% %hh=fspecial('average',10);
% bs=imfilter(bird(:,:,1),hh);
%
% md=linspace(xxlim(1),xxlim(2),201);
% rd=linspace(yylim(1),yylim(2),201);
% %% plot with contours
%  load('MyColormaps','avijet_bird');
%  bartag='Sample$\%$';
% %md=8:0.02:12;
% %rd=0:0.02:15;
%
% [mdMat,rdMat]=meshgrid(10.^md,rd);
%
% sigMat=mdMat./(2*pi.*rdMat.^2);
%
% hf1=figure;
% %[CB,hcB]=contour(md,rd,bs);
% %set(hcB,'ShowText','on','TextStep',get(hcB,'LevelStep')*2)
% contourf(md,rd,bs,0:0.001:0.05,'LineColor','none');colormap(avijet_bird)
% caxis([min(min(bs)) max(max(bs))]);
% hold on
% [C1,hc1]=contour(md,rd,log10(sigMat.*sigma_factor('half')),4:0.5:10);
% set(hc1,'ShowText','on','TextStep',get(hc1,'LevelStep')*2)
%
% bar=colorbar;
% set(get(bar,'Title'),'String',bartag,'Fontsize',12,'Interpreter','latex');
% set(gca,'Fontsize',14,'box','on','Ydir','Normal')
% xlabelmine('$M_s\,[\mathrm{M_\odot}]$',14)
% ylabelmine('$R_d\,[\mathrm{kpc}]$',14)
% titlemine('$\Sigma_{\mathrm{eff}}$')
%
% hf2=figure;
% contourf(md,rd,bs,0:0.001:0.05,'LineColor','none');colormap(avijet_bird)
% caxis([min(min(bs)) max(max(bs))]);
% hold on
%
% [C2,hc2]=contour(md,rd,log10(sigMat.*sigma_factor('5090')),4:0.5:10);
% set(hc2,'ShowText','on','TextStep',get(hc2,'LevelStep')*2)
%
% bar=colorbar;
% set(get(bar,'Title'),'String',bartag,'Fontsize',12,'Interpreter','latex');
% set(gca,'Fontsize',14,'box','on')
% xlabelmine('$M_s\,[\mathrm{M_\odot}]$',14)
% ylabelmine('$R_d\,[\mathrm{kpc}]$',14)
% titlemine('$\Sigma_{50-90}$')
%
% hf3=figure;
% contourf(md,rd,bs,0:0.001:0.05,'LineColor','none');colormap(avijet_bird)
% caxis([min(min(bs)) max(max(bs))]);
% hold on
% [C3,hc3]=contour(md,rd,log10(sigMat),4:0.5:10);
% set(hc3,'ShowText','on','TextStep',get(hc3,'LevelStep')*2)
%
% bar=colorbar;
% set(get(bar,'Title'),'String',bartag,'Fontsize',12,'Interpreter','latex');
% set(gca,'Fontsize',14,'box','on')
% xlabelmine('$M_s\,[\mathrm{M_\odot}]$',14)
% ylabelmine('$R_d\,[\mathrm{kpc}]$',14)
% titlemine('$\Sigma_s$')
