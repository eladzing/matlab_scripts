%% make mock catalog
units

global cosmoStruct

%nSample=1;
ngal=100;

%% set Galaxy parameter values

% we start with small mass range: 
msRange=[0.5 1].*1e10;
ms=msRange(1)+diff(msRange).*rand(ngal,1);

fgsRange=[0.05 1];
fgs=fgsRange(1)+diff(fgsRange).*rand(ngal,1);

betaRange=[0.4 1];
beta=betaRange(1)+diff(betaRange).*rand(ngal,1);

fbsRange=[0 1];
fbs=fbsRange(1)+diff(fbsRange).*rand(ngal,1);

xiRange=[0.5 4];
xi=xiRange(1)+diff(xiRange).*rand(ngal,1);

Mv=halomass_from_moster(ms.*(1+fbs),'noshow');
cv=cvir_Mvir(Mv,0,'random');
lambda=lambda_prime(ms);

% calculate Rd 
res=rscale_mmw_array(ms,'fg',fgs,'beta',beta,'fb',fbs,'xi',xi,...
    'Mv',Mv,'cv',cv,'lambda',lambda,'noshow');

% for i=1:ngal
% gal(i)=GALAXY('ms',ms(i),'rd',res.rd(i),'fgs',fgs(i),'beta',beta(i),'fbs',fbs(i),'xi',xi(i),'Mh',Mv(i),'cv',cv(i));
% end

%% set Host parameters 
Mh=1e14;
cvh=cvir_Mvir(Mh,0,'random');
fgh=0.15;
host=NFW('mv',Mh,'cc',cvh,'fg',fgh);


%R0=1.5.*host.Rvir;
%rp=logspace(-2,0,1000).*R0;
%rhoICM=gasDensity(host,rp,'kpc');
%rp0=rp;        


%% set orbit parameters 
y0=0;
x0=1.5.*host.Rvir;
vTotat=host.Vvir; %.*rand(ngal,1);
%vAngle=pi/2.*rand(ngal,1); % in radians
vx0= -1.*vTotat;%.*cos(vAngle);
vy0=0; % vTotat.*sin(vAngle);

tOrbit=x0./host.Vvir;
tmax=1.3*tOrbit;
dt=tOrbit/1e4;
dtmin=dt/1e3;

%orbits
    orbIC.x=x0;
    orbIC.y=y0;
    orbIC.vx=vx0;%(i);
    orbIC.vy=vy0;%(i);
    
    orb=orbits.rk4_orbitIntegration(orbIC,dt,tmax,dtmin,@orbits.rhs_nfw,host);




prc0=10;
step=10;
fprintf('Running evolution \n')

for i=1:ngal
    %gals
    
    prc=i/ngal*100;
    if prc>=prc0
        fprintf('completed %s %% of catalog \n',prc0)
        prc0=prc0+step;
        toc
    end
    gal(i)=GALAXY('ms',ms(i),'rd',res.rd(i),'fgs',fgs(i),'beta',beta(i),'fbs',fbs(i),'xi',xi(i),'Mh',Mv(i),'cv',cv(i));
    
    galRes(i)=galEvolutionMachineLight(gal(i),host,orb,'nobar');
    
end























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
