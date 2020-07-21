% global DEFAULT_MATFILE_DIR
% name='orbParam_mean_vs_crit.mat';
% load([DEFAULT_MATFILE_DIR '/' name],'orbParam')
% fprintf('saving to: %s \n',[DEFAULT_MATFILE_DIR '/' name])


%% plot individual hosts/mass ratio 
cc=brewermap(8,'Set1');
figure('color','w','position',[-1143  384 1144 945])


subplot(3,3,7)

vv1=abs(orbParam.mh1e12.mr05.vr);
vv2=abs(orbParam.mh5e12.mr05.vr);

vv3=cat(2,vv1,vv2);


vp=radialVelocity_pdf(2e12,0.25);

len=50;

[vvh1, bs1, xxlimH1]= histogram1d(vv1,ones(size(vv1)),'bins',len);
xx0=linspace(xxlimH1(1),xxlimH1(2),len+1);
xx1=xx0(1:end-1)+0.5.*diff(xx0);


[vvh2, bs2, xxlimH2]= histogram1d(vv2,ones(size(vv2)),'bins',len);
xx0=linspace(xxlimH2(1),xxlimH2(2),len+1);
xx2=xx0(1:end-1)+0.5.*diff(xx0);

[vvh3, bs3, xxlimH3]= histogram1d(vv3,ones(size(vv3)),'bins',len);
xx0=linspace(xxlimH3(1),xxlimH3(2),len+1);
xx3=xx0(1:end-1)+0.5.*diff(xx0);

plot(xx1,vvh1(:,1)./length(vv1)./bs1,'--','color',cc(1,:))
hold on 
plot(xx2,vvh2(:,1)./length(vv2)./bs2,'--','color',cc(2,:))
plot(xx3,vvh3(:,1)./length(vv3)./bs3,'-k')
plot(vp.vv,vp.pdf./trapz(vp.vv,vp.pdf),'m')

grid;xlim([0 1]);ylim([0 3.25]);
xlabelmine('$V_\mathrm{r}/V$');
ylabelmine('$P(V_\mathrm{r}/V)$');
set(gca,'fontsize',14)
titlemine('$10^{12},\, <0.5$');

%%

subplot(3,3,4)

vv1=abs(orbParam.mh5e12.mr005.vr);
vp=radialVelocity_pdf(5e12,0.025);

len=50;

[vvh1, bs1, xxlimH1]= histogram1d(vv1,ones(size(vv1)),'bins',len);
xx0=linspace(xxlimH1(1),xxlimH1(2),len+1);
xx1=xx0(1:end-1)+0.5.*diff(xx0);

plot(xx1,vvh1(:,1)./length(vv1)./bs1,'-k')
hold on 

plot(vp.vv,vp.pdf./trapz(vp.vv,vp.pdf),'m')

grid;xlim([0 1]);ylim([0 3.25]);
%xlabelmine('$V_\mathrm{r}/V$');
ylabelmine('$P(V_\mathrm{r}/V)$');
set(gca,'fontsize',14)
titlemine('$10^{12},\, <0.05$');

%% 

subplot(3,3,8)

vv1=abs(orbParam.mh1e13.mr05.vr);
vv2=abs(orbParam.mh5e13.mr05.vr);

vv3=cat(2,vv1,vv2);

vp=radialVelocity_pdf(2e13,0.25);

len=50;

[vvh1, bs1, xxlimH1]= histogram1d(vv1,ones(size(vv1)),'bins',len);
xx0=linspace(xxlimH1(1),xxlimH1(2),len+1);
xx1=xx0(1:end-1)+0.5.*diff(xx0);


[vvh2, bs2, xxlimH2]= histogram1d(vv2,ones(size(vv2)),'bins',len);
xx0=linspace(xxlimH2(1),xxlimH2(2),len+1);
xx2=xx0(1:end-1)+0.5.*diff(xx0);

[vvh3, bs3, xxlimH3]= histogram1d(vv3,ones(size(vv3)),'bins',len);
xx0=linspace(xxlimH3(1),xxlimH3(2),len+1);
xx3=xx0(1:end-1)+0.5.*diff(xx0);

plot(xx1,vvh1(:,1)./length(vv1)./bs1,'--','color',cc(1,:))
hold on 
plot(xx2,vvh2(:,1)./length(vv2)./bs2,'--','color',cc(2,:))
plot(xx3,vvh3(:,1)./length(vv3)./bs3,'-k')
plot(vp.vv,vp.pdf./trapz(vp.vv,vp.pdf),'m')

grid;xlim([0 1]);ylim([0 3.25]);
xlabelmine('$V_\mathrm{r}/V$');
%ylabelmine('$P(V_\mathrm{r}/V)$');
set(gca,'fontsize',14)
titlemine('$10^{13},\, <0.5$');

%% 

subplot(3,3,2)

vv1=abs(orbParam.mh5e13.mr0005.vr);

vp=radialVelocity_pdf(5e13,0.0025);
len=50;

[vvh1, bs1, xxlimH1]= histogram1d(vv1,ones(size(vv1)),'bins',len);
xx0=linspace(xxlimH1(1),xxlimH1(2),len+1);
xx1=xx0(1:end-1)+0.5.*diff(xx0);
h=[];
h(1)=plot(xx1,vvh1(:,1)./length(vv1)./bs1,'-k',...
    'DisplayName','$V(R_\mathrm{200,c})$');
hold on 

h(2)=plot(vp.vv,vp.pdf./trapz(vp.vv,vp.pdf),'m',...
    'DisplayName','PDF');

hl=legend(h);
set(hl,'Interpreter','latex','Location','NorthEast','fontsize',10)

grid;xlim([0 1]);ylim([0 3.25]);
%xlabelmine('$V_\mathrm{r}/V$');
ylabelmine('$P(V_\mathrm{r}/V)$');
set(gca,'fontsize',14)
titlemine('$10^{13},\, <0.005$');

%% 

subplot(3,3,5)

%vv1=abs(orbParam.mh1e13.mr005.vv.*orbParam.mh1e13.v200c./orbParam.mh1e13.host.Vvir;
%vv2=abs(orbParam.mh5e13.mr005.vv.*orbParam.mh5e13.v200c./orbParam.mh5e13.host.Vvir;
vv1=abs(orbParam.mh1e13.mr005.vr);
vv2=abs(orbParam.mh5e13.mr005.vr);
vv3=cat(2,vv1,vv2);

vp=radialVelocity_pdf(2e13,0.025);

len=50;

[vvh1, bs1, xxlimH1]= histogram1d(vv1,ones(size(vv1)),'bins',len);
xx0=linspace(xxlimH1(1),xxlimH1(2),len+1);
xx1=xx0(1:end-1)+0.5.*diff(xx0);


[vvh2, bs2, xxlimH2]= histogram1d(vv2,ones(size(vv2)),'bins',len);
xx0=linspace(xxlimH2(1),xxlimH2(2),len+1);
xx2=xx0(1:end-1)+0.5.*diff(xx0);

[vvh3, bs3, xxlimH3]= histogram1d(vv3,ones(size(vv3)),'bins',len);
xx0=linspace(xxlimH3(1),xxlimH3(2),len+1);
xx3=xx0(1:end-1)+0.5.*diff(xx0);

plot(xx1,vvh1(:,1)./length(vv1)./bs1,'--','color',cc(1,:))
hold on 
plot(xx2,vvh2(:,1)./length(vv2)./bs2,'--','color',cc(2,:))
plot(xx3,vvh3(:,1)./length(vv3)./bs3,'-k')
plot(vp.vv,vp.pdf./trapz(vp.vv,vp.pdf),'m')

grid;xlim([0 1]);ylim([0 3.25]);
%xlabelmine('$V_\mathrm{r}/V$');
%ylabelmine('$P(V_\mathrm{r}/V)$');
set(gca,'fontsize',14)
titlemine('$10^{13},\, <0.05$');


%% plot individual hosts/mass ratio 

subplot(3,3,3)

vv1=abs(orbParam.mh1e14.mr0005.vr);
vv2=abs(orbParam.mh5e14.mr0005.vr);
vv3=abs(orbParam.mh1e15.mr0005.vr);

vv4=cat(2,vv1,vv2,vv3);
vp=radialVelocity_pdf(2e14,0.0025);

len=50;

[vvh1, bs1, xxlimH1]= histogram1d(vv1,ones(size(vv1)),'bins',len);
xx0=linspace(xxlimH1(1),xxlimH1(2),len+1);
xx1=xx0(1:end-1)+0.5.*diff(xx0);

[vvh2, bs2, xxlimH2]= histogram1d(vv2,ones(size(vv2)),'bins',len);
xx0=linspace(xxlimH2(1),xxlimH2(2),len+1);
xx2=xx0(1:end-1)+0.5.*diff(xx0);

[vvh3, bs3, xxlimH3]= histogram1d(vv3,ones(size(vv3)),'bins',len);
xx0=linspace(xxlimH3(1),xxlimH3(2),len+1);
xx3=xx0(1:end-1)+0.5.*diff(xx0);

[vvh4, bs4, xxlimH4]= histogram1d(vv4,ones(size(vv4)),'bins',len);
xx0=linspace(xxlimH4(1),xxlimH4(2),len+1);
xx4=xx0(1:end-1)+0.5.*diff(xx0);

plot(xx1,vvh1(:,1)./length(vv1)./bs1,'--','color',cc(1,:))
hold on 
plot(xx2,vvh2(:,1)./length(vv2)./bs2,'--','color',cc(2,:))
plot(xx3,vvh3(:,1)./length(vv3)./bs3,'--','color',cc(3,:))
plot(xx4,vvh4(:,1)./length(vv4)./bs4,'-k')
plot(vp.vv,vp.pdf./trapz(vp.vv,vp.pdf),'m')


grid;xlim([0 1]);ylim([0 3.25]);
%xlabelmine('$V_\mathrm{r}/V$');
%ylabelmine('$P(V_\mathrm{r}/V)$');
set(gca,'fontsize',14)
titlemine('$10^{14},\, <0.005$');

%%


subplot(3,3,6)

vv1=abs(orbParam.mh1e14.mr005.vr);
vv2=abs(orbParam.mh5e14.mr005.vr);
vv3=abs(orbParam.mh1e15.mr005.vr);

vv4=cat(2,vv1,vv2,vv3);

vp=radialVelocity_pdf(2e14,0.025);

len=50;

[vvh1, bs1, xxlimH1]= histogram1d(vv1,ones(size(vv1)),'bins',len);
xx0=linspace(xxlimH1(1),xxlimH1(2),len+1);
xx1=xx0(1:end-1)+0.5.*diff(xx0);

[vvh2, bs2, xxlimH2]= histogram1d(vv2,ones(size(vv2)),'bins',len);
xx0=linspace(xxlimH2(1),xxlimH2(2),len+1);
xx2=xx0(1:end-1)+0.5.*diff(xx0);

[vvh3, bs3, xxlimH3]= histogram1d(vv3,ones(size(vv3)),'bins',len);
xx0=linspace(xxlimH3(1),xxlimH3(2),len+1);
xx3=xx0(1:end-1)+0.5.*diff(xx0);

[vvh4, bs4, xxlimH4]= histogram1d(vv4,ones(size(vv4)),'bins',len);
xx0=linspace(xxlimH4(1),xxlimH4(2),len+1);
xx4=xx0(1:end-1)+0.5.*diff(xx0);

plot(xx1,vvh1(:,1)./length(vv1)./bs1,'--','color',cc(1,:))
hold on 
plot(xx2,vvh2(:,1)./length(vv2)./bs2,'--','color',cc(2,:))
plot(xx3,vvh3(:,1)./length(vv3)./bs3,'--','color',cc(3,:))
plot(xx4,vvh4(:,1)./length(vv4)./bs4,'-k')
plot(vp.vv,vp.pdf./trapz(vp.vv,vp.pdf),'m')

grid;xlim([0 1]);ylim([0 3.25]);
%xlabelmine('$V_\mathrm{r}/V$');
%ylabelmine('$P(V_\mathrm{r}/V)$');
set(gca,'fontsize',14)
titlemine('$10^{14},\, <0.05$');

%%


subplot(3,3,9)

vv1=abs(orbParam.mh1e14.mr05.vr);
vv2=abs(orbParam.mh5e14.mr05.vr);
vv3=abs(orbParam.mh1e15.mr05.vr);

vv4=cat(2,vv1,vv2,vv3);

vp=radialVelocity_pdf(2e14,0.25);

len=50;

[vvh1, bs1, xxlimH1]= histogram1d(vv1,ones(size(vv1)),'bins',len);
xx0=linspace(xxlimH1(1),xxlimH1(2),len+1);
xx1=xx0(1:end-1)+0.5.*diff(xx0);

[vvh2, bs2, xxlimH2]= histogram1d(vv2,ones(size(vv2)),'bins',len);
xx0=linspace(xxlimH2(1),xxlimH2(2),len+1);
xx2=xx0(1:end-1)+0.5.*diff(xx0);

[vvh3, bs3, xxlimH3]= histogram1d(vv3,ones(size(vv3)),'bins',len);
xx0=linspace(xxlimH3(1),xxlimH3(2),len+1);
xx3=xx0(1:end-1)+0.5.*diff(xx0);

[vvh4, bs4, xxlimH4]= histogram1d(vv4,ones(size(vv4)),'bins',len);
xx0=linspace(xxlimH4(1),xxlimH4(2),len+1);
xx4=xx0(1:end-1)+0.5.*diff(xx0);

plot(xx1,vvh1(:,1)./length(vv1)./bs1,'--','color',cc(1,:))
hold on 
plot(xx2,vvh2(:,1)./length(vv2)./bs2,'--','color',cc(2,:))
plot(xx3,vvh3(:,1)./length(vv3)./bs3,'--','color',cc(3,:))
plot(xx4,vvh4(:,1)./length(vv4)./bs4,'-k')
plot(vp.vv,vp.pdf./trapz(vp.vv,vp.pdf),'m')

grid;xlim([0 1]);ylim([0 3.25]);
xlabelmine('$V_\mathrm{r}/V$');
%ylabelmine('$P(V_\mathrm{r}/V)$');
set(gca,'fontsize',14)
titlemine('$10^{14},\, <0.5$');


        
            