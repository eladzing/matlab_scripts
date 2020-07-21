%% plotting future sfr estimates. 


global DEFAULT_MATFILE_DIR

load([DEFAULT_MATFILE_DIR '/futureSFR_tff_snp99_TNG100.mat'])
load([DEFAULT_MATFILE_DIR '/fofs_subs_TNG100_z0.mat'])

%% create central mask 
subsInfo=illustris.infrastructure.build_sub_fof_connection(subs,fofs);
centralMask=subsInfo.isCentral & sfrStruct.galMask;

galMass=sfrStruct.galMass(centralMask);

sfr=sfrStruct.inGal.sfr(centralMask);
sfr2=subs.SubhaloSFRinRad(centralMask);



sfEff=0.01;

massG100=sfEff.*sfrStruct.inGal.fMs100MyrCool(centralMask);
massC100=sfEff.*sfrStruct.inCGM.fMs100MyrCool(centralMask);
massO100=sfEff.*sfrStruct.inOut.fMs100MyrCool(centralMask);
ts100=100.*1e6;

massG200=sfEff.*sfrStruct.inGal.fMs200MyrCool(centralMask);
massC200=sfEff.*sfrStruct.inCGM.fMs200MyrCool(centralMask);
massO200=sfEff.*sfrStruct.inOut.fMs200MyrCool(centralMask);
ts200=200.*1e6;

massG500=sfEff.*sfrStruct.inGal.fMs500MyrCool(centralMask);
massC500=sfEff.*sfrStruct.inCGM.fMs500MyrCool(centralMask);
massO500=sfEff.*sfrStruct.inOut.fMs500MyrCool(centralMask);
ts500=500.*1e6;

massG1=sfEff.*sfrStruct.inGal.fMs1GyrCool(centralMask);
massC1=sfEff.*sfrStruct.inCGM.fMs1GyrCool(centralMask);
massO1=sfEff.*sfrStruct.inOut.fMs1GyrCool(centralMask);
ts1G=1e9;

sfr100(1,:)=(massG100)./ts100;    %         ./(mass+massG100+massC100+sfrStruct.inGal.fMs100MyrSFR(centralMask))+1e-16;
sfr100(2,:)=(massG100+massC100)./ts100;
sfr100(3,:)=(massG100+massC100+massO100)./ts100;

sfr200(1,:)=(massG200)./ts200;    
sfr200(2,:)=(massG200+massC200)./ts200;
sfr200(3,:)=(massG200+massC200+massO200)./ts200;

sfr500(1,:)=(massG500)./ts500;    
sfr500(2,:)=(massG500+massC500)./ts500;
sfr500(3,:)=(massG500+massC500+massO500)./ts500;

sfr1(1,:)=(massG1)./ts1G;    
sfr1(2,:)=(massG1+massC1)./ts1G;
sfr1(3,:)=(massG1+massC1+massO1)./ts1G;

base=10^0.5*1e-17;
scat=0.5;

ssfr0=base.*10.^(scat.*rand(1,sum(centralMask)));

ssfr100(1,:)=sfr100(1,:)./(galMass+massG100)+ssfr0;
ssfr100(2,:)=sfr100(2,:)./(galMass+massG100+massC100)+ssfr0;
ssfr100(3,:)=sfr100(3,:)./(galMass+massG100+massC100+massO100)+ssfr0;

ssfr200(1,:)=sfr200(1,:)./(galMass+massG200)+ssfr0;
ssfr200(2,:)=sfr200(2,:)./(galMass+massG200+massC200)+ssfr0;
ssfr200(3,:)=sfr200(3,:)./(galMass+massG200+massC200+massO200)+ssfr0;

ssfr500(1,:)=sfr500(1,:)./(galMass+massG500)+ssfr0;
ssfr500(2,:)=sfr500(2,:)./(galMass+massG500+massC500)+ssfr0;
ssfr500(3,:)=sfr500(3,:)./(galMass+massG500+massC500+massO500)+ssfr0;

ssfr1(1,:)=sfr1(1,:)./(galMass+massG1)+ssfr0;
ssfr1(2,:)=sfr1(2,:)./(galMass+massG1+massC1)+ssfr0;
ssfr1(3,:)=sfr1(3,:)./(galMass+massG1+massC1+massO1)+ssfr0;




%% plotting 




figure 
loglog(sfr,sfr100(1,:),'.')
hold on
loglog(sfr,sfr100(2,:),'.')
loglog(sfr,sfr100(3,:),'.')
loglog([1e-4 1e2],[1e-4 1e2],'k--')

% figure
% h=[];
% h(1)=loglog(ssfr100(1,:),ssfr100(2,:)./ssfr100(1,:),'.','displayname','G+C/G');
% hold on
% h(2)=loglog(ssfr100(2,:),ssfr100(3,:)./ssfr100(2,:),'.','displayname','G+C+O/G+C');
% 
% hl=legend(h);
% set(hl,'Interpreter','latex','fontsize',14)
% 
% xlabelmine('sSFR')
% ylabelmine('sSFR Ratio')
% grid 


figure
h=[];
h(1)=loglog(ssfr100(1,:),ssfr100(2,:),'.','displayname','G+C vs.\@ G');
hold on
h(2)=loglog(ssfr100(2,:),ssfr100(3,:),'.','displayname','G+C+O vs.\@ G+C');
loglog([1e-17 1e-8],[1e-11 1e-11],'--k')
loglog([1e-11 1e-11],[1e-17 1e-8],'--k')
loglog([1e-17 1e-8],[1e-17 1e-8],'--k')
hl=legend(h);
set(hl,'Interpreter','latex','fontsize',14,'location','SouthEast')
xlim([1e-17 5e-9])
ylim([1e-17 5e-9])
xlabelmine('sSFR');
ylabelmine('sSFR');
titlemine('$t_s=100\,\mathrm{Myr}$')

figure
h=[];
h(1)=loglog(ssfr200(1,:),ssfr200(2,:),'.','displayname','G+C vs.\@ G');
hold on
h(2)=loglog(ssfr200(2,:),ssfr200(3,:),'.','displayname','G+C+O vs.\@ G+C');
loglog([1e-17 1e-8],[1e-11 1e-11],'--k')
loglog([1e-11 1e-11],[1e-17 1e-8],'--k')
loglog([1e-17 1e-8],[1e-17 1e-8],'--k')
hl=legend(h);
set(hl,'Interpreter','latex','fontsize',14,'location','SouthEast')
xlim([1e-17 5e-9])
ylim([1e-17 5e-9])
xlabelmine('sSFR');
ylabelmine('sSFR');
titlemine('$t_s=200\,\mathrm{Myr}$')

figure
h=[];
h(1)=loglog(ssfr500(1,:),ssfr500(2,:),'.','displayname','G+C vs.\@ G');
hold on
h(2)=loglog(ssfr500(2,:),ssfr500(3,:),'.','displayname','G+C+O vs.\@ G+C');
loglog([1e-17 1e-8],[1e-11 1e-11],'--k')
loglog([1e-11 1e-11],[1e-17 1e-8],'--k')
loglog([1e-17 1e-8],[1e-17 1e-8],'--k')
hl=legend(h);
set(hl,'Interpreter','latex','fontsize',14,'location','SouthEast')
xlim([1e-17 5e-9])
ylim([1e-17 5e-9])
xlabelmine('sSFR');
ylabelmine('sSFR');
titlemine('$t_s=500\,\mathrm{Myr}$')

figure
h=[];
h(1)=loglog(ssfr1(1,:),ssfr1(2,:),'.','displayname','G+C vs.\@ G');
hold on
h(2)=loglog(ssfr1(2,:),ssfr1(3,:),'.','displayname','G+C+O vs.\@ G+C');
loglog([1e-17 1e-8],[1e-11 1e-11],'--k')
loglog([1e-11 1e-11],[1e-17 1e-8],'--k')
loglog([1e-17 1e-8],[1e-17 1e-8],'--k')
hl=legend(h);
set(hl,'Interpreter','latex','fontsize',14,'location','SouthEast')
xlim([1e-17 5e-9])
ylim([1e-17 5e-9])
xlabelmine('sSFR');
ylabelmine('sSFR');
titlemine('$t_s=1\,\mathrm{Gyr}$')




filt=fspecial('disk',6);
xdata=log10(galMass);
ydata100=log10(ssfr100(2,:));
ydata200=log10(ssfr200(2,:));
ydata500=log10(ssfr500(2,:));
ydata1=log10(ssfr1(2,:));
    %obal popCont
popCont100=plot_population_contour(xdata,ydata100,'smooth',filt,'noplot');
popCont200=plot_population_contour(xdata,ydata200,'smooth',filt,'noplot');
popCont500=plot_population_contour(xdata,ydata500,'smooth',filt,'noplot');
popCont1=plot_population_contour(xdata,ydata1,'smooth',filt,'noplot');

 cc=brewermap(8,'Set1');

figure
plot(xdata,ydata100,'.','color',cc(1,:))
hold on 

contour(popCont100.xx,popCont100.yy,popCont100.popContour,'ShowText','off','LineColor',[0 0 0],...
    'LevelList',20:20:100,'Fill','off','linestyle','-');
grid
ylim([-16.5 -8.5])

plot([9 13],[-11 -11],':k')
xlabelmine('stellar mass (old)');
ylabelmine('sSFR projected');
titlemine('$t_s=100\,\mathrm{Myr}$');

figure
plot(xdata,ydata200,'.','color',cc(1,:))
hold on 

contour(popCont200.xx,popCont200.yy,popCont200.popContour,'ShowText','off','LineColor',[0 0 0],...
    'LevelList',20:20:100,'Fill','off','linestyle','-');
grid
ylim([-16.5 -8.5])

plot([9 13],[-11 -11],':k')
xlabelmine('stellar mass (old)');
ylabelmine('sSFR projected');
titlemine('$t_s=200\,\mathrm{Myr}$');

figure
plot(xdata,ydata500,'.','color',cc(1,:))
hold on 

contour(popCont500.xx,popCont500.yy,popCont500.popContour,'ShowText','off','LineColor',[0 0 0],...
    'LevelList',20:20:100,'Fill','off','linestyle','-');
grid
ylim([-16.5 -8.5])

plot([9 13],[-11 -11],':k')
xlabelmine('stellar mass (old)');
ylabelmine('sSFR projected');
titlemine('$t_s=500\,\mathrm{Myr}$');

figure
plot(xdata,ydata1,'.','color',cc(1,:))
hold on 

contour(popCont1.xx,popCont1.yy,popCont1.popContour,'ShowText','off','LineColor',[0 0 0],...
    'LevelList',20:20:100,'Fill','off','linestyle','-');
grid
ylim([-16.5 -8.5])

plot([9 13],[-11 -11],':k')
xlabelmine('stellar mass (old)');
ylabelmine('sSFR projected');
titlemine('$t_s=1\,\mathrm{Gyr}$');




contour(popCont200.xx,popCont200.yy,popCont200.popContour,'ShowText','off','LineColor',cc(2,:),...
    'LevelList',25:25:100,'Fill','off','linestyle','-');
contour(popCont500.xx,popCont500.yy,popCont500.popContour,'ShowText','off','LineColor',cc(3,:),...
    'LevelList',25:25:100,'Fill','off','linestyle','-');
contour(popCont1.xx,popCont1.yy,popCont1.popContour,'ShowText','off','LineColor',cc(4,:),...
    'LevelList',25:25:100,'Fill','off','linestyle','-');

figure 
loglog(galMass,ssfr100(3,:),'.')
hold on
loglog(galMass,ssfr200(3,:),'.')
loglog(galMass,ssfr500(3,:),'.')
loglog(galMass,ssfr1(3,:),'.')
% loglog(mass,ssfr100(2,:),'r.')
% loglog(mass,ssfr100(3,:),'g.')
loglog([1e9 1e13],[1e-11 1e-11],':k')
ylim([1e-17 1e-8])

%sfr100(1)=(massG100+massC100)/1e8)./(mass+massG100+massC100+sfrStruct.inGal.fMs100MyrSFR(centralMask))+1e-16;