%% plotting future sfr estimates.
snap=99;
global simDisplayName

if readFlag
    
    global DEFAULT_MATFILE_DIR
    
    load([DEFAULT_MATFILE_DIR '/futureSFR_tff_snp99_' simDisplayName '.mat'])
    
    %load([DEFAULT_MATFILE_DIR '/fofs_subs_TNG100_z0.mat'])
    
    fofs=illustris.groupcat.loadHalos(bp,snap);
    subs=illustris.groupcat.loadSubhalos(bp,snap);
    
end
%% create central mask
subsInfo=illustris.infrastructure.build_sub_fof_connection(subs,fofs);
spbMask = mk_splashback_mask('time',5,'both',0.1);


centralMask=subsInfo.isCentral & sfrStruct.galMask & ~spbMask;

galMass=sfrStruct.galMass(centralMask);

sfr=sfrStruct.inGal.sfr(centralMask);
sfrCGM=sfrStruct.inCGM.sfr(centralMask);
sfrOut=sfrStruct.inOut.sfr(centralMask);
sfrSub=subs.SubhaloSFRinRad(centralMask);
  
ssfr=sfr./galMass ;
mask=ssfr==0;
base=10^0.5*1e-17;scat=0.5;
ssfr(mask)=base.*10.^(scat.*rand(1,sum(mask)));


sfEff=1;
ttag=10;

massG100=sfEff.*sfrStruct.inGal.fMs10(1,centralMask);
massC100=sfEff.*sfrStruct.inCGM.fMs10(1,centralMask);
massO100=sfEff.*sfrStruct.inOut.fMs10(1,centralMask);
ts100=sfrStruct.ts(1).*1e9; % in Gyr

massG500=sfEff.*sfrStruct.inGal.fMs10(2,centralMask);
massC500=sfEff.*sfrStruct.inCGM.fMs10(2,centralMask);
massO500=sfEff.*sfrStruct.inOut.fMs10(2,centralMask);
ts500=sfrStruct.ts(2).*1e9; % in Gyr

massG1=sfEff.*sfrStruct.inGal.fMs10(3,centralMask);
massC1=sfEff.*sfrStruct.inCGM.fMs10(3,centralMask);
massO1=sfEff.*sfrStruct.inOut.fMs10(3,centralMask);
ts1=sfrStruct.ts(3).*1e9; % in Gyr

massG16=sfEff.*sfrStruct.inGal.fMs10(4,centralMask);
massC16=sfEff.*sfrStruct.inCGM.fMs10(4,centralMask);
massO16=sfEff.*sfrStruct.inOut.fMs10(4,centralMask);
ts16=sfrStruct.ts(4).*1e9; % in Gyr

massG2=sfEff.*sfrStruct.inGal.fMs10(5,centralMask);
massC2=sfEff.*sfrStruct.inCGM.fMs10(5,centralMask);
massO2=sfEff.*sfrStruct.inOut.fMs10(5,centralMask);
ts2=sfrStruct.ts(5).*1e9; % in Gyr

massG5=sfEff.*sfrStruct.inGal.fMs10(6,centralMask);
massC5=sfEff.*sfrStruct.inCGM.fMs10(6,centralMask);
massO5=sfEff.*sfrStruct.inOut.fMs10(6,centralMask);
ts5=sfrStruct.ts(6).*1e9; % in Gyr

delMass100=massG100+massC100+massO100;
delMass500=massG500+massC500+massO500;
delMass1=massG1+massC1+massO1;
delMass16=massG16+massC16+massO16;
delMass2=massG2+massC2+massO2;
delMass5=massG5+massC5+massO5;

sfr100(1,:)=(massG100)./ts100;    %         ./(mass+massG100+massC100+sfrStruct.inGal.fMs100MyrSFR(centralMask))+1e-16;
sfr100(2,:)=(massG100+massC100)./ts100;
sfr100(3,:)=(massG100+massC100+massO100)./ts100;

sfr500(1,:)=(massG500)./ts500;
sfr500(2,:)=(massG500+massC500)./ts500;
sfr500(3,:)=(massG500+massC500+massO500)./ts500;

sfr1(1,:)=(massG1)./ts1;
sfr1(2,:)=(massG1+massC1)./ts1;
sfr1(3,:)=(massG1+massC1+massO1)./ts1;

sfr16(1,:)=(massG16)./ts16;
sfr16(2,:)=(massG16+massC16)./ts16;
sfr16(3,:)=(massG16+massC16+massO16)./ts16;

sfr2(1,:)=(massG2)./ts2;
sfr2(2,:)=(massG2+massC2)./ts2;
sfr2(3,:)=(massG2+massC2+massO2)./ts2;

sfr5(1,:)=(massG5)./ts5;
sfr5(2,:)=(massG5+massC5)./ts5;
sfr5(3,:)=(massG5+massC5+massO5)./ts5;




%base=10^0.5*1e-17;
%scat=0.5;

ssfr0=base.*10.^(scat.*rand(1,sum(centralMask)));

ssfr100(1,:)=sfr100(1,:)./(galMass+massG100)+ssfr0;
ssfr100(2,:)=sfr100(2,:)./(galMass+massG100+massC100)+ssfr0;
ssfr100(3,:)=sfr100(3,:)./(galMass+massG100+massC100+massO100)+ssfr0;

ssfr500(1,:)=sfr500(1,:)./(galMass+massG500)+ssfr0;
ssfr500(2,:)=sfr500(2,:)./(galMass+massG500+massC500)+ssfr0;
ssfr500(3,:)=sfr500(3,:)./(galMass+massG500+massC500+massO500)+ssfr0;

ssfr1(1,:)=sfr1(1,:)./(galMass+massG1)+ssfr0;
ssfr1(2,:)=sfr1(2,:)./(galMass+massG1+massC1)+ssfr0;
ssfr1(3,:)=sfr1(3,:)./(galMass+massG1+massC1+massO1)+ssfr0;

ssfr16(1,:)=sfr16(1,:)./(galMass+massG16)+ssfr0;
ssfr16(2,:)=sfr16(2,:)./(galMass+massG16+massC16)+ssfr0;
ssfr16(3,:)=sfr16(3,:)./(galMass+massG16+massC16+massO16)+ssfr0;

ssfr2(1,:)=sfr2(1,:)./(galMass+massG2)+ssfr0;
ssfr2(2,:)=sfr2(2,:)./(galMass+massG2+massC2)+ssfr0;
ssfr2(3,:)=sfr2(3,:)./(galMass+massG2+massC2+massO2)+ssfr0;

ssfr5(1,:)=sfr5(1,:)./(galMass+massG5)+ssfr0;
ssfr5(2,:)=sfr5(2,:)./(galMass+massG5+massC5)+ssfr0;
ssfr5(3,:)=sfr5(3,:)./(galMass+massG5+massC5+massO5)+ssfr0;


%% plotting




% figure
% loglog(sfr,sfr100(1,:),'.')
% hold on
% loglog(sfr,sfr100(2,:),'.')
% loglog(sfr,sfr100(3,:),'.')
% loglog([1e-4 1e2],[1e-4 1e2],'k--')

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
titlemine('$t_s=100\,\mathrm{Myr}$');

name=sprintf('ssfrFuture_tf%i_comp_ts%s_%s',ttag,'100M',simDisplayName);
printout_fig(gcf,name);



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
titlemine('$t_s=500\,\mathrm{Myr}$');
name=sprintf('ssfrFuture_tf%i_comp_ts%s_%s',ttag,'500M',simDisplayName);
printout_fig(gcf,name);


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
titlemine('$t_s=1\,\mathrm{Gyr}$');
name=sprintf('ssfrFuture_tf%i_comp_ts%s_%s',ttag,'1G',simDisplayName);
printout_fig(gcf,name);

figure
h=[];
h(1)=loglog(ssfr16(1,:),ssfr16(2,:),'.','displayname','G+C vs.\@ G');
hold on
h(2)=loglog(ssfr16(2,:),ssfr16(3,:),'.','displayname','G+C+O vs.\@ G+C');
loglog([1e-17 1e-8],[1e-11 1e-11],'--k')
loglog([1e-11 1e-11],[1e-17 1e-8],'--k')
loglog([1e-17 1e-8],[1e-17 1e-8],'--k')
hl=legend(h);
set(hl,'Interpreter','latex','fontsize',14,'location','SouthEast')
xlim([1e-17 5e-9])
ylim([1e-17 5e-9])
xlabelmine('sSFR');
ylabelmine('sSFR');
titlemine('$t_s=1.6\,\mathrm{Gyr}$');
name=sprintf('ssfrFuture_tf%i_comp_ts%s_%s',ttag,'16G',simDisplayName);
printout_fig(gcf,name);



figure
h=[];
h(1)=loglog(ssfr2(1,:),ssfr2(2,:),'.','displayname','G+C vs.\@ G');
hold on
h(2)=loglog(ssfr2(2,:),ssfr2(3,:),'.','displayname','G+C+O vs.\@ G+C');
loglog([1e-17 1e-8],[1e-11 1e-11],'--k')
loglog([1e-11 1e-11],[1e-17 1e-8],'--k')
loglog([1e-17 1e-8],[1e-17 1e-8],'--k')
hl=legend(h);
set(hl,'Interpreter','latex','fontsize',14,'location','SouthEast')
xlim([1e-17 5e-9])
ylim([1e-17 5e-9])
xlabelmine('sSFR');
ylabelmine('sSFR');
titlemine('$t_s=2\,\mathrm{Gyr}$')
name=sprintf('ssfrFuture_tf%i_comp_ts%s_%s',ttag,'2G',simDisplayName);
printout_fig(gcf,name);










filt=fspecial('disk',6);
xdata=log10(galMass);
ydata0=log10(ssfr);
ydata100=log10(ssfr100(3,:));
ydata500=log10(ssfr500(3,:));
ydata1=log10(ssfr1(3,:));
ydata16=log10(ssfr16(3,:));
ydata2=log10(ssfr2(3,:));
ydata5=log10(ssfr5(3,:));
%obal popCont
popCont0=plot_population_contour(xdata,ydata0,'smooth',filt,'noplot');
popCont100=plot_population_contour(xdata,ydata100,'smooth',filt,'noplot');
popCont500=plot_population_contour(xdata,ydata500,'smooth',filt,'noplot');
popCont1=plot_population_contour(xdata,ydata1,'smooth',filt,'noplot');
popCont16=plot_population_contour(xdata,ydata16,'smooth',filt,'noplot');
popCont2=plot_population_contour(xdata,ydata2,'smooth',filt,'noplot');
popCont5=plot_population_contour(xdata,ydata5,'smooth',filt,'noplot');

cc=brewermap(8,'Set1');

figure
set(gcf,'color','w')

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
name=sprintf('ssfrFuture_tf%i_ts%s_%s',ttag,'100M',simDisplayName);
printout_fig(gcf,name);




figure
set(gcf,'color','w')

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
name=sprintf('ssfrFuture_tf%i_ts%s_%s',ttag,'500M',simDisplayName);
printout_fig(gcf,name);




figure
set(gcf,'color','w')

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
name=sprintf('ssfrFuture_tf%i_ts%s_%s',ttag,'1G',simDisplayName);
printout_fig(gcf,name);


figure
set(gcf,'color','w')

plot(xdata,ydata16,'.','color',cc(1,:))
hold on

contour(popCont16.xx,popCont16.yy,popCont16.popContour,'ShowText','off','LineColor',[0 0 0],...
    'LevelList',20:20:100,'Fill','off','linestyle','-');
grid
ylim([-16.5 -8.5])

plot([9 13],[-11 -11],':k')
xlabelmine('stellar mass (old)');
ylabelmine('sSFR projected');
titlemine('$t_s=1.6\,\mathrm{Gyr}$');
name=sprintf('ssfrFuture_tf%i_ts%s_%s',ttag,'16G',simDisplayName);
printout_fig(gcf,name);


figure
set(gcf,'color','w')

plot(xdata,ydata2,'.','color',cc(1,:))
hold on

contour(popCont2.xx,popCont2.yy,popCont2.popContour,'ShowText','off','LineColor',[0 0 0],...
    'LevelList',20:20:100,'Fill','off','linestyle','-');
grid
ylim([-16.5 -8.5])

plot([9 13],[-11 -11],':k')
xlabelmine('stellar mass (old)');
ylabelmine('sSFR projected');
titlemine('$t_s=2\,\mathrm{Gyr}$');
name=sprintf('ssfrFuture_tf%i_ts%s_%s',ttag,'2G',simDisplayName);
printout_fig(gcf,name);


figure
set(gcf,'color','w')

plot(xdata,ydata0,'.','color',cc(2,:))
hold on
plot(xdata,ydata2,'.','color',cc(1,:))

contour(popCont2.xx,popCont2.yy,popCont2.popContour,'ShowText','off','LineColor',[0 0 0],...
    'LevelList',20:20:100,'Fill','off','linestyle','-');
grid
ylim([-16.5 -8.5])

plot([9 13],[-11 -11],':k')
xlabelmine('stellar mass (old)');
ylabelmine('sSFR projected');
titlemine('$t_s=2\,\mathrm{Gyr}$');
name=sprintf('ssfrFuture_currentSSFR_tf%i_ts%s_%s',ttag,'2G',simDisplayName);
printout_fig(gcf,name);



figure
set(gcf,'color','w')
plot(xdata,ydata5,'.','color',cc(1,:))
hold on

contour(popCont5.xx,popCont5.yy,popCont5.popContour,'ShowText','off','LineColor',[0 0 0],...
    'LevelList',20:20:100,'Fill','off','linestyle','-');
grid
ylim([-16.5 -8.5])

plot([9 13],[-11 -11],':k')
xlabelmine('stellar mass (old)');
ylabelmine('sSFR projected');
titlemine('$t_s=5\,\mathrm{Gyr}$');
name=sprintf('ssfrFuture_tf%i_ts%s_%s',ttag,'5G',simDisplayName);
printout_fig(gcf,name);
%% 

xdata=log10(galMass);
ydelM100=log10(delMass100(3,:));
ydelM500=log10(delMass500(3,:));
ydelM1=log10(delMass1(3,:));
ydelM16=log10(delMass16(3,:));
ydelM2=log10(delMass2(3,:));
ydelM5=log10(delMass5(3,:));
%obal popCont
popContDM100=plot_population_contour(xdata,ydelM100,'smooth',filt,'noplot');
popContDM500=plot_population_contour(xdata,ydelM500,'smooth',filt,'noplot');
popContDM1=plot_population_contour(xdata,ydelM1,'smooth',filt,'noplot');
popContDM16=plot_population_contour(xdata,ydelM16,'smooth',filt,'noplot');
popContDM2=plot_population_contour(xdata,ydelM2,'smooth',filt,'noplot');
popContDM5=plot_population_contour(xdata,ydelM5,'smooth',filt,'noplot');


figure
set(gcf,'color','w')
mask=delMass100==0;
add=zeros(size(delMass100));

base=1e4;scat=0.5;
add(mask)=base.*10.^(scat.*rand(1,sum(mask)));
loglog(galMass,(delMass100+add),'.')
xlabelmine('stellar mass (old)');
ylabelmine('$\Delta$Mass');
titlemine('$t_s=100\,\mathrm{Myr}$');
name=sprintf('smassFuture_tf%i_ts%s_%s',ttag,'100M',simDisplayName);
printout_fig(gcf,name);

figure
set(gcf,'color','w')
mask=delMass500==0;
add=zeros(size(delMass500));

base=1e4;scat=0.5;
add(mask)=base.*10.^(scat.*rand(1,sum(mask)));
loglog(galMass,(delMass500+add),'.')
xlabelmine('stellar mass (old)');
ylabelmine('$\Delta$Mass');
titlemine('$t_s=500\,\mathrm{Myr}$');
name=sprintf('smassFuture_tf%i_ts%s_%s',ttag,'500M',simDisplayName);
printout_fig(gcf,name,'subdir','futureSsfr');

figure
set(gcf,'color','w')
mask=delMass1==0;
add=zeros(size(delMass1));

base=1e4;scat=0.5;
add(mask)=base.*10.^(scat.*rand(1,sum(mask)));
loglog(galMass,(delMass1+add),'.')
xlabelmine('stellar mass (old)');
ylabelmine('$\Delta$Mass');
titlemine('$t_s=1\,\mathrm{Gyr}$');
name=sprintf('smassFuture_tf%i_ts%s_%s',ttag,'1G',simDisplayName);
printout_fig(gcf,name);

figure
set(gcf,'color','w')
mask=delMass16==0;
add=zeros(size(delMass16));

base=1e4;scat=0.5;
add(mask)=base.*10.^(scat.*rand(1,sum(mask)));
loglog(galMass,(delMass16+add),'.')
xlabelmine('stellar mass (old)');
ylabelmine('$\Delta$Mass');
titlemine('$t_s=1.6\,\mathrm{Gyr}$');
name=sprintf('smassFuture_tf%i_ts%s_%s',ttag,'16G',simDisplayName);
printout_fig(gcf,name);

figure
set(gcf,'color','w')
mask=delMass2==0;
add=zeros(size(delMass2));

base=1e4;scat=0.5;
add(mask)=base.*10.^(scat.*rand(1,sum(mask)));
loglog(galMass,(delMass2+add),'.')
xlabelmine('stellar mass (old)');
ylabelmine('$\Delta$Mass');
titlemine('$t_s=2\,\mathrm{Gyr}$');
name=sprintf('smassFuture_tf%i_ts%s_%s',ttag,'2G',simDisplayName);
printout_fig(gcf,name);


figure
set(gcf,'color','w')
mask=delMass5==0;
add=zeros(size(delMass5));

base=1e4;scat=0.5;
add(mask)=base.*10.^(scat.*rand(1,sum(mask)));
loglog(galMass,(delMass5+add),'.')
xlabelmine('stellar mass (old)');
ylabelmine('$\Delta$Mass');
titlemine('$t_s=5\,\mathrm{Gyr}$');
name=sprintf('smassFuture_tf%i_ts%s_%s',ttag,'5G',simDisplayName);
printout_fig(gcf,name);


%
%
% contour(popCont200.xx,popCont200.yy,popCont200.popContour,'ShowText','off','LineColor',cc(2,:),...
%     'LevelList',25:25:100,'Fill','off','linestyle','-');
% contour(popCont500.xx,popCont500.yy,popCont500.popContour,'ShowText','off','LineColor',cc(3,:),...
%     'LevelList',25:25:100,'Fill','off','linestyle','-');
% contour(popCont1.xx,popCont1.yy,popCont1.popContour,'ShowText','off','LineColor',cc(4,:),...
%     'LevelList',25:25:100,'Fill','off','linestyle','-');
%
% figure
% loglog(galMass,ssfr100(3,:),'.')
% hold on
% loglog(galMass,ssfr200(3,:),'.')
% loglog(galMass,ssfr500(3,:),'.')
% loglog(galMass,ssfr1(3,:),'.')
% % loglog(mass,ssfr100(2,:),'r.')
% % loglog(mass,ssfr100(3,:),'g.')
% loglog([1e9 1e13],[1e-11 1e-11],':k')
% ylim([1e-17 1e-8])
%
% %sfr100(1)=(massG100+massC100)/1e8)./(mass+massG100+massC100+sfrStruct.inGal.fMs100MyrSFR(centralMask))+1e-16;