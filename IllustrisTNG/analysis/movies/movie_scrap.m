%% plotting many things 

zr=galhist.zred;

len=zr;

colorSet=brewermap(9,'Set1');
grey=colorSet(9,:);
greyType=':';

ii=34;

%% ssft & gas mass 
figure 
h=[];
yyaxis left
semilogy(zr,galhist.ssfr,greyType,'color',grey)
hold on
semilogy(zr(ii:end),galhist.ssfr(ii:end),'-','color',colorSet(2,:),'linewidth',1.5)
semilogy(zr(ii),galhist.ssfr(ii),'d','color',colorSet(2,:),...
    'linewidth',1.5,'markersize',10,'MarkerFaceColor',colorSet(2,:));

ylabelmine('ssfr $[\mathrm{yr^{-1}}]$')
ylim([1e-17 1e-8])
yyaxis right

semilogy(zr,galhist.stellarMass./1e10,greyType,'color',grey) 
hold on
h(1)=semilogy(zr(ii:end),galhist.stellarMass(ii:end)./1e10,'-',...
    'color',colorSet(1,:),'linewidth',1.5,'DisplayName','Stellar');
semilogy(zr(ii),galhist.stellarMass(ii)./1e10,'d','color',colorSet(1,:),...
    'linewidth',1.5,'markersize',10,'MarkerFaceColor',colorSet(1,:));
    

semilogy(zr,galhist.inGal.gasMass(2,:)./1e10,greyType,'color',grey) 
hold on
h(2)=semilogy(zr(ii:end),galhist.inGal.gasMass(2,ii:end)./1e10,'-',...
    'color',colorSet(3,:),'linewidth',1.5,'DisplayName','Gas');
semilogy(zr(ii),galhist.inGal.gasMass(2,ii)./1e10,'d','color',colorSet(3,:),...
    'linewidth',1.5,'markersize',10,'MarkerFaceColor',colorSet(3,:));
hl=legend(h);
set(hl,'Interpreter','latex','Fontsize',14);



ylabelmine('Gas / Stellar Mass $[\mathrm{10^{10}M_\odot}]$')

xlabelmine('redshift')
set(gca,'Fontsize',14)

%% gas components figure 

nanMask=galhist.inGal.gasMass(2,:)==0;

sfg=galhist.inGal.sfGas.mass./galhist.inGal.gasMass(2,:);sfg(nanMask)=0;
cdn=galhist.inGal.coldDenseGas.mass./galhist.inGal.gasMass(2,:);cdn(nanMask)=0;
cdi=galhist.inGal.coldDiluteGas.mass./galhist.inGal.gasMass(2,:);cdi(nanMask)=0;
wrm=(galhist.inGal.warmHotGas.mass-galhist.inGal.hotGas.mass)./galhist.inGal.gasMass(2,:);wrm(nanMask)=0;
hot=galhist.inGal.hotGas.mass./galhist.inGal.gasMass(2,:);hot(nanMask)=0;

figure 
h=[];

yyaxis left 

plot(zr,sfg,greyType,'color',grey,'linewidth',1);
hold on
plot(zr,cdn,greyType,'color',grey,'linewidth',1);
plot(zr,cdi,greyType,'color',grey,'linewidth',1);
plot(zr,wrm,greyType,'color',grey,'linewidth',1);
plot(zr,hot,greyType,'color',grey,'linewidth',1);



h(1)=plot(zr(ii:end),sfg(ii:end),'-','color',colorSet(4,:),'linewidth',1.5,'DisplayName','SF');
h(2)=plot(zr(ii:end),cdn(ii:end),'-','color',colorSet(2,:),'linewidth',1.5,'DisplayName','Cold Dense');
h(3)=plot(zr(ii:end),cdi(ii:end),'-','color',colorSet(3,:),'linewidth',1.5,'DisplayName','Cold Dilute');
h(4)=plot(zr(ii:end),wrm(ii:end),'-','color',colorSet(5,:),'linewidth',1.5,'DisplayName','Warm');
h(5)=plot(zr(ii:end),hot(ii:end),'-','color',colorSet(1,:),'linewidth',1.5,'DisplayName','Hot');

plot(zr(ii),sfg(ii),'d','MarkerFaceColor',colorSet(4,:),'markersize',8)
plot(zr(ii),cdn(ii),'d','MarkerFaceColor',colorSet(2,:),'markersize',8)
plot(zr(ii),cdi(ii),'d','MarkerFaceColor',colorSet(3,:),'markersize',8)
plot(zr(ii),wrm(ii),'d','MarkerFaceColor',colorSet(5,:),'markersize',8)
plot(zr(ii),hot(ii),'d','MarkerFaceColor',colorSet(1,:),'markersize',8)



ylabelmine('Mass Fraction')

yyaxis right

h(6)=plot(zr,galhist.inGal.gasMass(2,:)./1e10,'k-.','linewidth',0.8,'DisplayName','gas mass');%,'color',grey) 
hold on 
plot(zr(ii),galhist.inGal.gasMass(2,ii)./1e10,'dk',...
    'markersize',10);%,'MarkerFaceColor',colorSet(3,:));
ylabelmine('Gas Mass $[\mathrm{10^{10}M_\odot}]$')


hl=legend(h);
set(hl,'Interpreter','latex','Fontsize',14);

xlabelmine('redshift')
set(gca,'Fontsize',14)


%% bh 

EQM=galhist.inGal.cumEngQM;
ERM=galhist.inGal.cumEngRM;



figure 
h=[];


%yyaxis left

semilogy(zr,EQM,greyType,'color',grey,'linewidth',1);
hold on
semilogy(zr,ERM,greyType,'color',grey,'linewidth',1);

h(1)=semilogy(zr(ii:end),EQM(ii:end),'-','color',colorSet(2,:),'linewidth',1.5,'DisplayName','HAM');
h(2)=semilogy(zr(ii:end),ERM(ii:end),'-','color',colorSet(1,:),'linewidth',1.5,'DisplayName','LAM');

semilogy(zr(ii),EQM(ii),'d','MarkerFaceColor',colorSet(2,:),'markersize',8)
semilogy(zr(ii),ERM(ii),'d','MarkerFaceColor',colorSet(1,:),'markersize',8)

hl=legend(h);
set(hl,'Interpreter','latex','Fontsize',14,'location','NorthEast');

ylabelmine('Cumulative Energy $[\mathrm{10^{53}\,ergs}]$')


%yyaxis right 

xlabelmine('redshift')
set(gca,'Fontsize',14)







% 
% 
% zr2=galhist.zredExt;
% efac=illustris.utils.get_factorByRedshift('EnergyDissipationUnit',zr2);
% 
% edGal=fliplr(cumsum(fliplr(galhist.inGal.EnergyDissipation(2,:).*efac)))./1e53;
% edCGM=fliplr(cumsum(fliplr(galhist.inCGM.EnergyDissipation(2,:).*efac)))./1e53;
% edCGM2=fliplr(cumsum(fliplr((galhist.inSub.EnergyDissipation(2,:)-galhist.inCGM.EnergyDissipation(2,:)).*efac)))./1e53;
% 
% 
% figure
% iii=4;
% h1=[];
% 
% semilogy(zr2,edGal,greyType,'color',grey,'linewidth',1);
% hold on
% semilogy(zr2,edCGM,greyType,'color',grey,'linewidth',1);
% semilogy(zr2,edCGM2,greyType,'color',grey,'linewidth',1);
% 
% h1(1)=semilogy(zr2(iii:end),edGal(iii:end),'-','color',colorSet(2,:),'linewidth',1.5,'DisplayName','Gal');
% h1(2)=semilogy(zr2(iii:end),edCGM(iii:end),'-','color',colorSet(1,:),'linewidth',1.5,'DisplayName','$\mathrm{CGM_{in}}$');
% h1(3)=semilogy(zr2(iii:end),edCGM2(iii:end),'-','color',colorSet(3,:),'linewidth',1.5,'DisplayName','$\mathrm{CGM_{out}}$');
% 
% semilogy(zr2(iii),edGal(iii),'d','MarkerFaceColor',colorSet(2,:),'markersize',8)
% semilogy(zr2(iii),edCGM(iii),'d','MarkerFaceColor',colorSet(1,:),'markersize',8)
% semilogy(zr2(iii),edCGM2(iii),'d','MarkerFaceColor',colorSet(3,:),'markersize',8)
% 
% hl1=legend(h1);
% set(hl1,'Interpreter','latex','Fontsize',14,'location','SouthEast');
% 
% ylabelmine('Cumulative Energy $[10^{53}\,ergs]$')


%% track in mass - ssfr plot 

global DEFAULT_MATFILE_DIR
load([DEFAULT_MATFILE_DIR '/cooling_times_z0_TNG100.mat'])
load([DEFAULT_MATFILE_DIR '/tng100_z0_fofs_subs.mat'])
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);

centralMask= subsInfo.isCentral(tCoolStruct.galMask);
ssfr = illustris.utils.calc_ssfr(subs);
ssfr=ssfr(tCoolStruct.galMask);
galMass=tCoolStruct.galMass(tCoolStruct.galMask);

ssfrLab='$\log \mathrm{sSFR},[\mathrm{yr^{-1}}] $';
mstarLab='$\log M_\star\,[\mathrm{M_\odot}]$';

filt=fspecial('disk',8);
xdata=log10(galMass(centralMask));
ydata=log10(ssfr(centralMask));
popCont=plot_population_contour(xdata,ydata,'smooth',filt,'noplot');
xl=[9 12.5];
yl=[-16.5 -8];

figure

ax1=axes;
plot(log10(galhist.stellarMass),log10(galhist.ssfr),':','color','k');

hold on 
plot(log10(galhist.stellarMass(i:end)),log10(galhist.ssfr(i:end)),'-','color',colorSet(2,:),'linewidth',1.5)
plot(log10(galhist.stellarMass(i)),log10(galhist.ssfr(i)),'d','color',colorSet(2,:),...
    'linewidth',1.5,'markersize',10,'MarkerFaceColor',colorSet(2,:));

xlim(xl);ylim(yl);
grid
xlabelmine(mstarLab,12);
ylabelmine(ssfrLab,12);
set(gca,'fontsize',12)


ax2=axes;
contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','of','LineColor',grey,...
    'LevelList',[5:10:95 100],'Fill','off','linestyle','-');
xlim(xl);ylim(yl);
set(ax2,'position',get(ax1,'position'));
linkaxes([ax1,ax2])
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];


%% cooling time 
% tcG=interp1(zpl2,galhist.inGal.meanTcMW(1,:),zpl,'pchip');
% tcC=interp1(zpl2,galhist.inCGM.meanTcMW(1,:),zpl,'pchip');
% tcS=interp1(zpl2,galhist.inSub.meanTcMW(1,:),zpl,'pchip');

figure 

h=[];
loglog(zpl,tcG,greyType,'color',grey)
hold on
loglog(zpl,tcC,greyType,'color',grey)
loglog(zpl,tcS,greyType,'color',grey)

h(1)=loglog(zpl(i:end),tcG(i:end),'-','color',colorSet(2,:),'DisplayName','Galaxy');
h(2)=loglog(zpl(i:end),tcC(i:end),'-','color',colorSet(3,:),'DisplayName','Inner');
h(3)=loglog(zpl(i:end),tcS(i:end),'-','color',colorSet(1,:),'DisplayName','Outer');

indx0=find(zpl2==zpl(i));
if ~isempty(indx0)
    indx=indx0;
end
loglog(zpl2(indx:end),galhist.inGal.meanTcMW(1,indx:end),...
    'd','color',colorSet(2,:),'linewidth',1.5,'markersize',10,...
    'MarkerFaceColor',colorSet(2,:));

loglog(zpl2(indx:end),galhist.inCGM.meanTcMW(1,indx:end),...
    'd','color',colorSet(3,:),'linewidth',1.5,'markersize',10,...
    'MarkerFaceColor',colorSet(3,:));

loglog(zpl2(indx:end),galhist.inSub.meanTcMW(1,indx:end),...
    'd','color',colorSet(1,:),'linewidth',1.5,'markersize',10,...
    'MarkerFaceColor',colorSet(1,:));

hl=legend(h);
set(hl,'Interpreter','latex','Fontsize',14,'location','NorthEast');


grid minor

ylabelmine('Cooling Time')


%% mach 

MaG2=galhist.inGal.meanMachEW(2,:);
MaI2=galhist.inCGM.meanMachEW(2,:);
MaS2=galhist.inSub.meanMachEW(2,:);

% MaG2=MaG2(~isnan(MaG2));zpG=zpl2(~isnan(MaG2));
% MaI2=MaI2(~isnan(MaI2));zpI=zpl2(~isnan(MaI2));
% MaS2=MaS2(~isnan(MaS2));zpS=zpl2(~isnan(MaI2));


mI=galhist.inCGM.EnergyDissipation(2,:);
mS=galhist.inSub.EnergyDissipation(2,:);

MaO2= (MaS2.*mS -   MaI2.*mI )./(mS-mI);  % find the value for outer halo 

MaG=interp1(zpl2,MaG2,zpl,'pchip');
MaI=interp1(zpl2,MaI2,zpl,'pchip');
MaO=interp1(zpl2,MaO2,zpl,'pchip');

clear MaS2 mI mS

figure 
loglog(zpl,MaG,greyType,'color',grey)
    hold on
    loglog(zpl,MaI,greyType,'color',grey)
    loglog(zpl,MaO,greyType,'color',grey)
    
    loglog(zpl(i:end),MaG(i:end),'-','color',colorSet(2,:));
    loglog(zpl(i:end),MaI(i:end),'-','color',colorSet(3,:));
    loglog(zpl(i:end),MaO(i:end),'-','color',colorSet(1,:));
    
    indx0=find(zpl2==zpl(i));
    if ~isempty(indx0)
        indx=indx0;
    end
    
    h=[];
    h(1)=loglog(zpl2(indx:end),MaG2(1,indx:end),...
        'd','color',colorSet(2,:),'linewidth',1.5,'markersize',10,...
        'MarkerFaceColor',colorSet(2,:),'DisplayName','Gal');
    
    h(2)=loglog(zpl2(indx:end),MaI2(1,indx:end),...
        'd','color',colorSet(3,:),'linewidth',1.5,'markersize',10,...
        'MarkerFaceColor',colorSet(3,:),'DisplayName','Inner');
    
    h(3)=loglog(zpl2(indx:end),MaO2(1,indx:end),...
        'd','color',colorSet(1,:),'linewidth',1.5,'markersize',10,...
        'MarkerFaceColor',colorSet(1,:),'DisplayName','Outer');
    
    hl=legend(h);
    set(hl,'Interpreter','latex','Fontsize',12,'location','NorthEast');

    
    grid minor
    ylabelmine('Mach')
    xlabelmine('$\log(1+z)$')
    



%% Ediss

EdG2=galhist.inGal.EnergyDissipation(2,:);
mG=galhist.inGal.gasMass(2,fullMask);
EdI2=galhist.inCGM.EnergyDissipation(2,:);
mI=galhist.inCGM.gasMass(2,fullMask);
EdS2=galhist.inSub.EnergyDissipation(2,:);
mS=galhist.inSub.gasMass(2,fullMask);

EdO2= EdS2- EdI2;  % find the value for outer halo 
mO=mS-mI;

EdG2=EdG2./mG;
EdI2=EdI2./mI;
EdO2=EdO2./mO;


EdG=interp1(zpl2,EdG2,zpl,'pchip');
EdI=interp1(zpl2,EdI2,zpl,'pchip');
EdO=interp1(zpl2,EdO2,zpl,'pchip');

clear EdS2 mI mS mG

figure 
loglog(zpl,EdG,greyType,'color',grey)
    hold on
    loglog(zpl,EdI,greyType,'color',grey)
    loglog(zpl,EdO,greyType,'color',grey)
    
    loglog(zpl(i:end),EdG(i:end),'-','color',colorSet(2,:));
    loglog(zpl(i:end),EdI(i:end),'-','color',colorSet(3,:));
    loglog(zpl(i:end),EdO(i:end),'-','color',colorSet(1,:));
    
    indx0=find(zpl2==zpl(i));
    if ~isempty(indx0)
        indx=indx0;
    end
    
    h=[];
    h(1)=loglog(zpl2(indx:end),EdG2(1,indx:end),...
        'd','color',colorSet(2,:),'linewidth',1.5,'markersize',10,...
        'MarkerFaceColor',colorSet(2,:),'DisplayName','Gal');
    
    h(2)=loglog(zpl2(indx:end),EdI2(1,indx:end),...
        'd','color',colorSet(3,:),'linewidth',1.5,'markersize',10,...
        'MarkerFaceColor',colorSet(3,:),'DisplayName','Inner');
    
    h(3)=loglog(zpl2(indx:end),EdO2(1,indx:end),...
        'd','color',colorSet(1,:),'linewidth',1.5,'markersize',10,...
        'MarkerFaceColor',colorSet(1,:),'DisplayName','Outer');
    
    hl=legend(h);
    set(hl,'Interpreter','latex','Fontsize',12,'location','NorthEast');

    
    grid 
    ylabelmine('$E_{\mathrm{dis}}/ M_{\mathrm{g}}\,[\mathrm{erg\,yr^{-1}\,M_\odot}]$')
    xlabelmine('$\log(1+z)$')
    




%% mach / E components 

  
    %gmass=galhist.inGal.gasMass(2,:);
    %nanMask=gmass==0;
    %name1={'sfg','cdn','cdi','wrm','hot'}
    %name2={'sfGas','coldDenseGas','coldDiluteGas','warmHotGas','hotGas'}
    
    comp.sfg2=galhist.inGal.sfGas.EnergyDissipation./...
        galhist.inGal.sfGas.mass(fullMask);
    comp.sfg2(isnan(comp.sfg2))=0;
    
    comp.cdn2=galhist.inGal.coldDenseGas.EnergyDissipation./...
        galhist.inGal.coldDenseGas.mass(fullMask);
    comp.cdn2(isnan(comp.cdn2))=0;
    
    comp.cdi2=galhist.inGal.coldDiluteGas.EnergyDissipation./...
        galhist.inGal.coldDiluteGas.mass(fullMask);
    comp.cdi2(isnan(comp.cdi2))=0;
    
    comp.wrm2=(galhist.inGal.warmHotGas.EnergyDissipation-...
        galhist.inGal.hotGas.EnergyDissipation)./...
        (galhist.inGal.warmHotGas.mass(fullMask)-...
        galhist.inGal.hotGas.mass(fullMask));
    comp.wrm2(isnan(comp.wrm2))=0;
    
    comp.hot2=galhist.inGal.hotGas.EnergyDissipation./...
        galhist.inGal.hotGas.mass(fullMask);
    comp.hot2(isnan(comp.hot2))=0;
    
    %galMovie_plot_gasMass(zpl,comp,i,'yleft')
    
    comp.sfg=interp1(zpl2,comp.sfg2,zpl,'pchip');
    comp.cdn=interp1(zpl2,comp.cdn2,zpl,'pchip');
    comp.cdi=interp1(zpl2,comp.cdi2,zpl,'pchip');
    comp.wrm=interp1(zpl2,comp.wrm2,zpl,'pchip');
    comp.hot=interp1(zpl2,comp.hot2,zpl,'pchip');
    


figure 
loglog(zpl,EdG,greyType,'color',grey)
    hold on
    loglog(zpl,EdI,greyType,'color',grey)
    loglog(zpl,EdO,greyType,'color',grey)
    
    loglog(zpl(i:end),EdG(i:end),'-','color',colorSet(2,:));
    loglog(zpl(i:end),EdI(i:end),'-','color',colorSet(3,:));
    loglog(zpl(i:end),EdO(i:end),'-','color',colorSet(1,:));
    
    indx0=find(zpl2==zpl(i));
    if ~isempty(indx0)
        indx=indx0;
    end
    
    h=[];
    h(1)=loglog(zpl2(indx:end),EdG2(1,indx:end),...
        'd','color',colorSet(2,:),'linewidth',1.5,'markersize',10,...
        'MarkerFaceColor',colorSet(2,:),'DisplayName','Gal');
    
    h(2)=loglog(zpl2(indx:end),EdI2(1,indx:end),...
        'd','color',colorSet(3,:),'linewidth',1.5,'markersize',10,...
        'MarkerFaceColor',colorSet(3,:),'DisplayName','Inner');
    
    h(3)=loglog(zpl2(indx:end),EdO2(1,indx:end),...
        'd','color',colorSet(1,:),'linewidth',1.5,'markersize',10,...
        'MarkerFaceColor',colorSet(1,:),'DisplayName','Outer');
    
    hl=legend(h);
    set(hl,'Interpreter','latex','Fontsize',12,'location','NorthEast');

    
    grid 
    ylabelmine('$E_{\mathrm{dis}}/ M_{\mathrm{g}}\,[\mathrm{erg\,yr^{-1}\,M_\odot}]$')
    xlabelmine('$\log(1+z)$')
    


