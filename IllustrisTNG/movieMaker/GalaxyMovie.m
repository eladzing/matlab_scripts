%% perliminary

if perlimFlag
    
    global cosmoStruct
    
    global DEFAULT_MATFILE_DIR
    load([DEFAULT_MATFILE_DIR 'cooling_times_z0_TNG100.mat'])
    load([DEFAULT_MATFILE_DIR 'tng100_z0_fofs_subs.mat'])
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
    
    colorSet=brewermap(9,'Set1');
    grey=colorSet(9,:);
    greyType=':';
    
    %% cooling time
    
    % find full snaps
    
    zred=galhist.zred;
    zpl=zred+1;
    zpl2=galhist.zredExt+1;
    fullMask=false(size(zred));
    
    for j=1:length(zred)
        fullMask(j)=any(zpl2==zpl(j));
    end
    
    
    tcG2=galhist.inGal.meanTcMW(1,:);
    tcI2=galhist.inCGM.meanTcMW(1,:);
    tcS2=galhist.inSub.meanTcMW(1,:);
    
    mI=galhist.inCGM.gasMass(1,fullMask);
    mS=galhist.inSub.gasMass(1,fullMask);
    
    tcO2= (tcS2.*mS -   tcI2.*mI )./(mS-mI);  % find the value for outer halo
    
    tcG=interp1(zpl2,tcG2,zpl,'pchip');
    tcI=interp1(zpl2,tcI2,zpl,'pchip');
    tcO=interp1(zpl2,tcO2,zpl,'pchip');
    
    clear tcS2 mI mS
end


%% prepare movie frames



%set(hf,'position',[220 605 1642 718]);


cnt=1;


for i=length(galhist.zred):-1:1
    hf=figure;
    set(hf,'position',[32  68  2504  1261]);%76  71  2151 1267]);
    
    %% birds
    bird1=squeeze(galhistBird.inGal.bird(:,:,i));
    bird2=squeeze(galhistBird.inCGMBird.bird(:,:,i));
    bird3=squeeze(galhistBird.inSubBird.bird(:,:,i)-galhistBird.inCGMBird.bird(:,:,i));
    
    subplot(4,4,[1 5])
    illustris.plots.plot_phaseDiagram(galhistBird.birdXlim,galhistBird.birdYlim,...
        bird1,'fig',hf,'caxis',[7.5 10]);
    titlemine(sprintf('Gal $z=%4.2f$',galhist.zred(i)));
    
    subplot(4,4,[2 6])
    illustris.plots.plot_phaseDiagram(galhistBird.birdXlim,galhistBird.birdYlim,...
        bird2,'fig',hf,'caxis',[7.5 10]);
    titlemine(sprintf('Inner Halo $z=%4.2f$',galhist.zred(i)));
    %fprintf('i= %i \n',i)
    
    subplot(4,4,[3 7])
    illustris.plots.plot_phaseDiagram(galhistBird.birdXlim,galhistBird.birdYlim,...
        bird3,'fig',hf,'caxis',[7.5 10]);
    titlemine(sprintf('Outer Halo $z=%4.2f$',galhist.zred(i)));
    
    %% track in mass - ssfr
    
    subplot(4,4,[4 8])
    
    ax1=gca;% axes;
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
    
    
    ax2= axes;
    contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','of','LineColor',grey,...
        'LevelList',[5:10:95 100],'Fill','off','linestyle','-');
    xlim(xl);ylim(yl);
    set(ax2,'position',get(ax1,'position'));
    linkaxes([ax1,ax2])
    ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];
    
    
    %% mass fractoions Gal
    
    
    
    % in Galaxy
    subplot(4,4,9)
    
    gmass=galhist.inGal.gasMass(2,:);
    nanMask=gmass==0;
    
    comp.gmass=gmass;
    comp.sfg=galhist.inGal.sfGas.mass./gmass;comp.sfg(nanMask)=0;
    comp.cdn=galhist.inGal.coldDenseGas.mass./gmass;comp.cdn(nanMask)=0;
    comp.cdi=galhist.inGal.coldDiluteGas.mass./gmass;comp.cdi(nanMask)=0;
    comp.wrm=(galhist.inGal.warmHotGas.mass-galhist.inGal.hotGas.mass)./gmass;comp.wrm(nanMask)=0;
    comp.hot=galhist.inGal.hotGas.mass./gmass;comp.hot(nanMask)=0;
    
    galMovie_plot_gasMass(zpl,comp,i,'yleft')
    
    
    % inner sub
    subplot(4,4,10)
    
    gmass=galhist.inCGM.gasMass(2,:);
    nanMask=gmass==0;
    
    comp.gmass=gmass;
    comp.sfg=galhist.inCGM.sfGas.mass./gmass;comp.sfg(nanMask)=0;
    comp.cdn=galhist.inCGM.coldDenseGas.mass./gmass;comp.cdn(nanMask)=0;
    comp.cdi=galhist.inCGM.coldDiluteGas.mass./gmass;comp.cdi(nanMask)=0;
    comp.wrm=(galhist.inCGM.warmHotGas.mass-galhist.inCGM.hotGas.mass)./gmass;comp.wrm(nanMask)=0;
    comp.hot=galhist.inCGM.hotGas.mass./gmass;comp.hot(nanMask)=0;
    
    galMovie_plot_gasMass(zpl,comp,i)
    
    % outer sub
    
    subplot(4,4,11)
    
    gmass=galhist.inSub.gasMass(2,:)-galhist.inCGM.gasMass(2,:);
    nanMask=gmass==0;
    
    comp.gmass=gmass;
    comp.sfg=(galhist.inSub.sfGas.mass-galhist.inCGM.sfGas.mass)./gmass;comp.sfg(nanMask)=0;
    comp.cdn=(galhist.inSub.coldDenseGas.mass-galhist.inCGM.coldDenseGas.mass)./gmass;comp.cdn(nanMask)=0;
    comp.cdi=(galhist.inSub.coldDiluteGas.mass-galhist.inCGM.coldDiluteGas.mass)./gmass;comp.cdi(nanMask)=0;
    comp.wrm=((galhist.inSub.warmHotGas.mass-galhist.inCGM.warmHotGas.mass)-...
        (galhist.inSub.hotGas.mass-galhist.inCGM.hotGas.mass))./gmass;comp.wrm(nanMask)=0;
    comp.hot=(galhist.inSub.hotGas.mass-galhist.inCGM.hotGas.mass)./gmass;comp.hot(nanMask)=0;
    
    galMovie_plot_gasMass(zpl,comp,i,'yright','legend')
    
    
    %% sfr & gas mass
    subplot(4,4,13)
    
    h=[];
    yyaxis left
    semilogx(zpl,log10(galhist.sfr),greyType,'color',grey)
    hold on
    semilogx(zpl(i:end),log10(galhist.sfr(i:end)),'-','color',colorSet(2,:),'linewidth',1.5)
    semilogx(zpl(i),log10(galhist.sfr(i)),'d','color',colorSet(2,:),...
        'linewidth',1.5,'markersize',10,'MarkerFaceColor',colorSet(2,:));
    
    ylabelmine(' log sfr $[\mathrm{M_\odot\,yr^{-1}}]$')
    %ylim([1e-14 1e-8])
    yyaxis right
    
    semilogx(zpl,log10(galhist.stellarMass./1e10),greyType,'color',grey)
    hold on
    h(1)=semilogx(zpl(i:end),log10(galhist.stellarMass(i:end)./1e10),'-',...
        'color',colorSet(1,:),'linewidth',1.5,'DisplayName','Stellar');
    semilogx(zpl(i),log10(galhist.stellarMass(i)./1e10),'d','color',colorSet(1,:),...
        'linewidth',1.5,'markersize',10,'MarkerFaceColor',colorSet(1,:));
    
    
    semilogx(zpl,log10(galhist.inGal.gasMass(2,:)./1e10),greyType,'color',grey)
    hold on
    h(2)=semilogx(zpl(i:end),log10(galhist.inGal.gasMass(2,i:end)./1e10),'-',...
        'color',colorSet(3,:),'linewidth',1.5,'DisplayName','Gas');
    semilogx(zpl(i),log10(galhist.inGal.gasMass(2,i)./1e10),'d','color',colorSet(3,:),...
        'linewidth',1.5,'markersize',10,'MarkerFaceColor',colorSet(3,:));
    hl=legend(h);
    set(hl,'Interpreter','latex','Fontsize',14);
    
    
    
    ylabelmine(' log Gas / Stellar Mass $[\mathrm{10^{10}M_\odot}]$')
    
    xlabelmine('$\log(1+z)$')
    set(gca,'Fontsize',14)
    
    
    %% BH stuff
    subplot(4,4,14)
    
    EQM=galhist.inGal.cumEngQM;
    ERM=galhist.inGal.cumEngRM;
    
    dedtQ=zeros(size(EQM));
    dedtR=dedtQ;
    
    % derivative
    time=redshift2time(zred,'cosmo',cosmoStruct);
    [dum,~]=derive1(EQM,time.age);
    dedtQ(2:end-1)=log10(dum);
    [dum,tim]=derive1(ERM,time.age);
    dedtR(2:end-1)=log10(dum);
    
    
    h=[];
    
    yyaxis left
    
    semilogx(zpl(2:end-1),dedtQ(2:end-1),greyType,'color',grey,'linewidth',1);
    hold on
    semilogx(zpl(2:end-1),dedtR(2:end-1),greyType,'color',grey,'linewidth',1);
    
    if i==1
        ii=2;
    elseif i==length(zpl)
        ii=length(zpl)-1;
    else
        ii=i;
    end
    
    
    h(1)=semilogx(zpl(ii:end),dedtQ(ii:end),'-','color',colorSet(2,:),'linewidth',1.5,'DisplayName','HAR');
    h(2)=semilogx(zpl(ii:end),dedtR(ii:end),'-','color',colorSet(1,:),'linewidth',1.5,'DisplayName','LAR');
    
    semilogx(zpl(ii),dedtQ(ii),'d','MarkerFaceColor',colorSet(2,:),'markersize',8)
    semilogx(zpl(ii),dedtR(ii),'d','MarkerFaceColor',colorSet(1,:),'markersize',8)
    
    ylabelmine('$\log \dot{E} [\mathrm{10^{53}\,ergs\,Gyr^{-1}}]$')
    
    
    yyaxis right
    EQM=log10(EQM);
    ERM=log10(ERM);
    
    semilogx(zpl,EQM,greyType,'color',grey,'linewidth',1);
    hold on
    semilogx(zpl,ERM,greyType,'color',grey,'linewidth',1);
    
    semilogx(zpl(i:end),EQM(i:end),'--','color',colorSet(2,:),'linewidth',1.5,'DisplayName','HAM');
    semilogx(zpl(i:end),ERM(i:end),'--','color',colorSet(1,:),'linewidth',1.5,'DisplayName','LAM');
    
    semilogx(zpl(i),EQM(i),'d','MarkerFaceColor',colorSet(2,:),'markersize',8)
    semilogx(zpl(i),ERM(i),'d','MarkerFaceColor',colorSet(1,:),'markersize',8)
    
    ylabelmine('log Cumulative Energy $[\mathrm{10^{53}\,ergs}]$')
    
    
    hl=legend(h);
    set(hl,'Interpreter','latex','Fontsize',14,'location','SouthWest');
    %yyaxis right
    
    
    
    xlabelmine('$\log(1+z)$')
    set(gca,'Fontsize',12)
    
    
    %%  cooling time
    subplot(4,4,16)
    loglog(zpl,tcG,greyType,'color',grey)
    hold on
    loglog(zpl,tcI,greyType,'color',grey)
    loglog(zpl,tcO,greyType,'color',grey)
    
    loglog(zpl(i:end),tcG(i:end),'-','color',colorSet(2,:));
    loglog(zpl(i:end),tcI(i:end),'-','color',colorSet(3,:));
    loglog(zpl(i:end),tcO(i:end),'-','color',colorSet(1,:));
    
    indx0=find(zpl2==zpl(i));
    if ~isempty(indx0)
        indx=indx0;
    end
    
    h=[];
    h(1)=loglog(zpl2(indx:end),tcG2(1,indx:end),...
        'd','color',colorSet(2,:),'linewidth',1.5,'markersize',10,...
        'MarkerFaceColor',colorSet(2,:),'DisplayName','Gal');
    
    h(2)=loglog(zpl2(indx:end),tcI2(1,indx:end),...
        'd','color',colorSet(3,:),'linewidth',1.5,'markersize',10,...
        'MarkerFaceColor',colorSet(3,:),'DisplayName','Inner');
    
    h(3)=loglog(zpl2(indx:end),tcO2(1,indx:end),...
        'd','color',colorSet(1,:),'linewidth',1.5,'markersize',10,...
        'MarkerFaceColor',colorSet(1,:),'DisplayName','Outer');
    
    hl=legend(h);
    set(hl,'Interpreter','latex','Fontsize',12,'location','SouthWest');
    
    
    grid
    ylabelmine('Cooling Time [Gyr]')
    xlabelmine('$\log(1+z)$')
    
    
    %% 
    
    
    
    %% add text 
    str1=sprintf('$z= %s',
    
    
    
    pos=[1700 500 200 200];
    MyBox = uicontrol('style','text');
    set(MyBox,'String','Here is a lot more information')
    set(MyBox,'Position',pos);
    %% save Frame
    
    
    F(cnt)=getframe(gcf);
    
    
    cnt=cnt+1;
    
    
    close(gcf)
    
    
end

%% write movie to file

global DEFAULT_MOVIE_DIR

movieName=[DEFAULT_MOVIE_DIR '/galtest2' ];

vid = VideoWriter(movieName);

set(vid,'FrameRate',2);
open(vid)

for i=1:length(F)
    
    writeVideo(vid,F(i));
end


close(vid)

%% convert to mp4

comm=sprintf('ffmpeg -i %s.avi -acodec libfaac -b:a 128k -vcodec mpeg4 -b:v 1200k -flags +aic+mv4 %s.mp4',movieName,movieName);

system(comm);
system(sprintf('rm -f %s.avi',movieName));





% set properties









%
% writeVideo(vid,F)