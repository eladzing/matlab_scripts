%% perliminary

fprintf('*** Make sure you have uploaded the correct history structures  *** \n')

fullMovieFlag=true;

movieFolder='subbox0_list';
baseName='cgmQuenched';


if perlimFlag
    
    
    
    global cosmoStruct
    %     global DEFAULT_MOVIE_DIR
    
    global DEFAULT_MATFILE_DIR
%    load([DEFAULT_MATFILE_DIR '/cooling_times_z0_TNG100.mat'])
%    load([DEFAULT_MATFILE_DIR '/tng100_z0_fofs_subs.mat'])
%    subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
    
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
    
    clear tCoolStruct;
    
    
    % filename for output
    global simDisplayName
    fname=sprintf('BCG_Movies_%s',simDisplayName);
    clear subs fofs subsInfo
    
    %bird stuff
    guide(1).tmp=[7 8];
    guide(1).ro=[0 NaN];
    guide(1).a=-1;
    guide(1).c=guide(1).ro(1)-guide(1).a.*guide(1).tmp(1);
    guide(1).ro(2)=guide(1).a.*guide(1).tmp(2)+guide(1).c;
    
    guide(2).tmp=[4 5];
    guide(2).ro=[-5 NaN];
    guide(2).a=1.5;
    guide(2).c=guide(2).ro(1)-guide(2).a.*guide(2).tmp(1);
    guide(2).ro(2)=guide(2).a.*guide(2).tmp(2)+guide(2).c;
    
    birdCaxis=[7.5 10];
    bartag='Mass';
    bmap='*Spectral';
    cmap=brewermap(256,bmap);cmap(1,:)=[1 1 1];
    
    % load matfile
    
    %     MGH=matfile([DEFAULT_MATFILE_DIR '/BCG_history_z0_TNG100.mat']);
    %     MBH=matfile([DEFAULT_MATFILE_DIR '/BCG_history_birds_z0_TNG100.mat']);
    %     done=MGH.done;
    perlimFlag=false;
    
end


%% get the right galaxy

galInds=list;    [328251 350183 351071 356906 364479] +1;
%galInds=find(done);
   

for k=2:length(galInds)
    galID= galInds(k);
    %galhist=galHistory(1,galID);
    %galhistBird=galBirdHistory(1,galID);
    galhist=galHistory(1,k);
    galhistBird=galBirdHistory(1,k);
    
%     galhist=MGH.galHistory(1,galID);
%     galhistBird=MBH.galBirdHistory(1,galID);
    
    fprintf('Preparing movie %i out of %i \n',k,length(galInds))
    
    
    
    %% cooling time
    
    % find full snaps
    
    zred=galhist.zred;
    zpl=zred+1;
    zpl2=galhist.zredExt+1;
    fullMask=false(size(zred));
    
    % get times
    zredTimes=redshift2time(zred,'cosmo',cosmoStruct);
    zredTimes=zredTimes.age;
    indx=length(galhist.zredExt); % initialize indx
    
    
    
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
    
    % energy dissipation
    mG=galhist.inGal.gasMass(2,fullMask)./1e10;
    mI=galhist.inCGM.gasMass(2,fullMask)./1e10;
    mO=(galhist.inSub.gasMass(2,fullMask)-mI)./1e10;
    
    EdG0=galhist.inGal.EnergyDissipation(2,:);
    EdG0=cat(1,EdG0,EdG0./mG);
    
    EdI0=galhist.inCGM.EnergyDissipation(2,:);
    EdI0=cat(1,EdI0,EdI0./mI);
    
    EdO0=galhist.inSub.EnergyDissipation(2,:)- EdI0;
    EdO0=cat(1,EdO0,EdO0./mO);
    
    EdG0(isnan(EdG0))=0;
    EdI0(isnan(EdI0))=0;
    EdO0(isnan(EdO0))=0;
    
    EdG=zeros(2,length(zpl));
    EdI=zeros(2,length(zpl));
    EdO=zeros(2,length(zpl));
    
    % energy dissipated
    for lm=1:2
        EdG(lm,:)=interp1(zpl2,EdG0(lm,:),zpl,'pchip');
        EdI(lm,:)=interp1(zpl2,EdI0(lm,:),zpl,'pchip');
        EdO(lm,:)=interp1(zpl2,EdO0(lm,:),zpl,'pchip');
        
    end
    
    clear mO mI mG
    
    
    % BH stuff
    EQM=galhist.inGal.cumEngQM;
    ERM=galhist.inGal.cumEngRM;
    
    dedtQ=zeros(size(EQM));
    dedtR=dedtQ;
    
    % derivative of energy
    
    [dum,~]=derive1(EQM,zredTimes);
    dedtQ(2:end-1)=log10(dum);
    [dum,tim]=derive1(ERM,zredTimes);
    dedtR(2:end-1)=log10(dum);
    
    EQM=log10(EQM);
    ERM=log10(ERM);
    clear dum
    
    % find virial stuff
    [rvir, mvir, tvir, vvir]=calculate_virials('mv',galhist.hostMass,...
        'zred',galhist.zred,'delv',200,'crit','cosmo',cosmoStruct);
    %% prepare movie frames
    
    cnt=1;
    if fullMovieFlag
        frameRange=length(galhist.zred):-1:1;
    else
        frameRange=1;
    end
    for i=frameRange  %length(galhist.zred):-1:1
        
        hf=figure('Visible','off');
        %hf=figure('Visible','on');
        
        hf.Position=[32  68  2504  1261];
        %76  71  2151 1267]);
        
        %% birds
        bird1=squeeze(galhistBird.inGal.bird(:,:,i));
        bird2=squeeze(galhistBird.inCGM.bird(:,:,i));
        bird3=squeeze(galhistBird.inSub.bird(:,:,i)-galhistBird.inCGM.bird(:,:,i));
        
        
        
        subplot(6,4,[1 5 9])
        %illustris.plots.plot_phaseDiagram(galhistBird.birdXlim,galhistBird.birdYlim,...
        %    bird1,'fig',hf,'caxis',[7.5 10]);
        
        imsc(galhistBird.birdYlim,galhistBird.birdXlim,...
            log10(bird1'),cmap);
        
        set(gca,'Ydir','normal','Fontsize',10)
        hold on
        for lm=1:2;plot(guide(lm).ro,guide(lm).tmp,'--k');end
        plot(galhistBird.birdYlim,log10(tvir(i)).*[1 1],'-.k')
        
        grid('on')
        
        colormap(cmap);
        hb=colorbar;barTitle(hb,bartag)
        caxis(birdCaxis)
        
        ylabelmine('$\log(T)\,[\mathrm{K}]$');
        xlabelmine('$\log(n)\,[\mathrm{cm^{-3}}]$');
        titlemine('Gal',20);    %titlemine(sprintf('Gal $z=%4.2f$',galhist.zred(i)));
        
        subplot(6,4,[2 6 10])
        %         illustris.plots.plot_phaseDiagram(galhistBird.birdXlim,galhistBird.birdYlim,...
        %             bird2,'fig',hf,'caxis',[7.5 10]);
        %titlemine(sprintf('Inner Halo $z=%4.2f$',galhist.zred(i)));
        imsc(galhistBird.birdYlim,galhistBird.birdXlim,...
            log10(bird2'),cmap);
        set(gca,'Ydir','normal','Fontsize',10)
        hold on
        for lm=1:2;plot(guide(lm).ro,guide(lm).tmp,'--k');end
        plot(galhistBird.birdYlim,log10(tvir(i)).*[1 1],'-.k')
        grid('on')
        
        colormap(cmap);
        hb=colorbar;barTitle(hb,bartag)
        caxis(birdCaxis)
        
        ylabelmine('$\log(T)\,[\mathrm{K}]$');
        xlabelmine('$\log(n)\,[\mathrm{cm^{-3}}]$');
        titlemine('Inner Halo',20);
        
        
        subplot(6,4,[3 7 11])
        %         illustris.plots.plot_phaseDiagram(galhistBird.birdXlim,galhistBird.birdYlim,...
        %             bird3,'fig',hf,'caxis',[7.5 10]);
        %titlemine(sprintf('Outer Halo $z=%4.2f$',galhist.zred(i)));
        imsc(galhistBird.birdYlim,galhistBird.birdXlim,...
            log10(bird3'),cmap);
        set(gca,'Ydir','normal','Fontsize',10)
        hold on
        for lm=1:2;plot(guide(lm).ro,guide(lm).tmp,'--k');end
        plot(galhistBird.birdYlim,log10(tvir(i)).*[1 1],'-.k')
        grid('on')
        
        colormap(cmap);
        hb=colorbar;barTitle(hb,bartag)
        caxis(birdCaxis)
        
        ylabelmine('$\log(T)\,[\mathrm{K}]$');
        xlabelmine('$\log(n)\,[\mathrm{cm^{-3}}]$');
        titlemine('Outer Halo',20);
        
        %% track in mass - ssfr
        
        subplot(6,4,[4 8 ])
        
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
        subplot(6,4,[13 17])
        
        gmass=galhist.inGal.gasMass(2,:);
        nanMask=gmass==0;
        
        comp.gmass=gmass;
        comp.sfg=galhist.inGal.sfGas.mass./gmass;comp.sfg(nanMask)=0;
        comp.cdn=galhist.inGal.coldDenseGas.mass./gmass;comp.cdn(nanMask)=0;
        comp.cdi=galhist.inGal.coldDiluteGas.mass./gmass;comp.cdi(nanMask)=0;
        comp.wrm=(galhist.inGal.warmHotGas.mass-galhist.inGal.hotGas.mass)./gmass;comp.wrm(nanMask)=0;
        comp.hot=galhist.inGal.hotGas.mass./gmass;comp.hot(nanMask)=0;
        
        galMovie_plot_gasMass(zpl,comp,i,'yleft','xl')
        
        
        % inner sub
        subplot(6,4,[14 18])
        
        gmass=galhist.inCGM.gasMass(2,:);
        nanMask=gmass==0;
        
        comp.gmass=gmass;
        comp.sfg=galhist.inCGM.sfGas.mass./gmass;comp.sfg(nanMask)=0;
        comp.cdn=galhist.inCGM.coldDenseGas.mass./gmass;comp.cdn(nanMask)=0;
        comp.cdi=galhist.inCGM.coldDiluteGas.mass./gmass;comp.cdi(nanMask)=0;
        comp.wrm=(galhist.inCGM.warmHotGas.mass-galhist.inCGM.hotGas.mass)./gmass;comp.wrm(nanMask)=0;
        comp.hot=galhist.inCGM.hotGas.mass./gmass;comp.hot(nanMask)=0;
        
        galMovie_plot_gasMass(zpl,comp,i,'xl')
        
        % outer sub
        
        subplot(6,4,[15 19])
        
        gmass=galhist.inSub.gasMass(2,:)-galhist.inCGM.gasMass(2,:);
        nanMask=gmass==0;
        
        comp.gmass=gmass;
        comp.sfg=(galhist.inSub.sfGas.mass-galhist.inCGM.sfGas.mass)./gmass;comp.sfg(nanMask)=0;
        comp.cdn=(galhist.inSub.coldDenseGas.mass-galhist.inCGM.coldDenseGas.mass)./gmass;comp.cdn(nanMask)=0;
        comp.cdi=(galhist.inSub.coldDiluteGas.mass-galhist.inCGM.coldDiluteGas.mass)./gmass;comp.cdi(nanMask)=0;
        comp.wrm=((galhist.inSub.warmHotGas.mass-galhist.inCGM.warmHotGas.mass)-...
            (galhist.inSub.hotGas.mass-galhist.inCGM.hotGas.mass))./gmass;comp.wrm(nanMask)=0;
        comp.hot=(galhist.inSub.hotGas.mass-galhist.inCGM.hotGas.mass)./gmass;comp.hot(nanMask)=0;
        
        galMovie_plot_gasMass(zpl,comp,i,'yright','legend','xl')
        
        
        %% sfr & gas mass
        subplot(6,4,[12 16])
        
        h=[];
        yyaxis left
        semilogx(zpl,log10(galhist.sfr),greyType,'color',grey)
        hold on
        h(1)=semilogx(zpl(i:end),log10(galhist.sfr(i:end)),'-','color',colorSet(2,:),'linewidth',1.5,...
            'DisplayName','SFR');
        semilogx(zpl(i),log10(galhist.sfr(i)),'d','color',colorSet(2,:),...
            'linewidth',1.5,'markersize',10,'MarkerFaceColor',colorSet(2,:));
        
        ylabelmine(' log sfr $[\mathrm{M_\odot\,yr^{-1}}]$')
        %ylim([1e-14 1e-8])
        yyaxis right
        
        semilogx(zpl,log10(galhist.stellarMass./1e10),greyType,'color',grey)
        hold on
        h(2)=semilogx(zpl(i:end),log10(galhist.stellarMass(i:end)./1e10),'-',...
            'color',colorSet(1,:),'linewidth',1.5,'DisplayName','Stellar');
        semilogx(zpl(i),log10(galhist.stellarMass(i)./1e10),'d','color',colorSet(1,:),...
            'linewidth',1.5,'markersize',10,'MarkerFaceColor',colorSet(1,:));
        
        
        semilogx(zpl,log10(galhist.inGal.gasMass(2,:)./1e10),greyType,'color',grey)
        hold on
        h(3)=semilogx(zpl(i:end),log10(galhist.inGal.gasMass(2,i:end)./1e10),'-',...
            'color',colorSet(3,:),'linewidth',1.5,'DisplayName','Gas');
        semilogx(zpl(i),log10(galhist.inGal.gasMass(2,i)./1e10),'d','color',colorSet(3,:),...
            'linewidth',1.5,'markersize',10,'MarkerFaceColor',colorSet(3,:));
        hl=legend(h);
        set(hl,'Interpreter','latex','Fontsize',14);
        
        
        
        ylabelmine(' log Gas / Stellar Mass $[\mathrm{10^{10}M_\odot}]$')
        
        xlabelmine('$\log(1+z)$')
        %set(gca,'Fontsize',12)
        
        
        %% BH stuff
        subplot(6,4,[20 24])
        
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
        
        h(1)=semilogx(zpl(ii:end-1),dedtQ(ii:end-1),'-','color',colorSet(2,:),'linewidth',1.5,'DisplayName','HAR');
        h(2)=semilogx(zpl(ii:end-1),dedtR(ii:end-1),'-','color',colorSet(1,:),'linewidth',1.5,'DisplayName','LAR');
        
        semilogx(zpl(ii),dedtQ(ii),'d','MarkerFaceColor',colorSet(2,:),'markersize',8)
        semilogx(zpl(ii),dedtR(ii),'d','MarkerFaceColor',colorSet(1,:),'markersize',8)
        
        ylabelmine('$\log \dot{E} [\mathrm{10^{53}\,ergs\,Gyr^{-1}}]$')
        
        
        yyaxis right
        
        
        
        semilogx(zpl,EQM,greyType,'color',grey,'linewidth',1);
        hold on
        semilogx(zpl,ERM,greyType,'color',grey,'linewidth',1);
        
        semilogx(zpl(i:end),EQM(i:end),'--','color',colorSet(2,:),'linewidth',1.5,'DisplayName','HAM');
        semilogx(zpl(i:end),ERM(i:end),'--','color',colorSet(1,:),'linewidth',1.5,'DisplayName','LAM');
        
        semilogx(zpl(i),EQM(i),'d','MarkerFaceColor',colorSet(2,:),'markersize',8)
        semilogx(zpl(i),ERM(i),'d','MarkerFaceColor',colorSet(1,:),'markersize',8)
        
        ylabelmine('log Cumulative Energy $[\mathrm{10^{53}\,ergs}]$')
        
        
        hl=legend(h);
        set(hl,'Interpreter','latex','Fontsize',14,'location','NorthEast');
        %yyaxis right
        
        
        
        xlabelmine('$\log(1+z)$')
        %set(gca,'Fontsize',12)
        
%         %%   Ediss
%         indx0=find(zpl2==zpl(i));
%         if ~isempty(indx0)
%             indx=indx0;
%         end
%         
%         subplot(6,4,15)
%         
%         yyaxis left
%         
%         semilogx(zpl,log10(EdG(1,:)),greyType,'color',grey)
%         hold on
%         semilogx(zpl,log10(EdI(1,:)),greyType,'color',grey)
%         semilogx(zpl,log10(EdO(1,:)),greyType,'color',grey)
%         
%         semilogx(zpl(i:end),log10(EdG(1,i:end)),'-','color',colorSet(2,:));
%         semilogx(zpl(i:end),log10(EdI(1,i:end)),'-','color',colorSet(3,:));
%         semilogx(zpl(i:end),log10(EdO(1,i:end)),'-','color',colorSet(1,:));
%         
%         
%         h=[];
%         h(1)=semilogx(zpl2(indx:end),log10(EdG0(1,indx:end)),...
%             'd','color',colorSet(2,:),'linewidth',1.5,'markersize',8,...
%             'MarkerFaceColor',colorSet(2,:),'DisplayName','Gal');
%         
%         h(2)=semilogx(zpl2(indx:end),log10(EdI0(1,indx:end)),...
%             'd','color',colorSet(3,:),'linewidth',1.5,'markersize',8,...
%             'MarkerFaceColor',colorSet(3,:),'DisplayName','Inner');
%         
%         h(3)=semilogx(zpl2(indx:end),log10(EdO0(1,indx:end)),...
%             'd','color',colorSet(1,:),'linewidth',1.5,'markersize',8,...
%             'MarkerFaceColor',colorSet(1,:),'DisplayName','Outer');
%         
%         
%         
%         %grid
%         ylabelmine('$E_{\mathrm{dis}}\,[10^{45} \mathrm{erg\,yr^{-1}}]$')
%         
%         yyaxis right
%         
%         semilogx(zpl,log10(EdG(2,:)),greyType,'color',grey)
%         hold on
%         semilogx(zpl,log10(EdI(2,:)),greyType,'color',grey)
%         semilogx(zpl,log10(EdO(2,:)),greyType,'color',grey)
%         
%         semilogx(zpl(i:end),log10(EdG(2,i:end)),'--','color',colorSet(2,:));
%         semilogx(zpl(i:end),log10(EdI(2,i:end)),'--','color',colorSet(3,:));
%         semilogx(zpl(i:end),log10(EdO(2,i:end)),'--','color',colorSet(1,:));
%         
%         
%         
%         semilogx(zpl2(indx:end),log10(EdG0(2,indx:end)),...
%             'd','color',colorSet(2,:),'linewidth',1.5,'markersize',4,...
%             'MarkerFaceColor',colorSet(2,:),'DisplayName','Gal');
%         
%         semilogx(zpl2(indx:end),log10(EdI0(2,indx:end)),...
%             'd','color',colorSet(3,:),'linewidth',1.5,'markersize',4,...
%             'MarkerFaceColor',colorSet(3,:),'DisplayName','Inner');
%         
%         semilogx(zpl2(indx:end),log10(EdO0(2,indx:end)),...
%             'd','color',colorSet(1,:),'linewidth',1.5,'markersize',4,...
%             'MarkerFaceColor',colorSet(1,:),'DisplayName','Outer');
%         
%         hl=legend(h);
%         set(hl,'Interpreter','latex','Fontsize',12,'location','SouthWest');
%         
%         
%         %grid
%         ylabelmine('$E_{\mathrm{dis}}/ M_{\mathrm{g}}\,[10^{35}\mathrm{erg\,yr^{-1}\,M_\odot^{-1}}]$')
%         
%         xlabelmine('$\log(1+z)$')
%         
%         
        
        
%         
%         %%  cooling time
%         subplot(6,4,16)
%         loglog(zpl,tcG,greyType,'color',grey)
%         hold on
%         loglog(zpl,tcI,greyType,'color',grey)
%         loglog(zpl,tcO,greyType,'color',grey)
%         
%         loglog(zpl(i:end),tcG(i:end),'-','color',colorSet(2,:));
%         loglog(zpl(i:end),tcI(i:end),'-','color',colorSet(3,:));
%         loglog(zpl(i:end),tcO(i:end),'-','color',colorSet(1,:));
%         
%         %         indx0=find(zpl2==zpl(i));
%         %         if ~isempty(indx0)
%         %             indx=indx0;
%         %         end
%         
%         h=[];
%         h(1)=loglog(zpl2(indx:end),tcG2(1,indx:end),...
%             'd','color',colorSet(2,:),'linewidth',1.5,'markersize',10,...
%             'MarkerFaceColor',colorSet(2,:),'DisplayName','Gal');
%         
%         h(2)=loglog(zpl2(indx:end),tcI2(1,indx:end),...
%             'd','color',colorSet(3,:),'linewidth',1.5,'markersize',10,...
%             'MarkerFaceColor',colorSet(3,:),'DisplayName','Inner');
%         
%         h(3)=loglog(zpl2(indx:end),tcO2(1,indx:end),...
%             'd','color',colorSet(1,:),'linewidth',1.5,'markersize',10,...
%             'MarkerFaceColor',colorSet(1,:),'DisplayName','Outer');
%         
%         hl=legend(h);
%         set(hl,'Interpreter','latex','Fontsize',12,'location','SouthWest');
%         
%         
%         grid
%         ylabelmine('Cooling Time [Gyr]')
%         xlabelmine('$\log(1+z)$')
%         
%         
%         
%         
%         
%         
%         
        %% add text
        
        rvTag='Mpc';rvfac=1;
        if log10(rvir(i))<0
            rvTag='kpc';
            rvfac=1e3;
        end
        
        
        
        str1=sprintf('$z=%4.2f$, ',zred(i));
        str2=sprintf('$\\mathrm{time}=%4.2f \\, \\mathrm{Gyr}$ \n',zredTimes(i));
        str3=sprintf('%s ,  ',mk_mvir_string(galhist.hostMass(i),'host'));
        str4=sprintf('$R_{\\mathrm{v}}=%4.2f \\,\\mathrm{%s}$ \n',rvir(i).*rvfac,rvTag);
        str5=sprintf('%s \n',mk_mvir_string(galhist.inGal.bhMassMax(i),'BH'));
        str6=''; %sprintf('$R_{\\mathrm{gal}}= %4.2f\\,\\mathrm{kpc},\\,R_{\\mathrm{1/2,g}}=%4.2f\\,\\mathrm{kpc}$',...
            %2.*galhist.halfMassRadType(5,i),galhist.halfMassRadType(1,i));
        
        str=[str1 str2 str3 str4 str5 str6];
        
        %pos=[200 500 300 300];
        pos=[0.12, 0, 0.6, 0.22];
        annot=annotation('textbox', pos, 'String', str,'Interpreter','latex','Fontsize',28,'linestyle','none');
        
        %% save Frame
        if fullMovieFlag
            
            F(cnt)=getframe(hf);
            
        end
        
        cnt=cnt+1;
        
        % printout last frame
        if i==1
            %fname=sprintf('movies/%s_%s/%s_id%i',movieFolder,simDisplayName,baseName,galID-1);
            fname=sprintf('movies/%s/%s_id%i',movieFolder,baseName,galID-1);
            
            %             fname=sprintf('movies/BCG_%s/BCG_id%i_lastFrame',simDisplayName,galID-1);
            
            printout_fig( hf,fname,'nopdf');
        end
        
        close(hf)
        
        
    end
    
    %% write movie to file
    if fullMovieFlag
        
        %name=sprintf('%s_%s/%s_id%i',movieFolder,simDisplayName,baseName,galID-1);
        name=sprintf('%s/%s_id%i',movieFolder,baseName,galID-1);
        %    name=sprintf('/BCG_TNG300/BCG_id%i',galID-1);
        
        galaxyMovieAssembly(F,name);
        
        %
        %     %% save to struct
        %
        %      BCG_movies(k).finalFrame=F(end);
        %      BCG_movies(k).id=galID-1;
        %
        %     if mod(k,5)==0
        %
        %         save([DEFAULT_MATFILE_DIR '/' fname ],'BCG_movies','done','-v7.3')
        %
        %     end
        %
        clear F
    end
end
%
%
% %% save to mat file
%
% save([DEFAULT_MATFILE_DIR '/' fname],'BCG_movies','done','-v7.3')
