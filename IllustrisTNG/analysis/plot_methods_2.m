%% plot method exploration compared to tng100
if perlimFlag
    
    load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/tng100_methodBackground.mat')
    load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/methodStructure_snp4.mat')
    
    yl=[-3 3];
    xl=[9 12.5];
    snap35=4;
    
    markerList={'x','o','^','s'};
    
    mstarLab='log Stellar Mass $[\mathrm{M_\odot}]$';
    tcLab='$\log t_\mathrm{cool}\,[\mathrm{Gyr}]$';
    fgLab='$\log M_\mathrm{gas}/M_\mathrm{Baryon}$';
    ssfrLab='$\log \mathrm{sSFR}\,[\mathrm{yr^{-1}}]$';
    densLab='$\log n\,[\mathrm{cm^{-3}}]$';
    entLab='$\log S\,[\mathrm{KeV\,cm^2}]$';
    tempLab='$\log T/T_{\mathrm{vir}} $';
    
    cmapB=brewermap(256,'Greys');
    
    %cc=brewermap(8,'Set1');
    
    method_list={'0000','2201','3000',...
        '3101','3102',...
        '3103','3104',...
        '3301','3302',...
        '3403','3404',...
        '3801','3802',...
        '2101'};
    nameList={'fiducial','no BH','no LAM',...
        'low LAM','high LAM',...
        'low HAM','high HAM',...
        'high $\chi$','low $\chi$',...   
        'high $M_{\mathrm{piv}}$','low $M_{\mathrm{piv}}$',...
        'less Bursty','more Bursty',...
        'no SN'};
    
    ptagList={'fid','noBH','noKM',...
        'loKM','hiKM',...
        'loQM','hiQM',...
        'hiChi','lowChi',...
        'hiMp','loMp',...
        'loBurst','hiBurst',...
        'noSN'};
    
    printoutDir='/home/zinger/workProjects/IllustrisTNG/printout/methodsComparison';
    
    compo={'gal','cgm','out'};
    compoName={'Gal','CGM','Out'};
    
    paramType={'Ssfr','Temp','Tc','Dens','fgs'};
end

plotList={'0000','2201','3000'};
%% plot galaxy
for kk=1:length(paramType)
    skipComp=false(size(compo));
    
    if strcmp(paramType{kk},'fgs')
        skipComp(2:3)=true;
    end
    
    
    
    for jj=1:3
        
        if skipComp(jj)
            continue
        end
        
        switch paramType{kk}
            case 'Ssfr'
                param=[compo{jj} paramType{kk}];
                barLab=ssfrLab;
                cax=[-13 -9];
                cmapD=brewermap(256,'RdBu');   
            case 'fgs'
                param=paramType{kk};
                barLab=fgLab;
              
                 cax=[-3.5 0];
                 cmapD=brewermap(256,'*YlGnBu');
            case 'Tc'
                param=[compo{jj} paramType{kk}];
                barLab=tcLab;
                 cax=[-2.5 1.5];
                 cmapD=brewermap(256,'YlOrRd');
            case 'Dens'
                param=[compo{jj} paramType{kk}];
                barLab=densLab;
                 cax=[-3 -0.5];
                 cmapD=brewermap(256,'*YlGnBu');
            case 'Temp'
                param=[compo{jj} paramType{kk}];
                barLab=tempLab;
                 cax=[-1.5 1.5];
                 cmapD=brewermap(256,'*Spectral');
        end
        
        switch compo{jj}
            case 'gal'
                bird=squeeze(galBird(:,:,1));
            case 'cgm'
                bird=squeeze(cgmBird(:,:,1));
            case 'out'
                bird=squeeze(outBird(:,:,1));
        end
        
        
        
        h=[];
        figure
        set(gcf,'position',[1432 421 1000 750],'Color','w')
        
        ax0=axes('Position',[0.13 0.11 0.75 0.815]);
        xlim(xl);ylim(yl);
        
        
        imagesc(xl,yl,log10(bird),'alphadata',0.75);
        set(gca,'Ydir','Normal','Fontsize',18)
        xfac=0.05; yfac=0.9;
        text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),compoName{jj},...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
            'Interpreter','latex','fontsize',28,'fontweight','bold')
        grid
        xlabelmine(mstarLab,20);
        ylabelmine(entLab,20);
        set(gca,'fontsize',20)
        
        %hold on
        ax1=axes;
        ptag='';
        for ii=1:length(plotList)
            xx=methodStruct.(['meth' plotList{ii}]).([compo{jj} 'Mass']);
            yy=methodStruct.(['meth' plotList{ii}]).([compo{jj} 'Ent']);
            
            cc=methodStruct.(['meth' plotList{ii}]).(param);
            
            methMask=strcmp(method_list,plotList{ii});
            dispName=nameList{methMask};
            switch markerList{ii}
                case {'x','+','.','*'}
                    h(ii)=scatter(xx,yy,30,cc,markerList{ii},...
                        'DisplayName',dispName,'linewidth',2);
                case {'o','d','s','^','<','>'}
                    h(ii)=scatter(xx,yy,30,cc,markerList{ii},'filled',...
                        'DisplayName',dispName);
            end
            if ii==1; hold on ;end
            ptag=[ptag ptagList{methMask} '_'];
        end
        
        hl=legend(h);
        set(hl,'Interpreter','latex','Fontsize',16,'location','SouthEast')
        
        xlim(xl);ylim(yl);
        set(ax1,'position',get(ax0,'position'));
        
        hb=colorbar('peer',ax1,'Position',...
            [0.90 0.11 0.025 0.82],'fontsize',20);
        barTitle(hb,barLab,'fontsize',20)
        caxis(cax);
        
        linkaxes([ax0,ax1])  %linkaxes([ax0,ax1,ax2])
        ax1.Visible = 'off';ax1.XTick = [];ax1.YTick = [];
        
        set(ax0,'colormap',cmapB)
        set(ax1,'colormap',cmapD)
        
        fname=['methComp_' ptag paramType{kk} '_' compoName{jj} '_' 'snp4'];
        
        printout_fig(gcf,fname,'dir',printoutDir,'v');
        
        close(gcf);
        
    end
end

    
    
