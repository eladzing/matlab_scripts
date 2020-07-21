% generate histories of centrals.

%global matFilePath
%global cosmoStruct
if ~skip2point
    global DRACOFLAG
    global simDisplayName
    %sim='100';
    snap=99; %z=0
    zmax=3.0;
    
    bp=illustris.set_env(sim,'draco');
    
    global illUnits
    
    if ~exist('readFlag','var')
        readFlag=true;
    end
    
    if readFlag
        fprintf(' *** Reading data *** \n');
        fofs=illustris.groupcat.loadHalos(bp,snap);
        subs=illustris.groupcat.loadSubhalos(bp,snap);
        readFlag=false;
        
        
    end
    
    
    massThresh=10^9; % threshold for *stellar* mass
    
    %% load subsinfo
    subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs,snap);
    units; % load general unit structure in cgs.
    
    
    massAllGals= subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit; % stellar mass within 2*rhalf
    massMask=massAllGals>massThresh;
    
    % this mask selects all galaxies with dm component & stars, above *stellar* mass limit whose host has virials parameters
    %galMask= subsInfo.isCentral & subsInfo.DM10perc & subsInfo.hasStars & massMask & subsInfo.hostHasVirial & subsInfo.hasGas ;%  & subs.SubhaloMassInRadType(1,:)>0;
    galMask=generateMask('subs',subs,'fofs',fofs,'central','gas','mass',massThresh);
    
    
    
    % define important parameters
    % load relevant data structure
    global DEFAULT_MATFILE_DIR
    load([DEFAULT_MATFILE_DIR '/cooling_times_z0_TNG100.mat'])
    
    
    galMass=tCoolStruct.galMass;  % galaxy stellar mass
    %gasMass=tCoolStruct.inGal.gasMass(1,:);
    gasEnt=tCoolStruct.inGal.meanEntMW(1,:);
    
    
    ssfr=illustris.utils.calc_ssfr(subs);
    
    %% select galaxies
    
    if selectFlag
        
        % create tree
        filt=fspecial('disk',8);
        xdata=log10(galMass(galMask));
        ydata=log10(gasEnt(galMask));
        %popCont=plot_population_contour(xdata,ydata,'smooth',filt,'noplot');
        
        minLev=6;
        splitParam=30;
        yl=[-3 2.5];
        xl=[9 12.5];
        cc=log10(ssfr(galMask));
        cax=[min(cc) max(cc)];
        
        cmap=brewermap(256,'Spectral');
        ssfrLab='$\log \mathrm{sSFR}_{\mathrm{Gal}},[\mathrm{yr^{-1}}] $';
        mstarLab='$\log M_\star\,[\mathrm{M_\odot}]$';
        entLab=sprintf('$\\log S_{\\mathrm{%s}} \\,[\\mathrm{KeV\\,cm^2}]$','Gal');
        
        tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
        celVal=points2tree(cc,tre,'median');
        
        plot2dTreeMap(celVal,tre,'cmap',cmap,'minmax',cax);
        hb=colorbar;barTitle(hb,ssfrLab,'fontsize',16)
        set(hb,'fontsize',14);
        colormap(cmap);caxis(cax);
        xlim(xl);ylim(yl);
        
        
        grid
        xlabelmine(mstarLab,18);
        ylabelmine(entLab,18);
        set(gca,'fontsize',14)
        
        [xPol,yPol]=ginput;
        
        hold on
        
        plot(xPol,yPol,'k-')
        %plot_massEnt('tree',xdata,ydata,log10(ssfr(galMask)),tre,ssfrLab,flipud(cmap),'inGal')
        
    else
        
        xPol=[10.0133,9.9614,10.0378,10.2059,10.3832,10.7347,10.8784,11.1321,11.0862,9.9950,9.9950];
        
        yPol=[2.1161,0.4951,0.0514,-0.0851,0.0685,0.5292,0.8250,1.2174,2.3663,2.3607,2.3607];
        % quenched polygon
        %xPol=[11.1499   10.0046   10.0046   10.1872   11.1005   11.1252];
        %yPol=[2.2345    2.2975    0.2811   -0.1600    1.0822    2.2345];
        
    end
    
    polyMask=inpolygon(log10(galMass),log10(gasEnt),xPol,yPol);
    groupMask=polyMask & galMask;
    %qGroupInd=find(groupMask);
    
end

%% generate values
fprintf(' *** Running over %s Galaxies *** \n',num2str(length(groupMask)));

step=5;
stepNext=5;
len=double(subs.count);
cnt=0;
nextP=0;

done=false(size(galMask));

for id=0:len-1
    
    perCent=floor((id+1)/len*100);
    
    if (perCent>=stepNext)     %       mod(perCent,5)==0 && perCent>=5)
        fprintf('%s %% of Galaxies done  \n',num2str(floor(perCent)))
        stepNext=stepNext+step;
    end
    
    
    if groupMask(id+1)
        cnt=cnt+1;
        
        prc=floor(cnt./sum(groupMask)*100);
        if (prc>=nextP)
            
            fprintf(' *** %i %% of histories completed  *** \n',prc);
            nextP=nextP+2;
        end
        
        [galHistory(id+1),galBirdHistory(id+1)]=get_galaxy_history_full(id,'zmax',zmax,'startSnap',snap);
        done(id+1)=true;
        
        
        if mod(cnt,50)==0 % print interim
            
            
            fname=sprintf('quenchedGroup_history_z%s_%s_id%i',num2str(illustris.utils.get_zred(snap)),simDisplayName,id);
            fnameBird=sprintf('quenchedGroup_history_birds_z%s_%s_id%i',num2str(illustris.utils.get_zred(snap)),simDisplayName,id);
            save([DEFAULT_MATFILE_DIR '/' fname],'galHistory','done','-v7.3')
            save([DEFAULT_MATFILE_DIR '/' fnameBird],'galBirdHistory','done','-v7.3')
            fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);
            fprintf(' *** Bird Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fnameBird]);
            
        end
    end
end
if ~done(len)
    ii=find(~done,1,'first');
    galHistory(len)=galHistory(ii);
    galBirdHistory(len)=galBirdHistory(ii);
end

%% organize output



if DRACOFLAG
    %global DEFAULT_MATFILE_DIR
    %global simDisplayName
    fname=sprintf('quenchedGroup_history_z%s_%s',num2str(illustris.utils.get_zred(snap)),simDisplayName);
    fnameBird=sprintf('quenchedGroup_history_birds_z%s_%s',num2str(illustris.utils.get_zred(snap)),simDisplayName);
    save([DEFAULT_MATFILE_DIR '/' fname],'galHistory','done','-v7.3')
    save([DEFAULT_MATFILE_DIR '/' fnameBird],'galBirdHistory','done','-v7.3')
    fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);
    fprintf(' *** Bird Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fnameBird]);
    
end
fprintf(' *** DONE!  *** \n');

