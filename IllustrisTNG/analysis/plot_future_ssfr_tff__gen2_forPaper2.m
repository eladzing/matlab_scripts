%% plotting future sfr estimates.
snap=99;
global simDisplayName
global DRACOFLAG

if ~DRACOFLAG
    tagFont=16;
    axFont=18;
    bigTagFont=18;
    lw=3;
else
    tagFont=30;
    axFont=30;
    bigTagFont=34;
    lw=5;
end

if readFlag
    
    global DEFAULT_MATFILE_DIR
    
    load([DEFAULT_MATFILE_DIR '/futureSFR_tff_gen2_snp99_' simDisplayName '.mat'])
    
    %load([DEFAULT_MATFILE_DIR '/fofs_subs_TNG100_z0.mat'])
    
    global DRACOFLAG
    if DRACOFLAG
        fofs=illustris.groupcat.loadHalos(bp,snap);
        subs=illustris.groupcat.loadSubhalos(bp,snap);
        
    else
        loadFofSubTNG100
    end
    
    
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



sfEffList=[0.01 0.05 0.1 0.5 1.0];
sfeTag={'001' '005' '01' '05' '1'};
ttagList=[1 10];

for k=1:length(ttagList)
    ttag=ttagList(k);
    
    for j=2 % :2 %length(sfEffList)
        
        sfEff=sfEffList(j);
        
        
        
        switch(ttag)
            case 1
                fieldname='fMs';
                tfTag='$t_\mathrm{c}< t_\mathrm{ff}$';
            case 10
                fieldname='fMs10';
                tfTag='$t_\mathrm{c}< 10t_\mathrm{ff}$';
            otherwise
                error('Illegal ttag')
        end
        
        %         massG100=sfEff.*sfrStruct.inGal.(fieldname)(1,centralMask);
        %         massC100=sfEff.*sfrStruct.inCGM.(fieldname)(1,centralMask);
        %         massO100=sfEff.*sfrStruct.inOut.(fieldname)(1,centralMask);
        %         ts100=sfrStruct.ts(1).*1e9; % in yr
        %
        %         massG500=sfEff.*sfrStruct.inGal.(fieldname)(2,centralMask);
        %         massC500=sfEff.*sfrStruct.inCGM.(fieldname)(2,centralMask);
        %         massO500=sfEff.*sfrStruct.inOut.(fieldname)(2,centralMask);
        %         ts500=sfrStruct.ts(2).*1e9; % in yr
        %
        massG1=sfrStruct.inGal.(fieldname)(1,centralMask);
        massC1=sfrStruct.inCGM.(fieldname)(1,centralMask);
        massO1=sfrStruct.inOut.(fieldname)(1,centralMask);
        ts1=sfrStruct.ts(1).*1e9; % in yr
        
        %         massG16=sfEff.*sfrStruct.inGal.(fieldname)(4,centralMask);
        %         massC16=sfEff.*sfrStruct.inCGM.(fieldname)(4,centralMask);
        %         massO16=sfEff.*sfrStruct.inOut.(fieldname)(4,centralMask);
        %         ts16=sfrStruct.ts(4).*1e9; % in yr
        %
        %         massG2=sfEff.*sfrStruct.inGal.(fieldname)(5,centralMask);
        %         massC2=sfEff.*sfrStruct.inCGM.(fieldname)(5,centralMask);
        %         massO2=sfEff.*sfrStruct.inOut.(fieldname)(5,centralMask);
        %         ts2=sfrStruct.ts(5).*1e9; % in yr
        %
        %         massG5=sfEff.*sfrStruct.inGal.(fieldname)(6,centralMask);
        %         massC5=sfEff.*sfrStruct.inCGM.(fieldname)(6,centralMask);
        %         massO5=sfEff.*sfrStruct.inOut.(fieldname)(6,centralMask);
        %         ts5=sfrStruct.ts(6).*1e9; % in yr
        
        % only cgm + out gas
        %         delMass100=massC100+massO100;
        %         delMass500=massC500+massO500;
        delMass1=massC1+massO1;
        %         delMass16=massC16+massO16;
        %         delMass2=massC2+massO2;
        %         delMass5=massC5+massO5;
        
        %% plotting
        
        
        
        
        filt=fspecial('disk',12);
        xdata=log10(galMass);
        ydata0=log10(ssfr);
        
        %         xdata100=log10(galMass+delMass100);
        %         xdata500=log10(galMass+delMass500);
        xdata1=log10(galMass+delMass1);
        %         xdata16=log10(galMass+delMass16);
        %         xdata2=log10(galMass+delMass2);
        %         xdata5=log10(galMass+delMass5);
        
        
        
        
        %% plotting
        
        fprintf('plotting \n')
        
        cc=brewermap(8,'Set1');
        cc2=brewermap(8,'Set2');
        cc3=brewermap(8,'Set3');
        ccp=brewermap(6,'Paired');
        
        qmask=ssfr<1e-11;
        
        %xdata=log10(galMass);
        %         for i=3  %1:6
        %             switch i
        %                 case 1
        %                     dm=delMass100;
        %                     xdataN=xdata100;
        %                     qlim=1e-11.*ts100;
        %                     %ssfrN=ssfr100(3,:);
        %                     %ydata=ydata100;
        %                     %popC=popCont100;
        %                     %tsTag='$t_\mathrm{s}=100\,\mathrm{Myr}$';
        %                     tsTag='$100\,\mathrm{Myr}$';
        %                     nameTag='100M';
        %                 case 2
        %                     dm=delMass500;
        %                     xdataN=xdata100;
        %                     qlim=1e-11.*ts500;
        % %                     ssfrN=ssfr100(3,:);
        % %                     ydata=ydata500;
        % %                     popC=popCont500;
        %                      tsTag='$500\,\mathrm{Myr}$';
        %                     nameTag='500M';
        %                 case 3
        dm=delMass1;
        xdataN=xdata1;
        qlim=1e-11.*ts1;
        %                     ssfrN=ssfr1(3,:);
        %                     ydata=ydata1;
        %                     popC=popCont1;
        tsTag='$1\,\mathrm{Gyr}$';
        nameTag='1G';
        %                 case 4
        %                     dm=delMass16;
        %                     xdataN=xdata16;
        %                     qlim=1e-11.*ts16;
        % %                     ssfrN=ssfr16(3,:);
        % %                     ydata=ydata16;
        % %                     popC=popCont16;
        %                     tsTag='$1.6\,\mathrm{Gyr}$';
        %                     nameTag='16G';
        %                 case 5
        %                     dm=delMass2;
        %                     xdataN=xdata2;
        %                     qlim=1e-11.*ts2;
        % %                     ssfrN=ssfr2(3,:);
        % %                     ydata=ydata2;
        % %                     popC=popCont2;
        %
        %                     tsTag='$2\,\mathrm{Gyr}$';
        %                     nameTag='2G';
        %                 case 6
        %                     dm=delMass5;
        %                     xdataN=xdata5;
        %                     qlim=1e-11.*ts5;
        % %                     ssfrN=ssfr5(3,:);
        % %                     ydata=ydata5;
        % %                     popC=popCont5;
        %                     tsTag='$5\,\mathrm{Gyr}$';
        %                     nameTag='5G';
        %             end
        titTag=sprintf('$\\epsilon_\\mathrm{sf}=%3.2g$, %s, %s',sfEff,tfTag,tsTag);
        
        
        % plot mass increase
        
        
        
        
        mask=dm==0;
        add=zeros(size(dm));
        
        mbase=1e-6.*galMass;mscat=0.5;
        add(mask)=mbase(mask).*10.^(mscat.*rand(1,sum(mask)));
        yd=dm+add;
        
        
        
        % fractional mass increase
        
        yds=yd./galMass;
        %popD=plot_population_contour(xdata,log10(yds),'smooth',filt,'noplot');
        
        popQ=plot_population_contour(xdata(qmask),log10(yds(qmask)),'smooth',filt,'noplot');
        popSF=plot_population_contour(xdata(~qmask),log10(yds(~qmask)),'smooth',filt,'noplot');
        
        figure
        %set(gcf,'color','w')
        set(gcf,'position',[1432 421 1000 750],'Color','w')
        
        contour(popSF.xx,popSF.yy,popSF.popContour,'ShowText','off','LineColor',cc(2,:),...
            'LevelList',[25 50 75 95],'Fill','off','linestyle','-','linewidth',lw);
        hold on
        contour(popQ.xx,popQ.yy,popQ.popContour,'ShowText','off','LineColor',cc(1,:),...
            'LevelList',[25 50 75 95],'Fill','off','linestyle','-','linewidth',lw);
        
        %             plot(xdata(~qmask),log10(yds(~qmask)),'.','color',cc(2,:))
        %             hold on
        %             plot(xdata(qmask),log10(yds(qmask)),'.','color',cc(1,:))
        %
        %             contour(popD.xx,popD.yy,popD.popContour,'ShowText','off','LineColor',[0 0 0],...
        %                 'LevelList',20:20:100,'Fill','off','linestyle','-');
        %
        plot([9 12.5],log10([qlim qlim]),'k--','linewidth',lw);
        grid
        xl=[9 12.5];
        yl=[-6 1];
        xlim(xl)
        ylim(yl)
        
        xlabelmine('log $M_\star$ (z=0)',tagFont);
        ylabelmine('$\log\,\Delta M_\mathrm{SF}/M_\star$',tagFont);
        
        % add stuff
        qAbove=sum(yds(qmask)>qlim)/sum(qmask).*100;
        qAbove10=sum(yds(qmask & galMass>1e10)>qlim)/sum(qmask& galMass>1e10).*100;
        sfAbove=sum(yds(~qmask)>qlim)/sum(~qmask).*100;
        strQ=sprintf('Q %3.1f\\%% / %3.1f\\%%',qAbove,qAbove10);
        strS=sprintf('SF  %3.1f\\%%',sfAbove);
        
        xfac=0.6; yfac=0.95;
        text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),strQ,...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
            'Interpreter','latex','fontsize',bigTagFont,'color',cc(1,:))
        xfac=0.05; yfac=0.95;
        text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),strS,...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
            'Interpreter','latex','fontsize',bigTagFont,'color',cc(2,:))
        
        
        str1=sprintf('$\\epsilon_\\mathrm{SF}=%s\\%% $',num2str(sfEff*100));
        
        xfac=0.8; yfac=0.7;
        text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),str1,...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
            'Interpreter','latex','fontsize',bigTagFont)
        xfac=0.8; yfac=0.8;
        text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),tsTag,...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
            'Interpreter','latex','fontsize',bigTagFont)
        
        %titlemine(titTag);
        set(gca,'Fontsize',axFont)
        if(printFlag)
            name=sprintf('smassFuture2_cgm_gen2_sfe%s_tf%i_ts%s_%s',sfeTag{j},ttag,nameTag,simDisplayName);
            % printout_fig(gcf,name,'subdir','futureSsfr','v');
            printout_fig(gcf,name,'v');
        end
    end
    %fprintf('go ahead, press something...\n')
    %pause
    
    
    
    %close all
end


