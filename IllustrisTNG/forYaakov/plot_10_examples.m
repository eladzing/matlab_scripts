
%% make a series of plots for an example subhalo
global illUnits
global cosmoStruct
subList=fieldnames(massHistStruct);

cmapStar=brewermap(256,'*Greys');

newStarThresh=0.1;
newStarTag='100\,\mathrm{Myr}';
newStarColor=[0.09,0.63,0.86];

rhog200=rho_crit(illUnits.zred,'cosmo',cosmoStruct).*illUnits.physUnits.Ms./(illUnits.physUnits.Mpc)^3./(cosmoStruct.muMass.*illUnits.physUnits.mp)*200.*cosmoStruct.Omb;

outdirBase="C:\Users\eladz\Documents\workProjects\IllustrisTNG\printout\forYakov\";
cl=[1e-7 1e-1];
cl2=[1e4 1e6];
%% run over subhalos
for i=2:length(subList)

    disp(i)

    sbh=massHistStruct.(subList{i});



    circ(1).radius=sbh.rgal./illUnits.lengthUnit;
    circ(1).type='-';

    circ(2).radius=sbh.rgas./illUnits.lengthUnit;
    circ(2).type='--';

    circ(3).radius=sbh.R200c./illUnits.lengthUnit;
    circ(3).type=':';



    nam=string(subList{i});%.extractAfter('_');
    outdir=outdirBase+nam +"\";

    dstat=exist(outdir,'dir');
    if dstat==0
        mkdir(outdir);
    elseif dstat~=7
        error('something weird with output directory: %s \n',outdir);
    end


    %%

    %stars=illustris.snapshot.loadSubset(bp,snap,illustris.partTypeNum('stars'),{'Coordinates','Masses','GFM_StellarFormationTime'});

    starMask=sbh.stars.GFM_StellarFormationTime>=0; % remove wind particles

    %% plot star maps
    hf=myFigure('pos',[388  69 1510   1214],'color','k');
    tt=tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

    bsize=2*ceil(1.1*sbh.rgas./illUnits.lengthUnit);

    nexttile(1)
    fs=illustris.plots.mkmapStars('star',sbh.stars,'type','mass','xz','mask',starMask,...
        'ng',256,'fig',hf,'axes',gca,'circ',circ,...%'clims',cl,...
        'newStars',newStarThresh,'newStarscolor',newStarColor,'boxsize',bsize,'output');

    nexttile(3)
    illustris.plots.mkmapStars('star',sbh.stars,'type','mass','xy','mask',starMask,...
        'ng',256,'fig',hf,'axes',gca,'circ',circ,...%'clims',cl,...
        'newStars',newStarThresh,'newStarscolor',newStarColor,'boxsize',bsize,'fullcube',fs.cubeStruct);


    nexttile(4)
    illustris.plots.mkmapStars('star',sbh.stars,'type','mass','yz','mask',starMask,...
        'ng',256,'fig',hf,'axes',gca,'circ',circ,...%'clims',cl,...
        'newStars',newStarThresh,'newStarscolor',newStarColor,'boxsize',bsize,'fullcube',fs.cubeStruct);
    set(gcf, 'InvertHardcopy', 'off')
    exportgraphics(gcf,outdir+"sub_"+ nam + "_stars.png","BackgroundColor","black")


    %% plot density maps
    hf=myFigure('pos',[388  69 1510   1214]);
    tt=tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

    nexttile(1)
    fs=illustris.plots.mkmapGas('gas',sbh.gas,'type','n','xz','ng',500,'fig',hf,'axes',gca,'circ',circ,'clims',cl,'output');

    nexttile(3)
    illustris.plots.mkmapGas('gas',sbh.gas,'fullcube',fs.cubeStruct,'type','n','xy','ng',500,'fig',hf,'axes',gca,'circ',circ,'clims',cl);

    nexttile(4)
    illustris.plots.mkmapGas('gas',sbh.gas,'fullcube',fs.cubeStruct,'type','n','yz','ng',500,'fig',hf,'axes',gca,'circ',circ,'clims',cl);

    exportgraphics(gcf,outdir+ nam + "_map1.png")

      %% plot temperature maps
    hf=myFigure('pos',[388  69 1510   1214]);
    tt=tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

    nexttile(1)
    fs=illustris.plots.mkmapGas('gas',sbh.gas,'type','temperature','xz','ng',500,'fig',hf,'axes',gca,'circ',circ,'clims',cl2,'output');

    nexttile(3)
    illustris.plots.mkmapGas('gas',sbh.gas,'fullcube',fs.cubeStruct,'type','temperature','xy','ng',500,'fig',hf,'axes',gca,'circ',circ,'clims',cl2);

    nexttile(4)
    illustris.plots.mkmapGas('gas',sbh.gas,'fullcube',fs.cubeStruct,'type','temperature','yz','ng',500,'fig',hf,'axes',gca,'circ',circ,'clims',cl2);

    exportgraphics(gcf,outdir+ nam + "_map2.png")

    %% plot vr maps
    id=nam.extractAfter('_').double;
    vcm=subs.SubhaloVel(:,id+1);


    hf=myFigure('pos',[388  69 1510   1214]);
    tt=tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

    nexttile(1)
    fs=illustris.plots.mkmapGas('gas',sbh.gas,'type','vr','vcm',vcm,'xz','ng',500,'fig',hf,'axes',gca,'circ',circ,'output','clims',65*[-1 1]);

    nexttile(3)
    illustris.plots.mkmapGas('gas',sbh.gas,'fullcube',fs.cubeStruct,'type','vr','vcm',vcm,'xy','ng',500,'fig',hf,'axes',gca,'circ',circ,'clims',65*[-1 1]);

    nexttile(4)
    illustris.plots.mkmapGas('gas',sbh.gas,'fullcube',fs.cubeStruct,'type','vr','vcm',vcm,'yz','ng',500,'fig',hf,'axes',gca,'circ',circ,'clims',65*[-1 1]);

    exportgraphics(gcf,outdir+ nam + "_vrMap.png")

    %% plot density maps
    % hf=myFigure('pos',[388  69 1510   1214]);
    % tt=tiledlayout(2,2,'TileSpacing','tight','Padding','tight');
    %
    % bsize=2*ceil(rr./100)*100;
    %
    % nexttile(1)
    % illustris.plots.mkmapGas('gas',gas,'type','n','xz','ng',500,'fig',hf,'axes',gca,'circ',circ,'boxsize',bsize,'clims',cl);
    %
    % nexttile(3)
    % illustris.plots.mkmapGas('gas',gas,'type','n','xy','ng',500,'fig',hf,'axes',gca,'circ',circ,'boxsize',bsize,'clims',cl);
    %
    % nexttile(4)
    % illustris.plots.mkmapGas('gas',gas,'type','n','yz','ng',500,'fig',hf,'axes',gca,'circ',circ,'boxsize',bsize,'clims',cl);
    %
    % exportgraphics(gcf,outdir+"sub_"+ nam + "_map2.png")

    %% only plot sf gas
    sfMask=sbh.gas.StarFormationRate>0;
    hf=myFigure('pos',[388  69 1510   1214]);
    tt=tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

    nexttile(1)
    fs=illustris.plots.mkmapGas('gas',sbh.gas,'type','n','xz','ng',500,'fig',hf,'axes',gca,'mask',sfMask,'circ',circ(1),'output');

    nexttile(3)
    illustris.plots.mkmapGas('gas',sbh.gas,'type','n','xy','ng',500,'fig',hf,'axes',gca,'mask',sfMask,'circ',circ(1),'fullcube',fs.cubeStruct);

    nexttile(4)
    illustris.plots.mkmapGas('gas',sbh.gas,'type','n','yz','ng',500,'fig',hf,'axes',gca,'mask',sfMask,'circ',circ(1),'fullcube',fs.cubeStruct);

    exportgraphics(gcf,outdir+"sub_"+ nam + "_mapSF.png")
    %%
    colors=brewermap(9,'Set1');
    colors=colors([1:5 7:9],:);
    %% plot mass histograms
    compNames=["Gal", "CGMin", "CGMout", "CGM", "CGMoutskirt", "Sub"];
    pars=["Density" "Temp" "Entropy" "Tcool" "Zmet"];
    virs=[rhog200 sbh.T200c  sbh.K200c cosmoStruct.tHubble];
    for k=1:2

        hf=myFigure('pos',[388  69 2004   1214]);
        tt=tiledlayout(2,3,'TileSpacing','tight','Padding','tight');


        for ii=1:length(pars)

            nexttile
            h=[];
            for j=1:length(compNames)
                histo=sbh.(compNames(j));
                hst=histo.(pars(ii)+"MassHist");
                if j~=length(compNames)
                    h(end+1)=plot(hst.xAxis,...
                        hst.hist,...
                        'color',colors(j,:),'linewidth',1.8,...
                        'DisplayName',compNames(j));
                else
                    h(end+1)=plot(hst.xAxis,...
                        hst.hist,'k-',...
                        'linewidth',1.8,...
                        'DisplayName',compNames(j));
                end
                if j==1;hold on; end
                ym(j)=min(hst.hist(hst.hist>0));
            end
            yl=get(gca,'Ylim');
            if k==2; yl(1)=min(ym);end
            ylim(yl)
            if ii<5
                plot(log10(virs(ii)).*[1 1],yl,'--k');
            end
            legend(h,'Interpreter','latex')
            if k==2; set(gca,'yscale','log'); end
            xlabelmine("log " + pars(ii))
            myAxis;
        end




        tt.YLabel.String='Mass $[\mathrm{M_\odot}]$';
        tt.YLabel.FontSize=22;
        tt.YLabel.Interpreter='latex';



        fname=outdir+"sub_"+ nam + "_hist";
        if k==2; fname=fname+ "_log"; end
        exportgraphics(gcf,fname + ".png")

    end

    %% phase diagrams

    % xl=sbh.phaseDiagram.xlim;
    % yl=sbh.phaseDiagram.ylim;

    hf=myFigure('pos',[ 89         112        1669         982]);

    tt=tiledlayout(2,3,'Padding','compact','TileSpacing','tight');
    tag=["$r<r_\mathrm{gal}$","$r_\mathrm{gal}<r<r_\mathrm{gas}$",...
        "$r_\mathrm{gas}<r<r_\mathrm{200,c}$",...
        "$r_\mathrm{gal}<r<r_\mathrm{200,c}$",...
        "$r>r_\mathrm{200,c}$","Everything!"];

    massLim=log10([0 max(max(sbh.Sub.phaseDiagram.bird(:,:,1)))]);
    for k=1:length(compNames)

        nexttile;


        xl=sbh.(compNames(k)).phaseDiagram.xlim;
        yl=sbh.(compNames(k)).phaseDiagram.ylim;


        str=sbh.(compNames(k)).phaseDiagram;

        illustris.plots.plot_phaseDiagram(yl,xl,str.bird(:,:,1),...%'clim',massLim,...
            'fig',hf,'axes',gca,'nolabs');
        hold on;
        plot(xl,log10(sbh.T200c).*[1 1],'--k')
        plot(log10(rhog200).*[1 1],yl,'--k');

        text(xl(1)+0.1*diff(xl),yl(1)+0.95*diff(yl),tag(k),'Interpreter','latex','FontSize',16);
        text(xl(1)+0.1,log10(sbh.T200c)+0.1,'$T_\mathrm{200,c}$','Interpreter','latex','FontSize',16);
        text(log10(rhog200)-0.5,yl(1)+0.1,'$200\,\rho_\mathrm{crit}\Omega_\mathrm{b}$','Interpreter','latex','FontSize',16);

        if k==length(compNames)
            nRef=rhog200.*[0.1 1 2];
            tRef=sbh.T200c;
            entConst=illustris.utils.calcEntropy(tRef,nRef(2)./illUnits.numberDensityFactor);

            tmpRef=(entConst.*illUnits.physUnits.ev.*1000).*nRef.^(2/3)./illUnits.physUnits.kb;
            plot(log10(nRef),log10(tmpRef),'--k','LineWidth',1.5)

            nRef=rhog200.*[1 15];
            tRef=sbh.T200c;
            prConst=illustris.utils.calcPressure(nRef(1),tRef);

            tmpRef=(prConst).*nRef.^(-1)./illUnits.physUnits.kb;
            plot(log10(nRef),log10(tmpRef),':k','LineWidth',1.5)

        end

        titlemine(compNames(k));
    end

    tt.YLabel.String='$\log(T)\,[\mathrm{K}]$';
    tt.YLabel.FontSize=22;
    tt.YLabel.Interpreter='latex';

    tt.XLabel.String='$\log(n)\,[\mathrm{cm^{-3}}]$';
    tt.XLabel.FontSize=22;
    tt.XLabel.Interpreter='latex';

    fname=outdir+"sub_"+ nam + "_birds";
    exportgraphics(gcf,fname+".png");


    %% rad birds


    hf=myFigure('pos',[317          64        2000        1032]);

    tt=tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
    pars=["radTemp" "radDensity" "radEntropy" "radTcool" "radZmet"];


    % rgal=sbh.rgal;
    % rgas=sbh.rgas;
    % r200=sbh.r200c;
    xm=log10(0.5.*sbh.rgal./sbh.R200c);
    for k=1:length(pars)

        nexttile;

        str=sbh.(pars(k));

        xl=str.xlim;
        yl=str.ylim;

        yl1=yl(1);
        yl2=yl(2);

        switch pars(k)
            case "radTemp"
                tag="log temperature [K]";
                virLine=sbh.T200c;
                virTag='$T_\mathrm{200,c}$';
            case "radDensity"
                tag="log number density $[\mathrm{cm^{-3}}]$";
                virLine=rhog200;
                yl2=-0.5;
                virTag='$200\,\rho_\mathrm{crit}\Omega_\mathrm{b}$';
            case "radEntropy"
                tag="log Entropy $[\mathrm{keV\, cm^2}]$";
                virLine=sbh.K200c;
                virTag='$K_\mathrm{200,c}$';
            case"radTcool"
                tag="$\log t_\mathrm{cool}\, [\mathrm{Gyr}]$";
                virLine=cosmoStruct.tHubble;
                virTag='$t_\mathrm{hubble}$';
                yl1=-4.5;yl2=3.5;

            case "radZmet"
                tag="log Metallicity $[\mathrm{Z_\odot}]$";
                %virLine=sbh.K200c;
                %virTag='$K_\mathrm{200,c}$';
        end


        illustris.plots.plot_param_vs_radius_2Dhistogram(xl,yl,str.bird(:,:,1),...
            'fig',hf,'axes',gca,...
            'xlab','none', ...
            'ylab',tag);
        hold on

        plot(log10(sbh.rgal./sbh.R200c).*[1 1],yl,'-.k','LineWidth',1.5);
        plot(log10(sbh.rgas./sbh.R200c).*[1 1],yl,'--k','LineWidth',1.5);
        plot(log10(sbh.R200c./sbh.R200c).*[1 1],yl,':k','LineWidth',1.5);
        if k<5
            plot(xl,log10(virLine).*[1 1],'--k','LineWidth',1.5);

            dd=(yl2-yl1)./100*3;
            text(xm+0.08,log10(virLine)+dd,virTag,'Interpreter','latex','FontSize',16);
        end

        xlim([xm  0.5]); %xl(2)])
        ylim([yl1 yl2])
        if k==1
            legend({'$r_\mathrm{gal}$','$r_\mathrm{gas}$','$r_\mathrm{200,c}$'},...
                'Interpreter','latex','FontSize',16)
        end



        %titlemine(compNames(k));
    end

    tt.XLabel.String='$\log\, r/\mathrm{R_{200,c}}$';
    tt.XLabel.FontSize=22;
    tt.XLabel.Interpreter='latex';

    % tt.YLabel.String='$\log(n)\,[\mathrm{cm^{-3}}]$';
    % tt.YLabel.FontSize=22;
    % tt.YLabel.Interpreter='latex';
    %
    fname=outdir+"sub_"+ nam + "_radBirds";
    exportgraphics(gcf,fname+".png");

%%

    close('all')
end


%% get limits for parameters
for i=2:length(subList)
    sbh=massHistStruct.(subList{i});
    compNames=["Gal", "CGMin", "CGMout", "CGM", "CGMoutskirt", "Sub"];
    for j=1:length(compNames)
        lims.tempLim.(compNames(j))(:,i-1)=sbh.(compNames(j)).TempMassHist.xAxis([1 end]);
        lims.tempLim.(compNames(j)+"_pd")(:,i-1)=sbh.(compNames(j)).phaseDiagram.ylim;

        lims.densLim.(compNames(j))(:,i-1)=sbh.(compNames(j)).DensityMassHist.xAxis([1 end]);
        lims.densLim.(compNames(j)+"_pd")(:,i-1)=sbh.(compNames(j)).phaseDiagram.xlim;

        lims.entLim.(compNames(j))(:,i-1)=sbh.(compNames(j)).EntropyMassHist.xAxis([1 end]);
        lims.tcLim.(compNames(j))(:,i-1)=sbh.(compNames(j)).TcoolMassHist.xAxis([1 end]);

        lims.zmetLim.(compNames(j))(:,i-1)=sbh.(compNames(j)).ZmetMassHist.xAxis([1 end]);
    end
end


%% plot limits.
cols=inferno(12);
dd=linspace(-0.2,0.2,10);

comp=["Gal" "CGMin" "CGMout" "CGM" "CGMoutskirt" "Sub"];

hf=myFigure('pos',[317          64        2000        1032]);

tt=tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

flds=fieldnames(lims);
for k=1:length(flds)
    nexttile

    for ii=1:10
        for j=1:length(comp)
            plot([j j]+dd(ii),lims.(flds{k}).(comp(j))(:,ii),'color',cols(ii,:))
            if ii==1;hold on;end
        end
    end
    xticks(1:6)
    xticklabels(comp)
    xlim([0.6 6.4])
    ylabelmine(flds{k});
    myAxis;

end

fname=outdirBase+"dwarf10_par_limits";
exportgraphics(gcf,fname+".png");


%% plot limits - t & rho comparison to phase diagram
cols=inferno(12);
dd=linspace(-0.2,0.2,10);

comp=["Gal" "CGMin" "CGMout" "CGM" "CGMoutskirt" "Sub"];

hf=myFigure('pos',[1406          73        1405        1220]);

tt=tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

flds=fieldnames(lims);

for k=1:4
    nexttile

    k1=round((k)/2);
    k2=2-mod(k,2);
    for ii=1:10
        for j=1:length(comp)
            compName=comp(j);
            if k2==2
                compName=compName+"_pd";
            end
            plot([j j]+dd(ii),lims.(flds{k1}).(compName)(:,ii),'color',cols(ii,:))
            if ii==1;hold on;end
        end
    end
    xticks(1:6)
    xticklabels(comp)
    xlim([0.6 6.4])
    ylabelmine(flds{k1});
    myAxis;

end

fname=outdirBase+"dwarf10_par_limits_pd";
exportgraphics(gcf,fname+".png");