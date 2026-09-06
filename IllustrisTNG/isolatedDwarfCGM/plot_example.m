
%% make a series of plots for an example subhalo


circ(1).radius=massHistStruct.SubHalo_605482.rgal./illUnits.lengthUnit;
circ(1).type='-';

circ(2).radius=massHistStruct.SubHalo_605482.rgas./illUnits.lengthUnit;
circ(2).type='--';

circ(3).radius=massHistStruct.SubHalo_605482.R200c./illUnits.lengthUnit;
circ(3).type=':';

cl=[1e-7 1e-1];
aa=fieldnames(massHistStruct);
nam=string(aa{3}).extractAfter('_');
outdir="C:\Users\eladz\Documents\workProjects\IllustrisTNG\printout\forYakov\";


rhog200=rho_crit(illUnits.zred,'cosmo',cosmoStruct).*illUnits.physUnits.Ms./(illUnits.physUnits.Mpc)^3./(cosmoStruct.muMass.*illUnits.physUnits.mp)*200.*cosmoStruct.Omb;

%%

%stars=illustris.snapshot.loadSubset(bp,snap,illustris.partTypeNum('stars'),{'Coordinates','Masses','GFM_StellarFormationTime'});
stars=illustris.snapshot.loadSubhalo(bp, snap, id, 'stars');
starMask= stars.GFM_StellarFormationTime>=0; % remove wind particles
stars.newCoord = illustris.utils.centerObject(stars.Coordinates,subs.SubhaloPos(:,id+1));
cmapStar=brewermap(256,'*Greys');

newStarThresh=0.1;
newStarTag='100\,\mathrm{Myr}';
newStarColor=[0.09,0.63,0.86];
%% plot star maps
hf=myFigure('pos',[388  69 1510   1214],'color','k');
tt=tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

bsize=50; %2*ceil(rhalfGas./10)*10;

nexttile(1)
illustris.plots.mkmapStars('star',stars,'type','mass','xz','mask',starMask,...
    'ng',256,'fig',hf,'axes',gca,'circ',circ,...%'clims',cl,...
    'newStars',newStarThresh,'newStarscolor',newStarColor,'boxsize',bsize);

nexttile(3)
illustris.plots.mkmapStars('star',stars,'type','mass','xy','mask',starMask,...
    'ng',256,'fig',hf,'axes',gca,'circ',circ,...%'clims',cl,...
    'newStars',newStarThresh,'newStarscolor',newStarColor,'boxsize',bsize);


nexttile(4)
illustris.plots.mkmapStars('star',stars,'type','mass','yz','mask',starMask,...
    'ng',256,'fig',hf,'axes',gca,'circ',circ,...%'clims',cl,...
    'newStars',newStarThresh,'newStarscolor',newStarColor,'boxsize',bsize);
set(gcf, 'InvertHardcopy', 'off')
exportgraphics(gcf,outdir+"sub_"+ nam + "_stars.png","BackgroundColor","black")




%% plot density maps
hf=myFigure('pos',[388  69 1510   1214]);
tt=tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

nexttile(1)
fs=illustris.plots.mkmapGas('gas',gas,'type','n','xz','ng',500,'fig',hf,'axes',gca,'circ',circ,'clims',cl,'output');

nexttile(3)
illustris.plots.mkmapGas('gas',gas,'fullcube',fs.cubeStruct,'type','n','xy','ng',500,'fig',hf,'axes',gca,'circ',circ,'clims',cl);

nexttile(4)
illustris.plots.mkmapGas('gas',gas,'fullcube',fs.cubeStruct,'type','n','yz','ng',500,'fig',hf,'axes',gca,'circ',circ,'clims',cl);

%exportgraphics(gcf,outdir+"sub_"+ nam + "_map1.png")
%% plot density maps
hf=myFigure('pos',[388  69 1510   1214]);
tt=tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

bsize=2*ceil(rr./100)*100;

nexttile(1)
illustris.plots.mkmapGas('gas',gas,'type','n','xz','ng',500,'fig',hf,'axes',gca,'circ',circ,'boxsize',bsize,'clims',cl);

nexttile(3)
illustris.plots.mkmapGas('gas',gas,'type','n','xy','ng',500,'fig',hf,'axes',gca,'circ',circ,'boxsize',bsize,'clims',cl);

nexttile(4)
illustris.plots.mkmapGas('gas',gas,'type','n','yz','ng',500,'fig',hf,'axes',gca,'circ',circ,'boxsize',bsize,'clims',cl);

exportgraphics(gcf,outdir+"sub_"+ nam + "_map2.png")

%% only plot sf gas

hf=myFigure('pos',[388  69 1510   1214]);
tt=tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

nexttile(1)
illustris.plots.mkmapGas('gas',gas,'type','n','xz','ng',500,'fig',hf,'axes',gca,'mask',sfMask,'circ',circ(1));

nexttile(3)
illustris.plots.mkmapGas('gas',gas,'type','n','xy','ng',500,'fig',hf,'axes',gca,'mask',sfMask,'circ',circ(1));

nexttile(4)
illustris.plots.mkmapGas('gas',gas,'type','n','yz','ng',500,'fig',hf,'axes',gca,'mask',sfMask,'circ',circ(1));

exportgraphics(gcf,outdir+"sub_"+ nam + "_map3.png")
%%
colors=brewermap(9,'Set1');
%%

for k=1:2

    hf=myFigure('pos',[388  69 1510   1214]);
    tt=tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

    compNames=["Gal", "CGMin", "CGMout", "CGM", "CGMoutskirt", "Sub"];
    nexttile
    h=[];
    for j=1:5
        histo=massHistStruct.SubHalo_605482.(compNames(j));
        h(end+1)=plot(histo.DensityMassHist.xAxis,...
            histo.DensityMassHist.hist,...
            'color',colors(j,:),'linewidth',1.8,...
            'DisplayName',compNames(j));
        if j==1;hold on; end
        ym(j)=min(histo.DensityMassHist.hist(histo.DensityMassHist.hist>0));
    end
    yl=get(gca,'Ylim');
    if k==2; yl(1)=min(ym);end
    ylim(yl)
    plot(log10(rhog200).*[1 1],yl,'--k');
    legend(h,'Interpreter','latex')
    if k==2; set(gca,'yscale','log'); end
    xlabelmine('log number density')
    myAxis;

    nexttile
    h=[];
    for j=1:5
        histo=massHistStruct.SubHalo_605482.(compNames(j));
        h(end+1)=plot(histo.TempMassHist.xAxis,...
            histo.TempMassHist.hist,...
            'color',colors(j,:),'linewidth',1.8,...
            'DisplayName',compNames(j));
        if j==1;hold on; end
        ym(j)=min(histo.TempMassHist.hist(histo.TempMassHist.hist>0));
    end
    yl=get(gca,'Ylim');
    if k==2; yl(1)=min(ym);end
    ylim(yl)
    plot(log10(massHistStruct.SubHalo_605482.T200c).*[1 1],yl,'--k');
    %legend(h,'Interpreter','latex')
    if k==2; set(gca,'yscale','log'); end

    xlabelmine('log Temperature density')
    myAxis;

    nexttile
    h=[];
    for j=1:5
        histo=massHistStruct.SubHalo_605482.(compNames(j));
        h(end+1)=plot(histo.EntropyMassHist.xAxis,...
            histo.EntropyMassHist.hist,...
            'color',colors(j,:),'linewidth',1.8,...
            'DisplayName',compNames(j));
        if j==1;hold on; end
        ym(j)=min(histo.EntropyMassHist.hist(histo.EntropyMassHist.hist>0));
    end
    yl=get(gca,'Ylim');
    if k==2; yl(1)=min(ym);end
    ylim(yl)
    plot(log10(massHistStruct.SubHalo_605482.K200c).*[1 1],yl,'--k');
    %legend(h,'Interpreter','latex')
    if k==2; set(gca,'yscale','log'); end
    xlabelmine('log entropy')
    myAxis;


    nexttile
    h=[];
    for j=1:5
        histo=massHistStruct.SubHalo_605482.(compNames(j));
        h(end+1)=plot(histo.TcoolMassHist.xAxis,...
            histo.TcoolMassHist.hist,...
            'color',colors(j,:),'linewidth',1.8,...
            'DisplayName',compNames(j));
        if j==1;hold on; end
        ym(j)=min(histo.TcoolMassHist.hist(histo.TcoolMassHist.hist>0));
    end
    yl=get(gca,'Ylim');
    if k==2; yl(1)=min(ym);end
    ylim(yl)
    %legend(h,'Interpreter','latex')
    if k==2; set(gca,'yscale','log'); end
    xlabelmine('log cooling time')
    myAxis;

    fname=outdir+"sub_"+ nam + "_hist";
    if k==2; fname=fname+ "_log"; end
    exportgraphics(gcf,fname + ".png")

end

%% phase diagrams

% xl=massHistStruct.SubHalo_605482.phaseDiagram.xlim;
% yl=massHistStruct.SubHalo_605482.phaseDiagram.ylim;

hf=myFigure('pos',[ 89         112        1669         982]);

tt=tiledlayout(2,3,'Padding','tight','TileSpacing','tight');
tag=["$r<r_\mathrm{gal}$","$r_\mathrm{gal}<r<r_\mathrm{gas}$",...
    "$r_\mathrm{gas}<r<r_\mathrm{200,c}$",...
    "$r_\mathrm{gal}<r<r_\mathrm{200,c}$",...
    "$r>r_\mathrm{200,c}$","Everything!"];
for k=1:length(compNames)

    nexttile;


    xl=massHistStruct.SubHalo_605482.(compNames(k)).phaseDiagram.xlim;
    yl=massHistStruct.SubHalo_605482.(compNames(k)).phaseDiagram.ylim;


    str=massHistStruct.SubHalo_605482.(compNames(k)).phaseDiagram;

    illustris.plots.plot_phaseDiagram(yl,xl,str.bird(:,:,1),...
        'fig',hf,'axes',gca,'nolabs');
    hold on;
    plot(xl,log10(massHistStruct.SubHalo_605482.T200c).*[1 1],'--k')
    plot(log10(rhog200).*[1 1],yl,'--k');
    text(-5.1,6.3,tag(k),'Interpreter','latex','FontSize',16);
    text(xl(1)+0.1,log10(massHistStruct.SubHalo_605482.T200c)+0.1,'$T_\mathrm{200,c}$','Interpreter','latex','FontSize',16);
    text(log10(rhog200)-0.5,yl(1)+0.1,'$200\,\rho_\mathrm{crit}\Omega_\mathrm{b}$','Interpreter','latex','FontSize',16);

    if k==length(compNames)
        nRef=rhog200.*[0.1 1 2];
        tRef=massHistStruct.SubHalo_605482.T200c;
        entConst=illustris.utils.calcEntropy(tRef,nRef(2)./illUnits.numberDensityFactor);
        
        tmpRef=(entConst.*illUnits.physUnits.ev.*1000).*nRef.^(2/3)./illUnits.physUnits.kb;
        plot(log10(nRef),log10(tmpRef),'--k','LineWidth',1.5)

        nRef=rhog200.*[1 15];
        tRef=massHistStruct.SubHalo_605482.T200c;
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


hf=myFigure('pos',[317          64        1173        1032]);

tt=tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
pars=["radTemp" "radDensity" "radEntropy" "radTcool"];


rgal=massHistStruct.SubHalo_605482.rgal;
rgas=massHistStruct.SubHalo_605482.rgas;
r200=massHistStruct.SubHalo_605482.r200c;
for k=1:length(pars)

    nexttile;

    str=massHistStruct.SubHalo_605482.(pars(k));

    xl=str.xlim;
    yl=str.ylim;

    yl1=yl(1);
    yl2=yl(2);

    switch pars(k)
        case "radTemp"
            tag="log temperature [K]";
            virLine=massHistStruct.SubHalo_605482.T200c;
            virTag='$T_\mathrm{200,c}$';
        case "radDensity"
            tag="log number density $[\mathrm{cm^{-3}}]$";
            virLine=rhog200;
            yl2=-0.5;
            virTag='$200\,\rho_\mathrm{crit}\Omega_\mathrm{b}$';
        case "radEntropy"
            tag="log Entropy $[\mathrm{keV\, cm^2}]$";
            virLine=massHistStruct.SubHalo_605482.K200c;
            virTag='$K_\mathrm{200,c}$';
        case"radTcool"
            tag="$\log t_\mathrm{cool}\, [\mathrm{Gyr}]$";
            virLine=cosmoStruct.tHubble;
            virTag='$t_\mathrm{hubble}$';
            yl1=-4.5;yl2=3.5;
    end


    illustris.plots.plot_param_vs_radius_2Dhistogram(xl,yl,str.bird(:,:,1),...
        'fig',hf,'axes',gca,...
        'xlab','none', ...
        'ylab',tag);
    hold on

    plot(log10(rgal./rv).*[1 1],yl,'-.k','LineWidth',1.5);
    plot(log10(rgas./rv).*[1 1],yl,'--k','LineWidth',1.5);
    plot(log10(rv./rv).*[1 1],yl,':k','LineWidth',1.5);
    plot(xl,log10(virLine).*[1 1],'--k','LineWidth',1.5);
    
    dd=(yl2-yl1)./100*3;
    text(-2+0.1,log10(virLine)+dd,virTag,'Interpreter','latex','FontSize',16);
    
    
    xlim([-2 xl(2)])
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





