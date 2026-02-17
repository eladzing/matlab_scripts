%% make a series of plots for an example subhalo 


circ(1).radius=rhalfStar;
circ(1).type='-';

circ(2).radius=rhalfGas;
circ(2).type='--';

circ(3).radius=r200c(id+1);
circ(3).type=':';

cl=[1e-7 1e-1];
aa=fieldnames(massHistStruct);
nam=string(aa{3}).extractAfter('_');
outdir="C:\Users\eladz\Documents\workProjects\IllustrisTNG\printout\forYakov\";

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
saveas(gcf,outdir+"sub_"+ nam + "_stars.png")




%% plot density maps 
hf=myFigure('pos',[388  69 1510   1214]);
tt=tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

nexttile(1)
illustris.plots.mkmapGas('gas',gas,'type','n','xz','ng',500,'fig',hf,'axes',gca,'circ',circ,'clims',cl);

nexttile(3)   
illustris.plots.mkmapGas('gas',gas,'type','n','xy','ng',500,'fig',hf,'axes',gca,'circ',circ,'clims',cl);

nexttile(4)   
illustris.plots.mkmapGas('gas',gas,'type','n','yz','ng',500,'fig',hf,'axes',gca,'circ',circ,'clims',cl);

saveas(gcf,outdir+"sub_"+ nam + "_map1.png")
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

saveas(gcf,outdir+"sub_"+ nam + "_map2.png")

%% only plot sf gas

hf=myFigure('pos',[388  69 1510   1214]);
tt=tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

nexttile(1)
illustris.plots.mkmapGas('gas',gas,'type','n','xz','ng',500,'fig',hf,'axes',gca,'mask',sfMask,'circ',circ(1));

nexttile(3)   
illustris.plots.mkmapGas('gas',gas,'type','n','xy','ng',500,'fig',hf,'axes',gca,'mask',sfMask,'circ',circ(1));

nexttile(4)   
illustris.plots.mkmapGas('gas',gas,'type','n','yz','ng',500,'fig',hf,'axes',gca,'mask',sfMask,'circ',circ(1));

saveas(gcf,outdir+"sub_"+ nam + "_map3.png")
%%
colors=brewermap(9,'Set1');
%%

for k=1:2

hf=myFigure('pos',[388  69 1510   1214]);
tt=tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

compNames=["Gal", "CGMin", "CGMout", "CGMall" "Sub"];

nexttile
h=[];
for j=1:4
    histo=massHistStruct.SubHalo_605482.(compNames(j));
    h(end+1)=plot(histo.DensityMassHist.xAxis,...
        histo.DensityMassHist.hist,...
        'color',colors(j,:),'linewidth',1.8,...
        'DisplayName',compNames(j));
    if j==1;hold on; end

end
legend(h,'Interpreter','latex')
if k==2; set(gca,'yscale','log'); end
xlabelmine('log number density')
myAxis;

nexttile
h=[];
for j=1:4
    histo=massHistStruct.SubHalo_605482.(compNames(j));
    h(end+1)=plot(histo.TempMassHist.xAxis,...
        histo.TempMassHist.hist,...
        'color',colors(j,:),'linewidth',1.8,...
        'DisplayName',compNames(j));
    if j==1;hold on; end

end
%legend(h,'Interpreter','latex')
if k==2; set(gca,'yscale','log'); end
xlabelmine('log Temperature density')
myAxis;

nexttile
h=[];
for j=1:4
    histo=massHistStruct.SubHalo_605482.(compNames(j));
    h(end+1)=plot(histo.EntropyMassHist.xAxis,...
        histo.EntropyMassHist.hist,...
        'color',colors(j,:),'linewidth',1.8,...
        'DisplayName',compNames(j));
    if j==1;hold on; end

end
%legend(h,'Interpreter','latex')
if k==2; set(gca,'yscale','log'); end
xlabelmine('log entropy')
myAxis;


nexttile
h=[];
for j=1:4
    histo=massHistStruct.SubHalo_605482.(compNames(j));
    h(end+1)=plot(histo.TcoolMassHist.xAxis,...
        histo.TcoolMassHist.hist,...
        'color',colors(j,:),'linewidth',1.8,...
        'DisplayName',compNames(j));
    if j==1;hold on; end

end
%legend(h,'Interpreter','latex')
if k==2; set(gca,'yscale','log'); end
xlabelmine('log cooling time')
myAxis;

fname=outdir+"sub_"+ nam + "_hist";
if k==2; fname=fname+ "_log"; end
saveas(gcf,fname + ".png")

end