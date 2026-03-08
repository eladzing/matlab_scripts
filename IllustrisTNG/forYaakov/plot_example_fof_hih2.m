%% plot HI H2 example 
id=39;
indx=id+1;
m200=fofs.Group_M_Crit200(indx).*illUnits.massUnit;

rv=fofs.Group_R_Crit200(indx);

circ(1).radius=rv;
circ(1).type=':';

outdir="C:\Users\eladz\Documents\workProjects\IllustrisTNG\printout\forYakov\";

centre=fofs.GroupPos(:,indx);

%% star maps 

starMask= stars.GFM_StellarFormationTime>=0; % remove wind particles
stars.newCoord = illustris.utils.centerObject(stars.Coordinates,centre);
cmapStar=brewermap(256,'*Greys');


gas.newCoord = illustris.utils.centerObject(gas.Coordinates,centre);

newStarThresh=0.1;
newStarTag='100\,\mathrm{Myr}';
newStarColor=[0.09,0.63,0.86];

%% plot star maps 
hf=myFigure('pos',[388  69 1510   1214],'color','k');
tt=tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

bsize=50; %2*ceil(rhalfGas./10)*10;

nexttile(1)
illustris.plots.mkmapStars('star',stars,'type','mass','xz','mask',starMask,...
    'ng',256,'fig',hf,'axes',gca,'circ',circ);%,...%'clims',cl,...
    %'newStars',newStarThresh,'newStarscolor',newStarColor);%,'boxsize',bsize);

nexttile(3)   
illustris.plots.mkmapStars('star',stars,'type','mass','xy','mask',starMask,...
    'ng',256,'fig',hf,'axes',gca,'circ',circ);%,...%'clims',cl,...
    %'newStars',newStarThresh,'newStarscolor',newStarColor);%,'boxsize',bsize);


nexttile(4)   
illustris.plots.mkmapStars('star',stars,'type','mass','yz','mask',starMask,...
    'ng',256,'fig',hf,'axes',gca,'circ',circ);%,...%'clims',cl,...
    %'newStars',newStarThresh,'newStarscolor',newStarColor);%,'boxsize',bsize);
set(gcf, 'InvertHardcopy', 'off')
exportgraphics(gcf,outdir+"fof_"+ num2str(id) + "_stars.png","BackgroundColor","black")


%% plot density maps 
hf=myFigure('pos',[388  69 1510   1214]);
tt=tiledlayout(2,2,'TileSpacing','tight','Padding','tight');
cl=[1e-7 1e-1];

nexttile(1)
illustris.plots.mkmapGas('gas',gas,'type','n','xz','ng',500,'fig',hf,'axes',gca,'circ',circ,'clims',cl);

nexttile(3)   
illustris.plots.mkmapGas('gas',gas,'type','n','xy','ng',500,'fig',hf,'axes',gca,'circ',circ,'clims',cl);

nexttile(4)   
illustris.plots.mkmapGas('gas',gas,'type','n','yz','ng',500,'fig',hf,'axes',gca,'circ',circ,'clims',cl);

exportgraphics(gcf,outdir+"fof_"+ num2str(id) + "_map1.png")


%% plot HI maps 
hf=myFigure('pos',[388  69 1510   1214]);
tt=tiledlayout(2,2,'TileSpacing','tight','Padding','tight');
cl=[1e-7 1e-1];

nexttile(1)
illustris.plots.mkmapGas('gas',gas,'type','hi','BR','xz','ng',500,'fig',hf,'axes',gca,'circ',circ);%,'clims',cl);

nexttile(3)   
illustris.plots.mkmapGas('gas',gas,'type','hi','BR','xy','ng',500,'fig',hf,'axes',gca,'circ',circ);%,'clims',cl);

nexttile(4)   
illustris.plots.mkmapGas('gas',gas,'type','hi','BR','yz','ng',500,'fig',hf,'axes',gca,'circ',circ);%,'clims',cl);

exportgraphics(gcf,outdir+"fof_"+ num2str(id) + "_map_hi.png")

%% plot HI maps 
hf=myFigure('pos',[388  69 1510   1214]);
tt=tiledlayout(2,2,'TileSpacing','tight','Padding','tight');
cl=[1e14 1e21];

nexttile(1)
illustris.plots.mkmapGas('gas',gas,'type','hicol','BR','xz','ng',500,'fig',hf,'axes',gca,'circ',circ,'clims',cl);

nexttile(3)   
illustris.plots.mkmapGas('gas',gas,'type','hicol','BR','xy','ng',500,'fig',hf,'axes',gca,'circ',circ,'clims',cl);

nexttile(4)   
illustris.plots.mkmapGas('gas',gas,'type','hicol','BR','yz','ng',500,'fig',hf,'axes',gca,'circ',circ,'clims',cl);

exportgraphics(gcf,outdir+"fof_"+ num2str(id) + "_map_hicol.png")
