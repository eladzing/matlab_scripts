%% plot method exploration compared to tng100

method_list={'0000','2201','3000','3101','3102','3103','3104',...
    '3301','3302','3403','3404','3801','3802','2101'};


load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/tng100_methodBackground.mat')
load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/methodStructure_snp4.mat') 

yl=[-3 3];
xl=[9 12.5];
snap35=4;
    


%% plot stuff
mstarLab='log Stellar Mass $[\mathrm{M_\odot}]$';
entLab='$\log K\,[\mathrm{KeV\,cm^2}]$';
cmap=brewermap(256,'Greys');
cc=brewermap(8,'Set1');


% galaxy
ii=0;
h=[];
figure

imagesc(xl,yl,log10(squeeze(galBird(:,:,1))),'alphadata',0.75);
set(gca,'Ydir','Normal','Fontsize',18)
colormap(cmap);

hold on

ii=ii+1;
xx=methodStruct.meth0000.xdataG;
yy=methodStruct.meth0000.ydataG;
h(ii)=plot(xx,yy,'x','color',cc(ii,:),'markersize',5,...
    'DisplayName','fiducial');

ii=ii+1;
xx=methodStruct.meth2201.xdataG;
yy=methodStruct.meth2201.ydataG;
h(ii)=plot(xx,yy,'o','color',cc(ii,:),'markersize',5,...
 'DisplayName','noBH');
% 
% ii=ii+1;
% xx=methodStruct.meth3000.xdataG;
% yy=methodStruct.meth3000.ydataG;
% h(ii)=plot(xx,yy,'o','color',cc(ii,:),'markersize',5,...
%  'DisplayName','noKM');

ii=ii+1;
xx=methodStruct.meth3101.xdataG;
yy=methodStruct.meth3101.ydataG;
h(ii)=plot(xx,yy,'o','color',cc(ii,:),'markersize',5,...
 'DisplayName','lowKM');

hl=legend(h);
set(hl,'Interpreter','latex','Fontsize',16,'location','SouthEast')

grid

xlabelmine(mstarLab);
ylabelmine(entLab);


% sfr=subs.SubhaloSFRinRad(tCoolStruct.galMask);  % sfr in galaxy
% ssfr=illustris.utils.calc_ssfr(subs);% sfr./galMass + 10^0.5*1e-17.*10.^(0.5.*rand(size(sfr)));
% ssfr=ssfr(tCoolStruct.galMask);
%
% ssfrThresh=1e-11;
% qMask=ssfr<=ssfrThresh;






