%% script for removing point sources from an x-ray map
% after plotting a given projection, allow a interactive choice of a region
% to extract. Then fill in the area by the average along a radial shell

%% get relevant information

%load('C:\Users\eladzing\Documents\matlab_scripts\cluster\mat_files\xrayMapsFull_2.mat')
load('/home/zinger/workProjects/matlab_scripts/cluster/mat_files/xrayMapsFull_2.mat')

% get cluster and projection

choice= chooseClusterProjection;

new_env(choice.cluster)

global ClusterIndex
global hub
boxx=8;
xMapStruct=xRayMap(ClusterIndex);

switch(lower(choice.proj))
    case 'xy'
        xMap0=xMapStruct.projXY;
        xLab='X';yLab='Y';
        printag='XY';
    case 'xz'
        xMap0=xMapStruct.projXZ;
        xLab='X';yLab='Z';
        printag='XZ';
    case 'yz'
        xMap0=xMapStruct.projYZ;
        xLab='Z';yLab='Y';
        printag='YZ';
end

if isempty(xMap0)
    error('%s does not yet have an xray map',choice.cluster)
end
r500=get_rvir(500);
r200=get_rvir(200);


cmap=brewermap(256,'*RdBu');
% set scale
xl=r500.*[-2.2 2.2].*hub;

figure;
imagesc([-4 4],[-4 4],log10(xMap0));
set(gca,'Ydir','Normal','Fontsize',14,'XLim',xl,'YLim',xl)
bar=colorbar;set(bar,'Fontsize',12)
barTitle(bar,'$\mathrm{log(counts\,arcsec^{-2}\,s^{-1})}$')
caxis([-12 -2 ])
draw_circle_boxx(gcf,r200,'white');      % draw rvir
draw_circle_boxx(gcf,r500,'white');      % draw rvir
draw_circle_boxx(gcf,2*r500,'white');      % draw rvir
xlabelmine(sprintf('$%s\\,[\\mathrm{Mpc/h}]$',xLab));
ylabelmine(sprintf('$%s\\,[\\mathrm{Mpc/h}]$',yLab));
titlemine(xMapStruct.cluster);
colormap(cmap);
axis square


%% generate position and radial grid
sz=size(xMap0);
cl=boxx/sz(1);
[meshX, meshY] = meshgrid(1:sz(1), 1:sz(2));

meshX = meshX - (sz(1)+1)/2;
meshY = meshY - (sz(2)+1)/2;
meshX = meshX * ((boxx)./sz(1));
meshY = meshY * ((boxx)./sz(2));



%% select region and threshold
extractFlag=true;
xMap=xMap0;
while(extractFlag)
    
    %warndlg('Choose area for extraction');
    
    
    %% plot map
    hf1=figure;
    imagesc([-4 4],[-4 4],log10(xMap));
    set(gca,'Ydir','Normal','Fontsize',14,'XLim',xl,'YLim',xl)
    bar=colorbar;set(bar,'Fontsize',12)
    barTitle(bar,'$\mathrm{log(counts\,arcsec^{-2}\,s^{-1})}$')
    caxis([-12 -2 ])
    draw_circle_boxx(gcf,r200,'white');      % draw rvir
    draw_circle_boxx(gcf,r500,'white');      % draw rvir
    draw_circle_boxx(gcf,2*r500,'white');      % draw rvir
    xlabelmine(sprintf('$%s\\,[\\mathrm{Mpc/h}]$',xLab));
    ylabelmine(sprintf('$%s\\,[\\mathrm{Mpc/h}]$',yLab));
    titlemine(xMapStruct.cluster);
    colormap(cmap);
    axis square
    
    [xp, yp] = getline(gcf); % uncomment this line for interactive selection
    line(xp,yp, 'Color','g');
    
    xlZoom=[min(xp) max(xp)];
    ylZoom=[min(yp) max(yp)];
    xlZoom=xlZoom+0.5.*[-1 1].*diff(xlZoom);
    ylZoom=ylZoom+0.5.*[-1 1].*diff(ylZoom);
    set(gca,'XLim',xlZoom,'YLim',ylZoom)
    
    %warndlg('Choose value for lower threshold')
    [xd,yd]=ginput(1);
    hold on
    plot(xd,yd,'w+','markersize',8)
    ix=ceil( (xd+boxx./2)./cl);
    iy=ceil( (yd+boxx./2)./cl);
    thresh=xMap(iy,ix);
    
    mask = inpolygon(meshX,meshY,xp,yp) & xMap>thresh;
    
    hf2=figure;
    imagesc([-4 4],[-4 4],log10(xMap.*~mask));
    
    set(gca,'Ydir','Normal','Fontsize',14,'XLim',xl,'YLim',xl)
    bar=colorbar;set(bar,'Fontsize',12)
    barTitle(bar,'$\mathrm{log(counts\,arcsec^{-2}\,s^{-1})}$')
    caxis([-12 -2 ])
    draw_circle_boxx(gcf,r200,'white');      % draw rvir
    draw_circle_boxx(gcf,r500,'white');      % draw rvir
    draw_circle_boxx(gcf,2*r500,'white');      % draw rvir
    xlabelmine(sprintf('$%s\\,[\\mathrm{Mpc/h}]$',xLab));
    ylabelmine(sprintf('$%s\\,[\\mathrm{Mpc/h}]$',yLab));
    titlemine(xMapStruct.cluster);
    colormap(cmap);
    axis square
    line(xp,yp, 'Color','g');
    
    
    %% find the mean along the radial shell
    radMap=sqrt(meshX.^2+meshY.^2);
    
    
    radRange=[min(radMap(mask)) max(radMap(mask))];
    
    radMask=radMap>=radRange(1) & radMap<=radRange(2);
    
    figure(hf2)
    hold on
    draw_circle_boxx(gcf,radRange(1)./hub,'red');
    draw_circle_boxx(gcf,radRange(2)./hub,'red');
    
    xMean=mean(xMap(radMask & ~mask));
    
    xMap2=xMap;
    xMap2(mask)=xMean;
    
    xMap=xMap2;
    
    button = questdlg('Extract another region?');
    
    extractFlag=strcmp(button,'Yes');
    close(hf1,hf2);
end

%% replot and save

hf3=figure;
imagesc([-4 4],[-4 4],log10(xMap2));

set(gca,'Ydir','Normal','Fontsize',14); %,'XLim',xl,'YLim',xl)
bar=colorbar;set(bar,'Fontsize',12)
barTitle(bar,'$\mathrm{log(counts\,arcsec^{-2}\,s^{-1})}$')
caxis([-12 -2 ])
draw_circle_boxx(gcf,r200,'white');      % draw rvir
draw_circle_boxx(gcf,r500,'white');      % draw rvir
draw_circle_boxx(gcf,2*r500,'white');      % draw rvir
xlabelmine(sprintf('$%s\\,[\\mathrm{Mpc/h}]$',xLab));
ylabelmine(sprintf('$%s\\,[\\mathrm{Mpc/h}]$',yLab));
titlemine(xMapStruct.cluster);
colormap(cmap);
axis square

xRayMapFilterd=xMapStruct;

switch(lower(choice.proj))
    case 'xy'
        xRayMapFilterd.projXY=xMap2;
    case 'xz'
        xRayMapFilterd.projXZ=xMap2;
    case 'yz'
        xRayMapFilterd.projYZ=xMap2;
end


%% plot and save

global NCELL
global aexpn

mapLen=8*NCELL;

cl=get_cellsize(1,'mpc');

cLen=ceil(min(2.1*r500/cl,mapLen/2));

indx=mapLen/2-cLen+1:mapLen/2+cLen;
xx=(indx-(mapLen/2+0.5)).*cl;


%     switch(lower(choice.proj))
%     case 'xy'
%         xMap0=xMapStruct.projXY;
%         xLab='X';yLab='Y';
%     case 'xz'
%         xMap0=xMapStruct.projXZ;
%         xLab='X';yLab='Z';
%     case 'yz'
%         xMap0=xMapStruct.projYZ;
%         xLab='Z';yLab='Y';
% end

proj0=xMap2(indx,indx);
%projYZ0=xRayMap(k).projYZ(indx,indx);
%projXZ0=xRayMap(k).projXZ(indx,indx);



%% degrade

degScale = 500;
[yg,xg]=meshgrid(xx,xx);

xxx=reshape(xg,[numel(xg) 1]);
yyy=reshape(yg,[numel(yg) 1]);

vvv=reshape(proj0,[numel(proj0) 1]);

[proj, binsize, xxlim,yylim]= remap_fine2coarse(xxx,yyy,vvv,cl,'len',500,'wt',cl.^2);
proj=proj(:,:,1)./prod(binsize);

xl=xxlim(1)+0.5.*binsize(1);

%% plot images
global DEFAULT_PRINTOUT_DIR
printoutdir=sprintf('%s/sz_data',DEFAULT_PRINTOUT_DIR);

map = brewermap(256,'*RdBu');
colormap(map);


figure
imagesc(xxlim*hub,yylim*hub,log10(proj));
set(gca,'Ydir','Normal','Fontsize',14)
bar=colorbar;set(bar,'Fontsize',12)
barTitle(bar,'$\mathrm{log(counts\,arcsec^{-2}\,s^{-1})}$')
caxis([-12 -2 ])
draw_circle_boxx(gcf,r200,'white');      % draw rvir
draw_circle_boxx(gcf,r500,'white');      % draw rvir
draw_circle_boxx(gcf,2*r500,'white');      % draw rvir
xlabelmine(sprintf('$%s\\,[\\mathrm{Mpc/h}]$',xLab));
ylabelmine(sprintf('$%s\\,[\\mathrm{Mpc/h}]$',yLab));
titlemine(xMapStruct.cluster);
colormap(map);
axis square
printout_fig(gcf,sprintf('%s_xrayMap%sFiltered_n%d_%s',xMapStruct.cluster,printag,degScale,aexpn),'subdir','sz_data','png')

%% remove center
str=sprintf('removal radius in pixels (%s kpc)',num2str(binsize(1).*1000));
answer=inputdlg(str);

cutRad=str2double(answer{1});
if cutRad>1
    len=length(proj);
    cen=len/2+0.5;
    [xx,yy]=meshgrid(1:500);
    rr=sqrt((xx-cen).^2+(yy-cen).^2);
    mask=rr<=cutRad;
    
    mask2=rr>=cutRad & rr<cutRad+1;
    
    proj(mask)=mean(mean(proj(mask2)));
    
    
    
    figure
    imagesc(xxlim*hub,yylim*hub,log10(proj));
    set(gca,'Ydir','Normal','Fontsize',14)
    bar=colorbar;set(bar,'Fontsize',12)
    barTitle(bar,'$\mathrm{log(counts\,arcsec^{-2}\,s^{-1})}$')
    caxis([-12 -2 ])
    draw_circle_boxx(gcf,r200,'white');      % draw rvir
    draw_circle_boxx(gcf,r500,'white');      % draw rvir
    draw_circle_boxx(gcf,2*r500,'white');      % draw rvir
    xlabelmine(sprintf('$%s\\,[\\mathrm{Mpc/h}]$',xLab));
    ylabelmine(sprintf('$%s\\,[\\mathrm{Mpc/h}]$',yLab));
    titlemine(xMapStruct.cluster);
    colormap(map);
    axis square
    printout_fig(gcf,sprintf('%s_xrayMap%sFiltered_n%d_%s',xMapStruct.cluster,printag,degScale,aexpn),'subdir','sz_data')
    
end

%% write to file


head=sprintf('## bottom/left (Mpc): %s / %s , n_elements: %i x %i , cellSize: %s kpc, units: counts/sec/arcsec^2',...
    num2str(xl),num2str(xl),degScale,degScale,num2str(1000*binsize(1)));
global HALO_PATH

fname=sprintf('%s/%s_xrayMapFiltered_%s_n%d_%s.dat',HALO_PATH,xMapStruct.cluster,printag,degScale,aexpn);
fid=fopen(fname,'w');
fprintf(fid,'%s \n',head);
fprintf(fid,'%12.5g \n',reshape(proj',[numel(proj) 1]));

fclose('all');
%close all;





% %% auxilary function fior dialog box
% function res = chooseDialog
% 
% 
% d = dialog('Position',[300 300 350 250],'Name','Select Cluster and Projection');
% txt1 = uicontrol('Parent',d,...
%     'Style','text',...
%     'Position',[20 180 210 40],...
%     'String','Select a cluster');
% 
% popup = uicontrol('Parent',d,...
%     'Style','popup',...
%     'Position',[100 170 100 25],...
%     'String',{'CL101';'CL102';'CL103';...
%     'CL104';'CL105';'CL106';'CL107';...
%     'CL3';'CL5';'CL6';'CL7';'CL9';...
%     'CL10';'CL11';'CL14';'CL24'},...
%     'Callback',@popup_callback1);
% 
% txt2 = uicontrol('Parent',d,...
%     'Style','text',...
%     'Position',[20 120 210 40],...
%     'String','Select a projection');
% 
% popup2 = uicontrol('Parent',d,...
%     'Style','popup',...
%     'Position',[100 100 100 25],...
%     'String',{'XY','YZ','XZ'},...
%     'Callback',@popup_callback2);
% 
% 
% 
% btn = uicontrol('Parent',d,...
%     'Position',[89 20 70 25],...
%     'String','Submit',...
%     'Callback','delete(gcf)');
% 
% %choice = 'CL101';
% 
% % Wait for d to close before running to completion
% uiwait(d);
% 
%     function popup_callback1(popup,~)
%         idx = popup.Value;
%         popup_items = popup.String;
%         % This code uses dot notation to get properties.
%         % Dot notation runs in R2014b and later.
%         % For R2014a and earlier:
%         % idx = get(popup,'Value');
%         % popup_items = get(popup,'String');
%         res.cluster = char(popup_items(idx,:));
%     end
% 
%     function popup_callback2(popup,~)
%         idx = popup.Value;
%         popup_items = popup.String;
%         % This code uses dot notation to get properties.
%         % Dot notation runs in R2014b and later.
%         % For R2014a and earlier:
%         % idx = get(popup,'Value');
%         % popup_items = get(popup,'String');
%         res.proj = char(popup_items(idx,:));
%     end
% 
% 
% end