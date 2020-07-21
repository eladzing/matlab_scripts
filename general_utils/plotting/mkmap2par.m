function res = mkmap2par(mainPar,secPar,colMap,mainParLim,secParLim)
%MKMAP2PAR plot colormap where 2 parameters set the color & brigtness
%   prepare a map in which the color(hue& saturation are set by a main
%   parameter and the brightness is set by a secondary parameter. This is
%   acheived by converting the single value map of the main parameter in an
%   RGB version and then into a HSV. The secondary parameter is linearly
%   scaled to the range [0,1] and then used to set the value (3 index) of
%   the HSV map. The result is converted back to an RGB image which can be
%   displayed. A 'colobar' map is also produced. 
%   Arguments: 
%     mainPar - an MxN map of the main parameter
%     secPar  -  an MxN map of the main parameter
%     colMap  - color map NC X 3 color map to be used for plotting
%     mainParLim,secParLim - limits for Main and Secondary Parameters
%                            respectevely(OPTIONAL)
%   OUTPUT: 
%      rgbOut   -  new RGB image. 
%      rgbColorPlane  - 2D 'colorbar' with the Main parameter plotted on
%      the x Axis and the Secondary parameter plotted on the y axis 
%      mainParLim,secParLimK  - limits for Main and Secondary Parameters
%                            respectevely


%% set up limits if not given explicitly 
if ~exist('mainParLim','var')
    mainParLim=[min(mainPar(isfinite(mainPar))) max(mainPar(isfinite(mainPar)))];
end
if ~exist('secParLim','var')
    secParLim=[min(secPar(isfinite(secPar))) max(secPar(isfinite(secPar)))];
end


%% convert main parameter to RGB and then to HSV
rgbCB=real2rgb(mainPar,colMap,mainParLim);
hsv=rgb2hsv(rgbCB);

% figure;image(rgbCB);
% colormap(colMap);colorbar;
% title('Main Par. Original')
% 
%  rgbSec=real2rgb(secPar,colMap,secParLim);
% % 
%  figure;image(rgbSec);
%  colormap(colMap);colorbar;
%  title('Sec Par. Original')
%  clear rgbSec

%% scale secondary parameter to values between 0 and 1 
secParScaled=(secPar-secParLim(1))./diff(secParLim);


%% scale brightness by secondary parameter 
%hsv=hsv0;
hsv(:,:,3)=hsv(:,:,3).*(secParScaled);
%hsv(:,:,3)=secParScaled;%    1-secParScaled;
rgbOut=hsv2rgb(hsv);
%figure;image(rgbOut);


%% prepare colorplane 
nn=256;
xx=linspace(mainParLim(1),mainParLim(2),nn);
yy=linspace(secParLim(1),secParLim(2),nn);

[X,Y]=meshgrid(xx,yy);

rgbCB=real2rgb(X,colMap,mainParLim);
hsvCB=rgb2hsv(rgbCB);
Yscaled=(Y-secParLim(1))./diff(secParLim);

%hsv=hsv0;
hsvCB(:,:,3)=hsvCB(:,:,3).*(Yscaled);
%hsv(:,:,3)=Yscaled;

rgbColorPlane=hsv2rgb(hsvCB);

res.rgbOut=rgbOut;
res.rgbColorPlane=rgbColorPlane;
res.mainParLim=mainParLim;
res.secParLim=secParLim;


% figure;imagesc(xx,yy,rgbColorPlane);
% set(gca,'Ydir','normal')
%% test different options 

% 
% hsv1=hsv;
% hsv1(:,:,2)=secParScaled;
% rgb1=hsv2rgb(hsv1);
% figure;image(rgb1);
% title('direct Sat scaling')
% 
% hsv2=hsv;
% hsv2(:,:,2)=hsv2(:,:,2).*(secParScaled);
% rgb2=hsv2rgb(hsv2);
% figure;image(rgb2);
% title('factor Sat scaling')
% 
% hsv3=hsv;
% hsv3(:,:,3)=(1-secParScaled);
% rgb3=hsv2rgb(hsv3);
% figure;image(rgb3);
% title('factor value scaling')
% 
% 
% hsv4=hsv;
% hsv4(:,:,2)=hsv4(:,:,2).*(1-secParScaled);
% hsv4(:,:,3)=hsv4(:,:,3).*(secParScaled);
% rgb4=hsv2rgb(hsv4);
% figure;image(rgb4);
% title('factor sat & value')
% 



end

