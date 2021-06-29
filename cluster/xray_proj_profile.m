%% extract profule from 2D Image. Use recursive division of cells into sub cells 

boxx=1; %in Mpc
global NCELL
global hub 
imaj=projI;

%% generate bins - should be larger than cell size
cellsize = (boxx/hub)/NCELL;
binFac=3; 
binSize0=binFac*cellsize;

rmin=0.0;
rmax=boxx/hub/2;

nBins=floor((rmax-rmin)/binSize0);
binSize=(rmax-rmin)/nBins;

binEdges=rmin:binSize:rmax;
bins=binEdges(1:end-1)+diff(binEdges)./2;

%% create  distance array 

cm=[0 0];
[meshX, meshY] = meshgrid(1:NCELL);

meshX = meshX - (NCELL+1)/2 -cm(1);
meshY = meshY - (NCELL+1)/2 -cm(2);
% Fix Units (to be in Mpc)
meshX = meshX * ((boxx/hub)/NCELL);
meshY = meshY * ((boxx/hub)/NCELL);


rcube=sqrt(meshX.^2+meshY.^2) ; % r cube in Mpc
%% divide image into bins 
profile=zeros(nBins,1);
for i=1:nBins 
    mask=rcube>=binEdges(i) & rcube<binEdges(i+1);
    
    profile(i)=sum(imaj(mask))/(pi *( binEdges(i+1)^2-binEdges(i)^2))*cellsize^2  ;
end

