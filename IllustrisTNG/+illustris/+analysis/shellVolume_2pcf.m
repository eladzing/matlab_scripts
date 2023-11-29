function res = shellVolume_2pcf(rProf,boxSize,vf)
%SHELLVOLUME_2PCF calculating the relevant shell volume for the calculation
%of the 2-point corllation function. 
%   to take advantge of the objects in a simulation box which are in the
%   regions beyond a sphere contained in the simulation cube, we need the
%   volume which intersects a shell of radius  greater than half the cube
%   side and the cube itself. Using a pre-prepared array from the
%   CUBE_IN_SHELL_FRACVOL function the actual volume is computed. 

if ~exist('boxSize','var')
    global LBox
    boxSize=LBox;
end

maxBoxRad=boxSize/2;
maxRad=sqrt(3)*maxBoxRad;
inRad=rProf(1:end-1);


totVolume= 4*pi/3*(rProf(2:end).^3-rProf(1:end-1).^3);
fracVol=ones(size(totVolume));

fracVol(inRad>=maxRad)=0;

%% 
inds=find(inRad>=maxBoxRad & inRad<=maxRad);


if ~isempty(inds)
    
    % load array 
    global DEFAULT_MATFILE_DIR
    load([DEFAULT_MATFILE_DIR '\fractionalShellVolume.mat'])
    
    %%  interpolate points 
    if inds(end)<length(rProf)
        inds=[inds inds(end)+1];
    end
    
    intGridTop=griddedInterpolant(shellFracVolume.rr,shellFracVolume.fracVolume,'makima');
    intGridBottom=griddedInterpolant(shellFracVolume.rr,shellFracVolume.totalVolume,'makima');
    
    nTop=intGridTop(rProf(inds)./maxBoxRad);
%     nTop=interp1(shellFracVolume.rr,shellFracVolume.fracVolume,rProf(inds)./maxBoxRad,'linear');
    nBottom=intGridBottom(rProf(inds)./maxBoxRad);
%     nBottom=interp1(shellFracVolume.rr,shellFracVolume.totalVolume,rProf(inds)./maxBoxRad,'linear');
        
    
        
    fracVol(inds(1:end-1))=diff(nTop)./diff(nBottom);

    fracVol(isnan(fracVol))=0;
end


res.shellVolume=totVolume.*fracVol;
res.fracVol=fracVol;
res.rProf=rProf;
res.rShell=0.5.*(rProf(2:end)+rProf(1:end-1));



end

