function res = cube_in_shell_fracVol(nSample0,rSample)
%CUBE_IN_SHELL_FRACVOL find the fractional volume of a shell of given
%radius and thickness which is found within a cube of given length
%   If a spherical shell is of radius greater than half the cube length but
%   smaller than sqrt(3)/2 of the cube length only part of the spherical
%   shell is found within the cube volume. This function computes the
%   fractional volume of the spherical shell within the cube volume.
%   Usefule for calculating the representational volume of for a 2-point
%   correlation function estimation. The volume is found by generating
%   sample points on the shells of increasing radius between r=1 and r=sqrt(3)and finding the fraction which is found
%   within the cube. The cube has a side 2 (-1 to 1) on all sides. The function ourputs two arrays:
%    intrsectVolume - the total number of points fonund in the intersecting
%                     volume cumulative from r to sqrt(3)
%    totalVolume   - the total number of points fonund in the sphere
%                     volume cumulative from r to sqrt(3)
%    rProf        - the radii of teh shells of points
%    To find the fractional volume of a shell intersecting a cube of radius
%    R and width dr you must find the relavent indices i1,i2 in the rProf
%    such that rProf(i1)=R and rProf(i2)-rProf(i1)=dr. The fractional voume is given by:
%    v_frac=diff(intrsectVolume([i2 i1])./diff(totalVolume([i2 i1]).
%
%    or eve better, by interpolation.

%    arguments :
%    nSamnple - no. of sample points on unit sphere in each azimuthat direction (default is 100 (10,000 points)
%    rsample - no. of sampling shells in the r- direction (defualt =100 ;
%
%% validate arguments

% if rad
% rmin=1;
% rmax=sqrt(3);

if ~exist('nSample0','var')
    nSample0=100;
end

if ~exist('rSample','var')
    rSample=100;
end


%%
rr=linspace(1,sqrt(3),rSample);
nSample=round(nSample0.*rr);
phase=(2.*pi./nSample).*rand(1,rSample);

xxx=[];
yyy=[];
zzz=[];

%% generate unit sphere sampleing
nsamp=nSample(1);
smplePnt=generateUnitSphereSampling(nsamp);
stp=10;
prc=stp;
for i=1:rSample
    
    if i>=(prc/100)*rSample
        fprintf('completed %i %% \n',prc);
        prc=prc+stp;
    end
    
    
    if nSample(i)>nsamp
        nsamp=nSample(i);
        smplePnt=generateUnitSphereSampling(nsamp);
    end
    %% generate points in shell and convert to cartesian
    
    % add a random shift in phi to get better coverage
    
    theta=smplePnt.theta;
    phi=smplePnt.phi+phase(i);
    
    
    %convert sample points to cartesian
    xx=rr(i).*sin(theta).*cos(phi);
    yy=rr(i).*sin(theta).*sin(phi);
    zz=rr(i).*cos(theta);
    
    xxx=cat(2,xxx,xx);
    yyy=cat(2,yyy,yy);
    zzz=cat(2,zzz,zz);
end

%% calculate fraction volume.
sampleRad=hypot3(xxx,yyy,zzz);
mskC=abs(xxx)<=1 & abs(yyy)<=1 & abs(zzz)<=1;
for i=1:rSample
    mskR=sampleRad>=rr(i); %  & sampleRad<rr(i+1);
    
    res.fracVolume(i)=sum(mskR & mskC);
    res.totalVolume(i)=sum(mskR);
end
res.rr=rr;

