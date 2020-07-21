function res = generateSphericalSampling(nSample,varargin)
%gENERATESPHERICALsAMPLING - generate sampling points to fill a sphere
% Sampling is preformed on a units sphere such that the points are
% uniformly distributed (see generateUnitSphereSamplping ). Points are then
% extended according to radial binning supplied in arguments.
%
% output is in cartesian coordinates. Optin exists for printing to file.

%% read arguments

% defualts
%nSample=100;

nbins=100;
rmin=0;
rmax=inf;
rRange=[];
rProf=[];
outFile='';
%cellsize=1;

%binFlag=false;  %% if true, use user supplied bins

%% get arguments
i=1;
while i<=length(varargin)
    
    switch lower(varargin{i})
        case{'rmin'}
            i=i+1;
            rmin=varargin{i};
        case{'rmax'}
            i=i+1;
            rmax=varargin{i};
        case{'range','rrange','rads'}
            i=i+1;
            rRange=varargin{i};
        case {'rbins','bins','nb','nbins'}
            i=i+1;
            nbins=varargin{i};
        case {'rprof'}
            i=i+1;
            rProf=varargin{i};
        case {'file','outputfile','output','outfile'}
            i=i+1;
            outFile=varargin{i};
            if isempty(outFile)
                error('generateSphericalSmapling - no file name given');
            end
        otherwise
            error('generateSphericalSmapling - Illegal argument %s; ',varargin{i});
    end
    i=i+1;
end


%% generate radial bins
if isempty(rProf)
    
    if isempty(rRange)
        if isinf(rmax)% set rmax if unset
            error('generateSphericalSmapling - must enter either radial array or rmax')
        end
        
        rRange=[rmin rmax];
    else
        if length(rRange)~=2
            error('generateSphericalSmapling - rRange must be a two value array');
        end
        
        if diff(rRange)<0
            tm=rRange(1);
            rRange(1)=rRange(2);
            rRange(2)=tm;
            clear tm
        end
    end
    rProf=rRange(1):diff(rRange)/(nbins-1):rRange(2);
end

%binSize=diff(rp);

%% generate unit sphere sampling
unitSphere=generateUnitSphereSampling(nSample);

%% extend to entire sphere

xx=zeros(length(rProf),length(unitSphere.area));
yy=xx;
zz=xx;
for i=1:length(rProf)
    
    xx(i,:)=rProf(i).*sin(unitSphere.theta).*cos(unitSphere.phi);
    yy(i,:)=rProf(i).*sin(unitSphere.theta).*sin(unitSphere.phi);
    zz(i,:)=rProf(i).*cos(unitSphere.theta);
end

res.rProf=rProf;
res.xs=xx;
res.ys=yy;
res.zs=zz;
%% print out result to file.

if ~isempty(outFile)
    xxx=reshape(xx',[numel(xx) 1]);
    yyy=reshape(yy',[numel(yy) 1]);
    zzz=reshape(zz',[numel(zz) 1]);
    rrr=hypot3(xxx,yyy,zzz);
    arr=cat(2,rrr,xxx,yyy,zzz);
    %arr=
    
    header='## radius x y z';
    fid=fopen(outFile,'w');
    fprintf(fid,'%s \n',header);
    fprintf(fid,'%12.5g %12.5g  %12.5g  %12.5g \n',arr');
    
    fclose(fid);
end
end

