function res = generate_halo_masses(nHalo,varargin)
%GENERATE_HALO_MASSES generate a random ensemble of halo masses (M_200,c) based on TNG100 Mass function 
%   given a mass range, randomly generate halo masses based on the TNG100
%   mass function (fofs) 

path='/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles';

massRange0=[10 15];
massR=[];
massType='m200c';
source='TNG300';
i=1;
while i<=length(varargin)
    switch lower(varargin{i})
        case {'massrange','range','masses'}
            i=i+1;
            massR=varargin{i};
        case {'m200m','mean200','200mean','mean','200m'}
            massType='m200m';
        case {'m200c','crit200','200crit','crit','200c'}
            massType='m200c';
        case {'m500c','500','500crit','500c'}
            massType='m500c';
        case {'100','tng100'}
            source='TNG100';
        case {'300','tng300'}
            source='TNG300';
        case {'draco','isaac'}
            path='/isaac/ptmp/gc/eladzing/IllustrisTNG/matFiles';
        otherwise
            error('GENERATE_HALO_MASSES - Illegal argument: %s',varargin{i});
    end
    i=i+1;
end



%% arrange mass range 

if isempty(massR)
    massRange=massRange0;
else
    
    if any(massR>1e6)
        massR=log10(massR);
    end
    
    if length(massR)==1
        massRange=[max(massRange0(1),massR) min(massRange0(2),massR)];
    elseif length(massR)==2
        if diff(massR<=0)
            massR=fliplr(massR);
        end
        massRange=massR; %[max(massRange0(1),massR(1)) min(massRange0(2),massR(2))];
    else
         error('GENERATE_HALO_MASSES - Illegal mass range: %s',massR);
    
    end
    
end




% % if looking for halos beyond 10^14 then use TNG300; 
% if any(massRange>14)
%     source='TNG300';
% end 


%% load TNG mass function 



name=sprintf('fofMasses_%s.mat',source);

switch massType
    case 'm200c'
        load([path '/' name],'M200c');
        masses=log10(double(M200c));
    case 'm200m'
        load([path '/' name],'M200m');
        masses=log10(double(M200m));
    case 'm500c'
        load([path '/' name],'M500c');
        masses=log10(double(M500c));
end


masses=masses(masses>=massRange(1) & masses<=massRange(2));

%% generate cdf 

[nc, edges]=histcounts(masses,100);

mm=0.5.*(edges(1:end-1)+edges(2:end));
cdf=cumtrapz(mm,nc)./trapz(mm,nc);


%% generate masses 

[cdf, mask] = unique(cdf);
mm=mm(mask);

rv=rand(1,nHalo);
mas=interp1(cdf,mm,rv);

%% enfors the mass range limits 
mas=mas(mas>massRange(1) & mas<massRange(2));

df=nHalo-length(mas);
cnt=0;
while df>0 || cnt>20
    cnt=cnt+1;
    rv=rand(1,df);
    mas1=interp1(cdf,mm,rv);
    
    mas=cat(2,mas,mas1);
    mas=mas(mas>massRange(1) & mas<massRange(2));
    
    df=nHalo-length(mas);
end

res=10.^mas;


end



% cdf1=cumsum(nc);
% cdf1=cdf1./cdf1(end);
% 
% cdf21=cumtrapz(edges(1:end-1),nc)./trapz(edges(1:end-1),nc);
% cdf22=cumtrapz(edges(2:end),nc)./trapz(edges(2:end),nc);

