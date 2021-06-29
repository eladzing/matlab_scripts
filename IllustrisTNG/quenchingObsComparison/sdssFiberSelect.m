function res = sdssFiberSelect(pos,mass,lenParam,mask)
%USDSSFIBERSELECT  - function to account for SDSS fiber effect on galaxy
%selection 
%   The SDSS fiber effect: if there are 2 gal's within 33'' of each other,
%   i.e., in the same fiber, only the more massive one is selected and the
%   other is ignored. This function takes a list of galaxies, sorts them by
%   mass and then goind done in mass, removes the less massive in pairs
%   which are closer than some given distance 
% arguments: 
% pos - should be projected, and normalized to box size! should be kpc, but in general could be
%       anything as long as it is consistent with lenght parameter 
% mass - mass of galaxies - units not important since it is a ranking
% lenParam - excising parameter (fiber size) should be consistent with
% position and normalized to box size 
% mask - optional mask for galaxy list 

if ~exist('mask','var')
    mask=true(size(mass));
end

%% sort by mass
indList=find(mask);
mass2=mass(mask);
pos2=pos(:,mask);

if length(lenParam)>1
    lpar=lenParam(mask);
else
    lpar=lenParam.*ones(1,sum(mask));
end

[mass2,ix]=sort(mass2,'descend');

indList2=indList(ix);

pos2=pos2(:,ix);

lpar=lpar(ix);

for i=1:length(mass2)
    
    if indList2(i)<0
        continue
    end
    % deal w/ periodic boundary 
    dx=abs(pos2(1,:)-pos2(1,i));
    dx=min(dx,1-dx);
    
    dy=abs(pos2(2,:)-pos2(2,i));
    dy=min(dy,1-dy);
    
    dist=hypot(dx,dy);
    
    inMask= dist<=lpar(i) & dist>0;
    
    indList2(inMask)=-1;
    
end
    
    
outList=sort(indList2);  
res.outList=outList(outList>0);
    
     
res.outMask=false(size(mask));
res.outMask(res.outList)=true;

end

