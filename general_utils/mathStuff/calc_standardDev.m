function res= calc_standardDev(pnts,wts )
%CALC_STANDARDDEV  Calculate the standard deviation for a bunch of points.
%   Calculation of the standard deviation. Also - weighted standard
%   deviation

np=length(pnts);
if np>1
        
    if ~exist('wts','var')
        wts=ones(size(pnts));
    end
    
    if any(wts<0)
        error('calc_standardDev - negative weights?')
    end
    
    nw=sum(wts~=0);
    
    if nw<2
        error('calc_standardDev - must have more than 1 non-zero weight.')
    end
    
    wmeen=sum(pnts.*wts)./sum(wts);
    
    res=sqrt(sum(wts.*(pnts-wmeen).^2)./((nw-1)/nw.*sum(wts)));
    
else
    res=0;
    
end

