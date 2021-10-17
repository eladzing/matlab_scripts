function res = calcMode(arr,nbins)
%CALCMODE _calculate the mode of non-integer distributions by constructing
%a histogram
%   for a given array of values build a histogram and identify the most
%   frequent value

%nbinsDef=20;



if ~exist('nbins','var')
    [nc,edges]=histcounts(arr);
else
     [nc,edges]=histcounts(arr,nbins);
end


[~,ix]=max(nc);

res=0.5.*(sum(edges(ix:ix+1)));




end

