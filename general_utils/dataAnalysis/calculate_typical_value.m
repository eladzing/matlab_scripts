function res =calculate_typical_value(arr,mass,nbins)
%CALCULATE TYPICAL VALUE  function for calculating the ypical
%value of a given array of vlues based on several metrics : 
% median, mean and mode  - explores value distribution 
% mass-weighted median, mean and mode - explores what is the typical value
% of the mass distribution. 

if ~exist('nbins','var')
    nbins=100;
end

% value distribution
res.arrMean=mean(arr);
res.arrMedian=median(arr);
res.arrMode=calcMode(arr,nbins);

% mass weighted 
res.arrMassMean=sum(arr.*mass)./sum(mass);

[xx, massDist, mus]=mk_mass_histogram(arr,mass,0.5,nbins);
res.arrMedianMass=mus(2);

[~,ix]=max(massDist);
res.arrModeMass=xx(ix);


end

