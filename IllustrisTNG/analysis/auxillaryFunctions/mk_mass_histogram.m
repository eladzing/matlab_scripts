function [xx, massDist, mus]=mk_mass_histogram(val,mass,qus,leng)
%MK_MASS_HISTOGRAM Auxilary function for GENERATE_GASPROPERTIES
%   create cumulative mass histogram according to a given parameter. calculate 'mass
%   quantiles'.

if ~exist('leng','var')
    leng=200;
end

[bird, ~, xxlim]= histogram1d(val,mass,'len',leng);
len=size(bird,1);
xx0=linspace(xxlim(1),xxlim(2),len+1);
xx=xx0(1:end-1)+0.5.*diff(xx0);

massDist=bird(:,1);
massDistC=cumsum(bird(:,1))./sum(mass);
%xx=xxlim(1)+0.5*binsize:binsize:xxlim(end)-0.5*binsize;

mus=zeros(size(qus));
for j=1:length(qus)
    ind=find(massDistC>qus(j),1,'first');
    if ind>1
        mus(j)=interp1(massDistC(ind-1:ind),xx(ind-1:ind),qus(j));
    else
        mus(j)=xx(1);
    end
end

end