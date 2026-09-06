function [xx, massDist, mus]=mk_mass_histogram(val,mass,varargin)
%MK_MASS_HISTOGRAM Auxilary function for GENERATE_GASPROPERTIES
%   create cumulative mass histogram according to a given parameter. calculate 'mass
%   quantiles'.

% defaults 
qus=[0.1 0.25 0.5 0.75 0.9];
leng=200;
xlim=[];

% parse arguments
i=1;
while i<=length(varargin)
    switch(lower(varargin{i}))
        case {'length','leng','len'}
            i=i+1;
            leng=varargin{i};
        case {'qus','quantiles','quants'}
            i=i+1;
            qus=varargin{i};
        case {'lims','xlims','xlim','lim'}
            i=i+1;
            xlim=varargin{i};
        otherwise
            error('%s - unknown argument:%s',current_function().upper,varargin{i})
    end
    i=i+1;
end

if isempty(xlim)
    [bird, ~, xxlim]= histogram1d(double(val),double(mass),'len',leng);
else
    [bird, ~, xxlim]= histogram1d(double(val),double(mass),'len',leng,'xlim',xlim);
end

len=size(bird,1);
xx0=linspace(xxlim(1),xxlim(2),len+1);
xx=xx0(1:end-1)+0.5.*diff(xx0);

massDist=bird(:,1);
massDistC=cumsum(bird(:,1))./sum(mass);
%xx=xxlim(1)+0.5*binsize:binsize:xxlim(end)-0.5*binsize;

%% Calculate the mass quantiles 
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