function md=md_from_moster(mstar,show)
%% Based on Moster et al. 2009, we match a stellar mass to halo mass

%mstar=10.^(8:0.1:12);

showFlag=false;

if exist('show','var')
    showFlag=strcmp(show,'show') | strcmp(show,'plot');
end


%% create SHM based on Mosteret al. values 
mh=10.^(7:0.01:15);
a=2;

lm1=11.884;
lm1Ps=lm1+a.*0.03;
lm1Ms=lm1-a.*0.023;

m1=10.^lm1;
m1Ps=10^lm1Ps;
m1Ms=10^lm1Ms;

mM0=0.0282;
mM0Ps=mM0+a.*0.000061;
mM0Ms=mM0-a.*0.000053;

beta=1.057;
betaPs=beta+a.*0.054;
betaMs=beta-a.*0.046;

gamma=0.556;
gammaPs=gamma+a.*0.01;
gammaMs=gamma-a.*0.004;


ms=mh.*2.*mM0.*( (mh./m1).^(-beta)+(mh./m1).^(gamma)).^-1;

msM=mh.*2.*mM0Ps.*( (mh./m1Ms).^(-betaMs)+(mh./m1Ps).^(gammaMs)).^-1;
msP=mh.*2.*mM0Ms.*( (mh./m1Ps).^(-betaPs)+(mh./m1Ms).^(gammaPs)).^-1;

%% find mean and range for given stellar mass

mHaloMean=interp1(ms,mh,mstar);
mHaloP=interp1(msP,mh,mstar);
mHaloM=interp1(msM,mh,mstar);

mean=log10(mHaloMean);
sig=max(abs(log10(mHaloP)-mean),abs(log10(mHaloM)-mean));
mHalo=10.^(mean+2.*sig.*randn(size(mstar)));

md=mstar./mHalo;


if showFlag

figure
loglog(ms,mh,'-b',msP,mh,'-r',msM,mh,'-g')
hold on 
loglog(mstar,mHalo,'.b')
xlim([1e8 1e12])

figure
semilogx(mstar,mstar./mHalo,'o')