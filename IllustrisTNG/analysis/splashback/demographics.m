%% We're looking at the demoraphics of "splashback galaxies - 
% centrals which were 

%% load relavent things 

bp=illustris.set_env(100,'nomount');
load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/fofs_subs_TNG100_z0.mat')


global simDisplayName
global DEFAULT_MATFILE_DIR
global illUnits

galMass=subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit;

snap=99;
massThresh=1e9;
galMask=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'centrals');
ssfr=illustris.utils.calc_ssfr(subs);

sfMask=ssfr>=1e-11;
qMask=~sfMask;

%% look back as a function of time 
ftype='Time';
fname=sprintf('%s/splashbackCatalog_by%s_%s',DEFAULT_MATFILE_DIR,ftype,simDisplayName);
load(fname,'splashbackTime');

% only by sat
binE=9:0.25:12.25;

for i=1:5
    sats= splashbackTime.isSat(1:i,:)==1;
    sats=any(sats,1);
    
    [hst(i).counts, hst(i).edges]=histcounts(log10(galMass(sats & galMask)),binE);
    [hst(i).Qcounts, hst(i).Qedges]=histcounts(log10(galMass(sats & galMask & qMask)),binE);
    
    satNum(i)=sum(sats);
    satNumQ(i)=sum(sats(qMask));
    %satNumSF(i)=sum(sats(sfMask));
    
    
    fars= splashbackTime.isFar_50prc(1:i,:)==1;
    fars=any(fars,1);
    
    [hstF(i).counts, hstF(i).edges]=histcounts(log10(galMass(fars & galMask)),binE);
    [hstF(i).Qcounts, hstF(i).Qedges]=histcounts(log10(galMass(fars & galMask & qMask)),binE);
    
    bothNum(i)=sum(sats & fars);
    bothNumQ(i)=sum(sats(qMask) &  fars(qMask));
    %farNumSF(i)=sum(sats(sfMask));
    
    [hstB(i).counts, hstB(i).edges]=histcounts(log10(galMass(sats & fars & galMask)),binE);
    [hstB(i).Qcounts, hstB(i).Qedges]=histcounts(log10(galMass(sats & fars & galMask & qMask)),binE);
    
    farNum(i)=sum(fars);
    farNumQ(i)=sum(fars(qMask));
    %farNumSF(i)=sum(sats(sfMask));
    
    
    
    
    
    
    
end

