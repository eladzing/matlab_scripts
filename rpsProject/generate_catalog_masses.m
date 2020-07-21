% bp=illustris.set_env(100);
% global DEFAULT_MATFILE_DIR
% load([DEFAULT_MATFILE_DIR '/fofs_subs_TNG100_z0.mat']

function res=generate_catalog_masses(groupMass,haloRange,msRange,Ngal)


mh=groupMass(log10(groupMass)>haloRange(1) & ...
log10(groupMass)<haloRange(2));

[~,ind]=sort(rand(size(mh)));

samp=mh(ind);

 

msSamp=stellarmass_from_moster(samp);

mask=log10(msSamp)>=msRange(1) & log10(msSamp)<=msRange(2);
stellarMass=msSamp(mask);
haloMass=samp(mask);

%[~,ind]=sort(rand(size(stellarMass)));
indx=(1:Ngal);
res.stellarMass=stellarMass(indx);
res.haloMass=haloMass(indx);

end
