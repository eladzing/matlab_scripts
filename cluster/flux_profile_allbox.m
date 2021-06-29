function [flprof rprof m200 mg200 r200 t200 v200]= flux_profile_allbox(halopath,smallbox)

%global FILE_FORMAT;
%FILE_FORMAT = sprintf('%s/%s', halopath, haloformat);
global FILE_FORMAT_SPHERE;
FILE_FORMAT_SPHERE = sprintf('%s/%s', halopath, '%s_sphere_%d.mat');

load(sprintf('%s/virial%d', halopath, smallbox));
full_ff = flux_sphere(smallbox);
M_DOT    = [];

for ridx = 1:(length(R_Profile)-1)
    ff = squeeze(full_ff(ridx,:,:));
    M_DOT(end+1) = sum(ff(:));
end

rp=R_Profile(1:(length(R_Profile)-1));

%%
for ii=smallbox:3
    boxx=2^ii
    load(sprintf('%s/virial%d', halopath, boxx));
    full_ff = flux_sphere(boxx);
    
    for ridx =(length(R_Profile)/2):(length(R_Profile)-1)
        ff = squeeze(full_ff(ridx,:,:));
        M_DOT(end+1) = sum(ff(:));
    end
    rp=cat(2,rp,R_Profile(length(R_Profile)/2:length(R_Profile)-1));
end

%%
rprof=rp;
flprof=M_DOT;

mg200=read_MGAS_Profile(halopath,RVIR);
m200=MVIR;
r200=RVIR;
t200=TVIR;
v200=VVIR
%m_norm=-560*(MVIR./1e13).^(0.15).*(MGAS./1e13)  %virial accretion rate in Msun/yr

