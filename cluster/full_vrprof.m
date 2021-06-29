function [roprof vrprof r_prof rvir]= full_vrprof(halopath,smallbox,bigbox)

%global FILE_FORMAT;
%FILE_FORMAT = sprintf('%s/%s', halopath, haloformat);
global FILE_FORMAT_SPHERE;
global aexpn;
FILE_FORMAT_SPHERE = sprintf('%s/%s', halopath, '%s_sphere_%d.mat');

boxx=smallbox
load(sprintf('%s/virial%d_%s', halopath, boxx,aexpn));
full_vr = Vr_sphere(boxx);
full_ro =  RHOG_sphere(boxx);
vrsh= [];
rop= [];

for ridx = 1:(length(R_Profile)-1)
    rosh=squeeze(full_ro(ridx,:,:));
    rovr=(squeeze(full_ro(ridx,:,:))).*(squeeze(full_vr(ridx,:,:)));
    vrsh(end+1) = sum(sum(rovr))./sum(sum(rosh));
    rop(end+1) = mean(mean(rosh));
end

rp=R_Profile(1:(length(R_Profile)-1));

%%
if smallbox<bigbox
for ii=(log2(smallbox)+1):log2(bigbox)
    boxx=2^ii
    load(sprintf('%s/virial%d_%s', halopath, boxx,aexpn));
    full_vr = Vr_sphere(boxx);
    full_ro =  RHOG_sphere(boxx);
    
    for ridx =(length(R_Profile)/2):(length(R_Profile)-1)
        rosh=squeeze(full_ro(ridx,:,:));
        rovr=(squeeze(full_ro(ridx,:,:))).*(squeeze(full_vr(ridx,:,:)));
        vrsh(end+1) = sum(sum(rovr))./sum(sum(rosh));
        rop(end+1) = mean(mean(rosh));
    end
    rp=cat(2,rp,R_Profile(length(R_Profile)/2:length(R_Profile)-1));
end

%ind = find(rp<(0.01*RVIR));
%l=length(rp);lb=ind(length(ind))+1; %dx=0.01*(rp(l)- rp(lb))./RVIR;

r_prof=rp;
vrprof=vrsh;
roprof=rop;
rvir=RVIR;
end
