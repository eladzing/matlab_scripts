function [full_ff full_ro full_ts full_ss rp m200 r200 t200 v200]=catspheres(smallbox,bigbox,vcmflag,rvcm,hubflag)

%global FILE_FORMAT_SPHERE;
%FILE_FORMAT_SPHERE = sprintf('%s/%s', halopath, '%s_sphere_%d.mat');
global HALO_PATH
%global FILE_FORMAT_SPHERE;
global aexpn
%load(sprintf('%s/virial%d', HALO_PATH, 8));

% if ~exist('rvcm')
%     rvcm=[0 1];
% end
% if ~exist('hubflag')
%     hf=0;
% else
%     hf=(strcmpi(hubflag,'hub'));    
% end

boxx=smallbox
load(sprintf('%s/virial%d_%s', HALO_PATH, boxx, aexpn));
ro = RHOG_sphere(boxx);
ts = T_sphere(boxx);
ss = ts./ro.^(2/3); %S_sphere(boxx);  %ts./ro.^(2/3);

% if strcmp(vcmflag,'vcm') %%find Vrcm within Rvir
%     rvc=rvcm.*get_rvir();
%     vrcm=Vrcm_SPHERE(boxx,rvc,[0,0,0],hubflag);
% else
%     vrcm=zeros(size(ro));
% end
% clear rvc;

ff = flux_sphere(boxx); %,vrcm,hubflag);
%vr = Vr_sphere(boxx)+hf.*hubble_sphere(boxx)-vrcm;



lastind = length(R_Profile)-1;
rp=R_Profile(1:lastind);

full_ff = ff(1:lastind,:,:);
full_ro = ro(1:lastind,:,:);
full_ts = ts(1:lastind,:,:);
%full_vr = vr(1:lastind,:,:);
full_ss = ss(1:lastind,:,:);

if smallbox<bigbox
for ii=(log2(smallbox)+1):log2(bigbox)
    boxx=2^ii
    load(sprintf('%s/virial%d_%s', HALO_PATH, boxx,aexpn));
    ff = flux_sphere(boxx);%,vrcm,hubflag);
    ro = RHOG_sphere(boxx);
    ts = T_sphere(boxx);
 %   vr = Vr_sphere(boxx)+hf.*hubble_sphere(boxx)-vrcm;
    ss = ts./ro.^(2/3); %S_sphere(boxx);%ts./ro.^(2/3);

    
    ind=find(R_Profile<rp(end));
    indx=ind(end)+1;
    lastind = length(R_Profile)-1;
    
    rp=cat(2,rp,R_Profile(indx:lastind));
    full_ff = cat(1,full_ff,ff(indx:lastind,:,:));
    full_ro = cat(1,full_ro,ro(indx:lastind,:,:));
    full_ts = cat(1,full_ts,ts(indx:lastind,:,:));
  %  full_vr = cat(1,full_vr,vr(indx:lastind,:,:));
    full_ss = cat(1,full_ss,ss(indx:lastind,:,:));
end

m200=get_mvir;
r200=get_rvir;
t200=get_tvir;
v200=get_vvir;

end
