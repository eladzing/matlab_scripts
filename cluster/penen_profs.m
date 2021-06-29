function [m_dot m_in m_hlfvr roprof vrprof vrin r_prof rvir mvir vvir]= penen_profs(smallbox,bigbox,vcmflag,rvcm,hubflag)

%global FILE_FORMAT;
%FILE_FORMAT = sprintf('%s/%s', halopath, haloformat);
global HALO_PATH
global FILE_FORMAT_SPHERE;
global aexpn;

if ~exist('rvcm','var')
    rvcm=[0 1];
end
if ~exist('hubflag')
    hf=0;
else
    hf=(strcmpi(hubflag,'hub'));    
end

boxx=smallbox

full_ro =  RHOG_sphere(boxx);
%full_ts =  T_sphere(boxx);
%full_s = full_ts./full_ro.^(2/3);

if strcmp(vcmflag,'vcm') %%find Vrcm within Rvir
    rvc=rvcm.*get_rvir();
    vrcm=Vrcm_SPHERE(boxx,rvc,[0,0,0],hubflag);
else
    vrcm=zeros(size(full_ro));
end
clear rvc;

load(sprintf('%s/virial%d_%s', HALO_PATH, boxx,aexpn));
full_ff = flux_sphere(boxx,vrcm,hubflag);
full_vr = Vr_sphere(boxx)+hf.*hubble_sphere(boxx)-vrcm;
%tsh=[];
vrsh=[];
%ssh=[];
rop=[];
md=[];
mdin=[];
mdvr=[];
vris=[];

vvi=get_vvir();
for ridx = 1:(length(R_Profile)-1)
    ff = squeeze(full_ff(ridx,:,:));
    vr = squeeze(full_vr(ridx,:,:));
    md(end+1)= sum(ff(:));
    mdin(end+1)= sum(ff(ff<0));
    mdvr(end+1)= sum(ff(vr<-0.5.*vvi));
    rosh=squeeze(full_ro(ridx,:,:));
    rovr=rosh.*(squeeze(full_vr(ridx,:,:)));
    %rots=rosh.*(squeeze(full_ts(ridx,:,:)));
    %ross=rosh.*(squeeze(full_s(ridx,:,:)));
    %tsh(end+1)= sum(sum(rots))./sum(sum(rosh));
    vrsh(end+1) = sum(sum(rovr))./sum(sum(rosh));
    vris(end+1) = sum(sum(rovr(vr<0)))./sum(sum(rosh(vr<0)));
    %ssh(end+1) = sum(sum(ross))./sum(sum(rosh));
    rop(end+1) = mean(mean(rosh));
end

rp=R_Profile(1:(length(R_Profile)-1));

%%
if smallbox<bigbox
for ii=(log2(smallbox)+1):log2(bigbox)
    boxx=2^ii
    load(sprintf('%s/virial%d_%s', HALO_PATH, boxx,aexpn));
    full_vr = Vr_sphere(boxx)+hf.*hubble_sphere(boxx)-vrcm;
    full_ff = flux_sphere(boxx,vrcm,hubflag);
    full_ro =  RHOG_sphere(boxx);
    %full_ts =  T_sphere(boxx);
    %full_s = full_ts./full_ro.^(2/3);
    for ridx =(length(R_Profile)/2):(length(R_Profile)-1)
        ff = squeeze(full_ff(ridx,:,:));
        vr = squeeze(full_vr(ridx,:,:));
        md(end+1)= sum(ff(:));
        mdin(end+1)= sum(ff(ff<0));
        mdvr(end+1)= sum(ff(vr<-0.5.*vvi));
        rosh=squeeze(full_ro(ridx,:,:));
        rovr=rosh.*(squeeze(full_vr(ridx,:,:)));
%         rots=rosh.*(squeeze(full_ts(ridx,:,:)));
%         ross=rosh.*(squeeze(full_s(ridx,:,:)));
%         tsh(end+1)= sum(sum(rots))./sum(sum(rosh));
        vrsh(end+1) = sum(sum(rovr))./sum(sum(rosh));
        vris(end+1) = sum(sum(rovr(vr<0)))./sum(sum(rosh(vr<0)));
%         ssh(end+1) = sum(sum(ross))./sum(sum(rosh));
         rop(end+1) = mean(mean(rosh));
%         md(end+1)= sum(ff(:));
    end
    rp=cat(2,rp,R_Profile(length(R_Profile)/2:length(R_Profile)-1));
end
end

if bigbox<8
    load(sprintf('%s/virial%d_%s', HALO_PATH, 8,aexpn));
end

%t_prof=tsh;
r_prof=rp;
%s_prof=ssh;
vrprof=vrsh;
vrin=vris;
roprof=rop;
rvir=RVIR;
mvir=MVIR;
%tvir=TVIR;
vvir=VVIR;
m_dot=md;
m_in=mdin;
m_hlfvr=mdvr;
end
