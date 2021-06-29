function [m_dot roprof vr_prof r_prof t_prof s_prof rvir mvir tvir vvir]= full_profs2(smallbox,bigbox,vcmflag,rvcm,hubflag)

%global FILE_FORMAT;
%FILE_FORMAT = sprintf('%s/%s', halopath, haloformat);
global HALO_PATH
global FILE_FORMAT_SPHERE;
global aexpn;

if ~exist('rvcm','var')
    rvcm=[0 1];
end
if ~exist('hubflag','var')
    hf=0;
else
    hf=(strcmpi(hubflag,'hub'));    
end

boxx=smallbox

full_ro =  RHOG_sphere(boxx);
full_ts =  T_sphere(boxx);
full_s = full_ts./full_ro.^(2/3);

if strcmp(vcmflag,'vcm') %%find Vrcm within Rvir
    rvc=rvcm.*get_rvir();
    vrcm=Vrcm_SPHERE(boxx,rvc,[0,0,0],hubflag);
else
    vrcm=zeros(size(full_ro));
end
clear rvc;

load(sprintf('%s/virial%d_%s', HALO_PATH, boxx,aexpn));
%vv=get_vvir(); 
full_ff = flux_sphere(boxx,vrcm,hubflag);
full_vr = Vr_sphere(boxx)+hf.*hubble_sphere(boxx)-vrcm;
tsh= [];
vrsh= [];
ssh= [];
rop= [];
md1= []; % flux
md2= []; % flux in
md3= []; % flux out
md4= []; % flux in strong

for ridx = 1:(length(R_Profile)-1)
    
    squee_vr=squeeze(full_vr(ridx,:,:));
    
    rosh=squeeze(full_ro(ridx,:,:));
    rovr=rosh.*(squee_vr);   %squeeze(full_vr(ridx,:,:)));
    rots=rosh.*(squeeze(full_ts(ridx,:,:)));
    ross=rosh.*(squeeze(full_s(ridx,:,:)));
    tsh(end+1)= sum(sum(rots))./sum(sum(rosh));
    vrsh(end+1) = sum(sum(rovr))./sum(sum(rosh));
    ssh(end+1) = sum(sum(ross))./sum(sum(rosh));
    rop(end+1) = mean(mean(rosh));
    
    inflow_mask = (squee_vr <= 0) ; 
    strong_inflow_mask = (squee_vr <= -0.5*VVIR); % If VVir isn't defined use vv (see above)
    
    ff = squeeze(full_ff(ridx,:,:));
    ff_in = ff(inflow_mask);
    ff_out = ff(~inflow_mask);
    ff_in_strong = ff(strong_inflow_mask);
    
    md1(end+1)= sum(ff(:));
    md2(end+1)= sum(ff_in(:));
    md3(end+1)= sum(ff_out(:));
    md4(end+1)= sum(ff_in_strong(:));
    
end

rp=R_Profile(1:(length(R_Profile)-1));

%%
if smallbox<bigbox
for ii=(log2(smallbox)+1):log2(bigbox)
    boxx=2^ii
    load(sprintf('%s/virial%d_%s', HALO_PATH, boxx,aexpn));
    full_ff = flux_sphere(boxx,vrcm,hubflag);
    full_vr = Vr_sphere(boxx)+hf.*hubble_sphere(boxx)-vrcm;
    full_ro =  RHOG_sphere(boxx);
    full_ts =  T_sphere(boxx);
    full_s = full_ts./full_ro.^(2/3);
    for ridx =(length(R_Profile)/2):(length(R_Profile)-1)
        squee_vr=squeeze(full_vr(ridx,:,:));
    
        rosh=squeeze(full_ro(ridx,:,:));
        rovr=rosh.*(squee_vr);   %squeeze(full_vr(ridx,:,:)));
        rots=rosh.*(squeeze(full_ts(ridx,:,:)));
        ross=rosh.*(squeeze(full_s(ridx,:,:)));
        tsh(end+1)= sum(sum(rots))./sum(sum(rosh));
        vrsh(end+1) = sum(sum(rovr))./sum(sum(rosh));
        ssh(end+1) = sum(sum(ross))./sum(sum(rosh));
        rop(end+1) = mean(mean(rosh));
    
        inflow_mask = (squee_vr <= 0) ; 
        strong_inflow_mask = (squee_vr <= -0.5*VVIR); % If VVir isn't defined use vv (see above)
    
        ff = squeeze(full_ff(ridx,:,:));
        ff_in = ff(inflow_mask);
        ff_out = ff(~inflow_mask);
        ff_in_strong = ff(strong_inflow_mask);
    
        md1(end+1)= sum(ff(:));
        md2(end+1)= sum(ff_in(:));
        md3(end+1)= sum(ff_out(:));
        md4(end+1)= sum(ff_in_strong(:));
            
    end
    rp=cat(2,rp,R_Profile(length(R_Profile)/2:length(R_Profile)-1));
end
end

if bigbox<8
    load(sprintf('%s/virial%d_%s', HALO_PATH, 8,aexpn));
end

t_prof=tsh;
r_prof=rp;
s_prof=ssh;
vr_prof=vrsh;
roprof=rop;
rvir=RVIR;
mvir=MVIR;
tvir=TVIR;
vvir=VVIR;

m_dot{1}=md1; % full flux
m_dot{2}=md2; % inflowing flux
m_dot{3}=md3; % outflowing
m_dot{4}=md4; % strong inflowing

end
