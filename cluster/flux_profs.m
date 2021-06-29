function [m_dot m_dot_in m_dot_out m_dot_in_strong vrprof vrprof_in vrprof_out vrprof_in_strong r_prof rvir mvir vvir]= flux_profs(smallbox,bigbox,vcmflag,rvcm,hubflag)

%global FILE_FORMAT;
%FILE_FORMAT = sprintf('%s/%s', halopath, haloformat);
global HALO_PATH
global FILE_FORMAT_SPHERE;
global aexpn;

if ~exist('rvcm')
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

%inflow_mask = full_vr <= 0 ; 
%strong_inflow_mask = full_vr <= -0.5*VVIR;


%ff_in = full_ff(inflow_mask);
%ff_out = full_ff(~inflow_mask);
%ff_in_strong = full_ff(strong_inflow_mask);

%vr_in = full_vr(inflow_mask);
%vr_out = full_vr(~inflow_mask);
%vr_in_strong = full_vr(strong_inflow_mask);

%tsh= [];
vrsh= [];
%ssh= [];
rop= [];
md1= [];
md2= [];
md3= [];
md4= [];

for ridx = 1:(length(R_Profile)-1)
    squee_vr=squeeze(full_vr(ridx,:,:));
    inflow_mask = (squee_vr <= 0) ; 
    strong_inflow_mask = (squee_vr <= -0.5*VVIR);
    
    ff = squeeze(full_ff(ridx,:,:));
    rosh=squeeze(full_ro(ridx,:,:));
    
    ff_in = ff(inflow_mask);
    ff_out = ff(~inflow_mask);
    ff_in_strong = ff(strong_inflow_mask);
    
    vr_in = squee_vr(inflow_mask);
    vr_out = squee_vr(~inflow_mask);
    vr_in_strong = squee_vr(strong_inflow_mask);
                   
    rovr=rosh.*(squeeze(full_vr(ridx,:,:)));
    %rots=rosh.*(squeeze(full_ts(ridx,:,:)));
    %ross=rosh.*(squeeze(full_s(ridx,:,:)));
    %tsh(end+1)= sum(sum(rots))./sum(sum(rosh));
    vrsh(end+1) = sum(sum(rovr))./sum(sum(rosh));
    %ssh(end+1) = sum(sum(ross))./sum(sum(rosh));
    rop(end+1) = mean(mean(rosh));
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
    %full_ts =  T_sphere(boxx);
    %full_s = full_ts./full_ro.^(2/3);
    for ridx =(length(R_Profile)/2):(length(R_Profile)-1)
        squee_vr=squeeze(full_vr(ridx,:,:));
        inflow_mask = (squee_vr <= 0) ; 
        strong_inflow_mask = (squee_vr <= -0.5*VVIR);
    
        ff = squeeze(full_ff(ridx,:,:));
        rosh=squeeze(full_ro(ridx,:,:));
    
        ff_in = ff(inflow_mask);
        ff_out = ff(~inflow_mask);
        ff_in_strong = ff(strong_inflow_mask);
    
           
        rovr=rosh.*(squeeze(full_vr(ridx,:,:)));
        vrsh(end+1) = sum(sum(rovr))./sum(sum(rosh));
        rop(end+1) = mean(mean(rosh));
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

m_dot=md1;
m_dot_in=md2;
m_dot_out=md3;
m_dot_in_strong=md4;
vrprof=vrsh;
vrprof_in=vrsh;
vrprof_out=vrsh;
vrprof_in_strong=vrsh;
r_prof=rp;
rvir=RVIR;
mvir=MVIR;
vvir=VVIR;

end
