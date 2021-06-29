%% generalized jeans equation 
% Generate the terms in the generalized Jeans equation 
list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

rp=0.01:0.001:4/0.7; %till size of biggest box
units; 
%rp=rp.*Mpc;
%smooth parameters
nSmoothBase=5;
%smoothLength=2*nSmooth+1;
smoothSigBase=2;

az='a1';

for j=1:length(list)
    tic
% for k=1:2
%     switch k
%         case 1 
%             az='a1';
%         case 2
%             az='a06';
%     end
%          

new_env(list(j),az)
%rp=r0.*get_rvir.*Mpc;

global NCELL
global hub

%% calculate gravity 

mtot=read_MTOT_Profile(rp);
[rog,~]=read_RHO_Profiles(rp);

grav=-G.*(mtot.*Ms).*(rog.*Ms./Mpc.^3)./(rp.*Mpc).^2; % in cgs

%% calculate pressure
tm=read_T_Profile(rp);
pressTherm=kb/(mp.*muMass).*(rog.*Ms./Mpc.^3).*tm;

%% calculate turbulent pressuer
boxs=[1 2 4 8];

for i=1:length(boxs)
    nSmooth=ceil(nSmoothBase./boxs(i));
    smoothSig=ceil(smoothSigBase./boxs(i));
    smoothLength=2*nSmooth+1;
    %vr=Vr_full(boxs(i)).*km;
    [vx,vy,vz]=get_velocities(boxs(i));
    ro=RHOG(boxs(i));
    ros=smooth3(ro,'gaussian',smoothLength,smoothSig);
    vx=vx.*km;vy=vy.*km;vz=vz.*km;
        
    vxs=smooth3(vx.*ro,'gaussian',smoothLength,smoothSig)./ros;
    vys=smooth3(vy.*ro,'gaussian',smoothLength,smoothSig)./ros;
    vzs=smooth3(vz.*ro,'gaussian',smoothLength,smoothSig)./ros;
     
    vx2s=smooth3(vx.^2.*ro,'gaussian',smoothLength,smoothSig)./ros;
    vy2s=smooth3(vy.^2.*ro,'gaussian',smoothLength,smoothSig)./ros;
    vz2s=smooth3(vz.^2.*ro,'gaussian',smoothLength,smoothSig)./ros;
     
%     ros=smooth3(ro,'box',smoothLength);%,smoothSig);
%     vrs=smooth3(vr.*ro,'box',smoothLength)./ros;%,smoothSig)./ros;
%     vr2s=smooth3(vr.^2.*ro,'box',smoothLength)./ros;%,smoothSig)./ros;
%      
    sigmaSqX=vx2s-vxs.^2;
    sigmaSqY=vy2s-vys.^2;
    sigmaSqZ=vz2s-vzs.^2;
    sigmaSq=sigmaSqX+sigmaSqY+sigmaSqZ;
    sigSq=MAKE_PROFILE_FROM_CUBE_DENSE(ro.*sigmaSq)./MAKE_PROFILE_FROM_CUBE_DENSE(ro);
        
    cl=0.5.*(boxs(i)./hub)/NCELL;
    
    sig(i).box=boxs(i);
    sig(i).sigProf=sigSq;
    sig(i).rp=linspace(cl/2,0.5.*boxs(i)./hub-cl/2,NCELL);
    
end
sigProf0=cat(2,sig(1).sigProf(1:end-1),sig(2).sigProf(128:end-1),sig(3).sigProf(128:end-1),sig(4).sigProf(128:end));
rpSig0=cat(2,sig(1).rp(1:end-1),sig(2).rp(128:end-1),sig(3).rp(128:end-1),sig(4).rp(128:end));
sigProf=interp1(rpSig0,sigProf0,rp);

pressTurb=(rog.*Ms./Mpc.^3).*sigProf;

[fTurb,~]=derive1(pressTurb,rp.*Mpc);
[fTherm,rn]=derive1(pressTherm,rp.*Mpc);


Jean(j).grav=grav(2:end-1);
Jean(j).fTherm=fTherm;
Jean(j).fTurb=fTurb;
Jean(j).press=pressTherm;
Jean(j).Turb=pressTurb;
Jean(j).rp=rp(2:end-1);
Jean(j).rn=rn;
Jean(j).sigProf=sigProf;
Jean(j).rog=rog;
Jean(j).rv=get_rvir;
Jean(j).sigma=sig;
Jean(j).cl=list(j);
toc
end
   
    
    
