function sigmaSqCube=mk_sigma_squared(box) % sigma square in (km/sec)^2

%units;
nSmoothBase=5;
%smoothLength=2*nSmooth+1;
smoothSigBase=2;
nSmooth=ceil(nSmoothBase./box);
smoothSig=ceil(smoothSigBase./box);
smoothLength=2*nSmooth+1;

[vx,vy,vz]=get_velocities(box);
ro=RHOG(box);

ros=smooth3(ro,'gaussian',smoothLength,smoothSig);
%vx=vx.*km;vy=vy.*km;vz=vz.*km;

vxs=smooth3(vx.*ro,'gaussian',smoothLength,smoothSig)./ros;
vys=smooth3(vy.*ro,'gaussian',smoothLength,smoothSig)./ros;
vzs=smooth3(vz.*ro,'gaussian',smoothLength,smoothSig)./ros;

vx2s=smooth3(vx.^2.*ro,'gaussian',smoothLength,smoothSig)./ros;
vy2s=smooth3(vy.^2.*ro,'gaussian',smoothLength,smoothSig)./ros;
vz2s=smooth3(vz.^2.*ro,'gaussian',smoothLength,smoothSig)./ros;

sigmaSqX=vx2s-vxs.^2;
sigmaSqY=vy2s-vys.^2;
sigmaSqZ=vz2s-vzs.^2;
sigmaSqCube=sigmaSqX+sigmaSqY+sigmaSqZ;
