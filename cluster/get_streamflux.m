function [rps,fls]=get_streamflux(clname,smallbox,bigbox)

global HALO_PATH;
streamfluxfile='%s/%s_streamflux_z0_b%d.dat';


boxx=smallbox;


fname=sprintf(streamfluxfile,HALO_PATH,clname,boxx);

[rps,fls]=textread(fname,'%f %f');

if smallbox<bigbox
for ii=(log2(smallbox)+1):log2(bigbox)
    boxx=2^ii;
    fname=sprintf(streamfluxfile,HALO_PATH,clname,boxx);
    [rp,fl]=textread(fname,'%f %f');
    rps=cat(1,rps,rp(ceil(length(rp)/2):length(rp)));
    fls=cat(1,fls,fl(ceil(length(fl)/2):length(fl)));
end
end