% build orbit bank 

mv0=5e14;
cv0=cvir_Mvir(mv0,0);
fg0=0.15;
    
% set host
    
host=NFW('mv',mv0,'cc',cv0,'fg',fg0);

tOrbit=2*pi*host.Rvir./host.Vvir;
tmax=15.34; 
dt=tOrbit/1e4;
dtmin=dt/1e3;   

fac=0:0.1:1;

IC.x=1.5.*host.Rvir;
IC.y=0;
IC.vx=0; %-host.Vvir;%./sqrt(2);

for i=1:length(fac) 
    IC.vy=fac(i).*vcirc(host,1.5.*host.Rvir,'kpc');%  
    
    
    orbitBank.orb(i)=orbits.rk4_orbitIntegration(IC,dt,tmax,dtmin,@orbits.rhs_nfw,host);
    orbitBank.IC(i)=IC;
end

    orbitBank.host=host;
    orbitBank.fac=fac;