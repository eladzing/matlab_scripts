clust=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

%a='a06';
rp=0.4:0.0001:0.7;
for i=1:length(clust)
    new_env(clust(i),a);
    global zred
    global hub
    rr=rp.*get_rvir;
    
    mp=read_MTOT_Profile(rr);
    
    rhom=mp./(4*pi/3.*rr.^3);
   
    ref=500.*rho_crit(zred,hub);
    %ref=deltavir(zred).*rho_mean(zred);%
    
    r500c(i)=interp1(rhom,rr,ref);
end
    
    

