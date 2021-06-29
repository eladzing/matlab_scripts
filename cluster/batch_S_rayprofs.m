boxx=8;

tc=S(boxx);
wt=RHOG(boxx);
len=8;
units;

profs=[];

rv=get_rvir();


for i=1:len
    [p1 p2 rp]=ray_profile('boxx',boxx,'data',tc,'wt',wt,'width',3);
    
    profs(end+1,:)=p1;
    profs(end+1,:)=p2;
    

    
end
tp=read_S_Profile(rp);

figure
loglog(rp./rv,profs.*f_ent,'--');
hold on
loglog(rp./rv,tp,'linewidth',2.5);hold off
xlim([0.1 2.5]);
%ylim([5e-3 3])
grid 
hold off
xlabelmine('$r/R_{vir}$')
ylabelmine('$S\,[\mathrm{KeV cm^2}]$')
titlemine('Ray Profiles for CL6')    