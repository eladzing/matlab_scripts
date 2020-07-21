boxx=8;

tc=T(boxx);
wt=RHOG(boxx);
len=8;


profs=[];


for i=1:len
    [p1 p2 rp]=ray_profile('boxx',boxx,'data',tc,'wt',wt,'width',3);
    
    profs(end+1,:)=p1;
    profs(end+1,:)=p2;
    

    
end
tp=read_T_Profile(rp);

figure
loglog(rp./rv,profs./tv,'--');
hold on
loglog(rp./rv,tp./tv,'linewidth',2.5);
xlim([0.1 2.5]);
ylim([5e-3 3])
grid 
hold off
xlabelmine('$r/R_{vir}$')
ylabelmine('T/T_{\mathrm{vir}}$')
titlemine('Ray Profiles for CL6')    