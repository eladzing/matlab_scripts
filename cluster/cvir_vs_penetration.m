%% compare cvir to penetration 

%% find cvir  - quick and dirty 

list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

r0=0.01:0.001:1;

for i=1:length(list)
    
    new_env(list(i))
    
    rv=get_rvir;
    
    rp=r0.*rv;
    [rog,rotot]=read_RHO_Profiles(rp);
    
    rotot2=smooth(rotot,30);
    
    
    figure
    loglog(rp./rv,rotot,rp./rv,rotot2)
    hold on
    
    legend('reg','smooth');
    
    xlabelmine('$r\,[\mathrm{Mpc}]$')
    ylabelmine('$\rho\,[\mathrm{M_\odot/Mpc^3}]$')
    
    [dlro1,lr1]=derive1(log10(rotot),log10(rp));
    %[dlro2,lr2]=derive1(log10(rotot-rog),log10(rp));
    [dlro2,lr2]=derive1(log10(rotot2'),log10(rp));
    dlro3=smooth(dlro1,50);
    figure
    %plot(lr1-log10(rv),dlro1,lr2-log10(rv),dlro2)
    plot(rp(2:end-1)./rv,dlro1,rp(2:end-1)./rv,dlro2,rp(2:end-1)./rv,dlro3')
    grid
    legend('reg','smooth','smooth2');
    xlabelmine('$\log(r)\,[\mathrm{Mpc}]$')
    ylabelmine('$\log(\rho)\,[\mathrm{M_\odot/Mpc^3}]$')
   
    pause
    close all
end
    
    
    
