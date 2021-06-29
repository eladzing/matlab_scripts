list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
aExp='a1';
pdir='C:\Users\eladzing\Documents\cluster\printout\pressure_gradient_maps';
for i=1:length(list)
    
    new_env(list(i),aExp)
    
    press=RHOG(1).*T(1)./1e20;
    [gx,gy,gz]=gradient(pp);
    dP=sqrt(gx.^2+gy.^2+gz.^2);
    
    mkmap(1,'data',dP,'log','print','printag','dpress','outputdir',pdir,'png');
    mkmap(2,'data',dP,'log','print','printag','dpress','outputdir',pdir,'png');
    
    close all
end
 
