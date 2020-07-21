%% plot a stripping as a function of cluster radius 

%% prepare satellite 

ms=10^10;
beta=1;
fg=0.1;
fb=0;
md=0.1/(1+fg+fb);
res=rscale_mmw_array(ms,'beta',beta,'fg',fg,'fb',fb,...
'md',md,'noshow');
sigma=ms./(2*pi*res.rd^2);


%% prepare cluster 

Mc=1e15;
fc=0.1;
r=0.01:0.01:100;
etap=3:-0.05:0.05 ;
for i=1:length(etap)  % in units of RV
    pval=rps_factor_expdisk('etap',etap(i),'sigma',sigma,'fd',fg,'fc',fc,'Mc',Mc);
   pv(i)=pval;
    f=disk_force_reduced(r,'fg',fg,'beta',beta,'fb',fb);
    
    if max(f)>pval
        rStr=interp1(f,r,pval);
        mStr(i)=exp_disk_mass(rStr,beta);
    else
        mStr(i)=0;
    end
    
end