%% draw cluster profiles for density, temperature & entropy


printoutdir='/home/eladzing/Ubuntu One/cluster/printout/profiles';

titletag='%s %s, $z=%d$';



list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
mask=ones(size(list));
%mask =[ 0   0   0   0   0   0   0  0 0 1 0 0  0  0  0  0];

ll=sum(mask);
roo=[];
too=[];
soo=[]; 
rroo=[];

%energyprofs=cell(ll,10);
   
%% global defineitions 
typ='csf';
aexp='a06';

%hubflag='hub';
%vcenterflag='defualt';
%center=[0,0,0];


% kb=1.38066e-16; %erg/K
% mu=0.5926; %mass per particle for primordial composition
% mp=1.672649e-24; % gram
% yr=3.155815e7; %  sec 
% pc=3.0856e18; % center
% km=1e5; %center
% Ms=1.989e33; %gr
% xi=0.5185185; 

%%numerical factors converting to units of M_sun*(km/sec)^2 per Mpc^3
% f_th=1.5.*kb./(mu.*mp)/km.^2 ;%E_th factor converts to units cited above
% f_pr=kb./(mu.*mp)/km.^2 ;%E_th factor converts to units cited above
% f_pv=kb./(mu.*mp.*km.*1.0e6.*pc).*yr; %units of M_sun*(km/sec)^2 per Mpc^3 per yr
% fu=km./(1e6.*pc).*yr;
% f_cl=(xi./(mu.*mp)).^2*yr*Ms/(1e6*pc)^3/km^2*1.e-22;

% in units of Rv 
%rsh=[0.1 1];


for id=1:length(list)
    if(~mask(id)) 
        continue
    end

    
    cluster=sprintf('CL%d',list(id))   

    %begin calculation 
    new_env(cluster,typ,aexp);

    %global HALO_PATH
    %global aexpn;
    global hub;
    global NCELL;   
    global zred;
    global CLUSTER;
        
    cellsize=1/hub/NCELL;
    rv=get_rvir();tv=get_tvir();
    %%ri=1:0.5:  (4/(1/NCELL));
    %%rr=ri.*cellsize;
    rrr=cellsize:cellsize/2:10/hub;
    rrr=rrr.*(1+zred);
    
    [roprof rotot_prof] = read_RHO_Profiles(rrr);
    tprof = read_T_Profile(rrr);
    sprof = read_S_Profile(rrr);
    
    roo(:,end+1)=roprof;
    soo(:,end+1)=sprof;
    too(:,end+1)=tprof;
    
    
    
    rr=rrr./rv;
    rroo(:,end+1)=rr;
    reslim=find_inner_reslim(1,1,3)./rv;
    
    xl=[reslim 5];
    
    name=sprintf('%s/%s_%s_profile_%s.%s',printoutdir,CLUSTER,'%s',aexp,'%s');
    
    figure
    loglog(rr,roprof,rr,rotot_prof-roprof,rr,rotot_prof,'linewidth',2)
    xlim(xl);
    grid
%     hold on
%     yl=ylim;
%     
%     loglog([rv rv],yl,'--')
%     loglog([4*rv 4*rv],yl,'--')
%     hold off
    h=legend('$\rho_{gas}$','$\rho_{DM}$','$\rho_{tot}$');
    set(h,'Interpreter','latex','FontSize',12);
    xlabelmine('$r\, [\mathrm{Mpc}]$')
    ylabelmine('$\rho\, [\mathrm{M_\odot / Mpc^3}]$')
    titlemine(sprintf(titletag,CLUSTER,'Density',zred));
    set(gca,'Fontsize',12) 
    
    %exportfig(gcf,sprintf(name,'rho','png'),'format','png');
    %exportfig(gcf,sprintf(name,'rho','eps'))
    
    
    
    figure
    loglog(rr,tprof,'linewidth',2)
    xlim(xl);
    grid
     hold on
%     yl=ylim;
     loglog(xl,[tv tv],'--','linewidth',1.5);
    xlim(xl);
%     loglog([3*rv 32*rv],yl,'--')
%     loglog([3*rv 3*rv],yl,'--')
%     hold off
    xlabelmine('$r\, [\mathrm{Mpc}]$')
    ylabelmine('$T\, [\mathrm{K}]$')
    titlemine(sprintf(titletag,CLUSTER,'Temperature',zred));
    set(gca,'Fontsize',12)
    %exportfig(gcf,sprintf(name,'tmp','png'),'format','png');
    %exportfig(gcf,sprintf(name,'tmp','eps'))
    
    
    
    
    figure
    loglog(rr,sprof,'linewidth',2);
    xlim(xl);
    grid
%     hold on
%     yl=ylim;
%     loglog([rv rv],yl,'--')
%     loglog([3*rv 3*rv],yl,'--')
%     hold off
%     set(h,'Interpreter','latex','FontSize',12);
    xlabelmine('$r\, [\mathrm{Mpc}]$')
    ylabelmine('$S\, [\mathrm{Kev\,cm^2}]$')
    titlemine(sprintf(titletag,CLUSTER,'Entropy',zred));
    set(gca,'Fontsize',12)
    %exportfig(gcf,sprintf(name,'ent','png'),'format','png');
    %exportfig(gcf,sprintf(name,'ent','eps'))
    
    close all
end

