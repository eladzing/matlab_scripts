%% draw cluster profiles for density, temperature & entropy

list=[101 102 103 104 105 106 107  3  5  6  7  9 10 11 14 24];
mask=[ 0   0   1   0   0   1   0   0  0  1  0  0  0  1  1  0];
mask=true(size(list))
%printoutdir='/home/eladzing/Ubuntu One/cluster/printout/profiles';

titletag='%s %s, $z=%d$';


ll=sum(mask);
roo=[];
too=[];
soo=[]; 
rroo=[];
ress=[];
%energyprofs=cell(ll,10);
   
%% global defineitions 
typ='csf';
aexp='a06';


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
    soo(:,end+1)=sprof./tv;
    too(:,end+1)=tprof./tv;
    
    
    
    rr=rrr./rv;
    rroo(:,end+1)=rr;
    reslim=find_inner_reslim(1,8)./rv;
    ress(end+1)=reslim;
    xl=[reslim 5];
    
    name=sprintf('%s/%s_%s_profile_%s.%s',printoutdir,CLUSTER,'%s',aexp,'%s');
    
end







name='%s/cledge_%s_profs_%s.%s';
xl=([max(ress) 5]);

%figure
    
    loglog(rroo,soo,'linewidth',2);
    ylim([1.0e-6 1.0e-3]);
    xlim(xl);
    grid
    hold on
    yl=ylim;
    loglog([1 1],yl,'--k','linewidth',1)
    hold off
      
    hl=legend('CL103','CL106','CL6','CL11','CL14');
    set(hl,'Interpreter','latex','FontSize',12,'location','NorthWest');
    xlabelmine('$ r/R_{vir}$')
    ylabelmine('$S/T_{vir}, [\mathrm{Kev\,cm^2\,K^{-1}}]$')
    titlemine(sprintf('%s  $z=%3.2g$','Entropy',zred));
    set(gca,'Fontsize',12)
    exportfig(gcf,sprintf(name,printoutdir,'ent',aexp,'png'),'format','png');
    exportfig(gcf,sprintf(name,printoutdir,'ent',aexp,'eps'))
    %pause
 figure
    
    loglog(rroo,too,'linewidth',2);
    ylim([1e-3 1e1]);   
    xlim(xl);
    grid
    hold on
    yl=ylim;
    loglog([1 1],yl,'--k','linewidth',1)
    loglog(xl,[1 1],'--k','linewidth',1)
    hold off
    
    hl=legend('CL103','CL106','CL6','CL11','CL14');
    set(hl,'Interpreter','latex','FontSize',12,'location','SouthWest');
    xlabelmine('$r/R_{vir}$')
    ylabelmine('$T/T_{vir}$')
    titlemine(sprintf('%s  $z=%3.2g$','Temperature',zred));
    set(gca,'Fontsize',12)
    exportfig(gcf,sprintf(name,printoutdir,'tmp',aexp,'png'),'format','png');
    exportfig(gcf,sprintf(name,printoutdir,'tmp',aexp,'eps'))