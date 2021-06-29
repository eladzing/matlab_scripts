%% draw cluster profiles for density, temperature & entropy

list=[101 102 103 104 105 106 107  3  5  6  7  9 10 11 14 24];
mask=[ 0   0   0   0   0   0   0   0  0  0  0  0  0  0  1  0];

%printoutdir='/home/eladzing/Ubuntu One/cluster/printout/profiles';

titletag='%s %s, $z=%d$';
load('C:\Users\eladzing\Documents\cluster\matlab\mat_files\shockedge.mat')
%plot_shockedge

ll=sum(mask);
roo=[];
too=[];
soo=[]; 
rroo=[];
ress=[];
%energyprofs=cell(ll,10);
   
%% global defineitions 
typ='csf';
aexp='a1';

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
    %roo(:,end+1)=roprof;
    %soo(:,end+1)=sprof./tv;
    %too(:,end+1)=tprof./tv;
    
    
    
    rr=rrr./rv;
    %rroo(:,end+1)=rr;
    reslim=find_inner_reslim(1,8)./rv;
    %ress(end+1)=reslim;
    xl=[0.1 5];
    
    name=sprintf('%s/%s_%s_cledge_profile_%s.%s','%s',CLUSTER,'%s',aexp,'%s');
    

%name='%s/cledge_%s_profs_%s.%s';
%xl=([max(ress) 5]);
switch aexp
    case 'a1'
        e1=shockedge_a1{id,4};%./rv_a1(i);
        e2=shockedge_a1{id,5};%./rv_a1(i);
        
        e1=e1(shockedgeMask1_a1(id,:));
        e2=e2(shockedgeMask2_a1(id,:));
            ed=cat(2,e1,e2);

    case 'a06'
        e1=shockedge_a06{id,4};%./rv_a1(i);
        e2=shockedge_a06{id,5};%./rv_a1(i);

        e1=e1(shockedgeMask1_a06(id,:));
        e2=e2(shockedgeMask2_a06(id,:));
        ed=cat(2,e1,e2);
end
    
    shk=ed;
    shk=shk(shk>0);
    shk=cat(1,shk,shk);


figure
    
    loglog(rr,sprof./(kb*tv./(1000.*ev)),'linewidth',2);
    %ylim([5.0e-6 1.0e-4]);
     ylim([5e1 1e4])
    xlim(xl);
    grid
    hold on
    loglog(shk./rv,ylim,'--k','linewidth',1.5)
    hold off
      
    %hl=legend('CL103','CL106','CL6','CL11','CL14');
    %set(hl,'Interpreter','latex','FontSize',12,'location','NorthWest');
    xlabelmine('$ r/R_{\mathrm{vir}}$')
    %ylabelmine('$S/T_{\mathrm{vir}}\,[\mathrm{Kev\,cm^2\,K^{-1}}]$')
    ylabelmine('$S/k_{\mathrm{B}}T_{\mathrm{vir}}\,[\mathrm{cm^2}]$');
    %titlemine(sprintf('%s %s  $z=%3.2g$',CLUSTER,'Entropy',zred));
    set(gca,'Fontsize',14)
    %exportfig(gcf,sprintf(name,printoutdir,'ent','png'),'format','png');
    %exportfig(gcf,sprintf(name,printoutdir,'ent','eps'))
    %pause
    printout_fig(gcf,sprintf('CL%s_entProf_%s',num2str(list(id)),aexp));
figure
     
    loglog(rr,tprof./tv,'linewidth',2);
    ylim([5e-3 5e0]);   
    xlim(xl);
    grid
    hold on
    loglog(shk,ylim,'--k','linewidth',1)
        
    hold off
    
    %hl=legend('CL103','CL106','CL6','CL11','CL14');
    %set(hl,'Interpreter','latex','FontSize',12,'location','SouthWest');
    xlabelmine('$r/R_{vir}$')
    ylabelmine('$T/T_{vir}$')
    titlemine(sprintf('%s %s $z=%3.2g$',CLUSTER,'Temperature',zred));
    set(gca,'Fontsize',12)
 %   exportfig(gcf,sprintf(name,printoutdir,'tmp','png'),'format','png');
 %   exportfig(gcf,sprintf(name,printoutdir,'tmp','eps'))

end