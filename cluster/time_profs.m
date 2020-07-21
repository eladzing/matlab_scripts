

set(0,'DefaultAxesColorOrder',[0 0 1;1 0 0;0 1 0])
global HALO_PATH
%%load(sprintf('%s/virial%d', HALO_PATH, 8));
[roprof vrprof r_prof rvir]= full_vrprof(HALO_PATH,smallbox,bigbox);

clear rvir;

[flprof rprof m200 mg200 r200 t200 v200]= flux_profile_allbox(HALO_PATH,smallbox);


thick=0.001;lth=10^thick;

RP=r200*0.001:thick:2.5*r200;

mdfull=interp1(rprof,flprof,RP,'spline');
vrp=interp1(r_prof,vrprof,RP,'spline');

clear flprof rprof roprof vrprof r_prof rvir;


Ggrav=4.43e-15;
zmet=1.;

Tp = read_T_Profile(RP);
[RHOG_Profile RHOTOT_Profile] = read_RHO_Profiles(RP);
[MG_Profile MSTAR_Profile MDM_Profile] = read_Mass_Profiles(RP);
Mtot= read_MTOT_Profile(RP);


[Mg200 ms200 M_dm200]=read_Mass_Profiles(r200);


if strcmp(plotflag,'pppppplot')
figure;
title(sprintf( '%s Profiles',clustername));
subplot(2,2,1);
loglog(RP/r200,Tp./t200);grid;xlim([RP(1) 3]);ylim([5e-2 6]);
xlabel('r/R_{vir}');ylabel('T/T_{vir}');title(sprintf('%s T Profile, T_{vir}=%sK',clustername,num2str(t200,'%1.2d')));


subplot(2,2,2);
rhovir=(3.*m200./(4.*pi.*r200.^3));
loglog(RP/r200,RHOG_Profile./rhovir);grid;xlim([RP(1) 3]);ylim([8e-3 1e4])
xlabel('r/R_{vir}');ylabel('\rho/\rho_{vir}');title(texlabel(sprintf('%s rho_{gas} Profile, rho_{vir}=%s M_{sun}/Mpc^3',clustername,num2str(rhovir,'%1.2d'))));
clear rhovir;

subplot(2,2,3);
loglog(RP./r200,MG_Profile/mg200);grid;xlim([RP(1) 3]);ylim([5e-3 3])
xlabel('r/R_{vir}');ylabel('M_{gas}(<r)/M_{gas}(<R_{vir})');title(sprintf('%s M_{gas}(<r) Profile, M_{gas}(R_{vir}=%s M_{sun}',clustername,num2str(Mg200,'%1.2d')));


subplot(2,2,4);
semilogx(RP./r200,vrp);grid;xlim([RP(1) 3]);ylim([-1000 200])
xlabel('r/R_{vir}');ylabel('V_r [km/sec]');title(sprintf('%s V_{r} Profile, V_{vir}=%s km/sec',clustername,num2str(v200,'%3.1f')));



if strcmp(pflag,'print')
        saveas(gcf,sprintf('%s/%s_profs.png',result_dir,clustername));
end
end
%%



cson=sqrt(1.52e8.*Tp);
tdyn1=(Ggrav.*Mtot.*RP.^(-3)).^(-1/2);
tdyn2=(4.*pi.*Ggrav.*RHOTOT_Profile).^(-1/2);

T6=Tp./1e6;

%figure;
%loglog(RP/RVIR,tdyn1,RP/RVIR,tdyn2);grid;axis tight
%xlabel('r/R_{vir}');ylabel('t_{dyn} [Gyr]');


lmb22=0.6.*(zmet/0.3).^0.7.*T6.^(-1)+0.02.*T6.^0.5;

tcool=2.61.*(RHOG_Profile.*6.77e-13).^(-1).*T6.*lmb22.^-1;
tcon=RP./abs(vrp).*(979.365);
tdyn3=RP./cson.*3.085e24./(3.15e16);
Mdp=2.*MG_Profile./(tcool.*1.e9);



%load(sprintf('%s/virial%d', halopath, 8));


%%plot data
%%plotflag='plot';
if strcmp(plotflag,'plot')

    figure;subplot(3,1,1:2);
    loglog(RP/r200,tcool,RP/r200,tcon,RP/r200,tdyn1,RP/r200,tdyn3);grid;xlim([RP(1) 3]);ylim([5e-3 7e3]) 
    set(gcf,'DefaultAxesColorOrder',[0 0 1;1 0 0;0 1 0;1 0 1])
    %xlabel('r/R_{vir}');
    ylabel(sprintf('t_{dyn},t_{cool},t_{con} [Gyr]'),'Fontsize',12);
    title(sprintf( '%s Dynamical, Contraction and Cooling Time Profile',clustername),'Fontsize',12);
    legend('t_{cool}','t_{con}', 't_{dyn}','r/c_s','Location','NorthWest');
    line([min(RP/r200) max(RP/r200)], [14 14],'Color','Black');
    subplot(3,1,3);
    loglog(RP/r200,tcon./tcool);grid;xlim([RP(1) 3]);ylim([5e-3 30])
    xlabel('r/R_{vir}','Fontsize',12);ylabel(sprintf('t_{con}/t_{cool}'),'Fontsize',12);
    set(gca,'fontsize',12)   
    
    if strcmp(pflag,'print')
        saveas(gcf,sprintf('%s/%s_times.png',result_dir,clustername));
    end
end
if strcmp(plotflag,'plot3333')
    %figure;
    %loglog(RP/r200,Tp./t200);grid;xlim([RP(1) 3]);
    %set(gcf,'DefaultAxesColorOrder',[0 0 1;1 0 0;0 1 0])
    %xlabel('r/R_{vir}');ylabel(sprintf('T/T_{vir}'));
    %title(sprintf( '%s Temperature Profile',clustername));
     
    %if strcmp(pflag,'print')
    %    saveas(gcf,sprintf('%s/%s_tprof.png',result_dir,clustername));
    %end

    figure;
    subplot(3,1,1:2)
    loglog(RP/r200,abs(mdfull),RP/r200,abs(Mdp));grid;xlim([RP(1) 3]);ylim([10 1e5])
    set(gcf,'DefaultAxesColorOrder',[0 0 1;1 0 0;0 1 0])
    %xlabel('r/R_{vir}','Fontsize',12);
    ylabel(sprintf('|dM_{gas}/dt|, |M_{gas}/t_{cool}| [M_{sun}/yr]'),'Fontsize',12);
    title(sprintf( '%s Mass Flux Profile',clustername),'Fontsize',12);
    legend('Mass flux', 'M_{gas}/t_{cool}','Location','SouthEast');
 
    subplot(3,1,3)
    loglog(RP/r200,abs(mdfull./Mdp));grid;xlim([RP(1) 3]);ylim([1e-2 1e2])
    xlabel('r/R_{vir}','Fontsize',12);ylabel(sprintf('|dM_{gas}/dt|/|M_{gas}/t_{cool}|'),'Fontsize',12);
       
    if strcmp(pflag,'print')
        saveas(gcf,sprintf('%s/%s_massflux2.png',result_dir,clustername));
    end
end
    
    
    
%clear all 
