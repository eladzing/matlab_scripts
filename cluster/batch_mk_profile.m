%% constructs profiles for all the halos.  

flist=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
%mask=[ 0   0   1   0   0   1   0   0  0  1  0  0  0  1  1  0];
mask=true(size(flist));
%mask(1:2)=true;

pflag='noprint';
cflist=[104 105 107 3 6 7 10 14;4 5 7 8 10 11 13 15];
ncflist=[101 102 103 106 5 9 11 24;1 2 3 6 9 12 14 16];

ll=sum(mask);
roo=[];
too=[];
soo=[];
rroo=[];
ress=[];
flxPr=[];
rPr=[];

%% global defineitions
%typ='csf';


profile=struct();

for k=1:2
    switch k
        case 2
            aexp='a06';
        case 1
            aexp='a1';
    end
    for id=1:length(flist)
        if(~mask(id))
            continue
        end
        
        cluster=sprintf('CL%d',flist(id))
        profile(k).cluster(id).name=cluster;
        
        %begin calculation
        new_env(cluster,aexp);
        
        global HALO_PATH
        global aexpn;
        global hub;
        global NCELL;
        global zred;
        %global CLUSTER;
        
        profile(k).zred=zred;
        
        cellsize=1/hub/NCELL; % comoving;
        rv=get_rvir();
        tv=get_tvir();
        mv=get_mvir();
        
        %rhovir=3.*mv./(4.*pi.*rv.^3);
        %fluxnorm=0.056.*(mv./1e13).^0.15;
        %svir=(tvir./rhovir.^(2/3));
        [Mgas, ~,~] = read_Mass_Profiles(rv);
        
        profile(k).cluster(id).mvir=mv;
        profile(k).cluster(id).rvir=rv;
        profile(k).cluster(id).tvir=tv;
        profile(k).cluster(id).mgas=Mgas;
        
        %% prepare profiles from profile files
        rrr=cellsize:cellsize/2:10/hub; %comoving
        rrr=rrr./(1+zred);  %proper
        
        [roprof, rotot_prof] = read_RHO_Profiles(rrr);
        tprof = read_T_Profile(rrr);
        sprof = read_S_Profile(rrr);
        
        roo(:,end+1)=roprof;
        soo(:,end+1)=sprof./tv;
        too(:,end+1)=tprof./tv;
        
        
        
        
        rr=rrr./rv;
        rroo(:,end+1)=rr;
        reslim=find_inner_reslim(1,8)./rv;
        ress(end+1)=reslim;
        %xl=[reslim 5];
        
        profile(k).cluster(id).rProf=rr;
        profile(k).cluster(id).rhoProf=roprof;
        profile(k).cluster(id).tmpProf=tprof./tv;
        profile(k).cluster(id).sProf=sprof./tv;
        profile(k).cluster(id).reslim=reslim;
        
        %% prepare flux
        fprof=zeros(4,NCELL);
        fprofIn=zeros(4,NCELL);
        fprofOut=zeros(4,NCELL);
        rp=zeros(4,NCELL);
        
        for i=1:4
            b=2^(i-1);
            
            % create flux
            %vr_sph=cart2sphere(Vr_full(b));
            %f=new_flux(b,'vr_sph',vr_sph); %   flux_sphere(b);
            
            %save_cube(f, HALO_PATH, sprintf('flux_sphere_%s_%d',aexpn,b));
            %save_cube(vr_sph, HALO_PATH, sprintf('Vr_sphere_%s_%d',aexpn,b));
            
            %f=flux_sphere(b);
            f=new_flux(b);
            inMask=f<0;
            fIn=f.*inMask;
            fOut=f.*~inMask;
            
            
            
            fprof(i,:)=sum(sum(f,2),3);
            fprofIn(i,:)=sum(sum(fIn,2),3);
            fprofOut(i,:)=sum(sum(fOut,2),3);
            clear f fIn fOut
            %prepare radius
            dr=(b/2)/NCELL;
            rp(i,:)=dr:dr:b/2;
        end
        rrp=cat(2,rp(1,1:end-1),rp(2,128:end-1),rp(3,128:end-1),rp(4,128:end-1));
        ffprop=cat(2,fprof(1,1:end-1),fprof(2,128:end-1),fprof(3,128:end-1),fprof(4,128:end-1));
        ffpropIn=cat(2,fprofIn(1,1:end-1),fprofIn(2,128:end-1),fprofIn(3,128:end-1),fprofIn(4,128:end-1));
        ffpropOut=cat(2,fprofOut(1,1:end-1),fprofOut(2,128:end-1),fprofOut(3,128:end-1),fprofOut(4,128:end-1));
        rrp=rrp./hub;
        
        rpr=(rrp(1):cellsize/2:4/hub);
        flxp=interp1(rrp,ffprop,rpr,'PCHIP')./Mgas.*1e9; % in units of 1/Gyr
        flxpIn=interp1(rrp,ffpropIn,rpr,'PCHIP')./Mgas.*1e9; % in units of 1/Gyr
        flxpOut=interp1(rrp,ffpropOut,rpr,'PCHIP')./Mgas.*1e9; % in units of 1/Gyr
        flxPr(:,end+1)=flxp;
        rPr(:,end+1)=rpr./(1+zred)./rv;
        
        profile(k).cluster(id).rProflux=rpr./(1+zred)./rv;
        profile(k).cluster(id).fluxProf=flxp;
        profile(k).cluster(id).fluxProfIn=flxpIn;
        profile(k).cluster(id).fluxProfOut=flxpOut;
        
        
        %name=sprintf('%s/%s_%s_profile_%s.%s',printoutdir,CLUSTER,'%s',aexp,'%s');
        
    end
end



