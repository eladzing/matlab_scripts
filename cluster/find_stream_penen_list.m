%% plot streams according to conditions
list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

pen=zeros(16,4,2);
%pen06=zeros(16,4);

for i=1:length(list)
    
    for k=1:2
        %k=1;
        switch k
            case 1
                new_env(list(i));
            case 2
                new_env(list(i),'a06')
        end
        
        units;
        
        %global SIMTYPE
        %global CLUSTER;
        global VCM
        %global aexpn;
        global NCELL;
        global zred
        global hub;
        
        RVIR=get_rvir();
        MVIR=get_mvir();
        %VVIR=get_vvir();
        TVIR=get_tvir();
        
        vcm = VCM;
        cm=[0,0,0];
        
        
        for bx=1:4;
            tic;
            boxx=2^(bx-1)
            
            
            
            %% read data
            
            % get velocity
            [hubX,hubY,hubZ] = hubble_flow(boxx,[0,0,0]);
            vx = Vx(boxx)+hubX-vcm(1);
            vy = Vy(boxx)+hubY-vcm(2);
            vz = Vz(boxx)+hubZ-vcm(3);
            
            vc=sqrt(vx.^2+vy.^2+vz.^2);
            
            % get positions
            [meshY, meshX, meshZ] = meshgrid(1:size(vy,1), 1:size(vx,2), 1:size(vz,3));
            %convert to center origin coordinates
            meshX = meshX - (size(vx,1)+1)/2 -cm(1);
            meshY = meshY - (size(vy,2)+1)/2 -cm(2);
            meshZ = meshZ - (size(vz,3)+1)/2 -cm(3);
            % Fix Units (to be in Mpc)
            meshX = meshX * ((boxx/hub)/NCELL);
            meshY = meshY * ((boxx/hub)/NCELL);
            meshZ = meshZ * ((boxx/hub)/NCELL);
            
            rcube=sqrt(meshX.^2+meshY.^2+meshZ.^2) ; % r cube in Mpc
            
            % read density
            rog=RHOG(boxx);
            [roShell,~]=read_RHO_Profiles(rcube);
            ros=RHOTOT(boxx)-RHODM(boxx)-rog;
            
            
            % read flux
            vr=(vx.*meshX+vy.*meshY+vz.*meshZ)./rcube ; %radial velocity
            
            fluxnorm=(0.1117.*(MVIR./1e15).^0.15*(1+zred)^2.25)./(4.*pi); % normalization for flux
            [Mg200 , ~, ~]=read_Mass_Profiles(RVIR);
            flux=rog.*vr.*rcube.^2./Mg200.*(km/Mpc*Gyr)./fluxnorm ; % normalized
            
            
            
            % read temperature
            tmp=T(boxx);
            tShell=read_T_Profile(rcube);
            
            % entropy
            ent=S(boxx);
            entShell=read_S_Profile(rcube);
            
            % mach
            gm=5/3;
            cso=(gm*kb/mm.*tShell).^0.5;  % ./1e5;  % so
            cso=cso./km;
            
            mach=vc./cso;
            
            %smooth
            smoothLen=100; %in kpc
            smLen=ceil(smoothLen./(1e3.*boxx./hub).*NCELL);
            smLen=smLen-1+mod(smLen,2);
            flux=smooth3(flux,'box',smLen);
            mach=smooth3(mach,'box',smLen); %/cso;
            ross=smooth3(ros,'box',smLen);
            ents=smooth3(ent,'box',smLen);
            
            %tmp=smooth3(tmp,'box',smLen); %/cso;
            
            %% set conditions
            fluxThresh=-3;
            fMask=flux<fluxThresh;
            
            roRatioThresh=1;
            roMask=rog./roShell>=roRatioThresh;
            
            tmpRatioThresh=1;
            tmpMask=tmp./tShell<tmpRatioThresh;
            
            entRatioThresh=1.0;
            entMask=ents./entShell<=entRatioThresh;
            
            machThresh=0.9;
            machMask=mach>machThresh;
            
            roStarThresh=0;
            rosMask=ross<=roStarThresh;
            %% build radius array
            
            radM=rcube( machMask & fMask & entMask & rosMask);
            if ~isempty(radM)
                pen(i,bx,k)=min(radM)./RVIR;
                %pen06(i,bx)=min(radM)./RVIR;
            else
               pen(i,bx,k)=1;
                %pen06(i,bx)=1;
            end
            toc
        end
    end
end

penStruct.penArr=pen;
penStruct.fluxThresh=fluxThresh;
penStruct.roRatioThresh=roRatioThresh;
penStruct.tmpTatioThresh=tmpRatioThresh;
penStruct.roStarThresh=roStarThresh;
penStruct.machThresh=machThresh;
penStruct.entThresh=entRatioThresh;
penStruct.smooth=smoothLen;
penStruct.note='ent mach and flux';