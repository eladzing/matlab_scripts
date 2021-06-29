%% plot streams according to conditions
list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

units

vMean10=zeros(16,2);
vMean25=zeros(16,2);
rv=zeros(16,2);
mv=zeros(16,2);
vv=zeros(16,2);
to10=vv;
to25=vv;
%pen06=zeros(16,4);

for i=1:length(list)
    
    for k=1:1 %2
        %k=1;
        switch k
            case 1
                new_env(list(i));
            case 2
                new_env(list(i),'a06')
        end
        
        
        
        %global SIMTYPE
        %global CLUSTER;
        global VCM
        %global aexpn;
        global NCELL;
        %global zred
        global hub;
        
        RVIR=get_rvir();
        MVIR=get_mvir();
        VVIR=get_vvir();
        %TVIR=get_tvir();
        
        vcm = VCM;
        cm=[0,0,0];
        
        rv(i,k)=RVIR;
        mv(i,k)=MVIR;
        vv(i,k)=VVIR;
        
        switch list(i)
            case{101,102,103}
                bxPar=2;
            otherwise
                bxPar=1;
        end
        
        for bx=1:bxPar
            

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
            
            % make mask
            
            if bx==1
                mask10=rcube<=0.1*RVIR;
                
                vv1=vc(mask10);
                ro1=rog(mask10);
                
                vMean10(i,k)=sum(vv1.*ro1)./sum(ro1);
            end
            
            if bxPar==2 && bx==1
                continue
            end
            
            mask25=rcube<=0.25*RVIR;
            
            
            vv2=vc(mask25);
            ro2=rog(mask25);
            
            vMean25(i,k)=sum(vv2.*ro2)./sum(ro2);
            
            
            
        end
    end
end
for k=1:2
    to10(:,k)=(0.1.*rv(:,k).*Mpc)./(vMean10(:,k).*km)./(yr);
    to25(:,k)=(0.25.*rv(:,k).*Mpc)./(vMean25(:,k).*km)./(yr);
end
