%% build orbit bank

%% define host
%   1e12  5e12  1e13  5e13  1e14  5e14  1e15 - 3e15/5e15
mh=[12.22 12.83 13.23 13.83 14.19 14.81 15.12];
kstart=[3  2    2     1     1     1     1];
hostTags={'1e12' '5e12' ...
    '1e13' '5e13' ...
    '1e14' '5e14' '1e15'};
for kk=1:length(mh)
    
    fprintf('Generating orbits for host of log(M)=%f \n',mh(kk))
    
    mhost=10.^mh(kk);
    cvh=cvir_Mvir(mhost,0);
    fgh=0.15;
    host=NFW('mv',mhost,'cc',cvh,'fg',fgh,'mean');
    
    %% generate velocities
    nOrb=1000;
    for k=kstart(kk):3
        
        switch(k)
            case(1)
                massRatio=0.0025;
                mrTag='0005';
            case 2
                massRatio=0.025;
                mrTag='005';
            case 3
                massRatio=0.25;
                mrTag='05';
        end
        
        
        vr=generateRadialVelocity(nOrb,mhost,massRatio);
        vv=generateTotallVelocity(nOrb,mhost,massRatio);
        
        %     %% enforce v<2Vvir limit and keep total number of orbits
        %     vr=vr(vv<=2);
        %     vv=vv(vv<=2);
        %
        %     df=nOrb-length(vv);
        %     cnt=0;
        %     while df>0 || cnt>20
        %         cnt=cnt+1;
        %
        %         vr1=generateRadialVelocity(df,mhost,massRatio);
        %         vv1=generateTotallVelocity(df,mhost,massRatio);
        %
        %         vv=cat(2,vv,vv1);
        %         vr=cat(2,vr,vr1);
        %
        %         df=nOrb-length(vv);
        %     end
        
        % find v_tangetial
        vt=vv.*sqrt(1-vr.^2);
        
        % create vx vy
        vx=-vr.*vv.*host.Vvir;
        vy=vt.*host.Vvir;
        
        %% build orbits
        y0=0;
        x0=host.Rvir;
        
        tOrbit=2*x0./host.Vvir;
        tmax=2*tOrbit;
        dt=tOrbit/2e4;
        dtmin=dt/1e3;
        
        %% calculate orbits
        orbIC.x=x0;
        orbIC.y=y0;
        
        step=10;
        prc0=10;
        nOrb=length(vx);
        
        fprintf('Calculating %s orbits for mass ratio %s \n',num2str(nOrb),mrTag)
        for i=1:nOrb
            
            prc=i./nOrb*100;
            if prc>=prc0
                fprintf('completed %s %% of orbits \n',num2str(prc0))
                prc0=prc0+step;
            end
            
            
            orbIC.vx=vx(i);
            orbIC.vy=vy(i);
            
            orb(i)=orbits.rk4_orbitIntegration(orbIC,dt,tmax,dtmin,@orbits.rhs_nfw,host);
        end
        
        global DEFAULT_MATFILE_DIR
        name=sprintf('orbBank%i_mean_hostMass%s_massRatio%s.mat',nOrb,hostTags{kk},mrTag);
        save([DEFAULT_MATFILE_DIR '/' name],'orb','-v7.3')
        fprintf('saving to: %s \n',[DEFAULT_MATFILE_DIR '/' name])
        
    end
    
end
