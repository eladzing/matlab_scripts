%% build orbit bank

% 1. start with known distribution of orbital parameters at r200,c
% 2. integrate *outwards* to r200,m
% 3. Use orbital parameters at r200,m and recalculate orbits inwards. 



setEnv_RPS

global cosmoStruct

%% define host based on 'mean' 
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
    
    % find mass at r200,c , define new host 
    rr=0.1:0.01:1;
    rom=meanDensity(host,rr);
    refc=200.*rho_crit(0,'cosmo',cosmoStruct)/1e9;
    r200c=interp1(rom,rr,refc);
    m200c=mass(host,r200c);
    
    % find rscale - define new ccvir
    rs=host.Rvir./host.cvir;
    cv2=r200c*host.Rvir/rs;
    
    host2=NFW('mv',m200c,'cc',cv2,'fg',fgh,'crit');
          
    %% generate velocities at r200,c
    nOrb=1500;
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
                
        % find v_tangetial
        vt=vv.*sqrt(1-vr.^2);
        
        % create vx vy, but going outwards 
        vx=vr.*vv.*host2.Vvir;  % vr is pointing outwards! 
        vy=-vt.*host2.Vvir;
        
        % set condition
        condition = ['rad(cnt) >= ' num2str(host.Rvir)];
        
        %% build orbits
        y0=0;
        x0=host2.Rvir;
        
        tOrbit=2*x0./host2.Vvir;
        tmax=2*tOrbit;
        dt=tOrbit/1e4;
        dtmin=dt/1e3;
        
        %% calculate orbits
        orbIC.x=x0;
        orbIC.y=y0;
        
        step=10;
        prc0=10;
        nOrb=length(vx);
        
        fprintf('Calculating %s *outgoing* orbits for mass ratio %s \n',num2str(nOrb),mrTag)
        for i=1:nOrb
            
            prc=i./nOrb*100;
            if prc>=prc0
                fprintf('completed %s %% of orbits \n',num2str(prc0))
                prc0=prc0+step;
            end
            
            
            orbIC.vx=vx(i);
            orbIC.vy=vy(i);
            
            orbO=orbits.rk4_orbitIntegration(orbIC,dt,tmax,dtmin,@orbits.rhs_nfw,host2,...
                'condition',condition);
            
           
            ind1=find(orbO.rad==max(orbO.rad),1,'first');
            vr=abs((orbO.vx.*orbO.x+orbO.vy.*orbO.y)./orbO.rad);
            ind2=find(vr==min(vr),1,'first');
            ind=min(ind1,ind2);
            
            newIC(i).x=orbO.x(ind);
            newIC(i).y=orbO.y(end);
            newIC(i).vx=-orbO.vx(ind);
            newIC(i).vy=-orbO.vy(ind);
%             
%             indx(i)=ind;
%             orbOut(i)=orbO;
        end
        
        
        %% generate new orbits 
        
        tOrbit=2*host2.Rvir./host2.Vvir;
        tmax=3.5*tOrbit;
        dt=tOrbit/1e4;
        dtmin=dt/1e3;
               
        
        step=10;
        prc0=10;
                
        fprintf('Calculating %s *ingoing* orbits for mass ratio %s \n',num2str(nOrb),mrTag)
        for i=1:nOrb
            
            prc=i./nOrb*100;
            if prc>=prc0
                fprintf('completed %s %% of orbits \n',num2str(prc0))
                prc0=prc0+step;
            end
            
            orb(i)=orbits.rk4_orbitIntegration(newIC(i),dt,tmax,dtmin,@orbits.rhs_nfw,host);
                
            
        end
        
        
        global DEFAULT_MATFILE_DIR
        name=sprintf('orbBank%i_mean_hostMass%s_massRatio%s.mat',nOrb,hostTags{kk},mrTag);
        save([DEFAULT_MATFILE_DIR '/' name],'orb','-v7.3')
        fprintf('saving to: %s \n',[DEFAULT_MATFILE_DIR '/' name])
        
    end
    
end
