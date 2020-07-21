
mh=[12.22 12.83 13.23 13.83 14.19 14.81 15.12];
kstart=[3  2    2     1     1     1     1];
hostTags={'1e12' '5e12' ...
    '1e13' '5e13' ...
    '1e14' '5e14' '1e15'};

units;
global DEFAULT_MATFILE_DIR
global cosmoStruct

for kk=1:length(mh)
    
    %fprintf('Generating orbits for host of log(M)=%f \n',mh(kk))
    
    mhost=10.^mh(kk);
    cvh=cvir_Mvir(mhost,0);
    fgh=0.15;
    host=NFW('mv',mhost,'cc',cvh,'fg',fgh,'mean');
    hostt(kk)=host;
    %%  Find, for given host what r_200,c and M_200,c is
    rr=0.5:0.001:1;
    roref=200*rho_crit(0,'cosmo',cosmoStruct)/1e9;
    rom=meanDensity(host,rr);
    
    r200c(kk)=interp1(rom,rr,roref).*host.Rvir;
    m200c(kk)=mass(host,r200c(kk)./host.Rvir);
    v200c(kk)=sqrt(Units.G.*(m200c(kk).*Units.Ms)./(r200c(kk).*Units.kpc))./Units.km; %in km/sec
    
    mhFld=['mh' hostTags{kk}];
    orbParam.(mhFld).host=host;
    orbParam.(mhFld).r200c=r200c(kk);
    orbParam.(mhFld).m200c=m200c(kk);
    orbParam.(mhFld).v200c=v200c(kk);
    
    
    
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
        
        %% load orbits.
        name=sprintf('orbBank1500_mean_hostMass%s_massRatio%s.mat',hostTags{kk},mrTag);
        load([DEFAULT_MATFILE_DIR '/' name],'orb')
        fprintf('loading from: %s \n',name)
        
        nOrb=length(orb);
        vv=zeros(1,nOrb);
        vr=vv;
        for i=1:nOrb
            
            ind=find(orb(i).rad<r200c(kk),1,'first');
            
            if ~isempty(ind)
                
                vorb=orb(i).vel(ind);
                vv(i)=vorb./v200c(kk);
                vrOrb=(orb(i).vx(ind)*orb(i).x(ind)+...
                    orb(i).vy(ind)*orb(i).y(ind))/orb(i).rad(ind);
                vr(i)=vrOrb./vorb;
                
            else
                vv(i)= NaN;
                vr(i)=NaN;
            end
            
        end
        
        mrFld=['mr' mrTag];
        orbParam.(mhFld).(mrFld).vv=vv;
        orbParam.(mhFld).(mrFld).vr=vr;
    end
    
end

global DEFAULT_MATFILE_DIR
name='orbParam_Rmean.mat';
save([DEFAULT_MATFILE_DIR '/' name],'orbParam','-v7.3')
fprintf('saving to: %s \n',[DEFAULT_MATFILE_DIR '/' name])














