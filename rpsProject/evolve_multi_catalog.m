%global DEFAULT_MATFILE_DIR
%load([DEFAULT_MATFILE_DIR '/catalog_masses.mat'])


%% set Host parameters
Mh=5e14;
cvh=cvir_Mvir(Mh,0,'random');
fgh=0.15;
host=NFW('mv',Mh,'cc',cvh,'fg',fgh,'mean');



%% set orbit parameters
if orbitFlag
    y0=0;
    x0=1.5.*host.Rvir;
    vTotat=host.Vvir; %.*rand(ngal,1);
    %vAngle=pi/2.*rand(ngal,1); % in radians
    vx0= -1.*vTotat.*cos(vAngle);
    vy0=vTotat.*sin(vAngle);
    
    tOrbit=x0./host.Vvir;
    tmax=2*tOrbit;
    dt=tOrbit/1e4;
    dtmin=dt/1e3;
    
    %orbits
    orbIC.x=x0;
    orbIC.y=y0;
    orbIC.vx=vx0;%(i);
    orbIC.vy=vy0;%(i);
    fprintf('Orbit Initial conditions: \n')
    orbIC 
    
    orb=orbits.rk4_orbitIntegration(orbIC,dt,tmax,dtmin,@orbits.rhs_nfw,host);
end

fprintf('Running evolution \n')


sz=size(cata);
betaRange=[1/2.5  1/2 1/1.5  1  1.5 2  2.5 Inf];
fgsRange=[0.05 0.2 0.35 0.5 0.65 0.8 0.95 Inf];


for k=1:sz(1)
    for j=1:sz(2)
        
        catg=cata(k,j);
        ngal=length(catg.Ms);
        fprintf('Running on catalog %i of %i, evoloving galaxies \n',...
            (k-1)*sz(2)+j,prod(sz))
        prc0=20;
        step=20;
        for i=1:ngal
            %gals
            
            prc=i/ngal*100;
            if prc>=prc0
                
                fprintf('completed %s %% of catalog \n',num2str(prc0))
                
                prc0=prc0+step;
                
            end
            gal=GALAXY('ms',catg.Ms(i),'rd',catg.rd(i),'fgs',catg.fgs(i),...
                'beta',catg.beta(i),'fbs',catg.fbs(i),'xi',catg.xi(i),...
                'Mh',catg.Mv(i),'cv',catg.cv(i));
            
            
            
            
            galRes(i)=galEvolutionMachineLight(gal,host,orb,'nobar');
            
        end
        catEvol(k,j).galRes=galRes;
        
        
        %% calculate quenched fraction
        
        iend=length(galRes(1).rpos);
        %iend=6599;
        edj=0:0.1:1.2;
        bCen=edj(1:end-1)+0.5.*diff(edj);
        
        qf1=zeros(1,length(bCen));
        qf2=zeros(1,length(bCen));
        cnt1=qf1;
        cnt2=qf2;
        for ii=1:length(galRes)
            rpp=galRes(ii).rpos(1:iend)./host.Rvir;
            rpp2=rpp.*generate_projectionFac(length(rpp))';
            bInd=discretize(rpp,edj);
            bInd2=discretize(rpp2,edj);
            
            
            qn=galRes(ii).ssfr(1:iend)<1e-11;
            onn=ones(size(qn));
            
            qff1=zeros(1,length(bCen));
            cntt1=qff1;
            qff2=qff1;
            cntt2=qff1;
            
            for jj=1:length(bCen)
                qff1(jj)=sum(qn(bInd==jj));
                cntt1(jj)=sum(onn(bInd==jj));
                
                qff2(jj)=sum(qn(bInd2==jj));
                cntt2(jj)=sum(onn(bInd2==jj));
            end
            
            qf1=qf1+qff1;
            qf2=qf2+qff2;
            cnt1=cnt1+cntt1;
            cnt2=cnt2+cntt2;
            
        end
        
        qFrac(k,j).qf=qf1;
        qFrac(k,j).qfProj=qf2;
        qFrac(k,j).cnt=cnt1;
        qFrac(k,j).cntProj=cnt2;
        qFrac(k,j).bCen=bCen;
        qFrac(k,j).fgs=fgsRange(k);
        qFrac(k,j).beta=betaRange(k);
        
        
        
        
        
    end
end













