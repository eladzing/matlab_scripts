


%% set Host parameters
Mh=1e14;
cvh=cvir_Mvir(Mh,0,'random');
fgh=0.15;
host=NFW('mv',Mh,'cc',cvh,'fg',fgh,'mean');



%% set orbit parameters
y0=0;
x0=1.5.*host.Rvir;
vTotat=host.Vvir; %.*rand(ngal,1);
%vAngle=pi/2.*rand(ngal,1); % in radians
vx0= -1.*vTotat; %.*cos(vAngle);
vy0=0; %vTotat.*sin(vAngle);

tOrbit=x0./host.Vvir;
tmax=1.3*tOrbit;
dt=tOrbit/1e4;
dtmin=dt/1e3;

%orbits
    orbIC.x=x0;
    orbIC.y=y0;
    orbIC.vx=vx0;%(i);
    orbIC.vy=vy0;%(i);

    orb=orbits.rk4_orbitIntegration(orbIC,dt,tmax,dtmin,@orbits.rhs_nfw,host);




prc0=10;
step=10;
fprintf('Running evolution \n')
tic
for i=1:ngal
    %gals
    
    prc=i/ngal*100;
    if prc>=prc0
       
        fprintf('completed %s %% of catalog - ',num2str(prc0))
        toc
        prc0=prc0+step;
        
    end
    gal(i)=GALAXY('ms',cata.Ms(i),'rd',cata.rd(i),'fgs',cata.fgs(i),...
        'beta',cata.beta(i),'fbs',cata.fbs(i),'xi',cata.xi(i),...
        'Mh',cata.Mv(i),'cv',cata.cv(i));
    
    
    
    
    galRes(i)=galEvolutionMachineLight(gal(i),host,orb,'nobar');
    
end











