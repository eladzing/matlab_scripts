global DEFAULT_MATFILE_DIR
%load([DEFAULT_MATFILE_DIR '/catalog_masses.mat'])

mh=       [12.22 12.83 13.23 13.83 14.19 14.81 15.12];
fghTNG100=[0.08 0.05  0.07   0.1    0.1  0.12    0.13];   
fghTNG300=[0.1  0.05  0.07   0.1    0.12  0.13   0.13 ];   

kstart=[1     1     1     1     1     1     1];
kend=  [2     4     4     5     5     5     5];

% full set: 
%kstart=[1     1     1     1     1     1     1];
%kend=  [2     4     4     5     5     5     5];



%kstart=[3  2    2     1     1     1     1];
hostTags={'1e12' '5e12' ...
    '1e13' '5e13' ...
    '1e14' '5e14' '1e15'};
tags={'9_95','95_10','10_105','105_11','11_115'};

if ~exist('simEmulate','var')
    error('please define a simulation to draw fgs from')
end


for hInd=1:length(mh)
    %% set Host parameters
    
    fprintf('Evolving galaxies in host of log(M)=%f \n',mh(hInd))
    
    mhost=10.^mh(hInd);
    cvh=cvir_Mvir(mhost,0);
    fgh=fghTNG300(hInd);
    host=NFW('mv',mhost,'cc',cvh,'fg',fgh,'mean');
    
    for kk=kstart(hInd):kend(hInd)
        
        %% load catalog
        
        %name=sprintf('catalog_%s_ng100_fgs%s.mat',tags{kk},simEmulate);
        name=sprintf('catalog_fg02_bet05_%s_ng100.mat',tags{kk});
        load([DEFAULT_MATFILE_DIR '/' name],'fullCatalog')
        fprintf('loaded: %s \n', name)
        
        cata=fullCatalog.cata; %(1,1);
        
       
        
        
        %% orbits
        
        fprintf('Arranging orbits evolution \n')
        
        satHalo=cata(1,1).Mv;
        
        massRat=satHalo./host.Mvir;
        
        galMask=massRat<0.25;
        galInds=find(galMask);
        
        mr=[0 0.005 0.05 0.25];
        mrTag={'null' '0005' '005' '05'};
        orbs=struct([]);
        inds=[];
        tic
        for k=2:length(mr)
            fprintf('Mass ratio %s \n',mrTag{k})
            mask=massRat>=mr(k-1) & massRat<mr(k);
            inds0=find(mask);
            
            if ~any(mask)
                continue
            end
            
            fname=sprintf('orbBank1500_mean_hostMass%s_massRatio%s.mat',...
                hostTags{hInd},mrTag{k});
            load([DEFAULT_MATFILE_DIR '/' fname]);
            
            
            oInd=1:length(orb);
            oInd=shuffleArray(oInd);
            
            oIndF=oInd(1:sum(mask));
            
            orbs0=orb(oIndF);
            clear orb
            
            orbs=cat(2,orbs,orbs0);
            inds=cat(1,inds,inds0);
        end
        [iii,sInd]=sort(inds);
        if ~all(iii==galInds)
            error('EVOLVE_CATALOG - index mismatch')
        end
        
        
        orbs=orbs(sInd);
        toc
        
        fprintf('Running evolution \n')
        
        ngal=length(satHalo); %length(galInds); %length(catg.Ms);
             
        %sz=size(cata);
        %sz=[1 1];
        %betaRange=[1/2.5  1/2 1/1.5  1  1.5 2  2.5 Inf];
        %fgsRange=[0.05 0.2 0.35 0.5 0.65 0.8 0.95 Inf];
        
        %betaRange=[1/2 Inf];
        %fgsRange=[0.2 Inf];
        
        %betaRange=[1/5 1/3 1/2 1/1.5 4/5 1.0 Inf];
       betaRange=0.5; 
        nCat=1;
        nBet=1; %[3 7];
        catCount=0;
        for k=nCat
            for j=1:length(nBet) %   1:sz(2)  %go over beta range 
		    catCount=catCount+1;
                catg=cata(k,nBet(j));
                  fprintf('Running on catalog %i of %i, evoloving galaxies \n',...
                    catCount,length(nCat).*length(nBet))
%                 fprintf('Running on catalog %i of %i, evoloving galaxies \n',...
%                     (k-1)*sz(2)+j,prod(sz))
                prc0=20;
                step=20;
                
                galRes(1,ngal)=struct('time',0,'rpos',0,'pram',0,'rhoICM',0,'sfr',0,'ssfr',0,'gasMass',0,'stellarMass',0);
                
                iend1Pass=zeros(1,ngal);
                    iendRv1=zeros(1,ngal);
                    iendApo=zeros(1,ngal);
                    orbInd=0;
                for i=1:ngal %running only over galaxies within mask 
                              
                    prc=i/ngal*100;
                    if prc>=prc0
                        
                        fprintf('completed %s %% of catalog \n',num2str(prc0))
                        
                        prc0=prc0+step;
                        
                    end
                                      
                    if galMask(i) % only run on galaxies of small enough mass ratios
                        
                        gal=GALAXY('ms',catg.Ms(i),'rd',catg.rd(i),'fgs',catg.fgs(i),...
                            'beta',catg.beta(i),'fbs',0,'xi',1,...
                            'Mh',catg.Mv(i),'cv',catg.cv(i));
                        
                        orbInd=orbInd+1;
                        orb=orbs(orbInd);
                        galRes(i)=galEvolutionMachineLight(gal,host,orb,'nobar','rpsAlfa',1.0);
                        
                        % find final index for orbit
                        iend1=find(orb.rad>host.Rvir,1,'first');
                        vra=orb.x.*orb.vx+orb.y.*orb.vy;
                        iend2=find(vra(1:end-1)>0 & vra(2:end)<0,1,'first');
                        
                        if isempty(iend1)
                            iend1=length(orb.rad);
                        end
                        
                        if isempty(iend2)
                            iend2=length(orb.rad);
                        end
                        
                        iend1Pass(i)=min(iend1,iend2);
                        iendRv1(i)=iend1;
                        iendApo(i)=iend2;
                        
                    end
                end
                catEvol(k,j).galRes=galRes;
                
                
                %% calculate quenched fraction by position
                
                % do 1 pass;
                
                edj=0:0.1:1;
                bCen=edj(1:end-1)+0.5.*diff(edj);
                
                qf1=zeros(1,length(bCen));
                qf2=zeros(1,length(bCen));
                cnt1=qf1;
                cnt2=qf2;
                for ii=1:ngal %length(galRes)
                    
                     if ~galMask(ii)
                        continue
                     end
                    
                    rpp=galRes(ii).rpos(1:iend1Pass(ii))./host.Rvir;
                    rpp2=rpp.*generate_projectionFac(length(rpp));
                    bInd=discretize(rpp,edj);
                    bInd2=discretize(rpp2,edj);
                    
                    
                    qn=galRes(ii).ssfr(1:iend1Pass(ii))<1e-11;
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
                
                qFrac.r1Pass(k,j).qf=qf1;
                qFrac.r1Pass(k,j).qfProj=qf2;
                qFrac.r1Pass(k,j).cnt=cnt1;
                qFrac.r1Pass(k,j).cntProj=cnt2;
                qFrac.r1Pass(k,j).bCen=bCen;
                qFrac.r1Pass(k,j).fgs=simEmulate;
                qFrac.r1Pass(k,j).beta=betaRange(j);
                
                
                % do till end pass;
                
                edj=0:0.1:1;
                bCen=edj(1:end-1)+0.5.*diff(edj);
                
                qf1=zeros(1,length(bCen));
                qf2=zeros(1,length(bCen));
                cnt1=qf1;
                cnt2=qf2;
                for ii=1:ngal %length(galRes)
                    
                     if ~galMask(ii)
                        continue
                    end
                    
                    rpp=galRes(ii).rpos./host.Rvir;
                    rpp2=rpp.*generate_projectionFac(length(rpp));
                    bInd=discretize(rpp,edj);
                    bInd2=discretize(rpp2,edj);
                    
                    
                    qn=galRes(ii).ssfr<1e-11;
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
                
                qFrac.rFull(k,j).qf=qf1;
                qFrac.rFull(k,j).qfProj=qf2;
                qFrac.rFull(k,j).cnt=cnt1;
                qFrac.rFull(k,j).cntProj=cnt2;
                qFrac.rFull(k,j).bCen=bCen;
                qFrac.rFull(k,j).fgs=simEmulate;
                qFrac.rFull(k,j).beta=betaRange(j);
                
                
                
                %% qf by stellar mass
                
                edj=9:0.5:11.5;
                
                bCen=edj(1:end-1)+0.5.*diff(edj);
                
                qf1=zeros(1,length(bCen));
                qf2=zeros(1,length(bCen));
                cnt1=qf1;
                cnt2=qf2;
                
                qfP1=zeros(1,length(bCen));
                qfP2=zeros(1,length(bCen));
                cntP1=qfP1;
                cntP2=qfP2;
                
                for ii=1:ngal %length(galRes)
                    
                     if ~galMask(ii)
                        continue
                     end
                    
                    % only count when galaxy is in 1Rvir
                    rp=galRes(ii).rpos./host.Rvir;
                    rpP=rp.*generate_projectionFac(length(rp));
                                          
                    rmask=rp<1;
                    rmaskP=rpP<1;
                                        
                    rmask1=rmask;
                    rmask1(iend1Pass(ii)+1:end)=false;
                    
                    rmaskP1=rmaskP;
                    rmaskP1(iend1Pass(ii)+1:end)=false;
                     
                    sMass=log10(galRes(ii).stellarMass);
                    sMass(sMass<9)=9.01;
                    
                    %             rpp=galRes(ii).rpos./host.Rvir;
                    %             rpp2=rpp.*generate_projectionFac(length(rpp));
                    bInd1=discretize(sMass(rmask),edj);
                    bInd2=discretize(sMass(rmask1),edj);
                    
                    bIndP1=discretize(sMass(rmaskP),edj);
                    bIndP2=discretize(sMass(rmaskP1),edj);
                                     
                    
                    qn1=galRes(ii).ssfr(rmask)<1e-11;
                    onn1=ones(size(qn1));
                    qn2=galRes(ii).ssfr(rmask1)<1e-11;
                    onn2=ones(size(qn2));
                    
                    qnP1=galRes(ii).ssfr(rmaskP)<1e-11;
                    onnP1=ones(size(qnP1));
                    qnP2=galRes(ii).ssfr(rmaskP1)<1e-11;
                    onnP2=ones(size(qnP2));
                    
                    
                    
                    qff1=zeros(1,length(bCen));
                    cntt1=qff1;
                    qff2=qff1;
                    cntt2=qff1;
                    
                    qffP1=zeros(1,length(bCen));
                    cnttP1=qffP1;
                    qffP2=qffP1;
                    cnttP2=qffP1;
                    
                    
                    
                    for jj=1:length(bCen)
                        qff1(jj)=sum(qn1(bInd1==jj));
                        cntt1(jj)=sum(onn1(bInd1==jj));
                        
                        qff2(jj)=sum(qn2(bInd2==jj));
                        cntt2(jj)=sum(onn2(bInd2==jj));
                        
                        % projected
                        qffP1(jj)=sum(qnP1(bIndP1==jj));
                        cnttP1(jj)=sum(onnP1(bIndP1==jj));
                        
                        qffP2(jj)=sum(qnP2(bIndP2==jj));
                        cnttP2(jj)=sum(onnP2(bIndP2==jj));
                    end
                    
                    qf1=qf1+qff1;
                    qf2=qf2+qff2;
                    cnt1=cnt1+cntt1;
                    cnt2=cnt2+cntt2;
                    
                    qfP1=qfP1+qffP1;
                    qfP2=qfP2+qffP2;
                    cntP1=cntP1+cnttP1;
                    cntP2=cntP2+cnttP2;
                    
                end
                
                qFrac.msFull(k,j).qf=qf1;
                qFrac.msFull(k,j).cnt=cnt1;
                qFrac.msFull(k,j).qfProj=qfP1;
                qFrac.msFull(k,j).cntProj=cntP1;
                                
                qFrac.msFull(k,j).bCen=bCen;
                qFrac.msFull(k,j).fgs=simEmulate;
                qFrac.msFull(k,j).beta=betaRange(j);
                
                
                qFrac.ms1Pass(k,j).qf=qf2;
                qFrac.ms1Pass(k,j).cnt=cnt2;
                qFrac.ms1Pass(k,j).qfProj=qfP2;
                qFrac.ms1Pass(k,j).cntProj=cntP2;
                
                qFrac.ms1Pass(k,j).bCen=bCen;
                qFrac.ms1Pass(k,j).fgs=simEmulate;
                qFrac.ms1Pass(k,j).beta=betaRange(j);
                
                
                
                %% %% qf by Mvsat/Mvhost
                
                edj=[1e-4 0.0005 0.001 0.005 0.01 0.05 0.1 0.25];
                
                bCen=edj(1:end-1)+0.5.*diff(edj);
                
                qf1=zeros(1,length(bCen));
                qf2=zeros(1,length(bCen));
                cnt1=qf1;
                cnt2=qf2;
                
                qfP1=zeros(1,length(bCen));
                qfP2=zeros(1,length(bCen));
                cntP1=qfP1;
                cntP2=qfP2;
                
                massR=massRat;
                massR(massR<edj(1))=edj(1);
                
                for ii=1:ngal %length(galRes)
                    
                     if ~galMask(ii)
                        continue
                     end
                    
                    % only count when galaxy is in 1Rvir
                    rp=galRes(ii).rpos./host.Rvir;
                    rpP=rp.*generate_projectionFac(length(rp));
                                          
                    rmask=rp<1;
                    rmaskP=rpP<1;
                                        
                    rmask1=rmask;
                    rmask1(iend1Pass(ii)+1:end)=false;
                    
                    rmaskP1=rmaskP;
                    rmaskP1(iend1Pass(ii)+1:end)=false;
                     
                    %massRat=catg(ii).Mv./host.Mvir;
                    %massRat(massRat<edj(1))=edj(1);
                    mr=massRat.*ones(size(galRes(ii).stellarMass));
                    
                    %             rpp=galRes(ii).rpos./host.Rvir;
                    %             rpp2=rpp.*generate_projectionFac(length(rpp));
                    bInd1=discretize(mr(rmask),edj);
                    bInd2=discretize(mr(rmask1),edj);
                    
                    bIndP1=discretize(mr(rmaskP),edj);
                    bIndP2=discretize(mr(rmaskP1),edj);
                    
                    qn1=galRes(ii).ssfr(rmask)<1e-11;
                    onn1=ones(size(qn1));
                    qn2=galRes(ii).ssfr(rmask1)<1e-11;
                    onn2=ones(size(qn2));
                    
                     qnP1=galRes(ii).ssfr(rmaskP)<1e-11;
                    onnP1=ones(size(qnP1));
                    qnP2=galRes(ii).ssfr(rmaskP1)<1e-11;
                    onnP2=ones(size(qnP2));
                    
                    qff1=zeros(1,length(bCen));
                    cntt1=qff1;
                    qff2=qff1;
                    cntt2=qff1;
                    
                    qffP1=zeros(1,length(bCen));
                    cnttP1=qffP1;
                    qffP2=qffP1;
                    cnttP2=qffP1;
                    
                    
                    
                    for jj=1:length(bCen)
                        qff1(jj)=sum(qn1(bInd1==jj));
                        cntt1(jj)=sum(onn1(bInd1==jj));
                        
                        qff2(jj)=sum(qn2(bInd2==jj));
                        cntt2(jj)=sum(onn2(bInd2==jj));
                        
                        % projected
                        qffP1(jj)=sum(qnP1(bIndP1==jj));
                        cnttP1(jj)=sum(onnP1(bIndP1==jj));
                        
                        qffP2(jj)=sum(qnP2(bIndP2==jj));
                        cnttP2(jj)=sum(onnP2(bIndP2==jj));
                    end
                    
                    qf1=qf1+qff1;
                    qf2=qf2+qff2;
                    cnt1=cnt1+cntt1;
                    cnt2=cnt2+cntt2;
                    
                    qfP1=qfP1+qffP1;
                    qfP2=qfP2+qffP2;
                    cntP1=cntP1+cnttP1;
                    cntP2=cntP2+cnttP2;
                    
                end
                
                qFrac.massRatFull(k,j).qf=qf1;
                qFrac.massRatFull(k,j).cnt=cnt1;
                qFrac.massRatFull(k,j).qfProj=qfP1;
                qFrac.massRatFull(k,j).cntProj=cntP1;
                                
                qFrac.massRatFull(k,j).bCen=bCen;
                qFrac.massRatFull(k,j).fgs=simEmulate;
                qFrac.massRatFull(k,j).beta=betaRange(j);
                
                
                qFrac.massRat1Pass(k,j).qf=qf2;
                qFrac.massRat1Pass(k,j).cnt=cnt2;
                qFrac.massRat1Pass(k,j).qfProj=qfP2;
                qFrac.massRat1Pass(k,j).cntProj=cntP2;
                
                qFrac.massRat1Pass(k,j).bCen=bCen;
                qFrac.massRat1Pass(k,j).fgs=simEmulate;
                qFrac.massRat1Pass(k,j).beta=betaRange(j);
                
                
                 %% %% qf by Mstellar/Mvhost
                
                edj=10.^(-7:1:-1);
                
                bCen=edj(1:end-1)+0.5.*diff(edj);
                
                qf1=zeros(1,length(bCen));
                qf2=zeros(1,length(bCen));
                cnt1=qf1;
                cnt2=qf2;
                
                qfP1=zeros(1,length(bCen));
                qfP2=zeros(1,length(bCen));
                cntP1=qfP1;
                cntP2=qfP2;
                
                %massR=massRat;
                %massR(massR<edj(1))=edj(1);
                
                for ii=1:ngal %length(galRes)
                    
                     if ~galMask(ii)
                        continue
                     end
                    
                    % only count when galaxy is in 1Rvir
                    rp=galRes(ii).rpos./host.Rvir;
                    rpP=rp.*generate_projectionFac(length(rp));
                                          
                    rmask=rp<1;
                    rmaskP=rpP<1;
                                        
                    rmask1=rmask;
                    rmask1(iend1Pass(ii)+1:end)=false;
                    
                    rmaskP1=rmaskP;
                    rmaskP1(iend1Pass(ii)+1:end)=false;
                     
                    mr=galRes(ii).stellarMass./host.Mvir;
                    
                    %             rpp=galRes(ii).rpos./host.Rvir;
                    %             rpp2=rpp.*generate_projectionFac(length(rpp));
                    bInd1=discretize(mr(rmask),edj);
                    bInd2=discretize(mr(rmask1),edj);
                    
                    bIndP1=discretize(mr(rmaskP),edj);
                    bIndP2=discretize(mr(rmaskP1),edj);
                    
                    qn1=galRes(ii).ssfr(rmask)<1e-11;
                    onn1=ones(size(qn1));
                    qn2=galRes(ii).ssfr(rmask1)<1e-11;
                    onn2=ones(size(qn2));
                    
                    qnP1=galRes(ii).ssfr(rmaskP)<1e-11;
                    onnP1=ones(size(qnP1));
                    qnP2=galRes(ii).ssfr(rmaskP1)<1e-11;
                    onnP2=ones(size(qnP2));
                    
                    qff1=zeros(1,length(bCen));
                    cntt1=qff1;
                    qff2=qff1;
                    cntt2=qff1;
                    
                    qffP1=zeros(1,length(bCen));
                    cnttP1=qffP1;
                    qffP2=qffP1;
                    cnttP2=qffP1;
                    
                    
                    
                    for jj=1:length(bCen)
                        qff1(jj)=sum(qn1(bInd1==jj));
                        cntt1(jj)=sum(onn1(bInd1==jj));
                        
                        qff2(jj)=sum(qn2(bInd2==jj));
                        cntt2(jj)=sum(onn2(bInd2==jj));
                        
                        % projected
                        qffP1(jj)=sum(qnP1(bIndP1==jj));
                        cnttP1(jj)=sum(onnP1(bIndP1==jj));
                        
                        qffP2(jj)=sum(qnP2(bIndP2==jj));
                        cnttP2(jj)=sum(onnP2(bIndP2==jj));
                    end
                    
                    qf1=qf1+qff1;
                    qf2=qf2+qff2;
                    cnt1=cnt1+cntt1;
                    cnt2=cnt2+cntt2;
                    
                    qfP1=qfP1+qffP1;
                    qfP2=qfP2+qffP2;
                    cntP1=cntP1+cnttP1;
                    cntP2=cntP2+cnttP2;
                    
                end
                
                qFrac.massRatSFull(k,j).qf=qf1;
                qFrac.massRatSFull(k,j).cnt=cnt1;
                qFrac.massRatSFull(k,j).qfProj=qfP1;
                qFrac.massRatSFull(k,j).cntProj=cntP1;
                                
                qFrac.massRatSFull(k,j).bCen=bCen;
                qFrac.massRatSFull(k,j).fgs=simEmulate;
                qFrac.massRatSFull(k,j).beta=betaRange(j);
                
                
                qFrac.massRatS1Pass(k,j).qf=qf2;
                qFrac.massRatS1Pass(k,j).cnt=cnt2;
                qFrac.massRatS1Pass(k,j).qfProj=qfP2;
                qFrac.massRatS1Pass(k,j).cntProj=cntP2;
                
                qFrac.massRatS1Pass(k,j).bCen=bCen;
                qFrac.massRatS1Pass(k,j).fgs=simEmulate;
                qFrac.massRatS1Pass(k,j).beta=betaRange(j);
                
                
            end
        end
        
        name=sprintf('catalogEvolved_Ms100_%s_host_%s_fgh300_fgs%s.mat',tags{kk},hostTags{hInd},simEmulate) ;
        save([DEFAULT_MATFILE_DIR '/' name],...
            'catEvol','galMask','iend1Pass','iendRv1','iendApo','-v7.3')
        fprintf('saving galRes to: %s \n',[DEFAULT_MATFILE_DIR '/' name])
        
        name=sprintf('qFrac_Ms100_%s_host%s_fgh300_fgs%s.mat',tags{kk},hostTags{hInd},simEmulate);
        save([DEFAULT_MATFILE_DIR '/' name],'qFrac','-v7.3')
        fprintf('saving qfrac to: %s \n',[DEFAULT_MATFILE_DIR '/' name])
        
    end
end







