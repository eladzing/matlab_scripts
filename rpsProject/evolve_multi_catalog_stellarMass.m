global DEFAULT_MATFILE_DIR
%load([DEFAULT_MATFILE_DIR '/catalog_masses.mat'])

mh=[12.22 12.83 13.23 13.83 14.19 14.81 15.12];
kend=[3  2    2     1     1     1     1];
%kstart=[3  2    2     1     1     1     1];
hostTags={'1e12' '5e12' ...
    '1e13' '5e13' ...
    '1e14' '5e14' '1e15'};
tags={'9_95','95_10','10_105','105_11','11_115'};

for jj=1:length(mh)
%% set Host parameters

fprintf('Evolving galaxies in host of log(M)=%f \n',mh(jj))
    
mhost=10.^mh(kk);
cvh=cvir_Mvir(mhost,0);
fgh=0.15;
host=NFW('mv',mhost,'cc',cvh,'fg',fgh,'crit');



for kk=
    
    %% load catalog
    
    name=sprintf('catalog_%s_ng100.mat',tags{kk});
    load([DEFAULT_MATFILE_DIR '/' name],'fullCatalog')
    fprintf('loaded: %s \n', name)
    
    cata=fullCatalog.cata;
    
    
    %% orbits
    if orbitFlag
        fprintf('Arranging orbits evolution \n')
        satHalo=cata.Mv;
        
        massRat=cata(1,1).Mv./host.Mvir;
        
        mr=[0 0.005 0.05 0.5];
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
            
            fname=sprintf('orbBank_hostMass142_massRatio%s.mat',mrTag{k});
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
        orbs=orbs(sInd);
        toc
    end
    fprintf('Running evolution \n')
 
    
    %sz=size(cata);
    sz=[1 1];
    %betaRange=[1/2.5  1/2 1/1.5  1  1.5 2  2.5 Inf];
    %fgsRange=[0.05 0.2 0.35 0.5 0.65 0.8 0.95 Inf];
    
    betaRange=[1/2 Inf];
    fgsRange=[0.2 Inf];
    
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
                    'beta',catg.beta(i),'fbs',0,'xi',1,...
                    'Mh',catg.Mv(i),'cv',catg.cv(i));
                
                
                %orbInd=oInd(i);
                
                galRes(i)=galEvolutionMachineLight(gal,host,orbs(i),'nobar');
                
                
                % find final index for orbit
                iend1=find(orbs(i).rad>host.Rvir,1,'first');
                vra=orbs(i).x.*orbs(i).vx+orbs(i).y.*orbs(i).vy;
                iend2=find(vra(1:end-1)>0 & vra(2:end)<0,1,'first');
                
                if isempty(iend1)
                    iend1=length(orbs(i).rad);
                end
                
                if isempty(iend2)
                    iend2=length(orbs(i).rad);
                end
                
                iend1Pass(i)=min(iend1,iend2);
                iendRv1(i)=iend1;
                iendApo(i)=iend2;
                
                
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
            for ii=1:length(galRes)
                
                
                rpp=galRes(ii).rpos(1:iend1Pass(ii))./host.Rvir;
                rpp2=rpp.*generate_projectionFac(length(rpp))';
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
            qFrac.r1Pass(k,j).fgs=fgsRange(k);
            qFrac.r1Pass(k,j).beta=betaRange(j);
            
            
            % do till end pass;
            
            edj=0:0.1:1;
            bCen=edj(1:end-1)+0.5.*diff(edj);
            
            qf1=zeros(1,length(bCen));
            qf2=zeros(1,length(bCen));
            cnt1=qf1;
            cnt2=qf2;
            for ii=1:length(galRes)
                
                
                rpp=galRes(ii).rpos./host.Rvir;
                rpp2=rpp.*generate_projectionFac(length(rpp))';
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
            qFrac.rFull(k,j).fgs=fgsRange(k);
            qFrac.rFull(k,j).beta=betaRange(j);
            
            
            
            %% qf by stellar mass
            
            edj=9:0.5:11.5;
            
            bCen=edj(1:end-1)+0.5.*diff(edj);
            
            qf1=zeros(1,length(bCen));
            qf2=zeros(1,length(bCen));
            cnt1=qf1;
            cnt2=qf2;
            
            for ii=1:length(galRes)
                
                
                sMass=log10(galRes(ii).stellarMass);
                sMass(sMass<9)=9.01;
                
                %             rpp=galRes(ii).rpos./host.Rvir;
                %             rpp2=rpp.*generate_projectionFac(length(rpp))';
                bInd1=discretize(sMass,edj);
                bInd2=discretize(sMass(1:iend1Pass(ii)),edj);
                
                %bInd2=discretize(rpp2,edj);
                
                qn1=galRes(ii).ssfr<1e-11;
                onn1=ones(size(qn1));
                qn2=galRes(ii).ssfr(1:iend1Pass(ii))<1e-11;
                onn2=ones(size(qn2));
                
                qff1=zeros(1,length(bCen));
                cntt1=qff1;
                qff2=qff1;
                cntt2=qff1;
                
                for jj=1:length(bCen)
                    qff1(jj)=sum(qn1(bInd1==jj));
                    cntt1(jj)=sum(onn1(bInd1==jj));
                    
                    qff2(jj)=sum(qn2(bInd2==jj));
                    cntt2(jj)=sum(onn2(bInd2==jj));
                end
                
                qf1=qf1+qff1;
                qf2=qf2+qff2;
                cnt1=cnt1+cntt1;
                cnt2=cnt2+cntt2;
                
            end
            
            qFrac.msFull(k,j).qf=qf1;
            qFrac.msFull(k,j).cnt=cnt1;
            qFrac.msFull(k,j).bCen=bCen;
            qFrac.msFull(k,j).fgs=fgsRange(k);
            qFrac.msFull(k,j).beta=betaRange(j);
            
            qFrac.ms1Pass(k,j).qf=qf2;
            qFrac.ms1Pass(k,j).cnt=cnt2;
            qFrac.ms1Pass(k,j).bCen=bCen;
            qFrac.ms1Pass(k,j).fgs=fgsRange(k);
            qFrac.ms1Pass(k,j).beta=betaRange(j);
            
            
            
            
            
        end
    end
    
    name=sprintf('catalogEvolved_Ms100_%s.mat',tags{kk}) ;
    save([DEFAULT_MATFILE_DIR '/' name],...
        'catEvol','iend1Pass','iendRv1','iendApo','-v7.3')
    fprintf('saving galRes to: %s \n',[DEFAULT_MATFILE_DIR '/' name])
    
    name=sprintf('qFrac_Ms100_%s.mat',tags{kk});
    save([DEFAULT_MATFILE_DIR '/' name],'qFrac','-v7.3')
    fprintf('saving qfrac to: %s \n',[DEFAULT_MATFILE_DIR '/' name])
    
    
    
    
    
end








