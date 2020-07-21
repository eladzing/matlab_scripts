units

global DEFAULT_MATFILE_DIR
global DRACOFLAG
%% generate halo masses
ng=100;

massType='mean';
simEmulate='TNG100';


betaRange=[1/5 1/3 1/2 1/1.5 4/5 1.0];
fgsRange=[0.05 0.1 0.2 0.3 0.4 0.5];

betaVals=0.5; %1/2;
fgsVals=1;


%9 - 9.5
mhRange(1,:)=log10([1.26e11 2.29e11]);
%9.5 -10
mhRange(2,:)=log10([2.29e11 4.455e11]);
%10 -10.5
mhRange(3,:)=log10([4.455e11 1.07e12]);
%10.5 -11
mhRange(4,:)=log10([1.07e12 5.56e12]);
%11 - 11.5
mhRange(5,:)=log10([5.56e12 6.8e13]);
tags={'9_95','95_10','10_105','105_11','11_115'};

for kk=1:5
    fprintf('working on halo mass range %s \n',tags{kk});
    satHaloRange=mhRange(kk,:);
    
    if ~DRACOFLAG
        satHalo=generate_halo_masses(ng,'range',satHaloRange,massType);
    else
        satHalo=generate_halo_masses(ng,'range',satHaloRange,massType,'draco');
    end
    
    stellarMass=stellarmass_from_moster(satHalo);
    
    nCats=1;
    %nCats=min(nCats,length(massCatalog));
    for k=1:nCats
        fprintf('beginning iteration %i of %i \n',k,nCats);
        tic
        
        ngal=length(stellarMass);
        
        %fgs=generateGasRatio(stellarMass,simEmulate);
        fgs=fgsVals.*ones(size(stellarMass));
        
        beta=betaVals.*ones(size(stellarMass));
        
        %betaRange=[1/2.5  1/2 1/1.5  1  1.5 2  2.5];
        %fgsRange=[0.05 0.2 0.35 0.5 0.65 0.8 0.95];
        
        
        fbs=zeros(size(stellarMass));
        xi=ones(size(stellarMass));
        
        cv=cvir_Mvir(satHalo,0,'random');
        lambda=lambda_prime(satHalo);
        
        % calculate Rd
        cnt=0;
        tot=(length(betaVals)+1)*(length(fgsVals)+1);
        for i=1:length(betaVals)+1
            
            
            for j=1:length(fgsVals)+1
                cnt=cnt+1;
                
                if i~=length(betaVals)+1
                    
                    beta=ones(size(stellarMass)).*betaVals(i);
                else
                    beta=betaRange(1)+diff(betaRange([1 end])).*rand(ngal,1);
                end
                
                
                if j~=length(fgsVals)+1
                    
                    fgs=ones(size(stellarMass)).*fgsVals(j);
                else
                    fgs=fgsRange(1)+diff(fgsRange([1 end])).*rand(ngal,1);
                end
                
                
                res=rscale_mmw_array(stellarMass,'fg',fgs,'beta',beta,'fb',fbs,'xi',xi,...
                    'Mv',satHalo,'cv',cv,'lambda',lambda,'noshow');
                %      res=rscale_mmw_array(stellarMass,'fg',fgs,'beta',beta,...
                %          'Mv',satHalo,'cv',cv,'lambda',lambda,'noshow');
                %
                cata(j,i).Ms=stellarMass;
                cata(j,i).rd=res.rd;
                cata(j,i).Mv=res.Mv;
                cata(j,i).md=res.md;
                %cata(j,i).BT=res.BT;
                cata(j,i).fbs=fbs;
                cata(j,i).fgs=fgs;
                cata(j,i).xi=xi;
                cata(j,i).cv=res.cv;
                cata(j,i).lambda=res.lambda;
                %cata(j,i).sigma=stellarMass./(2*pi*res.rd.^2);
                %cata(j,i).sigmaeff=stellarMass./(2*pi*res.rd.^2).*sigma_factor('half');
                %cata(j,i).sigma5090=stellarMass./(2*pi*res.rd.^2).*sigma_factor('5090');
                cata(j,i).beta=beta;
                % cata(j,i).qmin=qmin;
                % cata(j,i).qsub=qsub;
                % cata(j,i).eps=eps;
                cata(j,i).lambdaMask=res.lambdaMask;
                
                fprintf('completed %i %% of catalogs \n',ceil(cnt/tot*100));
            end
        end
        
        fullCatalog(k).cata=cata;
        toc
    end
    
    
    name=sprintf('catalog_%s_ng100.mat',tags{kk});
    save([DEFAULT_MATFILE_DIR '/' name],'fullCatalog','-v7.3')
    fprintf('saving to: %s \n',[DEFAULT_MATFILE_DIR '/' name])
    
end