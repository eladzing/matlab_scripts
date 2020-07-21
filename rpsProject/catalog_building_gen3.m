%% This builds galaxy catalogs:
% first choose halo mass range which corresponds to stellar mass range
% stellar mass set by halo mass
% fgs set by stellar mass (using relation from TNG100/TNG300 set by
% simEmulate
% generate several catalogs for varying beta values

units

global DEFAULT_MATFILE_DIR
global DRACOFLAG
%% generate halo masses
ng=100;
%snap=99;
massType='mean';
%simEmulate='TNG100';

if exist('simEmulate','var')
   fprintf('fgs Emulating by %s \n',simEmulate)
else 
    error('CATALOG_BUILDING_GEN3 - simEmulate not defined. ');
end

%betaRange=[1/5 1/3 1/2 1/1.5 4/5 1.0];
betaRange=1/2;
%fgsRange=[0.05 0.1 0.2 0.3 0.4 0.5];

%betaVals=1/2;
%fgsVals=0.2;


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
        
        fgs=generateGasRatio(stellarMass,simEmulate,'snap',snap);
        
        
        beta=ones(size(stellarMass));
        
        %betaRange=[1/2.5  1/2 1/1.5  1  1.5 2  2.5];
        %fgsRange=[0.05 0.2 0.35 0.5 0.65 0.8 0.95];
        
        
        fbs=zeros(size(stellarMass));
        xi=ones(size(stellarMass));
        
        cv=cvir_Mvir(satHalo,0,'random');
        lambda=lambda_prime(satHalo);
        
        % calculate Rd
        %cnt=0;
        %tot=(length(betaVals)+1)*(length(fgsVals)+1);
        for i=1:length(betaRange)+1
            
            
            %             for j=1:length(fgsVals)+1
            %                 cnt=cnt+1;
            %
            if i~=length(betaRange)+1
                
                beta=ones(size(stellarMass)).*betaRange(i);
            else
                beta=betaRange(1)+diff(betaRange([1 end])).*rand(ngal,1);
            end
            
            
            %                 if j~=length(fgsVals)+1
            %
            %                     fgs=ones(size(stellarMass)).*fgsVals(j);
            %                 else
            %                     fgs=fgsRange(1)+diff(fgsRange([1 end])).*rand(ngal,1);
            %                 end
            
            
            res=rscale_mmw_array(stellarMass,'fg',fgs,'beta',beta,'fb',fbs,'xi',xi,...
                'Mv',satHalo,'cv',cv,'lambda',lambda,'noshow');
            %      res=rscale_mmw_array(stellarMass,'fg',fgs,'beta',beta,...
            %          'Mv',satHalo,'cv',cv,'lambda',lambda,'noshow');
            %
            cata(i).Ms=stellarMass;
            cata(i).rd=res.rd;
            cata(i).Mv=res.Mv;
            cata(i).md=res.md;
            %cata,i).BT=res.BT;
            cata(i).fbs=fbs;
            cata(i).fgs=fgs;
            cata(i).xi=xi;
            cata(i).cv=res.cv;
            cata(i).lambda=res.lambda;
            %cata(j,i).sigma=stellarMass./(2*pi*res.rd.^2);
            %cata(j,i).sigmaeff=stellarMass./(2*pi*res.rd.^2).*sigma_factor('half');
            %cata(j,i).sigma5090=stellarMass./(2*pi*res.rd.^2).*sigma_factor('5090');
            cata(i).beta=beta;
            % cata(j,i).qmin=qmin;
            % cata(j,i).qsub=qsub;
            % cata(j,i).eps=eps;
            cata(i).lambdaMask=res.lambdaMask;
            
            fprintf('completed %i %% of catalogs \n',ceil(i/(length(betaRange)+1)*100));
            
        end
        
        fullCatalog(k).cata=cata;
        toc
    end
    
    
    name=sprintf('catalog_%s_ng100_fgs%s_snp%s.mat',tags{kk},simEmulate,num2str(snap));
    save([DEFAULT_MATFILE_DIR '/' name],'fullCatalog','-v7.3')
    fprintf('saving to: %s \n',[DEFAULT_MATFILE_DIR '/' name])
    
end