units

global DEFAULT_MATFILE_DIR
ng=100;
if readFlag
    load([DEFAULT_MATFILE_DIR '/catalog_masses.mat'])
end

nCats=1;
nCats=min(nCats,length(massCatalog));
for k=1:nCats
    fprintf('beginning iteration %i of %i \n',k,nCats);
    tic
    [~,ind]=sort(rand(size(massCatalog(k).stellarMass)));
    ind=ind(1:ng);
    haloMass=double(massCatalog(k).haloMass(ind));
    stellarMass=double(massCatalog(k).stellarMass(ind));
    
    ngal=length(stellarMass);
    
    %fgs=ones(size(stellarMass));
    beta=ones(size(stellarMass));
    
    %betaRange=[1/2.5  1/2 1/1.5  1  1.5 2  2.5];
    %fgsRange=[0.05 0.2 0.35 0.5 0.65 0.8 0.95];
    
    %betaRange=[1/5 1/3 1/2 1/1.5 4/5 1.0]; 
    %fgsRange=[0.05 0.1 0.2 0.3 0.4 0.5];
    
    
    fbs=zeros(size(stellarMass));
    xi=ones(size(stellarMass));
    
    cv=cvir_Mvir(haloMass,0,'random');
    lambda=lambda_prime(haloMass);
    
    % calculate Rd
    cnt=0;
    tot=(length(betaRange)+1)*(length(fgsRange)+1);
    for i=1%:length(betaRange)+1
        
        
        for j=1%:length(fgsRange)+1
            cnt=cnt+1;
            
            if i~=length(betaRange)+1
                
                beta=ones(size(stellarMass)).*betaRange(i);
            else
                beta=betaRange(1)+diff(betaRange([1 end])).*rand(ngal,1);
            end
            
            
            if j~=length(fgsRange)+1
                
                fgs=ones(size(stellarMass)).*fgsRange(j);
            else
                fgs=fgsRange(1)+diff(fgsRange([1 end])).*rand(ngal,1);
            end
            
            
            res=rscale_mmw_array(stellarMass,'fg',fgs,'beta',beta,'fb',fbs,'xi',xi,...
                'Mv',haloMass,'cv',cv,'lambda',lambda,'noshow');
            %      res=rscale_mmw_array(stellarMass,'fg',fgs,'beta',beta,...
            %          'Mv',haloMass,'cv',cv,'lambda',lambda,'noshow');
            %
            cata(j,i).Ms=stellarMass;
            cata(j,i).rd=res.rd;
            cata(j,i).Mv=res.Mv;
            cata(j,i).md=res.md;
            cata(j,i).BT=res.BT;
            cata(j,i).fbs=fbs;
            cata(j,i).fgs=fgs;
            cata(j,i).xi=xi;
            cata(j,i).cv=res.cv;
            cata(j,i).lambda=res.lambda;
            cata(j,i).sigma=stellarMass./(2*pi*res.rd.^2);
            cata(j,i).sigmaeff=stellarMass./(2*pi*res.rd.^2).*sigma_factor('half');
            cata(j,i).sigma5090=stellarMass./(2*pi*res.rd.^2).*sigma_factor('5090');
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
