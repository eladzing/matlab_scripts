%% claculate the average density in the same slab used to find CF's.
%%
%%ll=65:192;
rhoVal=1e-3; %sqrt(0.7)*

for i=1:length(cfList)
    new_env(cfList(i).cluster)
    
    rhoRad=find_radius_of_density(rhoVal,'n').*0.7;
    
    r1=sqrt(sum(cfList(i).cf1.^2,2));
    
    if ~isempty(r1)
        for j=1:length(r1)
            
            ll=125:132;
            if r1(j) <=0.5
                boxx=1;
                
            elseif r1(j)<=1
                boxx=2;
            elseif r1(j)<=2
                boxx=4;
            elseif r1(j)<=4
                boxx=8;
                ll=121:136;
            else
                error('r1 all fucked up at %i, %g',j,r2(j))
            end
            
            nn=RHOGN(boxx);
            
            switch cfList(i).cf1Prj(j)
                case 2
                    nn=squeeze(sum(nn(ll,:,:),1))./length(ll);
                case 3
                    nn=squeeze(sum(nn(:,ll,:),2))./length(ll);
                    nn=nn';
                case 1
                    nn=squeeze(sum(nn(:,:,ll),3))./length(ll);
                otherwise
                    error('cf1Prj fucked up at %i',j)
            end
            
            ind1=ceil((cfList(i).cf1(j,1)+boxx/2)./boxx.*256);
            ind2=ceil((cfList(i).cf1(j,2)+boxx/2)./boxx.*256);
            
            val= max(max(nn(ind1-2:ind1+2,ind2-2:ind2+2)));
            
            cfList(i).cf1Count(j,1)=val>rhoVal;
            cfList(i).cf1Count(j,2)=r1(j)<=rhoRad;
        end
    end
    
    r2=sqrt(sum(cfList(i).cf2.^2,2));
    if ~isempty(r2)
        for j=1:length(r2)
            ll=125:132;
            if r2(j) <=0.5
                boxx=1;
            elseif r2(j)<=1
                boxx=2;
            elseif r2(j)<=2
                boxx=4;
            elseif r2(j)<=4
                boxx=8;
                ll=121:136;
            else
                error('r2 all fucked up at %i, %g',j,r2(j))
            end
            
            nn=RHOGN(boxx);
            
            switch cfList(i).cf2Prj(j)
                case 1
                    nn=squeeze(sum(nn(ll,:,:),1))./length(ll);
                case 2
                    nn=squeeze(sum(nn(:,ll,:),2))./length(ll);
                    nn=nn';
                case 3
                    nn=squeeze(sum(nn(:,:,ll),3))./length(ll);
                otherwise
                    error('cf1Prj fucked up at %i',j)
            end
            
            ind1=ceil((cfList(i).cf2(j,1)+boxx/2)./boxx.*256);
            ind2=ceil((cfList(i).cf2(j,2)+boxx/2)./boxx.*256);
            
            
            val=max(max(nn(ind1-2:ind1+2,ind2-2:ind2+2)));
            
            cfList(i).cf2Count(j,1)=val>rhoVal;
            cfList(i).cf2Count(j,2)=r2(j)<=rhoRad;
        end
    end
end


%% count results 
for i=1:length(cfList)
    if ~isempty(cfList(i).cf1Count)
        obs1(i)=sum(cfList(i).cf1Count(:,1));
        obs1P(i)=sum(cfList(i).cf1Count(:,2));
        
        tot1(i)=length(cfList(i).cf1Count(:,1));
    end
    
    if ~isempty(cfList(i).cf2Count)
        obs2(i)=sum(cfList(i).cf2Count(:,1));
        obs2P(i)=sum(cfList(i).cf2Count(:,2));
        tot2(i)=length(cfList(i).cf2Count(:,1));
    end
        
end
    