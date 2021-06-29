
% prepare edge definitions 
plot_shockedge

list=[101 102 103 104 105 106 107  3  5  6  7  9 10 11 14 24];
%mask=[ 0   0   1   0   0   1   0   0  0  1  0  0  0  1  1  0];
%mask=[ 0   0   0   0   0   0   0   0  0  1  0  0  0  1  1  0];
mask=ones(size(list));
%list=[103 106 6 10 11 14];
%mask=[ 1   1  1  1  1  1]; 
%boxx=8;

mas_a1=zeros(length(list),4);
mas_a06=zeros(length(list),4);



for i=1:length(list);
%    continue % THIS IS HERE SINCE THESE PICS ARE DONE
    if~mask(i) 
        continue;
    end
    
  % z=0 
      
    new_env(sprintf('CL%d',list(i)),'csf','a1');
    global hub
    global NCELL 

   
   for boxx=8;  
       vol=boxx./hub/NCELL;
    
       rog=RHOG(boxx); 
       rodm=RHODM(boxx);
       SVIR=get_tvir.*get_rvir.^2./(3.*get_mvir./(4.*pi)).^(2./3); %normalization for entropy
       ent=S(boxx)./SVIR;
       rc=mk_rcube(boxx,rog);
      
       edge=max(edge_a1(i,:).*rv_a1(i));
       
       bmask= (rc<edge & ent>=1);
       
       clear rc ent
       mg=sum(sum(sum(rog.*bmask.*vol)));
       md=sum(sum(sum(rodm.*bmask.*vol)));
       mas_a1(i,:)=[mg md mg+md get_mvir];
       
       
    end

% z=0.6
       
    new_env(sprintf('CL%d',list(i)),'csf','a06');
    %global hub
    %global NCELL 

   
   for boxx=8;  
       vol=boxx./hub/NCELL;
    
       rog=RHOG(boxx); 
       rodm=RHODM(boxx);
       SVIR=get_tvir.*get_rvir.^2./(3.*get_mvir./(4.*pi)).^(2./3); %normalization for entropy
       ent=S(boxx)./SVIR;
       rc=mk_rcube(boxx,rog);
      
       edge=max(edge_a06(i,:).*rv_a06(i));
       
       bmask= (rc<edge & ent>=1);
       
       clear rc ent
       mg=sum(sum(sum(rog.*bmask.*vol)));
       md=sum(sum(sum(rodm.*bmask.*vol)));
       mas_a06(i,:)=[mg md mg+md get_mvir];
       
       
   end
end
