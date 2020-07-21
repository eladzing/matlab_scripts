function [v_x v_y]=mk_vfield(Vxx,Vyy,Vzz,weight,slind) 
  
  Vxx=Vxx.*weight;Vyy=Vyy.*weight;Vzz=Vzz.*weight;
  v_x=zeros([size(Vxx,1),size(Vxx,2),3]);
  v_y=zeros([size(Vxx,1),size(Vxx,2),3]);
  
  for projection = 1:3
        
    switch projection
      case 1
        v_x(:,:,projection) = squeeze(mean_nan(Vzz(slind,:,:),projection)./ mean_nan(weight(slind,:,:),projection));
        v_y(:,:,projection)  = squeeze(mean_nan(Vyy(slind,:,:),projection)./ mean_nan(weight(slind,:,:),projection));
            
      case 2
        v_y(:,:,projection)  = transpose(squeeze(mean_nan(Vzz(:,slind,:),projection)./ mean_nan(weight(:,slind,:),projection)));
        v_x(:,:,projection)  = transpose(squeeze(mean_nan(Vxx(:,slind,:),projection)./ mean_nan(weight(:,slind,:),projection)));
        
      case 3
        v_y(:,:,projection)  = transpose(squeeze(mean_nan(Vyy(:,:,slind),projection)./ mean_nan(weight(:,:,slind),projection)));
        v_x(:,:,projection)  = transpose(squeeze(mean_nan(Vxx(:,:,slind),projection)./ mean_nan(weight(:,:,slind),projection)));
        
    end
      
  end
    
    
end
  

