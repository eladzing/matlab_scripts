for i=[1 2 4 8];
    
     %[pgs,pdm,pdot]=Pot_temp_all(i);
     figure

     
     pp=-1.*log10(abs(pdm));
     mn=min(min(min(pp)));
     mx=max(max(max(pp)));
         
     
     for j=1:256
     
     %subplot(2,2,1)
     
     imagesc(squeeze(pp(:,:,j)))
     colorbar
     caxis([mn mx]) 
     %subplot(2,2,2)
     %imagesc(squeeze(pgs(:,:,j)))
     
     %subplot(2,2,3)
     %imagesc(squeeze(pdot(:,:,j)))
     
     pause
     end
     close gcf
     pause
end