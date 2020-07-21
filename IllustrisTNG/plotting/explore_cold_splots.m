wcu=cubeStr.cube.*cubeStr.weights;


crc=-1.*crStr.cube./1e-22;

ro=roStr.cube;
eng=engStr.cube;

figure('position',[1101 -234  1446 1176]);

for i=1:256
   % titlemine(sprintf('%i',i));
    
    subplot(2,2,1)
    imagesc(transpose(log10(squeeze(cr2(:,i,:)))))
    caxis([6 13]);% caxis([-10 -4])
    titlemine(sprintf('tc %i',i));
    
    subplot(2,2,2)
    imagesc(transpose(log10(squeeze(crc(:,i,:)))))
    caxis([-6 0])
    titlemine(sprintf('cooling rate %i',i));
    
    subplot(2,2,3)
    imagesc(transpose(log10(squeeze(ro(:,i,:)))))
    caxis([-8 -3])
    titlemine(sprintf('$\rho$ %i',i));
    
    subplot(2,2,4)
    imagesc(transpose(log10(squeeze(eng(:,i,:)))))
    caxis([2 5])
    titlemine(sprintf('Energy %i',i));
    
    
    pause
end

