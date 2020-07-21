

function slice=mk_slice(cube,weight,slind)
 % creates 3 slices for x,y,z projections. slices can be weighted. 
 %  slind is the index range for the slice
 
  cube=cube.*weight;
  
  slice=zeros([size(cube,1),size(cube,2),3]);
  for projection = 1:3
    switch projection
      case 1
        slice(:,:,projection) = squeeze(sum(cube(slind,:,:),projection)./sum(weight(slind,:,:),projection));
      case 2
        slice(:,:,projection) = transpose(squeeze(sum(cube(:,slind,:),projection)./sum(weight(:,slind,:),projection)));
      case 3
        slice(:,:,projection) = transpose(squeeze(sum(cube(:,:,slind),projection)./sum(weight(:,:,slind),projection)));
    end
  end
  
end
