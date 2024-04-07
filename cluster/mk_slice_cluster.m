

function slice=mk_slice_cluster(cube,weight,slind,type)
% creates 3 slices for x,y,z projections. slices can be weighted.
%  slind is the index range for the slice

if ~exist('type','var')
    type='avg';
end


cube=cube.*weight;

slice=zeros([size(cube,1),size(cube,2),3]);

switch lower(type)
    case 'avg'
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
    case 'sum'
        for projection = 1:3
            switch projection
                case 1
                    slice(:,:,projection) = squeeze(sum(cube(slind,:,:),projection));
                case 2
                    slice(:,:,projection) = transpose(squeeze(sum(cube(:,slind,:),projection)));
                case 3
                    slice(:,:,projection) = transpose(squeeze(sum(cube(:,:,slind),projection)));
            end
        end
    otherwise
        error('%s - Illegal type: %s',current_function().upper,type);
end

end
