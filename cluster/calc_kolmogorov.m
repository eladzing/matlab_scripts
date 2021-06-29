function [bins_x bins] = calc_kolmogorov(small)

%%
MPc = 1;

%%
%[vx, vy, vz, v1,v2,v3] = V_Vcm(MPc);

%vel_abs = (vx.^2+vy.^2+vz.^2);

%%
%small = vel_abs(97:160,97:160,97:160);
%small = vel_abs;
%small = vel_abs(range1,range2,range3);

%%
%small = zeros(size(small));
%small((1:8)*8,(1:8)*8,(1:8)*8) = 1;


small_len = length(small)
h = 0.7;
L = small_len/256*MPc/h;
%%
%%small = vel_abs;
sk = fftn(small);
sk2 = (abs(sk));
sk2s = fftshift(sk2);
%imshow(log(squeeze(sum(sk2s(32:33,:,:),1))),[])
%imshow(log(squeeze(sum(sk2s(:,32:33,:),2))),[])
%imshow(log(squeeze(sum(sk2s(:,:,32:33),3))),[])

%%
%sk2sr = cart2sphere_new(sk2s);
% loglog(squeeze(sum(sum(sk2sr,2),3)))
% yy = (squeeze(sum(sum(sk2sr,2),3)));
% xx = 1:length(small);
% xx = xx(:)
% yy = yy(:);
% yy2 = yy.*(xx.^2)
% ttt = (yy2(2:end)-yy2(1:end-1)).*((xx(2:end)+xx(1:end-1))/2)./((yy2(2:end)+yy2(1:end-1))/2);
% semilogx(ttt)


%%
grid_vec = (1:small_len)-((small_len+1)/2);
[xxx,yyy,zzz] = meshgrid(grid_vec,grid_vec,grid_vec);
radii = sqrt(xxx.^2+yyy.^2+zzz.^2);
%imshow((squeeze(sum(radii(32:33,:,:),1))),[])

bins = zeros(1,small_len*2);
bin_counts = bins;

for idx = 1:length(radii(:))
    rad = round(radii(idx));
    %bins(rad)= bins(rad) + sk2s(idx);
    bins(rad)= bins(rad) + sk2s(idx).^2;
    bin_counts(rad) = bin_counts(rad) + 1;
end

bins = bins./bin_counts;
bins_x = ((1:length(bins))-1)/L;