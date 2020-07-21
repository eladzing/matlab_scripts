function plot_bird2slice(MPSec, thick, rho_range, temp_range, mask)
%thick = half of the thickness

dilute = 4; %dilute parameter for vector field

%load data
[Vxx Vyy Vzz VcmX VcmY VcmZ] = V_Vcm(MPSec);
rhog = RHOG(MPSec);
temp = T(MPSec);

if ~exist('mask')
    mask = bitand(bitand(rhog>rho_range(1),rhog<rho_range(2)), bitand(temp>temp_range(1),temp<temp_range(2)));
end

Vzz(~mask) = NaN;
Vxx(~mask) = NaN;
Vyy(~mask) = NaN;

%figure;

for projection = 1:3
    
    rgbcube = zeros([256,256,3]);
    switch projection
        case 1
            rho_channel = log10(squeeze(sum(rhog(129-thick:128+thick,:,:),projection)/(thick*2)));
            temp_channel = log10(squeeze(sum(temp(129-thick:128+thick,:,:),projection)/(thick*2)));
            subplot(2,2,4);
        case 2
            rho_channel = transpose(log10(squeeze(sum(rhog(:,129-thick:128+thick,:),projection)/(thick*2))));
            temp_channel = transpose(log10(squeeze(sum(temp(:,129-thick:128+thick,:),projection)/(thick*2))));
            subplot(2,2,1);
        case 3
            rho_channel = transpose(log10(squeeze(sum(rhog(:,:,129-thick:128+thick),projection)/(thick*2))));
            temp_channel = transpose(log10(squeeze(sum(temp(:,:,129-thick:128+thick),projection)/(thick*2))));
            subplot(2,2,3);
    end

    hold on;

    rho_channel = rho_channel - min(rho_channel(:));
    rho_channel = rho_channel ./ max(rho_channel(:));
    temp_channel = temp_channel - min(temp_channel(:));
    temp_channel = temp_channel ./ max(temp_channel(:));

    rgbcube(:,:,1) = rho_channel;
    rgbcube(:,:,3) = temp_channel;

    diluted_len = length(1:dilute:256);
    diluted_jump = round(256/(diluted_len-1));
    [xx yy] = meshgrid(1:diluted_jump:256,1:diluted_jump:256);

    %figure;
    plot(1:256)
    hold on
    imagesc(-MPSec/2:MPSec/2,-MPSec/2:MPSec/2,rgbcube);

    switch projection
        case 1
            v_x = squeeze(mean_nan(Vzz(129-thick:128+thick,1:dilute:256,1:dilute:256),projection));
            v_y = squeeze(mean_nan(Vyy(129-thick:128+thick,1:dilute:256,1:dilute:256),projection));
            counts = sum(mask(129-thick:128+thick,:,:),projection) > 0;
            [yyy xxx] = ind2sub([256 256], find(counts));
            xlabel('Z');ylabel('Y');title(sprintf('ZY slice (%dMpc, thickness=%d)',MPSec,thick*2));
        case 2
            v_y = transpose(squeeze(mean_nan(Vzz(1:dilute:256,129-thick:128+thick,1:dilute:256),projection)));
            v_x = transpose(squeeze(mean_nan(Vxx(1:dilute:256,129-thick:128+thick,1:dilute:256),projection)));
            counts = sum(mask(:,129-thick:128+thick,:),projection) > 0;
            [xxx yyy] = ind2sub([256 256], find(counts));
            xlabel('X');ylabel('Z');title(sprintf('XZ slice (%dMpc, thickness=%d)',MPSec,thick*2));
        case 3
            v_y = transpose(squeeze(mean_nan(Vyy(1:dilute:256,1:dilute:256,129-thick:128+thick),projection)));
            v_x = transpose(squeeze(mean_nan(Vxx(1:dilute:256,1:dilute:256,129-thick:128+thick),projection)));
            counts = sum(mask(:,:,129-thick:128+thick),projection) > 0;
            [xxx yyy] = ind2sub([256 256], find(counts));
            xlabel('X');ylabel('Y');title(sprintf('XY slice (%dMpc, thickness=%d)',MPSec,thick*2));
    end

    
    plot(xxx,yyy,'.y','MarkerSize',1);

    good_idxs = find(~isnan(v_x));
    quiver(xx(good_idxs),yy(good_idxs),v_x(good_idxs),v_y(good_idxs), 'g');
    %%
    xlim([0 256]);ylim([0 256]);
end
