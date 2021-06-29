function [pcube_profile masked_profile notmasked_profile] = calc_masked_profile(MPSec, pcube, mask, wcube)
%pcube - cube in polar coordinates
%mask - mask to the cube
%wcube - wheights cube

if (~exist('wcube'))
    pcube_profile = squeeze(sum(sum(pcube,3),2))'/(256^2);

    pcube1 = pcube;
    pcube1(~mask) = 0;
    masked_profile = squeeze(sum(sum(pcube1,3),2))'./(sum(sum(mask,3),2))';

    pcube1 = pcube;
    pcube1(mask) = 0;
    notmasked_profile = squeeze(sum(sum(pcube1,3),2))'./(sum(sum(~mask,3),2))';
else
    pcube_profile = squeeze(sum(sum(pcube.*wcube,3),2))'./squeeze(sum(sum(wcube,3),2))';

    pcube1 = pcube.*wcube;
    pcube1(~mask) = 0;
    wcube1 = wcube;
    wcube1(~mask) = 0;
    masked_profile = squeeze(sum(sum(pcube1,3),2))'./(sum(sum(wcube1,3),2))';

    pcube1 = pcube.*wcube;
    pcube1(mask) = 0;
    wcube1 = wcube;
    wcube1(mask) = 0;
    notmasked_profile = squeeze(sum(sum(pcube1,3),2))'./(sum(sum(wcube1,3),2))';
end