function [ output_args ] = plot_phaseDiagram_polygons(hf,polys)
%UNTITLED plot the polygons Masks on the phase diagram 
%   Detailed explanation goes here

figure(hf)
patch(polys.coldDense_polygonMask(:,1),polys.coldDense_polygonMask(:,2),'y','FaceAlpha',0.15)
hold on 
patch(polys.coldDilute_polygonMask(:,1),polys.coldDilute_polygonMask(:,2),'m','FaceAlpha',0.15)
patch(polys.sf_polygonMask(:,1),polys.sf_polygonMask(:,2),'g','FaceAlpha',0.15)
patch(polys.warmHot_polygonMask(:,1),polys.warmHot_polygonMask(:,2),'b','FaceAlpha',0.15)
patch(polys.hotTemp_polygonMask(:,1),polys.hotTemp_polygonMask(:,2),'b','FaceAlpha',0.15)

end