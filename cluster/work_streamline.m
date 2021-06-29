%%
close all
%global FILE_FORMAT;
%FILE_FORMAT = '../data/%s_a1.001L%dMpc.dat';

MPSec = 8;
[Vxx Vyy Vzz VcmX VcmY VcmZ] = V_Vcm(MPSec);

[Vrr cosa] = Vr(MPSec);
%%
Vabs = sqrt(Vxx.^2 + Vyy.^2 + Vzz.^2);
mask = bitand(cosa<0, Vabs>200);

Vxx(~mask) = 0;
Vyy(~mask) = 0;
Vzz(~mask) = 0;

DILUTE = 50;
[startX1 startY1 startZ1] = meshgrid([1 256], [1:DILUTE:256], [1:DILUTE:256]);
[startX2 startY2 startZ2] = meshgrid([1:DILUTE:256], [1 256], [1:DILUTE:256]);
[startX3 startY3 startZ3] = meshgrid([1:DILUTE:256], [1:DILUTE:256], [1 256]);

%streamline(Vxx, Vyy, Vzz, [startX1(:) startX2(:) startX3(:)], [startY1(:) startY2(:) startY3(:)], [startZ1(:) startZ2(:) startZ3(:)])

%%
[meshX, meshY, meshZ] = meshgrid(1:size(Vxx,1), size(Vyy,2):-1:1, 1:size(Vzz,3));
verts = stream3(meshX, meshY, meshZ, Vxx, Vyy, Vzz, [startX1(:) startX2(:) startX3(:)], [startY1(:) startY2(:) startY3(:)], [startZ1(:) startZ2(:) startZ3(:)]);
sl = streamline(verts);

iverts = interpstreamspeed(meshX,meshY,meshZ, Vxx,Vyy,Vzz, verts,.025);
axis tight; view(30,30); daspect([1 1 .125])
camproj perspective; camva(8)
set(gca,'DrawMode','fast')
box on
streamparticles(iverts,35,'animate',10,'ParticleAlignment','on')