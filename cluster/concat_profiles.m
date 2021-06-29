%
function [rpr prbig] = concat_profiles(pr1,pr2,pr4,pr8)  %clustername,typ,aexp,cm,cen_fac) 
% mix together profiles from different boxes into one profile with bins the
% size of the smallest grid box (1/256 Mpc)
% profiles are assumed to be of units ##/Mpc i.e. summed over a shell and
% divided by the width of the shell

%if ~exist('cm')
%    cm = [0,0,0];
%end

%[q1 q2 q4 q8]=cooling_profs(clustername,typ,aexp,cm,cen_fac); 

%global hub;
global NCELL;

prbig=zeros(1,8*NCELL);

ri=1:NCELL;

%% method 1: interpolate values 
%r1=(ri-0.5).*(1./2./NCELL);
r2=(ri-0.5).*(2./2./NCELL);
r4=(ri-0.5).*(4./2./NCELL);
r8=(ri-0.5).*(8./2./NCELL);

rb2=1:2*NCELL;rb4=1:4*NCELL;rb8=1:8*NCELL;

rb2=(rb2-0.5).*(1./2./NCELL);
rb4=(rb4-0.5).*(1./2./NCELL);
rb8=(rb8-0.5).*(1./2./NCELL);

prbig2=interp1(r2,pr2./2,rb2,'cubic');
prbig4=interp1(r4,pr4./4,rb4,'cubic');
prbig8=interp1(r8,pr8./8,rb8,'cubic');

rpr=rb8;
clear r2 r4 r8 rb2 rb4 rb8    

prbig(1:NCELL)=pr1;
prbig(NCELL+1:2*NCELL)=prbig2(NCELL+1:2*NCELL);
prbig(2*NCELL+1:4*NCELL)=prbig4(2*NCELL+1:4*NCELL);
prbig(4*NCELL+1:8*NCELL)=prbig8(4*NCELL+1:8*NCELL);

%% method 2: simply divide each shell into identical value cells 

% qbig8=zeros(1,8*NCELL);
% qbig4=zeros(1,4*NCELL);
% qbig2=zeros(1,2*NCELL);
% indx2=0; indx4=0; indx8=0;
% for i=ri %NCELL/2+1:NCELL
%         
%     qv=q2(i)./2;
%     for k=1:2
%         indx2=indx2+1;
%         qbig2(indx2)=qv;
%     end
%     
%     qv=q4(i)./4;
%     for k=1:4
%         indx4=indx4+1;
%         qbig4(indx4)=qv;
%     end
%     
%     qv=q8(i)./8;
%     for k=1:8
%         indx8=indx8+1;
%         qbig8(indx8)=qv;
%     end
%     
% end





