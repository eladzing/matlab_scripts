
function [rq qbig] = concat_qprofs(q1,q2,q4,q8)  %clustername,typ,aexp,cm,cen_fac) 
%calculate qprofs and attach them at the ends. 
%if ~exist('cm')
%    cm = [0,0,0];
%end

%[q1 q2 q4 q8]=cooling_profs(clustername,typ,aexp,cm,cen_fac); 

%global hub;
global NCELL;

qbig=zeros(1,8*NCELL);

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

qbig2=interp1(r2,q2./2,rb2,'cubic');
qbig4=interp1(r4,q4./4,rb4,'cubic');
qbig8=interp1(r8,q8./8,rb8,'cubic');

rq=rb8;
clear r2 r4 r8 rb2 rb4 rb8    

qbig(1:NCELL)=q1;
qbig(NCELL+1:2*NCELL)=qbig2(NCELL+1:2*NCELL);
qbig(2*NCELL+1:4*NCELL)=qbig4(2*NCELL+1:4*NCELL);
qbig(4*NCELL+1:8*NCELL)=qbig8(4*NCELL+1:8*NCELL);

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





