%% Dissipation Energy Calculation
% We try to find the dissipation term by calculating all 
% other energy components, and their instantaneous change 
% by u grad E - this assumes a steady state: partital 
% derivative of E is zero.
% this time all values are per unit mass.

%% Calculate energy components. 
% all components are per unit volume
% units:  M_sun*(km/sec)^2 per Mpc^3 

list=[101 102 103 104 105 106 107  3  5  6  7  9 10 11 14 24];
%mask=[ 1   0   1   0   1   1   0   1  1  1  0  1  1  1  1  1];
mask=false(size(list));

mask(10)=true;

innermask=true(256,256,256);
innermask(65:192,65:192,65:192)=false;


tag1='test_dis_energy';
% pos=zeros(4,6);
% neg=zeros(4,6);
% 
% nc=zeros(4,6);
% avg=zeros(4,4,6);%1-total, 2 pos, 3 neg, 4 abs 
% med=zeros(4,4,6);
% rms=zeros(1,4,6);
% stdv=zeros(4,4,6);

units;
global syst
for i=1:length(list);
    if~mask(i) 
        continue;
    end
    
   new_env(sprintf('CL%d',list(i)),'csf','a1',syst);
   

  % qdiss_array_np=cell(6,1);
   
global CLUSTER;
%global VCM
global NCELL;
%global zred
global hub;


%diss=zeros(NCELL,NCELL,NCELL);
qdiss=zeros(NCELL,NCELL,NCELL);%,4);




for b=[2 4];
    bmask=true(256,256,256);
    boxx=2^(b-1);
    boxx
%     if boxx==1
%        bmask=true(256,256,256);
%        cl=ceil(0.05*get_rvir()/(boxx/hub/NCELL));
%        bmask(129-cl:128+cl,129-cl:128+cl,129-cl:128+cl)=false; 
%     else
%         bmask=innermask;
%     end
    
    
    cellsize=boxx/hub/NCELL;
    vol=cellsize^3;
    
[Vxx Vyy Vzz] = get_velocities(boxx,'defualt');
divu=calc_divV_single(Vxx,Vyy,Vzz,boxx);

%rog=RHOG(boxx);
[pdm,pgs,~,ro]=poteng(boxx);
tp=T(boxx); 

if pergramflag
    ro=1;
end

mas=ro.*vol;

%ro=ones(size(ro));

Pr=ro.*tp.*f_pr;
ekin=0.5.*ro.*(Vxx.^2+Vyy.^2+Vzz.^2).*bmask;
egrav=(pdm+pgs).*ro.*bmask;
eth=ro.*tp.*f_th.*bmask;
clear pdm pgs tp 

%qcube=cooling_cube(boxx)./ro;%.*bmask;  %% per gram
qcube=cooling_cube(boxx);%.*bmask;
pdv=-1.*Pr.*divu.*fu.*bmask; % comes out positive for convergent flow 

clear Pr divu
    
%% calculate u grad 



 [gKx gKy gKz]=gradient(ekin,cellsize);
 deltaK=(Vxx.*gKx + Vyy.*gKy + Vzz.*gKz).*fu;
 clear gKx gKy gKz 
 
 [gPx gPy gPz]=gradient(egrav,cellsize);
 deltaP=(Vxx.*gPx + Vyy.*gPy + Vzz.*gPz).*fu;
 clear gPx gPy gPz 
 
 [gTx gTy gTz]=gradient(eth,cellsize);
 deltaT=(Vxx.*gTx + Vyy.*gTy + Vzz.*gTz).*fu;
 clear gTx gTy gTz Vxx Vyy Vzz

normfac=get_mvir()/get_rvir()*f_eng;
%mkmap('data',ekin,'clims',[0 6],'proj',[1 0 0],'normalize',normfac','title','Kin','boxx',boxx,'log','on','bartag','$\log E\,[\mathrm{G M_{vir}/R_{vir}}]$');
%mkmap('data',abs(egrav),'clims',[0 6],'proj',[1 0 0],'normalize',normfac','title','Pot','boxx',boxx,'log','on','bartag','$\log E\,[\mathrm{G M_{vir}/R_{vir}}]$');
%mkmap('data',eth,'clims',[0 6],'proj',[1 0 0],'normalize',normfac','title','Therm','boxx',boxx,'log','on','bartag','$\log E\,[\mathrm{G M_{vir}/R_{vir}}]$');

%mkmap('data',deltaK,'proj',[1 0 0],'normalize',normfac','title','$\Delta$Kin','boxx',boxx,'log','on','bartag','$\log E\,[\mathrm{G M_{vir}/R_{vir}/Mpc^3/yr}]$');
%mkmap('data',deltaP,'proj',[1 0 0],'normalize',normfac','title','$\Delta$Pot','boxx',boxx,'log','on','bartag','$\log E\,[\mathrm{G M_{vir}/R_{vir}/Mpc^3/yr}]$');
%mkmap('data',deltaT,'proj',[1 0 0],'normalize',normfac','title','$\Delta$Therm','boxx',boxx,'log','on','bartag','$\log E\,[\mathrm{G M_{vir}/R_{vir}/Mpc^3/yr}]$');

%mkmap('data',qcube,'proj',[1 0 0],'normalize',normfac','title','Q','boxx',boxx,'log','on','bartag','$\log E\,[\mathrm{G M_{vir}/R_{vir}/Mpc^3}/yr]$');
%mkmap('data',pdv,'proj',[1 0 0],'normalize',normfac','title','PdV','boxx',boxx,'log','on','bartag','$\log E\,[\mathrm{G M_{vir}/R_{vir}/Mpc^3}/yr]$');
qdiss=deltaK+deltaP+pdv-qcube;

map_qdiss
end
end


%clear ekin egrav eth 
% 
% 
% nc(b,:)=sum(sum(sum(bmask)));
% 
% for k=1:6
%     switch k
%         case 1
%             dcub=qdiss(:,:,:,b);
%         case 2  
%             dcub=deltaK;
%         case 3
%             dcub=deltaP;
%         case 4  
%             dcub=deltaT;
%         case 5
%             dcub=pdv;
%         case 6  
%             dcub=qcube;
%     end
%             
%     
% pos(b,k)=sum(sum(sum(dcub>0)))./nc(b);
% neg(b,k)=sum(sum(sum(dcub<0)))./nc(b);
% 
% %l=dcub(bmask==1);
% %lp=dcub(dcub>0);
% %ln=dcub(dcub<0);
% 
% clear dcub
% 
% % avg(1,b,k)=mean(l);
% % avg(2,b,k)=mean(lp);
% % avg(3,b,k)=mean(ln);
% % avg(4,b,k)=mean(abs(l));
% % 
% % med(1,b,k)=median(l);
% % med(2,b,k)=median(lp);
% % med(3,b,k)=median(ln);
% % med(4,b,k)=median(abs(l));
% % 
% % rms(b,k)=sqrt(mean(l.^2));
% % 
% % stdv(1,b,k)=sqrt(var(l));
% % stdv(2,b,k)=sqrt(var(lp));
% % stdv(3,b,k)=sqrt(var(ln));
% % stdv(4,b,k)=sqrt(var(abs(l)));
% % 
% % clear l lp ln
% 
% end
% clear deltaK deltaP deltaT pdv qcube qdiss
% 
% %qdiss(:,:,:,b)=diss;
% 
% 
% end
% 
% %qdiss_array{i,1}=CLUSTER;
% qdiss_array_np{1}=pos;
% qdiss_array_np{2}=neg;
% qdiss_array_np{3}=nc;
% %qdiss_array{4}=rms;
% %qdiss_array{5}=stdv;
% %qdiss_array{6}='name q_dis(4box) avg med rms std'; 
% 
% save(sprintf('mat_files/qdiss_np_%s_%s.mat','a1',CLUSTER),'qdiss_array_np');
% 
% clear qdiss_array_np
% 
% end
% 
