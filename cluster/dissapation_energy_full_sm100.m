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


tag1={'testeng_dis' 'testeng_dis_prgr'};
% pos=zeros(4,6);
% neg=zeros(4,6);
%
% nc=zeros(4,6);
% avg=zeros(4,4,6);%1-total, 2 pos, 3 neg, 4 abs
% med=zeros(4,4,6);
% rms=zeros(1,4,6);
% stdv=zeros(4,4,6);

units;
% smoothing parameters: 
smoothtype='box';
%smoothscale=0.1; % 150 kpc

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
    
    normfac(1)=get_mvir()/get_rvir()*f_eng;
    normfac(2)=normfac(1)/get_mvir()/Ms;
    
    %diss=zeros(NCELL,NCELL,NCELL);
    qdiss=zeros(NCELL,NCELL,NCELL,2);%,4);
    ekin=zeros(NCELL,NCELL,NCELL,2);
    egrav=zeros(NCELL,NCELL,NCELL,2);
    eth=zeros(NCELL,NCELL,NCELL,2);
    qcube=zeros(NCELL,NCELL,NCELL,2);
    pdv=zeros(NCELL,NCELL,NCELL,2);
    deltaK=zeros(NCELL,NCELL,NCELL,2);
    deltaP=zeros(NCELL,NCELL,NCELL,2);
    deltaT=zeros(NCELL,NCELL,NCELL,2);
    
    for b=[2 4];
        bmask=true(256,256,256);
        boxx=2^(b-1)
        
        %     if boxx==1
        %        bmask=true(256,256,256);
        %        cl=ceil(0.05*get_rvir()/(boxx/hub/NCELL));
        %        bmask(129-cl:128+cl,129-cl:128+cl,129-cl:128+cl)=false;
        %     else
        %         bmask=innermask;
        %     end
        
        
        cellsize=boxx/hub/NCELL;
        vol=cellsize^3;
        smscale=ceil(smoothscale/cellsize);
    smscale=smscale+1-mod(smscale,2);
    
        [Vxx Vyy Vzz] = get_velocities(boxx,'defualt');
        Vxx=smooth3(Vxx,smoothtype,smscale);
    Vyy=smooth3(Vyy,smoothtype,smscale);
    Vzz=smooth3(Vzz,smoothtype,smscale);
        
        
        
        divu=calc_divV_single(Vxx,Vyy,Vzz,boxx);
        
       
        [pdm,pgs,~,ro]=poteng(boxx);
    pot=smooth3(pdm+pgs,smoothtype,smscale);
    %pgs=smooth3(pgs,smoothtype,smscale);
    ro=smooth3(ro,smoothtype,smscale);
    clear ro pdm pgs
         tp=smooth3(T(boxx),smoothtype,smscale);
      
        
        mass=ro./vol;
        
        Pr=ro.*tp.*f_pr;
        ekin(:,:,:,1)=0.5.*ro.*(Vxx.^2+Vyy.^2+Vzz.^2).*bmask;
        egrav(:,:,:,1)=(pdm+pgs).*ro.*bmask;
        eth(:,:,:,1)=ro.*tp.*f_th.*bmask;
        clear pdm pgs tp
        
        qcube=cooling_cube(boxx);%.*bmask;
    %qcube=cooling_cube(boxx)./ro;%.*bmask;   %% per gram
    qcube=smooth3(qcube,smoothtype,smscale);
        
        
        pdv(:,:,:,1)=-1.*Pr.*divu.*fu.*bmask; % comes out positive for convergent flow
        clear Pr divu
        
        ekin(:,:,:,2)=ekin(:,:,:,1)./mass;
        egrav(:,:,:,2)=egrav(:,:,:,1)./mass;
        eth(:,:,:,2)=eth(:,:,:,1)./mass;
        qcube(:,:,:,2)=qcube(:,:,:,1)./mass;
        pdv(:,:,:,2)= pdv(:,:,:,1)./mass;
        
        
        %% calculate u grad
        for kk=1:2
            [gKx gKy gKz]=gradient(ekin(:,:,:,kk),cellsize);
            deltaK(:,:,:,kk)=(Vxx.*gKx + Vyy.*gKy + Vzz.*gKz).*fu;
            clear gKx gKy gKz
            
            [gPx gPy gPz]=gradient(egrav(:,:,:,kk),cellsize);
            deltaP(:,:,:,kk)=(Vxx.*gPx + Vyy.*gPy + Vzz.*gPz).*fu;
            clear gPx gPy gPz
            
            [gTx gTy gTz]=gradient(eth(:,:,:,kk),cellsize);
            deltaT(:,:,:,kk)=(Vxx.*gTx + Vyy.*gTy + Vzz.*gTz).*fu;
            clear gTx gTy gTz Vxx Vyy Vzz
            
            qdiss(:,:,:,kk)=deltaK(:,:,:,kk)+deltaP(:,:,:,kk)+pdv(:,:,:,kk)-qcube(:,:,:,kk);
            %map_qdiss
            
            save('mat_files/diss_test_full.m',...
                normfac,qdiss,ekin,egrav,eth,qcube,pdv,...
                deltaK,deltaP,deltaT)      
            
            
        end
        
        
        %mkmap('data',ekin,'clims',[0 6],'proj',[1 0 0],'normalize',normfac','title','Kin','boxx',boxx,'log','on','bartag','$\log E\,[\mathrm{G M_{vir}/R_{vir}}]$');
        %mkmap('data',abs(egrav),'clims',[0 6],'proj',[1 0 0],'normalize',normfac','title','Pot','boxx',boxx,'log','on','bartag','$\log E\,[\mathrm{G M_{vir}/R_{vir}}]$');
        %mkmap('data',eth,'clims',[0 6],'proj',[1 0 0],'normalize',normfac','title','Therm','boxx',boxx,'log','on','bartag','$\log E\,[\mathrm{G M_{vir}/R_{vir}}]$');
        
        %mkmap('data',deltaK,'proj',[1 0 0],'normalize',normfac','title','$\Delta$Kin','boxx',boxx,'log','on','bartag','$\log E\,[\mathrm{G M_{vir}/R_{vir}/Mpc^3/yr}]$');
        %mkmap('data',deltaP,'proj',[1 0 0],'normalize',normfac','title','$\Delta$Pot','boxx',boxx,'log','on','bartag','$\log E\,[\mathrm{G M_{vir}/R_{vir}/Mpc^3/yr}]$');
        %mkmap('data',deltaT,'proj',[1 0 0],'normalize',normfac','title','$\Delta$Therm','boxx',boxx,'log','on','bartag','$\log E\,[\mathrm{G M_{vir}/R_{vir}/Mpc^3/yr}]$');
        
        %mkmap('data',qcube,'proj',[1 0 0],'normalize',normfac','title','Q','boxx',boxx,'log','on','bartag','$\log E\,[\mathrm{G M_{vir}/R_{vir}/Mpc^3}/yr]$');
        %mkmap('data',pdv,'proj',[1 0 0],'normalize',normfac','title','PdV','boxx',boxx,'log','on','bartag','$\log E\,[\mathrm{G M_{vir}/R_{vir}/Mpc^3}/yr]$');
        %qdiss=deltaK+deltaP+pdv-qcube;
        
        %map_qdiss
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
