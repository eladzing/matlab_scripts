%% Creating median profiles for emmisivity calculations

% load data and find median profiles
new_env(7);
units;

r500=get_rvir(500);
r200=get_rvir(200);
global hub

if createFlag
    for b=1:4
        
        boxx=2^(b-1);
        fprintf('*** reading boxx %s ***\n',num2str(boxx));
        
        rho=RHOG(boxx).*(Units.Ms/Units.Mpc^3/(Units.muMass.*Units.mp));
        tmp=T(boxx);
        zmet=ZIa(boxx)+ZII(boxx);
        
        rhoSphere=cart2sphere(rho);
        tmpSphere=cart2sphere(tmp);
        zmetSphere=cart2sphere(zmet);
        
        len=length(rhoSphere);
        rhoMean=zeros(1,len);
        rhoMed=rhoMean;
        rhoStd=rhoMean;
        
        tmpMean=zeros(1,len);
        tmpMed=tmpMean;
        tmpStd=tmpMean;
        
        zmetMean=zeros(1,len);
        zmetMed=zmetMean;
        zmetStd=zmetMean;
%         
%         gamMean=zeros(1,len);
%         gamMed=gamMean;
%         gamStd=gamMean;
        
        
        for i=1:len
            arr=rhoSphere(i,:,:);
            arr=reshape(arr,[numel(arr) 1]);
            
            
            
            rhoMean(i)=mean(arr);
            rhoMed(i)=median(arr);
            rhoStd(i)=std(arr);
            
            arr=tmpSphere(i,:,:);
            arr=reshape(arr,[numel(arr) 1]);
            
            tmpMean(i)=mean(arr);
            tmpMed(i)=median(arr);
            tmpStd(i)=std(arr);
            
                      
            arr=zmetSphere(i,:,:);
            arr=reshape(arr,[numel(arr) 1]);
            
            zmetMeanW(i)=sum(sum(zmetSphere(i,:,:).*rhoSphere(i,:,:)))./sum(sum(rhoSphere(i,:,:)));            
            zmetMean(i)=mean(arr);
            zmetMed(i)=median(arr);
            zmetStd(i)=std(arr);
            
%             arr=log10(tmpSphere(i,:,:))./log10(rhoSphere(i,:,:))+1;
%             arr=reshape(arr,[numel(arr) 1]);
%             
% %             gamMean(i)=mean(arr);
%             gamMed(i)=median(arr);
%             gamStd(i)=std(arr);
        end
        
        bin=boxx/hub/2/len.*1000; %in kpc
        rOut=(1:len).*bin;
        rIn=rOut-bin;
        rb=(0.5.*(rIn+rOut));
        shellVol=4*pi/3*(rOut.^3-rIn.^3);
        
        emsStruct(b).box=boxx;
        emsStruct(b).rhoMean=rhoMean;
        emsStruct(b).rhoMed=rhoMed;
        emsStruct(b).rhoStd=rhoStd;
        
        emsStruct(b).tmpMean=tmpMean;
        emsStruct(b).tmpMed=tmpMed;
        emsStruct(b).tmpStd=tmpStd;
        
        emsStruct(b).zmetMeanW=zmetMeanW;
        emsStruct(b).zmetMean=zmetMean;
        emsStruct(b).zmetMed=zmetMed;
        emsStruct(b).zmetStd=zmetStd;
        
%         emsStruct(b).gamMean=gamMean;
%         emsStruct(b).gamMed=gamMed;
%         emsStruct(b).gamStd=gamStd;
        
        emsStruct(b).bin=bin;
        emsStruct(b).rb=rb;
        emsStruct(b).shellVol=shellVol;
    end
    
end

%% tie stuff together
rb0=emsStruct(1).rb(1:end-1);

rhoMean0=emsStruct(1).rhoMean(1:end-1);
rhoMed0=emsStruct(1).rhoMed(1:end-1);
rhoStd0=emsStruct(1).rhoStd(1:end-1);

tmpMean0=emsStruct(1).tmpMean(1:end-1);
tmpMed0=emsStruct(1).tmpMed(1:end-1);
tmpStd0=emsStruct(1).tmpStd(1:end-1);

zmetMean0=emsStruct(1).zmetMean(1:end-1);
zmetMeanW0=emsStruct(1).zmetMeanW(1:end-1);
zmetMed0=emsStruct(1).zmetMed(1:end-1);
zmetStd0=emsStruct(1).zmetStd(1:end-1);

% gamMean0=emsStruct(1).gamMean(1:end-1);
% gamMed0=emsStruct(1).gamMed(1:end-1);
% gamStd0=emsStruct(1).gamStd(1:end-1);

shellVol0=emsStruct(1).shellVol(1:end-1);

for i=2:4
    
    rb0=cat(2,rb0,emsStruct(i).rb(128:end-1));
    
    rhoMean0=cat(2,rhoMean0,emsStruct(i).rhoMean(128:end-1));
    rhoMed0=cat(2,rhoMed0,emsStruct(i).rhoMed(128:end-1));
    rhoStd0=cat(2,rhoStd0,emsStruct(i).rhoStd(128:end-1));
    
    tmpMean0=cat(2,tmpMean0,emsStruct(i).tmpMean(128:end-1));
    tmpMed0=cat(2,tmpMed0,emsStruct(i).tmpMed(128:end-1));
    tmpStd0=cat(2,tmpStd0,emsStruct(i).tmpStd(128:end-1));
    
    zmetMean0=cat(2,zmetMean0,emsStruct(i).zmetMean(128:end-1));
    zmetMeanW0=cat(2,zmetMeanW0,emsStruct(i).zmetMeanW(128:end-1));
    zmetMed0=cat(2,zmetMed0,emsStruct(i).zmetMed(128:end-1));
    zmetStd0=cat(2,zmetStd0,emsStruct(i).zmetStd(128:end-1));
%     
%     gamMean0=cat(2,gamMean0,emsStruct(i).gamMean(128:end-1));
%     gamMed0=cat(2,gamMed0,emsStruct(i).gamMed(128:end-1));
%     gamStd0=cat(2,gamStd0,emsStruct(i).gamStd(128:end-1));
    
    shellVol0=cat(2,shellVol0,emsStruct(i).shellVol(128:end-1));
    
end

rp=50:1:4e3;

rhoMean=interp1(rb0,rhoMean0,rp,'spline');
rhoMed=interp1(rb0,rhoMed0,rp,'spline');
rhoStd=interp1(rb0,rhoStd0,rp,'spline');

tmpMean=interp1(rb0,tmpMean0,rp,'spline');
tmpMed=interp1(rb0,tmpMed0,rp,'spline');
tmpStd=interp1(rb0,tmpStd0,rp,'spline');

zmetMean=interp1(rb0,zmetMean0,rp,'spline');
zmetMeanW=interp1(rb0,zmetMeanW0,rp,'spline');
zmetMed=interp1(rb0,zmetMed0,rp,'spline');
zmetStd=interp1(rb0,zmetStd0,rp,'spline');

% gamMean=interp1(rb0,gamMean0,rp,'spline');
% gamMed=interp1(rb0,gamMed0,rp,'spline');
% gamStd=interp1(rb0,gamStd0,rp,'spline');

shellVol=interp1(rb0,shellVol0,rp,'spline');

%% plot

cm=brewermap(6,'Set1');


%% density 
figure


h(1)=loglog(rp,rhoMean,'Color',cm(1,:),'DisplayName','Mean');
hold on
loglog(rp,rhoMean-0.5.*rhoStd,'--','Color',cm(1,:));
loglog(rp,rhoMean+0.5.*rhoStd,'--','Color',cm(1,:));
h(2)=loglog(rp,rhoMed,'Color',cm(2,:),'DisplayName','Median');

h(3)=loglog(r500.*1000.*[1 1],[1e-10 1e10],'--k','linewidth',1.8,...
    'DisplayName','$R_{500}$');

h(4)=loglog(r200.*1000.*[1 1],[1e-10 1e10],':k','linewidth',1.8,...
    'DisplayName','$R_{200}$');


%h(2)=loglog(emsStruct.rb,emsStruct(1).rhoMed,'DisplayName','Med');
xlim([50 4e3])
ylim([1e-7 1e-1])
hl=legend(h);
set(hl,'Interpreter','latex','Fontsize',12,'Location','SouthWest')
set(gca,'fontsize',12);
xlabelmine('$r\,[\mathrm{kpc}]$'); 
ylabelmine('$n\,[\mathrm{cm^{-3}}]$');
grid
titlemine('Density')

%% temperature 

figure


h(1)=loglog(rp,tmpMean,'Color',cm(1,:),'DisplayName','Mean');
hold on
loglog(rp,tmpMean-0.5.*tmpStd,'--','Color',cm(1,:));
loglog(rp,tmpMean+0.5.*tmpStd,'--','Color',cm(1,:));
h(2)=loglog(rp,tmpMed,'Color',cm(2,:),'DisplayName','Median');

h(3)=loglog(r500.*1000.*[1 1],[1e-10 1e10],'--k','linewidth',1.8,...
    'DisplayName','$R_{500}$');

h(4)=loglog(r200.*1000.*[1 1],[1e-10 1e10],':k','linewidth',1.8,...
    'DisplayName','$R_{200}$');


%h(2)=loglog(emsStruct.rb,emsStruct(1).rhoMed,'DisplayName','Med');
xlim([50 4e3])
ylim([5e6 1e8])
hl=legend(h);
set(hl,'Interpreter','latex','Fontsize',12,'Location','SouthWest')
set(gca,'fontsize',12);
xlabelmine('$r\,[\mathrm{kpc}]$'); 
ylabelmine('$T\,[\mathrm{K}]$');
grid
titlemine('Temperature')


%% metallicity 
figure


h(1)=loglog(rp,zmetMean,'Color',cm(1,:),'DisplayName','Mean');
hold on
loglog(rp,zmetMean-0.5.*zmetStd,'--','Color',cm(1,:));
loglog(rp,zmetMean+0.5.*zmetStd,'--','Color',cm(1,:));
h(2)=loglog(rp,zmetMed,'Color',cm(2,:),'DisplayName','Median');
h(3)=loglog(rp,zmetMeanW,'Color',cm(3,:),'DisplayName','Weighted Mean');

h(4)=loglog(r500.*1000.*[1 1],[1e-10 1e10],'--k','linewidth',1.8,...
    'DisplayName','$R_{500}$');

h(5)=loglog(r200.*1000.*[1 1],[1e-10 1e10],':k','linewidth',1.8,...
    'DisplayName','$R_{200}$');


%h(2)=loglog(emsStruct.rb,emsStruct(1).rhoMed,'DisplayName','Med');
xlim([50 4e3])
ylim([5e-3 0.2])
hl=legend(h);
set(hl,'Interpreter','latex','Fontsize',12,'Location','SouthWest')
set(gca,'fontsize',12);
xlabelmine('$r\,[\mathrm{kpc}]$'); 
ylabelmine('$Z\,[\mathrm{Z_\odot}]$');
grid
titlemine('Metallicity')

%% polytropic index 
gamMean=derive1(log10(tmpMean),log10(rhoMean))+1;
gamMed=derive1(log10(tmpMed),log10(rhoMed))+1;

rpp=rp(2:end-1);


figure


h(1)=semilogx(rpp,gamMed,'Color',cm(1,:),'DisplayName','Mean');
hold on
h(2)=semilogx(rpp,smooth(gamMed,50)+2,'Color',cm(2,:),'DisplayName','Median');




h(3)=semilogx(r500.*1000.*[1 1],[1e-10 1e10],'--k','linewidth',1.8,...
    'DisplayName','$R_{500}$');

h(4)=semilogx(r200.*1000.*[1 1],[1e-10 1e10],':k','linewidth',1.8,...
    'DisplayName','$R_{200}$');


%h(2)=loglog(emsStruct.rb,emsStruct(1).rhoMed,'DisplayName','Med');
xlim([50 4e3])
ylim([0.6 1.6])
hl=legend(h);
set(hl,'Interpreter','latex','Fontsize',12,'Location','SouthWest')
set(gca,'fontsize',12);
xlabelmine('$r\,[\mathrm{kpc}]$'); 
ylabelmine('$\gamma$');
grid
titlemine('polytropic index')



%% write to file 
global DEFAULT_PRINTOUT_DIR
global CLUSTER

%% density 
fname=[DEFAULT_PRINTOUT_DIR '/' CLUSTER '_mean_median_density_prof.txt'];
fid=fopen(fname,'w');

arr=cat(1,rp,rhoMean,rhoMed);

fprintf(fid,'columns:radius in kpc, mean of density values in cm^-3, median density  \n');
fprintf(fid,'%12.5e %12.5e %12.5e \n',arr);

fclose(fid);

%% density 
fname=[DEFAULT_PRINTOUT_DIR '/' CLUSTER '_mean_median_temperature_prof.txt'];
fid=fopen(fname,'w');

arr=cat(1,rp,tmpMean,tmpMed);

fprintf(fid,'columns:radius in kpc, mean of temperature values in K, median Temperature  \n');
fprintf(fid,'%12.5e %12.5e %12.5e \n',arr);

fclose(fid);

%% density 
fname=[DEFAULT_PRINTOUT_DIR '/' CLUSTER '_mean_median_metallicity_prof.txt'];
fid=fopen(fname,'w');

arr=cat(1,rp,zmetMean,zmetMeanW,zmetMed);

fprintf(fid,'columns:radius in kpc, mean of metallictiy values in solar metallicity, mass-weighted mean of metallictiy values, median metallicity  \n');
fprintf(fid,'%12.5e %12.5e %12.5e %12.5e \n',arr);

fclose(fid);

%% polytropic index
fname=[DEFAULT_PRINTOUT_DIR '/' CLUSTER '_mean_median_polytropicIndex_prof.txt'];
fid=fopen(fname,'w');

arr=cat(1,rpp,gamMed,smooth(gamMed,50)');

fprintf(fid,'columns:radius in kpc, polytropic index of based on median profiles, polytropic index smoothed over 50 kpc scale  \n');
fprintf(fid,'%12.5e %12.5e %12.5e \n',arr);

fclose(fid);
