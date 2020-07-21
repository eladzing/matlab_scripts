list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
units ;
cf_ncf_list

boxx=1;
readFlag=false;

tTO=zeros(length(list),3);
vMean=zeros(length(list),2);
vMed=zeros(length(list),2);
sk=zeros(size(list));
mvir=zeros(size(list));
vvir=zeros(size(list));
rvir=zeros(size(list));
tvir=zeros(size(list));
for i=1:length(list)
    new_env(list(i))
    global VCM
    global CLUSTER
    
    if readFlag
        vx=Vx(boxx)-VCM(1);
        vy=Vy(boxx)-VCM(2);
        vz=Vz(boxx)-VCM(3);
        vv=sqrt(vx.^2+vy.^2+vz.^2);
        ro=RHOG(boxx);
    end
    
    %% get mask
    rc=mk_radius_cube(boxx);
    rv=get_rvir;
    rvir(i)=rv;
    mvir(i)=get_mvir;
    vvir(i)=get_vvir;
    tvir(i)=get_tvir;
    
    continue
    mask=rc<=0.2.*rv & rc>0.01.*rv;
    
    %% smooth
    
    vs1=smooth3(vv);
    %vs2=smooth3(vv,'box',7);
    %vs3=smooth3(vv,'gaussian',7,3);
    
    vHist=vv(mask);
    %vHist1=vs1(mask);
    %vHist2=vs2(mask);
    %vHist3=vs3(mask);
    
    global RELAXED
    if RELAXED
        tag='relaxed';
    else
        tag='unrelaxed';
    end

        
    figure
    %subplot(2,2,1)
    hist(vHist,100);
    xlabelmine('$v\,[\mathrm{km/sec}]$')
    titlemine(sprintf('%s, %s',CLUSTER,tag))
    %subplot(2,2,2)
    %hist(vHist1,100);
    %subplot(2,2,3)
    %hist(vHist2,100)
    %subplot(2,2,4)
    %hist(vHist3,100)
    
    vMean(i,1)=mean(vHist);
    vMean(i,2)=sum(sum(sum(ro(mask).*vv(mask))))/sum(sum(sum(ro(mask))));
    vMed(i)=median(vHist);
    sk(i)=skewness(vHist);
    
    tTO(i,1)=((0.2*rv*Mpc)/(vMean(i,1)*km))/yr/1e9;
    tTO(i,2)=((0.2*rv*Mpc)/(vMean(i,2)*km))/yr/1e9;
    tTO(i,3)=((0.2*rv*Mpc)/(vMed(i)*km))/yr/1e9;
    %tTO1=((0.2*rv*Mpc)/(mean(vHist1)*km))/yr/1e9
    %tTO2=((0.2*rv*Mpc)/(mean(vHist2)*km))/yr/1e9
    %tTO3=((0.2*rv*Mpc)/(mean(vHist3)*km))/yr/1e9
    
end