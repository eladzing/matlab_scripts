% plot results of galaxy catalog 

%iend=find(diff(galRes(1).rpos)>0,1,'first');
iend=length(galRes(1).time);
edj=0:0.1:1.5;
bCen=edj(1:end-1)+0.5.*diff(edj);

 %qf=zeros(size(galRes(1).ssfr(1:iend)));
 qf=zeros(1,length(bCen));
 cnt=qf;
for i=1:length(galRes)
    rpp=galRes(i).rpos(1:iend)./host.Rvir;
    rpp2=rpp.*generate_projectionFac(length(rpp))';
    bInd=discretize(rpp,edj);
    bInd2=discretize(rpp2,edj);

 qn=galRes(i).ssfr(1:iend)<1e-11;

 qf(bInd)=qf(bInd)+qn;
 cnt(bInd)=cnt(bInd)+2*ones(size(qn));

end


qf1=zeros(1,length(bCen));
qf2=zeros(1,length(bCen));
cnt1=qf1;
cnt2=qf2;
for i=1:length(galRes)
    rpp=galRes(i).rpos(1:iend)./host.Rvir;
    rpp2=rpp.*generate_projectionFac(length(rpp))';
    bInd=discretize(rpp,edj);
    bInd2=discretize(rpp2,edj);
    
        
    qn=galRes(i).ssfr(1:iend)<1e-11;
    onn=ones(size(qn));
    
    qff1=zeros(1,length(bCen));
    cntt1=qff1;
    qff2=qff1;
    cntt2=qff1;
    
    for j=1:length(bCen)
        qff1(j)=sum(qn(bInd==j));
        cntt1(j)=sum(onn(bInd==j));
        
        qff2(j)=sum(qn(bInd2==j));
        cntt2(j)=sum(onn(bInd2==j));
    end
    
    qf1=qf1+qff1;
    qf2=qf2+qff2;
    cnt1=cnt1+cntt1;
    cnt2=cnt2+cntt2;
    
end



figure


    semilogy(,...
        galRes(i).stellarMass(1:iend)./galRes(i).stellarMass(1));
    hold on
end
 qf=zeros(size(rp));  
for i=1:length(galRes)

    qn=galRes(i).ssfr(1:iend)<1e-11;
    
    qf=qf+qn;
end
    
rp=galRes(1).rpos(1:iend)./host.Rvir;

