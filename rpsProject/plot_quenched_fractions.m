% plot results of galaxy catalog 

%iend=find(diff(galRes(1).rpos)>0,1,'first');
iend=length(galRes(1).rpos);
%iend=6599;
edj=0:0.1:1.2;
bCen=edj(1:end-1)+0.5.*diff(edj);

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

lw=1.5;
h=[];
figure
h(1)=plot(bCen,qf1./cnt1,'bo-',...
    'linewidth',lw,'Displayname','non-projected');
hold on
h(2)=plot(bCen,qf2./cnt2,'rx-',...
       'linewidth',lw,'Displayname','projected');

h(3)=plot(Wetzel12_fig15(2,:),Wetzel12_fig15(1,:),'sk',...
'DisplayName','Wetzel+12','markersize',8);

ylim([0 1])
grid

hl=legend(h);
set(hl,'Interpreter','latex','fontsize',14,'Location','NorthEast');

xlabelmine('$r/R_{\mathrm{200,m}}$',16);
ylabelmine('quenched fraction',16);
set(gca,'fontsize',14)



