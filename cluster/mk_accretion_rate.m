%% calculate approximate accretion rate
load('C:\Users\eladzing\Documents\cluster\matlab\mat_files\shockedge.mat');
list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

mv1=zeros(size(list));
mv2=zeros(size(list));
edgMax=zeros(size(list));
edgMin=zeros(size(list));
edgMean=zeros(size(list));
edgMid=zeros(size(list));


a1=1;
a2=0.626;
for i=1:length(list)
    
    new_env(list(i),'a1')
    
    m1=get_mvir();
    m2=get_mvir();
    rv=get_rvir();
    mv1(i)=m1;
    fac=1; %m2/m1
    new_env(list(i),'a06')
    
    mv2(i)=get_mvir();
    
    e1=shockedge_a1{i,4}./rv;%./rv_a1(i);
    e2=shockedge_a1{i,5}./rv;%./rv_a1(i);
    e1=e1(shockedgeMask1_a1(i,:));
    e2=e2(shockedgeMask2_a1(i,:));
    edg=cat(2,e1,e2).*fac;%(m2/m1);
    
    edgMax(i)=max(edg);
    edgMin(i)=min(edg);
    edgMean(i)=mean(edg);
    edgMid(i)=0.5.*(max(edg)+min(edg));
    
end

gamma1=log10(mv1./mv2)./log10(a1/a2);

[gamma,ind]=sort(gamma1);

h=[];
plot(gamma,edgMax(ind),'-o',gamma,edgMin(ind),'-+',...
    gamma,edgMean(ind),'-d',gamma,edgMid(ind),'-s')

hl=legend('Max','Min','Mean','Mid');
set(hl,'Location','SouthEast','Interpreter','latex','Fontsize',14);
xlabelmine('$\Gamma$');
ylabelmine('Edge');


