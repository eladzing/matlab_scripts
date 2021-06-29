edge=[];
shockedge=shockedge_a1;
for i=1:16
    edge(i,:)=shockedge{i,4};
    edge2(i,:)=shockedge{i,5};
    rv(i)=shockedge{i,2}';
    mv(i)=shockedge{i,3}';
    if(length(edge)<2)
        edge=cat(2,edge,edge);
    end
    if(length(edge2)<2)
        edge2=cat(2,edge2,edge2);
    end
    edger(i,:)=edge(i,:)./rv(i);
    edger2(i,:)=edge2(i,:)./rv(i);
    
    rang(i)=diff(edger(i,:));
    rang2(i)=diff(edger2(i,:));
    
end
figure
errorbar(mv./1e14,mean(edger,2),rang,'.')
xlabelmine('$M_{vir}\,[10^{14}\mathrm{M_\odot}]$')
ylabelmine('$R_{edge}/R_{vir}$')
titlemine('Edge by $\nabla S$, $z=0.6$')
%exportfig(gcf,'edge1_vs_mv_a06.png','format','png');
%exportfig(gcf,'edge1_vs_mv_a06.eps')

figure 
errorbar(mv./1e14,mean(edger2,2),rang,'.')
xlabelmine('$M_{vir}\,[10^{14}\mathrm{M_\odot}]$')
ylabelmine('$R_{edge}/R_{vir}$')
titlemine('Edge by $S_{max}$, $z=0.6$')
%exportfig(gcf,'edge2_vs_mv_a06.png','format','png');
%exportfig(gcf,'edge2_vs_mv_a06.eps')
