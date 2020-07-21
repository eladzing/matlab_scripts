



cmap=brewermap(8,'Dark2');

tag={'$\log\,M_\mathrm{h}=11.9-12.4$',...
    '$\log\,M_\mathrm{h}=12.4-12.8$',...
    '$\log\,M_\mathrm{h}=12.8-13.3$',...
    '$\log\,M_\mathrm{h}=13.3-13.8$',...
    '$\log\,M_\mathrm{h}=113.8-14.3$',...
    '$\log\,M_\mathrm{h}=14.3-14.8$',...
    'Centrals'};


% tag={'$\log\,M_\mathrm{h}=12-12.5$',...
%     '$\log\,M_\mathrm{h}=12.5-13$',...
%     '$\log\,M_\mathrm{h}=13-13.5$',...
%     '$\log\,M_\mathrm{h}=13.5-13$',...
%     '$\log\,M_\mathrm{h}=14-14.5$',...
%     '$\log\,M_\mathrm{h}=14.5-15$',...
%     'Centrals'};

fld={'mh12' 'mh125' 'mh13' 'mh135' 'mh14' 'mh145' 'centrals'};

figure
h=[];
for k=1:length(fld)-1

    
ms=fig3a.(sprintf('ms_qf_%s',fld{k}))(2,:);
qf=fig3a.(sprintf('ms_qf_%s',fld{k}))(1,:);
ypos=abs(fig3a.(sprintf('ms_qf_%s_errP',fld{k}))(1,:)-qf);
yneg=abs(fig3a.(sprintf('ms_qf_%s_errM',fld{k}))(1,:)-qf);

h(k)=errorbar(ms,qf,yneg,ypos,':','color',cmap(k,:),...
    'Displayname',tag{k});
if k==1;hold on;end
end

hl=legend(h);
set(hl,'Fontsize',14,'interpreter','latex')

xlabelmine('log Stellar Mass $[\mathrm{M_\odot}]$');
ylabelmine('Quenched Fration');


%%

fld={'ms97' 'ms101' 'ms105' 'ms109'};
figure
h=[];
for k=1:length(fld)

    
ms=fig3b.(sprintf('mh_qf_%s',fld{k}))(2,:);
qf=fig3b.(sprintf('mh_qf_%s',fld{k}))(1,:);
ypos=abs(fig3b.(sprintf('mh_qf_%s_errP',fld{k}))(1,:)-qf);
yneg=abs(fig3b.(sprintf('mh_qf_%s_errM',fld{k}))(1,:)-qf);

h(k)=errorbar(ms,qf,yneg,ypos,'color',cmap(k,:),...
    'Displayname',tag{k});
if k==1;hold on;end
end

hl=legend(h);
set(hl,'Fontsize',14,'interpreter','latex')

xlabelmine('log halo Mass $[\mathrm{M_\odot}]$');
ylabelmine('Quenched Fration');




%%
fld={'mh14','mh13','mh12'};
figure
h=[];
for k=1:length(fld)

    
ms=fig4a.(sprintf('rp_qf_ms97_%s',fld{k}))(2,:);
qf=fig4a.(sprintf('rp_qf_ms97_%s',fld{k}))(1,:);
ypos=abs(fig4a.(sprintf('rp_qf_ms97_%s_errP',fld{k}))(1,:)-qf);
yneg=abs(fig4a.(sprintf('rp_qf_ms97_%s_errM',fld{k}))(1,:)-qf);

h(k)=errorbar(ms,qf,yneg,ypos,'color',cmap(k,:),...
    'Displayname',tag{k});
if k==1;hold on;end
end

hl=legend(h);
set(hl,'Fontsize',14,'interpreter','latex')

xlabelmine('log rpos $[\mathrm{R_{vir}}]$');
ylabelmine('Quenched Fration');

%%

fld={'ms109','ms105','ms101','ms97'};
figure
h=[];
for k=1:length(fld)

    
ms=fig4b.(sprintf('rp_qf_mh12_%s',fld{k}))(2,:);
qf=fig4b.(sprintf('rp_qf_mh12_%s',fld{k}))(1,:);
ypos=abs(fig4b.(sprintf('rp_qf_mh12_%s_errP',fld{k}))(1,:)-qf);
yneg=abs(fig4b.(sprintf('rp_qf_mh12_%s_errM',fld{k}))(1,:)-qf);

h(k)=errorbar(ms,qf,yneg,ypos,'color',cmap(k,:),...
    'Displayname',tag{k});
if k==1;hold on;end
end

hl=legend(h);
set(hl,'Fontsize',14,'interpreter','latex')

xlabelmine('log rpos $[\mathrm{R_{vir}}]$');
ylabelmine('Quenched Fration');


%%
fig3b.mh_qf_ms97=mh_qf_ms97;
fig3b.mh_qf_ms97_errM=mh_qf_ms97_errM(1,:);
fig3b.mh_qf_ms97_errP=mh_qf_ms97_errP(1,:);
fig3b.mh_qf_ms101=mh_qf_ms101;
fig3b.mh_qf_ms101_errM=mh_qf_ms101_errM(1,:);
fig3b.mh_qf_ms101_errP=mh_qf_ms101_errP(1,:);
fig3b.mh_qf_ms105=mh_qf_ms105;
fig3b.mh_qf_ms105_errM=mh_qf_ms105_errM(1,:);
fig3b.mh_qf_ms105_errP=mh_qf_ms105_errP(1,:);
fig3b.mh_qf_ms109=mh_qf_ms109;
fig3b.mh_qf_ms109_errM=mh_qf_ms109_errM(1,:);
fig3b.mh_qf_ms109_errP=mh_qf_ms109_errP(1,:);
fig3b.msGroups={'9.7-10.1',...
    '10.1-10.5','10.5-10.9','10.9-11.3'};

fig3a.mhGroups={'12-12.5',...
    '12.5-13',...
    '13-13.5',...
    '13.5-13',...
    '14-14.5',...
    '14.5-15',...
    'Centrals'};


fig3a.ms_qf_centrals=ms_qf_centrals;
fig3a.ms_qf_centrals_errM=ms_qf_centrals_errM(1,:);
fig3a.ms_qf_centrals_errP=ms_qf_centrals_errP(1,:);
fig3a.ms_qf_mh12=ms_qf_mh12;
fig3a.ms_qf_mh12_errM=ms_qf_mh12_errM(1,:);
fig3a.ms_qf_mh12_errP=ms_qf_mh12_errP(1,:);
fig3a.ms_qf_mh125=ms_qf_mh125;
fig3a.ms_qf_mh125_errM=ms_qf_mh125_errM(1,:);
fig3a.ms_qf_mh125_errP=ms_qf_mh125_errP(1,:);
fig3a.ms_qf_mh13=ms_qf_mh13;
fig3a.ms_qf_mh13_errM=ms_qf_mh13_errM(1,:);
fig3a.ms_qf_mh13_errP=ms_qf_mh13_errP(1,:);
fig3a.ms_qf_mh135=ms_qf_mh135;
fig3a.ms_qf_mh135_errM=ms_qf_mh135_errM(1,:);
fig3a.ms_qf_mh135_errP=ms_qf_mh135_errP(1,:);
fig3a.ms_qf_mh14=ms_qf_mh14;
fig3a.ms_qf_mh14_errM=ms_qf_mh14_errM(1,:);
fig3a.ms_qf_mh14_errP=ms_qf_mh14_errP(1,:);
fig3a.ms_qf_mh145=ms_qf_mh145;
fig3a.ms_qf_mh145_errM=ms_qf_mh145_errM(1,:);
fig3a.ms_qf_mh145_errP=ms_qf_mh145_errP(1,:);


fig4a.rp_qf_ms97_mh12=rp_qf_ms97_mh12;
fig4a.rp_qf_ms97_mh12_errM=rp_qf_ms97_mh12_errM;
fig4a.rp_qf_ms97_mh12_errP =rp_qf_ms97_mh12_errP;
fig4a.rp_qf_ms97_mh13=rp_qf_ms97_mh13;
fig4a.rp_qf_ms97_mh13_errM=rp_qf_ms97_mh13_errM;
fig4a.rp_qf_ms97_mh13_errP =rp_qf_ms97_mh13_errP;
fig4a.rp_qf_ms97_mh14=rp_qf_ms97_mh14;
fig4a.rp_qf_ms97_mh14_errM=rp_qf_ms97_mh14_errM;
fig4a.rp_qf_ms97_mh14_errP =rp_qf_ms97_mh14_errP;
fig4a.msGroup='9.7-10.5';
fig4a.mhGroups={'12.3-13.2','13.2-14.2','14.1-15'};

fig4b.rp_qf_mh12_ms97=rp_qf_mh12_ms97;
fig4b.rp_qf_mh12_ms97_errM=rp_qf_mh12_ms97_errM;
fig4b.rp_qf_mh12_ms97_errP=rp_qf_mh12_ms97_errP;
fig4b.rp_qf_mh12_ms101=rp_qf_mh12_ms101;
fig4b.rp_qf_mh12_ms101_errM=rp_qf_mh12_ms101_errM;
fig4b.rp_qf_mh12_ms101_errP=rp_qf_mh12_ms101_errP;
fig4b.rp_qf_mh12_ms105=rp_qf_mh12_ms105;
fig4b.rp_qf_mh12_ms105_errM=rp_qf_mh12_ms105_errM;
fig4b.rp_qf_mh12_ms105_errP=rp_qf_mh12_ms105_errP;
fig4b.rp_qf_mh12_ms109=rp_qf_mh12_ms109;
fig4b.rp_qf_mh12_ms109_errM=rp_qf_mh12_ms109_errM;
fig4b.rp_qf_mh12_ms109_errP=rp_qf_mh12_ms109_errP;
fig4b.mhGroup='>12.5';
fig4b.msGroups={'9.7-10.1','10.1-10.5','10.5-10.9','10.9-10.11'};