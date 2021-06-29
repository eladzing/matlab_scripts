%% exact disk rps model
clk0=clock;
Mc=1e15;
sigma=1e7./sigma_factor('50-90');

r=0.01:0.01:60; % in units of the scale radius


%% exploring Msat and rd
alfa=1; % fudge factor
fc=0.1;   % cluster gas fraction
%fd=0.1;  % disk gas fraction
%Mc=1e15;  % cluster mass
delv = 337; % delta_vir
zred=0;

nn=100;

bt_lim=log10([0.1 10]);
fd_lim=[0.02 0.5];

lbt=bt_lim(1):diff(bt_lim)/(nn-1):bt_lim(2);
beta=10.^lbt;
fd=fd_lim(1):diff(fd_lim)/(nn-1):fd_lim(2);

%% calculate pval coeeficients except for fd dependence
etap=1;  % position in cluster in units of Rv
pval0=rps_factor_expdisk('alpha',alfa,'fc',fc,'fd',1,'Mc',Mc,...
    'sigma',sigma,'delv',delv,'zred',zred,'etap',etap);

%% generate pval matrix
[~, fdm]=meshgrid(beta,fd);

pval=pval0./fdm;
clear bta fdm

mstrip=zeros(nn);
rpos90=zeros(nn);



for i=1:length(fd)
    for j=1:length(lbt)
        
        
        mf1=exp_disk_mass(r,beta(j));
        
        fg1=disk_force_reduced(r,'beta',beta(j),'fg',fd(i));
        
        [fgmax,ind1]=max(fg1);
        fg=fg1(ind1:end);
        mf=mf1(ind1:end);
        
        
        fg_mf01=interp1(mf(mf<0.5),fg(mf<0.5),0.1); % value of fg for mf=0.1
        
        if (pval(i,j)<=fgmax)
            mstrip(i,j)=interp1(fg,mf,pval(i,j));
        else
            mstrip(i,j)=0;
        end
        
        
        %calculate map of position in cluster for 90% stripping
        rpos90(i,j)=log10(1./sqrt(fg_mf01./pval(i,j)));
        
        
    end
end

% mass rfaction map
figure
imagesc(lbt,fd,100.*(1-mstrip))
set(gca,'Ydir','normal');
bar=colorbar;
set(get(bar,'Title'),'String','$\%M_{strip}$','Fontsize',12,'Interpreter','latex');
caxis([0 100])

ylabelmine('$f_{g,disk}$')
xlabelmine('$\log(\beta)$')

printout_fig(gcf,'rpsdisk_mstrip_alfabeta','v')
%
%rpos map

figure
imagesc(lbt,fd,rpos90)
set(gca,'Ydir','normal');
bar1=colorbar;
set(get(bar1,'Title'),'String','$r_{pos}\,[R_c]$','Fontsize',12,'Interpreter','latex');
caxis([-2 1])

ylabelmine('$f_{g,disk}$')
xlabelmine('$\log(\beta)$')

printout_fig(gcf,'rpsdisk_rpos_alfabeta','v')

clk=clock-clk0


