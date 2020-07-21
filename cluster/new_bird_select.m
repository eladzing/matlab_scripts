

%% create parameter plot which defines selection

global NCELL
global hub

par1=T(boxx);lim1=[2 9];
par2=RHOG(boxx);lim2=[7 17];

len=[500 500];

vol= (boxx/hub/NCELL)^3;
mass=par2.*vol;

[bird binx biny]=basic_bird(log10(par1),log10(par2),mass,ones(size(mass)),lim1,lim2,len);

clear mass

masshist=squeeze(bird(:,:,1)./sum(sum(bird(:,:,1))));

%% draw histogram

figure;
subplot(2,2,2);

imagesc(lim1,lim2,log10(masshist));
load('MyColormaps','avijet_bird')
set(gcf,'Colormap',avijet_bird)
set(gcf,'Colormap',avijet_bird)
set(gca,'Ydir','normal');
xlabel('log T','Fontsize',14);
ylabel('log \rho','Fontsize',14);
set(gca,'fontsize',12)
bar=colorbar;
set(get(bar,'Title'),'String','$m/M_{tot}$','Fontsize',12,'Interpreter','latex');

%% select points to map
[xv, yv] = getline(gcf); % uncomment this line for interactive selection
line(xv,yv, 'Color','g');

par2_log = log10(par2);
par1_log = log10(par1);

mask = inpolygon(par1_log,par2_log,xv,yv); 

%% plot the maps

% plot xz

subplot(2,2,1)
mkmap_bare(par2.*mask,'weight',mask,'thick',NCELL/2,'projection','zx','box',boxx,'circles',get_rvir())

% plot xy
subplot(2,2,3)
mkmap_bare(par2.*mask,'weight',mask,'thick',NCELL/2,'projection','xy','box',boxx,'circles',get_rvir())

% plot xz
subplot(2,2,4)
mkmap_bare(par2.*mask,'weight',mask,'thick',NCELL/2,'projection','yz','box',boxx,'circles',get_rvir())
%bar=colorbar;
%set(get(bar,'Title'),'String','$\rho$','Fontsize',12,'Interpreter','latex');






