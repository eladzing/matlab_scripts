

xlen =100;
ylen=100;

pln=zeros(ylen,xlen);

fgs=linspace(-5,0.5,ylen);
smass=linspace(9,11.5,xlen);
%% fit parameters 

ft = read_fgs_gaussian_fit_files('TNG100',99,200);

% mu1= [-0.1986 -0.2985 -0.3917 -0.4902 -0.5952 -0.6436 -0.6338 -1.3622 -1.7616 -1.2761];      %-1.2869 -1.0945];
% sig1= [0.2844 0.2636 0.2594 0.2630 0.2973 0.3063 0.3042 0.5076 0.0419 0.4343];               %0.3019 0.6537];
% mu2= [-4.0717 -4.7911 -5.0096 -4.1010 -3.1496 -3.1897 -2.6386 -2.4671 -1.5513 -1.9031];      %-1.2225 7.0520];
% sig2= [0.2416 0.2381 0.2526 0.2488 0.6008 0.5578 0.7994 0.3275 0.4220 0.2450];               %0.9018 0.2094];
% a1= [345.4941 290.3323 174.3282 160.9878 224.9213 125.2811 41.8338 22.5904 154.2026 136.4970]; %8.7071 3.9586];
% a2= [0 0 0 0 9.3505 50.4207 40.1317 19.8226 264.0177 -22.2368];                                %6.7216 0];
% 
% binEdge=[9 9.2500 9.5000 9.7500 10 10.2500 10.5000 10.7500 11 11.2500 11.5000];          


%smass=9:0.05:11.5;
binInd=discretize(smass,ft.binEdge);

for i=1:length(binInd)
    
    ii=binInd(i);

   yy=ft.a1(ii).*exp(-0.5.*((fgs-ft.mu1(ii))./ft.sig1(ii)).^2) + ...
        ft.a2(ii).*exp(-0.5.*((fgs-ft.mu2(ii))./ft.sig2(ii)).^2);
     pln(:,i)=yy./sum(yy);
end


% hf=figure('color','w','position',[823 392 1038 904]);
% for i=1:10
%     inds=find(smass>=binEdge(i) & smass<=binEdge(i+1));
%     s=surf(smass(inds),fgs,pln(:,inds),'facecolor',cc(i,:));
%     s.EdgeAlpha=0.2;
%     if i<5
%         s.EdgeColor='w';
%     else
%         s.EdgeColor='k';
%     end
%     
%     s.MeshStyle='column';
%     if i==1;hold on;end
% end


cc=brewermap(10,'*Spectral');%(magma(10));%

pln0=zeros(size(pln))-0.0001;

hf=figure('color','w','position',[823 392 1038 904]);
maxx=max(max(pln));
for i=1:10
    
    
    inds=find(smass>=ft.binEdge(i) & smass<=ft.binEdge(i+1));
    pln1=pln0;
    pln1(:,inds)=pln1(:,inds)+pln(:,inds)+0.0001;
    
    switch i
        case 1
            pln1(:,1)=-0.0001;
             s=surf(smass,fgs,pln1,'facecolor',cc(i,:));
            %s.FaceAlpha=0.75;
%         case 10           
%             s=surf(smass(inds),fgs,pln1(:,inds),'facecolor',cc(i,:));
        otherwise
            s=surf(smass,fgs,pln1,'facecolor',cc(i,:));
            
    end
    s.EdgeAlpha=0.3;
    
    switch(i)
        case {2,3,4,5,6,7,8}
            s.EdgeColor='k';
        otherwise
            s.EdgeColor='w';
    end
    
    s.MeshStyle='column';
    if i==1;hold on;end
end

hl=xlabelmine('log Stellar Mass $[\mathrm{M_\odot}]$',18);
set(hl,'Rotation',22);
%ylabelmine('$\log\,f_\mathrm{gs}$',16);
ylabelmine('log gas/stars mass ratio',18);
zlabelmine('PDF',18);
set(gca,'fontsize',16)

xlim([9 11.5])
ylim([-5 0.5])
zlim([-0.00005 1.058*maxx])
%view(gca,fliplr([64.043849716932456,10.156841357681486]));
    %-38.9431034482758 13.5248304084903]);
grid(gca,'on');
% Set the remaining axes properties
set(gca,'FontSize',14,'YTick',[-5 -4 -3 -2 -1 0]);

%cmap=brewermap(256,'YlOrRd');
%cmap=viridis(256);colormap(flipud(cmap))
%colormap(cmap)
