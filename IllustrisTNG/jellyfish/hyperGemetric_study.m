NN=20;

kkk=[0.05 0.2 0.25 0.4 0.5 0.6 0.75 0.8 0.95].*NN;

cmap=brewermap(256,'reds');
cmap(1,:)=[0.65 0.65 0.65];
myFigure;
t=tiledlayout(3,3);

h=[];
for k=1:length(kkk)
    KK=kkk(k);
    
    %mapp=-10.*ones(NN+1,NN);
    mapp=nan(NN+1,NN);
    nn=1:NN;
    
    for i=1:length(nn)
        nrm=nchoosek(NN,nn(i));
        
        kk=0:min(nn(i),KK);
        
        for j=1:length(kk)
            try
                mapp(j,i)=nchoosek(KK,kk(j))*nchoosek(NN-KK,nn(i)-kk(j))/nrm;
            catch
                %fprintf('nn=%i, kk=%i \n',nn(i),kk(j));
            end
            
        end
        
    end
    
    h(end+1)=nexttile;
    imagesc(nn,(0:NN),mapp);
    caxis([0 1])
    colormap(cmap)
    %colorbar
    set(gca,'YDir','normal')
    hold on
    plot([0 NN],[0 NN],':k')
    plot(nn,nn.*kkk(k)./NN,'--k')
    grid
    
    
    titlemine(['K=' num2str(kkk(k))],10);
    
    if kkk(k)/NN==0.5
        mapp50=mapp;
    end
end
hb=colorbar(h(end));
hb.Layout.Tile='east';
ylabel(t,'k','Interpreter','latex','fontsize',14);
xlabel(t,'n','Interpreter','latex','fontsize',14);
 title(t,['N=' num2str(NN)],'Interpreter','latex')  
t.Padding='compact';
t.TileSpacing='compact';

%% settign NN and nn constant and changing K and k 

myFigure
t=tiledlayout(2,4);
cmap=brewermap(256,'reds');
cmap(1,:)=[0.65 0.65 0.65];
h=[];
for NN=30:10:100

    nn=20;

KK=0:NN;
kk=0:nn;


mapp=nan(length(kk),length(KK));

nrm=nchoosek(NN,nn);

for i=1:length(KK)
    
    NNKK=NN-KK(i);
    
    for j=1:length(kk)
        try
            mapp(j,i)=nchoosek(KK(i),kk(j))*nchoosek(NNKK,nn-kk(j))/nrm;
        catch
            %fprintf('nn=%i, kk=%i \n',nn(i),kk(j));
        end
        
    end
    
end


h(end+1)=nexttile;
imagesc(KK./NN,kk./nn,mapp);
caxis([0 1])
colormap(cmap)
set(gca,'YDir','normal')
hold on
plot([0 1],[0 1],':k')
grid
titlemine(['N=' num2str(NN)],10);
end

hb=colorbar(h(end));
hb.Layout.Tile='east';

colorbar
set(gca,'YDir','normal')
hold on


ylabel(t,'k','Interpreter','latex');
xlabel(t,'K','Interpreter','latex');
 title(t,['n=' num2str(nn)],'Interpreter','latex')  
t.Padding='compact';
t.TileSpacing='compact';