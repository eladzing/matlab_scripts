%% read HI / H_2 data from Benedikt Deimer
% profiles

name='hih2_galaxy';
illustris.utils.read_catalog(name,'snap',99)

models=fliplr({'L08','GK11','K13','GD14','S14'});

%% select the right galaxies based on the
msk=3;
rgal=subs.SubhaloHalfmassRadType(illustris.partTypeNum('stars')+1,msk).*illUnits.lengthUnit;


%% generate a cumulative mass profile

%generalities
rp=hiStruct.profile_bins(:,msk);
area=hiStruct.profile_bins_area(:,msk);
vol=hiStruct.profile_bins_volume(:,msk);

%% total gas
gasDens2d=hiStruct.profile_gas_rho_2d(:,msk);
gasDens3d=hiStruct.profile_gas_rho_3d(:,msk);
gasDensMap=hiStruct.profile_gas_rho_map(:,msk);

gasMass2d=cumsum(gasDens2d.*area);
gasMass3d=cumsum(gasDens3d.*vol);
gasMassMap=cumsum(gasDensMap.*area);

%% neutral gas
neutRho2d=hiStruct.profile_f_neutral_H_2d(:,msk).*gasDens2d;
neutRho3d=hiStruct.profile_f_neutral_H_3d(:,msk).*gasDens3d;
neutRhoMap=hiStruct.profile_f_neutral_H_map(:,msk).*gasDensMap;

neutMass2d=cumsum(neutRho2d.*area);
neutMass3d=cumsum(neutRho3d.*vol);
neutMassMap=cumsum(neutRhoMap.*area);

%% molecular gas and atomic gas 
for k=1:length(models)
    prfName=['profile_f_mol_' models{k}];
    
    molRhoMap=hiStruct.([prfName '_map'])(:,msk).*neutRhoMap;
    molMassMap.(models{k})=cumsum(molRhoMap.*area);
    
    atmRhoMap=(1-hiStruct.([prfName '_map'])(:,msk)).*neutRhoMap;
    atmMassMap.(models{k})=cumsum(atmRhoMap.*area);
    
    if ~strcmp(models{k},'L08')
        molRho2d=hiStruct.([prfName '_3d'])(:,msk).*neutRho2d;
        molRho3d=hiStruct.([prfName '_2d'])(:,msk).*neutRho3d;
                
        molMass2d.(models{k})=cumsum(molRho2d.*area);
        molMass3d.(models{k})=cumsum(molRho3d.*vol);
        
        atmRho2d=(1-hiStruct.([prfName '_2d'])(:,msk)).*neutRho2d;
        atmRho3d=(1-hiStruct.([prfName '_3d'])(:,msk)).*neutRho3d;
        
        atmMass2d.(models{k})=cumsum(atmRho2d.*area);
        atmMass3d.(models{k})=cumsum(atmRho3d.*vol);
        
    end
end

%% plot atomic gas 

figure('color','w','position',[1090 331 766 602])
colors=brewermap(8,'Set1');
h=[];
for k=1:length(models)
        
    h(k)=plot(rp,atmMassMap.(models{k}),'-','color',colors(k,:),'linewidth',1.5,...
        'Displayname',models{k});
    hold on
    
    
    if ~strcmp(models{k},'L08')
        if k==1
            h(length(models)+1)=plot(rp,atmMass2d.(models{k}),'--','color',colors(k,:),'linewidth',1.5,...
                'DisplayName','2D');
           h(length(models)+2)=plot(rp,atmMass3d.(models{k}),':','color',colors(k,:),'linewidth',1.5,...
               'DisplayName','3D');
        else
            plot(rp,atmMass2d.(models{k}),'--','color',colors(k,:),'linewidth',1.5,...
                'DisplayName','2D');
            plot(rp,atmMass3d.(models{k}),':','color',colors(k,:),'linewidth',1.5,...
               'DisplayName','3D');
            
        end
        
        
    end
    
end
h(end+1)=plot(rp,neutMassMap,'-k','linewidth',1,...
                'DisplayName','NeutMap');
h(end+1)=plot(rp,neutMass2d,'--k','linewidth',1,...
                'DisplayName','Neut2D');
h(end+1)=plot(rp,neutMass3d,':k','linewidth',1,...
                'DisplayName','Neut3D');

plot([rgal rgal],[0 1e10],'-.k','linewidth',1.5)
            

            
set(gca,'fontsize',14)
hl=legend(h,'Interpreter','latex','fontsize',14,'location','Northwest');

xlabelmine('radius [kpc]');
ylabelmine('$M_\mathrm{HI}(<r)\,[\mathrm{M_\odot}]$');

% repeat only smaller 
figure('color','w','position',[1090 331 766 602])
colors=brewermap(8,'Set1');
h=[];
for k=1:length(models)
        
    h(k)=plot(rp,atmMassMap.(models{k}),'-','color',colors(k,:),'linewidth',1.5,...
        'Displayname',models{k});
    hold on
    
    
    if ~strcmp(models{k},'L08')
        if k==1
            h(length(models)+1)=plot(rp,atmMass2d.(models{k}),'--','color',colors(k,:),'linewidth',1.5,...
                'DisplayName','2D');
           h(length(models)+2)=plot(rp,atmMass3d.(models{k}),':','color',colors(k,:),'linewidth',1.5,...
               'DisplayName','3D');
        else
            plot(rp,atmMass2d.(models{k}),'--','color',colors(k,:),'linewidth',1.5,...
                'DisplayName','2D');
            plot(rp,atmMass3d.(models{k}),':','color',colors(k,:),'linewidth',1.5,...
               'DisplayName','3D');
            
        end
        
        
    end
    
end
h(end+1)=plot(rp,neutMassMap,'-k','linewidth',1,...
                'DisplayName','NeutMap');
h(end+1)=plot(rp,neutMass2d,'--k','linewidth',1,...
                'DisplayName','Neut2D');
h(end+1)=plot(rp,neutMass3d,':k','linewidth',1,...
                'DisplayName','Neut3D');

plot([rgal rgal],[0 1e10],'-.k','linewidth',1.5)
            
xlim([0 100])
ylim([0 2e9])
            
set(gca,'fontsize',14)
hl=legend(h,'Interpreter','latex','fontsize',14,'location','Northwest');

xlabelmine('radius [kpc]');
ylabelmine('$M_\mathrm{HI}(<r)\,[\mathrm{M_\odot}]$');








%% plot Molecular gas 

figure('color','w')
colors=brewermap(8,'Set1');
h=[];
for k=1:length(models)
        
    h(k)=plot(rp,molMassMap.(models{k}),'-','color',colors(k,:),'linewidth',1.5,...
        'Displayname',models{k});
    hold on
    
    
    if ~strcmp(models{k},'L08')
        if k==1
            h(length(models)+1)=plot(rp,molMass2d.(models{k}),'--','color',colors(k,:),'linewidth',1.5,...
                'DisplayName','2D');
           h(length(models)+2)=plot(rp,molMass3d.(models{k}),':','color',colors(k,:),'linewidth',1.5,...
               'DisplayName','3D');
        else
            plot(rp,molMass2d.(models{k}),'--','color',colors(k,:),'linewidth',1.5,...
                'DisplayName','2D');
            plot(rp,molMass3d.(models{k}),':','color',colors(k,:),'linewidth',1.5,...
               'DisplayName','3D');
            
        end
        
        
    end
    
end
h(end+1)=plot(rp,neutMassMap,'-k','linewidth',1,...
                'DisplayName','NeutMap');
h(end+1)=plot(rp,neutMass2d,'--k','linewidth',1,...
                'DisplayName','Neut2D');
h(end+1)=plot(rp,neutMass3d,':k','linewidth',1,...
                'DisplayName','Neut3D');

set(gca,'fontsize',14)
hl=legend(h,'Interpreter','latex','fontsize',14,'location','Northwest');

xlabelmine('radius [kpc]');
ylabelmine('$M_\mathrm{H_2}(<r)\,[\mathrm{M_\odot}]$');









%% interpolate to galactic radius