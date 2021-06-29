%% create a movie of hammer plots of flux from out to in
global NCELL
global hub
global zred

maxInd=NCELL*5/2-1;

[~,mesh_phi,mesh_theta] = sphere_grid(8);

% thik is thickness of shell
thik=0;
step=2*thik+1;
inner=30;

outInd=maxInd-thik;
inInd=inner+thik;

readcube=false(1,4);
fluxnorm=(0.056.*(get_mvir./1e13).^0.15*(1+zred)^2.25)./(4.*pi); % normalization for flux
clim=[-10 10];


nframe=outInd-inInd+1;
F(nframe) = struct('cdata',[],'colormap',[]);
k=0;
for i=outInd:-2:inInd
    k=k+1;
    if i<NCELL
        boxx=1;ii=1;
        if ~readcube(ii)
            cu=flux_sphere(boxx);
            cus=S_sphere(boxx);
            readcube(ii)=true;
        end
        
        ind=i;
        
    elseif i <NCELL*3/2
        boxx=2;ii=2;
        if ~readcube(ii)
            cu=flux_sphere(boxx);
            cus=S_sphere(boxx);
            readcube(ii)=true;
        end
        
        ind=i-NCELL/2;
        
    elseif i < NCELL*2
        boxx=4;ii=3;
        if ~readcube(ii)
            cu=flux_sphere(boxx);
            cus=S_sphere(boxx);
            readcube(ii)=true;
        end
        ind=i-NCELL;
        
    else
        boxx=8;ii=4;
        if ~readcube(ii)
            cu=flux_sphere(boxx);
            cus=S_sphere(boxx);
            readcube(ii)=true;
        end
        ind=i-NCELL*3/2;
        
    end
    
    
    
    shell=squeeze(sum(cu(ind-thik:ind+thik,:,:),1))./fluxnorm;
    shells=squeeze(sum(cus(ind-thik:ind+thik,:,:),1));
    mphi=squeeze(mesh_phi(ind,:,:));
    mthet=squeeze(mesh_theta(ind,:,:));
    rad=ind/NCELL*0.5*boxx/(1+zred)/hub/get_rvir;
    
   
    figure('Position',[100 100 600 700])
   
    hs1=subplot('Position',[0.03 0.52 0.94 0.45]);
    plot_hammer_shell(shell,mthet,mphi,'nc',500);
    titlemine('Mass Flux')
    colorbar
    caxis(clim)
    text(0.82,0.93,sprintf('r/R_v=%3.3g',rad),'units','normalized')
    axis fill
    
    
    hs2=subplot('Position',[0.03 0.03 0.94 0.45]);
    plot_hammer_shell(log10(shells),mthet,mphi,'nc',500);
    titlemine('Entropy')
    colorbar
    caxis([-1.5 0])
    axis fill
    
    F(k)=getframe(gcf);
   
    close(gcf)
    
    
end
%  [h,w,p]=size(F(1).cdata);
%  hf=figure;
%  set(hf,'Position',[150 150 w h]);
 
 
