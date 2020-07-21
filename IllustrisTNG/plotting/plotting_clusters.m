%% Script to make nice pictures of clusters
global illUnits
units;

snap=99; %z=0;
%% tng 100
bp=illustris.set_env('300');
%load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/tng100_z0_fofs_subs.mat')
%if readFlag
    fofs=illustris.groupcat.loadHalos(bp,snap);
%end
m200=fofs.Group_M_Mean200.*illUnits.massUnit;

m200Mask=m200>1e14;
idList=find(m200Mask)-1;

[~,ii]=sort(m200(idList+1),'descend');




%% plot
cnt=0;
for i=ii(1:10)
    cnt=cnt+1;
    fprintf('plotting no. %s out of %s \n',num2str(cnt),num2str(10))
    id=idList(i);
    
    %if readFlag
        gas=illustris.snapshot.loadHalo(bp, 99,id, 'gas');
    %end
    
    gas.newCoord=illustris.utils.centerObject(gas.Coordinates,fofs.GroupPos(:,id+1));
    gas=illustris.utils.addEntropy(gas);
    
    circ(1).radius=fofs.Group_R_Crit200(id+1);
    circ(1).width=1.5;
    circ(2).radius=fofs.Group_R_Crit500(id+1);
    circ(2).width=1.5;
    
    boxx=2.1.*fofs.Group_R_Crit200(id+1);
    vcm=fofs.GroupVel(:,id+1);
    ng=800;
    thk=fofs.Group_R_Crit500(id+1);
    %vDil=12;
    name=sprintf('fof%s_noVel',num2str(id));
    bmap='*PuOr';
    
    
    cubD=cell2grid(gas.newCoord,gas.Density,(gas.Masses./gas.Density).^(1/3),'box',boxx,'ext','ngrid',ng);
    
    
    illustris.plots.mkmapGas('gas',gas,'type','n','clims',[-5 -3],...
        'box',boxx,'ng',ng,'circ',circ,'brewer',bmap,'nobackground','savefig',name)
%         'vfield','vfieldcolor','k','dilute',vDil,'vcm',vcm,...
%         'print',name)
    
    illustris.plots.mkmapGas('gas',gas,'type','temperature','clims',[5 7.5],...
        'box',boxx,'ng',ng,'circ',circ,'brewer',bmap,'nobackground','savefig',name)
    %         'vfield','vfieldcolor','k','dilute',vDil,'vcm',vcm,...
%         'print',name)
    
    illustris.plots.mkmapGas('gas',gas,'type','entropy','clims',[1.5 3],...
        'box',boxx,'ng',ng,'circ',circ,'brewer',bmap,'nobackground','savefig',name)
%         'vfield','vfieldcolor','k','dilute',vDil,'vcm',vcm,...
%         'print',name)
    
end