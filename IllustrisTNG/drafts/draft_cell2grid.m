function res=cell2grid(coords,vals,varargin)...weights,Ngrid)
%% put cells into uniform grid 

len=length(vals);


%% defaults 
Ngrid=256;
wt=ones(1,len);
cellNum=len;


%% parse arguments 
i=1
while i<length(varargin)
    switch lower(varargin{i})
        case{'nc','ncell','ngrid'
               






%% find box side length

m1=abs(max(coords,[],2));
m2=abs(min(coords,[],2));
ll=max(cat(1,m1,m2)); %in ckpc/h;

boxSide=ceil(2.*ll); %in kpc 
boxSide=boxSide+10*boxSide/Ngrid; %increase boxsize so smoothed cells won't go out of bounds
clear m1 m2 ll

%% find location of sim cells in grid 
indX=ceil((coords(1,:)./boxSide+0.5).*Ngrid);
indY=ceil((coords(2,:)./boxSide+0.5).*Ngrid);
indZ=ceil((coords(3,:)./boxSide+0.5).*Ngrid);

cl=(gasStruct.Masses./gasStruct.Density).^(1/3); % "radius" of sim cells
gcl=ceil(cl./(boxSide/Ngrid));     %"radius" of sim cells in grid units

gcl=gcl+(1-mod(gcl,2)); % fix values to be only odd 

gind=(gcl-1)/2;   % extract no. of cells below and above the center cell;
%% set value and weights

value=gasStruct.Masses.*massUnit;

vv=value./(gcl.^3); % divide value by volume of grid cells - only for extensive values

% build grid 

cube=zeros(Ngrid,Ngrid,Ngrid);


tic
for i=1:cellNum
    
    
    cube(indX(i)-gind(i):indX(i)+gind(i),...
        indY(i)-gind(i):indY(i)+gind(i),...
        indZ(i)-gind(i):indZ(i)+gind(i))=cube(indX(i)-gind(i):indX(i)+gind(i),...
        indY(i)-gind(i):indY(i)+gind(i),...
        indZ(i)-gind(i):indZ(i)+gind(i))+vv(i);
end
   toc; 