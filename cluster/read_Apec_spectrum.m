%% read spectrum file for a single cell of the Apec code. 

% open file 
function res=read_Apec_spectrum
path='C:\Users\eladzing\Documents\cluster\xray_code\Apec\oneCell';
filterSpec=[path,'\*.dat'];

file=uigetfile(filterSpec);

fid=fopen([path,'/',file]);

% read temperature and abundance
res.temp=fscanf(fid,'%e',1); 
res.abun=fscanf(fid,'%e',1);
res.zred_emit=fscanf(fid,'%e',1);

% read
[spec,count]=fscanf(fid,'%e');
res.spec=reshape(spec,[2 count/2])';