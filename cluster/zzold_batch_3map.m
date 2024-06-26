%list1=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
list1=[3 7];
%projs=[1 1 1; 1 1 1];
tag3='no';
result_dir='/home/carina/eladzing_a/cold_flows/tit3/printout/';
%result_dir='/home/eladzing/work/sshfs/titan3/cold_flows/tit3/printout/';

pflag='print';


%rmfac=0.0;rxfac=2;
%tmfac=-0.5;

zone=[0.02 0.2 1];
nthick=(1./256).*4; 

lims=zeros([5 4 2]);

lims(1,:,:)=[0.2 0.7;0.2 0.7;0.2 0.7;0.2 0.7]; % entropy
lims(2,:,:)=[0 1.1;0 1.1;0 1.1;0 1.1]; % entropy gradient
lims(3,:,:)=[-2 2;-2 2;-2 2;-2 2]; % Mach 
lims(4,:,:)=[-3 3;-3 3;-3 3;-3 3]; % flux
lims(5,:,:)=[-0.5 0.5;-0.5 0.5;-0.5 0.5;-0.5 0.5]; % radial entropy gradient
%lims(4,:,:)=[-10 10;-10 10;-10 10;-10 10]; % flux

for id=1:length(list1)
    %halopath=sprintf(hpath,list1(id));
    clustername=sprintf('CL%d',list1(id))
    new_env(clustername);
    plotproj=[1 1 1];
    %plotproj=projs(id,:);
    %stack{id,1}=clustername;
    switch list1(id)
        case {101,102,103,104,105,106,107,5}
            smallbox=2;bigbox=2;
        otherwise
            smallbox=1; bigbox=1;    
    end
    %bigbox=1;
    for ii=log2(smallbox):log2(bigbox)
        thick=nthick.*2.^ii;
        boxx=2^ii  
        s_lim=squeeze(lims(1,ii+1,:));
        gs_lim=squeeze(lims(2,ii+1,:));
        m_lim=squeeze(lims(3,ii+1,:));
        f_lim=squeeze(lims(4,ii+1,:)); 
        gsr_lim=squeeze(lims(5,ii+1,:)); 
        mappin
    end

    %close all
end


  
