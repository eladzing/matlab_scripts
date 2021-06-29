%% explolre beta & fg dependence of stripping

%fg=0.1:0.1:0.3;%    0:0.01:0.25;
%beta=0.1:0.1:0.5; %       0.1:0.05:1;
fg=0.01:0.0025:0.25;
beta=0.2:0.01:1;
BT=0;%.2;

md=0.2;
mass=10.^10;
lambda=0.035;
len=length(fg)*length(beta);

ffg=reshape(repmat(fg,[length(beta) 1]),[len 1]);
bbeta=repmat(beta,[1 length(fg)])';
mmass=mass.*ones(size(ffg));
BT=BT.*ones(size(ffg));
ffb=BT./(1-BT);

mmass=mmass./(1+ffg+ffb);  % set total mass as constant

if Flag
    tic
    res=rscale_mmw_array(mmass,'lambda',lambda,'md',md,...
        'fg',ffg,'beta',bbeta,'BT',BT(i));
    toc
end
rd1=reshape(res.rd,[length(beta) length(fg)]);
si=mmass./(2*pi*res.rd.^2);
sig1=reshape(si,[length(beta) length(fg)]);

hf1=figure;
imagesc(fg,beta,log10(sig1))
set(gca,'Ydir','normal')
colorbar
xlabelmine('$f_g$')
ylabelmine('$\beta$')
titlemine('$\Sigma_s$')

%% calculate rps factor
pval=rps_factor_expdisk('etap',1,'sigma',mmass./(2*pi*res.rd.^2),'fd',ffg,...
    'Mc',1e15);
ppval=reshape(pval,[length(beta) length(fg)]);

hf2=figure;
imagesc(fg,beta,log10(ppval))
set(gca,'Ydir','normal')
colorbar
xlabelmine('$f_g$')
ylabelmine('$\beta$')
titlemine('$P_{rps}$')


r=0.01:0.01:100;
mStr=zeros(size(ffg));
for i=1:length(ffg)
    
    f=disk_force_reduced(r,'fg',ffg(i),'beta',bbeta(i),'fb',ffb(i));
    
    if max(f)>pval(i)
        rStr=interp1(f,r,pval(i));
        mStr(i)=exp_disk_mass(rStr,bbeta(i));
    else
        mStr(i)=0;
    end
    
end

mStr2d=reshape(mStr,[length(beta) length(fg)]);

hf3=figure;
imagesc(fg,beta,1-mStr2d)
set(gca,'Ydir','normal')
colorbar
xlabelmine('$f_g$')
ylabelmine('$\beta$')
titlemine('Mass stripped')
