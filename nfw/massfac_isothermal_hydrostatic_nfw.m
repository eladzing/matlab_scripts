

mm=10.^(8:0.001:16);

cv=cvir_Mvir(mm,0);

aa=nfwA(1,cv);

xx=0.0001:0.0001:1;

for i=1:length(mm)
    bfac(i)=( trapz(xx,xx.^2.*exp(2/aa(i).*log(cv(i).*xx+1)./xx))).^-1;
end