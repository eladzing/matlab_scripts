load(sprintf('%s/virial%d', hpath, 8));
full_mgas = mgas_sphere(8);
glen= length(R_Profile)-1; %  ceil(256*RVIR/R_Profile(length(R_Profile)));

M_gas = [];
MGAS =[];
%M_gas = zeros(1,glen1);
for ridx = 1:glen
    if ridx < 2  
        mg = squeeze(full_mgas(ridx,:,:)).*(R_Profile(ridx));
    else
        mg = squeeze(full_mgas(ridx,:,:)).*(R_Profile(ridx)-R_Profile(ridx-1));
    end
    M_gas(end+1)=sum(mg(:));
            
end



MVIR
m_norm=560*(MVIR./1e13).^0.15.*(MGAS./1e13)  %virial accretion rate in Msun/yr