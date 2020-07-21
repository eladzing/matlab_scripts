function ro_nfw=nfw_profile(r,ro_s,rs)

ro_nfw=ro_s./((r./rs).*(1+r./rs).^2);