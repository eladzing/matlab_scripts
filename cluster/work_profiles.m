global FILE_FORMAT;
FILE_FORMAT = '../data/%s_a1.001L%dMpc.dat';

SPHERE_RES = 50*2; %*16;
INTERP = 'linear';
MPSec = 8;
h = 0.7;

RRHOG = RHOG(MPSec);
RRHOTOT = RHOTOT(MPSec);
TT = T(MPSec);
SS = S(MPSec);

RHOG_Profile   = [];
RHOTOT_Profile = [];
T_Profile = [];
S_Profile = [];
R_Profile = 5:5:125;

[Sx Sy Sz] = uni_sphere(SPHERE_RES);
n = length(Sx(:));

for R = R_Profile
    disp(num2str(R))
    
    RSx = R*Sx+(128.5);
    RSy = R*Sy+(128.5);
    RSz = R*Sz+(128.5);

    %extensive properties (proportional to the size of the system)
    current_RHOG   = interp3(RRHOG,   RSx, RSy, RSz, INTERP);
    current_RHOTOT = interp3(RRHOTOT, RSx, RSy, RSz, INTERP);
    current_T      = interp3(TT,      RSx, RSy, RSz, INTERP);
    current_S      = interp3(SS,      RSx, RSy, RSz, INTERP);
    
    ds = 4*pi*((R*(MPSec*h)/256)^2)/n;
    
    RHOG_Profile(end+1)   = sum(sum(current_RHOG))*ds;
    RHOTOT_Profile(end+1) = sum(sum(current_RHOTOT))*ds;
    S_Profile(end+1)      = sum(sum(current_S))*ds;
    
    T_Profile(end+1)      = sum(sum(current_T.*current_RHOG))/sum(sum(current_RHOG));
end

% Fix R units
R_Profile = R_Profile*(MPSec*h)/256;

%we would like to integrate the Rhos to we multiply by dr
dr = R_Profile(2)-R_Profile(1);

%%% WARNING: We miss the inner sphere R<5, profile is shifted by constant
RHOG_Profile = cumsum(dr*RHOG_Profile);
RHOTOT_Profile = cumsum(dr*RHOTOT_Profile);


% Plots!
figure(1); semilogy(R_Profile, T_Profile, '*-'); title('T Profile'); xlabel('R (MPc)'); ylabel('log(T(R))'); saveas(gcf, getresultsdir(sprintf('T Profile (%d MPc).png',MPSec)));
figure(2); plot(R_Profile, S_Profile, '*-'); title('Entropy Profile'); xlabel('R (MPc)'); ylabel('Entropy'); saveas(gcf, getresultsdir(sprintf('S Profile (%d MPc).png',MPSec)));
figure(3); plot(R_Profile, RHOG_Profile, '*-'); title('Total Gas Mass Profile'); xlabel('R (MPc)'); ylabel('Total Gas Mass (Solar Masses)'); saveas(gcf, getresultsdir(sprintf('Cumulative Gas Mass Profile (%d MPc).png',MPSec)));
figure(4); plot(R_Profile, RHOTOT_Profile, '*-'); title('Total Mass Profile'); xlabel('R (MPc)'); ylabel('Total Mass (Solar Masses)'); saveas(gcf, getresultsdir(sprintf('Cumulative Total Mass Profile (%d MPc).png',MPSec)));
% TODO

figure(5); loglog(R_Profile, T_Profile, '*-'); title('T Profile'); xlabel('log(R) (MPc)'); ylabel('log(T(R))'); saveas(gcf, getresultsdir(sprintf('T Profile (%d MPc) loglog.png',MPSec)));
figure(6); loglog(R_Profile, S_Profile, '*-'); title('Entropy Profile'); xlabel('log(R) (MPc)'); ylabel('Entropy'); saveas(gcf, getresultsdir(sprintf('S Profile (%d MPc) loglog.png',MPSec)));
figure(7); loglog(R_Profile, RHOG_Profile, '*-'); title('Total Gas Mass Profile'); xlabel('log(R) (MPc)'); ylabel('Total Gas Mass (Solar Masses)'); saveas(gcf, getresultsdir(sprintf('Cumulative Gas Mass Profile (%d MPc) loglog.png',MPSec)));
figure(8); loglog(R_Profile, RHOTOT_Profile, '*-'); title('Total Mass Profile'); xlabel('log(R) (MPc)'); ylabel('Total Mass (Solar Masses)'); saveas(gcf, getresultsdir(sprintf('Cumulative Total Mass Profile (%d MPc) loglog.png',MPSec)));

