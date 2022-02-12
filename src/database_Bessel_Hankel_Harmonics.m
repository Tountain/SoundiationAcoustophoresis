function [db_filename_BHH] = database_Bessel_Hankel_Harmonics(N, rr_visual, theta_visual, phi_visual)
%%
% this function is used to build the database of all repeat value of
% Spherical Bessel, Hankel and Harmonics functions used in
% "partial_wave_expansion.m".
%
% As the function "partial_wave_expansion.m" is the core function to obtain
% pressure contour. Therefore, the speed of program mainly depends on the
% running speed of "partial_wave_expansion.m" so depends on
% Bessel, Hankel, Harmonics functions. So, in order to avoid repeating
% calling "sBessel.m", "sHankel.m" and "sHarmonics.m", we build a database
% here.
%%

parameters_names;

if exist([db_filename_BHH, '.mat']) ~= 0        % if already exist the database, 
    return;                                     % then do not create again for saving time.
end

kr_visual = fluid_k * rr_visual;
kr_len = length(kr_visual);
theta_len = length(theta_visual);
phi_len = length(phi_visual);

db_Bessel = zeros(N+1, kr_len);
db_Hankel = zeros(N+1, kr_len);
db_Harmonics = zeros(N+1, 2*N+1, theta_len, phi_len);


%% build database
% projecting relations are similar with "database_beam_scattering_coeffs.m"

for nn = 0:N
    for kr_ii = 1:kr_len
        db_Bessel(nn+1,kr_ii) = sBessel(nn,kr_visual(kr_ii),1);
        db_Hankel(nn+1,kr_ii) = sHankel(nn,kr_visual(kr_ii),1);
    end
end

for nn = 0:N
    indices_1 = nn + 1;
    % parfor mm = -nn:nn   % Matlab2016b can work, but Matlab2010a can not.
    for mm = -nn:nn
        for theta_jj = 1:theta_len
            for phi_kk = 1:phi_len
                db_Harmonics(indices_1,indices_1+mm,theta_jj,phi_kk) = ...
                    sHarmonics(nn,mm,theta_visual(theta_jj),phi_visual(phi_kk));
            end
        end
    end
    fprintf('Bessel, Hankel and Harmonics Functions Database Preparing %d%% \n',round(100*nn/N));
end


%% save these databases by meaningful name

save([db_filename_BHH, '.mat'], 'db_Bessel', 'db_Hankel', 'db_Harmonics');
