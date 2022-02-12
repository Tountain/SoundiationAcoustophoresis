function [N, db_bs_coeff, db_s_coeff, db_Bessel, db_Hankel, db_Harmonics, ...
        db_b_s_coeff, db_rs_in_coeff, db_rs_out_coeff, particle_involve, layer_range] = ...
        loading_preparing_database(multi_particle, rr_visual, theta_visual, phi_visual)
%%
% this function is used to:
%   1.loading "beam-shape coefficient" and "scattering coefficient" for
%     single particle system; (if multi_particle == 0)
%   2.preparing the "interaction beam-shape coefficient of probe particle"
%     and "scattering coefficient of probe particle" for multi-particle
%     system. (if multi_particle ~= 0)
%%

% build the database, if already exist, then it will jump to next sentense.
[db_filename] = database_beam_scattering_coeffs();      
load([db_filename, '.mat'], 'N', 'db_bs_coeff', 'db_s_coeff');

% build in "pressure_contour.m" or "pressure_contour_vertical.m"
[db_filename_BHH] = database_Bessel_Hankel_Harmonics(N, rr_visual, theta_visual, phi_visual);
load([db_filename_BHH, '.mat'], 'db_Bessel', 'db_Hankel', 'db_Harmonics');

if multi_particle == 0
    db_b_s_coeff = 0;
    db_rs_in_coeff = 0;
    db_rs_out_coeff = 0;
    particle_involve = 0;
    layer_range = 0;
end
if multi_particle ~= 0      % for multi-particle system 
    [db_filename_ts] = database_translation_coeffs(N);
    load([db_filename_ts, '.mat'], 'Snmvu_1', 'Snmvu_2', 'r_ij');
    [db_filename_ibs] = database_interaction_beam_rescattering_coeffs(N, db_s_coeff, Snmvu_1, Snmvu_2, r_ij); 
    load([db_filename_ibs, '.mat'], 'db_b_s_coeff', 'db_rs_in_coeff', 'db_rs_out_coeff', 'particle_involve', 'layer_range');
end

%%