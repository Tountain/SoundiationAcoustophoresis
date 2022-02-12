function absorption_factor()
%%
% This function is used to calculate the attenuation coefficient of mixture
% medium due to scattering attenuation. (temporally, ignore the viscosity)
%
% Ref. Wei_ronghao@1965@ActaPhyiscaSinica and Sewell@1910@RoyalSociety.
%%

%% define the host fluid parameters and particles' parameters

global absp_factor

% assume that no attenuation first and then based on "absp_factor == 0" to
% determine mixture's attenuation "absp_factor"    
absp_factor = 0;

parameters;

particles_num_per_unit_volume = particles_volume_fraction / (4/3 * pi * particle_radius^3);
particles_occupy_per_meter = particles_num_per_unit_volume^(1/3) * (2*particle_radius)

fluid_k_real = omega / fluid_c;
ratio_coeff = particles_num_per_unit_volume / fluid_k_real^2;


%% temporally build a database for "beam-shape and scattering coefficient"
%% for NO attenuation mixture

% [db_filename] = database_beam_scattering_coeffs();      
% load([db_filename, '.mat'], 'N', 'db_bs_coeff', 'db_s_coeff');
if multi_particle == 0
    [db_filename] = database_beam_scattering_coeffs();      
    load([db_filename, '.mat'], 'N', 'db_bs_coeff', 'db_s_coeff');
else
    [db_filename] = database_beam_scattering_coeffs();      
    load([db_filename, '.mat'], 'N', 'db_bs_coeff', 'db_s_coeff');
    [db_filename_ts] = database_translation_coeffs(N);
    load([db_filename_ts, '.mat'], 'Snmvu_1', 'Snmvu_2', 'r_ij');
    [db_filename_ibs] = database_interaction_beam_rescattering_coeffs(N, db_s_coeff, Snmvu_1, Snmvu_2, r_ij); 
    % load([db_filename_ibs, '.mat'], 'db_ibs_coeff');    % over_consider_scattering_diffusion(rescattering_has_been_over_considered_multiple_times) 
    % db_bs_coeff = db_ibs_coeff(:,:,1);      % "1" means probe particle
    % db_s_coeff = db_s_coeff(:,1);
    load([db_filename_ibs, '.mat'], 'db_b_s_coeff', 'db_rs_out_coeff', 'layer_range');
    for nn = 0:N
        for mm = -nn:nn
            db_bs_coeff(nn+1, nn+mm+1) = (db_b_s_coeff(nn+1, nn+mm+1, 1) + ...
                    db_rs_out_coeff(nn+1, nn+mm+1, length(layer_range))) ./ db_s_coeff(nn+1, 1);
        end
    end
    db_s_coeff = db_s_coeff(:,1);
end


%% calculate the attenuation factor "absp_factor"

ratio_energy_dissipation = 0;
ratio_energy_dissipation_vis = 0;
for nn = 0:N
    for mm = -nn:nn
    
        Attenuation_coeff = db_bs_coeff(nn+1, nn+mm+1) * db_s_coeff(nn+1, 1);
        ratio_energy_dissipation = ratio_energy_dissipation + (Attenuation_coeff * conj(Attenuation_coeff));
        
        ratio_energy_dissipation_vis = ratio_energy_dissipation_vis + ...
                        real((db_bs_coeff(nn+1, nn+mm+1) * conj(Attenuation_coeff)) + ...
                             (Attenuation_coeff * conj(Attenuation_coeff)));
        
    end
end
ratio_energy_dissipation = ratio_coeff * ratio_energy_dissipation;
ratio_energy_dissipation_vis = ratio_coeff * ratio_energy_dissipation_vis;

% NOTE: revised the absorption factor "absp_factor" for mixture fluid.
% As "absp_factor" define as global variable, so this change will affect
% the whole program after this assignment finished, PARTICULARLY in code
% "parameters.m" 
absp_factor = ratio_energy_dissipation / (2) / p_inlet^2;
absp_factor_vis = ratio_energy_dissipation_vis / (2) / p_inlet^2;


%% delete the temporally database
absp_factor
absp_factor_vis                 % not right!! even viscosity miu ~= 0.
absp_factor / fluid_k_real;

% delete([db_filename, '.mat']);
if multi_particle == 0
    delete([db_filename, '.mat']);
else
    delete([db_filename, '.mat']);
    delete([db_filename_ts, '.mat']);
    delete([db_filename_ibs, '.mat']);
end


%%
