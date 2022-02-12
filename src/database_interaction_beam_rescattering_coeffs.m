function [db_filename_ibs] = ...
    database_interaction_beam_rescattering_coeffs(N, db_s_coeff, Snmvu_1, Snmvu_2, r_ij)
%%
% build the database for "interaction beam-shape coeffcient of all
% particle" and the "rescattering coefficient" "rs_nm" for calculate the
% far field pressure distribution of multi-particle system to avoid
% repeating solving the eigen-matrix in "interaction_beam_shape_coeff.m".
%
% NOTE: this database is created only for multi-particle system
% (multi_particle == 1).
%%

parameters_names;

if exist([db_filename_ibs, '.mat']) ~= 0        % if already exist the database, 
    return;                                     % then do not create again for saving time.
end

%% solving eigen-matrix for interaction beam-shape coeffcient of all
%% particle

BB = multi_particle_eigen_matrix(multi_particle, Snmvu_2);         
A = multi_particle_beam_shape_coeff(Snmvu_1);
db_ibs_coeff = interaction_beam_shape_coeff(BB, A, N);
% NOTE: "db_ibs_coeff" only uses for multi-particle force&torque.
%       i.e. "multi_particle_force_torque ~= 0"

% if we just want to obtain the radition force or torque of multi-particle
% system (then "multi_particle_force_torque ~= 0" is TURE), we do not care
% the farfield (outer field) pressure or velocity potential. Therefore, we
% do not need to carry on left codes.

% if multi_particle_force_torque ~= 0
%     save([db_filename_ibs, '.mat'], 'db_ibs_coeff');
%     return;
% end


%% through "snm,j = bnm,j * sn,j", getting "snm,j" and save in "db_b_s_coeff"

particles_Cartesian_data;

% scattering coefficient for multi-particle system
db_b_s_coeff = zeros(size(db_ibs_coeff));   
for ii = 1:particle_number
    for nn = 0:N
        for mm = -nn:nn
            db_b_s_coeff(nn+1, nn+mm+1, ii) = db_ibs_coeff(nn+1, nn+mm+1, ii) * db_s_coeff(nn+1, ii);
        end
    end
end


%% average distance of particles away probe particle
%% outer rescattering coefficient

r_iL = r_ij(:,1);
% "r_layer" means the distance from the closest to farest of all particles
% "particle_order" means the order of which particle closer to probe particle
[r_order, particle_order] = sort(r_iL);

% "particle_involve(layer)" means inner particles' number in terms of "layer"
% "layer_range(layer)" means width of "layer" of "r" between "layer_range(layer-1)" to "layer_range(layer)" 
[particle_involve, layer_range] = relation_of_layer_and_involve_particles(r_order);

db_rs_in_coeff = zeros(N+1,2*N+1,length(layer_range));
db_rs_out_coeff = zeros(N+1,2*N+1,length(layer_range));
for layer = 1:length(layer_range)
    for nn = 0:N
        for mm = -nn:nn
            [db_rs_in_coeff(nn+1,nn+mm+1,layer), db_rs_out_coeff(nn+1,nn+mm+1,layer)] = ...
                rescattering_coefficient(nn, mm, particle_involve(layer), particle_order, db_b_s_coeff, Snmvu_1, Snmvu_2, N);
        end
    end
    fprintf('Rescattering Coefficients Database Preparing %d%% \n', ...
            round(100*layer/length(layer_range)));
end

% ======================== OLD VERSION ======================    
% r_kL = zeros(particle_number, 1);
% theta_kL = zeros(particle_number, 1);
% phi_kL = zeros(particle_number, 1);
% for kk = 2:particle_number
%     [r_kL(kk), theta_kL(kk), phi_kL(kk)] = coords_system_relative_positions(kk, 1);
% end
% r_avg = sum(r_kL)/(particle_number-1);      % except of "L" itself
% kr_kL = fluid_k * r_kL;
% 
% for nn = 0:N
%     for mm = -nn:nn
%         rs_nm_L = 0;
%         for kk = 2:particle_number
%             rs_nm_L = rs_nm_L + p_rs_nmj(db_ibs_coeff, db_s_coeff, kr_kL, theta_kL, phi_kL, N, kk, nn, mm);
%         end
%         db_rs_coeff(nn+1, nn+mm+1) = rs_nm_L;
%     end
%     fprintf('Rescattering Coefficients Database Preparing %d%% \n', ...
%             round(100*nn/N));
% end
% ======================== OLD VERSION ======================    


%% save the database by meaningful name

% save([db_filename_ibs, '.mat'], 'db_ibs_coeff', 'db_rs_coeff', 'r_kL');
save([db_filename_ibs, '.mat'], 'db_ibs_coeff', 'db_b_s_coeff', 'db_rs_in_coeff', 'db_rs_out_coeff', 'particle_involve', 'layer_range');

%%