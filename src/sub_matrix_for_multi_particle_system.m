function B_ij = sub_matrix_for_multi_particle_system(ii, jj, Snmvu_2)
%%
% This function is used to calculate the sub-matrix of "Bij", which are
% size of "(N+1)^2 * (N+1)^2" (Here "N" as Truncation number).
%
% Truncation number N is loaded from "database_beam_scattering_coeffs.m".
% database maximum size (maximum Truncation number): db_size_nn = 20
% (seeing "parameters.m");
%%

parameters_names;

if multi_particle == 0
    disp('Single particle system do not need to run this function!');
    return;
end

[db_filename] = database_beam_scattering_coeffs();      
load([db_filename, '.mat'], 'N', 'db_bs_coeff', 'db_s_coeff');

B_ij = zeros((N+1)^2, (N+1)^2);

if ii == jj
    B_ij = -eye((N+1)^2);
    return;
else
%     [r_ji, theta_ji, phi_ji] = coords_system_relative_positions(jj, ii);
%     kr_ji = fluid_k * r_ji;
    
    sub_ii = 0;
    for nn = 0:N    
        for mm = -nn:nn
            sub_ii = sub_ii + 1;            % row
            sub_jj = 0;
            for nu = 0:N    
                for mu = -nu:nu
                    sub_jj = sub_jj + 1;    % colume
%                     B_ij(sub_ii, sub_jj) = db_s_coeff(nu+1, jj) * ...
%                         Snmvu_coeff(nu, mu, nn, mm, kr_ji, theta_ji, phi_ji, 2);
                    B_ij(sub_ii, sub_jj) = db_s_coeff(nu+1, jj) * ...
                        Snmvu_2(nu+1, nu+mu+1, nn+1, nn+mm+1, jj, ii);
                end
            end
        end
    end
end

%%