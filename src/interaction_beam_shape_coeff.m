function db_ibs_coeff = interaction_beam_shape_coeff(BB, A, N)
%%
% this function is used to obtain the interaction beam-shape coefficient
% of probe particle "bnm,L" by "B = inv(BB) * A".
% 
% NOTE: in order to have a convinient interface with "pressure_contour.m",
% "pressure_contour_vertical.m", "radiation_force_based_NumIntPres.m",
% "radiation_force_based_Analyses.m", "radiation_torque_based_Analyses.m"
% and so on, we need to rearrange the solution of "B = inv(BB) * A" using
% database projecting rule (seeing "database_beam_scattering_coeffs.m" for
% details).
%%

particles_Cartesian_data;


%% solution of interaction beam-shape coefficient B

B = BB \ A;             % equavalue B = inv(BB) * A

if rcond(BB) < 10^(-10)
    fprintf('Warning: the condition number of the coefficient matrix BB is too high.\n')
end

upper = length(B) / particle_number;
% B_L = B(1 : upper);     % interaction beam-shape coefficient for probe particle


%% rearrange the "B_L" based on the database projecting rule

% db_ibs_coeff_L = zeros(N+1, 2*N+1);     % suffix "_L" means "probe particle"
db_ibs_coeff = zeros(N+1, 2*N+1, particle_number);
ii = 0;
for jj = 1:particle_number
    for nn = 0:N
        for mm = -nn:nn
            ii = ii + 1;
            % db_ibs_coeff_L(nn+1, nn+mm+1) = B_L(ii);
            db_ibs_coeff(nn+1, nn+mm+1, jj) = B(ii);
        end
    end
end

%%