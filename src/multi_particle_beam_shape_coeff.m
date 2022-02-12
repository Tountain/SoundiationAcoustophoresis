function A = multi_particle_beam_shape_coeff(Snmvu_1)
%%
% This function is used to establish external beam-shape coefficient of all
% particles.
%
% Firstly, loading the external beam-shape coefficient for probe particle
% "L";
% Secondly, using translation relation among the probe coordinates and
% other source coordinates, we can determine all source particles'
% beam-shape coefficients.
%
% NOTE: before using the beam-shape coefficient for probe particle "anm,L", we
% need to rearrenge the indexes' rule, i.e. obtaining all numbers and
% rearrange them as "1D" matrix following the order: 
%       n = 0, m = 0;
%       n = 1, m = -1, 0, 1;
%       ...
%       n = N, m = -N, -N+1, ... , 0, ..., N-1, N;
% where "N" is Truncation number.
%%

parameters_names;

if multi_particle == 0
    disp('Single particle system do not need to run this function!');
    return;
end

[db_filename] = database_beam_scattering_coeffs();      
load([db_filename, '.mat'], 'N', 'db_bs_coeff', 'db_s_coeff');

%% rearrange the loading beam-shape coefficient "db_bs_coeff" as 1D matrix

AA_L = zeros((N+1)^2, 1);       % number (N+1)^2 is total elements number of "anm,L"
ii = 0;
for nn = 0:N
    for mm = -nn:nn
        ii = ii + 1;
        AA_L(ii) = db_bs_coeff(nn+1, nn+mm+1);
    end
end


%% using the translation relation for the beam-shape coefficient of other
%% source particles (seeing P.A.Martin@2006@Theorem 3.26 and J.H@2016@IEEE)

particles_Cartesian_data;
A = zeros((N+1)^2, particle_number);
A(:, 1) = AA_L;
for qq = 2:particle_number
%     [r_1q, theta_1q, phi_1q] = coords_system_relative_positions(1, qq);   % number "1" for the whole program represents the probe particle
%     kr_1q = fluid_k * r_1q;
    ii = 0;
    for nn = 0:N
        for mm = -nn:nn
            sum = 0;
            ii = ii + 1;
            for nu = 0:N
                for mu = -nu:nu
%                     sum = sum + Snmvu_coeff(nu, mu, nn, mm, kr_1q, theta_1q, phi_1q, 1) * ...
%                             db_bs_coeff(nu+1, nu+mu+1);
                    sum = sum + Snmvu_1(nu+1, nu+mu+1, nn+1, nn+mm+1, 1, qq) * ...
                            db_bs_coeff(nu+1, nu+mu+1);
                end
            end
            A(ii, qq) = sum;
        end
    end
    fprintf('All Particles'' Beam-Shape Coefficient Preparing %d%% \n', ...
            round(100*(qq-1)/(particle_number-1)));
end

A = reshape(A, (N+1)^2 * particle_number, 1);


%%