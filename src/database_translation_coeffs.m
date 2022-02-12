function [db_filename_ts] = database_translation_coeffs(N)
%%
% build the database for 2 kinds of translation coefficients (first and
% second types). Therefore, saving the time for repeating calculating these
% two types of translation coefficients.
% 
% NOTE: this database is created only for multi-particle system
% (multi_particle == 1).
%%

parameters_names;

if exist([db_filename_ts, '.mat']) ~= 0        % if already exist the database, 
    return;                                     % then do not create again for saving time.
end


%% particles' relative position

particles_Cartesian_data;

% "ij" means particle "i" toward particle "j"
r_ij = zeros(particle_number, particle_number);
theta_ij = zeros(particle_number, particle_number);
phi_ij = zeros(particle_number, particle_number);
% number "jj == 1" represents the probe particle "L"
for ii = 1:particle_number
    for jj = 1:particle_number
        if ii == jj
            continue;
        end
        [r_ij(ii,jj), theta_ij(ii,jj), phi_ij(ii,jj)] = ...
            coords_system_relative_positions(ii,jj);
    end
end
kr_ij = fluid_k * r_ij;


%% the fisrt and second type translation coefficients

Snmvu_1 = zeros(N+1, 2*N+1, N+1, 2*N+1, particle_number, particle_number);
Snmvu_2 = zeros(N+1, 2*N+1, N+1, 2*N+1, particle_number, particle_number);

% % ================== Single thread version ==================
% for nn = 0:N
%     for mm = -nn:nn
%         for nu = 0:N
%             for mu = -nu:nu
%                 for ii = 1:particle_number
%                     for jj = 1:particle_number
%                         if ii == jj
%                             Snmvu_1(nn+1, nn+mm+1, nu+1, nu+mu+1, ii, jj) = 0;
%                             Snmvu_2(nn+1, nn+mm+1, nu+1, nu+mu+1, ii, jj) = 0;
%                         else
%                             [Snmvu_1(nn+1, nn+mm+1, nu+1, nu+mu+1, ii, jj), Snmvu_2(nn+1, nn+mm+1, nu+1, nu+mu+1, ii, jj)] = ...
%                                 Snmvu_coeff(nn, mm, nu, mu, kr_ij(ii, jj), theta_ij(ii,jj), phi_ij(ii,jj));
%                         end
%                     end
%                 end
%             end
%         end
%     end
%     fprintf('Translation Coefficients Database Preparing %d%% \n', ...
%             round(100*nn/N));
% end
% % ================== Single thread version ==================

% % ============= Parallel version 1 for speed up =============
% for nn = 0:N
%     for mm = -nn:nn
%         for nu = 0:N
%             for mu = -nu:nu
%                 for ii = 1:particle_number
%                     indices_1 = nn+1;
%                     indices_2 = nn+mm+1;
%                     indices_3 = nu+1;
%                     indices_4 = nu+mu+1;
%                     parfor jj = 1:particle_number
%
%                         if ii == jj
%                             Snmvu_1(indices_1, indices_2, indices_3, indices_4, ii, jj) = 0;
%                             Snmvu_2(indices_1, indices_2, indices_3, indices_4, ii, jj) = 0;
%                         else
%                             [Snmvu_1(indices_1, indices_2, indices_3, indices_4, ii, jj), Snmvu_2(indices_1, indices_2, indices_3, indices_4, ii, jj)] = ...
%                                Snmvu_coeff(nn, mm, nu, mu, kr_ij(ii, jj), theta_ij(ii,jj), phi_ij(ii,jj));
%                         end                      
%
%                     end
%                 end
%             end
%         end
%     end
%     fprintf('Translation Coefficients Database Preparing %d%% \n', ...
%             round(100*nn/N));
% end
% % ============== Parallel version 1 for speed up ==============

% % ============== Parallel version 2 for speed up ==============
% for nn = 0:N
%     for mm = -nn:nn
%         for nu = 0:N
%             for mu = -nu:nu
%                 indices_1 = nn+1;
%                 indices_2 = nn+mm+1;
%                 indices_3 = nu+1;
%                 indices_4 = nu+mu+1;
%                 particle_number = particle_number;
%                 parfor ii = 1:particle_number
%
%                     temp_Snmvu_1 = zeros(1, particle_number);
%                     temp_Snmvu_2 = zeros(1, particle_number);
%                     for jj = 1:particle_number
%                         if ii == jj
%                             temp_Snmvu_1(jj) = 0;
%                             temp_Snmvu_2(jj) = 0;
%                         else
%                             [temp_Snmvu_1(jj), temp_Snmvu_2(jj)] = ...
%                                 Snmvu_coeff(nn, mm, nu, mu, kr_ij(ii, jj), theta_ij(ii,jj), phi_ij(ii,jj));
%                         end
%                     end 
%                     Snmvu_1(indices_1, indices_2, indices_3, indices_4, ii, :) = temp_Snmvu_1;
%                     Snmvu_2(indices_1, indices_2, indices_3, indices_4, ii, :) = temp_Snmvu_2;
%                     
%                 end
%             end
%         end
%     end
%     fprintf('Translation Coefficients Database Preparing %d%% \n', ...
%             round(100*nn/N));
% end
% % ============== Parallel version 2 for speed up ==============

% % ============= Parallel version 4 for speed up =============
% for nn = 0:N
%     for mm = -nn:nn
%         indices_1 = nn+1;
%         indices_2 = nn+mm+1;
%         particle_number = particle_number;
%         parfor nu = 0:N
%             
%             temp_Snmvu_1 = zeros(2*N+1, particle_number, particle_number);
%             temp_Snmvu_2 = zeros(2*N+1, particle_number,
%             particle_number);
%             for mu = -nu:nu
%                 indices_4 = nu+mu+1;
%                 
%                 for ii = 1:particle_number
%                     for jj = 1:particle_number
%                         if ii == jj
%                             temp_Snmvu_1(indices_4, ii, jj) = 0;
%                             temp_Snmvu_2(indices_4, ii, jj) = 0;
%                         else
%                             [temp_Snmvu_1(indices_4, ii, jj), temp_Snmvu_2(indices_4, ii, jj)] = ...
%                                 Snmvu_coeff(nn, mm, nu, mu, kr_ij(ii, jj), theta_ij(ii,jj), phi_ij(ii,jj));
%                         end
%                     end
%                 end
%             end
%             Snmvu_1(indices_1, indices_2, nu+1, :, :, :) = temp_Snmvu_1;
%             Snmvu_2(indices_1, indices_2, nu+1, :, :, :) = temp_Snmvu_2;
%         
%         end
%     end
%     fprintf('Translation Coefficients Database Preparing %d%% \n', ...
%             round(100*nn/N));
% end
% % ============= Parallel version 4 for speed up =============

% ============= Parallel version 5 for speed up =============
for nn = 0:N
    indices_1 = nn+1;
    particle_number = particle_number;
    parfor mm = -nn:nn
        
        temp_Snmvu_1 = zeros(N+1, 2*N+1, particle_number, particle_number);
        temp_Snmvu_2 = zeros(N+1, 2*N+1, particle_number, particle_number);   
        for nu = 0:N
            for mu = -nu:nu
                indices_3 = nu+1;
                indices_4 = nu+mu+1;
                for ii = 1:particle_number
                    for jj = 1:particle_number
                        if ii == jj
                            temp_Snmvu_1(indices_3, indices_4, ii, jj) = 0;
                            temp_Snmvu_2(indices_3, indices_4, ii, jj) = 0;
                        else
                            [temp_Snmvu_1(indices_3, indices_4, ii, jj), temp_Snmvu_2(indices_3, indices_4, ii, jj)] = ...
                                Snmvu_coeff(nn, mm, nu, mu, kr_ij(ii, jj), theta_ij(ii,jj), phi_ij(ii,jj));
                        end
                    end
                end
            end
        end
        Snmvu_1(indices_1, indices_1 + mm, :, :, :, :) = temp_Snmvu_1;
        Snmvu_2(indices_1, indices_1 + mm, :, :, :, :) = temp_Snmvu_2;
        
    end
    fprintf('Translation Coefficients Database Preparing %d%% \n', ...
            round(100*nn/N));
end
% ============= Parallel version 5 for speed up =============

% % ============= Parallel version 6 for speed up =============
% particle_number = particle_number;
% 
% parfor nn = 0:N
%     
%     temp_Snmvu_1 = zeros(2*N+1, N+1, 2*N+1, particle_number, particle_number);
%     temp_Snmvu_2 = zeros(2*N+1, N+1, 2*N+1, particle_number, particle_number);
%     for mm = -nn:nn
%         indices_2 = nn+mm+1;
%         for nu = 0:N
%             for mu = -nu:nu
%                 indices_3 = nu+1;
%                 indices_4 = nu+mu+1;
%                 for ii = 1:particle_number
%                     for jj = 1:particle_number
%                         if ii == jj
%                             temp_Snmvu_1(indices_2, indices_3, indices_4, ii, jj) = 0;
%                             temp_Snmvu_2(indices_2, indices_3, indices_4, ii, jj) = 0;
%                         else
%                             [temp_Snmvu_1(indices_2, indices_3, indices_4, ii, jj), temp_Snmvu_2(indices_2, indices_3, indices_4, ii, jj)] = ...
%                                 Snmvu_coeff(nn, mm, nu, mu, kr_ij(ii, jj), theta_ij(ii,jj), phi_ij(ii,jj));
%                         end
%                     end
%                 end
%             end
%         end
%     end
%     Snmvu_1(nn+1, :, :, :, :, :) = temp_Snmvu_1;
%     Snmvu_2(nn+1, :, :, :, :, :) = temp_Snmvu_2;
%     fprintf('Translation Coefficients Database Preparing %d%% \n', ...
%             round(100*nn/N));
% 
% end
% % ============= Parallel version 6 for speed up =============


%% save the database by meaningful name

save([db_filename_ts, '.mat'], 'Snmvu_1', 'Snmvu_2', 'r_ij');

%%