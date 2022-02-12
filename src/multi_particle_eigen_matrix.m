function BB = multi_particle_eigen_matrix(multi_particle, Snmvu_2)
%%
% this function is used to establish eigen-matrix for multi-particle system
% by arranging all sub-matrix "B_ij" obtaining from sub-function
% "sub_matrix_for_multi_particle_system.m".
%%

if multi_particle == 0
    disp('Single particle system do not need to run this function!');
    return;
end

particles_Cartesian_data;

BB = [];
for ii = 1:particle_number
    BB_row = [];
    for jj = 1:particle_number
        B_ij = sub_matrix_for_multi_particle_system(ii, jj, Snmvu_2);
        BB_row = [BB_row, -B_ij];
        fprintf('Eigen-matrix Preparing %d%% \n', ...
            round(100*((ii-1)/particle_number + 1/particle_number * jj/particle_number)));
    end
    BB = [BB; BB_row];
end
        
%%