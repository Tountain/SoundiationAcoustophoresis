function [forces_matrix_matlab, torques_matrix_matlab] = readin_the_force_torque()
%%
% this function is used to read in the saved file for forces and torques in
% different incidence angles 'theta_inc =
% [0,pi/12,2*pi/12,3*pi/12,4*pi/12,5*pi/12,6*pi/12]' and length of mapping
% coefficients 'length(Cn)=[2,3,4,5]' (for sphere, ellipsoid, cone and
% diamond).
%%

parameters;

theta_inc = [0, pi/12, 2*pi/12, 3*pi/12, 4*pi/12, 5*pi/12, 6*pi/12];
length_Cn = [2, 3, 4, 5];
% freq = [1.00 1.50 2.00 2.50 3.00 3.50 4.00 4.50 5.00 5.50 6.00 6.50, 7.00 7.50 8.00]*10^6;
freq = [1.50 2.00 2.50 3.00 3.50 4.00 4.50 5.00 5.50 6.00 6.50, 7.00 7.50 8.00]*10^6;

forces_matrix_matlab = [];
torques_matrix_matlab = [];

for ii = 1 : length(theta_inc)              % outer loop: theta_inc
    for jj = 1 : length(length_Cn)          % middle loop: length_Cn
        for kk = 1 : length(freq)           % inner loop: freq
            f_ka = 2*pi*freq(kk)*particle_radius / fluid_c;
            loading_file_name = ['Arb', num2str(roundn(theta_inc(ii)*180/pi, -1)), ...
                '_0_irregular_1_Cn_', num2str(round(length_Cn(jj))), ...
                '_Pin_1Pa_ka', num2str(roundn(f_ka, -2)), ...
                '_plain_a50um_freq', num2str(round(10^-3*freq(kk))) , ...
                'kHz_X1(pr)0_Y1(pr)0_Z1(pr)0_olive_oil_range20_GridRes100A.mat'];
            dir_file_forces = ['.\data_forces_water\', loading_file_name];
            dir_file_torques = ['.\data_torques_water\', loading_file_name];
            load(dir_file_forces);  % for Frad_x, Frad_y, Frad_z
            load(dir_file_torques); % for Torque_x, Torque_y, Torque_z
            forces_matrix_matlab = [forces_matrix_matlab; Frad_x, Frad_y, Frad_z];
            torques_matrix_matlab = [torques_matrix_matlab; Torque_x, Torque_y, Torque_z];
        end
    end
    fprintf('Read in finished for incidence angle theta_inc = %d.\n', roundn(theta_inc(ii)*180/pi, -1));
end

% saving the 'forces_matrix' and 'torques_matrix' into file 'data.xlsx' in
% sheet 'original data'

matrix_matlab = [forces_matrix_matlab, torques_matrix_matlab];
xlswrite('data_matlab.xlsx', matrix_matlab, 'Water Original data');


%%