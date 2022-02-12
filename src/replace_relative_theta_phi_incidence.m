function replace_relative_theta_phi_incidence(theta_inc)
%%
% This function is used to change the incidence information 'theta_inc' and
% 'phi_inc' in the script "parameters.m". 
% Considering that the 'theta_inc' is defined in [0,pi], while the input
% 'theta_inc' is varied between [0,2pi]. In order to fix out this
% controversy, if the 'theta_inc' = [pi.2pi], we have to set 
%       'theta_inc = 2*pi - theta_inc'; and 'phi_inc = pi'.
%%

%% fix the controversy

if theta_inc <= pi && theta_inc >= 0
    phi_inc = 0;
elseif theta_inc <= 2*pi && theta_inc > pi
    theta_inc = 2*pi - theta_inc;
    phi_inc = pi;
elseif theta_inc < 0 && theta_inc >= -pi
    theta_inc = - theta_inc;
    phi_inc = pi;
end

%% Parameters need to computation in this script

file_parameters = 'parameters.m';

%% Read in the function files

fid_parameters = fopen(file_parameters, 'r');

ii = 0;
while feof(fid_parameters) ~= 1 
    ii = ii + 1;
    line = fgetl(fid_parameters);
    if isempty(strfind(line, 'theta_inc = ')) ~= 1
        line_num_theta_inc = ii;
    end
    if isempty(strfind(line, 'phi_inc = ')) ~= 1
        line_num_phi_inc = ii;
    end
end

fclose(fid_parameters);

%% Change specific parameters for the files

fid_parameters = fopen(file_parameters, 'r');

ii = 0;
content_parameters = {};
while feof(fid_parameters) ~= 1 
    ii = ii + 1;
    line = fgetl(fid_parameters);
    if ii == line_num_theta_inc
        line = ['        theta_inc = ', num2str(theta_inc), ';                  % (theta_inc, phi_inc) = (pi/2, 0) for along X-axis'];
    end
    if ii == line_num_phi_inc
        line = ['        phi_inc = ', num2str(phi_inc), ';                    % (theta_inc, phi_inc) = (pi/2, pi/2) for along Y-axis'];
    end
    content_parameters{ii} = line;
end

fclose(fid_parameters);

%% Write back to the files

fid_parameters = fopen(file_parameters, 'w');

for ii = 1:length(content_parameters)
    fprintf(fid_parameters, '%s\n', content_parameters{ii});
end

fclose(fid_parameters);

fclose all;
clear all;      % !!!!!!!!!

%%