function replace_position_theta(position, rotation)
%%
% This function is used to change the incidence information 'theta_inc' and
% 'phi_inc' in the script "parameters.m". 
% Considering that the 'theta_inc' is defined in [0,pi], while the input
% 'theta_inc' is varied between [0,2pi]. In order to fix out this
% controversy, if the 'theta_inc' = [pi.2pi], we have to set 
%       'theta_inc = 2*pi - theta_inc'; and 'phi_inc = pi'.
%%

%% Parameters need to computation in this script

file_parameters = 'parameters.m';

%% Read in the function files

fid_parameters = fopen(file_parameters, 'r');

ii = 0;
while feof(fid_parameters) ~= 1 
    ii = ii + 1;
    line = fgetl(fid_parameters);
    if isempty(strfind(line, 'deviationX = ')) ~= 1
        line_num_deviationX = ii;
    end
    if isempty(strfind(line, 'deviationY = ')) ~= 1
        line_num_deviationY = ii;
    end
    if isempty(strfind(line, 'deviationZ = ')) ~= 1
        line_num_deviationZ = ii;
    end
    if isempty(strfind(line, 'theta_x = ')) ~= 1
        line_num_theta_x = ii;
    end
    if isempty(strfind(line, 'theta_y = ')) ~= 1
        line_num_theta_y = ii;
    end
    if isempty(strfind(line, 'theta_z = ')) ~= 1
        line_num_theta_z = ii;
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
    if ii == line_num_deviationX
        line = ['deviationX = ', num2str(position(1)), ';'];
    end
    if ii == line_num_deviationY
        line = ['deviationY = ', num2str(position(2)), ';'];
    end
    if ii == line_num_deviationZ
        line = ['deviationZ = ', num2str(position(3)), ';'];
    end
    if ii == line_num_theta_x
        line = ['    theta_x = ', num2str(rotation(1)), ';                % particle counter-clock wise rotation (position Tx) along x-axis for positive ''theta_x'''];
    end
    if ii == line_num_theta_y
        line = ['    theta_y = ', num2str(rotation(2)), ';                % particle counter-clock wise rotation (position Ty) along y-axis for positive ''theta_y'''];
    end
    if ii == line_num_theta_z
        line = ['    theta_z = ', num2str(rotation(3)), ';                % particle counter-clock wise rotation (position Tz) along z-axis for positive ''theta_z'''];
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