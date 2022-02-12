function replace_parameters_for_par_initial_orientation(theta_X, theta_Y, theta_Z)
%%
% This function is used to change the information 'theta_X, _Y, _Z' in the
% script "parameters.m".   
% This is a sub-function of GUI function of 'main_interface.m'
%%

%% Parameters need to computation in this script

file_parameters = 'parameters.m';

%% Read in the function files

fid_parameters = fopen(file_parameters, 'r');

ii = 0;
while feof(fid_parameters) ~= 1 
    ii = ii + 1;
    line = fgetl(fid_parameters);
    if isempty(strfind(line, 'theta_x = ')) ~= 1
        line_num_thetaX = ii;
    end
    if isempty(strfind(line, 'theta_y = ')) ~= 1
        line_num_thetaY = ii;
    end
    if isempty(strfind(line, 'theta_z = ')) ~= 1
        line_num_thetaZ = ii;
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
    if ii == line_num_thetaX
        line = theta_X;
    end
    if ii == line_num_thetaY
        line = theta_Y;
    end
    if ii == line_num_thetaZ
        line = theta_Z;
    end
    content_parameters{ii} = line;
end

fclose(fid_parameters);

%% Write back to the files

% str_Cn = Cn(find(Cn == '[') : find(Cn == ']'));

fid_parameters = fopen(file_parameters, 'w');
 
for ii = 1:length(content_parameters)
    fprintf(fid_parameters, '%s\n', content_parameters{ii});
end

fclose(fid_parameters);
 
fclose all;
clear all;      % !!!!!!!!!

%%