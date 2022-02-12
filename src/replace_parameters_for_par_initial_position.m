function replace_parameters_for_par_initial_position(positionX, positionY, positionZ)
%%
% This function is used to change the information 'deviationX, Y, Z' in the
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
    if isempty(strfind(line, 'deviationX = ')) ~= 1
        line_num_deviationX = ii;
    end
    if isempty(strfind(line, 'deviationY = ')) ~= 1
        line_num_deviationY = ii;
    end
    if isempty(strfind(line, 'deviationZ = ')) ~= 1
        line_num_deviationZ = ii;
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
        line = positionX;
    end
    if ii == line_num_deviationY
        line = positionY;
    end
    if ii == line_num_deviationZ
        line = positionZ;
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