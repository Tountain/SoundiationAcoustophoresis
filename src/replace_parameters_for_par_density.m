function replace_parameters_for_par_density(par_density)
%%
% This function is used to change the information 'fluid' in the
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
    if isempty(strfind(line, 'particle_rho = ')) ~= 1
        line_num_par_density = ii;
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
    if ii == line_num_par_density
        line = par_density;
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