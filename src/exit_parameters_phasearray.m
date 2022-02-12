function exit_parameters_phasearray()
%%
% This function is used to reset the script "parameters.m".   
% This is a sub-function of GUI function of 'main_interface.m'
%%

%% Parameters need to computation in this script

file_parameters = 'parameters_reset.m';
file_phase_array = 'phase_array_beam_shape_coeff_reset.m';

%% Change specific parameters for the files

fid_parameters = fopen(file_parameters, 'r');
fid_phase_array = fopen(file_phase_array, 'r');

ii = 0;
content_parameters = {};
while feof(fid_parameters) ~= 1 
    ii = ii + 1;
    line = fgetl(fid_parameters);
    content_parameters{ii} = line;
end

ii = 0;
content_phase_array = {};
while feof(fid_phase_array) ~= 1 
    ii = ii + 1;
    line = fgetl(fid_phase_array);
    content_phase_array{ii} = line;
end

fclose(fid_parameters);
fclose(fid_phase_array);

%% Write back to the files

file_parameters = 'parameters.m';
file_phase_array = 'phase_array_beam_shape_coeff.m';

fid_parameters = fopen(file_parameters, 'w');
for ii = 1:length(content_parameters)
    fprintf(fid_parameters, '%s\n', content_parameters{ii});
end

fid_phase_array = fopen(file_phase_array, 'w');
for ii = 1:length(content_phase_array)
    fprintf(fid_phase_array, '%s\n', content_phase_array{ii});
end

fclose(fid_parameters);
fclose(fid_phase_array);

fclose all;
clear all;      % !!!!!!!!!

%%