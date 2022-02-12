function reset_parameters_phasearray()
%%
% This function is used to reset the script "parameters.m".   
% This is a sub-function of GUI function of 'main_interface.m'
%%

%% Parameters need to computation in this script

file_parameters = 'parameters_reset.m';
file_phase_array = 'phase_array_beam_shape_coeff_reset.m';

%% Read in the function files (do not reset popup menus)

fid_parameters = fopen(file_parameters, 'r');

ii = 0;
while feof(fid_parameters) ~= 1 
    ii = ii + 1;
    line = fgetl(fid_parameters);
    if isempty(strfind(line, 'fluid = ''')) ~= 1
        line_num_fluid = ii;
    end
    if isempty(strfind(line, 'BC = ''')) ~= 1
        line_num_BC = ii;
    end
    if isempty(strfind(line, '    Cn = ')) ~= 1
        line_num_Cn = ii;
    end
    if isempty(strfind(line, 'wave_type = ''')) ~= 1 && isempty(strfind(line, '''plain'', ''zero-Bessel'' or ''non-zero-Bessel'' and others')) ~= 1
        line_num_wavetype = ii;
    end
    if isempty(strfind(line, 'particle_radius = ')) ~= 1
        line_num_par_radius = ii;
    end
end

fclose(fid_parameters);


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


fid_parameters = fopen(file_parameters, 'r');       % keep the pop up manu features
ii = 0;
content_parameters2 = {};
while feof(fid_parameters) ~= 1 
    ii = ii + 1;
    line = fgetl(fid_parameters);
    content_parameters2{ii} = line;
end
fclose(fid_parameters);


fid_parameters = fopen(file_parameters, 'w');
for ii = 1:length(content_parameters)
    if ii ~= line_num_fluid && ii ~= line_num_BC && ii ~= line_num_Cn && ii ~= line_num_wavetype && ii ~= line_num_par_radius
        fprintf(fid_parameters, '%s\n', content_parameters{ii});
    else
        fprintf(fid_parameters, '%s\n', content_parameters2{ii});   % keep the pop up manu features
    end
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