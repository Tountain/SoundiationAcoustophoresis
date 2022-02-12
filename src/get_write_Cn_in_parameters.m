function get_write_Cn_in_parameters(handles)
%%
% This function is used to get the 'Cn' from Edit text in GUI and write it
% into "parameters.m". 
% For user-specified geometry.
%%

str_Cn = get(handles.mapping_coeff, 'String');
str_Cn = ['    Cn = ', str_Cn, ';'];

%% replace Cn

replace_parameters_for_Cn(str_Cn);

%%
