function varargout = main_interface(varargin)
% MAIN_INTERFACE M-file for main_interface.fig
%      MAIN_INTERFACE, by itself, creates a new MAIN_INTERFACE or raises the existing
%      singleton*.
%
%      H = MAIN_INTERFACE returns the handle to a new MAIN_INTERFACE or the handle to
%      the existing singleton*.
%
%      MAIN_INTERFACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN_INTERFACE.M with the given input arguments.
%
%      MAIN_INTERFACE('Property','Value',...) creates a new MAIN_INTERFACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main_interface_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_interface_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main_interface

% Last Modified by GUIDE v2.5 02-Jan-2022 14:20:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_interface_OpeningFcn, ...
                   'gui_OutputFcn',  @main_interface_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before main_interface is made visible.
function main_interface_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main_interface (see VARARGIN)

% Choose default command line output for main_interface
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes main_interface wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = main_interface_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in reset_data.
function reset_data_Callback(hObject, eventdata, handles)
% hObject    handle to reset_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% reset the popup menus??
% handles.irregular.value = 1;

% set particle parameters
set(handles.par_density, 'String', '15');
% set(handles.par_radius, 'String', '0.002');
set(handles.initial_position, 'String', '[0 0 0]');
set(handles.initial_orientation, 'String', '[0 0 0]');

% set medium parameters
% set(handles.fluid_rho, 'String', '1.224');
% set(handles.fluid_c, 'String', '340');
% set(handles.fluid_vis, 'String', '0.0000185');

% set wavefront parameters
set(handles.freq, 'String', '40000');
set(handles.integration_radius, 'String', '0.005');

set(handles.p0, 'String', '1');

set(handles.trans_num, 'String', '5');
set(handles.trans_radius, 'String', '0.005');
set(handles.v0, 'String', '1.5');
set(handles.dt, 'String', '0.02');
set(handles.trans_position, 'String', '[0 0 0; 0.01 0 0; -0.01 0 0; 0 0.01 0; 0 -0.01 0]');
set(handles.phase_delay, 'String', '[0 0 0 0 0]');
set(handles.amp_delay, 'String', '[1 1 1 1 1]');

% set radiation force and torque
set(handles.rad_force, 'String', '');
set(handles.rad_torque, 'String', '');

% set acoustophoresis parameters
set(handles.delta_t, 'String', '0.0001');
set(handles.t_end, 'String', '0.1');

% set all static texts as blank
set(handles.wavefield_status, 'String', '');
set(handles.force_torque_status, 'String', '');
set(handles.current_time, 'String', '');
set(handles.acoustophoresis_status, 'String', '');

reset_parameters_phasearray();
clc;clear;


% --- Executes on button press in reset_exit.
function reset_exit_Callback(hObject, eventdata, handles)
% hObject    handle to reset_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

exit_parameters_phasearray();
delete(handles.figure1);
clc;clear;


% --- Executes on selection change in boundary_conditions.
function boundary_conditions_Callback(hObject, eventdata, handles)
% hObject    handle to boundary_conditions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns boundary_conditions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from boundary_conditions

str = get(hObject, 'String');
val = get(hObject, 'Value');

switch str{val}
    case 'Sound-hard B.C.'
        % change the 'irregular' variable in 'parameters'
        BC = 'BC = ''rigid'';               % Boundary Conditions "rigid" or "soft".';
        replace_parameters_for_BC(BC);
	case 'Sound-soft B.C.'
        % readin 'parameters' for 'Cn'
        BC = 'BC = ''soft'';               % Boundary Conditions "rigid" or "soft".';
        replace_parameters_for_BC(BC);
    otherwise
        BC = 'BC = ''rigid'';               % Boundary Conditions "rigid" or "soft".';
        replace_parameters_for_BC(BC);
end



% --- Executes during object creation, after setting all properties.
function boundary_conditions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to boundary_conditions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in medium_type.
function medium_type_Callback(hObject, eventdata, handles)
% hObject    handle to medium_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns medium_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from medium_type

str = get(hObject, 'String');
val = get(hObject, 'Value');

switch str{val}
    case 'Air'
        % change the 'fluid' variable in 'parameters'
        medium = 'fluid = ''air'';               % medium types "air" or "water" and others';
        replace_parameters_for_fluid(medium);
        set(handles.fluid_rho, 'String', '1.224');
        set(handles.fluid_c, 'String', '340');
        set(handles.fluid_vis, 'String', '0.0000185');
	case 'Water'
        % change the 'fluid' variable in 'parameters'
        medium = 'fluid = ''water'';               % medium types "air" or "water" and others';
        replace_parameters_for_fluid(medium);
        set(handles.fluid_rho, 'String', '1000');
        set(handles.fluid_c, 'String', '1500');
        set(handles.fluid_vis, 'String', '0.00101');
    otherwise
        medium = 'fluid = ''air'';               % medium types "air" or "water" and others';
        replace_parameters_for_fluid(medium);
        set(handles.fluid_rho, 'String', '1.224');
        set(handles.fluid_c, 'String', '340');
        set(handles.fluid_vis, 'String', '0.0000185');
end



% --- Executes during object creation, after setting all properties.
function medium_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to medium_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fluid_rho_Callback(hObject, eventdata, handles)
% hObject    handle to fluid_rho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fluid_rho as text
%        str2double(get(hObject,'String')) returns contents of fluid_rho as a double


% --- Executes during object creation, after setting all properties.
function fluid_rho_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fluid_rho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fluid_c_Callback(hObject, eventdata, handles)
% hObject    handle to fluid_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fluid_c as text
%        str2double(get(hObject,'String')) returns contents of fluid_c as a double


% --- Executes during object creation, after setting all properties.
function fluid_c_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fluid_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fluid_vis_Callback(hObject, eventdata, handles)
% hObject    handle to fluid_vis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fluid_vis as text
%        str2double(get(hObject,'String')) returns contents of fluid_vis as a double


% --- Executes during object creation, after setting all properties.
function fluid_vis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fluid_vis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function freq_Callback(hObject, eventdata, handles)
% hObject    handle to freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of freq as text
%        str2double(get(hObject,'String')) returns contents of freq as a double


% --- Executes during object creation, after setting all properties.
function freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in wavefront_selection.
function wavefront_selection_Callback(hObject, eventdata, handles)
% hObject    handle to wavefront_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns wavefront_selection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from wavefront_selection

str = get(hObject, 'String');
val = get(hObject, 'Value');

switch str{val}
    case 'Plane travelling wave'
        % change the 'irregular' variable in 'parameters'
        wave_type = 'wave_type = ''plain'';                   % ''plain'', ''zero-Bessel'' or ''non-zero-Bessel'' and others';
        replace_parameters_for_wavetype(wave_type);
	case 'Transducer array (circular oscillator)'
        % readin 'parameters' for 'Cn'
        wave_type = 'wave_type = ''phase_array_transducer'';                   % ''plain'', ''zero-Bessel'' or ''non-zero-Bessel'' and others';
        replace_parameters_for_wavetype(wave_type);
    otherwise
        wave_type = 'wave_type = ''plain'';                   % ''plain'', ''zero-Bessel'' or ''non-zero-Bessel'' and others';
        replace_parameters_for_wavetype(wave_type);
end


% --- Executes during object creation, after setting all properties.
function wavefront_selection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wavefront_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in wavefield_planetype.
function wavefield_planetype_Callback(hObject, eventdata, handles)
% hObject    handle to wavefield_planetype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.wavefield_status, 'String', 'Busy');
pause(0.001);

% replace particle information in "parameters.m"
str_particle_radius = get(handles.par_radius, 'String');
par_radius = ['    particle_radius = ',str_particle_radius,';                   % particle radial mean-radius (a); 30 um;'];
replace_parameters_for_par_radius(par_radius);

str_par_initial_position = get(handles.initial_position, 'String');
par_initial_position = str2num(str_par_initial_position);
positionX = ['deviationX = ',num2str(par_initial_position(1)),';'];
positionY = ['deviationY = ',num2str(par_initial_position(2)),';'];
positionZ = ['deviationZ = ',num2str(par_initial_position(3)),';'];
replace_parameters_for_par_initial_position(positionX, positionY, positionZ);

str_par_initial_orientation = get(handles.initial_orientation, 'String');
par_initial_orientation = str2num(str_par_initial_orientation);
theta_X = ['    theta_x = ',num2str(par_initial_orientation(1)),';                % particle counter-clock wise rotation (position Tx) along x-axis for positive ''theta_x'''];
theta_Y = ['    theta_y = ',num2str(par_initial_orientation(2)),';                % particle counter-clock wise rotation (position Ty) along y-axis for positive ''theta_y'''];
theta_Z = ['    theta_z = ',num2str(par_initial_orientation(3)),';                % particle counter-clock wise rotation (position Tz) along z-axis for positive ''theta_z'''];
replace_parameters_for_par_initial_orientation(theta_X, theta_Y, theta_Z);

% replace the frequency in "parameters.m"
str_freq = get(handles.freq, 'String');
freq = ['freq = ',str_freq,';               % 1 MHz'];
replace_parameters_for_frequency(freq);

str_integration_radius = get(handles.integration_radius, 'String');
integration_radius = ['b = ',str_integration_radius,';           % integral radial distance (b) for beam-shape coeffcient, EMPIRICAL: MUST follow  b > a and b < inter_dist BUT not b >> a'];
replace_parameters_for_integration_radius(integration_radius);

% get the current wavefront type and write the "parameters.m" and/or
% "phase_array_beam_shape_coeff.m"
wavefront_str = get(handles.wavefront_selection, 'String');
wavefront_val = get(handles.wavefront_selection, 'Value');
switch wavefront_str{wavefront_val}
    case 'Plane travelling wave'
        % replace wavefront parameters 'p_inlet' into "parameters.m" for plane wave
        str_p0 = get(handles.p0, 'String');
        p_inlet = ['p_inlet = ',str_p0,';               % initial amplitude of proba particle [Pa] if no attenuation effects'];
        replace_parameters_for_p_inlet(p_inlet);
    case 'Transducer array (circular oscillator)'
        % replace wavefront parameters into "parameters.m" and
        % "phase_array_beam_shape_coeff.m" for transducer array
        str_trans_num = get(handles.trans_num, 'String');
        trans_num = ['transducer_number = ',str_trans_num,';          % NOTE: keep same in "total_effects_all_plain_wave_component.m". for phase array'];
        replace_parameters_for_trans_num(trans_num);

        str_trans_radius = get(handles.trans_radius, 'String');
        trans_radius = ['trans_radius = ',str_trans_radius,';           % NOTE: keep same in "total_effects_all_plain_wave_component.m". radius of transducer is 5mm'];
        replace_parameters_for_trans_radius(trans_radius);

        str_trans_v0 = get(handles.v0, 'String');
        trans_v0 = ['v0 = ',str_trans_v0,';                                                                     % transducer vibration velocity amplitude, m/s, assuming v0=1 m/s'];
        replace_parameters_for_trans_v0(trans_v0);

        str_trans_dt = get(handles.dt, 'String');
        trans_dt = ['inter_dist = ',str_trans_dt,';                                                         % distance between the transducer and the particle center (i.e., the origin of the coordinate system)'];
        replace_parameters_for_trans_dt(trans_dt);

        % replace transducer parameters into "phase_array_beam_shape_coeff.m"
        str_trans_position = get(handles.trans_position, 'String');
        transducer = ['transducer = ',str_trans_position,';'];
        replace_phase_array_for_trans_position(transducer);

        str_trans_phase_delay = get(handles.phase_delay, 'String');
        phi_delay = ['phi_delay = ',str_trans_phase_delay,';                                              % rad'];
        replace_phase_array_for_trans_phi_delay(phi_delay);

        str_trans_amp_delay = get(handles.amp_delay, 'String');
        amp_delay = ['A_delay = ',str_trans_amp_delay,';'];
        replace_phase_array_for_trans_amp_delay(amp_delay);
    otherwise
        % replace wavefront parameters 'p_inlet' into "parameters.m" for plane wave
        str_p0 = get(handles.p0, 'String');
        p_inlet = ['p_inlet = ',str_p0,';               % initial amplitude of proba particle [Pa] if no attenuation effects'];
        replace_parameters_for_p_inlet(p_inlet);
end

% particle rotation status should be [0,0,0] in "parameters.m", since we
% just want to visualize the wavefield along z-axis.
% switch wavefront_str{wavefront_val}
%     case 'Plane travelling wave'
%         par_initial_orientation = str2num([0 0 0]);
%         theta_X = ['    theta_x = ',num2str(par_initial_orientation(1)),';                % particle counter-clock wise rotation (position Tx) along x-axis for positive ''theta_x'''];
%         theta_Y = ['    theta_y = ',num2str(par_initial_orientation(2)),';                % particle counter-clock wise rotation (position Ty) along y-axis for positive ''theta_y'''];
%         theta_Z = ['    theta_z = ',num2str(par_initial_orientation(3)),';                % particle counter-clock wise rotation (position Tz) along z-axis for positive ''theta_z'''];
%         replace_parameters_for_par_initial_orientation(theta_X, theta_Y, theta_Z);
%     case 'Transducer array (circular oscillator)'
%         par_initial_orientation = str2num([0 0 0]);
%         theta_X = ['    theta_x = ',num2str(par_initial_orientation(1)),';                % particle counter-clock wise rotation (position Tx) along x-axis for positive ''theta_x'''];
%         theta_Y = ['    theta_y = ',num2str(par_initial_orientation(2)),';                % particle counter-clock wise rotation (position Ty) along y-axis for positive ''theta_y'''];
%         theta_Z = ['    theta_z = ',num2str(par_initial_orientation(3)),';                % particle counter-clock wise rotation (position Tz) along z-axis for positive ''theta_z'''];
%         replace_parameters_for_par_initial_orientation(theta_X, theta_Y, theta_Z);
%     otherwise
%         par_initial_orientation = str2num([0 0 0]);
%         theta_X = ['    theta_x = ',num2str(par_initial_orientation(1)),';                % particle counter-clock wise rotation (position Tx) along x-axis for positive ''theta_x'''];
%         theta_Y = ['    theta_y = ',num2str(par_initial_orientation(2)),';                % particle counter-clock wise rotation (position Ty) along y-axis for positive ''theta_y'''];
%         theta_Z = ['    theta_z = ',num2str(par_initial_orientation(3)),';                % particle counter-clock wise rotation (position Tz) along z-axis for positive ''theta_z'''];
%         replace_parameters_for_par_initial_orientation(theta_X, theta_Y, theta_Z);
% end

% visualize the wavefield on xOz plane
cla(handles.axes3);
axes(handles.axes3);
visualize_incident_wavefield();
set(handles.wavefield_status, 'String', '');

clc;


% --- Executes on selection change in irregular.
function irregular_Callback(hObject, eventdata, handles)
% hObject    handle to irregular (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns irregular contents as cell array
%        contents{get(hObject,'Value')} returns selected item from irregular

str = get(hObject, 'String');
val = get(hObject, 'Value');
str_particle_radius = get(handles.par_radius, 'String');

switch str{val}
    case 'Spherical particle'
        % readin 'parameters' for 'Cn' and 'particle_radius'
        Cn = '    Cn = [particle_radius, 0];                    % spherical';
        par_radius = ['    particle_radius = ',str_particle_radius,';                   % particle radial mean-radius (a); 30 um;'];
        replace_parameters_for_par_radius(par_radius);
        replace_parameters_for_Cn(Cn);
        parameters;
        set(handles.mapping_coeff, 'String', ['[', num2str(Cn), ']']);
	case 'Ellipsoidal particle'
        % readin 'parameters' for 'Cn' and 'particle_radius'
        Cn = '    Cn = [particle_radius, 0, particle_radius/5];   % spheroidal';
        par_radius = ['    particle_radius = ',str_particle_radius,';                   % particle radial mean-radius (a); 30 um;'];
        replace_parameters_for_par_radius(par_radius);
        replace_parameters_for_Cn(Cn);
        parameters;
        set(handles.mapping_coeff, 'String', ['[', num2str(Cn), ']']);
    case 'Cone particle'
        % readin 'parameters' for 'Cn' and 'particle_radius'
        Cn = '    Cn = [particle_radius, 0, 0, particle_radius/8];   % trangular-cone';
        par_radius = ['    particle_radius = ',str_particle_radius,';                   % particle radial mean-radius (a); 30 um;'];
        replace_parameters_for_par_radius(par_radius);
        replace_parameters_for_Cn(Cn);
        parameters;
        set(handles.mapping_coeff, 'String', ['[', num2str(Cn), ']']);
    case 'Diamond particle'
        % readin 'parameters' for 'Cn' and 'particle_radius'
        Cn = '    Cn = [particle_radius, 0, 0, 0, particle_radius/10];   % diamond';
        par_radius = ['    particle_radius = ',str_particle_radius,';                   % particle radial mean-radius (a); 30 um;'];
        replace_parameters_for_par_radius(par_radius);
        replace_parameters_for_Cn(Cn);
        parameters;
        set(handles.mapping_coeff, 'String', ['[', num2str(Cn), ']']);
    case 'Others (user-specified ''Cn'')'
        % self specified the 'Cn'
        par_radius = ['    particle_radius = ',str_particle_radius,';                   % particle radial mean-radius (a); 30 um;'];
        replace_parameters_for_par_radius(par_radius);
        set(handles.mapping_coeff, 'String', ['[', str_particle_radius, '  ']);
    otherwise
        % readin 'parameters' for 'Cn' and 'particle_radius'
        Cn = '    Cn = [particle_radius, 0];                    % spherical';
        par_radius = ['    particle_radius = ',str_particle_radius,';                   % particle radial mean-radius (a); 30 um;'];
        replace_parameters_for_par_radius(par_radius);
        replace_parameters_for_Cn(Cn);
        parameters;
        set(handles.mapping_coeff, 'String', ['[', num2str(Cn), ']']);
end



% --- Executes during object creation, after setting all properties.
function irregular_CreateFcn(hObject, eventdata, handles)
% hObject    handle to irregular (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mapping_coeff_Callback(hObject, eventdata, handles)
% hObject    handle to mapping_coeff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mapping_coeff as text
%        str2double(get(hObject,'String')) returns contents of mapping_coeff as a double


% --- Executes during object creation, after setting all properties.
function mapping_coeff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mapping_coeff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function par_density_Callback(hObject, eventdata, handles)
% hObject    handle to par_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of par_density as text
%        str2double(get(hObject,'String')) returns contents of par_density as a double


% --- Executes during object creation, after setting all properties.
function par_density_CreateFcn(hObject, eventdata, handles)
% hObject    handle to par_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in par_visual.
function par_visual_Callback(hObject, eventdata, handles)
% hObject    handle to par_visual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% replace mapping coefficients into "parameters.m"
str_Cn = get(handles.mapping_coeff, 'String');
Cn = ['    Cn = ', str_Cn, ';'];
replace_parameters_for_Cn(Cn);

% replace particle parameters into "parameters.m"
str_par_density = get(handles.par_density, 'String');
par_density = ['particle_rho = ',str_par_density,';'];
replace_parameters_for_par_density(par_density);

str_particle_radius = get(handles.par_radius, 'String');
par_radius = ['    particle_radius = ',str_particle_radius,';                   % particle radial mean-radius (a); 30 um;'];
replace_parameters_for_par_radius(par_radius);

str_par_initial_position = get(handles.initial_position, 'String');
par_initial_position = str2num(str_par_initial_position);
positionX = ['deviationX = ',num2str(par_initial_position(1)),';'];
positionY = ['deviationY = ',num2str(par_initial_position(2)),';'];
positionZ = ['deviationZ = ',num2str(par_initial_position(3)),';'];
replace_parameters_for_par_initial_position(positionX, positionY, positionZ);

str_par_initial_orientation = get(handles.initial_orientation, 'String');
par_initial_orientation = str2num(str_par_initial_orientation);
theta_X = ['    theta_x = ',num2str(par_initial_orientation(1)),';                % particle counter-clock wise rotation (position Tx) along x-axis for positive ''theta_x'''];
theta_Y = ['    theta_y = ',num2str(par_initial_orientation(2)),';                % particle counter-clock wise rotation (position Ty) along y-axis for positive ''theta_y'''];
theta_Z = ['    theta_z = ',num2str(par_initial_orientation(3)),';                % particle counter-clock wise rotation (position Tz) along z-axis for positive ''theta_z'''];
replace_parameters_for_par_initial_orientation(theta_X, theta_Y, theta_Z);


% visualize the particle in 'axes1'
cla(handles.axes1);
axes(handles.axes1);
establish_3D_model(100);


function par_radius_Callback(hObject, eventdata, handles)
% hObject    handle to par_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of par_radius as text
%        str2double(get(hObject,'String')) returns contents of par_radius as a double


% --- Executes during object creation, after setting all properties.
function par_radius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to par_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function initial_position_Callback(hObject, eventdata, handles)
% hObject    handle to initial_position (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of initial_position as text
%        str2double(get(hObject,'String')) returns contents of initial_position as a double


% --- Executes during object creation, after setting all properties.
function initial_position_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initial_position (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function initial_orientation_Callback(hObject, eventdata, handles)
% hObject    handle to initial_orientation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of initial_orientation as text
%        str2double(get(hObject,'String')) returns contents of initial_orientation as a double


% --- Executes during object creation, after setting all properties.
function initial_orientation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initial_orientation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function delta_t_Callback(hObject, eventdata, handles)
% hObject    handle to delta_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delta_t as text
%        str2double(get(hObject,'String')) returns contents of delta_t as a double


% --- Executes during object creation, after setting all properties.
function delta_t_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delta_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t_end_Callback(hObject, eventdata, handles)
% hObject    handle to t_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_end as text
%        str2double(get(hObject,'String')) returns contents of t_end as a double


% --- Executes during object creation, after setting all properties.
function t_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rad_force_Callback(hObject, eventdata, handles)
% hObject    handle to rad_force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rad_force as text
%        str2double(get(hObject,'String')) returns contents of rad_force as a double


% --- Executes during object creation, after setting all properties.
function rad_force_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rad_force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rad_torque_Callback(hObject, eventdata, handles)
% hObject    handle to rad_torque (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rad_torque as text
%        str2double(get(hObject,'String')) returns contents of rad_torque as a double


% --- Executes during object creation, after setting all properties.
function rad_torque_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rad_torque (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_num_Callback(hObject, eventdata, handles)
% hObject    handle to trans_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_num as text
%        str2double(get(hObject,'String')) returns contents of trans_num as a double


% --- Executes during object creation, after setting all properties.
function trans_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_radius_Callback(hObject, eventdata, handles)
% hObject    handle to trans_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_radius as text
%        str2double(get(hObject,'String')) returns contents of trans_radius as a double


% --- Executes during object creation, after setting all properties.
function trans_radius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function v0_Callback(hObject, eventdata, handles)
% hObject    handle to v0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of v0 as text
%        str2double(get(hObject,'String')) returns contents of v0 as a double


% --- Executes during object creation, after setting all properties.
function v0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dt_Callback(hObject, eventdata, handles)
% hObject    handle to dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dt as text
%        str2double(get(hObject,'String')) returns contents of dt as a double


% --- Executes during object creation, after setting all properties.
function dt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_position_Callback(hObject, eventdata, handles)
% hObject    handle to trans_position (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_position as text
%        str2double(get(hObject,'String')) returns contents of trans_position as a double


% --- Executes during object creation, after setting all properties.
function trans_position_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_position (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function phase_delay_Callback(hObject, eventdata, handles)
% hObject    handle to phase_delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phase_delay as text
%        str2double(get(hObject,'String')) returns contents of phase_delay as a double


% --- Executes during object creation, after setting all properties.
function phase_delay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phase_delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function amp_delay_Callback(hObject, eventdata, handles)
% hObject    handle to amp_delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of amp_delay as text
%        str2double(get(hObject,'String')) returns contents of amp_delay as a double


% --- Executes during object creation, after setting all properties.
function amp_delay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to amp_delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in trans_visual.
function trans_visual_Callback(hObject, eventdata, handles)
% hObject    handle to trans_visual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% replace mapping coefficients into "parameters.m"
str_Cn = get(handles.mapping_coeff, 'String');
Cn = ['    Cn = ', str_Cn, ';'];
replace_parameters_for_Cn(Cn);

% replace particle parameters into "parameters.m"
str_par_density = get(handles.par_density, 'String');
par_density = ['particle_rho = ',str_par_density,';'];
replace_parameters_for_par_density(par_density);

str_particle_radius = get(handles.par_radius, 'String');
par_radius = ['    particle_radius = ',str_particle_radius,';                   % particle radial mean-radius (a); 30 um;'];
replace_parameters_for_par_radius(par_radius);

str_par_initial_position = get(handles.initial_position, 'String');
par_initial_position = str2num(str_par_initial_position);
positionX = ['deviationX = ',num2str(par_initial_position(1)),';'];
positionY = ['deviationY = ',num2str(par_initial_position(2)),';'];
positionZ = ['deviationZ = ',num2str(par_initial_position(3)),';'];
replace_parameters_for_par_initial_position(positionX, positionY, positionZ);

str_par_initial_orientation = get(handles.initial_orientation, 'String');
par_initial_orientation = str2num(str_par_initial_orientation);
theta_X = ['    theta_x = ',num2str(par_initial_orientation(1)),';                % particle counter-clock wise rotation (position Tx) along x-axis for positive ''theta_x'''];
theta_Y = ['    theta_y = ',num2str(par_initial_orientation(2)),';                % particle counter-clock wise rotation (position Ty) along y-axis for positive ''theta_y'''];
theta_Z = ['    theta_z = ',num2str(par_initial_orientation(3)),';                % particle counter-clock wise rotation (position Tz) along z-axis for positive ''theta_z'''];
replace_parameters_for_par_initial_orientation(theta_X, theta_Y, theta_Z);

% visualize unrotated and untranslated particle in axes1
cla(handles.axes1);
axes(handles.axes1);
[X,Y,Z] = establish_3D_model(100);


% rotate and translate the particle based on 'theta_rotation'
theta_rotation = str2num(get(handles.initial_orientation, 'String'));
translation = str2num(get(handles.initial_position, 'String'));
Rx=[1 0 0;
    0 cos(-theta_rotation(1)) -sin(-theta_rotation(1)); 
    0 sin(-theta_rotation(1)) cos(-theta_rotation(1))];
Ry=[cos(-theta_rotation(2)) 0 sin(-theta_rotation(2));
    0 1 0; 
    -sin(-theta_rotation(2)) 0 cos(-theta_rotation(2))];
Rz=[cos(-theta_rotation(3)) -sin(-theta_rotation(3)) 0; 
    sin(-theta_rotation(3)) cos(-theta_rotation(3)) 0; 
    0 0 1];
Rxyz = Rx * Ry * Rz;
[r, c] = size(X);
for ii = 1 : r
    for jj = 1 : c
    
        T = [X(ii,jj), Y(ii,jj), Z(ii,jj)] * Rxyz;
        X(ii,jj) = T(1);
        Y(ii,jj) = T(2);
        Z(ii,jj) = T(3);
        X(ii,jj) = X(ii,jj) + translation(1);
        Y(ii,jj) = Y(ii,jj) + translation(2);
        Z(ii,jj) = Z(ii,jj) + translation(3);
        
    end
end

% visualize the rotated and translated particle based on GUI 'initial orientation'
cla(handles.axes2);
axes(handles.axes2);
col = ones(r,c,3);
col(:,:,1) = 0.5;
col(:,:,2) = 0.5;
col(:,:,3) = 0.5; % 'col' for gray color
h = surf(X, Y, Z, col);
set(h,'edgecolor','none');



% replace transducer parameters into "parameters.m"
str_trans_num = get(handles.trans_num, 'String');
trans_num = ['transducer_number = ',str_trans_num,';          % NOTE: keep same in "total_effects_all_plain_wave_component.m". for phase array'];
replace_parameters_for_trans_num(trans_num);

str_trans_radius = get(handles.trans_radius, 'String');
trans_radius = ['trans_radius = ',str_trans_radius,';           % NOTE: keep same in "total_effects_all_plain_wave_component.m". radius of transducer is 5mm'];
replace_parameters_for_trans_radius(trans_radius);

str_trans_v0 = get(handles.v0, 'String');
trans_v0 = ['v0 = ',str_trans_v0,';                                                                     % transducer vibration velocity amplitude, m/s, assuming v0=1 m/s'];
replace_parameters_for_trans_v0(trans_v0);

str_trans_dt = get(handles.dt, 'String');
trans_dt = ['inter_dist = ',str_trans_dt,';                                                         % distance between the transducer and the particle center (i.e., the origin of the coordinate system)'];
replace_parameters_for_trans_dt(trans_dt);

% replace transducer parameters into "phase_array_beam_shape_coeff.m"
str_trans_position = get(handles.trans_position, 'String');
transducer = ['transducer = ',str_trans_position,';'];
replace_phase_array_for_trans_position(transducer);

str_trans_phase_delay = get(handles.phase_delay, 'String');
phi_delay = ['phi_delay = ',str_trans_phase_delay,';                                              % rad'];
replace_phase_array_for_trans_phi_delay(phi_delay);

str_trans_amp_delay = get(handles.amp_delay, 'String');
amp_delay = ['A_delay = ',str_trans_amp_delay,';'];
replace_phase_array_for_trans_amp_delay(amp_delay);


% visualize the transducer based on GUI 'transducer positions'
transducer = str2num(get(handles.trans_position, 'String'));
visualize_transducers(transducer);

xlabel('\rm{\fontname{Times new roman}\it{x}{\rm{-axis [mm]}}}');
ylabel('\rm{\fontname{Times new roman}\it{y}{\rm{-axis [mm]}}}');
zlabel('\rm{\fontname{Times new roman}\it{z}{\rm{-axis [mm]}}}');
set(gca,'XTick', [-0.01, 0, 0.01]);
set(gca,'XTickLabel', [-10, 0, 10]);
set(gca,'YTick', [-0.01, 0, 0.01]);
set(gca,'YTickLabel', [-10, 0, 10]);
set(gca,'ZTick', [-0.06, -0.04, -0.02, -0.00, 0.02]);
set(gca,'ZTickLabel', [-60, -40, -20, 0, 20]);
set(gca, 'FontName', 'Times new roman');
% set(get(gca,'XLabel'),'FontSize',figure_FontSize);
% set(get(gca,'YLabel'),'FontSize',figure_FontSize);
% set(get(gca,'ZLabel'),'FontSize',figure_FontSize);
% set(findobj('Fontsize',10),'fontsize', figure_FontSize-5);
axis image; camlight; lighting gouraud;
axis([-0.015, 0.015, -0.015, 0.015, -0.06, 0.02]);





function p0_Callback(hObject, eventdata, handles)
% hObject    handle to p0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p0 as text
%        str2double(get(hObject,'String')) returns contents of p0 as a double


% --- Executes during object creation, after setting all properties.
function p0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calculation_dynamics.
function calculation_dynamics_Callback(hObject, eventdata, handles)
% hObject    handle to calculation_dynamics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% replace the frequency in "parameters.m"
str_freq = get(handles.freq, 'String');
freq = ['freq = ',str_freq,';               % 1 MHz'];
replace_parameters_for_frequency(freq);

str_integration_radius = get(handles.integration_radius, 'String');
integration_radius = ['b = ',str_integration_radius,';           % integral radial distance (b) for beam-shape coeffcient, EMPIRICAL: MUST follow  b > a and b < inter_dist BUT not b >> a'];
replace_parameters_for_integration_radius(integration_radius);

% replace wavefront parameters into "parameters.m" and "phase_array_beam_shape_coeff.m"
wavefront_str = get(handles.wavefront_selection, 'String');
wavefront_val = get(handles.wavefront_selection, 'Value');
switch wavefront_str{wavefront_val}
    case 'Plane travelling wave'
        % replace wavefront parameters 'p_inlet' into "parameters.m" for plane wave
        str_p0 = get(handles.p0, 'String');
        p_inlet = ['p_inlet = ',str_p0,';               % initial amplitude of proba particle [Pa] if no attenuation effects'];
        replace_parameters_for_p_inlet(p_inlet);
    case 'Transducer array (circular oscillator)'
        % replace wavefront parameters into "parameters.m" and
        % "phase_array_beam_shape_coeff.m" for transducer array
        str_trans_num = get(handles.trans_num, 'String');
        trans_num = ['transducer_number = ',str_trans_num,';          % NOTE: keep same in "total_effects_all_plain_wave_component.m". for phase array'];
        replace_parameters_for_trans_num(trans_num);

        str_trans_radius = get(handles.trans_radius, 'String');
        trans_radius = ['trans_radius = ',str_trans_radius,';           % NOTE: keep same in "total_effects_all_plain_wave_component.m". radius of transducer is 5mm'];
        replace_parameters_for_trans_radius(trans_radius);

        str_trans_v0 = get(handles.v0, 'String');
        trans_v0 = ['v0 = ',str_trans_v0,';                                                                     % transducer vibration velocity amplitude, m/s, assuming v0=1 m/s'];
        replace_parameters_for_trans_v0(trans_v0);

        str_trans_dt = get(handles.dt, 'String');
        trans_dt = ['inter_dist = ',str_trans_dt,';                                                         % distance between the transducer and the particle center (i.e., the origin of the coordinate system)'];
        replace_parameters_for_trans_dt(trans_dt);

        % replace transducer parameters into "phase_array_beam_shape_coeff.m"
        str_trans_position = get(handles.trans_position, 'String');
        transducer = ['transducer = ',str_trans_position,';'];
        replace_phase_array_for_trans_position(transducer);

        str_trans_phase_delay = get(handles.phase_delay, 'String');
        phi_delay = ['phi_delay = ',str_trans_phase_delay,';                                              % rad'];
        replace_phase_array_for_trans_phi_delay(phi_delay);

        str_trans_amp_delay = get(handles.amp_delay, 'String');
        amp_delay = ['A_delay = ',str_trans_amp_delay,';'];
        replace_phase_array_for_trans_amp_delay(amp_delay);
    otherwise
        % replace wavefront parameters 'p_inlet' into "parameters.m" for plane wave
        str_p0 = get(handles.p0, 'String');
        p_inlet = ['p_inlet = ',str_p0,';               % initial amplitude of proba particle [Pa] if no attenuation effects'];
        replace_parameters_for_p_inlet(p_inlet);
end

% replace mapping coefficients into "parameters.m"
str_Cn = get(handles.mapping_coeff, 'String');
Cn = ['    Cn = ', str_Cn, ';'];
replace_parameters_for_Cn(Cn);

% replace particle parameters into "parameters.m"
str_par_density = get(handles.par_density, 'String');
par_density = ['particle_rho = ',str_par_density,';'];
replace_parameters_for_par_density(par_density);

str_particle_radius = get(handles.par_radius, 'String');
par_radius = ['    particle_radius = ',str_particle_radius,';                   % particle radial mean-radius (a); 30 um;'];
replace_parameters_for_par_radius(par_radius);

str_par_initial_position = get(handles.initial_position, 'String');
par_initial_position = str2num(str_par_initial_position);
positionX = ['deviationX = ',num2str(par_initial_position(1)),';'];
positionY = ['deviationY = ',num2str(par_initial_position(2)),';'];
positionZ = ['deviationZ = ',num2str(par_initial_position(3)),';'];
replace_parameters_for_par_initial_position(positionX, positionY, positionZ);

str_par_initial_orientation = get(handles.initial_orientation, 'String');
par_initial_orientation = str2num(str_par_initial_orientation);
theta_X = ['    theta_x = ',num2str(par_initial_orientation(1)),';                % particle counter-clock wise rotation (position Tx) along x-axis for positive ''theta_x'''];
theta_Y = ['    theta_y = ',num2str(par_initial_orientation(2)),';                % particle counter-clock wise rotation (position Ty) along y-axis for positive ''theta_y'''];
theta_Z = ['    theta_z = ',num2str(par_initial_orientation(3)),';                % particle counter-clock wise rotation (position Tz) along z-axis for positive ''theta_z'''];
replace_parameters_for_par_initial_orientation(theta_X, theta_Y, theta_Z);

% run dynamic function and visualize the dynamics graphs in GUI for
% Transducer array ONLY.
X_vec = (par_initial_position(1));
Y_vec = (par_initial_position(2));
Z_vec = (par_initial_position(3));
theta_x = (par_initial_orientation(1));
theta_y = (par_initial_orientation(2));
theta_z = (par_initial_orientation(3));
switch wavefront_str{wavefront_val}
    case 'Plane travelling wave'
        set(handles.acoustophoresis_status, 'String', 'Array Only!');
    case 'Transducer array (circular oscillator)'
        set(handles.acoustophoresis_status, 'String', 'Busy');pause(0.001);
        delta_t = str2double(get(handles.delta_t, 'String'));
        ending_t = str2double(get(handles.t_end, 'String'));
        [time, X_position, Y_position, Z_position, X_rotation, Y_rotation, Z_rotation] = ...
            dynamics_in_array(X_vec, Y_vec, Z_vec, theta_x, theta_y, theta_z, delta_t, ending_t, handles);
        set(handles.acoustophoresis_status, 'String', '');
%         handles.time = time;
%         handles.X_position = X_position;
%         handles.Y_position = Y_position;
%         handles.Z_position = Z_position;
%         handles.X_rotation = X_rotation;
%         handles.Y_rotation = Y_rotation;
%         handles.Z_rotation = Z_rotation;
    otherwise
        set(handles.acoustophoresis_status, 'String', 'Array Only!');
end




% % % % % Saving Data % % % % %

switch wavefront_str{wavefront_val}
    
    case 'Transducer array (circular oscillator)'
        
        % get mapping coefficients into "parameters.m"
        str_geometry = get(handles.irregular, 'String');
        geometry_val = get(handles.irregular, 'Value');
        str_Cn = get(handles.mapping_coeff, 'String');

        % get particle parameters into "parameters.m"
        str_par_density = get(handles.par_density, 'String');
        str_particle_radius = get(handles.par_radius, 'String');
        str_par_initial_position = get(handles.initial_position, 'String');
        str_par_initial_orientation = get(handles.initial_orientation, 'String');

        % get boundary condition
        str_BC = get(handles.boundary_conditions, 'String');
        BC_val = get(handles.boundary_conditions, 'Value');

        % get medium
        str_medium = get(handles.medium_type, 'String');
        medium_val = get(handles.medium_type, 'Value');
        str_fluid_rho = get(handles.fluid_rho, 'String');
        str_fluid_c = get(handles.fluid_c, 'String');
        str_fluid_vis = get(handles.fluid_vis, 'String');

        % get the frequency in "parameters.m"
        str_freq = get(handles.freq, 'String');

        str_integration_radius = get(handles.integration_radius, 'String');

        % get wavefront parameters into "parameters.m" and "phase_array_beam_shape_coeff.m"
        wavefront_str = get(handles.wavefront_selection, 'String');
        wavefront_val = get(handles.wavefront_selection, 'Value');
        
        % get wavefront parameters into "parameters.m" and
        % "phase_array_beam_shape_coeff.m" for transducer array
        str_trans_num = get(handles.trans_num, 'String');
        str_trans_radius = get(handles.trans_radius, 'String');
        str_trans_v0 = get(handles.v0, 'String');
        str_trans_dt = get(handles.dt, 'String');

        % get transducer parameters into "phase_array_beam_shape_coeff.m"
        str_trans_position = get(handles.trans_position, 'String');
        str_trans_phase_delay = get(handles.phase_delay, 'String');
        str_trans_amp_delay = get(handles.amp_delay, 'String');
            
        % acoustophorestic data
%         time = handles.time;
%         X_position = handles.X_position;
%         Y_position = handles.Y_position;
%         Z_position = handles.Z_position;
%         X_rotation = handles.X_rotation;
%         Y_rotation = handles.Y_rotation;
%         Z_rotation = handles.Z_rotation;


        % Write to file
        % [f_name,fpath,filterindex] = uiputfile({'.txt','Text File (*.txt)'},'Save data','MyFilename');
        f_name=(get(handles.f_name,'String'));

        % Check if file exists and warn user that it will be overwritten
        if(fopen(f_name, 'r') ~= -1)
            if strcmp(questdlg('This file already exists. Overwrite?', 'Overwrite Saved Simulation?', 'Yes', 'No', 'No'), 'No')
                return;
            end
        end

        fid = fopen(f_name, 'w');

        fprintf(fid, '%%Computational parameters\n\n');

        fprintf(fid, 'Geometry = %s;\n', str_geometry{geometry_val});
        fprintf(fid, '    Cn = %s;\n', str_Cn);
        fprintf(fid, '    Particle density = %s [kg/m^3];\n', str_par_density);
        fprintf(fid, '    Particle radius = %s [m];\n', str_particle_radius);
        fprintf(fid, '    Initial position = %s [m];\n', str_par_initial_position);
        fprintf(fid, '    Initial orientation = %s [rad];\n\n', str_par_initial_orientation);

        fprintf(fid, 'Boundary condition = %s;\n\n', str_BC{BC_val});

        fprintf(fid, 'Medium = %s;\n', str_medium{medium_val});
        fprintf(fid, '    Fluid density = %s [kg/m^3];\n', str_fluid_rho);
        fprintf(fid, '    Fluid sound speed = %s [m/s];\n', str_fluid_c);
        fprintf(fid, '    Fluid dynamic viscosity = %s [Pa s];\n\n', str_fluid_vis);

        fprintf(fid, 'Wave type = %s;\n', wavefront_str{wavefront_val});
        fprintf(fid, '    Frequency = %s [Hz];\n', str_freq);
        fprintf(fid, '    Integration radius = %s [m];\n', str_integration_radius);
        fprintf(fid, '    Transducer number = %s;\n', str_trans_num);
        fprintf(fid, '    Transducer radius = %s [m];\n', str_trans_radius);
        fprintf(fid, '    Transducer vibration radial velocity = %s [m/s];\n', str_trans_v0);
        fprintf(fid, '    Interdistance between the particle and the array = %s [m];\n', str_trans_dt);
        fprintf(fid, '    Transducer position matrix = %s [m];\n', str_trans_position);
        fprintf(fid, '    Transducer phase delay = %s [rad];\n', str_trans_phase_delay);
        fprintf(fid, '    Transducer relative amplitude delay = %s;\n\n', str_trans_amp_delay);

        
        fprintf(fid, '\n\n%%Dynamic data:\n\n');
        
        fprintf(fid, 'Time [ms]     X_position [m]     Y_position [m]     Z_position [m]     X_angle [rad]     Y_angle [rad]     Z_angle [rad]\n');
        for ii = 1 : length(time)
            fprintf(fid, '%f  %f  %f  %f  %f  %f  %f  \n', time(ii)*1000, X_position(ii), Y_position(ii), Z_position(ii), X_rotation(ii), Y_rotation(ii), Z_rotation(ii));
        end

        fclose(fid);
        
    otherwise
        
end


% --- Executes on button press in calculation_force_torque.
function calculation_force_torque_Callback(hObject, eventdata, handles)
% hObject    handle to calculation_force_torque (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% replace the frequency in "parameters.m"
str_freq = get(handles.freq, 'String');
freq = ['freq = ',str_freq,';               % 1 MHz'];
replace_parameters_for_frequency(freq);

str_integration_radius = get(handles.integration_radius, 'String');
integration_radius = ['b = ',str_integration_radius,';           % integral radial distance (b) for beam-shape coeffcient, EMPIRICAL: MUST follow  b > a and b < inter_dist BUT not b >> a'];
replace_parameters_for_integration_radius(integration_radius);

% replace wavefront parameters into "parameters.m" and "phase_array_beam_shape_coeff.m"
wavefront_str = get(handles.wavefront_selection, 'String');
wavefront_val = get(handles.wavefront_selection, 'Value');
switch wavefront_str{wavefront_val}
    case 'Plane travelling wave'
        % replace wavefront parameters 'p_inlet' into "parameters.m" for plane wave
        str_p0 = get(handles.p0, 'String');
        p_inlet = ['p_inlet = ',str_p0,';               % initial amplitude of proba particle [Pa] if no attenuation effects'];
        replace_parameters_for_p_inlet(p_inlet);
    case 'Transducer array (circular oscillator)'
        % replace wavefront parameters into "parameters.m" and
        % "phase_array_beam_shape_coeff.m" for transducer array
        str_trans_num = get(handles.trans_num, 'String');
        trans_num = ['transducer_number = ',str_trans_num,';          % NOTE: keep same in "total_effects_all_plain_wave_component.m". for phase array'];
        replace_parameters_for_trans_num(trans_num);

        str_trans_radius = get(handles.trans_radius, 'String');
        trans_radius = ['trans_radius = ',str_trans_radius,';           % NOTE: keep same in "total_effects_all_plain_wave_component.m". radius of transducer is 5mm'];
        replace_parameters_for_trans_radius(trans_radius);

        str_trans_v0 = get(handles.v0, 'String');
        trans_v0 = ['v0 = ',str_trans_v0,';                                                                     % transducer vibration velocity amplitude, m/s, assuming v0=1 m/s'];
        replace_parameters_for_trans_v0(trans_v0);

        str_trans_dt = get(handles.dt, 'String');
        trans_dt = ['inter_dist = ',str_trans_dt,';                                                         % distance between the transducer and the particle center (i.e., the origin of the coordinate system)'];
        replace_parameters_for_trans_dt(trans_dt);

        % replace transducer parameters into "phase_array_beam_shape_coeff.m"
        str_trans_position = get(handles.trans_position, 'String');
        transducer = ['transducer = ',str_trans_position,';'];
        replace_phase_array_for_trans_position(transducer);

        str_trans_phase_delay = get(handles.phase_delay, 'String');
        phi_delay = ['phi_delay = ',str_trans_phase_delay,';                                              % rad'];
        replace_phase_array_for_trans_phi_delay(phi_delay);

        str_trans_amp_delay = get(handles.amp_delay, 'String');
        amp_delay = ['A_delay = ',str_trans_amp_delay,';'];
        replace_phase_array_for_trans_amp_delay(amp_delay);
    otherwise
        % replace wavefront parameters 'p_inlet' into "parameters.m" for plane wave
        str_p0 = get(handles.p0, 'String');
        p_inlet = ['p_inlet = ',str_p0,';               % initial amplitude of proba particle [Pa] if no attenuation effects'];
        replace_parameters_for_p_inlet(p_inlet);
end


% replace mapping coefficients into "parameters.m"
str_Cn = get(handles.mapping_coeff, 'String');
Cn = ['    Cn = ', str_Cn, ';'];
replace_parameters_for_Cn(Cn);

% replace particle parameters into "parameters.m"
str_par_density = get(handles.par_density, 'String');
par_density = ['particle_rho = ',str_par_density,';'];
replace_parameters_for_par_density(par_density);

str_particle_radius = get(handles.par_radius, 'String');
par_radius = ['    particle_radius = ',str_particle_radius,';                   % particle radial mean-radius (a); 30 um;'];
replace_parameters_for_par_radius(par_radius);

str_par_initial_position = get(handles.initial_position, 'String');
par_initial_position = str2num(str_par_initial_position);
positionX = ['deviationX = ',num2str(par_initial_position(1)),';'];
positionY = ['deviationY = ',num2str(par_initial_position(2)),';'];
positionZ = ['deviationZ = ',num2str(par_initial_position(3)),';'];
replace_parameters_for_par_initial_position(positionX, positionY, positionZ);

str_par_initial_orientation = get(handles.initial_orientation, 'String');
par_initial_orientation = str2num(str_par_initial_orientation);
theta_X = ['    theta_x = ',num2str(par_initial_orientation(1)),';                % particle counter-clock wise rotation (position Tx) along x-axis for positive ''theta_x'''];
theta_Y = ['    theta_y = ',num2str(par_initial_orientation(2)),';                % particle counter-clock wise rotation (position Ty) along y-axis for positive ''theta_y'''];
theta_Z = ['    theta_z = ',num2str(par_initial_orientation(3)),';                % particle counter-clock wise rotation (position Tz) along z-axis for positive ''theta_z'''];
replace_parameters_for_par_initial_orientation(theta_X, theta_Y, theta_Z);


% call the functions "radiation_force_based_Analyses.m" and
% "radiation_torque_based_Analyses.m" for the radiation force and torque.
set(handles.force_torque_status, 'String', 'Busy');
pause(0.001);
[Frad_x,Frad_y,Frad_z,~] = radiation_force_based_Analyses();
[Torque_x,Torque_y,Torque_z,~] = radiation_torque_based_Analyses();
set(handles.force_torque_status, 'String', '');

% visualize the force and torque in GUI
set(handles.rad_force, 'String', ['[', num2str(Frad_x), '  ', num2str(Frad_y), '  ', num2str(Frad_z), ']']);
set(handles.rad_torque, 'String', ['[', num2str(Torque_x), '  ', num2str(Torque_y), '  ', num2str(Torque_z), ']']);

clc;






% --- Executes on button press in save_figs.
function save_figs_Callback(hObject, eventdata, handles)
% hObject    handle to save_figs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str_fig = get(handles.select_figures, 'String');
val_fig = get(handles.select_figures, 'Value');

switch str_fig{val_fig}
    case 'figure 1'
        [fname, fpath, filterindex] = uiputfile({'.jpeg','JPEG (*.jpeg)';'.png','Portable Network Graphics (*.png)';'.tiff','Tagged Image File Format (*.tiff)';'.eps','EPS File (*.eps)'},'Save Plot','particle');
        if fname ~= 0
            % Create an invisible figure containing just the pressure plot and save that
            h = figure(1);
            set(h,'Visible','Off');
            s = subplot(1,1,1);
            copyobj(allchild(handles.axes1), s);
%             pix = getframe(handles.axes1);
            switch(filterindex)
                case 1
                    saveas(s, fullfile(fpath, fname), 'jpeg');
%                     imwrite(pix.cdata, fullfile(fpath, fname), 'jpeg');
                case 2
                    saveas(s, fullfile(fpath, fname), 'png');
                case 3
                    saveas(s, fullfile(fpath, fname), 'tiffn');
                case 4
                    saveas(s, fullfile(fpath, fname), 'eps');
            end
            delete(h);
        end
	case 'figure 2'
        [fname, fpath, filterindex] = uiputfile({'.jpeg','JPEG (*.jpeg)';'.png','Portable Network Graphics (*.png)';'.tiff','Tagged Image File Format (*.tiff)'},'Save Plot','transducer-particle_system');
        if fname ~= 0
            % Create an invisible figure containing just the pressure plot and save that
            h = figure;
            set(h,'Visible','Off');
            s = subplot(1,1,1);
            copyobj(allchild(handles.axes2), s);
            switch(filterindex)
                case 1
                    saveas(s, fullfile(fpath, fname), 'jpeg');
                case 2
                    saveas(s, fullfile(fpath, fname), 'png');
                case 3
                    saveas(s, fullfile(fpath, fname), 'tiffn');
                case 4
                    saveas(s, fullfile(fpath, fname), 'eps');
            end
            delete(h);
        end
    case 'figure 3'
        [fname, fpath, filterindex] = uiputfile({'.jpeg','JPEG (*.jpeg)';'.png','Portable Network Graphics (*.png)';'.tiff','Tagged Image File Format (*.tiff)'},'Save Plot','pressure-wavefield_(xOz plane)');
        if fname ~= 0
            % Create an invisible figure containing just the pressure plot and save that
            h = figure;
            set(h,'Visible','Off');
            s = subplot(1,1,1);
            copyobj(allchild(handles.axes3), s);
            switch(filterindex)
                case 1
                    saveas(s, fullfile(fpath, fname), 'jpeg');
                case 2
                    saveas(s, fullfile(fpath, fname), 'png');
                case 3
                    saveas(s, fullfile(fpath, fname), 'tiffn');
                case 4
                    saveas(s, fullfile(fpath, fname), 'eps');        
            end
            delete(h);
        end
    case 'figure 4'
        [fname, fpath, filterindex] = uiputfile({'.jpeg','JPEG (*.jpeg)';'.png','Portable Network Graphics (*.png)';'.tiff','Tagged Image File Format (*.tiff)'},'Save Plot','Particle_Xposition');
        if fname ~= 0
            % Create an invisible figure containing just the pressure plot and save that
            h = figure;
            set(h,'Visible','Off');
            s = subplot(1,1,1);
            copyobj(allchild(handles.axes4), s);
            switch(filterindex)
                case 1
                    saveas(s, fullfile(fpath, fname), 'jpeg');
                case 2
                    saveas(s, fullfile(fpath, fname), 'png');
                case 3
                    saveas(s, fullfile(fpath, fname), 'tiffn');
                case 4
                    saveas(s, fullfile(fpath, fname), 'eps');        
            end
            delete(h);
        end
    case 'figure 5'
        [fname, fpath, filterindex] = uiputfile({'.jpeg','JPEG (*.jpeg)';'.png','Portable Network Graphics (*.png)';'.tiff','Tagged Image File Format (*.tiff)'},'Save Plot','Particle_Yposition');
        if fname ~= 0
            % Create an invisible figure containing just the pressure plot and save that
            h = figure;
            set(h,'Visible','Off');
            s = subplot(1,1,1);
            copyobj(allchild(handles.axes5), s);
            switch(filterindex)
                case 1
                    saveas(s, fullfile(fpath, fname), 'jpeg');
                case 2
                    saveas(s, fullfile(fpath, fname), 'png');
                case 3
                    saveas(s, fullfile(fpath, fname), 'tiffn');
                case 4
                    saveas(s, fullfile(fpath, fname), 'eps');        
            end
            delete(h);
        end
    case 'figure 6'
        [fname, fpath, filterindex] = uiputfile({'.jpeg','JPEG (*.jpeg)';'.png','Portable Network Graphics (*.png)';'.tiff','Tagged Image File Format (*.tiff)'},'Save Plot','Particle_Zposition');
        if fname ~= 0
            % Create an invisible figure containing just the pressure plot and save that
            h = figure;
            set(h,'Visible','Off');
            s = subplot(1,1,1);
            copyobj(allchild(handles.axes6), s);
            switch(filterindex)
                case 1
                    saveas(s, fullfile(fpath, fname), 'jpeg');
                case 2
                    saveas(s, fullfile(fpath, fname), 'png');
                case 3
                    saveas(s, fullfile(fpath, fname), 'tiffn');
                case 4
                    saveas(s, fullfile(fpath, fname), 'eps');        
            end
            delete(h);
        end
    case 'figure 7'
        [fname, fpath, filterindex] = uiputfile({'.jpeg','JPEG (*.jpeg)';'.png','Portable Network Graphics (*.png)';'.tiff','Tagged Image File Format (*.tiff)'},'Save Plot','Particle_Xrotation');
        if fname ~= 0
            % Create an invisible figure containing just the pressure plot and save that
            h = figure;
            set(h,'Visible','Off');
            s = subplot(1,1,1);
            copyobj(allchild(handles.axes7), s);
            switch(filterindex)
                case 1
                    saveas(s, fullfile(fpath, fname), 'jpeg');
                case 2
                    saveas(s, fullfile(fpath, fname), 'png');
                case 3
                    saveas(s, fullfile(fpath, fname), 'tiffn');
                case 4
                    saveas(s, fullfile(fpath, fname), 'eps');        
            end
            delete(h);
        end
    case 'figure 8'
        [fname, fpath, filterindex] = uiputfile({'.jpeg','JPEG (*.jpeg)';'.png','Portable Network Graphics (*.png)';'.tiff','Tagged Image File Format (*.tiff)'},'Save Plot','Particle_Yrotation');
        if fname ~= 0
            % Create an invisible figure containing just the pressure plot and save that
            h = figure;
            set(h,'Visible','Off');
            s = subplot(1,1,1);
            copyobj(allchild(handles.axes8), s);
            switch(filterindex)
                case 1
                    saveas(s, fullfile(fpath, fname), 'jpeg');
                case 2
                    saveas(s, fullfile(fpath, fname), 'png');
                case 3
                    saveas(s, fullfile(fpath, fname), 'tiffn');
                case 4
                    saveas(s, fullfile(fpath, fname), 'eps');        
            end
            delete(h);
        end
    case 'figure 9'
        [fname, fpath, filterindex] = uiputfile({'.jpeg','JPEG (*.jpeg)';'.png','Portable Network Graphics (*.png)';'.tiff','Tagged Image File Format (*.tiff)'},'Save Plot','Particle_Zrotation');
        if fname ~= 0
            % Create an invisible figure containing just the pressure plot and save that
            h = figure;
            set(h,'Visible','Off');
            s = subplot(1,1,1);
            copyobj(allchild(handles.axes9), s);
            switch(filterindex)
                case 1
                    saveas(s, fullfile(fpath, fname), 'jpeg');
                case 2
                    saveas(s, fullfile(fpath, fname), 'png');
                case 3
                    saveas(s, fullfile(fpath, fname), 'tiffn');
                case 4
                    saveas(s, fullfile(fpath, fname), 'eps');        
            end
            delete(h);
        end
    case 'figure 10'
        [fname, fpath, filterindex] = uiputfile({'.jpeg','JPEG (*.jpeg)';'.png','Portable Network Graphics (*.png)';'.tiff','Tagged Image File Format (*.tiff)'},'Save Plot','Particle-trajactory');
        if fname ~= 0
            % Create an invisible figure containing just the pressure plot and save that
            h = figure;
            set(h,'Visible','Off');
            s = subplot(1,1,1);
            copyobj(allchild(handles.axes10), s);
            switch(filterindex)
                case 1
                    saveas(s, fullfile(fpath, fname), 'jpeg');
                case 2
                    saveas(s, fullfile(fpath, fname), 'png');
                case 3
                    saveas(s, fullfile(fpath, fname), 'tiffn');
                case 4
                    saveas(s, fullfile(fpath, fname), 'eps');        
            end
            delete(h);
        end
    otherwise
        
end

clc;


% --- Executes on selection change in select_figures.
function select_figures_Callback(hObject, eventdata, handles)
% hObject    handle to select_figures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns select_figures contents as cell array
%        contents{get(hObject,'Value')} returns selected item from select_figures
        
        
% --- Executes during object creation, after setting all properties.
function select_figures_CreateFcn(hObject, eventdata, handles)
% hObject    handle to select_figures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function current_time_Callback(hObject, eventdata, handles)
% hObject    handle to current_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of current_time as text
%        str2double(get(hObject,'String')) returns contents of current_time as a double


% --- Executes during object creation, after setting all properties.
function current_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to current_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_dynamic_data.
function save_dynamic_data_Callback(hObject, eventdata, handles)
% hObject    handle to save_dynamic_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% get mapping coefficients into "parameters.m"
str_geometry = get(handles.irregular, 'String');
geometry_val = get(handles.irregular, 'Value');
str_Cn = get(handles.mapping_coeff, 'String');

% get particle parameters into "parameters.m"
str_par_density = get(handles.par_density, 'String');
str_particle_radius = get(handles.par_radius, 'String');
str_par_initial_position = get(handles.initial_position, 'String');
str_par_initial_orientation = get(handles.initial_orientation, 'String');

% get boundary condition
str_BC = get(handles.boundary_conditions, 'String');
BC_val = get(handles.boundary_conditions, 'Value');

% get medium
str_medium = get(handles.medium_type, 'String');
medium_val = get(handles.medium_type, 'Value');
str_fluid_rho = get(handles.fluid_rho, 'String');
str_fluid_c = get(handles.fluid_c, 'String');
str_fluid_vis = get(handles.fluid_vis, 'String');

% get the frequency in "parameters.m"
str_freq = get(handles.freq, 'String');

str_integration_radius = get(handles.integration_radius, 'String');

% get wavefront parameters into "parameters.m" and "phase_array_beam_shape_coeff.m"
wavefront_str = get(handles.wavefront_selection, 'String');
wavefront_val = get(handles.wavefront_selection, 'Value');
switch wavefront_str{wavefront_val}
    case 'Plane travelling wave'
        % get wavefront parameters 'p_inlet' into "parameters.m" for
        % plane wave
    case 'Transducer array (circular oscillator)'
        % get wavefront parameters into "parameters.m" and
        % "phase_array_beam_shape_coeff.m" for transducer array
        str_trans_num = get(handles.trans_num, 'String');
        str_trans_radius = get(handles.trans_radius, 'String');
        str_trans_v0 = get(handles.v0, 'String');
        str_trans_dt = get(handles.dt, 'String');
        
        % get transducer parameters into "phase_array_beam_shape_coeff.m"
        str_trans_position = get(handles.trans_position, 'String');
        str_trans_phase_delay = get(handles.phase_delay, 'String');
        str_trans_amp_delay = get(handles.amp_delay, 'String');
    otherwise
        % get wavefront parameters 'p_inlet' into "parameters.m" for
        % plane wave
end

% acoustophorestic data
time = handles.time;
X_position = handles.X_position;
Y_position = handles.Y_position;
Z_position = handles.Z_position;
X_rotation = handles.X_rotation;
Y_rotation = handles.Y_rotation;
Z_rotation = handles.Z_rotation;


% Write to file
% [f_name,fpath,filterindex] = uiputfile({'.txt','Text File (*.txt)'},'Save data','MyFilename');
f_name=(get(handles.f_name,'String'));

% Check if file exists and warn user that it will be overwritten
if(fopen(f_name, 'r') ~= -1)
    if strcmp(questdlg('This file already exists. Overwrite?','Overwrite Saved Simulation?','Yes','No','No'), 'No')
        return;
    end
end

fid = fopen(f_name, 'w');

fprintf(fid, '%%Computational parameters\n\n');

fprintf(fid, 'Geometry = %s;\n', str_geometry{geometry_val});
fprintf(fid, '    Cn = %s;\n', str_Cn);
fprintf(fid, '    Particle density = %s [kg/m^3];\n', str_par_density);
fprintf(fid, '    Particle radius = %s [m];\n', str_particle_radius);
fprintf(fid, '    Initial position = %s [m];\n', str_par_initial_position);
fprintf(fid, '    Initial orientation = %s [rad];\n\n', str_par_initial_orientation);

fprintf(fid, 'Boundary condition = %s;\n\n', str_BC{BC_val});

fprintf(fid, 'Medium = %s;\n', str_medium{medium_val});
fprintf(fid, '    Fluid density = %s [kg/m^3];\n', str_fluid_rho);
fprintf(fid, '    Fluid sound speed = %s [m/s];\n', str_fluid_c);
fprintf(fid, '    Fluid dynamic viscosity = %s [Pa s];\n\n', str_fluid_vis);

fprintf(fid, 'Wave type = %s;\n', wavefront_str{wavefront_val});
fprintf(fid, '    Frequency = %s [Hz];\n', str_freq);
fprintf(fid, '    Integration radius = %s [m];\n', str_integration_radius);
fprintf(fid, '    Transducer number = %s;\n', str_trans_num);
fprintf(fid, '    Transducer radius = %s [m];\n', str_trans_radius);
fprintf(fid, '    Transducer vibration radial velocity = %s [m/s];\n', str_trans_v0);
fprintf(fid, '    Interdistance between the particle and the array = %s [m];\n', str_trans_dt);
fprintf(fid, '    Transducer position matrix = %s [m];\n', str_trans_position);
fprintf(fid, '    Transducer phase delay = %s [rad];\n', str_trans_phase_delay);
fprintf(fid, '    Transducer relative amplitude delay = %s;\n\n', str_trans_amp_delay);

fprintf(fid, '\n\n%%Dynamic data:\n\n');

fprintf(fid, 'Time [ms]     X_position [m]     Y_position [m]     Z_position [m]     X_angle [rad]     Y_angle [rad]     Z_angle [rad]\n');
for ii = 1 : length(time)
    fprintf(fid, '%f  %f  %f  %f  %f  %f  %f  \n', time(ii)*1000, X_position(ii), Y_position(ii), Z_position(ii), X_rotation(ii), Y_rotation(ii), Z_rotation(ii));
end

fclose(fid);


        
        
% --- Executes during object creation, after setting all properties.
function uipanel10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function integration_radius_Callback(hObject, eventdata, handles)
% hObject    handle to integration_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of integration_radius as text
%        str2double(get(hObject,'String')) returns contents of integration_radius as a double


% --- Executes during object creation, after setting all properties.
function integration_radius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to integration_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f_name_Callback(hObject, eventdata, handles)
% hObject    handle to f_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f_name as text
%        str2double(get(hObject,'String')) returns contents of f_name as a double


% --- Executes during object creation, after setting all properties.
function f_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
