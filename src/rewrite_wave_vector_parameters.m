function rewrite_wave_vector_parameters(irregular_body, k_Fourier, theta_Fourier, phi_Fourier, trans_radius, transducer_number, fluid_c, p0nm)
%%
% This function is used to write the wave vector information (k_Fourier,
% theta_Fourier, phi_Fourier) in the parameter control file 'parameters.m'.
% 
% change:
%  1. 'freq' in function "parameters.m"; 
%  2. 
%%

freq = fluid_c * k_Fourier / (2*pi);
file_parameters = 'parameters.m';

%% Read in the function files

fid_parameters = fopen(file_parameters, 'r');

ii = 0;
while feof(fid_parameters) ~= 1
    ii = ii + 1;
    line = fgetl(fid_parameters);
    if isempty(strfind(line, 'irregular_body = ')) ~= 1 && (isempty(strfind(line, '%irregular_body = ')) == 1 && isempty(strfind(line, '% irregular_body = ')) == 1)
        line_num_irreg = ii;
    end
    if isempty(strfind(line, 'freq = ')) ~= 1 && (isempty(strfind(line, '%freq = ')) == 1 && isempty(strfind(line, '% freq = ')) == 1)
        line_num_freq = ii;
    end
    if isempty(strfind(line, 'p_inlet = ')) ~= 1 && (isempty(strfind(line, '%p_inlet = ')) == 1 && isempty(strfind(line, '% p_inlet = ')) == 1)
        line_num_p_inlet = ii;
    end
    if isempty(strfind(line, 'theta_Fourier = ')) ~= 1 && (isempty(strfind(line, '%theta_Fourier = ')) == 1 && isempty(strfind(line, '% theta_Fourier = ')) == 1)
        line_num_theta_Fourier = ii;
    end
    if isempty(strfind(line, 'phi_Fourier = ')) ~= 1 && (isempty(strfind(line, '%phi_Fourier = ')) == 1 && isempty(strfind(line, '% phi_Fourier = ')) == 1)
        line_num_phi_Fourier = ii;
    end
    if isempty(strfind(line, 'trans_radius = ')) ~= 1 && (isempty(strfind(line, '%trans_radius = ')) == 1 && isempty(strfind(line, '% trans_radius = ')) == 1)
        line_num_trans_radius = ii;
    end
    if isempty(strfind(line, 'transducer_number = ')) ~= 1 && (isempty(strfind(line, '%transducer_number = ')) == 1 && isempty(strfind(line, '% transducer_number = ')) == 1)
        line_num_transducer_number = ii;
    end
end

%% Change specific parameters for the files

fid_parameters = fopen(file_parameters, 'r');

ii = 0;
content_parameters = {};
while feof(fid_parameters) ~= 1 
    ii = ii + 1;
    line = fgetl(fid_parameters);
    if ii == line_num_irreg
        line = ['irregular_body = ', num2str(irregular_body), ';                 % 1: the particle is axisymmetric arbitrary nonspherical object; 0: the particle is spherical object.'];
    end
    if ii == line_num_freq
        line = ['freq = ', num2str(freq), ';               % 1 MHz'];
    end
    if ii == line_num_p_inlet
        line = ['p_inlet = ', num2str(p0nm), ';               % initial amplitude of proba particle [Pa] if no attenuation effects'];
    end
    if ii == line_num_theta_Fourier
        line = ['theta_Fourier = ', num2str(theta_Fourier), ';'];
    end
    if ii == line_num_phi_Fourier
        line = ['phi_Fourier = ', num2str(phi_Fourier), ';'];
    end
    if ii == line_num_trans_radius
        line = ['trans_radius = ', num2str(trans_radius), ';       % NOTE: keep same in "total_effects_all_plain_wave_component.m". radius of transducer is 5mm'];
    end
    if ii == line_num_transducer_number
        line = ['transducer_number = ', num2str(transducer_number), ';          % NOTE: keep same in "total_effects_all_plain_wave_component.m". for phase array'];
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






%% ============Standard model of file "parameters.m"============== 
% %%
% % save the parameters, which uses to control the properties of the system.
% %%
% 
% %%
% % ================ fundimental parameters ==================
% 
% particle_radius = 2000*10^-6;   % particle radius (a); 30 um;
% b = 20 * particle_radius;       % integral radial distance (b) for beam-shape coeffcient, EMPIRICAL: MUST follow  b > a BUT not b >> a
% R = 1000 * particle_radius;     % farfield integral distance (R) MUST follow  R >> a
% 
% freq = 44901.6703;               % 1 MHz
% omega = freq*2*pi;
% 
% %sound_intensity = 33*10^3;      % w/m2
% 
% % ==========================================================
% %%
% 
% %%
% % ============ Database&Grid control parameters ============
% 
% multi_particle = 0;          % 1: multi-particle system; 0: single-particle system "pressure_contour.m" and "radiation_force_based_Analyses.m"
% multi_particle_force_torque = 1*multi_particle;     % function valid only if 'multi_particle == 1'
%                              % 2: obtain force&torque of all particle of multi-particle system; 
%                              %      Use for "multi_radiation_force_based_Analyses.m"
%                              % 1: obtain force&torque of singe probe particle of multi-particle system; 
%                              %      Use for "radiation_force_based_Analyses.m" and "radiation_torque_based_Analyses.m" 
%                              % 0: obtain pressure field of multi-particle system (if multi_particle == 1).
%                              %      Use for "pressure_contour.m" and "pressure_contour_vertical.m" 
% 
% vector_amination_visualization = 0;                 % function valid for both 'multi_particle == 0' and 'multi_particle == 1'
%                              % 1: for visualizing the force vectors or making animation for particle's moving; 
%                              %      Use for "visual_radiation_force_vector.m" and "animation.m" 
%                              % 0: for visualizing the pressure contour or getting single position's forces and torques. 
%                              %      Use for "pressure_contour.m", "pressure_contour_vertical.m", "radiation_force_based_Analyses.m" and "radiation_torque_based_Analyses.m" 
% 
% range_r_coeff = 16;          % the range enlarge coefficient of radius
% grid_resolution = 100;       % the grid notes along r-axis and theta-axis
% 
% if multi_particle == 0       % database maximum size (larger than Truncation requirement) for outsidest root
%     db_size_nn = 30;         % for single-particle system, maximum size of database larger than 30 will waste the computational resource. 
% else
%     db_size_nn = 25;         % for multi-particle system, maximum size of database larger than 20 will out of MATLAB computational memory.
% end
% 
% % ==========================================================
% %%
% 
% %% 
% % =================== particle properties ==================
% 
% BC = 'rigid';               % Boundary Conditions "rigid" or "compressible".
% %BC = 'olive_oil';
% %BC = 'benzene';
% %BC = 'polyurethane';
% 
% [particle_rho, particle_k] = material_property(BC, freq);
% 
% % ==========================================================
% %%
% 
% %%
% % =================== medium properties ====================
% 
% fluid = 'air';                  % medium types "air" or "water" and others
% %fluid = 'mixture';              % "mixture" means considering attenuation  
% %fluid = 'water';
% 
% %======= define the Mie scattering crystals' moving velocity ======
% if multi_particle == 1
%     velocity = [0, 0, 0] * particle_radius/1;    % particle radius per second
% end
% %==================================================================
% 
% % for "mixture (attenuation)" fluid only, the ratio of particles' volume to the whole
% % calculating domain volume
% particles_volume_fraction = 0.064;   
% % (Details of parameters setting refers "absorption_manual.txt")
% % =========== AIR, rigid, 100um, 1MHz (kr=1.85) ============
% % particles_volume_fraction = 0.065  -> particles_occupy_per_meter = 50.0%, Dis = 4.0r cr N = 11 
% % particles_volume_fraction = 0.049  -> particles_occupy_per_meter = 45.5%, Dis = 4.4r    N = 12  
% % particles_volume_fraction = 0.043  -> particles_occupy_per_meter = 43.5%, Dis = 4.6r cr  
% % particles_volume_fraction = 0.034  -> particles_occupy_per_meter = 40.0%, Dis = 5.0r    N = 13
% % particles_volume_fraction = 0.025  -> particles_occupy_per_meter = 36.4%, Dis = 5.5r    N = 14  
% % particles_volume_fraction = 0.023  -> particles_occupy_per_meter = 35.1%, Dis = 5.7r cr  
% % particles_volume_fraction = 0.019  -> particles_occupy_per_meter = 33.3%, Dis = 6.0r    N = 15 
% % particles_volume_fraction = 0.012  -> particles_occupy_per_meter = 28.6%, Dis = 7.0r    N = 18 
% % particles_volume_fraction = 0.008  -> particles_occupy_per_meter = 25.0%, Dis = 8.0r cr N = 20 
% % =========== AIR, rigid, 100um, 1MHz (kr=1.85) ============   
% 
% [fluid_c, fluid_k, fluid_rho, fluid_viscosity] = fluid_property(fluid, omega);
% 
% f_ka = real(fluid_k)*particle_radius;
% 
% % ==========================================================
% %%
% 
% %%
% % ===== compressible particle BC for scattering coeffs======
% 
% gama = (fluid_rho * particle_k) / (particle_rho * fluid_k);
% 
% 
% % ==========================================================
% %%
% 
% %%
% % ========= incident progress Directions&Deviation =========
% 
% % wave_type = 'plain';                   % 'plain', 'zero-Bessel' or 'non-zero-Bessel' and others  a
% % wave_type = 'standing_plain';
% % wave_type = 'zero-Bessel';
% % wave_type = 'non-zero-Bessel';
% wave_type = 'single_transducer';       % actually, it is revision of zero Bessel beam.
% % wave_type = 'phase_array_transducer';
% 
% if strcmp(wave_type, 'standing_plain') == 1
%     %direction = 'X';
%     %direction = 'Y';
%     direction = 'Z';
%     %direction = 'XY'; (functionless)
%     %direction = 'XZ';
%     
%     %================= attenuation adjustion ===============
%     % as determination of attenuation coefficient of mixture
%     % medium is based on single direction plain 
%     % (Wei_ronghao@1965@). therefore, for standing plain wave
%     % or multi-plaine waves, we need to turn back to single
%     % plain wave for attenuation coefficient "absp_factor". 
%     % once "absp_factor ~=0" means the coefficient already be
%     % determined, then below code will not run.
%     global absp_factor          % "absorption_factor.m"
%     if strcmp(fluid, 'mixture') == 1 && absp_factor == 0
%         wave_type = 'plain';    % temporally change "standing_plain" as "plain"
%         direction = 'Z';
%     end
%     %================= attenuation adjustion ===============
%     
%     symbol = 'positive';                % DON'T CHANGE! for two-direction (X&Z-axis) standing wave, we define symbol always "positive", then saving file looks order-well
% else
%     %direction = 'X';
%     %direction = 'Y';
%     direction = 'Z';
%     
%     symbol = 'positive';
%     %symbol = 'negative';
% end
% 
% if vector_amination_visualization == 1 || multi_particle_force_torque == 2
%     global X_vec Y_vec Z_vec        % for "visual_radiation_force_vector.m", "animation.m" and "multi_radiation_force_based_Analyses.m"
%     deviationX = X_vec;             % particle deviation from wave input positions (if "XZ" case, specific for "X")
%     deviationY = Y_vec;             % particle deviation from wave input positions (if "XY" case, specific for "Y")??
%     deviationZ = Z_vec;             % particle deviation from wave input positions (if "XZ" case, specific for "Z")
% else
%     deviationX = 0*particle_radius;             % particle deviation from wave input positions (if "XZ" case, specific for "X")
%     deviationY = 0*particle_radius;             % particle deviation from wave input positions (if "XY" case, specific for "Y")
%     deviationZ = 0*particle_radius;             % particle deviation from wave input positions (if "XZ" case, specific for "Z")
% end
% 
% [dir_sign, x_translation, y_translation, z_translation] = ...
%     wave_direction_particle_deviation(wave_type, direction, symbol, deviationX, deviationY, deviationZ);
% 
% % ==========================================================
% %%
% 
% %%
% % ==== incident progress wave Amplitude&Phase of Micros ====
% 
% %p_inlet = round(sqrt(2 * sound_intensity * fluid_c * fluid_rho));      % initial amplitude given based on input intensity
% p_inlet = 2000;               % initial amplitude of proba particle [Pa] if no attenuation effects
% phase_inlet = 0;            % initial phase of proba particle [rad]
% 
% % adjusting amplitude of proba particle [Pa] to keep microphones' amplitude
% % consistent no matter exist attenuation or not.
% [p_0] = proba_particle_input_amplitude(p_inlet, fluid, fluid_k, particle_radius, range_r_coeff);
% 
% % ==========================================================
% %%
% 
% %%
% % =========== incident progress wave Type input ============
% % The 'p_i' will be used in "beam_shape_coeff.m" for the numerical
% % beam-shape coefficient. However, for the transducer scenarios, as the
% % beam-shape coefficient of the wave functions can not be numerical
% % decomposite independently with the position vector vec{r}. So, for the
% % wave functions other than 'plain', 'zero-Bessel' and 'non-zero-Bessel'
% % types, we have to use Fourier expansion to simplify the wave functions as
% % a series of plain waves, and each of them can be express by
% % position-independent beam-shape coefficient. 
% % Therefore, for complex wave function, the below codes are useless, while
% % a new counterpart function "G_i.m" is used to obtain the transerve shape
% % function for further Fourier expansion.
% 
% BN = 0;
% BETA = 0;
% X0 = 0;
% Y0 = 0;
% Z0 = 0;
% trans_radius = 5 * 10^-3;       % radius of transducer is 5mm
% transducer_number = 3;          % for phase array
% theta_Fourier = 0.47162;
% phi_Fourier = 3.1416;
% [p_i] = ...
%     wave_function(p_0, wave_type, direction, dir_sign, ...
%     fluid_k, x_translation, y_translation, z_translation, trans_radius);
% 
% if (strcmp(wave_type, 'plain') == 1 || strcmp(wave_type, 'standing_plain') == 1) ...
%         && BETA ~= 0 && X0 ~= 0 && Y0 ~= 0 && Z0 ~= 0
%     error('No wave deviations for plain wave!\n');
% end
% if strcmp(wave_type, 'zero-Bessel') == 1 && BN ~= 0
%     error('Bn == 0 for zero Bessel beam!\n');
% end
% if strcmp(wave_type, 'non-zero-Bessel') == 1 && BN < 1
%     error('Bn >= 1 for non-zero Bessel beam!\n');
% end
% 
% % ==========================================================
% %%
%% ============Standard model of file "parameters.m"============== 