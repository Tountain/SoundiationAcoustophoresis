%%
% save the parameters, which uses to control the properties of the system.
%%

%%
% ================ fundimental parameters ==================

irregular_body = 1;                 % 1: the particle is axisymmetric arbitrary nonspherical object; 0: the particle is spherical object.

if irregular_body == 0
    particle_radius = 2*10^-3;   % particle radius (a); 30 um;
    Cn = [particle_radius, 0];
    theta_x = 0;
    theta_y = 0;
    theta_z = 0;
    theta_rotation = [theta_x, theta_y, theta_z];
elseif irregular_body == 1
    particle_radius = 0.002;                   % particle radial mean-radius (a); 30 um;
    %particle_major_length = 2*particle_radius * 2;    % the major length (along axisymmetric axis) of the irregular particle
    Cn = [0.002           0           0     0.00025];
    theta_x = 3.8452e-011;                % particle counter-clock wise rotation (position Tx) along x-axis for positive 'theta_x'
    theta_y = 0.79747;                % particle counter-clock wise rotation (position Ty) along y-axis for positive 'theta_y'
    theta_z = -2.5943e-012;                % particle counter-clock wise rotation (position Tz) along z-axis for positive 'theta_z'
    theta_rotation = [theta_x, theta_y, theta_z];
    % 'theta_rotation' follow right-hand principle to determine the
    % particle rotating direction.
end
b = 0.005;           % integral radial distance (b) for beam-shape coeffcient, EMPIRICAL: MUST follow  b > a and b < inter_dist BUT not b >> a
% R = 1000 * particle_radius;         % farfield integral distance (R) MUST follow  R >> a

freq = 40000;               % 1 MHz
omega = freq*2*pi;

%sound_intensity = 33*10^3;         % w/m2

% transducers parameters used in calculation
v0 = 1.5;                                                                     % transducer vibration velocity amplitude, m/s, assuming v0=1 m/s
inter_dist = 0.02;                                                         % distance between the transducer and the particle center (i.e., the origin of the coordinate system)

% ==========================================================
%%

%%
% ============ Database&Grid control parameters ============

multi_particle = 0;          % 1: multi-particle system; 0: single-particle system "pressure_contour.m" and "radiation_force_based_Analyses.m"
multi_particle_force_torque = 0*multi_particle;     % function valid only if 'multi_particle == 1'
                             % 2: obtain force&torque of all particle of multi-particle system; 
                             %      Use for "multi_radiation_force_based_Analyses.m"
                             % 1: obtain force&torque of singe probe particle of multi-particle system; 
                             %      Use for "radiation_force_based_Analyses.m" and "radiation_torque_based_Analyses.m" 
                             % 0: obtain pressure field of multi-particle system (if multi_particle == 1).
                             %      Use for "pressure_contour.m" and "pressure_contour_vertical.m" 

vector_amination_visualization = 0;                 % function valid for both 'multi_particle == 0' and 'multi_particle == 1'
                             % 1: for visualizing the force vectors or making animation for particle's moving; 
                             %      Use for "visual_radiation_force_vector.m" and "animation.m" 
                             % 0: for visualizing the pressure contour or getting single position's forces and torques. 
                             %      Use for "pressure_contour.m", "pressure_contour_vertical.m", "radiation_force_based_Analyses.m" and "radiation_torque_based_Analyses.m" 

range_r_coeff = 10;          % the range enlarge coefficient of radius
grid_resolution = 100;       % the grid notes along r-axis and theta-axis

if multi_particle == 0       % database maximum size (larger than Truncation requirement) for outsidest root
    db_size_nn = 8;         % for single-particle system, maximum size of database larger than 30 will waste the computational resource. 
else
    db_size_nn = 25;         % for multi-particle system, maximum size of database larger than 20 will out of MATLAB computational memory.
end

% ==========================================================
%%

%% 
% =================== particle properties ==================

BC = 'rigid';               % Boundary Conditions "rigid" or "soft".

% [particle_rho, particle_k] = material_property(BC, freq);
particle_rho = 15;
particle_k = 0;
particle_mass = 4/3 * pi * particle_radius^3 * particle_rho;
motion_inertia = (2/5) * particle_mass * particle_radius^3;

% ==========================================================
%%

%%
% =================== medium properties ====================

fluid = 'air';               % medium types "air" or "water" and others
%fluid= 'mixture';              % "mixture" means considering attenuation  

%======= define the Mie scattering crystals' moving velocity ======
if multi_particle == 1
    velocity = [0, 0, 0] * particle_radius/1;    % particle radius per second
end
%==================================================================

% for "mixture (attenuation)" fluid only, the ratio of particles' volume to the whole
% calculating domain volume
particles_volume_fraction = 0.064;   
% (Details of parameters setting refers "absorption_manual.txt")
% =========== AIR, rigid, 100um, 1MHz (kr=1.85) ============
% particles_volume_fraction = 0.065  -> particles_occupy_per_meter = 50.0%, Dis = 4.0r cr N = 11 
% particles_volume_fraction = 0.049  -> particles_occupy_per_meter = 45.5%, Dis = 4.4r    N = 12  
% particles_volume_fraction = 0.043  -> particles_occupy_per_meter = 43.5%, Dis = 4.6r cr  
% particles_volume_fraction = 0.034  -> particles_occupy_per_meter = 40.0%, Dis = 5.0r    N = 13
% particles_volume_fraction = 0.025  -> particles_occupy_per_meter = 36.4%, Dis = 5.5r    N = 14  
% particles_volume_fraction = 0.023  -> particles_occupy_per_meter = 35.1%, Dis = 5.7r cr  
% particles_volume_fraction = 0.019  -> particles_occupy_per_meter = 33.3%, Dis = 6.0r    N = 15 
% particles_volume_fraction = 0.012  -> particles_occupy_per_meter = 28.6%, Dis = 7.0r    N = 18 
% particles_volume_fraction = 0.008  -> particles_occupy_per_meter = 25.0%, Dis = 8.0r cr N = 20 
% =========== AIR, rigid, 100um, 1MHz (kr=1.85) ============   

[fluid_c, fluid_k, fluid_rho, fluid_viscosity] = fluid_property(fluid, omega);

f_ka = real(fluid_k)*particle_radius;

% ==========================================================
%%

%%
% ===== compressible particle BC for scattering coeffs======

% gama = (fluid_rho * particle_k) / (particle_rho * fluid_k);


% ==========================================================
%%

%%
% ========= incident progress Directions&Deviation =========

wave_type = 'phase_array_transducer';                   % 'plain', 'zero-Bessel' or 'non-zero-Bessel' and others
% wave_type= 'standing_plain';
% wave_type= 'zero-Bessel';
% wave_type= 'non-zero-Bessel';
% wave_type= 'single_transducer';                      % not used, since it can be replaced by setting 'transducer_number=1' and wave_type = 'phase_array_transducer'
% wave_type= 'phase_array_transducer';                 % based on addition theorem
% % wave_type= 'phase_array_transducer2';                  % not based on addition theorem, call 'array_wave_function.m'
% wave_type= 'single_transducer_standing';             % not used
% wave_type= 'phase_array_transducer_standing';        % not used

if strcmp(wave_type, 'standing_plain') == 1
    %direction = 'X';
    %direction = 'Y';
    direction = 'Z';
    %direction = 'XY'; (functionless)
    %direction = 'XZ';
    
%     symbol = 'positive';                % DON'T CHANGE! for two-direction (X&Z-axis) standing wave, we define symbol always "positive", then saving file looks order-well
    symbol = 'arbitrary';               % for arbitrary incidence of standing wave at potential well ONLY.
    if strcmp(symbol, 'arbitrary') == 1 % 'direction' must be 'Z'
        direction = 'Z';                % (theta_inc, phi_inc) = (0, 0) for along Z-axis
        theta_inc = 0+eps;                  % (theta_inc, phi_inc) = (pi/2, 0) for along X-axis
        phi_inc = 0+eps;                    % (theta_inc, phi_inc) = (pi/2, pi/2) for along Y-axis
        if theta_inc < 0                % theta_inc in [0, pi]
            theta_inc = theta_inc+pi;
        elseif theta_inc > pi
            theta_inc = theta_inc-pi;
        end
    end
    
    %================= attenuation adjustion ===============
    % as determination of attenuation coefficient of mixture
    % medium is based on single direction plain 
    % (Wei_ronghao@1965@). therefore, for standing plain wave
    % or multi-plaine waves, we need to turn back to single
    % plain wave for attenuation coefficient "absp_factor". 
    % once "absp_factor ~=0" means the coefficient already be
    % determined, then below code will not run.
    global absp_factor          % "absorption_factor.m"
    if strcmp(fluid, 'mixture') == 1 && absp_factor == 0
        wave_type = 'plain';    % temporally change "standing_plain" as "plain"
        direction = 'Z';
    end
    %================= attenuation adjustion ===============
    
elseif strcmp(wave_type, 'plain') == 1
    %direction = 'X';
    %direction = 'Y';
    direction = 'Z';
    
    symbol = 'positive';
    %symbol = 'negative';
    %symbol = 'arbitrary';               % for arbitrary incidence of plane wave mode ONLY.
    if strcmp(symbol, 'arbitrary') == 1 % 'direction' must be 'Z'
        direction = 'Z';                % (theta_inc, phi_inc) = (0, 0) for along Z-axis
        theta_inc = 0+eps;                  % (theta_inc, phi_inc) = (pi/2, 0) for along X-axis
        phi_inc = 0+eps;                    % (theta_inc, phi_inc) = (pi/2, pi/2) for along Y-axis
        if theta_inc < 0                % theta_inc in [0, pi]
            theta_inc = theta_inc+pi;
        elseif theta_inc > pi
            theta_inc = theta_inc-pi;
        end
    end
elseif strcmp(wave_type, 'phase_array_transducer') == 1 || strcmp(wave_type, 'phase_array_transducer2') == 1
    direction = 'Z';
    symbol = 'positive';
else  
    %direction = 'X';
    %direction = 'Y';
    direction = 'Z';
    
    symbol = 'positive';
    %symbol = 'negative';
end

deviationX = 0.0073383;
deviationY = -1.0163e-013;
deviationZ = 0.016263;

[dir_sign, x_translation, y_translation, z_translation] = ...
    wave_direction_particle_deviation(wave_type, direction, symbol, deviationX, deviationY, deviationZ);

% ==========================================================
%%

%%
% ==== incident progress wave Amplitude&Phase of Micros ====

p_inlet = 1;               % initial amplitude of proba particle [Pa] if no attenuation effects
phase_inlet = 0;            % initial phase of proba particle [rad]

% adjusting amplitude of proba particle [Pa] to keep microphones' amplitude
% consistent no matter exist attenuation or not.
[p_0] = proba_particle_input_amplitude(p_inlet, fluid, fluid_k, particle_radius, range_r_coeff);

% ==========================================================
%%

%%
% =========== incident progress wave Type input ============
% The 'p_i' will be used in "beam_shape_coeff.m" for the numerical
% beam-shape coefficient. However, for the transducer scenarios, as the
% beam-shape coefficient of the wave functions can not be numerical
% decomposite independently with the position vector vec{r}. So, for the
% wave functions other than 'plain', 'zero-Bessel' and 'non-zero-Bessel'
% types, we have to use Fourier expansion to simplify the wave functions as
% a series of plain waves, and each of them can be express by
% position-independent beam-shape coefficient. 
% Therefore, for complex wave function, the below codes are useless, while
% a new counterpart function "G_i.m" is used to obtain the transerve shape
% function for further Fourier expansion.

% Bessel beam parameters
BN = 0;
BETA = 0;
X0 = 0;
Y0 = 0;
Z0 = 0;

% transducer array parameters
trans_radius = 0.005;           % NOTE: keep same in "total_effects_all_plain_wave_component.m". radius of transducer is 5mm
transducer_number = 5;          % NOTE: keep same in "total_effects_all_plain_wave_component.m". for phase array
theta_Fourier = 0;         % used for plane wave expansion feature 
phi_Fourier = 0;           % used for plane wave expansion feature 
Pr_0 = - fluid_rho * fluid_c * fluid_k * trans_radius^2 * v0 * 1i / 2;      % the transducer output power for 'Dims_Fourier = 3'.  

if (strcmp(wave_type, 'phase_array_transducer') == 1 || strcmp(wave_type, 'phase_array_transducer2') == 1) && b > inter_dist / 1.1 
    error('Integral radial distance (b) for beam-shape coeffcient must not include the source terms!\n');
end

% wavefront function
[p_i] = ...
    wave_function(p_0, wave_type, direction, dir_sign, ...
    fluid_k, x_translation, y_translation, z_translation, trans_radius, transducer_number, Pr_0, inter_dist, theta_rotation);

if (strcmp(wave_type, 'plain') == 1 || strcmp(wave_type, 'standing_plain') == 1) ...
        && BETA ~= 0 && X0 ~= 0 && Y0 ~= 0 && Z0 ~= 0
    error('No wave deviations for plain wave!\n');
end
if strcmp(wave_type, 'zero-Bessel') == 1 && BN ~= 0
    error('Bn == 0 for zero Bessel beam!\n');
end
if strcmp(wave_type, 'non-zero-Bessel') == 1 && BN < 1
    error('Bn >= 1 for non-zero Bessel beam!\n');
end

% ==========================================================
%%
