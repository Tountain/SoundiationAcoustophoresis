function [fluid_c, fluid_k, fluid_rho, fluid_viscosity] = fluid_property(fluid, omega)
%%
% this function is used to define the fluid property, such as dentisy,
% wavenumber, attenuation or absorption, sound speed and so on. 
%
% sub-function of "parameters.m"
%%

if strcmp(fluid, 'air') == 1
    fluid_c = 340;              % sound speed; m/s
    fluid_k = omega / fluid_c;  % wave number
    fluid_rho = 1.224;          % density; kg/m3
    fluid_viscosity = 18.5 * 10^-6;        % Pa*s
end
if strcmp(fluid, 'mixture') == 1
    % this factor will be zero (zero is requirement for determine the value
    % "absp_factor") for the first time before running
    % "absorption_factor.m", and will be revised after running
    % "absorption_factor.m".   
    global absp_factor          % "absorption_factor.m"
    
    fluid_c = 340;              % sound speed; m/s
    %absp_factor = abs(0.1 * omega / fluid_c);   % From my derivation of attenuation in mixture
    %absp_factor = absorption_factor(omega, fluid_c, particle_radius, particle_rho);
    fluid_k = omega / fluid_c + absp_factor * 1i;  % wave number
    fluid_rho = 1.224;           % density; kg/m3
    fluid_viscosity = 18.5 * 10^-6;        % Pa*s
end
if strcmp(fluid, 'water') == 1
    fluid_c = 1480;             % sound speed; m/s
    fluid_k = omega / fluid_c;  % wave number
    %fluid_k = omega / sqrt(fluid_c^2 + 4/3 * omega * 10^-6 * 1i);
    fluid_rho = 1000;           % density; kg/m3
    fluid_viscosity = 1.01 * 10^-3;      % dynamic viscosity Pa*s
end


%%