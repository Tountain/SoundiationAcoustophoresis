function [particle_rho, particle_k] = material_property(BC, freq)
%%
% this function is used to define the material (particles) property so that
% decide the boundary conditions and the scattering coefficient in
% "scattering_coefficient.m".
%
% sub-function of "parameters.m"
%%

omega = freq*2*pi;

if strcmp(BC, 'rigid') == 1
    particle_k = 0;
    particle_rho = 1500;    % water as example
    %gama = 0;
end
if strcmp(BC, 'olive_oil') == 1      % take "olive oil" as example [REF. G.T.S@2014(JASA)]
    particle_rho = 915.8;                               % density; kg/m3
    particle_c = 1464;                                  % sound speed of inside particle; m/s
    absorption = 4.10*10^(-14);                         % T.L.Szabo,?Time domain wave equations for lossy media obeying a frequency power law,? JASA,96,491?500(1994)
    miu = 2;
    ultra_absorp = absorption*(freq)^miu;               % Eq.(8) at G.T.S@2014@JASA
    particle_k = omega/particle_c + ultra_absorp*1i;     % inner wave number of particle
    %gama = (fluid_rho * particle_k) / (particle_rho * fluid_k);
end
if strcmp(BC, 'benzene') == 1        % take "benzene" as example  [REF. G.T.S@2014(JASA)]
    particle_rho = 870;                                 % density; kg/m3
    particle_c = 1295;                                  % sound speed of inside particle; m/s
    absorption = 2.21*10^(-14);                         % T.L.Szabo,?Time domain wave equations for lossy media obeying a frequency power law,? JASA,96,491?500(1994)
    miu = 2;
    ultra_absorp = absorption*(freq)^miu;               % Eq.(8) at G.T.S@2014@JASA
    particle_k = omega/particle_c + ultra_absorp*1i;     % inner wave number of particle
    %gama = (fluid_rho * particle_k) / (particle_rho * fluid_k);
end
if strcmp(BC, 'polyurethane') == 1   % take "polyurethane" as example  [REF. G.T.S@2012(IEEE)]
    particle_rho = 915.8;                               % density; kg/m3
    particle_c = 1468;                                  % sound speed of inside particle; m/s
    absorption = 0;                        
    miu = 0;
    ultra_absorp = 1.49;                                % Np*MHz/m
    particle_k = omega/particle_c + ultra_absorp*1i;     % inner wave number of particle
    %gama = (fluid_rho * particle_k) / (particle_rho * fluid_k);
end

%%