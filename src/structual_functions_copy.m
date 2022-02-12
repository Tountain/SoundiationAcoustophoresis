function [Gam, Gam_partial_u, Amd, Amd_partial_u, Rau, Rau_partial_u] = ...
                structual_functions(Cn, n, m)
%%
% this function is used to numerically evaluate the first structual
% function 'Gam(nn+1, mm+nn+1)' and its radial derivation
% 'Gam_partial_u(nn+1, mm+nn+1)'; the second structual function 'Amd(nn+1,
% mm+nn+1)' and its radial derivation 'Amd_partial_u(nn+1, mm+nn+1)'; and
% the third structual function 'Rau(nn+1, mm+nn+1)' and its radial
% derivation 'Rau_partial_u(nn+1, mm+nn+1)'.
%
% E.g., (refer to my note on 2021.4.28, point 8)
%     Cn=[0.005, 0]; (spherical)
%     Cn=[0.005, 0, 0.002]; (spheroidal)
%     Cn=[0.005, 0, 0, 0.001]; (trangular-cone)
%     Cn=[0.005, 0, 0, 0, 0.0005]; (diamond)
%%

parameters;
% freq = 40000;
% fluid_c = 340;
% particle_c = 1500;
% fluid_k = 2*pi*freq/fluid_c;
% particle_k = 2*pi*freq/particle_c;
% fluid_rho = 1.224;              % density of air
% particle_rho = 2000;               % density of particle

%% preparation of mapping functions 

coeff = sqrt(((2*n+1)/(4*pi)) * (factorial(n-m)/factorial(n+m)));
theta_w = linspace(0, pi, 1000);

[g, g_partial_u, f, f_partial_u] = mapping_functions(Cn, theta_w);
r0 = sqrt(g.^2 + f.^2);
theta0 = acos(g ./ r0);
kr0 = fluid_k * r0;
krp = particle_k * r0;

%% evaluation of structual functions

% abbreviations
r_partial_u0 = (g .* g_partial_u + f .* f_partial_u) ./ r0;
sin_theta0 = (g_partial_u .* r0 - r_partial_u0 .* g) ./ (r0.^2);
ortho_func = sLegendre(n, m, theta_w) .* sin(theta_w);
% ortho_func = sLegendre(n, m, theta0) .* sin(theta0);
% ortho_func = sLegendre(n, m, theta_w);


% the first structual function and its radial derivation
Gam_elements = fluid_rho * coeff * sBessel(n, kr0, 1) .* sLegendre(n, m, theta0) .* ortho_func;
Gam = trapz(theta_w, Gam_elements);

Gam_partial_u_elements = coeff * (fluid_k * r_partial_u0 .* sBessel_partial(n, kr0, 1) .* sLegendre(n, m, theta0) ...
                + sBessel(n, kr0, 1) .* sLegendre_partial(n, m, theta0) .* sin_theta0) .* ortho_func ./ sqrt(g_partial_u.^2 + f_partial_u.^2);
Gam_partial_u = trapz(theta_w, Gam_partial_u_elements);


% the second structual function and its radial derivation
Amd_elements = fluid_rho * coeff * sHankel(n, kr0, 1) .* sLegendre(n, m, theta0) .* ortho_func;
Amd = trapz(theta_w, Amd_elements);

Amd_partial_u_elements = coeff * (fluid_k * r_partial_u0 .* sHankel_partial(n, kr0, 1) .* sLegendre(n, m, theta0) ...
                + sHankel(n, kr0, 1) .* sLegendre_partial(n, m, theta0) .* sin_theta0) .* ortho_func ./ sqrt(g_partial_u.^2 + f_partial_u.^2);
Amd_partial_u = trapz(theta_w, Amd_partial_u_elements);


% the third structual function and its radial derivation

% for rigid B.C., as the 'particle_k' is set to zero, making the 'sBessel'
% and 'sBessel_partial' getting 'NaN' results. In order to avoid these
% errors, we intentionally set the 'Rau' and 'Rau_partial_u'. 
if strcmp(BC, 'rigid') == 1
    Rau = 1000;                    % any values except to zero
    Rau_partial_u = 0;
else
    Rau_elements = particle_rho * coeff * sBessel(n, krp, 1) .* sLegendre(n, m, theta0) .* ortho_func;
    Rau = trapz(theta_w, Rau_elements);

    Rau_partial_u_elements = coeff * (particle_k * r_partial_u0 .* sBessel_partial(n, krp, 1) .* sLegendre(n, m, theta0) ...
                    + sBessel(n, krp, 1) .* sLegendre_partial(n, m, theta0) .* sin_theta0) .* ortho_func ./ sqrt(g_partial_u.^2 + f_partial_u.^2);
    Rau_partial_u = trapz(theta_w, Rau_partial_u_elements);
end

%%