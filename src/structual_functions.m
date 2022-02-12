function [Gam, Gam_partial_u, Amd, Amd_partial_u, Rau, Rau_partial_u] = ...
                structual_functions(n_s, m_s)
%%
% this function is used to numerically evaluate the first structual
% function 'Gam' and its radial derivation 'Gam_partial_u'; the second
% structual function 'Amd' and its radial derivation 'Amd_partial_u'; and 
% the third structual function 'Rau' and its radial derivation
% 'Rau_partial_u'. 
%
% 'n_s, m_s' are the orthogonal indice of the eigenfunctions degrees.
%
% The return variables should be column vectors, which save the
% corresponding structural functions for different index 'n'. (Note that
% index 'm=m_s' due to the orthogonality of azimuthal direction)
%
% Mapping coefficient 'Cn': (refer to my note on 2021.4.28, point 8)
%     Cn=[0.005, 0]; (spherical)
%     Cn=[0.005, 0, 0.002]; (spheroidal)
%     Cn=[0.005, 0, 0, 0.001]; (trangular-cone)
%     Cn=[0.005, 0, 0, 0, 0.0005]; (diamond)
%
% Refer to my note: Point 10, Conformal mappring for irregular body, page
% 11 and 12.
%%

parameters;
delta_u0 = 0.00001;
% freq = 40000;
% fluid_c = 340;
% particle_c = 1500;
% fluid_k = 2*pi*freq/fluid_c;
% particle_k = 2*pi*freq/particle_c;
% fluid_rho = 1.224;              % density of air
% particle_rho = 2000;               % density of particle

mm = m_s;        % orthogonality requirement along the azimuthal angular coordinate

%% preparation of mapping functions 

coeff = @(nn,mm) sqrt(((2*nn+1)/(4*pi)) * (factorial(nn-mm)/factorial(nn+mm)));
theta_w = linspace(0, pi, 1000);

[g, dev_g, g_partial_u, f, dev_f, f_partial_u] = mapping_functions(Cn, theta_w, delta_u0);
r0 = sqrt(g.^2 + f.^2);
theta0 = acos(g ./ r0);
kr0 = fluid_k * r0;
krp = particle_k * r0;

dev_r0 = sqrt(dev_g.^2 + dev_f.^2);
dev_theta0 = acos(dev_g ./ dev_r0);
dev_kr0 = fluid_k * dev_r0;
dev_krp = particle_k * dev_r0;

%% evaluation of structual functions

%% abbreviations
r_partial_u0 = (g .* g_partial_u + f .* f_partial_u) ./ r0;
sin_theta0 = (g_partial_u .* r0 - r_partial_u0 .* g) ./ (r0.^2);
% theta0_partial_u0 = (g_partial_u .* r0 - r_partial_u0 .* g) ./ (r0.^2);
% sin_theta0 = sin(theta0);   % sin_theta0 = sqrt((f.^2) ./ (r0.^2));
ortho_func = sLegendre(n_s, m_s, theta_w) .* sin(theta_w);
% ortho_func = sLegendre(nn, mm, theta0) .* sin(theta0);
% ortho_func = sLegendre(nn, mm, theta_w);


%% the first structual function and its radial derivation
% Gam_elements = fluid_rho * coeff * sBessel(n, kr0, 1) .* sLegendre(n, m, theta0) .* ortho_func;
% Gam = trapz(theta_w, Gam_elements);
if strcmp(BC, 'rigid') ~= 1
    Gam = zeros(db_size_nn+1, 1);
    fluid_rho = fluid_rho;
    for nn = 0 : db_size_nn
        if nn < abs(mm)
            Gam(nn+1, 1) = 0;
        else
            Gam_elements = feval(coeff, nn, mm) * fluid_rho * sBessel(nn, kr0, 1) .* sLegendre(nn, mm, theta0) .* ortho_func;
            Gam(nn+1, 1) = trapz(theta_w, Gam_elements);
        end
    end
else
    Gam = 0;
end

% Gam_partial_u_elements = coeff * (fluid_k * r_partial_u0 .* sBessel_partial(n, kr0, 1) .* sLegendre(n, m, theta0) ...
%                 + sBessel(n, kr0, 1) .* sLegendre_partial(n, m, theta0) .* sin_theta0) .* ortho_func ./ sqrt(g_partial_u.^2 + f_partial_u.^2);
% Gam_partial_u = trapz(theta_w, Gam_partial_u_elements);

% The summation index 'nn' is truncated by 'db_size_nn' (defined in
% "parameters.m").
% if strcmp(BC, 'rigid') == 1
%     Gam_partial_u = zeros(db_size_nn+1, 1);
%     fluid_k = fluid_k;
%     for nn = 0 : db_size_nn
%         if nn < abs(mm)
%             Gam_partial_u(nn+1, 1) = 0;
%         else
%             Gam_partial_u_elements = feval(coeff, nn, mm) * (fluid_k * r_partial_u0 .* sBessel_partial(nn, kr0, 1) .* sLegendre(nn, mm, theta0) ...
%                         + sBessel(nn, kr0, 1) .* sLegendre_partial(nn, mm, theta0) .* sin_theta0) .* ortho_func ./ sqrt(g_partial_u.^2 + f_partial_u.^2);
%     %         Gam_partial_u_elements = coeff(nn, mm) * (fluid_k * r_partial_u0 .* sBessel_partial(nn, kr0, 1) .* sLegendre(nn, mm, theta0) ...
%     %                     + sBessel(nn, kr0, 1) .* sLegendre_partial(nn, mm, theta0) .* (-sin_theta0) .* theta0_partial_u0) ... 
%     %                     .* ortho_func ./ sqrt(g_partial_u.^2 + f_partial_u.^2);
%             Gam_partial_u(nn+1, 1) = trapz(theta_w, Gam_partial_u_elements);
%         end
%     end
% else
%     Gam_partial_u = 0;
% end

if strcmp(BC, 'rigid') == 1
    Gam_partial_u = zeros(db_size_nn+1, 1);
    fluid_k = fluid_k;
    for nn = 0 : db_size_nn
        if nn < abs(mm)
            Gam_partial_u(nn+1, 1) = 0;
        else
            Gam_elements = feval(coeff, nn, mm) * fluid_rho * sBessel(nn, kr0, 1) .* sLegendre(nn, mm, theta0) .* ortho_func;
            dev_Gam_elements = feval(coeff, nn, mm) * fluid_rho * sBessel(nn, dev_kr0, 1) .* sLegendre(nn, mm, dev_theta0) .* ortho_func;
            Gam_partial_u_elements = (dev_Gam_elements - Gam_elements) / delta_u0;
            Gam_partial_u(nn+1, 1) = trapz(theta_w, Gam_partial_u_elements);
        end
    end
else
    Gam_partial_u = 0;
end


%% the second structual function and its radial derivation
% Amd_elements = fluid_rho * coeff * sHankel(n, kr0, 1) .* sLegendre(n, m, theta0) .* ortho_func;
% Amd = trapz(theta_w, Amd_elements);
if strcmp(BC, 'rigid') ~= 1
    Amd = zeros(db_size_nn+1, 1);
    for nn = 0 : db_size_nn
        if nn < abs(mm)
            Amd(nn+1, 1) = 0;
        else
            Amd_elements = feval(coeff, nn, mm) * fluid_rho * sHankel(nn, kr0, 1) .* sLegendre(nn, mm, theta0) .* ortho_func;
            Amd(nn+1, 1) = trapz(theta_w, Amd_elements);
        end
    end
else
    Amd = 0;
end

% Amd_partial_u_elements = coeff * (fluid_k * r_partial_u0 .* sHankel_partial(n, kr0, 1) .* sLegendre(n, m, theta0) ...
%                 + sHankel(n, kr0, 1) .* sLegendre_partial(n, m, theta0) .* sin_theta0) .* ortho_func ./ sqrt(g_partial_u.^2 + f_partial_u.^2);
% Amd_partial_u = trapz(theta_w, Amd_partial_u_elements);

% The summation index 'nn' is truncated by 'db_size_nn' (defined in
% "parameters.m").
% if strcmp(BC, 'rigid') == 1
%     Amd_partial_u = zeros(db_size_nn + 1, 1);
%     for nn = 0 : db_size_nn
%         if nn < abs(mm)
%             Amd_partial_u(nn+1, 1) = 0;
%         else
%             Amd_partial_u_elements = feval(coeff, nn, mm) * (fluid_k * r_partial_u0 .* sHankel_partial(nn, kr0, 1) .* sLegendre(nn, mm, theta0) ...
%                         + sHankel(nn, kr0, 1) .* sLegendre_partial(nn, mm, theta0) .* sin_theta0) .* ortho_func ./ sqrt(g_partial_u.^2 + f_partial_u.^2);
%     %         Amd_partial_u_elements = coeff(nn, mm) * (fluid_k * r_partial_u0 .* sHankel_partial(nn, kr0, 1) .* sLegendre(nn, mm, theta0) ...
%     %                     + sHankel(nn, kr0, 1) .* sLegendre_partial(nn, mm, theta0) .* (-sin_theta0) .* theta0_partial_u0) ...
%     %                     .* ortho_func ./ sqrt(g_partial_u.^2 + f_partial_u.^2);
%             Amd_partial_u(nn+1, 1) = trapz(theta_w, Amd_partial_u_elements);
%         end
%     end
% else
%     Amd_partial_u = 0;
% end

if strcmp(BC, 'rigid') == 1
    Amd_partial_u = zeros(db_size_nn+1, 1);
    for nn = 0 : db_size_nn
        if nn < abs(mm)
            Amd_partial_u(nn+1, 1) = 0;
        else
            Amd_elements = feval(coeff, nn, mm) * fluid_rho * sHankel(nn, kr0, 1) .* sLegendre(nn, mm, theta0) .* ortho_func;
            dev_Amd_elements = feval(coeff, nn, mm) * fluid_rho * sHankel(nn, dev_kr0, 1) .* sLegendre(nn, mm, dev_theta0) .* ortho_func;
            Amd_partial_u_elements = (dev_Amd_elements - Amd_elements) / delta_u0;
            Amd_partial_u(nn+1, 1) = trapz(theta_w, Amd_partial_u_elements);
        end
    end
else
    Amd_partial_u = 0;
end


%% the third structual function and its radial derivation

% for rigid B.C., as the 'particle_k' is set to zero, making the 'sBessel'
% and 'sBessel_partial' getting 'NaN' results. In order to avoid these
% errors, we intentionally set the 'Rau' and 'Rau_partial_u'. 
if strcmp(BC, 'rigid') == 1
    Rau = 1000;                    % any values except to zero
    Rau_partial_u = 0;
else
%     Rau_elements = particle_rho * coeff * sBessel(nn, krp, 1) .* sLegendre(nn, mm, theta0) .* ortho_func;
%     Rau = trapz(theta_w, Rau_elements);
% 
%     Rau_partial_u_elements = coeff * (particle_k * r_partial_u0 .* sBessel_partial(nn, krp, 1) .* sLegendre(nn, mm, theta0) ...
%                     + sBessel(nn, krp, 1) .* sLegendre_partial(nn, mm, theta0) .* sin_theta0) .* ortho_func ./ sqrt(g_partial_u.^2 + f_partial_u.^2);
%     Rau_partial_u = trapz(theta_w, Rau_partial_u_elements);
    Rau = 1000;                    % any values except to zero
    Rau_partial_u = 0;
end

%%




% function [Gam, Gam_partial_u, Amd, Amd_partial_u, Rau, Rau_partial_u] = ...
%                 structual_functions(n_s, m_s)
% %%
% % this function is used to numerically evaluate the first structual
% % function 'Gam' and its radial derivation 'Gam_partial_u'; the second
% % structual function 'Amd' and its radial derivation 'Amd_partial_u'; and 
% % the third structual function 'Rau' and its radial derivation
% % 'Rau_partial_u'. 
% %
% % 'n_s, m_s' are the orthogonal indice of the eigenfunctions degrees.
% %
% % The return variables should be column vectors, which save the
% % corresponding structural functions for different index 'n'. (Note that
% % index 'm=m_s' due to the orthogonality of azimuthal direction)
% %
% % Mapping coefficient 'Cn': (refer to my note on 2021.4.28, point 8)
% %     Cn=[0.005, 0]; (spherical)
% %     Cn=[0.005, 0, 0.002]; (spheroidal)
% %     Cn=[0.005, 0, 0, 0.001]; (trangular-cone)
% %     Cn=[0.005, 0, 0, 0, 0.0005]; (diamond)
% %
% % Refer to my note: Point 10, Conformal mappring for irregular body, page
% % 11 and 12.
% %%
% 
% parameters;
% % freq = 40000;
% % fluid_c = 340;
% % particle_c = 1500;
% % fluid_k = 2*pi*freq/fluid_c;
% % particle_k = 2*pi*freq/particle_c;
% % fluid_rho = 1.224;              % density of air
% % particle_rho = 2000;               % density of particle
% 
% mm = m_s;        % orthogonality requirement along the azimuthal angular coordinate
% 
% %% preparation of mapping functions 
% 
% coeff = @(nn,mm) sqrt(((2*nn+1)/(4*pi)) * (factorial(nn-mm)/factorial(nn+mm)));
% theta_w = linspace(0, pi, 1000);
% 
% [g, g_partial_u, f, f_partial_u] = mapping_functions(Cn, theta_w);
% r0 = sqrt(g.^2 + f.^2);
% theta0 = acos(g ./ r0);
% kr0 = fluid_k * r0;
% krp = particle_k * r0;
% 
% 
% %% evaluation of structual functions
% 
% %% abbreviations
% r_partial_u0 = (g .* g_partial_u + f .* f_partial_u) ./ r0;
% sin_theta0 = (g_partial_u .* r0 - r_partial_u0 .* g) ./ (r0.^2);
% % theta0_partial_u0 = (g_partial_u .* r0 - r_partial_u0 .* g) ./ (r0.^2);
% % sin_theta0 = sin(theta0);   % sin_theta0 = sqrt((f.^2) ./ (r0.^2));
% ortho_func = sLegendre(n_s, m_s, theta_w) .* sin(theta_w);
% % ortho_func = sLegendre(nn, mm, theta0) .* sin(theta0);
% % ortho_func = sLegendre(nn, mm, theta_w);
% 
% 
% %% the first structual function and its radial derivation
% % Gam_elements = fluid_rho * coeff * sBessel(n, kr0, 1) .* sLegendre(n, m, theta0) .* ortho_func;
% % Gam = trapz(theta_w, Gam_elements);
% Gam = zeros(db_size_nn+1, 1);
% for nn = 0 : db_size_nn
%     if nn < abs(mm)
%         Gam(nn+1, 1) = 0;
%     else
%         Gam_elements = coeff(nn, mm) * fluid_rho * sBessel(nn, kr0, 1) .* sLegendre(nn, mm, theta0) .* ortho_func;
%         Gam(nn+1, 1) = trapz(theta_w, Gam_elements);
%     end
% end
% 
% % Gam_partial_u_elements = coeff * (fluid_k * r_partial_u0 .* sBessel_partial(n, kr0, 1) .* sLegendre(n, m, theta0) ...
% %                 + sBessel(n, kr0, 1) .* sLegendre_partial(n, m, theta0) .* sin_theta0) .* ortho_func ./ sqrt(g_partial_u.^2 + f_partial_u.^2);
% % Gam_partial_u = trapz(theta_w, Gam_partial_u_elements);
% 
% % The summation index 'nn' is truncated by 'db_size_nn' (defined in
% % "parameters.m").
% Gam_partial_u = zeros(db_size_nn+1, 1);
% for nn = 0 : db_size_nn
%     if nn < abs(mm)
%         Gam_partial_u(nn+1, 1) = 0;
%     else
%         Gam_partial_u_elements = coeff(nn, mm) * (fluid_k * r_partial_u0 .* sBessel_partial(nn, kr0, 1) .* sLegendre(nn, mm, theta0) ...
%                     + sBessel(nn, kr0, 1) .* sLegendre_partial(nn, mm, theta0) .* sin_theta0) .* ortho_func ./ sqrt(g_partial_u.^2 + f_partial_u.^2);
% %         Gam_partial_u_elements = coeff(nn, mm) * (fluid_k * r_partial_u0 .* sBessel_partial(nn, kr0, 1) .* sLegendre(nn, mm, theta0) ...
% %                     + sBessel(nn, kr0, 1) .* sLegendre_partial(nn, mm, theta0) .* (-sin_theta0) .* theta0_partial_u0) ... 
% %                     .* ortho_func ./ sqrt(g_partial_u.^2 + f_partial_u.^2);
%         Gam_partial_u(nn+1, 1) = trapz(theta_w, Gam_partial_u_elements);
%     end
% end
% 
% 
% %% the second structual function and its radial derivation
% % Amd_elements = fluid_rho * coeff * sHankel(n, kr0, 1) .* sLegendre(n, m, theta0) .* ortho_func;
% % Amd = trapz(theta_w, Amd_elements);
% Amd = zeros(db_size_nn+1, 1);
% for nn = 0 : db_size_nn
%     if nn < abs(mm)
%         Amd(nn+1, 1) = 0;
%     else
%         Amd_elements = coeff(nn, mm) * fluid_rho * sHankel(nn, kr0, 1) .* sLegendre(nn, mm, theta0) .* ortho_func;
%         Amd(nn+1, 1) = trapz(theta_w, Amd_elements);
%     end
% end
% 
% % Amd_partial_u_elements = coeff * (fluid_k * r_partial_u0 .* sHankel_partial(n, kr0, 1) .* sLegendre(n, m, theta0) ...
% %                 + sHankel(n, kr0, 1) .* sLegendre_partial(n, m, theta0) .* sin_theta0) .* ortho_func ./ sqrt(g_partial_u.^2 + f_partial_u.^2);
% % Amd_partial_u = trapz(theta_w, Amd_partial_u_elements);
% 
% % The summation index 'nn' is truncated by 'db_size_nn' (defined in
% % "parameters.m").
% Amd_partial_u = zeros(db_size_nn + 1, 1);
% for nn = 0 : db_size_nn
%     if nn < abs(mm)
%         Amd_partial_u(nn+1, 1) = 0;
%     else
%         Amd_partial_u_elements = coeff(nn, mm) * (fluid_k * r_partial_u0 .* sHankel_partial(nn, kr0, 1) .* sLegendre(nn, mm, theta0) ...
%                     + sHankel(nn, kr0, 1) .* sLegendre_partial(nn, mm, theta0) .* sin_theta0) .* ortho_func ./ sqrt(g_partial_u.^2 + f_partial_u.^2);
% %         Amd_partial_u_elements = coeff(nn, mm) * (fluid_k * r_partial_u0 .* sHankel_partial(nn, kr0, 1) .* sLegendre(nn, mm, theta0) ...
% %                     + sHankel(nn, kr0, 1) .* sLegendre_partial(nn, mm, theta0) .* (-sin_theta0) .* theta0_partial_u0) ...
% %                     .* ortho_func ./ sqrt(g_partial_u.^2 + f_partial_u.^2);
%         Amd_partial_u(nn+1, 1) = trapz(theta_w, Amd_partial_u_elements);
%     end
% end
% 
% %% the third structual function and its radial derivation
% 
% % for rigid B.C., as the 'particle_k' is set to zero, making the 'sBessel'
% % and 'sBessel_partial' getting 'NaN' results. In order to avoid these
% % errors, we intentionally set the 'Rau' and 'Rau_partial_u'. 
% if strcmp(BC, 'rigid') == 1
%     Rau = 1000;                    % any values except to zero
%     Rau_partial_u = 0;
% else
% %     Rau_elements = particle_rho * coeff * sBessel(nn, krp, 1) .* sLegendre(nn, mm, theta0) .* ortho_func;
% %     Rau = trapz(theta_w, Rau_elements);
% % 
% %     Rau_partial_u_elements = coeff * (particle_k * r_partial_u0 .* sBessel_partial(nn, krp, 1) .* sLegendre(nn, mm, theta0) ...
% %                     + sBessel(nn, krp, 1) .* sLegendre_partial(nn, mm, theta0) .* sin_theta0) .* ortho_func ./ sqrt(g_partial_u.^2 + f_partial_u.^2);
% %     Rau_partial_u = trapz(theta_w, Rau_partial_u_elements);
%     Rau = 1000;                    % any values except to zero
%     Rau_partial_u = 0;
% end
% 
% %%