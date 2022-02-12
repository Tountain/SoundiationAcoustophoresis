function [Anm, anm] = beam_shape_coeff(n, m)
%%
% this function is used to determine the beam shape coefficient "Anm" from
% complex incident pressure amplitude "p_i".
% numerical integral based on "trapz".
% Input: 
%        "p_i" -- is a handle function type.
%       eg. p_i = @(r,theta,phi) exp(i * r * f(phi') * g(theta))   
%                (Here "phi'" MUST before "theta")
%                where, "phi" and "theta" are row vectors, so "p_i" will be
%                a matrix of size of "length(phi)*length(theta)".
%%

% Based on paper "Glauber_T.Silva@2011@Off-Axis Scattering_of an Ultrasound 
% Bessel Beam by a Sphere" Eq.(8)

% p_i = @(r,theta,phi) ones(length(phi),1)*exp(-i*fluid_k*r*cos(theta));   
% plain wave
% p_i = @(r,beta,theta,phi) ones(length(phi),1)*(exp(i*fluid_k*r*cos(theta)*cos(beta)).*besselj(0,r*sin(theta)*sin(beta)));
% zero-order Bessel Beam; WHERE "beta" is the half-cone angle.

% kb = 100;                        % dimensionaless radial distance, MUST follow kb=k*b >> k*a
parameters;


theta = linspace(0, pi, 100*2);   % Empirical Relation: the interval of 
phi = linspace(0, 2*pi, 100*4);   % "theta" and "phi" should larger than "kb".

conj_Ynm = conj(sHarmonics(n, m, theta, phi));
inner_int = zeros(length(phi), 1);

if strcmp(wave_type, 'plain') == 1
    p_b = p_i(b, theta, phi);
end
if strcmp(wave_type, 'standing_plain') == 1
    p_b = p_i(b, theta, phi);
end
if strcmp(wave_type, 'zero-Bessel') == 1
    p_b = p_i(b, BETA, X0, Y0, Z0, theta, phi);
end
if strcmp(wave_type, 'non-zero-Bessel') == 1
    p_b = p_i(b, BN, BETA, X0, Y0, Z0, theta, phi);
end
if strcmp(wave_type, 'single_transducer') == 1 || strcmp(wave_type, 'single_transducer_standing') == 1    % circular piston wave function
    p_b = p_i(b, theta, phi);
%     p_b = p_i(b, BETA, X0, Y0, Z0, theta, phi);
end
if strcmp(wave_type, 'phase_array_transducer') == 1 || strcmp(wave_type, 'phase_array_transducer2') == 1 || strcmp(wave_type, 'phase_array_transducer_standing') == 1     % circular phase array wave function
    p_b = p_i(b, theta, phi);
end

% =============================================================
% there is a theoratical module for (standing) plain wave to improve
% robustness of farfield results.
% seeing my derivation of Page 7-(a) ~ 7-(b)
% if strcmp(wave_type, 'plain') == 1
%     if strcmp(direction, 'X') == 1
%         if strcmp(symbol, 'positive') == 1
%             Anm = p_0*4*pi * (1i)^(n) * conj(sHarmonics(n,m,pi/2,0)) * exp(1i*fluid_k*x_translation);
%             anm = Anm / ((2*n + 1) * 1i^n); 
%         end
%         if strcmp(symbol, 'negative') == 1
%             Anm = p_0*4*pi * (-1i)^(n) * conj(sHarmonics(n,m,pi/2,0)) * exp(-1i*fluid_k*x_translation);
%             anm = Anm / ((2*n + 1) * 1i^n); 
%         end
%         return;
%     end
%     if strcmp(direction, 'Y') == 1
%         if strcmp(symbol, 'positive') == 1
%             Anm = p_0*4*pi * (1i)^(n) * conj(sHarmonics(n,m,pi/2,pi/2)) * exp(1i*fluid_k*y_translation);
%             anm = Anm / ((2*n + 1) * 1i^n); 
%         end
%         if strcmp(symbol, 'negative') == 1
%             Anm = p_0*4*pi * (-1i)^(n) * conj(sHarmonics(n,m,pi/2,pi/2)) * exp(-1i*fluid_k*y_translation);
%             anm = Anm / ((2*n + 1) * 1i^n); 
%         end
%         return;
%     end
%     if strcmp(direction, 'Z') == 1
%         if strcmp(symbol, 'positive') == 1
%             Anm = p_0*4*pi * (1i)^(n) * conj(sHarmonics(n,m,0,0)) * exp(1i*fluid_k*z_translation);
%             anm = Anm / ((2*n + 1) * 1i^n); 
%         end
%         if strcmp(symbol, 'negative') == 1
%             Anm = p_0*4*pi * (-1i)^(n) * conj(sHarmonics(n,m,0,0)) * exp(-1i*fluid_k*z_translation);
%             anm = Anm / ((2*n + 1) * 1i^n); 
%         end
%         if strcmp(symbol, 'arbitrary') == 1
%             Anm = p_0*4*pi * (1i)^(n) * conj(sHarmonics(n,m,theta_inc,phi_inc)) * ...
%                 exp(1i*fluid_k*(sin(theta_inc)*cos(phi_inc)*x_translation + ...
%                 sin(theta_inc)*sin(phi_inc)*y_translation + cos(theta_inc)*z_translation));
%             anm = Anm / ((2*n + 1) * 1i^n); 
%         end
%         return;
%     end
% end
if strcmp(wave_type, 'standing_plain') == 1
    if strcmp(direction, 'X') == 1
        Anm = 2*p_0*4*pi * cos(fluid_k*x_translation + n*pi/2) * conj(sHarmonics(n,m,pi/2,0));
        anm = Anm / ((2*n + 1) * 1i^n);
        return;
    end
    if strcmp(direction, 'Y') == 1
        Anm = 2*p_0*4*pi * cos(fluid_k*y_translation + n*pi/2) * conj(sHarmonics(n,m,pi/2,pi/2));
        anm = Anm / ((2*n + 1) * 1i^n);
        return;
    end 
    if strcmp(direction, 'Z') == 1
        if strcmp(symbol, 'arbitrary') == 1         % particle placed in the center (in the acoustic radiation-force potential well)
            Anm = 2*p_0*4*pi * cos(fluid_k*(2*pi/fluid_k/4) + n*pi/2) * conj(sHarmonics(n,m,theta_inc,phi_inc));
            anm = Anm / ((2*n + 1) * 1i^n);
        else
            Anm = 2*p_0*4*pi * cos(fluid_k*z_translation + n*pi/2) * conj(sHarmonics(n,m,0,0));
            anm = Anm / ((2*n + 1) * 1i^n);
        end
        return;
    end
    if strcmp(direction, 'XZ') == 1
        Anm = 2*p_0*4*pi * (cos(fluid_k*x_translation + n*pi/2) * conj(sHarmonics(n,m,pi/2,0)) + ...
            cos(fluid_k*z_translation + n*pi/2) * conj(sHarmonics(n,m,0,0)));
        anm = Anm / ((2*n + 1) * 1i^n);
        return;
    end
end

% ================= General Plain Wave for each wave component decomposite from Fourier Expansion ===================
% if strcmp(wave_type, 'single_transducer') == 1          % Ref. notebook. 2020.9.5. and Eq.(6-175) of E.G.William
%     if strcmp(direction, 'Z') == 1
%         if strcmp(symbol, 'positive') == 1
%             Anm = p_0*4*pi * (1i)^(n) * conj(sHarmonics(n,m,theta_Fourier,phi_Fourier)) * ...
%                 exp(1i*fluid_k*(sin(theta_Fourier)*cos(phi_Fourier)*x_translation + ...
%                 sin(theta_Fourier)*sin(phi_Fourier)*y_translation + cos(theta_Fourier)*z_translation));
%             anm = Anm / ((2*n + 1) * 1i^n); 
%         end
%         if strcmp(symbol, 'negative') == 1
%             Anm = p_0*4*pi * (-1i)^(n) * conj(sHarmonics(n,m,theta_Fourier,phi_Fourier)) * ...
%                 exp(-1i*fluid_k*(sin(theta_Fourier)*cos(phi_Fourier)*x_translation + ...
%                 sin(theta_Fourier)*sin(phi_Fourier)*y_translation + cos(theta_Fourier)*z_translation));
%             anm = Anm / ((2*n + 1) * 1i^n); 
%         end
%         return;
%     end
% end
% if strcmp(wave_type, 'single_transducer_standing') == 1   %!!!
%     if strcmp(direction, 'Z') == 1
%         if strcmp(symbol, 'positive') == 1
%             Anm = 2*p_0*4*pi * cos(fluid_k*(sin(theta_Fourier)*cos(phi_Fourier)*x_translation + ...
%                 sin(theta_Fourier)*sin(phi_Fourier)*y_translation + cos(theta_Fourier)*z_translation) + n*pi/2) * ...
%                 conj(sHarmonics(n,m,theta_Fourier,phi_Fourier));
%             anm = Anm / ((2*n + 1) * 1i^n); 
%         end
%         return;
%     end
% end
% 
% 
% if strcmp(wave_type, 'phase_array_transducer') == 1    % Ref. notebook. 2020.9.5. and Eq.(6-175) of E.G.William
%     if strcmp(direction, 'Z') == 1
%         if strcmp(symbol, 'positive') == 1
%             Anm = p_0*4*pi * (1i)^(n) * conj(sHarmonics(n,m,theta_Fourier,phi_Fourier)) * ...
%                 exp(1i*fluid_k*(sin(theta_Fourier)*cos(phi_Fourier)*x_translation + ...
%                 sin(theta_Fourier)*sin(phi_Fourier)*y_translation + cos(theta_Fourier)*z_translation));
%             anm = Anm / ((2*n + 1) * 1i^n); 
%         end
%         if strcmp(symbol, 'negative') == 1
%             Anm = p_0*4*pi * (-1i)^(n) * conj(sHarmonics(n,m,theta_Fourier,phi_Fourier)) * ...
%                 exp(-1i*fluid_k*(sin(theta_Fourier)*cos(phi_Fourier)*x_translation + ...
%                 sin(theta_Fourier)*sin(phi_Fourier)*y_translation + cos(theta_Fourier)*z_translation));
%             anm = Anm / ((2*n + 1) * 1i^n); 
%         end
%         return;
%     end
% end
% if strcmp(wave_type, 'phase_array_transducer_standing') == 1   %!!!
%     if strcmp(direction, 'Z') == 1
%         if strcmp(symbol, 'positive') == 1
%             Anm = 2*p_0*4*pi * cos(fluid_k*(sin(theta_Fourier)*cos(phi_Fourier)*x_translation + ...
%                 sin(theta_Fourier)*sin(phi_Fourier)*y_translation + cos(theta_Fourier)*z_translation) + n*pi/2) * ...
%                 conj(sHarmonics(n,m,theta_Fourier,phi_Fourier));
%             anm = Anm / ((2*n + 1) * 1i^n); 
%         end
%         return;
%     end
% end
% ===================================================================================================================


for ii = 1:length(phi)
    obj_func = p_b(ii, :) .* conj_Ynm(ii, :) .* sin(theta); 
    inner_int(ii) = trapz(theta, obj_func);
end
integral = trapz(phi, inner_int);

fluid_kb = fluid_k * b;
anm = integral / ((2*n + 1) * 1i^n * (sBessel(n, fluid_kb, 1)+eps));   % G.T.S@2011@Off-axis, Eq.(8)
Anm = anm * (2*n + 1) * 1i^n;                                    % more general explaination

%% 
