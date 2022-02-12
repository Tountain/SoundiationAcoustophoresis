function [p_i] = array_wave_function(wave_type, fluid_k, x_translation, y_translation, z_translation, ...
                        trans_radius, transducer_number, Pr_0, inter_dist, Rxyz)
%%
% This function is used to sum all the wave function from all
% transdcuers together, so only valid for feature
% "wave_type='phase_array_transducer'" .
%%

%% getting transducer position information

[phase_delay_vector, transducer] = phase_array_beam_shape_coeff(wave_type, 0, 0, fluid_k, trans_radius, transducer_number, 0);
transducer_number = size(transducer, 1);

%% the array wave function

% p_trans = cell(transducer_number, 1);
p_i = @(r, theta, phi) 0;
for ii = 1 : transducer_number
    % Cartesian system
    x = @(r, theta, phi) r*cos(phi')*sin(theta) - transducer(ii,1);
    y = @(r, theta, phi) r*sin(phi')*sin(theta) - transducer(ii,2);
    z = @(r, theta, phi) r*ones(length(phi),1)*cos(theta) - transducer(ii,3);

    x_rot = @(r, theta, phi) x(r, theta, phi) * Rxyz(1,1) + y(r, theta, phi) * Rxyz(2,1) + z(r, theta, phi) * Rxyz(3,1);
    y_rot = @(r, theta, phi) x(r, theta, phi) * Rxyz(1,2) + y(r, theta, phi) * Rxyz(2,2) + z(r, theta, phi) * Rxyz(3,2);
    z_rot = @(r, theta, phi) x(r, theta, phi) * Rxyz(1,3) + y(r, theta, phi) * Rxyz(2,3) + z(r, theta, phi) * Rxyz(3,3);
    % making translations
    x_rot_t = @(r, theta, phi) x_rot(r, theta, phi) + x_translation;
    y_rot_t = @(r, theta, phi) y_rot(r, theta, phi) + y_translation;
    z_rot_t = @(r, theta, phi) z_rot(r, theta, phi) + inter_dist + z_translation;

    r_rot_t = @(r, theta, phi) (x_rot_t(r, theta, phi).^2 + y_rot_t(r, theta, phi).^2 + z_rot_t(r, theta, phi).^2).^(0.5);
    theta_rot_t = @(r, theta, phi) acos(z_rot_t(r, theta, phi) ./ r_rot_t(r, theta, phi));
    p_trans = @(r, theta, phi) phase_delay_vector(ii) * Pr_0 * 2 * ...          % single piston-like wave function
        besselj(1, fluid_k * trans_radius * sin(theta_rot_t(r, theta, phi))) ./ ...
        (fluid_k * trans_radius * sin(theta_rot_t(r, theta, phi) + eps)) .* ...
        exp(1i * fluid_k * (r_rot_t(r, theta, phi))) ./ (r_rot_t(r, theta, phi) + eps);
    
    p_i = @(r, theta, phi) p_i(r, theta, phi) + p_trans(r, theta, phi);
end  


%% the array wave function performs rotation transformation

% % Cartesian system
% x = @(r, theta, phi) r*cos(phi')*sin(theta);
% y = @(r, theta, phi) r*sin(phi')*sin(theta);
% z = @(r, theta, phi) r*ones(length(phi),1)*cos(theta);
% % making rotation
% x_rot = @(r, theta, phi) x(r, theta, phi) * Rxyz(1,1) + y(r, theta, phi) * Rxyz(2,1) + z(r, theta, phi) * Rxyz(3,1);
% y_rot = @(r, theta, phi) x(r, theta, phi) * Rxyz(1,2) + y(r, theta, phi) * Rxyz(2,2) + z(r, theta, phi) * Rxyz(3,2);
% z_rot = @(r, theta, phi) x(r, theta, phi) * Rxyz(1,3) + y(r, theta, phi) * Rxyz(2,3) + z(r, theta, phi) * Rxyz(3,3);
% % making translations
% x_rot_t = @(r, theta, phi) x_rot(r, theta, phi) + x_translation;
% y_rot_t = @(r, theta, phi) y_rot(r, theta, phi) + y_translation;
% z_rot_t = @(r, theta, phi) z_rot(r, theta, phi) + inter_dist + z_translation;
% 
% r_rot_t = @(r, theta, phi) (x_rot_t(r, theta, phi).^2 + y_rot_t(r, theta, phi).^2 + z_rot_t(r, theta, phi).^2).^(0.5);
% theta_rot_t = @(r, theta, phi) acos(z_rot_t(r, theta, phi) ./ r_rot_t(r, theta, phi));
% p_i = @(r, theta, phi) Pr_0 * 2 * ...          % single piston-like wave function
%     besselj(1, fluid_k * trans_radius * sin(theta_rot_t(r, theta, phi))) ./ ...
%     (fluid_k * trans_radius * sin(theta_rot_t(r, theta, phi) + eps)) .* ...
%     exp(1i * fluid_k * (r_rot_t(r, theta, phi))) ./ (r_rot_t(r, theta, phi) + eps);


%%