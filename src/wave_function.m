function [p_i] = ...
    wave_function(p_0, wave_type, direction, dir_sign, ...
    fluid_k, x_translation, y_translation, z_translation, trans_radius, transducer_number, Pr_0, inter_dist, theta_rotation)
%%
% this function is used to give the expression of incident progress wave
% based on specific variable "wave_type".
%
% sub-function of "parameters.m"
%%

p_i = 0;         % wave function, no zero return for 'plain', 'standing_plain', 'zero-Bessel' and 'non-zero-Bessel'.
% G_i = 0;         % transerve shape function, no zero return for 'single_transducer' and 'phase_array_transducer'.

%%

% plain wave
if strcmp(wave_type, 'plain') == 1
    if strcmp(direction, 'X') == 1          % x-direaction.                     
        p_i = @(r,theta,phi) p_0*exp(dir_sign * 1i*fluid_k*(r*cos(phi')*sin(theta) + x_translation));
    end
    if strcmp(direction, 'Y') == 1          % y-direaction.                                                        
        p_i = @(r,theta,phi) p_0*exp(dir_sign * 1i*fluid_k*(r*sin(phi')*sin(theta) + y_translation));
    end
    if strcmp(direction, 'Z') == 1          % z-direaction.
        p_i = @(r,theta,phi) p_0*ones(length(phi),1)*exp(dir_sign * 1i*fluid_k*(r*cos(theta) + z_translation));
    end
end

%%
% plain standing wave  
if strcmp(wave_type, 'standing_plain') == 1
    if strcmp(direction, 'X') == 1          % x-direaction.                     
        p_i = @(r,theta,phi) p_0*exp(1i*fluid_k*(r*cos(phi')*sin(theta) + x_translation)) + ...
                            p_0*exp(-1i*fluid_k*(r*cos(phi')*sin(theta) + x_translation));
    end
    if strcmp(direction, 'Y') == 1          % y-direaction.                                                        
        p_i = @(r,theta,phi) p_0*exp(1i*fluid_k*(r*sin(phi')*sin(theta) + y_translation)) + ...
                            p_0*exp(-1i*fluid_k*(r*sin(phi')*sin(theta) + y_translation));
    end
    if strcmp(direction, 'Z') == 1          % z-direaction.
        p_i = @(r,theta,phi) p_0*ones(length(phi),1)*exp(1i*fluid_k*(r*cos(theta) + z_translation)) + ...
                            p_0*ones(length(phi),1)*exp(-1i*fluid_k*(r*cos(theta) + z_translation));
    end
    if strcmp(direction, 'XZ') == 1         % plain two-direction (X&Z-axis) standing wave 
        p_i = @(r,theta,phi) p_0*exp(1i*fluid_k*(r*cos(phi')*sin(theta) + x_translation)) + ...
                            p_0*exp(-1i*fluid_k*(r*cos(phi')*sin(theta) + x_translation)) + ...
                            p_0*ones(length(phi),1)*exp(1i*fluid_k*(r*cos(theta) + z_translation)) + ...
                            p_0*ones(length(phi),1)*exp(-1i*fluid_k*(r*cos(theta) + z_translation));
    end
end
                    
%%
% zero-order Bessel Beam; WHERE "beta" is the half-cone angle and "(x0,y0)" is the launch position of the Bessel beam.
if strcmp(wave_type, 'zero-Bessel') == 1
    if strcmp(direction, 'X') == 1          % x-direaction.
        p_i = @(r,beta,x0,y0,z0,theta,phi) p_0*exp(dir_sign * 1i*fluid_k* (r*cos(phi')*sin(theta) + x_translation) *cos(beta)) .* ...                         
            besselj(0,fluid_k * sqrt(((r*sin(phi')*sin(theta) + y_translation)-y0).^2 + ((r*ones(length(phi),1)*cos(theta) + z_translation)-z0).^2) * sin(beta));
    end
    if strcmp(direction, 'Y') == 1          % y-direaction.
        p_i = @(r,beta,x0,y0,z0,theta,phi) p_0*exp(dir_sign * 1i*fluid_k* (r*sin(phi')*sin(theta) + y_translation) *cos(beta)) .* ...                         
            besselj(0,fluid_k * sqrt(((r*cos(phi')*sin(theta) + x_translation)-x0).^2 + ((r*ones(length(phi),1)*cos(theta) + z_translation)-z0).^2) * sin(beta));
    end
    if strcmp(direction, 'Z') == 1          % z-direaction.
        p_i = @(r,beta,x0,y0,z0,theta,phi) p_0*(ones(length(phi),1)*exp(dir_sign * 1i*fluid_k* (r*cos(theta) + z_translation) *cos(beta))) .* ...           
            besselj(0,fluid_k * sqrt(((r*cos(phi')*sin(theta) + x_translation)-x0).^2+((r*sin(phi')*sin(theta) + y_translation)-y0).^2) * sin(beta));
    end
end

%%
% other Bessel Beam; WHERE "Bn" is the order of Bessel beam, "beta" is the half-cone angle and "(x0,y0)" is the launch position of the Bessel beam.
if strcmp(wave_type, 'non-zero-Bessel') == 1
    if strcmp(direction, 'X') == 1          % x-direaction.
        p_i = @(r,Bn,beta,x0,y0,z0,theta,phi) p_0*exp(dir_sign * 1i*fluid_k*r*cos(phi')*sin(theta)*cos(beta)) .* ...                   
            besselj(Bn,fluid_k* sqrt((r*sin(phi')*sin(theta)-y0).^2 + (r*ones(length(phi),1)*cos(theta)-z0).^2) * sin(beta)) .* ...
            exp(1i*Bn*atan(dir_sign * (r*ones(length(phi),1)*cos(theta)-z0+eps) ./ (r*sin(phi')*sin(theta)-y0+eps)));
    end
    if strcmp(direction, 'Y') == 1          % y-direaction.
        p_i = @(r,Bn,beta,x0,y0,z0,theta,phi) p_0*exp(dir_sign * 1i*fluid_k*r*sin(phi')*sin(theta)*cos(beta)) .* ...                   
            besselj(Bn,fluid_k* sqrt((r*cos(phi')*sin(theta)-x0).^2 + (r*ones(length(phi),1)*cos(theta)-z0).^2) * sin(beta)) .* ...
            exp(1i*Bn*atan(dir_sign * (r*cos(phi')*sin(theta)-x0+eps) ./ (r*ones(length(phi),1)*cos(theta)-z0+eps)));
    end
    if strcmp(direction, 'Z') == 1          % z-direaction.
        p_i = @(r,Bn,beta,x0,y0,z0,theta,phi) p_0*(ones(length(phi),1)*exp(dir_sign * 1i*fluid_k*r*cos(theta)*cos(beta))) .* ...         
            besselj(Bn,fluid_k* sqrt((r*cos(phi')*sin(theta)-x0).^2 + (r*sin(phi')*sin(theta)-y0).^2) * sin(beta)) .* ...
            exp(1i*Bn*atan(dir_sign * (r*sin(phi')*sin(theta)-y0+eps) ./ (r*cos(phi')*sin(theta)-x0+eps)));
    end
end

%% Single transducer incident wave
% trans_radius = 5 * 10^-3;       % radius of transducer is 5mm
% cos_beta = sqrt(1 - (1/(fluid_k * trans_radius))^2);

% G_i(xx,yy):                                  o    y       (vertical view of OXYZ system)
%   [(xx1, yy1) (xx1, yy2) ... (xx1, yyn)       |¡ª¡ª¡ª>
%    (xx2, yy1) (xx2, yy2) ... (xx2, yyn)     x |
%         .                                     v
%         .
%         .
%    (xxn, yy1) (xxn, yy2) ... (xxn, yyn)]
if strcmp(wave_type, 'single_transducer') == 1 || strcmp(wave_type, 'single_transducer_standing') == 1
    if strcmp(direction, 'X') == 1          % x-direaction.
        error('x: please setting the wave propagates along +z-axes direction.\n');
    end
    if strcmp(direction, 'Y') == 1          % y-direaction.
        error('y: please setting the wave propagates along +z-axes direction.\n');
    end
    if strcmp(direction, 'Z') == 1          % z-direaction.
%         p_i = @(r, theta, phi) p_0 * besselj(0, fluid_k * trans_radius * ...
%             ones(length(phi),1) * sin(theta)) .* ...         
%             exp(1i * fluid_k * (r)) ./ ...
%             (r);
%         p_i = @(r, theta, phi) p_0 * ones(length(phi),1) * exp(1i * fluid_k * (r * cos(theta) + z_translation) * cos_beta) .* ...           % 0th Bessel beam 
%             besselj(0, sqrt(((r*cos(phi')*sin(theta) + x_translation)).^2+((r*sin(phi')*sin(theta) + y_translation)).^2) / trans_radius);
%         p_i = @(r,theta,phi) p_0*exp(1i*fluid_k*(r*ones(length(phi),1)*cos(theta) + z_translation)+1i*fluid_k*(r*cos(phi')*sin(theta) + x_translation));
%         G_i = @(xx, yy) p_0 * ...          % transerve shape function
%             besselj(0, 2.5 * sqrt((xx + x_translation).^2 + (yy + y_translation).^2) / trans_radius) .* ...
%             step_fun(sqrt((xx + x_translation).^2 + (yy + y_translation).^2) / trans_radius);
        r_t = @(r, theta, phi) ((r*cos(phi')*sin(theta)).^2 + (r*sin(phi')*sin(theta)).^2 + (r*ones(length(phi),1)*cos(theta) + inter_dist).^2).^(0.5);
        theta_t = @(r, theta, phi) acos((r*ones(length(phi),1)*cos(theta) + inter_dist) ./ r_t(r, theta, phi));
        p_i = @(r, theta, phi) Pr_0 * 2 * ...          % transerve shape function
            besselj(1, fluid_k * trans_radius * sin(theta_t(r, theta, phi))) ./ ...
            (fluid_k * trans_radius * sin(theta_t(r, theta, phi) + eps)) .* ...
            exp(1i * fluid_k * (r_t(r, theta, phi))) ./ (r_t(r, theta, phi) + eps);
    end
end

%% Phase array incident wave (based on addition theorem)

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
% dist = [0,0,inter_dist] * Rxyz;

% trans_radius = 5 * 10^-3;       % radius of transducer is 5mm
% cos_beta = sqrt(1 - (1/(fluid_k * trans_radius))^2);

if strcmp(wave_type, 'phase_array_transducer') == 1 || strcmp(wave_type, 'phase_array_transducer_standing') == 1
    if strcmp(direction, 'X') == 1          % x-direaction.
        error('x: please setting the wave propagates along +z-axes direction.\n');      
    end
    if strcmp(direction, 'Y') == 1          % y-direaction.
        error('y: please setting the wave propagates along +z-axes direction.\n');
    end
    if strcmp(direction, 'Z') == 1          % z-direaction.
        % Cartesian system
        x = @(r, theta, phi) r*cos(phi')*sin(theta);
        y = @(r, theta, phi) r*sin(phi')*sin(theta);
        z = @(r, theta, phi) r*ones(length(phi),1)*cos(theta);
        % making rotation
        x_rot = @(r, theta, phi) x(r, theta, phi) * Rxyz(1,1) + y(r, theta, phi) * Rxyz(2,1) + z(r, theta, phi) * Rxyz(3,1);
        y_rot = @(r, theta, phi) x(r, theta, phi) * Rxyz(1,2) + y(r, theta, phi) * Rxyz(2,2) + z(r, theta, phi) * Rxyz(3,2);
        z_rot = @(r, theta, phi) x(r, theta, phi) * Rxyz(1,3) + y(r, theta, phi) * Rxyz(2,3) + z(r, theta, phi) * Rxyz(3,3);
        % making translations
        x_rot_t = @(r, theta, phi) x_rot(r, theta, phi) + x_translation;
        y_rot_t = @(r, theta, phi) y_rot(r, theta, phi) + y_translation;
        z_rot_t = @(r, theta, phi) z_rot(r, theta, phi) + inter_dist + z_translation;
%         p_i = @(r, theta, phi) p_0 * ones(length(phi),1) * exp(1i * fluid_k * (r * cos(theta) + z_translation) * cos_beta) .* ...
%             besselj(0, sqrt(((r*cos(phi')*sin(theta) + x_translation)).^2+((r*sin(phi')*sin(theta) + y_translation)).^2) / trans_radius);
        r_rot_t = @(r, theta, phi) (x_rot_t(r, theta, phi).^2 + y_rot_t(r, theta, phi).^2 + z_rot_t(r, theta, phi).^2).^(0.5);
        theta_rot_t = @(r, theta, phi) acos(z_rot_t(r, theta, phi) ./ r_rot_t(r, theta, phi));
        p_i = @(r, theta, phi) Pr_0 * 2 * ...          % single piston-like wave function
            besselj(1, fluid_k * trans_radius * sin(theta_rot_t(r, theta, phi))) ./ ...
            (fluid_k * trans_radius * sin(theta_rot_t(r, theta, phi) + eps)) .* ...
            exp(1i * fluid_k * (r_rot_t(r, theta, phi))) ./ (r_rot_t(r, theta, phi) + eps);
%         p_i = @(r, theta, phi) p_0;
    end
end

% z = r*cos(theta)
% y = r*sin(theta)*sin(phi)
% x = r*sin(theta)*cos(phi)
% 
% theta = acos(z/r)
% r = x^2 + y^2 + z^2

%% plain wave with rotation transformation

Rx=Rx;
Ry=Ry;
Rz=Rz;
Rxyz = Rxyz;

if strcmp(wave_type, 'plain') == 1
    if strcmp(direction, 'Z') == 1          % z-direaction.
        % Cartesian system
        x = @(r, theta, phi) r*cos(phi')*sin(theta);
        y = @(r, theta, phi) r*sin(phi')*sin(theta);
        z = @(r, theta, phi) r*ones(length(phi),1)*cos(theta);
        % making rotation
%         x_rot = @(r, theta, phi) x(r, theta, phi) * Rxyz(1,1) + y(r, theta, phi) * Rxyz(2,1) + z(r, theta, phi) * Rxyz(3,1);
%         y_rot = @(r, theta, phi) x(r, theta, phi) * Rxyz(1,2) + y(r, theta, phi) * Rxyz(2,2) + z(r, theta, phi) * Rxyz(3,2);
        z_rot = @(r, theta, phi) x(r, theta, phi) * Rxyz(1,3) + y(r, theta, phi) * Rxyz(2,3) + z(r, theta, phi) * Rxyz(3,3);
        % making translations
%         x_rot_t = @(r, theta, phi) x_rot(r, theta, phi) + x_translation;
%         y_rot_t = @(r, theta, phi) y_rot(r, theta, phi) + y_translation;
        z_rot_t = @(r, theta, phi) z_rot(r, theta, phi) + inter_dist + z_translation;
        
        p_i = @(r,theta,phi) p_0*exp(dir_sign * 1i*fluid_k*z_rot_t(r, theta, phi));
    end
end


%% Phase array incident wave (NOT based on addition theorem)

Rx=Rx;
Ry=Ry;
Rz=Rz;
Rxyz = Rxyz;
% dist = [0,0,inter_dist] * Rxyz;

% trans_radius = 5 * 10^-3;       % radius of transducer is 5mm
% cos_beta = sqrt(1 - (1/(fluid_k * trans_radius))^2);

if strcmp(wave_type, 'phase_array_transducer2') == 1
    if strcmp(direction, 'X') == 1          % x-direaction.
        error('x: please setting the wave propagates along +z-axes direction.\n');      
    end
    if strcmp(direction, 'Y') == 1          % y-direaction.
        error('y: please setting the wave propagates along +z-axes direction.\n');
    end
    if strcmp(direction, 'Z') == 1          % z-direaction.
        p_i = array_wave_function(wave_type, fluid_k, x_translation, y_translation, z_translation, ...
                    trans_radius, transducer_number, Pr_0, inter_dist, Rxyz);
    end
end

%%

%% saving the old version

% % plain wave  
% if strcmp(wave_type, 'plain') == 1
%     BN = 0;     
%     BETA = 0;
%     X0 = 0;
%     Y0 = 0;
%     Z0 = 0;
% %     p_0 = 10;   % initial amplitude for plane wave
%     if strcmp(direction, 'X') == 1          % x-direaction.                     
%         p_i = @(r,theta,phi) p_0*exp(dir_sign * 1i*fluid_k*(r*cos(phi')*sin(theta) + x_translation));
%     end
%     if strcmp(direction, 'Y') == 1          % y-direaction.                                                        
%         p_i = @(r,theta,phi) p_0*exp(dir_sign * 1i*fluid_k*(r*sin(phi')*sin(theta) + y_translation));
%     end
%     if strcmp(direction, 'Z') == 1          % z-direaction.
%         p_i = @(r,theta,phi) p_0*ones(length(phi),1)*exp(dir_sign * 1i*fluid_k*(r*cos(theta) + z_translation));
%     end
% end
% 
% % plain standing wave  
% if strcmp(wave_type, 'standing_plain') == 1
%     if strcmp(direction, 'XZ') ~= 1
%         BN = 0; 
%         BETA = 0;
%         X0 = 0;
%         Y0 = 0;
%         Z0 = 0;
%         if strcmp(direction, 'X') == 1          % x-direaction.                     
%             p_i = @(r,theta,phi) p_0*exp(1i*fluid_k*(r*cos(phi')*sin(theta) + x_translation)) + ...
%                                 p_0*exp(-1i*fluid_k*(r*cos(phi')*sin(theta) + x_translation));
%         end
%         if strcmp(direction, 'Y') == 1          % y-direaction.                                                        
%             p_i = @(r,theta,phi) p_0*exp(1i*fluid_k*(r*sin(phi')*sin(theta) + y_translation)) + ...
%                                 p_0*exp(-1i*fluid_k*(r*sin(phi')*sin(theta) + y_translation));
%         end
%         if strcmp(direction, 'Z') == 1          % z-direaction.
%             p_i = @(r,theta,phi) p_0*ones(length(phi),1)*exp(1i*fluid_k*(r*cos(theta) + z_translation)) + ...
%                                 p_0*ones(length(phi),1)*exp(-1i*fluid_k*(r*cos(theta) + z_translation));
%         end
%         
%     else
%         % plain two-direction (X&Z-axis) standing wave 
%         direction = 'XZ';       % for two-direction (X&Z-axis) standing wave, we define direction always "XZ", then saving file looks order-well
%         BN = 0; 
%         BETA = 0;
%         X0 = 0;
%         Y0 = 0;
%         Z0 = 0;
%         p_i = @(r,theta,phi) p_0*exp(1i*fluid_k*(r*cos(phi')*sin(theta) + x_translation)) + ...
%                             p_0*exp(-1i*fluid_k*(r*cos(phi')*sin(theta) + x_translation)) + ...
%                              p_0*ones(length(phi),1)*exp(1i*fluid_k*(r*cos(theta) + z_translation)) + ...
%                             p_0*ones(length(phi),1)*exp(-1i*fluid_k*(r*cos(theta) + z_translation));
%     end
% end
%                     
% 
% % zero-order Bessel Beam; WHERE "beta" is the half-cone angle and "(x0,y0)" is the launch position of the Bessel beam.
% if strcmp(wave_type, 'zero-Bessel') == 1
%     BN = 0; 
%     BETA = pi/12;
%     X0 = particle_radius;
%     Y0 = 0;
%     Z0 = 0;
% %     p_0 = 10; % initial amplitude for Bessel beam
%     if strcmp(direction, 'X') == 1          % x-direaction.
%         p_i = @(r,beta,x0,y0,z0,theta,phi) p_0*exp(dir_sign * 1i*fluid_k* (r*cos(phi')*sin(theta) + x_translation) *cos(beta)) .* ...                         
%             besselj(0,fluid_k * sqrt(((r*sin(phi')*sin(theta) + y_translation)-y0).^2 + ((r*ones(length(phi),1)*cos(theta) + z_translation)-z0).^2) * sin(beta));
%     end
%     if strcmp(direction, 'Y') == 1          % y-direaction.
%         p_i = @(r,beta,x0,y0,z0,theta,phi) p_0*exp(dir_sign * 1i*fluid_k* (r*sin(phi')*sin(theta) + y_translation) *cos(beta)) .* ...                         
%             besselj(0,fluid_k * sqrt(((r*cos(phi')*sin(theta) + x_translation)-x0).^2 + ((r*ones(length(phi),1)*cos(theta) + z_translation)-z0).^2) * sin(beta));
%     end
%     if strcmp(direction, 'Z') == 1          % z-direaction.
%         p_i = @(r,beta,x0,y0,z0,theta,phi) p_0*(ones(length(phi),1)*exp(dir_sign * 1i*fluid_k* (r*cos(theta) + z_translation) *cos(beta))) .* ...           
%             besselj(0,fluid_k * sqrt(((r*cos(phi')*sin(theta) + x_translation)-x0).^2+((r*sin(phi')*sin(theta) + y_translation)-y0).^2) * sin(beta));
%     end
% end
% 
% % other Bessel Beam; WHERE "Bn" is the order of Bessel beam, "beta" is the half-cone angle and "(x0,y0)" is the launch position of the Bessel beam.
% if strcmp(wave_type, 'non-zero-Bessel') == 1
%     BN = 1;     % BN >= 1
%     BETA = pi/12;
%     X0 = 0;
%     Y0 = 0;
%     Z0 = 0;
% %     p_0 = 10; % initial amplitude for Bessel beam
%     if strcmp(direction, 'X') == 1          % x-direaction.
%         p_i = @(r,Bn,beta,x0,y0,z0,theta,phi) p_0*exp(dir_sign * 1i*fluid_k*r*cos(phi')*sin(theta)*cos(beta)) .* ...                   
%             besselj(Bn,fluid_k* sqrt((r*sin(phi')*sin(theta)-y0).^2 + (r*ones(length(phi),1)*cos(theta)-z0).^2) * sin(beta)) .* ...
%             exp(1i*Bn*atan(dir_sign * (r*ones(length(phi),1)*cos(theta)-z0+eps) ./ (r*sin(phi')*sin(theta)-y0+eps)));
%     end
%     if strcmp(direction, 'Y') == 1          % y-direaction.
%         p_i = @(r,Bn,beta,x0,y0,z0,theta,phi) p_0*exp(dir_sign * 1i*fluid_k*r*sin(phi')*sin(theta)*cos(beta)) .* ...                   
%             besselj(Bn,fluid_k* sqrt((r*cos(phi')*sin(theta)-x0).^2 + (r*ones(length(phi),1)*cos(theta)-z0).^2) * sin(beta)) .* ...
%             exp(1i*Bn*atan(dir_sign * (r*cos(phi')*sin(theta)-x0+eps) ./ (r*ones(length(phi),1)*cos(theta)-z0+eps)));
%     end
%     if strcmp(direction, 'Z') == 1          % z-direaction.
%         p_i = @(r,Bn,beta,x0,y0,z0,theta,phi) p_0*(ones(length(phi),1)*exp(dir_sign * 1i*fluid_k*r*cos(theta)*cos(beta))) .* ...         
%             besselj(Bn,fluid_k* sqrt((r*cos(phi')*sin(theta)-x0).^2 + (r*sin(phi')*sin(theta)-y0).^2) * sin(beta)) .* ...
%             exp(1i*Bn*atan(dir_sign * (r*sin(phi')*sin(theta)-y0+eps) ./ (r*cos(phi')*sin(theta)-x0+eps)));
%     end
% end