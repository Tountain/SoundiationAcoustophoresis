function [g, dev_g, g_partial_u, f, dev_f, f_partial_u] = mapping_functions(Cn, theta_w, delta_u0)
%%
% This function is used to calculate the mapping functions of 'g' and its
% radial deviation 'g_partial_u' as well as 'f' and its radial deviation
% 'f_partial_u' on the body surface 'u=0' at different polar coordinate
% 'theta_w', using the obtained mapping coefficients 'Cn' of the specific
% pattern.  
%
% E.g., (refer to my note on 2021.4.28, point 8)
%     Cn=[0.005, 0]; (spherical)
%     Cn=[0.005, 0, 0.002]; (spheroidal)
%     Cn=[0.005, 0, 0, 0.001]; (trangular-cone)
%     Cn=[0.005, 0, 0, 0, 0.0005]; (diamond)
%%

% mapping functions 'g' and 'f'

T = 2*pi;
u0 = 0;
rho0 = u0 + theta_w*1i;

G_u0 = Cn(1) * exp(2*pi/T * rho0);
for nn = 2 : length(Cn)
    G_u0 = G_u0 + Cn(nn) * exp(- 2*pi/T*(nn-2) * rho0);
end

g = real(G_u0);
f = imag(G_u0);

% deviations of 'g' and 'f' along radial coordinate 'u'

G_partial_u0 = Cn(1) * exp(2*pi/T * rho0);
for nn = 2 : length(Cn)
    G_partial_u0 = G_partial_u0 + Cn(nn) * (2-nn) * exp(- 2*pi/T*(nn-2) * rho0);
end

g_partial_u = real(G_partial_u0);
f_partial_u = imag(G_partial_u0);

% dev_g and dev_f
u0 = 0 + delta_u0;

rho0 = u0 + theta_w*1i;

G_u0 = Cn(1) * exp(2*pi/T * rho0);
for nn = 2 : length(Cn)
    G_u0 = G_u0 + Cn(nn) * exp(- 2*pi/T*(nn-2) * rho0);
end

dev_g = real(G_u0);
dev_f = imag(G_u0);


%% rotation the object on fog-plane

% theta_rotation = pi/6;
% if theta_rotation ~= 0;
%     f_rot = (f*cos(theta_rotation) - g*sin(theta_rotation));
%     g_rot = (f*sin(theta_rotation) + g*cos(theta_rotation));
% 
%     g = g_rot;
%     f = f_rot;
% 
%     dev_f_rot = (dev_f*cos(theta_rotation) - dev_g*sin(theta_rotation));
%     dev_g_rot = (dev_f*sin(theta_rotation) + dev_g*cos(theta_rotation));
% 
%     dev_g = dev_g_rot;
%     dev_f = dev_f_rot;
% end


%% checking for spheroidal cases

% g_checking = (Cn(1) + Cn(3)) .* cos(theta_w);
% f_checking = (Cn(1) - Cn(3)) .* sin(theta_w);
% g_checking_u = (Cn(1) - Cn(3)) .* cos(theta_w);
% f_checking_u = (Cn(1) + Cn(3)) .* sin(theta_w);
% g_checking;

%%