function [X,Y,Z] = establish_3D_model(sample_number)
%%
% This function is used to establish the 3D geometric model using the
% mapping coefficients 'Cn' and the mapping functions 'g' and 'f'. The data
% should be acceptable and compatible with COMSOL.
%
%
% The data saving tables 'X, Y, Z' for 3D geometric model:
%   X = [(w1,v1)  (w1,v2)  ... (w1, vn)      Y = [(w1,v1)  (w1,v2)  ... (w1, vn)
%        (w2,v1)  (w2,v2)  ... (w2, vn)           (w2,v1)  (w2,v2)  ... (w2, vn)
%          ...      ...           ...               ...      ...           ...
%        (wn,v1)  (wn,v2)  ... (wn, vn)];         (wn,v1)  (wn,v2)  ... (wn, vn)];         
%
%   Z = [(w1,v1)  (w1,v2)  ... (w1, vn)
%        (w2,v1)  (w2,v2)  ... (w2, vn)
%          ...      ...           ...
%        (wn,v1)  (wn,v2)  ... (wn, vn)];
%%

% sample_number = 100;

%% getting the mapping functions 'g' and 'f'

parameters;         % getting mapping coefficients 'Cn'.

if irregular_body == 0
    Cn = [particle_radius, 0];
end

T = 2*pi;
uu1 = 0;
ww1 = linspace(0, T/2, sample_number);
vv1 = linspace(0, T, sample_number);

rho1 = uu1 + ww1*1i;
G_u = Cn(1) * exp(2*pi/T * rho1);
for nn = 2 : length(Cn)
    G_u = G_u + Cn(nn) * exp(- 2*pi/T*(nn-2) * rho1);
end

g_data = real(G_u);
f_data = imag(G_u);

%% visualize the 2D slice surface on xOz or fOg plane

% figure;
% figure_FontSize = 18;
% box on;
% plot(f_data, g_data, 'b', 'linewidth', 3);
% xlabel('\rm{\fontname{Times new roman}\it{x {\rm{or}} f} {\rm{[m]}}}');
% ylabel('\rm{\fontname{Times new roman}\it{z {\rm{or}} g} {\rm{[m]}}}');
% set(get(gca,'XLabel'),'FontSize',figure_FontSize);
% set(get(gca,'YLabel'),'FontSize',figure_FontSize);
% set(findobj('Fontsize',10),'fontsize', figure_FontSize-2);
% axis equal;

%% rotating the slice along the azimuthal coordinate from 0 to 2pi to establish the 3D model

X = zeros(sample_number, sample_number);
Y = zeros(sample_number, sample_number);
Z = zeros(sample_number, sample_number);
for ii_w = 1 : sample_number                    % new coordinate system mapping back to Cartesian coordinate system
    for ii_v = 1: sample_number
        X(ii_w, ii_v) = f_data(ii_w) * cos(vv1(ii_v));
        Y(ii_w, ii_v) = f_data(ii_w) * sin(vv1(ii_v));
        Z(ii_w, ii_v) = g_data(ii_w);
    end
end

% rotate the particle based on 'theta_rotation'
% Rx=[1 0 0;
%     0 cos(-theta_rotation(1)) -sin(-theta_rotation(1)); 
%     0 sin(-theta_rotation(1)) cos(-theta_rotation(1))];
% Ry=[cos(-theta_rotation(2)) 0 sin(-theta_rotation(2));
%     0 1 0; 
%     -sin(-theta_rotation(2)) 0 cos(-theta_rotation(2))];
% Rz=[cos(-theta_rotation(3)) -sin(-theta_rotation(3)) 0; 
%     sin(-theta_rotation(3)) cos(-theta_rotation(3)) 0; 
%     0 0 1];
% Rxyz = Rx * Ry * Rz;
% [r, c] = size(X);
% for ii = 1 : r
%     for jj = 1 : c
%     
%         T = [X(ii,jj), Y(ii,jj), Z(ii,jj)] * Rxyz;
%         X(ii,jj) = T(1);
%         Y(ii,jj) = T(2);
%         Z(ii,jj) = T(3);
%         
%     end
% end

% figure;
[M,N] = size(X);
col = ones(M,N,3);
col(:,:,1) = 0.5;
col(:,:,2) = 0.5;
col(:,:,3) = 0.5; % 'col' for gray color
h = surf(X, Y, Z, col);
xlabel('\rm{\fontname{Times new roman}\it{x} {\rm{[m]}}}');
ylabel('\rm{\fontname{Times new roman}\it{y} {\rm{[m]}}}');
zlabel('\rm{\fontname{Times new roman}\it{z} {\rm{[m]}}}');
axis equal; camlight; lighting gouraud;
set(h, 'edgecolor', 'none');
set(gca, 'FontName', 'Times new roman');


% [x_min, x_max, y_min, y_max] = mat2txt('irregular_body_data', X, Y, Z);
% fprintf('The minimum and maximum value along x-axis are %f and %f, along y-axis are %f and %f.\n', x_min, x_max, y_min, y_max);

% surf2stl('irregular_body_data_spheroidal_0.2(1mm)', X, Y, Z, 'binary');
% surf2stl('irregular_body_data_cone_0.125(1mm)', X, Y, Z, 'binary');
% surf2stl('irregular_body_data_diamond_0.1(1mm)', X, Y, Z, 'binary');
% surf2stl('irregular_body_data_RBC_(-0.5-0.1-0.03)(1mm)', X, Y, Z, 'binary');
surf2stl('particle_data', X, Y, Z, 'binary');

%%