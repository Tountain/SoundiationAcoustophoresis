function [r_body, theta_body] = body_surface(sample_number)
%%
% this function is used to defined the surface coordiate 'r_body' of the
% irregular body on the spherical coordiate system.
% Note that: 
%       The 'r_body(theta_body)' corresponds the locus of points of the
%       body surface projecting on the XoZ plane. 
%%

parameters;

if irregular_body == 0
    error('The particle is not a irregular object!\n');
end


%% the surface function of the irregular body

% discrete data given in Cartesian coordinate system
x_series = [6 10 12 11 9 7 5 3];            % data input

x_start = 0;
x_end = 0;
x = [x_start, x_series, x_end];

z_start = particle_major_length/2;
z_end = -particle_major_length/2;
z1 = linspace(z_start, z_end, length(x_series)+2);
z2 = linspace(z_end, z_start, length(x_series)+2);
z = [z1, z2(2:end)];

% make sure the average radius of the body is 'particle_radius'
zz = linspace(z_start, z_end, sample_number);
xx = spline(z1, x, zz);
ratio = particle_radius / mean(xx);
x_series = x_series * ratio;        

x = [x_start, x_series, x_end];
for ii = 0:length(x_series)-1
    x = [x, -x_series(end-ii)];
end
x = [x, x_start];
 
% figure;
% plot(z, x, 'o');

% interpolation the data using 'spline' function
r = sqrt(x.^2 + z.^2);
theta1 = acos(z(1:length(x_series)+2) ./ r(1:length(x_series)+2));
theta = theta1;
for ii = 1:length(theta1)-1
    theta = [theta, pi + (pi-theta1(end-ii))];
end

theta_body = linspace(0, 2*pi, sample_number);
r_body = spline(theta, r, theta_body);

% visulize the body surface in xoz plane
% figure; box on;
% figure_FontSize = 18;
% plot(z, x, 'o', r_body.*cos(theta_body), r_body.*sin(theta_body), 'r');
% xlabel('\rm{\fontname{Times new roman}\it{z}\rm{-axis [m]}}');
% ylabel('\rm{\fontname{Times new roman}\it{x}\rm{-axis [m]}}');
% set(get(gca,'XLabel'),'FontSize',figure_FontSize);
% set(get(gca,'YLabel'),'FontSize',figure_FontSize);
% set(findobj('Fontsize',10),'fontsize', figure_FontSize-2);
% %axis([-particle_major_length/2, particle_major_length/2, 0, max(r_body.*sin(theta_body))]);
% axis equal;

%%

% zz = linspace(z_start, z_end, sample_number);
% xx = spline(z, x, zz);
% ratio = particle_radius / mean(xx);
% xx = xx * ratio;        % make sure the average radius of the body is 'particle_radius'

% translate the data to polar coordinate system
% r_body = sqrt(xx.^2 + zz.^2);
% theta_body = acos(zz ./ r_body);

% visulize the body surface in xoz plane
% figure;
% plot(x * ratio, z, 'o', r_body.*sin(theta_body), r_body.*cos(theta_body), 'r');

%%

