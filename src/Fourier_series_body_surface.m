function [avg_r, An, Bn, theta_body, Err] = Fourier_series_body_surface(N, sample_number)
%%
% Expanding radial function 'r_body' in a Fourier series relative to the
% polar angle 'theta_body', and return the Fourier series coefficients. 
% Refs: 
%    Eq. (5) in D.T.DiPerna@1994@(JASA);
%    Eq. (21) in D.B.Reeder@2003@(JASA).
%%

parameters;

if irregular_body == 0
    error('The particle is not a irregular object!\n');
end

%% Fourier expansion factor 'gnm'

T = 2*pi;		% 'T' is set to the period of the body of polar angle (0 to 2*pi), actually, only [0,pi] is meaningful.

[r_body, theta_body] = body_surface(sample_number);
a0 = r_body .* cos(2*pi/T * 0 * theta_body);
A0 = trapz(theta_body, a0) * 2/T;
avg_r = A0/2;
An = zeros(1, N);
Bn = zeros(1, N);
for nn = 1 : N
    an = r_body .* cos(2*pi/T * nn * theta_body);
    bn = r_body .* sin(2*pi/T * nn * theta_body);
    An(nn) = trapz(theta_body, an) * 2/T;
    Bn(nn) = trapz(theta_body, bn) * 2/T;
end

%% error

valid_theta_sample = sample_number/2;               % 'theta' only meaningful in [0, pi].
theta_body_Fourier = theta_body(1:valid_theta_sample);
r_body_Fourier = ones(1, valid_theta_sample)*avg_r;
for nn = 1 : N
    r_body_Fourier = r_body_Fourier + An(nn) * cos(2*pi/T*nn*theta_body_Fourier) + Bn(nn) * sin(2*pi/T*nn*theta_body_Fourier);
end

Err = sqrt(sum((r_body_Fourier - r_body(1:valid_theta_sample)).^2) / valid_theta_sample);        % RMSE

%% visualize

% figure_FontSize = 18;
% 
% figure;
% hold on; box on;
% h1 = plot(r_body(1:sample_number/20:end).*cos(theta_body(1:sample_number/20:end)), r_body(1:sample_number/20:end).*sin(theta_body(1:sample_number/20:end)), 'ob'); 
% h2 = plot(r_body_Fourier.*cos(theta_body_Fourier), r_body_Fourier.*sin(theta_body_Fourier), '--r', 'linewidth', 2);
% plot([-particle_major_length/2,particle_major_length/2],[0,0],'-.k');
% legend([h1, h2], '\rm{\fontname{Times new roman}Locus of irregular body}', '\rm{\fontname{Times new roman}Fitting of Fourier series}')
% title(['\it{\fontname{Times new roman}E_{\rm{RMS}}} \rm{\fontname{Times new roman} = } ', num2str(roundn(Err*100,-3)), '%'], ...
%     'FontSize', figure_FontSize);
% xlabel('\rm{\fontname{Times new roman}\it{z}\rm{-axis [m]}}');
% ylabel('\rm{\fontname{Times new roman}\it{x}\rm{-axis [m]}}');
% set(get(gca,'XLabel'),'FontSize',figure_FontSize);
% set(get(gca,'YLabel'),'FontSize',figure_FontSize);
% set(findobj('Fontsize',10),'fontsize', figure_FontSize-2);
% axis([-particle_major_length/2, particle_major_length/2, 0, max(r_body.*sin(theta_body))]);
% axis equal;
% hold off;

%%