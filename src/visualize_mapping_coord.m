function visualize_mapping_coord(sample_number)
%%
% This function is used to visualize the exterior mapping coordinate system
% based on the mapping coefficients 'Cn'.
%
% E.g., (refer to my note on 2021.4.28, point 8)
%     Cn=[0.005, 0]; (spherical)
%     Cn=[0.005, 0, 0.002]; (spheroidal)
%     Cn=[0.005, 0, 0, 0.001]; (trangular-cone)
%     Cn=[0.005, 0, 0, 0, 0.0005]; (diamond)
%
% Note: sample_number = 100;
%%

parameters;         % getting mapping coefficients 'Cn'.

if irregular_body == 0
    Cn = [particle_radius, 0];
end

T = 2*pi;
% Cn = Cn/2;
uu1 = linspace(0,1.2,13);
ww1 = linspace(0, T, sample_number);

G_u = zeros(length(uu1), sample_number);
for ii = 1 : length(uu1)
    rho1 = uu1(ii) + ww1*1i;
    G_u(ii, :) = Cn(1) * exp(2*pi/T * rho1);
    for nn = 2 : length(Cn)
        G_u(ii, :) = G_u(ii, :) + Cn(nn) * exp(- 2*pi/T*(nn-2) * rho1);
    end
end

uu2 = linspace(0,1.3,27);
ww2 = linspace(1*T, 2*T, 21);
G_w = zeros(length(ww2), length(uu2));
for jj = 1 : length(ww2)
    rho2 = uu2 + ww2(jj)*1i;
    G_w(jj, :) = Cn(1) * exp(2*pi/T * rho2);
    for nn = 2 : length(Cn)
        G_w(jj, :) = G_w(jj, :) + Cn(nn) * exp(- 2*pi/T*(nn-2) * rho2);
    end
end

figure;
figure_FontSize = 18;
hold on; box on;
for ii = 1 : length(uu1)
    if ii == 1
        plot(imag(G_u(ii, :)), real(G_u(ii, :)), 'b', 'linewidth', 3);
    else
        plot(imag(G_u(ii, :)), real(G_u(ii, :)), 'k', 'linewidth', 1);
    end
end
for jj = 1 : length(ww2)
    plot(imag(G_w(jj, :)), real(G_w(jj, :)), 'r', 'linewidth', 1);
end
% plot([-max(max(abs(G_u(3,:)))),max(max(abs(G_u(3,:))))],[0,0],'-.k');
% title(['\rm{\fontname{Times new roman}Sample number=}', num2str(sample_number), ...
%     '\rm{\fontname{Times new roman}; L=}', num2str(L), ...
%     '\rm{\fontname{Times new roman}; J=}', num2str(J), '.'], ...
%     'FontSize', figure_FontSize);
xlabel('\rm{\fontname{Times new roman}\it{x {\rm{or}} f} {\rm{[m]}}}');
ylabel('\rm{\fontname{Times new roman}\it{z {\rm{or}} g} {\rm{[m]}}}');
set(get(gca,'XLabel'),'FontSize',figure_FontSize);
set(get(gca,'YLabel'),'FontSize',figure_FontSize);
set(findobj('Fontsize',10),'fontsize', figure_FontSize-2);
axis equal;
axis([-0.92,0.92,-0.6,0.6]*2*10^-4);
hold off;

%%
