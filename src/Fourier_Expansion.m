function [gnm, kn_F, km_F, kp_F, Err] = Fourier_Expansion(NM, trans_radius, transducer_number, fluid_rho, fluid_c, freq, Pr_0, inter_dist, Dims_Fourier)
%%
%  Calculate the transerve shape function 'gnm' and the wave vector (kn_F,
%  km_F). 
%%

% parameters;             % for transducers' information, and 'fluid_c'
% if (strcmp(wave_type, 'single_transducer') || strcmp(wave_type, 'phase_array_transducer')) ~= 1
%     error('Fourier Expansion is useful for wave types of ''single_transducer'' and ''phase_array_transducer''.\n');
% end

%% Fourier expansion factor 'gnm'

if Dims_Fourier == 2
    
    % T = (trans_radius + 2*trans_radius * ceil(sqrt(transducer_number))) * 2;
    T = (trans_radius + 2*trans_radius * ceil(sqrt(3))) * 2;		% T = 0.05 mm, enough for 'transducer_number <= 9'
    x = linspace(-T, T, 100);
    y = x;

    N = NM(1);
    M = NM(2);
    kn_F = zeros(1, 2*N+1);
    km_F = zeros(1, 2*M+1);
    kp_F = 0;
    for n_F = -N : N
        for m_F = -M : M
            kn_F(n_F + N + 1) = 2*pi * n_F / T;
            km_F(m_F + M + 1) = 2*pi * m_F / T;
        end
    end

    gnm = zeros(N*2+1, M*2+1);
    for n_F = -N : N
        indice_1 = n_F+N+1;
        indice_2 = M+1;
        parfor m_F = -M : M

            first_int = zeros(length(x), length(y));
            second_int = zeros(1,length(y));
            
            for jj = 1:length(y)
                for ii = 1:length(x)
                    first_int(ii, jj) = G_i([x(ii), y(jj)], trans_radius, fluid_c, freq, Pr_0, inter_dist, Dims_Fourier) * exp(-1i*(kn_F(indice_1)*x(ii) + km_F(indice_2+m_F)*y(jj)));     % integration elements
                end
                second_int(jj) = trapz(x, first_int(:, jj));                        % first intergation along x variable
            end

            gnm(indice_1, indice_2+m_F) = trapz(y, second_int) / (T*T);             % second intergation along y variable

        end
        fprintf('The Fourier expansion is prepared at %d%%.\n', round(100*(n_F+N+1)/(2*N+1)));
    end
    
elseif Dims_Fourier == 3
    
    % T = (trans_radius + 2*trans_radius * ceil(sqrt(transducer_number))) * 2;
    T = (trans_radius + 2*trans_radius * ceil(sqrt(3))) * 2;        % T = 0.05 mm, enough for 'transducer_number <= 9'
    x = linspace(-T, T, 60);
    Tx = abs(max(x) - min(x));
    y = linspace(-T, T, 60);
    Ty = abs(max(y) - min(y));
    z = linspace(-T, T, 60);              % z should at least range from the transducer surface to the center of the particle
    Tz = abs(max(z) - min(z));
%     x = linspace(-inter_dist-2*trans_radius, inter_dist+2*trans_radius, 40);
%     y = linspace(-inter_dist-2*trans_radius, inter_dist+2*trans_radius, 40);
%     z = linspace(0*trans_radius, inter_dist+2*trans_radius, 20);
    
    N = NM(1);
    M = NM(2);
    P = NM(3);
    kn_F = zeros(1, 2*N+1);
    km_F = zeros(1, 2*M+1);
    kp_F = zeros(1, 2*P+1);
    for n_F = -N : N
        for m_F = -M : M
            for p_F = -P : P
                kn_F(n_F + N + 1) = 2*pi * n_F / Tx;
                km_F(m_F + M + 1) = 2*pi * m_F / Ty;
                kp_F(p_F + P + 1) = 2*pi * p_F / Tz;
            end
        end
    end

    gnm = zeros(N*2+1, M*2+1, P*2+1);
    for n_F = -N : N
        for m_F = -M : M
            indice_1 = n_F+N+1;
            kn_F_indice_1 = kn_F(indice_1);
            indice_2 = m_F+M+1;
            km_F_indice_2 = km_F(indice_2);
            indice_3 = P+1;
            parfor p_F = -P : P
                
                first_int = zeros(length(x), length(y), length(z));
                second_int = zeros(length(y), length(z));
                third_int = zeros(1, length(z));
                
%                 if indice_1 == 0 && indice_2 == 0 && (indice_3+p_F == 0)    % vec(k) = 0 means no wave.
%                     continue;
%                 else
                    for kk = 1:length(z)
                        for jj = 1:length(y)
                            for ii = 1:length(x)
                                first_int(ii, jj, kk) = G_i([x(ii), y(jj), z(kk)], trans_radius, fluid_c, freq, Pr_0, inter_dist, Dims_Fourier) * ...
                                    exp(-1i*(kn_F_indice_1*x(ii) + km_F_indice_2*y(jj) + kp_F(indice_3+p_F)*z(kk)));     % integration elements
                            end
                            second_int(jj, kk) = trapz(x, first_int(:, jj, kk));                        % first intergation along x variable
                        end
                        third_int(kk) = trapz(y, second_int(:, kk));                             % second intergation along x variable
                    end

                    gnm(indice_1, indice_2, indice_3+p_F) = trapz(z, third_int) / (Tx*Ty*Tz);             % third intergation along y variable
%                 end
                
            end
        end
        fprintf('The Fourier expansion is prepared at %d%%.\n', round(100*(n_F+N+1)/(2*N+1)));
    end
    
end

%% polar direction and graphs on xoz-plane

% if Dims_Fourier == 3
%     
%     k = 2*pi*freq / fluid_c; 
% %     v0 = - Pr_0 / fluid_rho / fluid_c / k / (trans_radius^2) / 1i * 2; 
% %     I0 = (fluid_rho * fluid_c * k^2 * v0^2 * pi*trans_radius^2) / (8 * pi * trans_radius^2);      % axial intensity; refer to Eq.(6.35) in Cheeke@2002@chapter 6
%     r_polar = inter_dist;
%     exp_rr = exp(1i * k * (r_polar));
%     theta_polar = linspace(0, 2*pi, 100);
%     x_polar = r_polar * sin(theta_polar);
%     y_polar = 0;
%     z_polar = r_polar * cos(theta_polar);
%     G_approx = zeros(1, length(theta_polar));
%     G_exact = G_approx;
%     for ii = 1:length(theta_polar)
%         for n_F = -N : N
%             for m_F = -M : M
%                 for p_F = -P : P
%                     G_approx(ii) = G_approx(ii) + gnm(n_F+N+1, m_F+M+1, p_F+P+1) * ...
%                         exp(1i*(kn_F(n_F + N + 1)*x_polar(ii) + km_F(m_F + M + 1)*y_polar + kp_F(p_F + P + 1)*z_polar(ii)));
%                 end
%             end
%         end
%         G_exact(ii) = G_i([x_polar(ii), y_polar, z_polar(ii)], trans_radius, fluid_c, freq, Pr_0, inter_dist, Dims_Fourier);
%     end
%     
%     p_amp_exact = G_exact / exp_rr;                             % complex pressure amplitude of exact wave
%     v_amp_exact = - p_amp_exact / (fluid_rho * fluid_c);        % momentum equation; complex velocity amplitude of exact wave
%     p_amp_approx = G_approx / exp_rr;                           % complex pressure amplitude of approximation wave
%     v_amp_approx = - p_amp_approx / (fluid_rho * fluid_c);      % momentum equation; complex velocity amplitude of exact wave
%     
%     % I_av = 0.5*Re(complex_amplitude_p * conj_complex_amplitude_v);
%     I_av_exact = 0.5 * real(p_amp_exact .* conj(v_amp_exact));   % time-average intensity of exact wave at position 'r_polar'
%     I_av_approx = 0.5 * real(p_amp_approx .* conj(v_amp_approx));% time-average intensity of approximation wave at position 'r_polar'
%     
%     figure(30);
%     polar(theta_polar, (I_av_exact), 'r');hold on;
%     polar(theta_polar, (I_av_approx), 'k--');
%     
% end


%% approximate result for original function 'G_i'

if Dims_Fourier == 2
    
    x0 = linspace(-trans_radius*2,trans_radius*2, 30);
    y0 = linspace(0,trans_radius*2, 5);
    G_approx = zeros(length(x0), length(y0));
    G_exact = G_approx;
    for ii = 1:length(x0)
        for jj = 1:length(y0)
            for n_F = -N : N
                for m_F = -M : M

                    G_approx(ii, jj) = G_approx(ii, jj) + gnm(n_F+N+1, m_F+M+1) * exp(1i*(kn_F(n_F + N + 1)*x0(ii) + km_F(m_F + M + 1)*y0(jj)));
                    
                end
            end
            G_exact(ii, jj) = G_i([x0(ii), y0(jj)], trans_radius, fluid_c, freq, Pr_0, inter_dist, Dims_Fourier);
        end
    end

elseif Dims_Fourier == 3
    
    x0 = linspace(-T/2, T/2, 50);
    y0 = linspace(-T/2, T/2, 50);
    z0 = linspace(-T, T, 21);
    G_approx = zeros(length(x0), length(y0), length(z0));
    G_exact = G_approx;
    for ii = 1:length(x0)
        for jj = 1:length(y0)
            for kk = 1:length(z0)
                for n_F = -N : N
                    for m_F = -M : M
                        for p_F = -P : P
                            G_approx(ii, jj, kk) = G_approx(ii, jj, kk) + gnm(n_F+N+1, m_F+M+1, p_F+P+1) * ...
                                exp(1i*(kn_F(n_F + N + 1)*x0(ii) + km_F(m_F + M + 1)*y0(jj) + kp_F(p_F + P + 1)*z0(kk)));
                        end
                    end
                end
                G_exact(ii, jj, kk) = G_i([x0(ii), y0(jj), z0(kk)], trans_radius, fluid_c, freq, Pr_0, inter_dist, Dims_Fourier);
            end
        end
    end
    
end

%% graph for approximate result 'G_approx' and original function 'G_i'

if Dims_Fourier == 2
    
    nMax1 = max(max(G_approx));
    nMax2 = max(max(G_exact));
    nMax = max([nMax1, nMax2]);
    nMin1 = min(min(G_approx));
    nMin2 = min(min(G_exact));
    nMin = min([nMin1, nMin2]);
    nG_approx = (G_approx - nMin) / (nMax - nMin);
    nG_exact = (G_exact - nMin) / (nMax - nMin);
    Err = sqrt(sum(sum((nG_approx - nG_exact).^2)) / (length(x0)*length(y0)));        % RMSE

    figure; hold on;grid on;
    resolution = 0.4;
    ii = 1;
    for jj = 1:length(y0)
        if (resolution*ii) <= 1
            color(ii,:) = [resolution*ii 0 0];
        elseif (resolution*ii) > 1 && (resolution*ii) <= 2
            color(ii,:) = [1 resolution*ii-1 0];
        elseif (resolution*ii) > 2 && (resolution*ii) <= 3
            color(ii,:) = [0 1 resolution*ii-2];
        else
            error('color out of range!\n');
        end
        plot3(x0(:), linspace(y0(jj),y0(jj),length(x0)), real(G_approx(:, jj)), 'o', 'color', color(ii,:));
        plot3(x0(:), linspace(y0(jj),y0(jj),length(x0)), real(G_exact(:, jj)), 'color', color(ii,:), 'linewidth', 1.5);
    %     figure; hold on;
    %     plot(x0(:), real(G_approx(:, jj)), 'bo');
    %     plot(x0(:), G_exact(:, jj), 'r', 'linewidth', 1.5);
        axis([-trans_radius*2,trans_radius*2,0,trans_radius*2,-0.2,1.1]);
        ii = ii + 1;
    end

    figure_FontSize = 30;

    title(['\it{\fontname{Times new roman}E_{\rm{RMS}}} \rm{\fontname{Times new roman} = } ', num2str(roundn(Err,-3))], ...
        'FontSize', figure_FontSize);
    xlabel('\rm{\fontname{Times new roman}\it{x}-axis}', 'rotation', 22);
    ylabel('\rm{\fontname{Times new roman}\it{y}-axis}', 'rotation', -20);
    zlabel('\rm{\fontname{Times new roman}\it{ G(x,y)}}');
    set(get(gca,'XLabel'),'FontSize',figure_FontSize);
    set(get(gca,'YLabel'),'FontSize',figure_FontSize);
    set(get(gca,'ZLabel'),'FontSize',figure_FontSize);
    set(findobj('Fontsize',10),'fontsize', figure_FontSize-8);

    view([1,1,1]);
    hold off;

elseif Dims_Fourier == 3
    
    G_approx = imag(G_approx);
    G_exact = imag(G_exact);
    nMax1 = max(max(max(G_approx)));
    nMax2 = max(max(max(G_exact)));
    nMax = max([nMax1, nMax2]);
    nMin1 = min(min(min(G_approx)));
    nMin2 = min(min(min(G_exact)));
    nMin = min([nMin1, nMin2]);
    nG_approx = (G_approx - nMin) / (nMax - nMin);
    nG_exact = (G_exact - nMin) / (nMax - nMin);
    Err = sqrt(sum(sum(sum((nG_approx - nG_exact).^2))) / (length(x0)*length(y0)*length(z0)));        % RMSE

    [X,Y] = meshgrid(x0,y0);
    for jj = 1:length(z0)
        figure(jj); figure_FontSize = 30;
%         plot3(X, Y, real(G_approx(:, :, jj)), 'ro'); hold on;
%         plot3(X, Y, real(G_exact(:, :, jj)), 'b*');
        subplot(1,2,1);pclr1 = pcolor(real(X), real(Y), -abs(real(G_approx(:, :, jj))));
        set(pclr1, 'LineStyle','none');
        xlabel('\rm{\fontname{Times new roman}{\it{x}-axis} [mm]}');
        ylabel('\rm{\fontname{Times new roman}{\it{y}-axis} [mm]}');
        set(get(gca,'XLabel'),'FontSize',figure_FontSize);
        set(get(gca,'YLabel'),'FontSize',figure_FontSize);
        set(findobj('Fontsize',10),'fontsize', figure_FontSize-8);
        set(gca,'XTick', [-0.01, -0.005, 0, 0.005, 0.01]);
        set(gca,'XTickLabel', [10, 5, 0, 5, 10]);
        set(gca,'YTick', [-0.01, -0.005, 0, 0.005, 0.01]);
        set(gca,'YTickLabel', [10, 5, 0, 5, 10]);
        set(gca, 'FontName', 'Times new roman');
%         axis equal
%         colorbar;
        caxis([-100,100]*0.5);
        
        subplot(1,2,2);pclr2 = pcolor(real(X), real(Y), -abs(real(G_exact(:, :, jj))));
        set(pclr2, 'LineStyle','none');
        xlabel('\rm{\fontname{Times new roman}{\it{x}-axis} [mm]}');
        ylabel('\rm{\fontname{Times new roman}{\it{y}-axis} [mm]}');
        set(get(gca,'XLabel'),'FontSize',figure_FontSize);
        set(get(gca,'YLabel'),'FontSize',figure_FontSize);
        set(findobj('Fontsize',10),'fontsize', figure_FontSize-8);
        set(gca,'XTick', [-0.01, -0.005, 0, 0.005, 0.01]);
        set(gca,'XTickLabel', [10, 5, 0, 5, 10]);
        set(gca,'YTick', [-0.01, -0.005, 0, 0.005, 0.01]);
        set(gca,'YTickLabel', [10, 5, 0, 5, 10]);
        set(gca, 'FontName', 'Times new roman');
%         axis equal
%         colorbar;
        caxis([-100,100]*0.5);
    end
    
    % isosurface
%     p_iso = linspace(20, 40, 5);
%     for jj = 1:length(p_iso)
%         figure(jj);
%         [X,Y,Z] = meshgrid(x0,y0,z0);
%         isosurface(X, Y, Z, real(G_approx));
%     end

end

%%