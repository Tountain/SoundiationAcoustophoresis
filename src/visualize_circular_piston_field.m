function visualize_circular_piston_field(p_0, freq, trans_radius, particle_radius)
%%
% This function is used to return a circular piston wave function (ref.
% Pierce@1999@Book Eqs. (5-5.1) and (5-5.4)) and Cheeke@2002@(6.32) 
% where x0, y0, z0 mean the relative position vector (i.e., \vec{Rr}) of
% the probe transducer ~Oq to the probe particle system Ol. Normally, for
% +z-axes direction, we set the probe transducer right beneath the probe
% particle (i.e., x_translation, y_translation) = (0, 0).   
%% solving

%trans_radius = 5 * 10^-3;       % radius of transducer is 5mm
%particle_radius = 2000 * 10^-6;
phi = 0;
%p_0 = 2000;
sound_speed = 340;
rho_0 = 1.224;
fluid_k = 2*pi * freq / sound_speed;

% the transducer output power
v0 = 1; % m/s ??
Pr_0 = - rho_0 * sound_speed * fluid_k * trans_radius^2 * v0 * 1i / 2;

range_r_coeff = 10;
grid_resolution = 100;
cos_beta = sqrt(1 - (1/(fluid_k * trans_radius))^2);

rr_visual = linspace(particle_radius, range_r_coeff*particle_radius, grid_resolution);   
theta_visual = linspace(0.01, pi, grid_resolution*2);

kk = 0;
for ii = 1:length(rr_visual)    % positive "phi" angle
    for jj = 1:length(theta_visual)
        kk = kk + 1;
        xx(kk) = rr_visual(ii)*cos(phi)*sin(theta_visual(jj));
        yy(kk) = rr_visual(ii)*sin(phi)*sin(theta_visual(jj));
        zz(kk) = rr_visual(ii)*cos(theta_visual(jj));
        p_i(kk) = Pr_0 * besselj(0, fluid_k * trans_radius * ...   % circuler piston wave 
            ones(length(phi),1) * sin(theta_visual(jj))) * ...         
            exp(1i * fluid_k * (rr_visual(ii))) / (rr_visual(ii));
        p_ii(kk) = Pr_0 * 2 * besselj(1, fluid_k * trans_radius * sin(theta_visual(jj))) / ...   % circuler piston wave 
            (fluid_k * trans_radius * sin(theta_visual(jj))) * ...         
            exp(1i * fluid_k * (rr_visual(ii))) / (rr_visual(ii));  
%         p_i(kk) = p_0 * besselj(0, rr_visual(ii)/trans_radius) * ...      
%             exp(1i * fluid_k * (rr_visual(ii)) * cos(theta_visual(jj)));
%         p_i(kk) = p_0 * besselj(0, rr_visual(ii) * sin(theta_visual(jj)) / trans_radius) * ...      % 0th Bessel beam        
%             exp(1i * fluid_k * (rr_visual(ii)) * cos(theta_visual(jj)) * cos_beta);
%         p_i(kk) = p_0 * besselj(0, 2.5 * rr_visual(ii) * sin(theta_visual(jj)) / trans_radius) * ...      % 0th Bessel beam (regional plain wave)        
%             (1-stepfun(rr_visual(ii) * sin(theta_visual(jj)) / trans_radius, 1)) * ...
%             exp(1i * fluid_k * (rr_visual(ii)) * cos(theta_visual(jj)));
%         p_i(kk) = p_0 * besselj(0, rr_visual(ii) * sin(theta_visual(jj)) / trans_radius) * ...       % 0th Bessel like beam (regional plain wave) 
%             exp(1i * fluid_k * (rr_visual(ii)) * cos(theta_visual(jj)));
    end
    fprintf('Finishing Calculation %f%%; \n', 100*ii/length(rr_visual)/2);
end
phi = phi + pi;     % negative "fixed_phi" angle
for ii = 1:length(rr_visual)    
    for jj = 1:length(theta_visual)
        kk = kk + 1;
        xx(kk) = rr_visual(ii)*cos(phi)*sin(theta_visual(jj));
        yy(kk) = rr_visual(ii)*sin(phi)*sin(theta_visual(jj));
        zz(kk) = rr_visual(ii)*cos(theta_visual(jj));
        p_i(kk) = Pr_0 * besselj(0, fluid_k * trans_radius * ...   % circuler piston wave 
            ones(length(phi),1) * sin(theta_visual(jj))) * ...         
            exp(1i * fluid_k * (rr_visual(ii))) / (rr_visual(ii));
        p_ii(kk) = Pr_0 * 2 * besselj(1, fluid_k * trans_radius * sin(theta_visual(jj))) / ...   % circuler piston wave 
            (fluid_k * trans_radius * sin(theta_visual(jj))) * ...         
            exp(1i * fluid_k * (rr_visual(ii))) / (rr_visual(ii)); 
%         p_i(kk) = p_0 * besselj(0, rr_visual(ii)/trans_radius) * ...      
%             exp(1i * fluid_k * (rr_visual(ii)) * cos(theta_visual(jj)));
%         p_i(kk) = p_0 * besselj(0, rr_visual(ii) * sin(theta_visual(jj)) / trans_radius) * ...      % 0th Bessel beam        
%             exp(1i * fluid_k * (rr_visual(ii)) * cos(theta_visual(jj)) * cos_beta);
%         p_i(kk) = p_0 * besselj(0, 2.5 * rr_visual(ii) * sin(theta_visual(jj)) / trans_radius) * ...      % 0th Bessel beam (regional plain wave)        
%             (1-stepfun(rr_visual(ii) * sin(theta_visual(jj)) / trans_radius, 1)) * ...
%             exp(1i * fluid_k * (rr_visual(ii)) * cos(theta_visual(jj)));
%         p_i(kk) = p_0 * besselj(0, rr_visual(ii) * sin(theta_visual(jj)) / trans_radius) * ...       % 0th Bessel like beam (regional plain wave)
%             exp(1i * fluid_k * (rr_visual(ii)) * cos(theta_visual(jj)));
    end
    fprintf('Finishing Calculation %f%%; \n', 100*ii/length(rr_visual)/2+50);
end

% p_i = p_i / max(abs(p_i));

%% visulizing

X = linspace(min(xx),max(xx),grid_resolution*2);
Z = linspace(min(zz),max(zz),grid_resolution*2);
[XX, ZZ, P1] = griddata(xx(1:10:end),zz(1:10:end),p_i(1:10:end),X',Z,'v4');
[~, ~, P2] = griddata(xx(1:10:end),zz(1:10:end),p_ii(1:10:end),X',Z,'v4');

% figure(1);
% pclr = pcolor(real(XX), real(ZZ), real(P1));
% set(pclr, 'LineStyle','none');
% %set(pclr, 'box', 'on');
% colorbar;
% caxis([-100,100]*0.5);

figure(2);
pclr = pcolor(real(XX), real(ZZ), real(P2));
set(pclr, 'LineStyle','none');
%set(pclr, 'box', 'on');
colorbar;
caxis([-100,100]*3);

%%