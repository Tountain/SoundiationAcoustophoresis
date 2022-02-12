function visualize_incident_wavefield()
%%
% this function is used to visualize the incident pressure data under
% different position of Cartesian coordinate.
%%

parameters_names;

%% partial wave expansion


%% numerical resolution and grid point information

rr_visual = linspace(particle_radius, range_r_coeff*particle_radius, grid_resolution);   
theta_visual = linspace(0, pi, grid_resolution*2);   
fixed_phi = 0;

global N db_bs_coeff db_s_coeff db_Bessel db_Hankel db_Harmonics db_b_s_coeff db_rs_in_coeff db_rs_out_coeff particle_involve layer_range
% these global variables use in "partial_wave_expansion.m"
[N, db_bs_coeff, db_s_coeff, db_Bessel, db_Hankel, db_Harmonics, ...
    db_b_s_coeff, db_rs_in_coeff, db_rs_out_coeff, particle_involve, layer_range] = ...
    loading_preparing_database(multi_particle, rr_visual, theta_visual, [fixed_phi,fixed_phi+pi]);
%[db_filename_BHH] = database_Bessel_Hankel_Harmonics(rr_visual, theta_visual, [fixed_phi,fixed_phi+pi]);

 
%% obtain the pressure data under different positions 

pi_v = zeros(1, 2*length(rr_visual) * length(theta_visual));        % these variables are used to save the 
xx = pi_v;
yy = pi_v;
zz = pi_v;
kk = 0;
for ii = 1:length(rr_visual)    % positive "fixed_phi" angle
    for jj = 1:length(theta_visual)
        
        kk = kk + 1;
        [pi_v(kk)] = partial_wave_expansion_incident(rr_visual(ii), ii, jj, 1);
        xx(kk) = rr_visual(ii)*cos(fixed_phi)*sin(theta_visual(jj));
        yy(kk) = rr_visual(ii)*sin(fixed_phi)*sin(theta_visual(jj));
        zz(kk) = rr_visual(ii)*cos(theta_visual(jj));
        
    end
    fprintf('Finishing Calculation %f%%; \n', 100*ii/length(rr_visual)/2);
end

fixed_phi = fixed_phi + pi;     % negative "fixed_phi" angle
for ii = 1:length(rr_visual)    
    for jj = 1:length(theta_visual)
        
        kk = kk + 1;
        [pi_v(kk)] = partial_wave_expansion_incident(rr_visual(ii), ii, jj, 2); 
        xx(kk) = rr_visual(ii)*cos(fixed_phi)*sin(theta_visual(jj));
        yy(kk) = rr_visual(ii)*sin(fixed_phi)*sin(theta_visual(jj));
        zz(kk) = rr_visual(ii)*cos(theta_visual(jj));
        
    end
    fprintf('Finishing Calculation %f%%; \n', 100*ii/length(rr_visual)/2+50);
end


%% visualize in xz-plain

X1 = linspace(min(xx),max(xx),grid_resolution);
Z1 = linspace(min(zz),max(zz),grid_resolution);
[XXI, ZZI, PI_V] = griddata(xx(1:10:end),zz(1:10:end),pi_v(1:10:end),X1',Z1,'v4');           

%contourf(XXT, ZZT, PT_V);
pclr = pcolor(real(XXI), real(ZZI), real(PI_V));
set(pclr, 'LineStyle','none');
%set(pclr, 'box', 'on');
colorbar;

title('\rm{\fontname{Times new roman}\rm{Computational coordinates (xOz Plane)}}');
xlabel('\rm{\fontname{Times new roman}\it{x}{\rm{-axis [m]}}}');
ylabel('\rm{\fontname{Times new roman}\it{z}{\rm{-axis [m]}}}');
set(gca, 'FontName', 'Times new roman');

delete([db_filename, '.mat']);
delete([db_filename_BHH, '.mat']);

%%


function [p_i_pwe] = partial_wave_expansion_incident(r, kr_ii, theta_jj, phi_kk)

%% single particle system

parameters;

global N db_bs_coeff db_Bessel db_Harmonics

% Incident wave expansion (pressure complex amplitude)     E.G.W@1999 Eq.(6.140)
p_i_pwe = 0;
for nn = 0:N
    for mm = -nn:nn

        p_i_pwe = p_i_pwe + db_bs_coeff(nn+1, mm+nn+1) * ... 
            db_Bessel(nn+1,kr_ii,1) * db_Harmonics(nn+1,nn+mm+1,theta_jj,phi_kk);

    end
end
