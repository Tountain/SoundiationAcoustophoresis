function visualize_polar_direction(sample_number)
%%
% This function is used to compare the polar direction between the model
% expansion method and its approximation wavefield (piston-like wavefield)
% of a single transduer. 
%%

parameters;

N = db_size_nn;

k = fluid_k;
phi_0 = 0;
phi_pi = pi;
theta_ap = asin(b/inter_dist);
theta = linspace (0.001, theta_ap, sample_number);

%% polar direction and graphs on r O theta-plane

% modal expansion method
for nn = 0:N
    for mm = -nn:nn
        db_bs_coeff(nn+1, mm+nn+1) = beam_shape_coeff(nn, mm);
    end
    fprintf('Beam-Shape and Scattering Coefficients Database Preparing %d%% \n',round(100*nn/N));
end

x = inter_dist * sin(theta);
z = inter_dist * cos(theta) - inter_dist;
rr = sqrt(x.^2 + z.^2);
theta_t = acos(z ./ (rr+eps));
p_i_pwe = zeros(sample_number, 2);
for ii = 1 : sample_number
    for nn = 0:N
        for mm = -nn:nn
            p_i_pwe(ii,1) = p_i_pwe(ii,1) + db_bs_coeff(nn+1, mm+nn+1) * ... 
                sBessel(nn,k*rr(ii),1) * sHarmonics(nn,mm,theta_t(ii),phi_0);
            p_i_pwe(ii,2) = p_i_pwe(ii,2) + db_bs_coeff(nn+1, mm+nn+1) * ... 
                sBessel(nn,k*rr(ii),1) * sHarmonics(nn,mm,theta_t(ii),phi_pi);
        end
    end
    fprintf('The Theta is %f.\n', theta(ii)*180/pi);
end

% piston-like wave function
p_ii = zeros(sample_number, 2);
for ii = 1 : sample_number
    p_ii(ii, 1) = Pr_0 * 2 * besselj(1, k * trans_radius * sin(theta(ii)) + eps) / ...   % circuler piston wave 
        (k * trans_radius * sin(theta(ii)) + eps) * ...         
        exp(1i * k * inter_dist) / (inter_dist + eps);  
end
p_ii(:, 2) = p_ii(:, 1);


% exp_rr = exp(1i * k * inter_dist);
p_amp_approx = p_i_pwe;                             % complex pressure amplitude of exact wave
v_amp_approx = - p_amp_approx / (fluid_rho * fluid_c);        % momentum equation; complex velocity amplitude of exact wave
p_amp_exact = p_ii;                           % complex pressure amplitude of approximation wave
v_amp_exact = - p_amp_exact / (fluid_rho * fluid_c);      % momentum equation; complex velocity amplitude of exact wave

% I_av = 0.5*Re(complex_amplitude_p * conj_complex_amplitude_v);
I_av_exact = 0.5 * real(p_amp_exact .* conj(v_amp_exact));   % time-average intensity of exact wave at position 'r_polar'
I_av_approx = 0.5 * real(p_amp_approx .* conj(v_amp_approx));% time-average intensity of approximation wave at position 'r_polar'
I_av_exact_norm = I_av_exact / I_av_exact(1,1);
I_av_approx_norm = I_av_approx / I_av_exact(1,1);

Err = sqrt(sum(sum((I_av_approx_norm - I_av_exact_norm).^2)) / (sample_number*2))        % RMSE
theta_ap
theta_ap = theta_ap*180/pi

% if b == 0.3*0.01
%     color1 = [1 0 0];
% elseif b == 0.4*0.01
%     color1 = [1 0.5 0.3];
% elseif b == 0.5*0.01
%     color1 = [1 0.27 0];
% elseif b == 0.6*0.01
%     color1 = [1 0.6 0.07];
% elseif b == 0.7*0.01
%     color1 = [1 0.84 0];
% elseif b == 0.8*0.01
%     color1 = [0.96 0.87 0.7];
% elseif b == 0.9*0.01
%     color1 = [0.64 0.58 0.5];
% else
%     color1 = [0.37 0.15 0.07];
% end

color1 = [1 0.27 0.];

figure(1);hold on;box on;
h1 = polar([-theta(end:-1:1) theta]'+pi/2, abs([I_av_exact_norm(end:-1:1,2);I_av_exact_norm(:,1)]));
hline1 = findobj(h1,'Type','line');
set(hline1, 'LineWidth', 3, 'color', [0.5 0.5 0.5]);

h2 = polar([-theta(end:-1:1) theta]'+pi/2, abs([I_av_approx_norm(end:-1:1,2);I_av_approx_norm(:,1)]), '--');
% polar(theta', real(p_amp_approx(:,1)), 'r');hold on;
% polar(theta', real(p_amp_exact(:,1)), 'k--');
hline2 = findobj(h2,'Type','line');
set(hline2, 'LineWidth', 3, 'color', color1);
axis([-0.15,0.15,0.8,1.1]);
% thetathicks(0:);

%%