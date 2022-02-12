function [Cn, EFn] = mapping_coeffs_NR(N, L, J, sample_number)
%%
% This function is used to: 
%     Firstly, determined the Fourier expansion coefficients 'En' and 'Fn',
% which connect the old and new polar angular coordinates using
% Newton-Raphson method.  
%     Secondly, based on the solved coefficients 'En' and 'Fn', we then
% evaluate the mapping coefficient 'Cn' used to map the old and new
% coordinate system.
% 
%
% Input variables:
%     N: the truncation number of Fourier expansion series for contour of
% length-wise surface in derivation of 'r_body' from 'theta_body'. 
%     L: the truncation number of Fourier expansion series for the old and
% new polar angular coordinate, the deviation of 'theta_body' from
% 'theta_map' on the surface contour of the length-wise slice. 
%     J: nonlinear equation numbers, J >= 2*L.
%     NN: the truncation number for mapping coefficient 'Cn'.  
%     sample_number: the sample pool to describe the locus of points of the
% length-wise slice of the irregular body.
%
% NOTE: for numerical convergence, the 'sample_number >= 10*max(N,L,NN)'
% condition should be satisied; 'sample_number' should be a even number.
% 
%
% Output variables: 
%     Cn: the mapping coefficients between old spherical coordinate system
% and new spherical coordinate system. 
%     En and Fn: the Fourier expansion coefficients for the derivation of
% 'theta_body' from 'theta_map' on the surface contour of the length-wise
% slice. 
%%

if J < 2*L
    error('Nonlinear equations number (J) should larger or equal to the number of the unknown variables (2*L)\n');
end

%% Obtain the Fourier expansion coefficients of the contour of length-wise slice surface

NN = J;

relaxation_factor = 1;  % introduce to control the deviation of 'EFn', and make sure the iterative adjustment does not jump too much and jump out the solution. 
% Rn_Max = 1000;
tolerance = 10^-7;      % the tolerance error
max_iteration = 1000;    % maximum iteration step
T = 2*pi;               % 'T' is set to the period of the body of polar angle (0 to 2*pi), actually, only [0,pi] is meaningful.

[avg_r, An, Bn, ~, ~] = Fourier_series_body_surface(N, sample_number);

Rn = 0.5 * (An + Bn*1i);

%% Mapping surface relationship in old and new coordinate systems
%% Part I: Solving the Fourier expansion coefficients 'En' and 'Fn'

% initial testing Fourier coefficients of 'En' and 'Fn'
En = 0.*ones(L, 1);
Fn = 0.*ones(L, 1);
EFn = [En; Fn];         % totally '2*L' unknown variables, thus '2*L' additional equations are required

% if length(EFn) < 2*L        % 'EFn' is the given initial Fourier coefficients, which can the convergent result for small 'L', and re-input for large 'L' cases.
%     EFn = [EFn; zeros(2*L-length(EFn), 1)];
% end

step = 1;
while step <= max_iteration
    
    int_j = zeros(2*L, 1);
    for jj = 1 : J-1
        int_j(jj) = integral_j(avg_r, Rn, T, EFn, jj+1, sample_number);             % 'jj' inputs to function "integral_j" should start from '2' to '2*L+1'
    end
    int_j;
    inv_J_EFn = Jacobian_matrix_nonlinear_system(avg_r, Rn, EFn, T, J, sample_number);     % Jacobian matrix
    
    EFn_t = EFn - relaxation_factor * inv_J_EFn * int_j;                          % Newton-Raphson method
%     EFn_t = EFn - real(relaxation_factor * inv_J_EFn * int_j);                        % Newton-Raphson method, make sure the EFn is real.
    
    % NOTE: if tolerance is too large, meaning the truncation number of
    % Fourier expansion series for the old and new polar angular coordinate
    % are too large (not satisfy: '2*L < sample_number'), which make the
    % system unstable, we need to reduce the 'L' (or increase the
    % 'sample_number') until finding a suitable 'L', in where 
    % '2*L < sample_number' (integration requirement).
    tol_cur = sum(abs(EFn_t - EFn));
    tol_int = sum(abs(int_j))/(J-1);
    if tol_cur > 10^6
        L = L - 1;
        EFn = [En(1:L); Fn(1:L)];
        fprintf('The Error = %f; and the truncation number ''L'' had been reduced to %d.\n', tol_cur, L);
        continue;
    else
        if tol_cur < tolerance
            EFn = EFn_t;
            fprintf('Iteration step: %d; Error = %f; integral = %f; Convergence!\n', step, tol_cur, tol_int);
            break;
        elseif step == max_iteration
            EFn = EFn_t;
            fprintf('Iteration step: %d; Error = %f; integral = %f; Reached the maximum iteration step!\n', step, tol_cur, tol_int);
        else
            EFn = EFn_t;
            fprintf('Iteration step: %d; Error = %f; integral = %f.\n', step, tol_cur, tol_int);
        end
        %EFn = real(EFn);        %
        step = step + 1;
    end
    
end

%% Part II: obtaining the mapping coefficient 'Cn'

int_j = zeros(NN+2, 1);
for nn = 1 : NN+2
    int_j(nn) = integral_j(avg_r, Rn, T, EFn, -nn+2, sample_number);            % 'jj' inputs to function "integral_j" should start from '-N' to '1'
end

Cn = int_j / T;             % Cn; n = -1, 0, 1, ..., N. half of them??

%% visualize the slice surface on the new coordinate system using 'Cn'

theta_body_Fourier = linspace(0, T, sample_number);
G_old1 = avg_r * exp(1i * 2*pi/T * theta_body_Fourier);
for nn = 1 : N
    G_old1 = G_old1 + conj(Rn(nn)) * exp(1i * (2*pi/T+nn) * theta_body_Fourier) + ...
        Rn(nn) * exp(1i * (2*pi/T-nn) * theta_body_Fourier);
end

u_map = 0;
theta_map_Fourier = linspace(0, T, sample_number*100);
G_new = Cn(1) * exp(1i * 2*pi/T * theta_map_Fourier);
for nn = 2 : length(Cn)
    G_new = G_new + Cn(nn) * exp(- 1i * 2*pi/T*(nn-2) * theta_map_Fourier);
end

theta_body_Fourier = theta_w(T, EFn, theta_map_Fourier, L);
G_old = avg_r * exp(1i * 2*pi/T * theta_body_Fourier);
for nn = 1 : N
    G_old = G_old + conj(Rn(nn)) * exp(1i * (2*pi/T+nn) * theta_body_Fourier) + ...
        Rn(nn) * exp(1i * (2*pi/T-nn) * theta_body_Fourier);
end

figure_FontSize = 18;

figure;
hold on; box on;
h1 = plot(real(G_old), imag(G_old), 'b', 'linewidth', 1);
h2 = plot(real(G_new), imag(G_new), 'k', 'linewidth', 2);
h3 = plot(real(G_old1), imag(G_old1), '--r', 'linewidth', 1);
title(['\rm{\fontname{Times new roman}Sample number=}', num2str(sample_number), ...
    '\rm{\fontname{Times new roman}; L=}', num2str(L), ...
    '\rm{\fontname{Times new roman}; J=}', num2str(J), '.'], ...
    'FontSize', figure_FontSize);
xlabel('\rm{\fontname{Times new roman}\it{z}}');
ylabel('\rm{\fontname{Times new roman}\it{x}}');
set(get(gca,'XLabel'),'FontSize',figure_FontSize);
set(get(gca,'YLabel'),'FontSize',figure_FontSize);
set(findobj('Fontsize',10),'fontsize', figure_FontSize-2);
axis equal;
hold off; 

%%