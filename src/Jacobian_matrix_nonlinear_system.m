function inv_J_EFn = Jacobian_matrix_nonlinear_system(avg_r, Rn, EFn, T, J, sample_number)
%%
% This function is used to calculate the Jacobian matrix of the nonlinear
% system for further processing the Newton-Raphson method.
%
%
% Inputs:
%     avg_r and Rn: the Fourier expansion coefficients related to the slice
% surface.
%     En and Fn: the Fourier expansion coefficients for the derivation of
% 'theta_body' from 'theta_map' on the surface contour of the length-wise
% slice. 
%     T: the period of the the body of polar angle (0 to 2*pi), actually,
% only [0,pi] is meaningful.
%     J: nonlinear equation numbers, J >= 2*L.
%     
%
% Output:
%     inv_J_EFn: the inverse of Jacobian matrix of the nonlinear system
% regarding to variable 'En' and 'Fn'.
%%

L = length(EFn) / 2;        % truncation number of Fourier series in derivation of 'theta_body_map' from 'theta_map'
N = length(Rn);             % truncation number of Fourier series in derivation of 'r_body' from 'theta_body'

theta_map = linspace(0, T, sample_number);


%% The integration of nonlinear function

% Jacobian matrix: single thread
% J_EFn = zeros(J, 2*L);
% for jj = 1 : J-1
%     for ll = 1 : 2*L
%         
%         if ll <= L
%             surface_complex_partial = avg_r * 1i * cos(ll * theta_map) .* exp(theta_w(T, EFn, theta_map, L) * 1i);
%             for nn = 1 : N
%                 surface_complex_partial = surface_complex_partial + ...
%                     conj(Rn(nn)) * 1i*(nn+1)*cos(ll * theta_map) .* exp((2*pi/T*nn+1) * theta_w(T, EFn, theta_map, L) * 1i) + ...
%                     Rn(nn) * 1i*(-nn+1)*cos(ll * theta_map) .* exp((-2*pi/T*nn+1) * theta_w(T, EFn, theta_map, L) * 1i);
%             end
%         elseif ll > L
%             surface_complex_partial = avg_r * 1i * sin((ll-L) * theta_map) .* exp(theta_w(T, EFn, theta_map, L) * 1i);
%             for nn = 1 : N
%                 surface_complex_partial = surface_complex_partial + ...
%                     conj(Rn(nn)) * 1i*(nn+1)*sin((ll-L) * theta_map) .* exp((2*pi/T*nn+1) * theta_w(T, EFn, theta_map, L) * 1i) + ...
%                     Rn(nn) * 1i*(-nn+1)*sin((ll-L) * theta_map) .* exp((-2*pi/T*nn+1) * theta_w(T, EFn, theta_map, L) * 1i);
%             end
%         end
%         J_EFn(jj, ll) = trapz(theta_map, exp(- 1i * theta_map * (jj+1)) .* surface_complex_partial);
%     
%     end
% end


% Jacobian matrix: multiple threads
J_EFn = zeros(J-1, 2*L);
parfor jj = 1 : J-1
    
    J_EFn_line = zeros(1, 2*L);         % parapool grobal
    surface_complex_partial = 0;
    for ll = 1 : 2*L
        
        if ll <= L
            surface_complex_partial = avg_r * 1i * cos(ll * theta_map) .* exp(theta_w(T, EFn, theta_map, L) * 1i);
            for nn = 1 : N
                surface_complex_partial = surface_complex_partial + ...
                    conj(Rn(nn)) * 1i*(nn+1)*cos(ll * theta_map) .* exp((2*pi/T*nn+1) * theta_w(T, EFn, theta_map, L) * 1i) + ...
                    Rn(nn) * 1i*(-nn+1)*cos(ll * theta_map) .* exp((-2*pi/T*nn+1) * theta_w(T, EFn, theta_map, L) * 1i);
            end
        elseif ll > L
            surface_complex_partial = avg_r * 1i * sin((ll-L) * theta_map) .* exp(theta_w(T, EFn, theta_map, L) * 1i);
            for nn = 1 : N
                surface_complex_partial = surface_complex_partial + ...
                    conj(Rn(nn)) * 1i*(nn+1)*sin((ll-L) * theta_map) .* exp((2*pi/T*nn+1) * theta_w(T, EFn, theta_map, L) * 1i) + ...
                    Rn(nn) * 1i*(-nn+1)*sin((ll-L) * theta_map) .* exp((-2*pi/T*nn+1) * theta_w(T, EFn, theta_map, L) * 1i);
            end
        end
        J_EFn_line(ll) = trapz(theta_map, exp(- 1i * theta_map * (jj+1)) .* surface_complex_partial);
        
    end
    J_EFn(jj, :) = J_EFn_line;
    
end

% inveres of the Jacobian matrix
if J == 2*L
    inv_J_EFn = inv(J_EFn);
else
    inv_J_EFn = (J_EFn'*J_EFn) \ (J_EFn');
end

%%