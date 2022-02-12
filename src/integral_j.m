function int_j = integral_j(avg_r, Rn, T, EFn, jj, sample_number)
%%
% This function is used to calculate the integration from 0 to 2*pi of the
% orthogonality relationship in terms of 'theta_map'.
%
%
% Inputs:
%     avg_r and Rn: the Fourier expansion coefficients related to the slice
% surface.
%     T: the period of the the body of polar angle (0 to 2*pi), actually,
% only [0,pi] is meaningful.
%     En and Fn: the Fourier expansion coefficients for the derivation of
% 'theta_body' from 'theta_map' on the surface contour of the length-wise
% slice. 
%     jj: the jj-th nonlinear equation of the system.
%
% Output:
%     int_j: the integration result.
%%

L = length(EFn) / 2;        % truncation number of Fourier series in derivation of 'theta_body_map' from 'theta_map'
N = length(Rn);             % truncation number of Fourier series in derivation of 'r_body' from 'theta_body'

theta_map = linspace(0, T, sample_number);

%% The integration of nonlinear function

surface_complex = avg_r * exp(theta_w(T, EFn, theta_map, L) * 1i);
for nn = 1 : N
    surface_complex = surface_complex + ...
        conj(Rn(nn)) * exp((2*pi/T*nn+1) * theta_w(T, EFn, theta_map, L) * 1i) + ...
        Rn(nn) * exp((-2*pi/T*nn+1) * theta_w(T, EFn, theta_map, L) * 1i);
end

elements = exp(- 1i * theta_map * jj) .* surface_complex;       % if delta(theta_map) is too large, then for jj>sample_number, the exp() jumps too much, and then the integration will invaild.
int_j = trapz(theta_map, elements);

%%