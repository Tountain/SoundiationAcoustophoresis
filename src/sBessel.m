function Bsl= sBessel(n,x,type)
%%
% Spherical Bessel Function
% n--order
% x--independent variable (SHOULD be row vector)
% type--the first kind or the second kind
%%

if type == 1    % Eq.(6.57)
    Bsl = (pi/2./x).^0.5.*besselj(n+0.5, x);
elseif type == 2
    Bsl = (pi/2./x).^0.5.*bessely(n+0.5, x);
end

% check table can be found in Eq.(6.60) and Eq.(6.61)
