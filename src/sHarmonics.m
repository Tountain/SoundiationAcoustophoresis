function Y_nm = sHarmonics(n,m,theta,phi)
%%
% Spherical Harmonics Function
% n--lower order
% m--upper order
% theta--polar angle, radian mesurement (SHOULD be row vector)   [0,pi]
% phi--azimuthal angle, radian mesurement (SHOULD be row vector) [0,2*pi]
% Y_nm-- row:length(phi), column:length(theta)
%%

if m >= 0           % p.190 -- definition
    
    P_nm = sLegendre(n,m,theta);        % 1 * length(theta)
    exp_phi = (exp(1i * m * phi'));      % length(phi) * 1
    coeff = sqrt(((2*n+1)/(4*pi)) * (factorial(n-m)/factorial(n+m)));
    
    Y_nm = coeff * exp_phi * P_nm;      % Y_nm: length(phi) * length(theta), m is positive
    
elseif m < 0        % Eq.(6.44)
    
    m = abs(m);     % transfer the negative m to positive m
    P_nm = sLegendre(n,m,theta);        % 1 * length(theta)
    exp_phi = (exp(1i * m * phi'));         % length(phi) * 1
    coeff = sqrt(((2*n+1)/(4*pi)) * (factorial(n-m)/factorial(n+m)));
    Y_nm = coeff * exp_phi * P_nm;      % Y_n(-m): length(phi) * length(theta)
    conj_Y_nm = conj(Y_nm);
    
    Y_nm = (-1)^(m) * conj_Y_nm;        % Y_n(-m): length(phi) * length(theta), m is positive
    
end

% check table can be found in Eq.(6.51) and Eq.(6.52)
