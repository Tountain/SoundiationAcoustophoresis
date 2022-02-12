function P_nm = sLegendre(n,m,theta)
%%
% Spherical Harmonics Function
% n--lower order
% m--upper order
% theta--polar angle, radian mesurement (SHOULD be row vector)  [0,pi]
%%

% 2021.5.6: making the function can work normally, if abs(m) > n.
if abs(m) > n
    P_nm = 0;
    return;
end

if m >= 0       % Eq.(6.29)
    
    P = legendre(n, cos(theta));  
    if length(theta) == 1
        P_nm = P(m+1);                       % P_nm(cos(theta)), m is positive
    else
        P_nm = P(m+1,:); 
    end
        
elseif m < 0    % Eq.(6.31)
    
    m = abs(m);                              % transfer the negative m to positive m
    coeff = factorial(n-m)/factorial(n+m);
    P = legendre(n, cos(theta));
    if length(theta) == 1
        P_nm = (-1)^(m) * coeff * P(m+1);    % P_n(-m)(cos(theta)), m is positive
    else
        P_nm = (-1)^(m) * coeff * P(m+1,:);  % P_n(-m)(cos(theta)), m is positive
    end
    
end

% check table can be found in Eq.(6.33), Eq.(6.34) and Eq.(6.35)