function Hk_n = sHankel(n,x,type)
%%
% Spherical Hankel Function
% n--order
% x--independent variable (SHOULD be row vector)
% type--the first kind or the second kind
%%

if type == 1    % Eq.(6.57)
    Hk_n = sBessel(n,x,1) + sBessel(n,x,2)*1i;
elseif type == 2
    Hk_n = sBessel(n,x,1) - sBessel(n,x,2)*1i;
end

% check table can be found in Eq.(6.62)
