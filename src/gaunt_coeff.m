function gaunt_v = gaunt_coeff(n, m, nu, mu, q)
%%
% this function is used to compute the Gaunt coefficient Based on
% definition Eq.(3.73) at P.A.MARTIN@2006@Multiple_Scattering.
%
% wigner 3-j symbol obtained by sub-function "wigner3j".
% W = wigner3j( J123, M123 ) 
%%

j123 = [n, nu, q];
m123 = [m, mu, -(m+mu)];

great_S = sqrt((2*n+1) * (2*nu+1) * (2*q+1) / (4*pi));
w1 = wigner3j(j123, [0, 0, 0]);
w2 = wigner3j(j123, m123);

gaunt_v = ((-1)^(m+mu)) * great_S * w1 * w2;

% check the correction of Gaunt coefficient can use Eq.(3.75)~(3.77) at
% P.A.MARTIN@2006@Multiple_Scattering
