function [Snmvu1, Snmvu2] = Snmvu_coeff(n, m, nu, mu, kb, theta, phi)
%%
% this function is used to compute the Translation coefficient "Snm_vu"
% Based on definition Eq.(3.80) for "type==1" and Eq.(3.87) for "type==2"
% at P.A.MARTIN@2006@Multiple_Scattering.
%
% Gaunt coefficient "gaunt value" and Spherical Harmonics "Ynm" are
% obtained by sub-function "gaunt_coeff" and "sHarmonics".
%%

coeff = 4*pi*1i^(nu-n);
q_lower = abs(n-nu);
q_upper = n+nu;

sum1 = 0;
sum2 = 0;

for qq = q_lower:q_upper
    if abs(mu-m)<=qq
        temp = conj(sHarmonics(qq,mu-m,theta,phi)) * gaunt_coeff(n,m,nu,-mu,qq);
        sum1 = sum1 + (1i)^qq * (-1)^m * sBessel(qq,kb,1) * temp;
        sum2 = sum2 + (1i)^qq * (-1)^m * sHankel(qq,kb,1) * temp;
    end
end
Snmvu1 = coeff * sum1;
Snmvu2 = coeff * sum2;


% % =======================OLD VERSION==========================
% 
% coeff = 4*pi*1i^(nu-n);
% q_lower = abs(n-nu);
% q_upper = n+nu;
% 
% sum = 0;
% if type == 1    
%     
%     for qq = q_lower:q_upper
%         if abs(mu-m)<=qq
%             sum = sum + (1i)^qq * (-1)^m * sBessel(qq,kb,1) * ...
%                 conj(sHarmonics(qq,mu-m,theta,phi)) * gaunt_coeff(n,m,nu,-mu,qq);
%         end
%     end
%     Snmvu = coeff * sum;
%     
% elseif type == 2
%     
%     for qq = q_lower:q_upper
%         if abs(mu-m)<=qq
%             sum = sum + (1i)^qq * (-1)^m * sHankel(qq,kb,1) * ...
%                 conj(sHarmonics(qq,mu-m,theta,phi)) * gaunt_coeff(n,m,nu,-mu,qq);
%         end
%     end
%     Snmvu = coeff * sum;
%     
% end
% % =======================OLD VERSION==========================


% check the correction of Translation coefficient "Snmvu" can use
% Eq.(3.82), Eq.(3.89), Eq.(3.93), Eq.(3.94), Eq.(3.95), Eq.(3.98),
% Eq.(3.99), Eq.(3.100) and Lamma 3.29 at
% P.A.MARTIN@2006@Multiple_Scattering
