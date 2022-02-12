function rs_nmj = p_rs_nmj(db_ibs_coeff, db_s_coeff, kr_kL, theta_kL, phi_kL, N, kk, nn, mm)
%%
% this function is used to calculate the inner function of rescattering
% wave of multi-particle system of outer rescattering field. 
% Call by function "partial_wave_expansion.m"
%%

rs_nmj = 0;
for nu = 0:N
    for mu = -nu:nu
        rs_nmj = rs_nmj + db_ibs_coeff(nu+1, nu+mu+1, kk) * db_s_coeff(nu+1, kk) * ...
            Snmvu_coeff(nu, mu, nn, mm, kr_kL(kk), theta_kL(kk), phi_kL(kk), 1);
    end
end


%%