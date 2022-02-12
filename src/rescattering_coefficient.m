function [db_rs_in_coeff, db_rs_out_coeff] = ...
                rescattering_coefficient(nn, mm, particle_involve_layer, ...
                particle_order, db_b_s_coeff, Snmvu_1, Snmvu_2, N)
%%
% this function calculate the rescattering coefficient "db_rs_in_coeff" and
% "db_rs_out_coeff" under different layers.
% each layer, we should have two rescattering coefficient (the inner and
% outer).
% these two coefficient are essential for calculate multi-paticle system's
% pressure field.
%%

particles_Cartesian_data;

db_rs_in_coeff = 0;
db_rs_out_coeff = 0;
for nu = 0:N
    for mu = -nu:nu
        % layer inner range
        for jj = particle_involve_layer+1:particle_number  
            db_rs_in_coeff = db_rs_in_coeff + db_b_s_coeff(nu+1, nu+mu+1, particle_order(jj)) * ...
                Snmvu_2(nu+1, nu+mu+1, nn+1, nn+mm+1, particle_order(jj), 1);
        end
        % layer outer range
        for jj = 2:particle_involve_layer   
            db_rs_out_coeff = db_rs_out_coeff + db_b_s_coeff(nu+1, nu+mu+1, particle_order(jj)) * ...
                Snmvu_1(nu+1, nu+mu+1, nn+1, nn+mm+1, particle_order(jj), 1);
        end
    end
end


%%
