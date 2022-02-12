function bs_nmi_q = transform_beam_shape_matrix(equ_A_j, kr_jq, theta_jq, phi_jq, N, ii)
%%
% build the layer transform beam-shape coefficient of the i-th transducer.
% 
% NOTE: this database is created only for phase array system
% (strcmp(wave_type, 'phase_array_transducer') == 1).
%%

bs_nmi_q = zeros(size(equ_A_j));
for nn = 0:N
    indices_1 = nn + 1;
    parfor mm = -nn:nn          % 'mm' is sliced variable
        temp_bs_nmi_q = 0;
        for nu = 0:N
            for mu = -nu:nu
                [Snmvu_1, ~] = Snmvu_coeff(nu, mu, nn, mm, kr_jq, theta_jq, phi_jq);
                temp_bs_nmi_q = temp_bs_nmi_q + ...
                    equ_A_j(nu+1, nu+mu+1) * Snmvu_1;
            end
        end
        bs_nmi_q(indices_1, indices_1 + mm) = temp_bs_nmi_q;
    end
    fprintf('Transducer-%d: Transform Beam-Shape Coefficient Preparing %d%% \n', ...
        ii, round(100*nn/N));
end

%%
