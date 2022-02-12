function translated_beam_shape_coeff()
%%
%
%
%%

%%

% the transform beam-shape coefficient 3D matrix: bs_nmi_q (ref. "database_translation_coeffs.m")
bs_nmi_q = cell(1, transducer_number);
for ii = 1:transducer_number
    if ii == 1                  % except to probe transducer 'ii==1'
        bs_nmi_q{ii} = equ_A{ii};
    else
        bs_nmi_q{ii} = transform_beam_shape_matrix(equ_A{ii}, kr_iq(ii), theta_iq(ii), phi_iq(ii), N, ii);
    end
    fprintf('Transform Beam-Shape Coefficients Database Preparing %d%% \n', ...
        round(100*ii/transducer_number));
end


%%