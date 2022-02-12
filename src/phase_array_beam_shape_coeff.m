function [equivalent_bs_q_coeff, transducer] = phase_array_beam_shape_coeff(wave_type, N, db_bs_coeff, fluid_k, trans_radius, transducer_number, theta_rotation)
%%
% This function is used to establish the equivalent beam-shape coefficient
% of the zero-th Bessel beam phase array w.r.t the probe transducer.
%
% Firstly, loading the phase delay vector 'phase_delay_vector' and the
% beam-shape coefficient of transducer from database
% "database_beam_scattering_coeffs.m", and setting all transducers share
% the same beam-shape coefficient of transducer 'db_bs_coeff'.
%
% Secondly, applied the 'phase_delay_vector'
% [A_delay*exp(1i*fluid_k*phi_delay)] to all transducers, then we
% get different launch wave modes from different transducers (the amplitude
% and initial phase). And the new produced coefficients
% ('phase_delay_vector' applied on the 'db_bs_coeff') are the beam-shape
% coefficients of different transducers.
%
% NOTE: for the probe transducer, we consider the corresponding component
% of 'phase_delay_vector' as reference, so that it is fixed as
% 'A_delay=p_0' and 'phi_delay=0'. 
%
% Thirdly, using the translation addition theroem to transform all the
% information of source transducers coordinate systems into the probe
% transducer coordinate system. And the comprehensive beam-shape
% coefficient is defined as the 'equivalent beam-shape coefficient' of the
% phase array.
%
% Finally, replace the 'beam-shape coefficient (db_bs_coeff)' save in the
% database "database_beam_scattering_coeffs.m" with the 'equivalent
% beam-shape coefficient (equivalent_bs_q_coeff)'. 
% i.e., 'db_bs_coeff = equivalent_bs_coeff'; and save 'db_bs_coeff' into
%       the database "database_beam_scattering_coeffs.m".
%
%%

% parameters;     % for 'fluid_k', 'trans_radius', 'transducer_number'

% the phase array information
transducer = [0, 0, 0; 0.01, 0, 0; -0.01, 0, 0; 0, 0.01, 0; 0, -0.01, 0];
 
if transducer_number ~= size(transducer, 1)
    error('The transducer number is wrong.\n');
end

% ultrasound transducer parameters 'phase_delay_vector'
A_delay = [1 1 1 1 1];
phi_delay = [0 0 0 0 0];                                              % rad
phase_delay_vector = A_delay .* exp(1i * -phi_delay);              % loading the 'phase_delay_vector'


if strcmp(wave_type, 'phase_array_transducer2') == 1    % DO NOT continue the below addition theorem
    equivalent_bs_q_coeff = phase_delay_vector;         % 
    return;
end


% rotation of relative positions among the transducers
% here 'theta_rotation' should take a opposite direction
Rx=[1 0 0;
    0 cos(theta_rotation(1)) -sin(theta_rotation(1)); 
    0 sin(theta_rotation(1)) cos(theta_rotation(1))];
Ry=[cos(theta_rotation(2)) 0 sin(theta_rotation(2));
    0 1 0; 
    -sin(theta_rotation(2)) 0 cos(theta_rotation(2))];
Rz=[cos(theta_rotation(3)) -sin(theta_rotation(3)) 0; 
    sin(theta_rotation(3)) cos(theta_rotation(3)) 0; 
    0 0 1];
Rxyz = Rx * Ry * Rz;
for ii = 1:transducer_number
    if ii == 1
        continue;
    end
    transducer(ii, :) = (transducer(ii, :) - transducer(1, :)) * Rxyz;
end


%% 2nd: applied the 'phase_delay_vector' to all transducers (i.e., db_bs_coeff) for cell 'equ_A'

% equ_A = zeros((N+1)^2, transducer_number);       % number (N+1)^2 is total elements number of "anm,L"
% ii = 0;
% for nn = 0:N
%     for mm = -nn:nn
%         ii = ii + 1;
%         equ_A(ii, :) = db_bs_coeff(nn+1, nn+mm+1);
%     end
% end

equ_A = cell(1, transducer_number);
for ii = 1 : transducer_number
    equ_A{ii} = db_bs_coeff * phase_delay_vector(ii);
end


%% 3rd: obtaining the 'equivalent beam-shape coefficient' of the phase array 

% the relative position of source transducers r.w.t the probe transducer (ref. "database_translation_coeffs.m") 
r_iq = zeros(transducer_number, 1);     % "iq" means particle "i" toward probe transducer "q"
theta_iq = zeros(transducer_number, 1);
phi_iq = zeros(transducer_number, 1);
% number "ii == 1" represents the probe transducer 'q'
for ii = 1:transducer_number
    if ii == 1
        continue;
    end
    [r_iq(ii,1), theta_iq(ii,1), phi_iq(ii,1)] = ...
        coords_system_relative_positions_general(transducer(ii, :), transducer(1, :));
end
kr_iq = fluid_k * r_iq;

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

% the equivalent beam-shape coefficient: equivalent_bs_coeff
% summing up all 2D matrix layers in the 3D matrix 'bs_nmi_q' in element by
% element.
equivalent_bs_q_coeff = zeros(size(bs_nmi_q{1}));
for ii = 1:transducer_number
    equivalent_bs_q_coeff = equivalent_bs_q_coeff + bs_nmi_q{ii};
end

%%
