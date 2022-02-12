function [db_filename] = database_beam_scattering_coeffs()
%%
% build the database for "beam_shape_coeffcient" and
% "scattering_coefficient" to avoid repeating call function "beam_shape_coeff" and
% "scattering_coefficient", which will cost tons of time.
% 
% database will automatically be created if current folder has not a
% corresponding database. (check the database file exists or not by the
% filename)
% database maximum size: db_size_nn = 30 (seeing "parameters.m");
%%

parameters_names;

if exist([db_filename, '.mat']) ~= 0        % if already exist the database, 
    return;                                 % then do not create again for saving time.
end

db_bs_coeff = zeros(db_size_nn+1, 2*db_size_nn+1);
if multi_particle == 0      % single particle system
    if irregular_body == 0      % single spherical object
        db_s_coeff = zeros(db_size_nn+1, 1);
    elseif irregular_body == 1  % single nonspherical object
        db_s_coeff = zeros(db_size_nn+1, 2*db_size_nn+1);
    end
else                        % multi-particle system
    particles_Cartesian_data;
    db_s_coeff = zeros(db_size_nn+1, particle_number);
end


%% determine the Truncation number of expansion terms: N

N = 0;          % Truncation number (based on "sn") for single-particle system
if multi_particle == 0
%     error = 1;
    % "cut-off theory" based on G.T.S@2011@IEEE p.299
%     while (error > 0.01) && (N < db_size_nn)
%         if irregular_body == 0
%             error = abs(db_s_coeff(N+1+1) / db_s_coeff(N+0+1));      
%             N = N + 1;
%         elseif irregular_body == 1
%             error = abs(sum(db_s_coeff(N+1+1, :)) / sum(db_s_coeff(N+0+1, :)));
%             N = N + 1;
%         end
%     end
    % get the small "N" from "N=ka+6" and "cut-off theory"
%     if N > round(fluid_k*particle_radius + 8)  
%         N = round(fluid_k*particle_radius + 8);
%     end
    N = db_size_nn;            % a larger "N", a "wider" plain wave
else
    % maximum "database size for beam-shape coeff, espacially plain wave,
    % as small database size will lose the plain wave information in
    % farfield, but for multi-particle system, we need farfield
    % information.
    N = db_size_nn;            % based on J.H@2016@IEEE, NOTE: the particle radius "ka<4"     
end

% avoid too many expansion terms, "db_size" is the maximum size of
% "beam_shape_coefficient" and "scattering_coefficient" in the database.
% if N > db_size_nn       
%     N = db_size_nn;
% end


%% database for "beam_shape_coeffcient": db_bs_coeff
% mapping relation of "db_bs_coeff" and "beam_shape_coeffcient"
%   1.first subscript for "scattering_coefficient" is "nn" 
%     but "nn+1" for "db_s_coeff"
%   2.second subscript for "scattering_coefficient" is "mm" 
%     but "mm+nn+1" for "db_s_coeff"
%
% particularly, beam_shape_coeffcient(NN, MM)
% then, db_bs_coeff(NN+1, NN+MM+1)
% where NN = nn + num_C, NN = mm + num_D
%
% NOTE: after finishing position translation(num_C or num_D), the new
%       position NN and MM should meet "abs(MM) <= abs(NN)", because for
%       "beam_shape_coeffcient", we require "abs(mm) <= abs(nn)".

for nn = 0:N
    for mm = -nn:nn
        
        db_bs_coeff(nn+1, mm+nn+1) = beam_shape_coeff(nn, mm);

    end
    fprintf('Beam-Shape and Scattering Coefficients Database Preparing %d%% \n',round(100*nn/N));
end

%% database for "equivalent beam_shape_coeffcient": equivalent_bs_coeff -> db_bs_coeff

if strcmp(wave_type, 'phase_array_transducer') == 1 
    [equivalent_bs_q_coeff, transducer] = phase_array_beam_shape_coeff(wave_type, N, db_bs_coeff, fluid_k, trans_radius, transducer_number, theta_rotation);
    db_bs_coeff = equivalent_bs_q_coeff;      % equivalent_bs_coeff -> db_bs_coeff
elseif strcmp(wave_type, 'phase_array_transducer2') == 1 
    [~, transducer] = phase_array_beam_shape_coeff(wave_type, N, db_bs_coeff, fluid_k, trans_radius, transducer_number, theta_rotation);
else
    transducer = [0,0,0];
end
fprintf('Equivalent Beam-Shape Coefficient is prepared well! \n');


%% database for "scattering_coefficient": db_s_coeff
% mapping relation of "db_s_coeff" and "scattering_coefficient"
%   subscript for "scattering_coefficient" is "nn" 
%   but "nn+1" for "db_s_coeff" 
%
% particularly, scattering_coefficient(NN)
% then, db_s_coeff(NN+1)
% where NN = nn + num_C

if multi_particle == 0      % single particle system
    
    if irregular_body == 0
        for nn = 0:db_size_nn   % single spherical object
            db_s_coeff(nn+1) = scattering_coefficient(nn); 
        end
    elseif irregular_body == 1  % single nonspherical object
        snm_ib = scattering_coeff_irregular_body(BC, db_bs_coeff);
        for nn = 0:db_size_nn
            for mm = -nn:nn
                db_s_coeff(nn+1, mm+nn+1) = snm_ib(nn+1, mm+nn+1);
            end
        end
    end
    
else                        % multi-particle system
    
    for nn = 0:db_size_nn   % multiple spherical objects
        db_s_coeff(nn+1, :) = scattering_coefficient(nn); 
    end
    
end

%% save these databases by meaningful name

% r = max(max(abs(db_s_coeff)))/max(max(abs(db_bs_coeff)))
% f_ka 

save([db_filename, '.mat'], 'N', 'db_bs_coeff', 'db_s_coeff', 'transducer');


%%