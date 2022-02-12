function [p_t_pwe, p_i_pwe, p_s_pwe, p_rs_pwe] = partial_wave_expansion(multi_particle, r, kr_ii, theta_jj, phi_kk)
%%
% this function is used to expansion the incident wave and scattering wave
% based on partial wave expansion method under the spherical coordinates.
%
% input data:
%    the Origin point locates in the center of the particle. and the sound
%    field spherical coordinate is related to this Origial point with 
%    "(r,theta,phi)==(kr_ii,theta_jj,phi_kk)"
%    ATTENTION: "(r,theta,phi)" all three input data support single number
%    only!
%
% output data:
%    p_t_pwe -- the partial wave expansion series of total complex 
%               amplitude pressure corresponding to point "(r,theta,phi)";
%    p_i_pwe -- the partial wave expansion series of incident complex 
%               amplitude pressure corresponding to point "(r,theta,phi)";
%    p_s_pwe -- the partial wave expansion series of scattering complex 
%               amplitude pressure corresponding to point "(r,theta,phi)".             
%%


%% Load matrix for "beam_shape_coeff" and "scattering_coefficient" 
%  (call the function "database_beam_scattering_coeffs" to build a connect 
%  matrix to avoid repeating call function "beam_shape_coeff" and
%  "scattering_coefficient") 
%% the number of expansion terms constant "N" also load by database.

% build the database, if already exist, then it will jump to next sentense.
% [db_filename, fluid_k, p_0] = database_beam_scattering_coeffs();      
% NOTE: the returning variables "fluid_k" and "p_0" of
% "database_beam_scattering_coeffs" are used to avoid repeating calling
% "parameters" and SO THAT saving time.

% build the database using interface
% [N, beam_shape, scattering] = loading_preparing_database(multi_particle);

global N db_bs_coeff db_s_coeff db_Bessel db_Hankel db_Harmonics db_b_s_coeff db_rs_in_coeff db_rs_out_coeff particle_involve layer_range


%% single particle system

if multi_particle == 0
    
    parameters;

    %% Incident wave expansion (pressure complex amplitude)     E.G.W@1999 Eq.(6.140)
    p_i_pwe = 0;
    for nn = 0:N
        for mm = -nn:nn

%             p_i_pwe = p_i_pwe + beam_shape_coeff(nn, mm) * ... 
%                 sBessel(nn,fluid_k*r,1) * sHarmonics(nn,mm,theta,phi);
%             p_i_pwe = p_i_pwe + beam_shape(nn+1, mm+nn+1) * ... 
%                 sBessel(nn,fluid_k*r,1) * sHarmonics(nn,mm,theta,phi);
%             p_i_pwe = p_i_pwe + db_bs_coeff(nn+1, mm+nn+1) * ... 
%                 sBessel(nn,fluid_k*r,1) * sHarmonics(nn,mm,theta,phi);
            p_i_pwe = p_i_pwe + db_bs_coeff(nn+1, mm+nn+1) * ... 
                db_Bessel(nn+1,kr_ii,1) * db_Harmonics(nn+1,nn+mm+1,theta_jj,phi_kk);
        
        end
    end
    %p_i_pwe = p_0 * p_i_pwe;

    %% scattering wave expansion (pressure complex amplitude)    E.G.W@1999 Eq.(6.92)
    p_s_pwe = 0;
    for nn = 0:N
        for mm = -nn:nn

            if irregular_body == 0          % single spherical object
%                 p_s_pwe = p_s_pwe + scattering_coefficient(nn) * ...
%                     beam_shape_coeff(nn, mm) * ... 
%                     sHankel(nn,fluid_k*r,1) * sHarmonics(nn,mm,theta,phi);
%                 p_s_pwe = p_s_pwe + scattering(nn+1) * ...
%                     beam_shape(nn+1, mm+nn+1) * ... 
%                     sHankel(nn,fluid_k*r,1) * sHarmonics(nn,mm,theta,phi);
%                 p_s_pwe = p_s_pwe + db_s_coeff(nn+1, 1) * ...
%                     db_bs_coeff(nn+1, mm+nn+1) * ... 
%                     sHankel(nn,fluid_k*r,1) * sHarmonics(nn,mm,theta,phi);
                p_s_pwe = p_s_pwe + db_s_coeff(nn+1, 1) * ...
                    db_bs_coeff(nn+1, mm+nn+1) * ... 
                    db_Hankel(nn+1,kr_ii,1) * db_Harmonics(nn+1,nn+mm+1,theta_jj,phi_kk);
            elseif irregular_body == 1      % single nonspherical object
                p_s_pwe = p_s_pwe + db_s_coeff(nn+1, nn+mm+1) * ...
                    db_bs_coeff(nn+1, mm+nn+1) * ... 
                    db_Hankel(nn+1,kr_ii,1) * db_Harmonics(nn+1,nn+mm+1,theta_jj,phi_kk);
            end

        end
    end
    %p_s_pwe = p_0 * p_s_pwe;
    p_rs_pwe = 0;   % No rescattering waves
    
end

%%


%% multi-particle system: scattering waves should consider other particles' rescattering effects  

if multi_particle ~= 0
    
    particles_Cartesian_data;
    
    p_i_pwe = 0;
    p_s_pwe = 0;
    p_rs_pwe = 0;
    
    %r_iL = r_ij(:,1);       % distances between particles and probe particle
    %r_layer = sort(r_iL);   % distances rank from closest to farest
    for layer = 1:length(layer_range)-1
        % innest range to sub-outest range
        if r >= layer_range(layer) && r < layer_range(layer+1)  
            
            for nn = 0:N
                for mm = -nn:nn
                    p_i_pwe = p_i_pwe + db_bs_coeff(nn+1, mm+nn+1) * ... 
                            db_Bessel(nn+1,kr_ii,1) * db_Harmonics(nn+1,nn+mm+1,theta_jj,phi_kk);
                    p_s_pwe = p_s_pwe + db_b_s_coeff(nn+1, mm+nn+1, 1) * ... 
                            db_Hankel(nn+1,kr_ii,1) * db_Harmonics(nn+1,nn+mm+1,theta_jj,phi_kk);
                    p_rs_pwe = p_rs_pwe + (db_rs_in_coeff(nn+1, nn+mm+1, layer) * db_Bessel(nn+1,kr_ii,1) +...
                            db_rs_out_coeff(nn+1, nn+mm+1, layer) * db_Hankel(nn+1,kr_ii,1)) * ...
                            db_Harmonics(nn+1,nn+mm+1,theta_jj,phi_kk);
                end
            end
            
        end
    end
    
    % outest range
    if r >= layer_range(end)
    
        for nn = 0:N
            for mm = -nn:nn
                p_i_pwe = p_i_pwe + db_bs_coeff(nn+1, mm+nn+1) * ... 
                        db_Bessel(nn+1,kr_ii,1) * db_Harmonics(nn+1,nn+mm+1,theta_jj,phi_kk);
                p_s_pwe = p_s_pwe + db_b_s_coeff(nn+1, mm+nn+1, 1) * ... 
                        db_Hankel(nn+1,kr_ii,1) * db_Harmonics(nn+1,nn+mm+1,theta_jj,phi_kk);
                p_rs_pwe = p_rs_pwe + (db_rs_in_coeff(nn+1, nn+mm+1, length(layer_range)) * db_Bessel(nn+1,kr_ii,1) +...
                        db_rs_out_coeff(nn+1, nn+mm+1, length(layer_range)) * db_Hankel(nn+1,kr_ii,1)) * ...
                        db_Harmonics(nn+1,nn+mm+1,theta_jj,phi_kk);
            end
        end
        
    end
    
    
% ======================== OLD VERSION ======================    
%     %% for innest sound field
%     if r <= r_avg  
%         p_in_pwe = 0;
%         for nn = 0:N
%             for mm = -nn:nn
% 
% %                 p_i_pwe = p_i_pwe + db_bs_coeff(nn+1, mm+nn+1) * ... 
% %                         sBessel(nn,fluid_k*r,1) * sHarmonics(nn,mm,theta,phi);
% %                 p_s_pwe = p_s_pwe + db_s_coeff(nn+1, 1) * ...
% %                         db_ibs_coeff(nn+1, mm+nn+1, 1) * ... 
% %                         sHankel(nn,fluid_k*r,1) * sHarmonics(nn,mm,theta,phi);
% %                 p_in_pwe = p_in_pwe + db_ibs_coeff(nn+1, mm+nn+1, 1) * ... 
% %                         sBessel(nn,fluid_k*r,1) * sHarmonics(nn,mm,theta,phi);
%                 p_i_pwe = p_i_pwe + db_bs_coeff(nn+1, mm+nn+1) * ... 
%                         db_Bessel(nn+1,kr_ii,1) * db_Harmonics(nn+1,nn+mm+1,theta_jj,phi_kk);
%                 p_s_pwe = p_s_pwe + db_s_coeff(nn+1, 1) * ...
%                         db_ibs_coeff(nn+1, mm+nn+1, 1) * ... 
%                         db_Hankel(nn+1,kr_ii,1) * db_Harmonics(nn+1,nn+mm+1,theta_jj,phi_kk);
%                 p_in_pwe = p_in_pwe + db_ibs_coeff(nn+1, mm+nn+1, 1) * ... 
%                         db_Bessel(nn+1,kr_ii,1) * db_Harmonics(nn+1,nn+mm+1,theta_jj,phi_kk);
%                 % for innest range, the "equivalent beam-shape coefficient"  
%             end
%         end
%         p_rs_pwe = p_in_pwe - p_i_pwe;
%     end
%     %% for outest sound field
%     if r > r_avg  
%         for nn = 0:N
%             for mm = -nn:nn
% 
% %                 p_i_pwe = p_i_pwe + db_bs_coeff(nn+1, mm+nn+1) * ... 
% %                         sBessel(nn,fluid_k*r,1) * sHarmonics(nn,mm,theta,phi);
% %                 p_s_pwe = p_s_pwe + db_s_coeff(nn+1, 1) * ...
% %                         db_ibs_coeff(nn+1, mm+nn+1, 1) * ... 
% %                         sHankel(nn,fluid_k*r,1) * sHarmonics(nn,mm,theta,phi);
% %                 p_rs_pwe = p_rs_pwe + db_rs_coeff(nn+1, nn+mm+1) * ...
% %                         sHankel(nn,fluid_k*r,1) * sHarmonics(nn,mm,theta,phi);
%                 p_i_pwe = p_i_pwe + db_bs_coeff(nn+1, mm+nn+1) * ... 
%                         db_Bessel(nn+1,kr_ii,1) * db_Harmonics(nn+1,nn+mm+1,theta_jj,phi_kk);
%                 p_s_pwe = p_s_pwe + db_s_coeff(nn+1, 1) * ...
%                         db_ibs_coeff(nn+1, mm+nn+1, 1) * ... 
%                         db_Hankel(nn+1,kr_ii,1) * db_Harmonics(nn+1,nn+mm+1,theta_jj,phi_kk);
%                 p_rs_pwe = p_rs_pwe + db_rs_coeff(nn+1, nn+mm+1) * ...
%                         db_Hankel(nn+1,kr_ii,1) * db_Harmonics(nn+1,nn+mm+1,theta_jj,phi_kk);
% 
%             end
%         end
%     end
% %     p_i_pwe = p_0 * p_i_pwe;   % DO NOT REPEAT CONSIDER "p_0"(beam-shape coefficient already included)
% %     p_s_pwe = p_0 * p_s_pwe;
% %     p_rs_pwe = p_0 * p_rs_pwe;
% ======================== OLD VERSION ======================   

end

%%

%% total wave expansion (pressure complex amplitude)

p_t_pwe = p_i_pwe + p_s_pwe + p_rs_pwe;


%%