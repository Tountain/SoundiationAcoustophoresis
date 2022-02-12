function [Frad_x,Frad_y,Frad_z,transducer] = radiation_force_based_Analyses()
%%
% this function is used to obtain the Radiation Force from the pressure
% distribution surrounding the particle.
%
% NOTE: the theory are based on lots of previous deriviation and some
% personal understandings. 
% In order to finally obtain the Radiation Force, major reference materials
% include (listed but not limited):    
%    P.J.W@1951(JASA); P.J.W@1957(JASA); L.P.Gorkov@1962; K.Yosioka@1955;
%    E.G.W@1999(book); J.D@1987(PRL); G.T.S@2011(IEEE); G.T.S@2013(IEEE); 
%    G.T.S@2011(JASA); G.T.S@2014(JASA);
%    (ABOVE refs. based on Inviscous Fluid!)
%%

testing_type = 4;

% elapse_time1 = datestr(now,'yyyy-mm-dd HH:MM:SS');
% disp(['Started Time: ',elapse_time1]);

global absp_factor
absp_factor = 0;
% parameters_names;
parameters;

% revised the attenuation (absorption) factor for mixture medium
if strcmp(fluid, 'mixture') == 1
    absorption_factor();        % "multi_particle_force_torque" should change to "0" for "attenuation coefficient".
    parameters;                 % revised parameters
end

% multi-particle system can not use this function 
% if multi_particle ~= 0
%     error('This function (for radiation forces) ONLY works for single particle system!');
% elseif testing_type ~= 4    % ONLY "test-4" can use for other incident wave directions 
if testing_type ~= 4
    if strcmp(direction, 'Z') ~= 1 || strcmp(symbol, 'positive') ~= 1
        error('This function (for radiation forces) ONLY works for Z-positive incident wave!');
    end
end

% build the database, if already exist, then it will jump to next sentense.
if multi_particle == 0
    [db_filename] = database_beam_scattering_coeffs();      
    load([db_filename, '.mat'], 'N', 'db_bs_coeff', 'db_s_coeff', 'transducer');
else
    [db_filename] = database_beam_scattering_coeffs();      
    load([db_filename, '.mat'], 'N', 'db_bs_coeff', 'db_s_coeff', 'transducer');
    [db_filename_ts] = database_translation_coeffs(N);
    load([db_filename_ts, '.mat'], 'Snmvu_1', 'Snmvu_2', 'r_ij');
    [db_filename_ibs] = database_interaction_beam_rescattering_coeffs(N, db_s_coeff, Snmvu_1, Snmvu_2, r_ij); 
    load([db_filename_ibs, '.mat'], 'db_ibs_coeff');
    
    db_bs_coeff = db_ibs_coeff(:,:,1);      % "1" means probe particle
    db_s_coeff = db_s_coeff(:,1);
end


%% obtain the force constant "force coeff", "AA", "BB", "CC", "DD", "EE" and "FF"

coeff_force = - (1)^2 / (2 * fluid_rho * fluid_c^2);

AA = @(nn,mm) (-1) * sqrt(((nn+mm-1) * (nn+mm)) / ((2*nn-1) * (2*nn+1)));
BB = @(nn,mm) sqrt(((nn-mm+2) * (nn-mm+1)) / ((2*nn+1) * (2*nn+3)));
CC = @(nn,mm) sqrt(((nn-mm-1) * (nn-mm)) / ((2*nn-1) * (2*nn+1)));
DD = @(nn,mm) (-1) * sqrt(((nn+mm+2) * (nn+mm+1)) / ((2*nn+1) * (2*nn+3)));
EE = @(nn,mm) sqrt(((nn-mm) * (nn+mm)) / ((2*nn-1) * (2*nn+1)));
FF = @(nn,mm) sqrt(((nn-mm+1) * (nn+mm+1)) / ((2*nn+1) * (2*nn+3)));


%% obtain the radiation force along X-axis

sum_1 = 0;
sum_2 = 0;
sum_3 = 0;
sum_4 = 0;
for nn = 0:(N-1)    % (n+1), (m+1)
    for mm = -nn:(nn)
        if (testing_type == 1)      % deriviation of myself for farfield approximately(a) 
            middle_var = -(1i^nn)*sin(fluid_k*R-nn*pi/2)*exp(-1i*fluid_k*R) + db_s_coeff(nn+1)*1i;       
        end
        if (testing_type == 2)      % deriviation of myself for farfield approximately(b) 
            middle_var = -(1i^nn)*sin(fluid_k*R-nn*pi/2)*exp(-1i*fluid_k*R) + db_s_coeff(nn+1)*1i;
        end
        if (testing_type == 3)      % deriviation of myself for exactly solution
            middle_var = (fluid_k*R)^2 * (sBessel(nn,fluid_k*R,1)*conj(sHankel(nn+1,fluid_k*R,1)) ...   
                + db_s_coeff(nn+1)*sHankel(nn,fluid_k*R,1)*conj(sHankel(nn+1,fluid_k*R,1)));
        end
        if (testing_type == 4)      % deriviation from G.T.S@2011@JASA
            if irregular_body == 0
                middle_var = 1 + db_s_coeff(nn+1);
            elseif irregular_body == 1
                middle_var = 1 + db_s_coeff(nn+1, nn+mm+1);
            end
        end
        
        n_plus_1 = nn + 1;      % adjust the projecting relation of database and beam-shape and scattering coefficents
        m_plus_1 = mm + 1;
        if abs(m_plus_1) <= n_plus_1    % make sure the adjusting position within the database
            if irregular_body == 0
                part_AA = conj(db_bs_coeff(n_plus_1+1,m_plus_1+n_plus_1+1)) * conj(db_s_coeff(n_plus_1+1)) * AA(nn+1,mm+1);
                sum_1 = sum_1 + db_bs_coeff(nn+1,mm+nn+1) * middle_var * part_AA;
            elseif irregular_body == 1
                part_AA = conj(db_bs_coeff(n_plus_1+1,m_plus_1+n_plus_1+1)) * conj(db_s_coeff(n_plus_1+1,m_plus_1+n_plus_1+1)) * AA(nn+1,mm+1);
                sum_1 = sum_1 + db_bs_coeff(nn+1,mm+nn+1) * middle_var * part_AA;
            end
        end
    
    end
end
for nn = 1:(N)      % (n-1), (m+1)
    for mm = -nn:(nn)
        if (testing_type == 1)      % deriviation of myself for farfield approximately(a) 
            middle_var = -(1i^nn)*sin(fluid_k*R-nn*pi/2)*exp(-1i*fluid_k*R) + db_s_coeff(nn+1)*1i;       
        end
        if (testing_type == 2)      % deriviation of myself for farfield approximately(b)
            middle_var = -(1i^nn)*sin(fluid_k*R-nn*pi/2)*exp(-1i*fluid_k*R) + db_s_coeff(nn+1)*1i;       
        end
        if (testing_type == 3)      % deriviation of myself for exactly solution
            middle_var = (fluid_k*R)^2 * (sBessel(nn,fluid_k*R,1)*conj(sHankel(nn-1,fluid_k*R,1)) ...   
                + db_s_coeff(nn+1)*sHankel(nn,fluid_k*R,1)*conj(sHankel(nn-1,fluid_k*R,1)));     
        end
        if (testing_type == 4)      % deriviation from G.T.S@2011@JASA
            if irregular_body == 0
                middle_var = 1 + db_s_coeff(nn+1);                                                        
            elseif irregular_body == 1
                middle_var = 1 + db_s_coeff(nn+1, nn+mm+1);
            end
        end
        
        n_minus_1 = nn - 1;
        m_plus_1 = mm + 1;
        if abs(m_plus_1) <= n_minus_1
            if irregular_body == 0
                part_BB = conj(db_bs_coeff(n_minus_1+1,m_plus_1+n_minus_1+1)) * conj(db_s_coeff(n_minus_1+1)) * BB(nn-1,mm+1);
                sum_2 = sum_2 + db_bs_coeff(nn+1,mm+nn+1) * middle_var * part_BB;
            elseif irregular_body == 1
                part_BB = conj(db_bs_coeff(n_minus_1+1,m_plus_1+n_minus_1+1)) * conj(db_s_coeff(n_minus_1+1,m_plus_1+n_minus_1+1)) * BB(nn-1,mm+1);
                sum_2 = sum_2 + db_bs_coeff(nn+1,mm+nn+1) * middle_var * part_BB;
            end
        end
    
    end
end
for nn = 0:(N-1)    % (n+1), (m-1)
    for mm = (-nn):(nn)
        if (testing_type == 1)      % deriviation of myself for farfield approximately(a) 
            middle_var = -(1i^nn)*sin(fluid_k*R-nn*pi/2)*exp(-1i*fluid_k*R) + db_s_coeff(nn+1)*1i;
        end
        if (testing_type == 2)      % deriviation of myself for farfield approximately(b)
            middle_var = -(1i^nn)*sin(fluid_k*R-nn*pi/2)*exp(-1i*fluid_k*R) + db_s_coeff(nn+1)*1i;
        end
        if (testing_type == 3)      % deriviation of myself for exactly solution
            middle_var = (fluid_k*R)^2 * (sBessel(nn,fluid_k*R,1)*conj(sHankel(nn+1,fluid_k*R,1)) ...   
                + db_s_coeff(nn+1)*sHankel(nn,fluid_k*R,1)*conj(sHankel(nn+1,fluid_k*R,1)));
        end
        if (testing_type == 4)      % deriviation from G.T.S@2011@JASA
            if irregular_body == 0
                middle_var = 1 + db_s_coeff(nn+1);                                                        
            elseif irregular_body == 1
                middle_var = 1 + db_s_coeff(nn+1, nn+mm+1);
            end
        end
        
        n_plus_1 = nn + 1;
        m_minus_1 = mm - 1;
        if abs(m_minus_1) <= n_plus_1
            if irregular_body == 0
                part_CC = conj(db_bs_coeff(n_plus_1+1,m_minus_1+n_plus_1+1)) * conj(db_s_coeff(n_plus_1+1)) * CC(nn+1,mm-1);
                sum_3 = sum_3 + db_bs_coeff(nn+1,mm+nn+1) * middle_var * part_CC;
            elseif irregular_body == 1
                part_CC = conj(db_bs_coeff(n_plus_1+1,m_minus_1+n_plus_1+1)) * conj(db_s_coeff(n_plus_1+1,m_minus_1+n_plus_1+1)) * CC(nn+1,mm-1);
                sum_3 = sum_3 + db_bs_coeff(nn+1,mm+nn+1) * middle_var * part_CC;
            end
        end
    
    end
end
for nn = 1:(N)      % (n-1), (m-1)
    for mm = (-nn):(nn)
        if (testing_type == 1)      % deriviation of myself for farfield approximately(a) 
            middle_var = -(1i^nn)*sin(fluid_k*R-nn*pi/2)*exp(-1i*fluid_k*R) + db_s_coeff(nn+1)*1i;
        end
        if (testing_type == 2)      % deriviation of myself for farfield approximately(b)
            middle_var = -(1i^nn)*sin(fluid_k*R-nn*pi/2)*exp(-1i*fluid_k*R) + db_s_coeff(nn+1)*1i;    
        end
        if (testing_type == 3)      % deriviation of myself for exactly solution
            middle_var = (fluid_k*R)^2 * (sBessel(nn,fluid_k*R,1)*conj(sHankel(nn-1,fluid_k*R,1)) ...   
            + db_s_coeff(nn+1)*sHankel(nn,fluid_k*R,1)*conj(sHankel(nn-1,fluid_k*R,1)));     
        end
        if (testing_type == 4)      % deriviation from G.T.S@2011@JASA
            if irregular_body == 0
                middle_var = 1 + db_s_coeff(nn+1);                                                        
            elseif irregular_body == 1
                middle_var = 1 + db_s_coeff(nn+1, nn+mm+1);
            end
        end
        
        n_minus_1 = nn - 1;
        m_minus_1 = mm - 1;
        if abs(m_minus_1) <= n_minus_1
            if irregular_body == 0
                part_DD = conj(db_bs_coeff(n_minus_1+1,n_minus_1+m_minus_1+1)) * conj(db_s_coeff(n_minus_1+1)) * DD(nn-1,mm-1);
                sum_4 = sum_4 + db_bs_coeff(nn+1,mm+nn+1) * middle_var * part_DD;
            elseif irregular_body == 1
                part_DD = conj(db_bs_coeff(n_minus_1+1,n_minus_1+m_minus_1+1)) * conj(db_s_coeff(n_minus_1+1,n_minus_1+m_minus_1+1)) * DD(nn-1,mm-1);
                sum_4 = sum_4 + db_bs_coeff(nn+1,mm+nn+1) * middle_var * part_DD;
            end
        end
    
    end
end

if (testing_type == 1)      % deriviation of myself for farfield approximately(a)
    Frad_x = coeff_force * 1/2 * real((sum_1 - sum_2 + sum_3 - sum_4) / (fluid_k^2));   
end
if (testing_type == 2)      % deriviation of myself for farfield approximately(b)
    Frad_x = coeff_force * 1/2 * real((sum_1 - sum_2 + sum_3 - sum_4) / (fluid_k^2));   
end
if (testing_type == 3)      % deriviation of myself for exactly solution
    Frad_x = coeff_force * 1/2 * real((sum_1 + sum_2 + sum_3 + sum_4) / (fluid_k^2));     
end
if (testing_type == 4)      % deriviation from G.T.S@2011@JASA
    Frad_x = coeff_force * 1/2 * real(1i * (sum_1 - sum_2 + sum_3 - sum_4) / (fluid_k^2));     
end


%% obtain the radiation force along Y-axis

if (testing_type == 1)      % deriviation of myself for farfield approximately(a)
    Frad_y = coeff_force * 1/2 * real((-1i)*(sum_1 - sum_2 - sum_3 + sum_4) / (fluid_k^2));
end
if (testing_type == 2)      % deriviation of myself for farfield approximately(b)
    Frad_y = coeff_force * 1/2 * real((-1i)*(sum_1 - sum_2 - sum_3 + sum_4) / (fluid_k^2));
end
if (testing_type == 3)      % deriviation of myself for exactly solution
    Frad_y = coeff_force * 1/2 * real((-1i)*(sum_1 + sum_2 - sum_3 - sum_4) / (fluid_k^2));
end
if (testing_type == 4)      % deriviation from G.T.S@2011@JASA
%     Frad_y = - coeff_force * 1/2 * real((sum_1 - sum_2 - sum_3 + sum_4) / (fluid_k^2)); 
    Frad_y = coeff_force * 1/2 * real((sum_1 - sum_2 - sum_3 + sum_4) / (fluid_k^2)); 
end


%% obtain the radiation force along Z-axis (parallel with incident wave)

sum_1 = 0;
sum_2 = 0;
sum_3 = 0;
for nn = 0:(N-1)     % (n+1), m
    for mm = -nn:(nn)
        if (testing_type == 1)      % deriviation of myself for farfield approximately(a)
            middle_var = db_s_coeff(nn+1)*1i;
        end
        if (testing_type == 2)      % deriviation of myself for farfield approximately(b)
            middle_var = -(1i^nn)*sin(fluid_k*R-nn*pi/2)*exp(-1i*fluid_k*R) + db_s_coeff(nn+1)*1i;   
        end
        if (testing_type == 3)      % deriviation of myself for exactly solution
            middle_var = (fluid_k*R)^2 * (sBessel(nn,fluid_k*R,1)*conj(sHankel(nn+1,fluid_k*R,1)) ...   
                + db_s_coeff(nn+1)*sHankel(nn,fluid_k*R,1)*conj(sHankel(nn+1,fluid_k*R,1)));     
        end
        if (testing_type == 4)      % deriviation from G.T.S@2011@JASA
            if irregular_body == 0
                middle_var = 1 + db_s_coeff(nn+1);
            elseif irregular_body == 1
                middle_var = 1 + db_s_coeff(nn+1, nn+mm+1);
            end
        end
            
        n_plus_1 = nn + 1;
        if abs(mm) <= n_plus_1
            if irregular_body == 0
                part_EE = conj(db_bs_coeff(n_plus_1+1,mm+n_plus_1+1)) * conj(db_s_coeff(n_plus_1+1)) * EE(nn+1,mm);
                sum_1 = sum_1 + db_bs_coeff(nn+1,mm+nn+1) * middle_var * part_EE;
            elseif irregular_body == 1
                part_EE = conj(db_bs_coeff(n_plus_1+1,mm+n_plus_1+1)) * conj(db_s_coeff(n_plus_1+1,mm+n_plus_1+1)) * EE(nn+1,mm);
                sum_1 = sum_1 + db_bs_coeff(nn+1,mm+nn+1) * middle_var * part_EE;
            end
        end
        
    end
end
for nn = 1:(N)      % (n-1), m
    for mm = -(nn):(nn)
        if (testing_type == 1)      % deriviation of myself for farfield approximately(a)
            middle_var = db_s_coeff(nn+1)*1i;
        end
        if (testing_type == 2)      % deriviation of myself for farfield approximately(b)
            middle_var = -(1i^nn)*sin(fluid_k*R-nn*pi/2)*exp(-1i*fluid_k*R) + db_s_coeff(nn+1)*1i;   
        end
        if (testing_type == 3)      % deriviation of myself for exactly solution
            middle_var = (fluid_k*R)^2 * (sBessel(nn,fluid_k*R,1)*conj(sHankel(nn-1,fluid_k*R,1)) ...   
                + db_s_coeff(nn+1)*sHankel(nn,fluid_k*R,1)*conj(sHankel(nn-1,fluid_k*R,1)));     
        end
        if (testing_type == 4)      % deriviation from G.T.S@2011@JASA
            if irregular_body == 0
                middle_var = 1 + db_s_coeff(nn+1);
            elseif irregular_body == 1
                middle_var = 1 + db_s_coeff(nn+1, nn+mm+1);
            end
        end
        
        n_minus_1 = nn - 1;
        if abs(mm) <= n_minus_1
            if irregular_body == 0
                part_FF = conj(db_bs_coeff(n_minus_1+1,mm+n_minus_1+1)) * conj(db_s_coeff(n_minus_1+1)) * FF(nn-1,mm);
                sum_2 = sum_2 + db_bs_coeff(nn+1,mm+nn+1) * middle_var * part_FF;
            elseif irregular_body == 1
                part_FF = conj(db_bs_coeff(n_minus_1+1,mm+n_minus_1+1)) * conj(db_s_coeff(n_minus_1+1,mm+n_minus_1+1)) * FF(nn-1,mm);
                sum_2 = sum_2 + db_bs_coeff(nn+1,mm+nn+1) * middle_var * part_FF;
            end
        end
        
    end
end
for nn = 0:(N)      % n, m
    for mm = (-nn):(nn)
        if (testing_type == 1)      % deriviation of myself for farfield approximately(a)
            middle_var = 1; 
        end
        if (testing_type == 2)      % deriviation of myself for farfield approximately(b)
            middle_var = (1i^(nn+1))*sin(fluid_k*R-nn*pi/2)*exp(-1i*fluid_k*R);                         
        end
        if (testing_type == 3)      % deriviation of myself for exactly solution
            middle_var = (fluid_k*R)^2 * (sBessel(nn,fluid_k*R,1)*conj(sHankel(nn,fluid_k*R,1)));       
        end
        if (testing_type == 4)      % deriviation from G.T.S@2011@JASA
            middle_var = 0;                                                                           
        end
        sum_3 = sum_3 + db_bs_coeff(nn+1,mm+nn+1) * conj(db_bs_coeff(nn+1,mm+nn+1)) * conj(db_s_coeff(nn+1)) * middle_var;
    end
end

if (testing_type == 1)      % deriviation of myself for farfield approximately(a)
    Frad_z = coeff_force * real((sum_1 - sum_2 +sum_3) / (fluid_k^2));
end
if (testing_type == 2)      % deriviation of myself for farfield approximately(b)
    Frad_z = coeff_force * real((sum_1 - sum_2 +sum_3) / (fluid_k^2));        
end
if (testing_type == 3)      % deriviation of myself for exactly solution
    Frad_z = coeff_force * real((sum_1 + sum_2 +sum_3) / (fluid_k^2));        
end
if (testing_type == 4)      % deriviation from G.T.S@2011@JASA
    Frad_z = coeff_force * real(1i * (sum_1 - sum_2) / (fluid_k^2));        
end


%% anti-rotation of the particle coordinate system for irregular particles

if irregular_body == 1
    Rx=inv([1 0 0;
            0 cos(theta_rotation(1)) -sin(theta_rotation(1)); 
            0 sin(theta_rotation(1)) cos(theta_rotation(1))]);
    Ry=inv([cos(theta_rotation(2)) 0 sin(theta_rotation(2));
            0 1 0; 
            -sin(theta_rotation(2)) 0 cos(theta_rotation(2))]);
    Rz=inv([cos(theta_rotation(3)) -sin(theta_rotation(3)) 0; 
            sin(theta_rotation(3)) cos(theta_rotation(3)) 0;
            0 0 1]);
    Rzyx = Rz * Ry * Rx;
%     Frad_anti = Rzyx_anti * [Frad_x; Frad_y; Frad_z];
    Frad_anti = [Frad_x, Frad_y, Frad_z] * Rzyx;
    Frad_x = Frad_anti(1);
    Frad_y = Frad_anti(2); 
    Frad_z = Frad_anti(3);
end


%% saving file  ('A.mat' for analyses results)

E0 = p_inlet^2 / (2 * fluid_rho * fluid_c^2);
pi_a2_E0 = (pi*particle_radius^2) * E0;
Y_x = Frad_x / pi_a2_E0;        % dimensionaless radiation force function x
Y_y = Frad_y / pi_a2_E0;        % dimensionaless radiation force function y
Y_z = Frad_z / pi_a2_E0;        % dimensionaless radiation force function z

% if (strcmp(wave_type, 'single_transducer') || strcmp(wave_type, 'phase_array_transducer') || strcmp(wave_type, 'phase_array_transducer2')) ~= 1
%     if multi_particle == 0
%         save([dir_file_forces, 'A.mat'], 'E0', 'pi_a2_E0', 'Frad_x', 'Frad_y', 'Frad_z');
%     else
%         particles_Cartesian_data;
%         dir_file_forces_multi = [dir_file_forces_multi, '_Coords('];
%         particle = particle / particle_radius;
%         for ii =1:particle_number
%             particle_code = [num2str(round(particle(ii,1))), num2str(round(particle(ii,2))), num2str(round(particle(ii,3)))];
%             dir_file_forces_multi = [dir_file_forces_multi, particle_code, ','];
%         end
%         dir_file_forces_multi(end) = ')';
%         save([dir_file_forces_multi, 'A.mat'], 'E0', 'pi_a2_E0', ...
%             'Frad_x', 'Frad_y', 'Frad_z', 'particle', 'multi_particle_radius');
%     end
% else        % saving for 'single_transducer' and 'phase_array_transducer' cases
%     if multi_particle == 0
%         save([dir_file_forces_trans, 'A.mat'], 'E0', 'pi_a2_E0', 'transducer', 'Frad_x', 'Frad_y', 'Frad_z');
%     else
%         particles_Cartesian_data;
%         dir_file_forces_trans = [dir_file_forces_trans, '_Coords('];
%         particle = particle / particle_radius;
%         for ii =1:particle_number
%             particle_code = [num2str(round(particle(ii,1))), num2str(round(particle(ii,2))), num2str(round(particle(ii,3)))];
%             dir_file_forces_trans = [dir_file_forces_trans, particle_code, ','];
%         end
%         dir_file_forces_trans(end) = ')';
%         save([dir_file_forces_trans, 'A.mat'], 'E0', 'pi_a2_E0', 'transducer', ...
%             'Frad_x', 'Frad_y', 'Frad_z', 'particle', 'multi_particle_radius');
%     end
% %     ii = 1;
% %     while 1
% %         dir = ['.\data_forces_transducer_', fluid, '\', num2str(ii)];
% %         if exist([dir, 'A.mat']) ~= 0        % if already exist the database, 
% %              ii = ii + 1;                    % then do not create again for saving time.
% %         else
% %             break;
% %         end
% %     end
% %     save([dir, 'A.mat'], 'E0', 'pi_a2_E0', 'Frad_x', 'Frad_y', 'Frad_z');
% end


%% clear database for beam and scattering coefficient

if multi_particle == 0
%     delete([db_filename, '.mat']);
else
    delete([db_filename, '.mat']);
    delete([db_filename_ts, '.mat']);
    delete([db_filename_ibs, '.mat']);
end
%delete([force_matrix_filename, '.mat']);  %% clear database for force matrix
%delete([db_filename_ibs, '.mat']);


%% record the elapse time

% elapse_time2 = datestr(now,'yyyy-mm-dd HH:MM:SS');
% disp(['Ending Time: ',elapse_time2]);
% disp([dir_file_forces, 'A.mat']);
% disp('');
