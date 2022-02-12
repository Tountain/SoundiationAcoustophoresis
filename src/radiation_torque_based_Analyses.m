function [Torque_x,Torque_y,Torque_z,transducer] = radiation_torque_based_Analyses()
%%
% this function is used to obtain the Radiation Torque from the pressure
% distribution surrounding the particle.
%
% NOTE: the theory are based on lots of previous deriviation and some
% personal understandings. 
% In order to finally obtain the Radiation Force, major reference materials
% include (listed but not limited):    
%    P.J.W@1951(JASA); P.J.W@1957(JASA); G.Maidanik@1958(JASA);
%    G.B.A@6th_version(book); E.G.W@1999(book); G.T.S@2011(JASA);
%    G.T.S@2012(EPL); G.T.S@2014(JASA).
%    (ABOVE refs. based on Inviscous Fluid!)
%
%%

testing_type = 4;

% elapse_time1 = datestr(now,'yyyy-mm-dd HH:MM:SS');
% disp(['Started Time: ',elapse_time1]);

% parameters_names;
parameters;

% build the database, if already exist, then it will jump to next sentense.
[db_filename] = database_beam_scattering_coeffs();      
load([db_filename, '.mat'], 'N', 'db_bs_coeff', 'db_s_coeff', 'transducer');


%% obtain the force constant "torque coeff", "BB_p" and "BB_n"

coeff_torque = - (1)^2 / (2 * fluid_rho * fluid_c^2);

BB_p = @(nn,mm) sqrt((nn-mm) * (nn+mm+1));
BB_n = @(nn,mm) sqrt((nn+mm) * (nn-mm+1));


%% obtain the radiation torque rotating with X-axis

sum_1 = 0;
sum_2 = 0;
sum_3 = 0;
sum_4 = 0;
for nn = 0:(N)    % (n), (m-1)
    for mm = -nn:(nn)
        
        if (testing_type == 1)      % deriviation of myself for exactly solution (a) 
            middle_var = (1i) * conj(db_bs_coeff(nn+1,mm+nn+1)) * conj(db_s_coeff(nn+1)) * ...
                        conj(sHankel_partial(nn,fluid_k*R,1)) * ...
                        (sBessel(nn,fluid_k*R,1) + db_s_coeff(nn+1)*sHankel(nn,fluid_k*R,1));
              
            % adjust the projecting relation of database and beam-shape and scattering coefficents
            m_minus_1 = mm - 1;
            if abs(m_minus_1) <= nn    % make sure the adjusting position within the database
                part_1 = (db_bs_coeff(nn+1,m_minus_1+nn+1)) * BB_p(nn,mm-1);
                sum_1 = sum_1 + R^2 * middle_var * part_1;
            end
        end
        
        if (testing_type == 2)      % deriviation of myself for exactly solution (b)
            middle_var = conj(db_s_coeff(nn+1)) * conj(sHankel_partial(nn,fluid_k*R,1)) * sBessel(nn,fluid_k*R,1) + ...
                        db_s_coeff(nn+1) * sHankel(nn,fluid_k*R,1) * sBessel_partial(nn,fluid_k*R,1) + ...
                        conj(db_s_coeff(nn+1)) * db_s_coeff(nn+1) * conj(sHankel_partial(nn,fluid_k*R,1)) * sHankel(nn,fluid_k*R,1);
            
            % adjust the projecting relation of database and beam-shape and scattering coefficents
            m_minus_1 = mm - 1;
            if abs(m_minus_1) <= nn    % make sure the adjusting position within the database
                part_1 = conj(db_bs_coeff(nn+1,m_minus_1+nn+1)) * BB_n(nn,mm);
                sum_1 = sum_1 + R^2 * (1i*db_bs_coeff(nn+1,mm+nn+1)) * middle_var * part_1;
            end
        end
        
        if (testing_type == 3)      % deriviation from G.T.S@2012@IEEE
            if irregular_body == 0
                middle_var = db_bs_coeff(nn+1,mm+nn+1) * db_s_coeff(nn+1);
            elseif irregular_body == 1
                middle_var = db_bs_coeff(nn+1,mm+nn+1) * db_s_coeff(nn+1, nn+mm+1);
            end
           
            % adjust the projecting relation of database and beam-shape and scattering coefficents
            m_minus_1 = mm - 1;
            if abs(m_minus_1) <= nn    % make sure the adjusting position within the database
                if irregular_body == 0
                    part_1 = (-1i * conj(db_bs_coeff(nn+1,m_minus_1+nn+1)) * ...
                            conj(db_s_coeff(nn+1)) - 1i*conj(db_bs_coeff(nn+1,m_minus_1+nn+1))) * ...
                            BB_n(nn,mm);
                    sum_1 = sum_1 + middle_var * part_1;
                elseif irregular_body == 1
                    part_1 = (-1i * conj(db_bs_coeff(nn+1,m_minus_1+nn+1)) * ...
                            conj(db_s_coeff(nn+1,m_minus_1+nn+1)) - 1i*conj(db_bs_coeff(nn+1,m_minus_1+nn+1))) * ...
                            BB_n(nn,mm);
                    sum_1 = sum_1 + middle_var * part_1;
                end
            end
        end
        
        if (testing_type == 4)      % deriviation from G.T.S@2012@IEEE
            if irregular_body == 0
                middle_var = (db_bs_coeff(nn+1,mm+nn+1)) + (db_bs_coeff(nn+1,mm+nn+1) * db_s_coeff(nn+1));
            elseif irregular_body == 1
                middle_var = (db_bs_coeff(nn+1,mm+nn+1)) + (db_bs_coeff(nn+1,mm+nn+1) * db_s_coeff(nn+1, nn+mm+1));
            end
           
            % adjust the projecting relation of database and beam-shape and scattering coefficents
            m_minus_1 = mm - 1;
            if abs(m_minus_1) <= nn    % make sure the adjusting position within the database
                if irregular_body == 0
                    part_1 = conj(db_bs_coeff(nn+1,m_minus_1+nn+1) * ...
                            db_s_coeff(nn+1)) * BB_n(nn,mm);
                    sum_1 = sum_1 + middle_var * part_1;
                elseif irregular_body == 1
                    part_1 = conj(db_bs_coeff(nn+1,m_minus_1+nn+1) * ...
                            db_s_coeff(nn+1,m_minus_1+nn+1)) * BB_n(nn,mm);
                    sum_1 = sum_1 + middle_var * part_1;
                end
            end
        end
        
        if (testing_type == 5)      % deriviation from myself
            if irregular_body == 0
                middle_var = conj(db_bs_coeff(nn+1,mm+nn+1)) + conj(db_bs_coeff(nn+1,mm+nn+1) * db_s_coeff(nn+1));
            elseif irregular_body == 1
                middle_var = conj(db_bs_coeff(nn+1,mm+nn+1)) + conj(db_bs_coeff(nn+1,mm+nn+1) * db_s_coeff(nn+1, nn+mm+1));
            end
           
            % adjust the projecting relation of database and beam-shape and scattering coefficents
            m_minus_1 = mm - 1;
            if abs(m_minus_1) <= nn    % make sure the adjusting position within the database
                if irregular_body == 0
                    part_1 = (db_bs_coeff(nn+1,m_minus_1+nn+1) * ...
                            db_s_coeff(nn+1)) * BB_p(nn,mm);
                    sum_1 = sum_1 + middle_var * part_1;
                elseif irregular_body == 1
                    part_1 = (db_bs_coeff(nn+1,m_minus_1+nn+1) * ...
                            db_s_coeff(nn+1,m_minus_1+nn+1)) * BB_p(nn,mm);
                    sum_1 = sum_1 + middle_var * part_1;
                end
            end
        end
        
    end
end

for nn = 0:(N)      % (n), (m+1)
    for mm = -nn:(nn)
        
        if (testing_type == 1)      % deriviation of myself for exactly solution (a)
            middle_var = (1i) * conj(db_bs_coeff(nn+1,mm+nn+1)) * conj(db_s_coeff(nn+1)) * ...
                        conj(sHankel_partial(nn,fluid_k*R,1)) * ...
                        (sBessel(nn,fluid_k*R,1) + db_s_coeff(nn+1)*sHankel(nn,fluid_k*R,1));

            m_plus_1 = mm + 1;
            if abs(m_plus_1) <= nn
                part_2 = (db_bs_coeff(nn+1,m_plus_1+nn+1)) * BB_n(nn,mm+1);
                sum_2 = sum_2 + R^2 * middle_var * part_2;
            end
        end
        
        if (testing_type == 2)      % deriviation of myself for exactly solution (b)
            middle_var = conj(db_s_coeff(nn+1)) * conj(sHankel_partial(nn,fluid_k*R,1)) * sBessel(nn,fluid_k*R,1) + ...
                        db_s_coeff(nn+1) * sHankel(nn,fluid_k*R,1) * sBessel_partial(nn,fluid_k*R,1) + ...
                        conj(db_s_coeff(nn+1)) * db_s_coeff(nn+1) * conj(sHankel_partial(nn,fluid_k*R,1)) * sHankel(nn,fluid_k*R,1);
            
            % adjust the projecting relation of database and beam-shape and scattering coefficents
            m_plus_1 = mm + 1;
            if abs(m_plus_1) <= nn    % make sure the adjusting position within the database
                part_2 = conj(db_bs_coeff(nn+1,m_plus_1+nn+1)) * BB_p(nn,mm);
                sum_2 = sum_2 + R^2 * (1i*db_bs_coeff(nn+1,mm+nn+1)) * middle_var * part_2;
            end
        end
        
        if (testing_type == 3)      % derivation from G.T.S@2012@IEEE
            if irregular_body == 0
                middle_var = db_bs_coeff(nn+1,mm+nn+1) * db_s_coeff(nn+1);
            elseif irregular_body == 1
                middle_var = db_bs_coeff(nn+1,mm+nn+1) * db_s_coeff(nn+1, nn+mm+1);
            end
           
            % adjust the projecting relation of database and beam-shape and scattering coefficents
            m_plus_1 = mm + 1;
            if abs(m_plus_1) <= nn    % make sure the adjusting position within the database
                if irregular_body == 0
                    part_2 = (-1i * conj(db_bs_coeff(nn+1,m_plus_1+nn+1)) * ...
                            conj(db_s_coeff(nn+1)) - 1i*conj(db_bs_coeff(nn+1,m_plus_1+nn+1))) * ...
                            BB_p(nn,mm);
                    sum_2 = sum_2 + middle_var * part_2;
                elseif irregular_body == 1
                    part_2 = (-1i * conj(db_bs_coeff(nn+1,m_plus_1+nn+1)) * ...
                            conj(db_s_coeff(nn+1,m_plus_1+nn+1)) - 1i*conj(db_bs_coeff(nn+1,m_plus_1+nn+1))) * ...
                            BB_p(nn,mm);
                    sum_2 = sum_2 + middle_var * part_2;
                end
            end
        end
        
        if (testing_type == 4)      % deriviation from G.T.S@2012@IEEE
            if irregular_body == 0
                middle_var = (db_bs_coeff(nn+1,mm+nn+1)) + (db_bs_coeff(nn+1,mm+nn+1) * db_s_coeff(nn+1));
            elseif irregular_body == 1
                middle_var = (db_bs_coeff(nn+1,mm+nn+1)) + (db_bs_coeff(nn+1,mm+nn+1) * db_s_coeff(nn+1, nn+mm+1));
            end
           
            % adjust the projecting relation of database and beam-shape and scattering coefficents
            m_plus_1 = mm + 1;
            if abs(m_plus_1) <= nn    % make sure the adjusting position within the database
                if irregular_body == 0
                    part_2 = conj(db_bs_coeff(nn+1,m_plus_1+nn+1) * ...
                            db_s_coeff(nn+1)) * BB_p(nn,mm);
                    sum_2 = sum_2 + middle_var * part_2;
                elseif irregular_body == 1
                    part_2 = conj(db_bs_coeff(nn+1,m_plus_1+nn+1) * ...
                            db_s_coeff(nn+1,m_plus_1+nn+1)) * BB_p(nn,mm);
                    sum_2 = sum_2 + middle_var * part_2;
                end
            end
        end
        
        if (testing_type == 5)      % derivation from myself
            if irregular_body == 0
                middle_var = conj(db_bs_coeff(nn+1,mm+nn+1)) + conj(db_bs_coeff(nn+1,mm+nn+1) * db_s_coeff(nn+1));
            elseif irregular_body == 1
                middle_var = conj(db_bs_coeff(nn+1,mm+nn+1)) + conj(db_bs_coeff(nn+1,mm+nn+1) * db_s_coeff(nn+1, nn+mm+1));
            end
           
            % adjust the projecting relation of database and beam-shape and scattering coefficents
            m_plus_1 = mm + 1;
            if abs(m_plus_1) <= nn    % make sure the adjusting position within the database
                if irregular_body == 0
                    part_2 = (db_bs_coeff(nn+1,m_plus_1+nn+1) * ...
                            db_s_coeff(nn+1)) * BB_n(nn,mm);
                    sum_2 = sum_2 + middle_var * part_2;
                elseif irregular_body == 1
                    part_2 = (db_bs_coeff(nn+1,m_plus_1+nn+1) * ...
                            db_s_coeff(nn+1,m_plus_1+nn+1)) * BB_n(nn,mm);
                    sum_2 = sum_2 + middle_var * part_2;
                end
            end
        end
        
    end
end

for nn = 0:(N)    	% (n), (m+1), conj
    for mm = (-nn):(nn)
        
        if (testing_type == 1)      % deriviation of myself for exactly solution (a)
            middle_var = (1i) * (db_bs_coeff(nn+1,mm+nn+1)) * (db_s_coeff(nn+1)) * ...
                        sHankel(nn,fluid_k*R,1) * sBessel_partial(nn,fluid_k*R,1);
        
            m_plus_1 = mm + 1;
            if abs(m_plus_1) <= nn
                part_3 = conj(db_bs_coeff(nn+1,m_plus_1+nn+1)) * BB_p(nn,mm);
                sum_3 = sum_3 + R^2 * middle_var * part_3;
            end
        end
        
    end
end

for nn = 0:(N)      % (n), (m-1), conj
    for mm = (-nn):(nn)
        
        if (testing_type == 1)      % deriviation of myself for exactly solution (a)
            middle_var = (1i) * (db_bs_coeff(nn+1,mm+nn+1)) * (db_s_coeff(nn+1)) * ...
                        sHankel(nn,fluid_k*R,1) * sBessel_partial(nn,fluid_k*R,1);
        
            m_minus_1 = mm - 1;
            if abs(m_minus_1) <= nn
                part_4 = conj(db_bs_coeff(nn+1,nn+m_minus_1+1)) * BB_n(nn,mm);
                sum_4 = sum_4 + R^2 * middle_var * part_4;
            end
        end
        
    end
end

if (testing_type == 1)      % deriviation of myself for exactly solution (a)
    Torque_x = coeff_torque / (fluid_k) * 1/2 * real(sum_1 + sum_2 + sum_3 + sum_4);   
end
if (testing_type == 2)      % deriviation of myself for exactly solution (b)
    Torque_x = coeff_torque / (fluid_k) * 1/2 * real(sum_2 + sum_1);   
end
if (testing_type == 3)      % deriviation from myself
    Torque_x = coeff_torque / (fluid_k)^3 * imag(sum_2 + sum_1);     
end
if (testing_type == 4)      % deriviation from G.T.S@2012@IEEE
%     Torque_x = - 1 / (pi*(particle_radius*fluid_k)^3) * 1/2 * real(sum_2 + sum_1);   
    Torque_x = coeff_torque / (fluid_k^3) * 1/2 * real(sum_2 + sum_1);
end
if (testing_type == 5)      % deriviation from myself
    Torque_x = - coeff_torque / (fluid_k^3) * real(sum_1 + sum_2);
end

%% obtain the radiation torque rotating with Y-axis

if (testing_type == 1)      % deriviation of myself for exactly solution (a)
    Torque_y = coeff_torque / (fluid_k) * 1/2 * real((sum_1 - sum_2 + sum_3 - sum_4)/(1i));
end
if (testing_type == 2)      % deriviation of myself for exactly solution (b)
	Torque_y = coeff_torque / (fluid_k) * 1/2 * real((sum_2 - sum_1)/(1i));
end
if (testing_type == 3)      % deriviation from myself
    Torque_y = - coeff_torque / (fluid_k)^3 * real(sum_2 - sum_1);     
end
if (testing_type == 4)      % deriviation from G.T.S@2012@IEEE
%     Torque_y = - 1 / (pi*(particle_radius*fluid_k)^3) * 1/2 * imag(sum_2 - sum_1);   
    Torque_y = coeff_torque / (fluid_k^3) * 1/2 * imag(sum_2 - sum_1);   
end
if (testing_type == 5)      % deriviation from myself
    Torque_y = - coeff_torque / (fluid_k^3) * imag(sum_1 - sum_2);
end

%% obtain the radiation torque rotating with Z-axis (parallel with incident wave)

sum_1 = 0;
sum_2 = 0;
sum_3 = 0;
for nn = 0:(N)     % (n), m
    for mm = -nn:(nn)
        
        if (testing_type == 1)      % deriviation of myself for exactly solution (a)
            middle_var = (1i*mm) * conj(db_bs_coeff(nn+1,mm+nn+1)) * ...
                        db_bs_coeff(nn+1,mm+nn+1);

            part_1 = conj(db_s_coeff(nn+1)) * conj(sHankel_partial(nn,fluid_k*R,1)) * ...
                    sBessel(nn,fluid_k*R,1);
            sum_1 = sum_1 + R^2 * middle_var * part_1;
        end
        
        if (testing_type == 2)      % deriviation of myself for exactly solution (b)
            middle_var = conj(db_s_coeff(nn+1)) * conj(sHankel_partial(nn,fluid_k*R,1)) * sBessel(nn,fluid_k*R,1) + ...
                        db_s_coeff(nn+1) * sHankel(nn,fluid_k*R,1) * sBessel_partial(nn,fluid_k*R,1) + ...
                        conj(db_s_coeff(nn+1)) * db_s_coeff(nn+1) * conj(sHankel_partial(nn,fluid_k*R,1)) * sHankel(nn,fluid_k*R,1);
            
            part_1 = conj(db_bs_coeff(nn+1,mm+nn+1));
            sum_1 = sum_1 + R^2 * (1i*mm*db_bs_coeff(nn+1,mm+nn+1)) * middle_var * part_1;
        end
        
        if (testing_type == 3)      % deriviation from myself
            if irregular_body == 0
                middle_var = mm * (db_bs_coeff(nn+1,mm+nn+1) * ...
                            db_s_coeff(nn+1));
                part_1 = (conj(db_bs_coeff(nn+1,mm+nn+1)) * ...
                            conj(db_s_coeff(nn+1)) - 1i*conj(db_bs_coeff(nn+1,mm+nn+1)));
                sum_1 = sum_1 + middle_var * part_1;
            elseif irregular_body == 1
                middle_var = mm * (db_bs_coeff(nn+1,mm+nn+1) * ...
                            db_s_coeff(nn+1, nn+mm+1));
                part_1 = (conj(db_bs_coeff(nn+1,mm+nn+1)) * ...
                            conj(db_s_coeff(nn+1,mm+nn+1)) - 1i*conj(db_bs_coeff(nn+1,mm+nn+1)));
                sum_1 = sum_1 + middle_var * part_1;
            end
        end
        
        if (testing_type == 4)      % deriviation from G.T.S@2012@IEEE
            if irregular_body == 0
                middle_var = mm * conj(db_bs_coeff(nn+1,mm+nn+1) * ...
                            db_s_coeff(nn+1));
                part_1 = db_bs_coeff(nn+1,mm+nn+1);
                sum_1 = sum_1 + middle_var * part_1;
            elseif irregular_body == 1
                middle_var = mm * conj(db_bs_coeff(nn+1,mm+nn+1) * ...
                            db_s_coeff(nn+1, nn+mm+1));
                part_1 = db_bs_coeff(nn+1,mm+nn+1);
                sum_1 = sum_1 + middle_var * part_1;
            end
        end
        
        if (testing_type == 5)      % deriviation from myself
            if irregular_body == 0
                middle_var = mm * (db_bs_coeff(nn+1,mm+nn+1) * ...
                            db_s_coeff(nn+1));
                part_1 = conj(db_bs_coeff(nn+1,mm+nn+1));
                sum_1 = sum_1 + middle_var * part_1;
            elseif irregular_body == 1
                middle_var = mm * (db_bs_coeff(nn+1,mm+nn+1) * ...
                            db_s_coeff(nn+1, nn+mm+1));
                part_1 = conj(db_bs_coeff(nn+1,mm+nn+1));
                sum_1 = sum_1 + middle_var * part_1;
            end
        end
        
    end
end

for nn = 0:(N)      % (n), m
    for mm = -(nn):(nn)
        
        if (testing_type == 1)      % deriviation of myself for exactly solution (a)
            middle_var = (1i*mm) * conj(db_bs_coeff(nn+1,mm+nn+1)) * ...
                        db_bs_coeff(nn+1,mm+nn+1);
        
            part_2 = (db_s_coeff(nn+1)) * (sHankel(nn,fluid_k*R,1)) * ...
                    sBessel_partial(nn,fluid_k*R,1);
            sum_2 = sum_2 + R^2 * middle_var * part_2;
        end
        
        if (testing_type == 4)      % deriviation from G.T.S@2012@IEEE
            if irregular_body == 0
                middle_var = mm * conj(db_bs_coeff(nn+1,mm+nn+1) * ...
                            db_s_coeff(nn+1));
                part_2 = db_bs_coeff(nn+1,mm+nn+1) * db_s_coeff(nn+1);
                sum_2 = sum_2 + middle_var * part_2;
            elseif irregular_body == 1
                middle_var = mm * conj(db_bs_coeff(nn+1,mm+nn+1) * ...
                            db_s_coeff(nn+1,mm+nn+1));
                part_2 = db_bs_coeff(nn+1,mm+nn+1) * db_s_coeff(nn+1,mm+nn+1);
                sum_2 = sum_2 + middle_var * part_2;
%                 middle_var = mm * db_bs_coeff(nn+1,mm+nn+1)^2 * ...
%                             conj(db_s_coeff(nn+1,mm+nn+1));
%                 part_2 = 1 + db_s_coeff(nn+1,mm+nn+1);
%                 sum_2 = sum_2 + middle_var * part_2;
            end
        end
        
        if (testing_type == 5)      % deriviation from myself
            if irregular_body == 0
                middle_var = mm * (db_bs_coeff(nn+1,mm+nn+1) * ...
                            db_s_coeff(nn+1));
                part_2 = conj(db_bs_coeff(nn+1,mm+nn+1) * db_s_coeff(nn+1));
                sum_2 = sum_2 + middle_var * part_2;
            elseif irregular_body == 1
                middle_var = mm * (db_bs_coeff(nn+1,mm+nn+1) * ...
                            db_s_coeff(nn+1,mm+nn+1));
                part_2 = conj(db_bs_coeff(nn+1,mm+nn+1) * db_s_coeff(nn+1,mm+nn+1));
                sum_2 = sum_2 + middle_var * part_2;
            end
        end
        
    end
end

for nn = 0:(N)      % (n), m
    for mm = (-nn):(nn)
        
        if (testing_type == 1)      % deriviation of myself for exactly solution (a)
            middle_var = (1i*mm) * conj(db_bs_coeff(nn+1,mm+nn+1)) * ...
                        db_bs_coeff(nn+1,mm+nn+1);
        
            part_3 = conj(db_s_coeff(nn+1)) * db_s_coeff(nn+1) * ...
                    (sHankel(nn,fluid_k*R,1)) * conj(sHankel_partial(nn,fluid_k*R,1));
            sum_3 = sum_3 + R^2 * middle_var * part_3;
        end
        
    end
end

if (testing_type == 1)      % deriviation of myself for exactly solution (a)
    Torque_z = coeff_torque / (fluid_k) * real(sum_1 + sum_2 +sum_3);
end
if (testing_type == 2)      % deriviation of myself for exactly solution (b)
	Torque_z = coeff_torque / (fluid_k) * real(sum_1);
end
if (testing_type == 3)      % deriviation of myself for exactly solution (b)
	Torque_z = 2 * coeff_torque / (fluid_k)^3 * imag(sum_1);
end
if (testing_type == 4)      % deriviation from G.T.S@2012@IEEE
%     Torque_z = - 1 / (pi*(particle_radius*fluid_k)^3) * real(sum_1 + sum_2);   
    Torque_z = coeff_torque / (fluid_k^3) * real(sum_1 + sum_2);  
end
if (testing_type == 5)      % deriviation from G.T.S@2012@IEEE
    Torque_z = - 2 * coeff_torque / (fluid_k^3) * real(sum_1 + sum_2);  
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
%     Torque_anti = Rzyx_anti * [Torque_x; Torque_y; Torque_z];
    Torque_anti = [Torque_x, Torque_y, Torque_z] * Rzyx;
    Torque_x = Torque_anti(1);
    Torque_y = Torque_anti(2);
    Torque_z = Torque_anti(3);
end


%% saving file  ('A.mat' for analyses results)

E0 = p_inlet^2 / (2 * fluid_rho * fluid_c^2);   % J/m3
pi_a3_E0 = (pi*particle_radius^3) * E0;     % J or N*m
% Torque_x = pi_a3_E0 * Torque_x;
% Torque_y = pi_a3_E0 * Torque_y;
% Torque_z = pi_a3_E0 * Torque_z;
N_x = Torque_x / pi_a3_E0;        % dimensionaless radiation torque function x
N_y = Torque_y / pi_a3_E0;        % dimensionaless radiation torque function y
N_z = Torque_z / pi_a3_E0;        % dimensionaless radiation torque function z
%[N_x,N_y,N_z]

% % save([dir_file_torques, 'A.mat'], 'E0', 'pi_a3_E0', 'Torque_x', 'Torque_y', 'Torque_z');
% if (strcmp(wave_type, 'single_transducer') || strcmp(wave_type, 'phase_array_transducer') || strcmp(wave_type, 'phase_array_transducer2')) ~= 1
%     if multi_particle == 0
%         save([dir_file_torques, 'A.mat'], 'E0', 'pi_a3_E0', 'Torque_x', 'Torque_y', 'Torque_z');
%     else
%         particles_Cartesian_data;
%         dir_file_torques_multi = [dir_file_torques_multi, '_Coords('];
%         particle = particle / particle_radius;
%         for ii =1:particle_number
%             particle_code = [num2str(round(particle(ii,1))), num2str(round(particle(ii,2))), num2str(round(particle(ii,3)))];
%             dir_file_torques_multi = [dir_file_torques_multi, particle_code, ','];
%         end
%         dir_file_torques_multi(end) = ')';
%         save([dir_file_torques_multi, 'A.mat'], 'E0', 'pi_a2_E0', ...
%             'Frad_x', 'Frad_y', 'Frad_z', 'particle', 'multi_particle_radius');
%     end
% else        % saving for 'single_transducer' and 'phase_array_transducer' cases
%     if multi_particle == 0
%         save([dir_file_torques_trans, 'A.mat'], 'E0', 'pi_a3_E0', 'transducer', 'Torque_x', 'Torque_y', 'Torque_z');
%     else
%         particles_Cartesian_data;
%         dir_file_torques_trans = [dir_file_torques_trans, '_Coords('];
%         particle = particle / particle_radius;
%         for ii =1:particle_number
%             particle_code = [num2str(round(particle(ii,1))), num2str(round(particle(ii,2))), num2str(round(particle(ii,3)))];
%             dir_file_torques_trans = [dir_file_torques_trans, particle_code, ','];
%         end
%         dir_file_torques_trans(end) = ')';
%         save([dir_file_torques_trans, 'A.mat'], 'E0', 'pi_a2_E0', 'transducer', ...
%             'Frad_x', 'Frad_y', 'Frad_z', 'particle', 'multi_particle_radius');
%     end
% %     ii = 1;
% %     while 1
% %         dir = ['.\data_torques_transducer_', fluid, '\', num2str(ii)];
% %         if exist([dir, 'A.mat']) ~= 0        % if already exist the database, 
% %              ii = ii + 1;                    % then do not create again for saving time.
% %         else
% %             break;
% %         end
% %     end
% %     save([dir, 'A.mat'], 'E0', 'pi_a3_E0', 'Torque_x', 'Torque_y', 'Torque_z');
% end

%% clear database for beam and scattering coefficient

delete([db_filename, '.mat']);
%delete([force_matrix_filename, '.mat']);  %% clear database for force matrix
%delete([db_filename_ibs, '.mat']);


%% record the elapse time

% elapse_time2 = datestr(now,'yyyy-mm-dd HH:MM:SS');
% disp(['Ending Time: ',elapse_time2]);
% disp([dir_file_torques, 'A.mat']);
% disp('');

