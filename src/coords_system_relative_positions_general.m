function [r_12, theta_12, phi_12] = coords_system_relative_positions_general(position_1, position_2)
%%
% this function is used to obtain the relative position of different
% particles' coordinates system.
% 
% Firstly, in this function, all particles position information (in term of
% a Visual Absolute Coordinates(VAC)) will be readin. 
% Secondly, based on two input particle indexes, quickly obtain the
% relative position information of each other. 
% Eg., if (indx_1, indx_2) = (1, 2), then 
%   "r_12" means distance between particle "1" and particle "2"; 
%   "theta_12" means polor angle (VAC) of vector "from 1 to 2";
%   "phi_12" means azimuthal angle (VAC) of vector "from 1 to 2";
% NOTE: All particles position information are described under Cartesian
% coordinates of the VAC.
%%

% acos(R) into [0,pi]; asin(R) into [-pi/2,pi/2]; R is vector;
% Therefore:
% 1.            abs(R) = sqrt(R_x^2+R_y^2+R_z^2);
% 2.            theta_R = acos(R_z/abs(R));
% 3.if theta_R ~= 0 or pi, 
%   then:
%       acos(R_x/abs(R)) into [0,pi/2] and asin(R_y/abs(R)) > 0;  
%               vector phi_R = acos(R)
%       acos(R_x/abs(R)) into [pi/2,pi] and asin(R_y/abs(R)) > 0; 
%               vector phi_R = acos(R)
%       acos(R_x/abs(R)) into [0,pi/2] and asin(R_y/abs(R)) < 0;  
%               vector phi_R = 2*pi - acos(R)
%       acos(R_x/abs(R)) into [pi/2,pi] and asin(R_y/abs(R)) > 0; 
%               vector phi_R = 2*pi - acos(R)

%% vector R_12, from 1 point to 2

relative_coords = position_2 - position_1;


%% distance "r_12", polor angle "theta_12"

r_12 = sqrt(sum(relative_coords.^2));
if r_12 == 0
    theta_12 = 0;
else
    theta_12 = acos(relative_coords(3) / r_12);
end

%% azimuthal angle "phi_12"

polor_dis = sqrt(relative_coords(1)^2 + relative_coords(2)^2);
if theta_12 == 0 || theta_12 == pi || polor_dis == 0
    phi_12 = 0;
else
    phi_1 = acos(relative_coords(1) / polor_dis);
    phi_2 = asin(relative_coords(2) / polor_dis);
    if (phi_1 >= 0 && phi_1 < pi/2) && phi_2 >= 0
        phi_12 = phi_1;
    end
    if (phi_1 >= pi/2 && phi_1 <= pi) && phi_2 >= 0
        phi_12 = phi_1;
    end
    if (phi_1 >= 0 && phi_1 < pi/2) && phi_2 < 0
        phi_12 = 2*pi - phi_1;
    end
    if (phi_1 >= pi/2 && phi_1 <= pi) && phi_2 < 0
        phi_12 = 2*pi - phi_1;
    end
end


%% End