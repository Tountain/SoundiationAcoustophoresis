function [compensation_x, compensation_y, compensation_z] = ...
            renew_coords_of_driving_crystals_with_time(velocity)
%%
% this function is used to renew the coordinate info of the Mie scattering
% crystals in different time and then transfer these new coordinates info
% to the 'particles_Cartesian_data.m' and then correct corresponding positions of
% the crystals structure. 
% We need to give the volecity of the crystal firstly.
%
% Temporarily, testing version only suits for 'uniform linear motion'.
% In future, we can change to arbitary movement.
%%

global time

%% uniform linear motion
compensation_x = velocity(1) * time;
compensation_y = velocity(2) * time;
compensation_z = velocity(3) * time;

%% arbitary movement (getting compensation by integral)

