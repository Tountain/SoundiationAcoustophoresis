function theta_body_map = theta_w(T, EFn, theta_map, L)
%%
% This function is used to calculate the mapping polar angular coordinate
% 'theta_body_map' using new polar angular coordinate based on Fourier
% series.  
%
%
% Inputs:
%     T: the period of the the body of polar angle (0 to 2*pi), actually,
% only [0,pi] is meaningful.
%     The Fourier expansion coefficients 'EFn', which includes 'En' for the
% first half data and 'Fn' for the rest half data, and new polar angular
% coordinate value 'theta_map'. 
%     L: the truncation number of Fourier expansion series for the old and
% new polar angular coordinate, the deviation of 'theta_body' from
% 'theta_map' on the surface contour of the length-wise slice. 
%
% Output:
%     The mapping polar angular coordinate 'theta_body_map'.
%%

En = EFn(1 : L);
Fn = EFn(L+1 : end);

theta_body_map = theta_map;
for ll = 1 : L
    theta_body_map = theta_body_map + En(ll) * cos(2*pi/T*ll*theta_map) + Fn(ll) * sin(2*pi/T*ll*theta_map);
end

%% limit the range of 'theta_body_map' inside [0, 2*pi] ?

%%