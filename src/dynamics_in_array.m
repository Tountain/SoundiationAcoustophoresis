function [time, X_position, Y_position, Z_position, X_rotation, Y_rotation, Z_rotation] = ...
            dynamics_in_array(X_vec, Y_vec, Z_vec, theta_x, theta_y, theta_z, delta_t, ending_t, handles)
%%
% This function is used to animate the translation and rotation state of
% various irregular bodies. 
%
% 1. we first give a incident angle, and call the function
% "radiation_force_based_Analyses.m" and
% "radiation_torque_based_Analyses.m" to obtain the radiation force and
% torque on the body. 
% 2. based on the radiation force and torque, as well as the drag force and
% torque to evaluate the angular velocity of the body. The drag force and
% torque are defined by a linear relation of drag coefficient and angular
% velocity. (For the first iterative step, the drag torque is set to zero).  
% As the translation and rotation speed are very slow, we assume that the
% moving and angular acceleration velocity is zero. So, the we determine
% the angular velocity by applying 'Frad + Fdrag = 0' and 'Trad + Tdrag =
% 0'. 
% 3. given a time step 'delta_t', based on the velocity and angular
% velocity obtained in above step to adjust the body's positions and
% postures.  
% 4. after adjusting the positions and postures of the body, repeat the
% above 2 to 4 step until the changes of velocity and angular velocity
% within acceptable level. 
%
%%

elapse_time1 = datestr(now,'yyyy-mm-dd HH:MM:SS');
disp(['Started Time: ',elapse_time1]);

parameters; 
Fg = particle_mass * 9.8;
% Fg = 0;
Stokes_Law = 1;

time_ratio = 10;        % rotational time step is 'time_ratio' times smaller than translational time step
error = +inf;
step = 1;
time = 0;
cur_time = 0;
max_step = ending_t / (delta_t);    % maximum simulated step


%% recording variables

% giving a initial position of probe particle "OL" 
% X_vec = 0.002;
% Y_vec = 0.002;
% Z_vec = 0*particle_radius;
% saving position information
X_position = X_vec;
Y_position = Y_vec;
Z_position = Z_vec;
% revise the X, Y and Z positions of the particle based on current forces and "delta_t"
% probe_velocity_x = 0;
% probe_velocity_y = 0;
% probe_velocity_z = 0;
delta_x = 0;
delta_y = 0;
delta_z = 0;
max_movement = 0;


% giving a initial posture of probe particle
% theta_x = 0;
% theta_y = 0;
% theta_z = 0;
% saving rotational angle
X_rotation = theta_x;
Y_rotation = theta_y;
Z_rotation = theta_z;
% revise the X, Y and Z postures of the particle based on current torques and "delta_t"
% rotational_velocity_x = 0;
% rotational_velocity_y = 0;
% rotational_velocity_z = 0;
delta_theta_x = 0;      % rotational angle influent by the drag torque and radiation torque within a time interval
delta_theta_y = 0;      % rotational angle influent by the drag torque and radiation torque within a time interval
delta_theta_z = 0;      % rotational angle influent by the drag torque and radiation torque within a time interval
max_rotation = 0;
% deviation_from_x_axis = 0;      % we define at the initial state, the body axisymmetric axis coincide with the x-axis.
% deviation_from_y_axis = 0;      % we define at the initial state, the body axisymmetric axis coincide with the y-axis.
% deviation_from_z_axis = 0;      % we define at the initial state, the body axisymmetric axis coincide with the z-axis.
% orientation_angle = [deviation_from_x_axis, deviation_from_y_axis, deviation_from_z_axis];


fprintf('Initial position point (%f, %f, %f).\n', X_vec, Y_vec, Z_vec);
fprintf('Initial incidence angle (%f, %f, %f).\n', theta_x, theta_y, theta_z);


%%
while error > 0.05 && step < max_step
    
    %%
    replace_position_theta([X_vec, Y_vec, Z_vec], [theta_x, theta_y, theta_z]);       % replace the position vector 'X_vec,Y_vec,Z_vec' and rotation angle 'theta_rotation' in "parameters.m"

    [Frad_x, Frad_y, Frad_z, transducer] = radiation_force_based_Analyses();

    Frad_z = Frad_z - Fg;       % considering the Gravity
    
    [Torque_x, Torque_y, Torque_z, ~] = radiation_torque_based_Analyses();

    %=================== the Stokes drag force is too big =============
    %     if (Frad_x^2+Frad_y^2+Frad_z^2) < (F_stokes_x^2+F_stokes_y^2+F_stokes_z^2)
    %         fprintf('The radiation forces can not move the probe particle, coz the viscosity is too high!\n');
    %         fprintf('please reduce the time step: ''delta_t'', and then try again!\n');
    %         break;
    %     end
    %==================================================================
    

    %% translation
    if Stokes_Law == 0
        a_x = (Frad_x) / particle_mass;
        a_y = (Frad_y) / particle_mass;
        a_z = (Frad_z) / particle_mass;
        delta_x = 0.5 * a_x * delta_t^2;
        delta_y = 0.5 * a_y * delta_t^2;
        delta_z = 0.5 * a_z * delta_t^2;
    elseif Stokes_Law == 1
        % assume that the Frad will balance with Fstokes in will short
        % distance, then will ignore this distance and consider the
        % whole physics behind is uniform velocity movement.
        probe_velocity_x = Frad_x / (6 * pi * fluid_viscosity * particle_radius);
        probe_velocity_y = Frad_y / (6 * pi * fluid_viscosity * particle_radius);
        probe_velocity_z = Frad_z / (6 * pi * fluid_viscosity * particle_radius);
        delta_x = probe_velocity_x * delta_t;
        delta_y = probe_velocity_y * delta_t;
        delta_z = probe_velocity_z * delta_t;
    end

    X_vec = X_vec + delta_x;
    Y_vec = Y_vec + delta_y;
    Z_vec = Z_vec + delta_z;
    if max_movement < abs(max([delta_x, delta_y, delta_z]))
        max_movement = abs(max([delta_x, delta_y, delta_z]));
    end
    
    error_translation = sqrt(delta_x^2 + delta_y^2 + delta_z^2) / max_movement;
    fprintf('Moving to point (%f, %f, %f) at Time %f s! Step = %d; Error = %f.\n', ...
                X_vec, Y_vec, Z_vec, cur_time, step, error_translation);        
    
	% recording the positions
	X_position = [X_position, X_vec];
    Y_position = [Y_position, Y_vec];
    Z_position = [Z_position, Z_vec];
    
    if (inter_dist+Z_vec) <= 0.012  % meaning the tweezer can not stably trap the object
        fprintf('The particle cannot be stably trapped in current setting.\n');
        break;
    end
    
    
    %% rotation
    if Stokes_Law == 0
        ar_x = (Torque_x) / motion_inertia;
        ar_y = (Torque_y) / motion_inertia;
        ar_z = (Torque_z) / motion_inertia;
        delta_theta_x = 0.5 * ar_x * delta_t^2;
        delta_theta_y = 0.5 * ar_y * delta_t^2;
        delta_theta_z = 0.5 * ar_z * delta_t^2;
    elseif Stokes_Law == 1
        rotational_velocity_x = Torque_x / (8 * pi * fluid_viscosity * particle_radius^3);
        rotational_velocity_y = Torque_y / (8 * pi * fluid_viscosity * particle_radius^3);
        rotational_velocity_z = Torque_z / (8 * pi * fluid_viscosity * particle_radius^3);
        delta_theta_x = rotational_velocity_x * delta_t;
        delta_theta_y = rotational_velocity_y * delta_t;
        delta_theta_z = rotational_velocity_z * delta_t;
    end
    
    %================ the relative incidence angle ====================
    theta_x = theta_x + delta_theta_x;
    theta_y = theta_y + delta_theta_y;
    theta_z = theta_z + delta_theta_z;
    % deviation_from_x_axis = deviation_from_x_axis + delta_theta_x;
    % deviation_from_y_axis = deviation_from_y_axis + delta_theta_y;
    % deviation_from_z_axis = deviation_from_z_axis + delta_theta_z;
    %==================================================================
    
    X_rotation = [X_rotation, theta_x];
    Y_rotation = [Y_rotation, theta_y];
    Z_rotation = [Z_rotation, theta_z];
%     orientation_angle = [orientation_angle, deviation_from_z_axis];
    
    if max_rotation < abs(max([delta_theta_x, delta_theta_y, delta_theta_z])) && step > 10
        max_rotation = abs(max([delta_theta_x, delta_theta_y, delta_theta_z]));
    end
    
    error_rotation = sqrt(delta_theta_x^2 + delta_theta_y^2 + delta_theta_z^2) / (max_rotation+eps);
    fprintf('Rotate to angle (%f, %f, %f) at Time %f s. Step = %d; Error = %f.\n', ...
                theta_x*180/pi, theta_y*180/pi, theta_z*180/pi, cur_time, step, error_rotation);        
    
    set(handles.current_time, 'String', num2str(cur_time));
    pause(0.001);
    
    %%
    
    cur_time = cur_time + delta_t;
    time = [time, cur_time];
    error = max(error_translation, error_rotation);
    step = step + 1;
    
end


%% visualizing the translation and rotation trajectories

color_sphere = [0, 0, 1];
color_elli = [0, 0.78, 0.55];
color_cone = [1, 0.38, 0.0];
color_diamond = [0.15, 0.15, 0.15];
figure_FontSize = 20;

% 3D
cla(handles.axes10);
axes(handles.axes10);
% figure(10);
box on; hold on;    % visualizing the trace with pressure contour

if length(Cn) == 2
    plot3(X_position, Y_position, Z_position, 'color', color_sphere, 'linewidth', 2);
elseif length(Cn) == 3
    plot3(X_position, Y_position, Z_position, 'color', color_elli, 'linewidth', 2);
elseif length(Cn) == 4
    plot3(X_position, Y_position, Z_position, 'color', color_cone, 'linewidth', 2);
elseif length(Cn) == 5
    plot3(X_position, Y_position, Z_position, 'color', color_diamond, 'linewidth', 2);
else
    plot3(X_position, Y_position, Z_position, 'r-', 'linewidth', 2);
end
% plot3(X_position, Y_position, Z_position, '.b');

% here, I edit "comet, comet3" and add a "pause" function
% comet(delta_t/10, X_position, Z_position, 0.1);
% comet3(X_position, Y_position, Z_position, 0.1);

%========= graph the initial and final position of the crystals ===========
% initial_final_position_of_crystals(fig_num, probe_initial_position)
%==========================================================================

ez = [0,0,0.005]*1;
particle_orientation = zeros(length(X_rotation),3);
for ii = 1:length(X_rotation)
    Rx=[1 0 0;
        0 cos(-X_rotation(ii)) -sin(-X_rotation(ii)); 
        0 sin(-X_rotation(ii)) cos(-X_rotation(ii))];
    Ry=[cos(-Y_rotation(ii)) 0 sin(-Y_rotation(ii));
        0 1 0; 
        -sin(-Y_rotation(ii)) 0 cos(-Y_rotation(ii))];
    Rz=[cos(-Z_rotation(ii)) -sin(-Z_rotation(ii)) 0; 
        sin(-Z_rotation(ii)) cos(-Z_rotation(ii)) 0; 
        0 0 1];
    Rxyz = Rx * Ry * Rz;
    particle_orientation(ii,:) = ez * Rxyz;
end

% quiver3(X_position(1:time_ratio:end)', Y_position(1:time_ratio:end)', Z_position(1:time_ratio:end)', ...
%     particle_orientation(1:time_ratio:end,1), particle_orientation(1:time_ratio:end,2), particle_orientation(1:time_ratio:end,3));
h = sum(ez)*0.03;
r1 = h/2;
r2 = r1 * 3;
color = linspace(0,1,length(X_rotation));         % 0 to 0.1 s 
time_ratio = ceil(length(X_rotation)/25);
for ii = 1:time_ratio:length(X_rotation)
    X1 = [X_position(ii), Y_position(ii), Z_position(ii)];
    V1 = [particle_orientation(ii, 1), particle_orientation(ii, 2), particle_orientation(ii, 3)];
    arrow3(X1, V1, r1, r2, h, [1, color(ii), 0]);
end

transducer = str2num(get(handles.trans_position, 'String'));
visualize_transducers(transducer);


xlabel('\rm{\fontname{Times new roman}\it{x}{\rm{-axis [mm]}}}');
ylabel('\rm{\fontname{Times new roman}\it{y}{\rm{-axis [mm]}}}');
zlabel('\rm{\fontname{Times new roman}\it{z}{\rm{-axis [mm]}}}');
set(gca,'XTick', [-0.01, 0, 0.01]);
set(gca,'XTickLabel', [-10, 0, 10]);
set(gca,'YTick', [-0.01, 0, 0.01]);
set(gca,'YTickLabel', [-10, 0, 10]);
set(gca,'ZTick', [-0.06, -0.04, -0.02, -0.00, 0.02, 0.04]);
set(gca,'ZTickLabel', [-60, -40, -20, 0, 20, 40]);
set(gca, 'FontName', 'Times new roman');
% set(get(gca,'XLabel'),'FontSize',figure_FontSize);
% set(get(gca,'YLabel'),'FontSize',figure_FontSize);
% set(get(gca,'ZLabel'),'FontSize',figure_FontSize);
% set(findobj('Fontsize',10),'fontsize', figure_FontSize-5);
view([1,0.3,0.3]); camlight; lighting gouraud;
axis image;
% rr = range_r_coeff * particle_radius;       % visible window range, 
axis([-0.015, 0.015, -0.015, 0.015, -0.06, 0.04]);

hold off;


% 2D: X position

cla(handles.axes4);
axes(handles.axes4);
% figure(4);
box on; hold on;   % visualizing the trace 

if length(Cn) == 2
    plot(time*1000, X_position(1:length(time))*1000, '-', 'color', color_sphere, 'linewidth', 2);
elseif length(Cn) == 3
    plot(time*1000, X_position(1:length(time))*1000, '-', 'color', color_elli, 'linewidth', 2);
elseif length(Cn) == 4
    plot(time*1000, X_position(1:length(time))*1000, '-', 'color', color_cone, 'linewidth', 2);
elseif length(Cn) == 5
    plot(time*1000, X_position(1:length(time))*1000, '-', 'color', color_diamond, 'linewidth', 2);
else
    plot(time*1000, X_position(1:length(time))*1000, 'r-', 'linewidth', 2);
end
    
xlabel('\rm{\fontname{Times new roman}\rm{Time }\rm{[ms]}}');
ylabel('\rm{\fontname{Times new roman}\rm{Position }{\it{x}}\rm{-axis [mm]}}');
set(gca, 'FontName', 'Times new roman');
% set(get(gca,'XLabel'),'FontSize',figure_FontSize);
% set(get(gca,'YLabel'),'FontSize',figure_FontSize);
% set(get(gca,'ZLabel'),'FontSize',figure_FontSize);
% set(findobj('Fontsize',10),'fontsize', figure_FontSize-5);
% set(gcf,'unit','centimeters','position',[10 3 9 9])

hold off;

% 2D: Y position

cla(handles.axes5);
axes(handles.axes5);
% figure(5);
box on; hold on;   % visualizing the trace 

if length(Cn) == 2
    plot(time*1000, Y_position(1:length(time))*1000, '-', 'color', color_sphere, 'linewidth', 2);
elseif length(Cn) == 3
    plot(time*1000, Y_position(1:length(time))*1000, '-', 'color', color_elli, 'linewidth', 2);
elseif length(Cn) == 4
    plot(time*1000, Y_position(1:length(time))*1000, '-', 'color', color_cone, 'linewidth', 2);
elseif length(Cn) == 5
    plot(time*1000, Y_position(1:length(time))*1000, '-', 'color', color_diamond, 'linewidth', 2);
else
    plot(time*1000, Y_position(1:length(time))*1000, 'r-', 'linewidth', 2);
end
    
xlabel('\rm{\fontname{Times new roman}\rm{Time }\rm{[ms]}}');
ylabel('\rm{\fontname{Times new roman}\rm{Position }{\it{y}}\rm{-axis [mm]}}');
set(gca, 'FontName', 'Times new roman');
% set(get(gca,'XLabel'),'FontSize',figure_FontSize);
% set(get(gca,'YLabel'),'FontSize',figure_FontSize);
% set(get(gca,'ZLabel'),'FontSize',figure_FontSize);
% set(findobj('Fontsize',10),'fontsize', figure_FontSize-5);
% set(gcf,'unit','centimeters','position',[10 3 9 9])

hold off;

% 2D: Z position

cla(handles.axes6);
axes(handles.axes6);
% figure(6);
box on; hold on;   % visualizing the trace 

if length(Cn) == 2
    plot(time*1000, Z_position(1:length(time))*1000, '-', 'color', color_sphere, 'linewidth', 2);
elseif length(Cn) == 3
    plot(time*1000, Z_position(1:length(time))*1000, '-', 'color', color_elli, 'linewidth', 2);
elseif length(Cn) == 4
    plot(time*1000, Z_position(1:length(time))*1000, '-', 'color', color_cone, 'linewidth', 2);
elseif length(Cn) == 5
    plot(time*1000, Z_position(1:length(time))*1000, '-', 'color', color_diamond, 'linewidth', 2);
else
    plot(time*1000, Z_position(1:length(time))*1000, 'r-', 'linewidth', 2);
end
    
xlabel('\rm{\fontname{Times new roman}\rm{Time }\rm{[ms]}}');
ylabel('\rm{\fontname{Times new roman}\rm{Position }{\it{z}}\rm{-axis [mm]}}');
set(gca, 'FontName', 'Times new roman');
% set(get(gca,'XLabel'),'FontSize',figure_FontSize);
% set(get(gca,'YLabel'),'FontSize',figure_FontSize);
% set(get(gca,'ZLabel'),'FontSize',figure_FontSize);
% set(findobj('Fontsize',10),'fontsize', figure_FontSize-5);
% set(gcf,'unit','centimeters','position',[10 3 9 9])

hold off;

% 2D: X orientation

cla(handles.axes7);
axes(handles.axes7);
% figure(7);
box on; hold on;   % visualizing the trace

if length(Cn) == 2
    plot(time*1000, X_rotation(1:length(time))*180/pi, '-', 'color', color_sphere, 'linewidth', 2);
elseif length(Cn) == 3
    plot(time*1000, X_rotation(1:length(time))*180/pi, '-', 'color', color_elli, 'linewidth', 2);
elseif length(Cn) == 4
    plot(time*1000, X_rotation(1:length(time))*180/pi, '-', 'color', color_cone, 'linewidth', 2);
elseif length(Cn) == 5
    plot(time*1000, X_rotation(1:length(time))*180/pi, '-', 'color', color_diamond, 'linewidth', 2);
else
    plot(time*1000, X_rotation(1:length(time))*180/pi, 'r-', 'linewidth', 2);
end

xlabel('\rm{\fontname{Times new roman}\rm{Time }\rm{[ms]}}');
ylabel('\rm{\fontname{Times new roman}\rm{Rotation angle }{\it{x}}\rm{-axis [бу]}}');
set(gca, 'FontName', 'Times new roman');
% set(get(gca,'XLabel'),'FontSize',figure_FontSize);
% set(get(gca,'YLabel'),'FontSize',figure_FontSize);
% set(get(gca,'ZLabel'),'FontSize',figure_FontSize);
% set(findobj('Fontsize',10),'fontsize', figure_FontSize-5);
% set(gcf,'unit','centimeters','position',[10 3 9 9])

hold off;

% 2D: Y orientation

cla(handles.axes8);
axes(handles.axes8);
% figure(8);
box on; hold on;   % visualizing the trace

if length(Cn) == 2
    plot(time*1000, Y_rotation(1:length(time))*180/pi, '-', 'color', color_sphere, 'linewidth', 2);
elseif length(Cn) == 3
    plot(time*1000, Y_rotation(1:length(time))*180/pi, '-', 'color', color_elli, 'linewidth', 2);
elseif length(Cn) == 4
    plot(time*1000, Y_rotation(1:length(time))*180/pi, '-', 'color', color_cone, 'linewidth', 2);
elseif length(Cn) == 5
    plot(time*1000, Y_rotation(1:length(time))*180/pi, '-', 'color', color_diamond, 'linewidth', 2);
else
    plot(time*1000, Y_rotation(1:length(time))*180/pi, 'r-', 'linewidth', 2);
end

xlabel('\rm{\fontname{Times new roman}\rm{Time }\rm{[ms]}}');
ylabel('\rm{\fontname{Times new roman}\rm{Rotation angle }{\it{y}}\rm{-axis [бу]}}');
set(gca, 'FontName', 'Times new roman');
% set(get(gca,'XLabel'),'FontSize',figure_FontSize);
% set(get(gca,'YLabel'),'FontSize',figure_FontSize);
% set(get(gca,'ZLabel'),'FontSize',figure_FontSize);
% set(findobj('Fontsize',10),'fontsize', figure_FontSize-5);
% set(gcf,'unit','centimeters','position',[10 3 9 9])

hold off;

% 2D: Z orientation

cla(handles.axes9);
axes(handles.axes9);
% figure(9);
box on; hold on;   % visualizing the trace

if length(Cn) == 2
    plot(time*1000, Z_rotation(1:length(time))*180/pi, '-', 'color', color_sphere, 'linewidth', 2);
elseif length(Cn) == 3
    plot(time*1000, Z_rotation(1:length(time))*180/pi, '-', 'color', color_elli, 'linewidth', 2);
elseif length(Cn) == 4
    plot(time*1000, Z_rotation(1:length(time))*180/pi, '-', 'color', color_cone, 'linewidth', 2);
elseif length(Cn) == 5
    plot(time*1000, Z_rotation(1:length(time))*180/pi, '-', 'color', color_diamond, 'linewidth', 2);
else
    plot(time*1000, Z_rotation(1:length(time))*180/pi, 'r-', 'linewidth', 2);
end

xlabel('\rm{\fontname{Times new roman}\rm{Time }\rm{[ms]}}');
ylabel('\rm{\fontname{Times new roman}\rm{Rotation angle }{\it{z}}\rm{-axis [бу]}}');
set(gca, 'FontName', 'Times new roman');
% set(get(gca,'XLabel'),'FontSize',figure_FontSize);
% set(get(gca,'YLabel'),'FontSize',figure_FontSize);
% set(get(gca,'ZLabel'),'FontSize',figure_FontSize);
% set(findobj('Fontsize',10),'fontsize', figure_FontSize-5);
% set(gcf,'unit','centimeters','position',[10 3 9 9])

hold off;

X_position = X_position(1:length(time));
Y_position = Y_position(1:length(time));
Z_position = Z_position(1:length(time));
X_rotation = X_rotation(1:length(time));
Y_rotation = Y_rotation(1:length(time));
Z_rotation = Z_rotation(1:length(time));

%%

elapse_time2 = datestr(now,'yyyy-mm-dd HH:MM:SS');
disp(['Started Time: ', elapse_time1]);
disp(['Ending Time: ',elapse_time2]);

%%
