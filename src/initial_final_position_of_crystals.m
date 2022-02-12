function initial_final_position_of_crystals(fig_num, probe_initial_position)
%%
% this function is used to visualize the initial and final position of the
% crystals.
%%

parameters;
particles_Cartesian_data;

global X_vec Y_vec Z_vec time
%time = 1.9;
%% determine the coordinate of initial and final position

distance = velocity * time;
initial_position_x = particle(2:end, 1) - (-X_vec) - distance(1);
initial_position_y = particle(2:end, 2) - (-Y_vec) - distance(2);
initial_position_z = particle(2:end, 3) - (-Z_vec) - distance(3);
initial_position = [initial_position_x, initial_position_y, initial_position_z];
final_position_x = particle(2:end, 1) - (-X_vec);
final_position_y = particle(2:end, 2) - (-Y_vec);
final_position_z = particle(2:end, 3) - (-Z_vec);
final_position = [final_position_x, final_position_y, final_position_z];

%% graph initial and final position of the cryticals and the probe particle

figure(fig_num);

hold on;
grid on;

% plot3(initial_position_x, initial_position_y, initial_position_z,'ro','MarkerFaceColor','k');
% plot3(final_position_x, final_position_y, final_position_z,'ro','MarkerFaceColor','k');

CO(:,:,1) = zeros(25); % color control
for ii = 2:particle_number
    [x y z]=sphere(25);
    surf(multi_particle_radius(ii)*x+(initial_position_x(ii-1)), ...
         multi_particle_radius(ii)*y+(initial_position_y(ii-1)), ...
         multi_particle_radius(ii)*z+(initial_position_z(ii-1)), CO); % ?????radius??
    surf(multi_particle_radius(ii)*x+(final_position_x(ii-1)), ...
         multi_particle_radius(ii)*y+(final_position_y(ii-1)), ...
         multi_particle_radius(ii)*z+(final_position_z(ii-1)), CO); % ?????radius??
end

[x y z]=sphere(10);
surf(particle_radius*x+(probe_initial_position(1)), ...
     particle_radius*y+(probe_initial_position(2)), ...
     particle_radius*z+(probe_initial_position(3))); % ?????radius??
surf(particle_radius*x+(X_vec), ...
     particle_radius*y+(Y_vec), ...
     particle_radius*z+(Z_vec)); % ?????radius??

plot3([0, distance(1)],[0, distance(2)],[0, distance(3)], 'm', 'linewidth', 2);
plot3(distance(1),distance(2),distance(3),'mp','MarkerSize',10,'MarkerFaceColor','k');
%annotation('arrow',initial_position,final_position,'LineStyle','-','color',[1 0 0]);
%annotation('arrow',initial_position,final_position,'LineStyle','-','color',[1 0 0],'HeadStyle','plain');
%annotation('arrow',initial_position,final_position,'LineStyle','-','color',[1 0 0],'HeadStyle','cback3');
%annotation('doublearrow',initial_position,final_position,'LineStyle','-','color',[1 0 0],'HeadStyle','cback3');
