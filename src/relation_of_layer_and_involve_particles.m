function [particle_involve, layer_range] = ...
            relation_of_layer_and_involve_particles(r_order)
%%
% this function is used to give the relations of corresponding "layer"
% involve how many particles inside from the probe particle to current
% "layer".
% Generally, when "layer" increase one, then the involving particle also
% increase one. 
%
% However, if the distances of some particles to probe particle are close,
% which these particles have same spherical distance to probe particle,
% then increase one "layer" will increase more than one particle!!
% Therefore, we can not simply regard the "layer number" equal to "particle
% number". We need a new relation of which "layer" involves how many
% "particles" so that function "rescattering_coefficient.m" can describe
% clearer the inner and outer coefficients.
%%

parameters;

delta_r = range_r_coeff / 10 * particle_radius / (grid_resolution);
minimum_effective_distance = 5*delta_r;

if minimum_effective_distance > particle_radius
    error('For multi-particle system, the visualizing resolution for pressure field is too low, please decrease the visualizing range!\n');
end

particle_involve = [];
layer_range = 0;
layer_involve = 0;
for ii = 1:length(r_order)-1
    layer_involve = layer_involve + 1;
    if (r_order(ii+1) - r_order(ii)) > minimum_effective_distance
        particle_involve = [particle_involve, layer_involve];
        layer_range = [layer_range, r_order(ii+1)];
    end
end
particle_involve = [particle_involve, length(r_order)];



%%