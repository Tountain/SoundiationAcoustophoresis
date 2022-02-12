function [p_0] = proba_particle_input_amplitude(p_inlet, fluid, fluid_k, particle_radius, range_r_coeff)
%%
% As the medium may exist the attenuation effects and we "solving frame" is
% fixed on the center of proba particle. If we do not adjust the proba
% particle input amplitude "p_0", then due to the attenuation effects, the
% microphones' input amplitude can not control as same value as no
% attenuation scenario. Therefore, we need to correct the proba particle
% input amplitude "p_0" so that keep consistent of microphones' input. 
% 
% sub-function of "parameters.m"
% correction principles can find from my derivation of Point 7.
%%

% 'mixture' means exist attenuation, then we need to make sure the inlet
% amplitude of incident wave is equal to "p_0", instead of particle's
% center position.
% if strcmp(fluid, 'mixture') == 1
%     if strcmp(direction, 'XZ') == 1
%         p_0x = p_inlet * exp(-abs(imag(fluid_k)) * (particle_radius * range_r_coeff));% + deviationX));
%         p_0z = p_inlet * exp(-abs(imag(fluid_k)) * (particle_radius * range_r_coeff));% + deviationZ));
%     else
%         p_0 = p_inlet * exp(-abs(imag(fluid_k)) * (particle_radius * range_r_coeff));% + deviation));
%     end
% else
%     if strcmp(direction, 'XZ') == 1
%         p_0x = p_inlet;
%         p_0z = p_inlet;
%     else
%         p_0 = p_inlet;
%     end
% end

if strcmp(fluid, 'mixture') == 1
    p_0 = p_inlet * exp(-abs(imag(fluid_k)) * (particle_radius * range_r_coeff));
    p_0;
else
    p_0 = p_inlet;
end