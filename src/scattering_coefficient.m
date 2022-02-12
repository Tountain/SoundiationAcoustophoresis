function [sn, tn] = scattering_coefficient(n)
%%
% this function calculate the scattering coefficient "sn" under
% different Boundary Conditions.
%%

parameters;

%% B.C. 1: rigid particles. Eq.(4) at G.T.S@2011@IEEE

if strcmp(BC, 'rigid') == 1
    
    if multi_particle == 0      % single particle system
        fluid_ka = particle_radius*fluid_k;

        sn = - (sBessel_partial(n,fluid_ka,1)) / (sHankel_partial(n,fluid_ka,1));
        tn = 0;
    else                        % multi-particle system
        particles_Cartesian_data;
        sn = zeros(1, particle_number);
        tn = zeros(1, particle_number);
        
        for ii = 1:particle_number
            fluid_ka_ii = multi_particle_radius(ii)*fluid_k;
            sn(ii) = - (sBessel_partial(n,fluid_ka_ii,1)) / (sHankel_partial(n,fluid_ka_ii,1));
            tn(ii) = 0;
        end
    end
    
end

%% B.C. 2: compressible particles. Eq.(7) at G.T.S@2014@JASA

if strcmp(BC, 'rigid') ~= 1
    
    if multi_particle == 0      % single particle system
        fluid_ka = particle_radius*fluid_k;
        particle_ka = particle_radius*particle_k;
        numerator = sBessel_partial(n,fluid_ka,1)*sBessel(n,particle_ka,1) - ...
                    gama*sBessel(n,fluid_ka,1)*sBessel_partial(n,particle_ka,1);
        demoninator = gama*sHankel(n,fluid_ka,1)*sBessel_partial(n,particle_ka,1) - ...
                      sHankel_partial(n,fluid_ka,1)*sBessel(n,particle_ka,1);

        sn = numerator / demoninator;
        tn = (sn * sHankel(n,fluid_ka,1) + sBessel(n,fluid_ka,1)) / sBessel(n,particle_ka,1);
    else                        % multi-particle system
        particles_Cartesian_data;
        sn = zeros(1, particle_number);
        tn = zeros(1, particle_number);
        
        for ii = 1:particle_number
            fluid_ka_ii = multi_particle_radius(ii)*fluid_k;
            particle_ka_ii = multi_particle_radius(ii)*multi_particle_k(ii);
            numerator_ii = sBessel_partial(n,fluid_ka_ii,1)*sBessel(n,particle_ka_ii,1) - ...
                           multi_gama(ii)*sBessel(n,fluid_ka_ii,1)*sBessel_partial(n,particle_ka_ii,1);
            demoninator_ii = multi_gama(ii)*sHankel(n,fluid_ka_ii,1)*sBessel_partial(n,particle_ka_ii,1) - ...
                             sHankel_partial(n,fluid_ka_ii,1)*sBessel(n,particle_ka_ii,1);

            sn(ii) = numerator_ii / demoninator_ii;
            tn(ii) = (sn(ii) * sHankel(n,fluid_ka_ii,1) + sBessel(n,fluid_ka_ii,1)) / sBessel(n,particle_ka_ii,1);
        end
    end
    
end

%% B.C. 3:

