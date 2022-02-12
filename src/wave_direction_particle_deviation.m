function [dir_sign, x_translation, y_translation, z_translation] = ...
    wave_direction_particle_deviation(wave_type, direction, symbol, deviationX, deviationY, deviationZ)
%%
% this function is used to reture the direction of incident wave and the
% probe particle's deviation from origin position under the "Frame coords
% (microphones' coords)".
%
% sub-function of "parameters.m"
%%

if (strcmp(wave_type, 'plain') == 1) || (strcmp(wave_type, 'standing_plain') == 1)
    if strcmp(direction, 'X') == 1
        if strcmp(symbol, 'positive') == 1
            dir_sign = 1;
        elseif strcmp(symbol, 'negative') == 1
            dir_sign = -1;
        end
        x_translation =  deviationX;     % X position of particle relative to incident wave
        y_translation =  0;             % Y position of particle relative to incident wave
        z_translation =  0;             % Z position of particle relative to incident wave
    end
    if strcmp(direction, 'Y') == 1
        if strcmp(symbol, 'positive') == 1
            dir_sign = 1;
        elseif strcmp(symbol, 'negative') == 1
            dir_sign = -1;
        end
        x_translation =  0;             % X position of particle relative to incident wave
        y_translation =  deviationY;     % Y position of particle relative to incident wave
        z_translation =  0;             % Z position of particle relative to incident wave
    end
    if strcmp(direction, 'Z') == 1
        if strcmp(symbol, 'positive') == 1 || strcmp(symbol, 'arbitrary') == 1
            dir_sign = 1;
        elseif strcmp(symbol, 'negative') == 1
            dir_sign = -1;
        end
        x_translation =  0;             % X position of particle relative to incident wave
        y_translation =  0;             % Y position of particle relative to incident wave
        z_translation =  deviationZ;     % Z position of particle relative to incident wave
    end
    if strcmp(direction, 'XZ') == 1
        dir_sign = 1;
        x_translation =  deviationX;    % X position of particle relative to incident wave
        y_translation =  0;             % Y position of particle relative to incident wave
        z_translation =  deviationZ;    % Z position of particle relative to incident wave
    end
end

if (strcmp(wave_type, 'zero-Bessel') == 1)
    if strcmp(symbol, 'positive') == 1
        dir_sign = 1;
    elseif strcmp(symbol, 'negative') == 1
        dir_sign = -1;
    end
    x_translation =  deviationX;         % X position of particle relative to incident wave
    y_translation =  deviationY;         % Y position of particle relative to incident wave
    z_translation =  deviationZ;         % Z position of particle relative to incident wave
end

if (strcmp(wave_type, 'non-zero-Bessel') == 1)
    if strcmp(symbol, 'positive') == 1
        dir_sign = 1;
    elseif strcmp(symbol, 'negative') == 1
        dir_sign = -1;
    end
    x_translation =  0;         % X position of particle relative to incident wave
    y_translation =  0;         % Y position of particle relative to incident wave
    z_translation =  0;         % Z position of particle relative to incident wave
end

% circular piston wave function
if (strcmp(wave_type, 'single_transducer') == 1 || strcmp(wave_type, 'single_transducer_standing') == 1)
    if strcmp(symbol, 'positive') == 1
        dir_sign = 1;
    elseif strcmp(symbol, 'negative') == 1
        dir_sign = -1;
    end
    x_translation =  deviationX;         % X position of particle relative to incident wave
    y_translation =  deviationY;         % Y position of particle relative to incident wave
    z_translation =  deviationZ;         % Z position of particle relative to incident wave
end

% circular piston wave function
if (strcmp(wave_type, 'phase_array_transducer') == 1 || strcmp(wave_type, 'phase_array_transducer2') == 1 || strcmp(wave_type, 'phase_array_transducer_standing') == 1 )
    if strcmp(symbol, 'positive') == 1
        dir_sign = 1;
    elseif strcmp(symbol, 'negative') == 1
        dir_sign = -1;
    end
    x_translation =  deviationX;         % X position of particle relative to incident wave
    y_translation =  deviationY;         % Y position of particle relative to incident wave
    z_translation =  deviationZ;         % Z position of particle relative to incident wave
end

%%