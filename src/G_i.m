function  Gxy = G_i(rr, trans_radius, fluid_c, freq, Pr_0, inter_dist, Dims_Fourier) 
%%
% This function is used to return a transerve shape function 'Gxy'. for
% further Fourier expansion into a plain wave serier.
% piston-like wave refers to Pierce@1999@Book Eqs. (5-5.1) and (5-5.4)) and
% Cheeke@2002@(6.32)
%
%
% The 'p_i' will be used in "beam_shape_coeff.m" for the numerical
% beam-shape coefficient. However, for the transducer scenarios, as the
% beam-shape coefficient of the wave functions can not be numerical
% decomposite independently with the position vector vec{r}. So, for the
% wave functions other than 'plain', 'zero-Bessel' and 'non-zero-Bessel'
% types, we have to use Fourier expansion to simplify the wave functions as
% a series of plain waves, and each of them can be express by
% position-independent beam-shape coefficient. 
% Therefore, for complex wave function, the below codes are useless, while
% a new counterpart function "G_i.m" is used to obtain the transerve shape
% function for further Fourier expansion.
%
%
% trans_radius = 5 * 10^-3;       % radius of transducer is 5mm
% cos_beta = sqrt(1 - (1/(fluid_k * trans_radius))^2);
% 
%
% rr = (xx,yy) or Dims_Fourier == 2
% G_i(xx,yy):                                 0     y       (vertical view of OXYZ system)
%   [(xx1, yy1) (xx1, yy2) ... (xx1, yyn)       |！！！>
%    (xx2, yy1) (xx2, yy2) ... (xx2, yyn)     x |
%         .                                     v
%         .
%         .
%    (xxn, yy1) (xxn, yy2) ... (xxn, yyn)]
%
%
% rr = (xx,yy,zz) or Dims_Fourier == 3
% G_i(xx,yy,zz==1):                                      0     y       (vertical view of OXYZ system)
%   [(xx1, yy1, 1) (xx1, yy2, 1) ... (xx1, yyn, 1)       |！！！>
%    (xx2, yy1, 1) (xx2, yy2, 1) ... (xx2, yyn, 1)     x |
%         .                                              v
%         .
%         .
%    (xxn, yy1, 1) (xxn, yy2, 1) ... (xxn, yyn), 1]
%
%       . . . . . . 
%
% G_i(xx,yy,zz==end):                                          0     y       (vertical view of OXYZ system)
%   [(xx1, yy1, end) (xx1, yy2, end) ... (xx1, yyn, end)       |！！！>      (Z-axis from inside to outside)  
%    (xx2, yy1, end) (xx2, yy2, end) ... (xx2, yyn, end)     x |
%         .                                                    v
%         .
%         .
%    (xxn, yy1, end) (xxn, yy2, end) ... (xxn, yyn, end)]
%
%
% NOTE: 'Gxy' should be unit pressure (pa), because the initial pressure
%       'p_0' have been involved in the calculation for beam-shape
%       coefficient in "beam_shape_coeff.m" or (phase array) Equivalent
%       beam-shape coefficient in "phase_array_beam_shape_coeff.m".
%
%%

if Dims_Fourier == 2
    
    % transerve shape function (bessel-like beam)
    Gxy = 1 * ...          
        besselj(0, 2.412 * sqrt((rr(1)).^2 + (rr(2)).^2) / trans_radius) .* ...
        step_fun(sqrt((rr(1)).^2 + (rr(2)).^2) / trans_radius);

elseif Dims_Fourier == 3
    
    fluid_k = 2*pi*freq / fluid_c;
    
    rr(1) = rr(1) + 0;
    rr(2) = rr(2) + 0;
    rr(3) = rr(3) + inter_dist;       % invalid this line for polar figure
    
    rr_dist = sqrt(rr(1)^2 + rr(2)^2 + rr(3)^2);
    theta = acos(rr(3)/rr_dist);                % theta = [0, pi];
    % spherical piston-like wave: refer to Cheeke@2002@(6.32)
%     Gxy = Pr_0 * 1 * ...
%         besselj(0, fluid_k * trans_radius * sin(theta)) * ...
%         exp(1i * fluid_k * (rr_dist)) / (rr_dist);
    Gxy = Pr_0 * 2 * besselj(1, fluid_k * trans_radius * sin(theta)) / ...
        (fluid_k * trans_radius * sqrt(rr(1)^2 + rr(2)^2 + eps)) * ...
        exp(1i * fluid_k * (rr_dist));
%     Gxy = Pr_0 * 2 * besselj(1, fluid_k * trans_radius * sin(theta+eps)) / ...
%         (fluid_k * trans_radius * sqrt(rr(1)^2 + rr(2)^2 + eps)); 
%     Gxy = rr(1) + rr(2) + rr(3);  % test

end


%% ================================================================
% 
% % parameters;             % for transducers' information, and 'fluid_c'
% % if (strcmp(wave_type, 'single_transducer') || strcmp(wave_type, 'phase_array_transducer')) ~= 1
% %     error('Fourier Expansion is useful for wave types of ''single_transducer'' and ''phase_array_transducer''.\n');
% % end
% 
% %% Single transducer incident wave
% 
% % Gxy = 2000 * ...          % transerve shape function
% %     besselj(0, 2.5 * sqrt((xx + 0).^2 + (yy + 0).^2) / 0.005) .* ...
% %     step_fun(sqrt((xx + 0).^2 + (yy + 0).^2) / 0.005);
% 
% if strcmp(wave_type, 'single_transducer') == 1
%     if strcmp(direction, 'X') == 1          % x-direaction.
%         error('x: please setting the wave propagates along +z-axes direction.\n');      
%     end
%     if strcmp(direction, 'Y') == 1          % y-direaction.
%         error('y: please setting the wave propagates along +z-axes direction.\n');
%     end
%     if strcmp(direction, 'Z') == 1          % z-direaction.
%         Gxy = 1 * ...          % transerve shape function
%             besselj(0, 2.5 * sqrt((xx).^2 + (yy).^2) / trans_radius) .* ...
%             step_fun(sqrt((xx).^2 + (yy).^2) / trans_radius);
%     end
% end
% 
% %% Phase array incident wave
% 
% % trans_radius = 5 * 10^-3;       % radius of transducer is 5mm
% % cos_beta = sqrt(1 - (1/(fluid_k * trans_radius))^2);
% 
% if strcmp(wave_type, 'phase_array_transducer') == 1
%     if strcmp(direction, 'X') == 1          % x-direaction.
%         error('x: please setting the wave propagates along +z-axes direction.\n');      
%     end
%     if strcmp(direction, 'Y') == 1          % y-direaction.
%         error('y: please setting the wave propagates along +z-axes direction.\n');
%     end
%     if strcmp(direction, 'Z') == 1          % z-direaction.
%         Gxy = 1 * ...           % transerve shape function
%             besselj(0, 2.5 * sqrt((xx).^2 + (yy).^2) / trans_radius) .* ...
%             step_fun(sqrt((xx).^2 + (yy).^2) / trans_radius);
%     end
% end

%%

