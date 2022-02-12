function visualize_transducers(transducer)
%%
% This function is used to visualize the transducers by cylinders in
% specific positions.
%%

parameters;
transducer_height = 0.003;
X2 = transducer;
X1 = X2;
X2(:, 3) = transducer(:, 3) - inter_dist;
X1(:, 3) = transducer(:, 3) - inter_dist - transducer_height;
color = [0.5, 0.5, 0.5];      %  grey color

for ii = 1 : transducer_number
    cylinder3(X1(ii,:), X2(ii,:), trans_radius, color);
end

%% 