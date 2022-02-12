function arrow3(X1, V1, r1, r2, h, color)
%%
% Plot a arrow graph, similar to quiver3, while the arrow is much better in
% visualization.
% given the vector position 'X1';
% given the vector direction 'V1';
% the radius of the arrow 'r1';
% the radius of the cone 'r2', r2 > r1;
% the height of the cone 'h';
% the color of the arrow.
%%

X2 = X1 + V1;
hold on;
cylinder3(X1, X2, r1, color);
X3 = (X2-X1) / norm(X2-X1) * h + X2;
cone3(X2,X3,r2,color);

